"""
      size_inner_tank(fuse_tank, t_cond::Vector{Float64}, ρfuel::Float64,
                        Rfuse::Float64, dRfuse::Float64, wfb::Float64, nfweb::Float64,
                        Wfuel::Float64)

`size_inner_tank` calculates the weight of the cryogenic fuel tank for a LH-fueled aircraft.

!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `fuse_tank::Struct`: structure with tank parameters.
      - `t_cond::Float64`: Vector with tank isulation layer thickness. Provided separately from fuse_tank as it changes during 
      non-linear solve process.

      **Outputs:**
      - `Wtank_total::Float64`: Total tank weight including fuel (N).
      - `l_cyl::Float64`: Length of the cylindrical portion (m).
      - `tskin::Float64`: Thickness of the tank's skin (m).
      - `Rtank_outer::Float64`: Outer radius of the tank (m).
      - `Vfuel::Float64`: Volume of fuel (m^3).
      - `Wtank::Float64`: Weight of the empty tank (N).
      - `Wfuel_tot::Float64`: Total weight of fuel (N).
      - `Winsul_sum::Float64`: Sum of insulation weight (N).
      - `t_head::Float64`: Thickness of the tank's head (m).
      - `Whead::Float64`: Weight of the tank's head (N).
      - `Wcyl::Float64`: Weight of the tank's cylinder (N).
      - `Wstiff::Float64`: Total stiffener weight (N)
      - `Winsul::Float64`: Weight of insulation (N).
      - `Shead_insul::Float64`: Insulated surface area of the head (m^2).
      - `l_tank::Float64`: Total length of the tank (m).

NOTE: Al alloy 2219 has been recommended as tank material (from H2 tank paper in OneNote)

See [here](@ref fueltanks).
"""
function size_inner_tank(fuse_tank, t_cond::Vector{Float64})

      #Unpack parameters in fuse_tank
      Rfuse = fuse_tank.Rfuse
      dRfuse = fuse_tank.dRfuse
      wfb = fuse_tank.wfb
      nfweb = fuse_tank.nfweb

      Wfuel = fuse_tank.Wfuelintank

      ρfuel = fuse_tank.rhofuel
      ρfgas = fuse_tank.rhofuelgas
      ftankadd = fuse_tank.ftankadd
      Δp = fuse_tank.ptank
      sigskin = fuse_tank.inner_material.UTS
      material_insul = fuse_tank.material_insul
      rhoskin = fuse_tank.inner_material.ρ
      clearance_fuse = fuse_tank.clearance_fuse
      AR = fuse_tank.ARtank
      ullage_frac = fuse_tank.ullage_frac
      weld_eff = fuse_tank.ew
      θ = fuse_tank.theta_inner

# Total thickness:
      thickness_insul = sum(t_cond)

      s_a = sigskin / 4 #Maximum allowable stress is 1/4 Ultimate tensile strength (Barron 1985, p. 359)
      
      Rtank_outer = Rfuse - thickness_insul - clearance_fuse

      tskin = Δp * (2 * Rtank_outer) / (2 * s_a * weld_eff + 0.8 * Δp) #(7.1) in Barron (1985)

      Rtank = Rtank_outer - tskin
      #tfweb = 2.0 * Δp * wfb  / ew
      Lhead = Rtank / AR       # eg. for a 2:1 ellipsoid majorax/minorax = 2/1 ⟹ R/Lhead = 2/1 
      
      K = (1/6) * (AR^2 + 2) # Aspect ratio of 2:1 for the head (# Barron pg 359) 
      t_head = Δp* (2*Rtank_outer) * K/ (2 * s_a * weld_eff + 2 * Δp * (K - 0.1)) #(7.2) in Barron (1985)

#--- Calculate length of cylindrical portion
      Wfuel_tot = Wfuel #Wfuel already includes the amount that boils off
      ρfmix = (1 - ullage_frac) * ρfuel + ullage_frac * ρfgas #Density of saturated mixture in tank
      Vfuel = Wfuel_tot / (gee * ρfmix) #Total tank volume taken by saturated mixture
      Vinternal = Vfuel  # required internal volume
      V_ellipsoid = 2π * (Rtank^3 / AR) / 3  # Half the vol of std ellipsoid = 1/2×(4π/3 ×(abc)) where a,b,c are the semi-axes length. Here a = R/AR, b=c=R
                                       # Also see: https://neutrium.net/equipment/volume-and-wetted-area-of-partially-filled-horizontal-vessels/
      V_cylinder = Vinternal - 2*V_ellipsoid
      l_cyl = V_cylinder / (π * (Rtank^2)) #required length of cylindrical portion

#--- tank cross-section geometric parameters
      wfblim = max( min( wfb , Rtank) , 0.0 )
      thetafb = asin(wfblim / Rtank)

#--- areas
      Scyl = (2.0*π+4.0*nfweb*thetafb)*Rtank*l_cyl + 2.0*dRfuse*l_cyl # Surface area of cylindrical part
      #Afweb = nfweb*(2.0*hfb+dRfuse)*tfweb
      #Atank = (π + nfweb*(2.0*thetafb + sin2t))*Rtank^2 + 2.0*Rtank*dRfuse #+ 2.0*(Rtank+nfweb*wfb)*dRfuse
      Shead = (2.0*π + 4.0*nfweb*thetafb)*Rtank^2* ( 0.333 + 0.667*(Lhead/Rtank)^1.6 )^0.625 # This form is better for insul thickness 
                                                                                          # but just as a note to reader this comes from  semi- oblate spheroid surf area is ≈ 2π×R²[1/3 + 2/3×(1/AR)^1.6]^(1/1.6)
#--- component volumes
      Vcyl  = Scyl*tskin    # volume of the metal in the cylindrical part
      #Vhead = Shead*tskin
      Vhead = Shead * t_head # volume of head

      Sinternal = Scyl + 2 * Shead

#--- weights and weight moments
      Whead = rhoskin*gee*Vhead
      Wcyl  = rhoskin*gee*Vcyl
      Winnertank = Wcyl + 2*Whead + Wfuel #Weight of inner vessel without stiffeners and supports, inc. fuel

      #stiffeners
      Nmain = 2 #Tanks typically have two main support rings
      Wmainstiff = stiffener_weight("inner", Winnertank / Nmain, Rtank_outer, 
                                s_a, rhoskin, θ) #Weight of one main stiffener

      Wstiff = Nmain * Wmainstiff
      Wtank = (Wcyl + 2*Whead + Wstiff) * (1 + ftankadd)

#--- insulation weight!
      N = length(t_cond)
      Vcyl_insul = zeros(N)
      Winsul = zeros(N)
      Shead_insul = zeros(N + 1) #add one for first (tank wall) surface 
      Vhead_insul = zeros(N)
      rho_insul = zeros(N)
      L = Lhead + tskin

      #Assemble array with layer densities
      for i = 1:N
            rho_insul[i] = insulation_density_calc(material_insul[i])
      end

      Ro = Ri = Rtank_outer # Start calculating insulation from the outer wall of the metal tank ∴Ri of insul = outer R of tank
      Shead_insul[1] = (2.0*π + 4.0*nfweb*thetafb)*(Ro)^2* ( 0.333 + 0.667*(L/Ro)^1.6 )^0.625
      
      for n in 1:N
            
            Ro = Ro + t_cond[n]
            L  = L  + t_cond[n]
            # println("AR ≈ $(Ro/L)")
            Vcyl_insul[n]  = (π * ( Ro^2 - Ri^2 ) * l_cyl)
            Shead_insul[n+1] = (2.0*π + 4.0*nfweb*thetafb)*(Ro)^2* ( 0.333 + 0.667*(L/Ro)^1.6 )^0.625

            Area_coeff = Shead_insul[n+1] / Ro^2 #coefficient that relates area and radius squared
            Vhead_insul[n] = ((Shead_insul[n] + Shead_insul[n+1])/2 - Area_coeff/(6) * t_cond[n]^2) * t_cond[n] #Closed-form solution
            
            Winsul[n] = (Vcyl_insul[n] + 2*Vhead_insul[n]) * rho_insul[n] * gee
            # println("AR = $(Ro/L)")
            Ri = Ro
      end

      Winsul_sum = sum(Winsul)
      Wtank = (Wtank + Winsul_sum)
      l_tank = l_cyl + 2*Lhead + 2*thickness_insul + 2*t_head #Total longitudinal length of the tank
#--- overall tank weight
      Wtank_total = Wtank + Wfuel_tot

return  Wtank_total, l_cyl, tskin, Rtank_outer, Vfuel, Wtank, Wfuel_tot, Winsul_sum, t_head, Whead, Wcyl, Wstiff, Winsul, Sinternal, Shead_insul, l_tank
end


"""
    size_outer_tank(fuse_tank, Winnertank, l_cyl, Ninterm)
This function sizes the outer vessel and calculates the weights of its components.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `fuse_tank::Struct`: structure with tank parameters.
    - `Winnertank::Float64`: weight of inner vessel and contents (N).
    - `l_cyl::Float64`: length of cylindrical portion of outer vessel (m).
    - `Ninterm::Float64`: optimum number of intermediate stiffener rings.

    **Outputs:**
    - `Wtank::Float64`: total weight of outer vessel (N).
    - `Wcyl::Float64`: weight of cylindrical portion of outer vessel (N).
    - `Whead::Float64`: weight of one elliptical outer-tank head (N).
    - `Wstiff::Float64`: total weight of stiffener material (N).
    - `S_outer::Float64`: surface area of outer vessel (m^2).
    - `Shead::Float64`: surface area of one outer vessel head (m^2).
    - `Scyl::Float64`: surface area of cylindrical portion of tank (m^2).
    - `t_cyl::Float64`: wall thickness of cylindrical portion of tank (m).
    - `t_head::Float64`: wall thickness of tank head (m). 
    - `l_tank::Float64`: Total length of the tank (m).
"""
function size_outer_tank(fuse_tank, Winnertank::Float64, l_cyl::Float64, Ninterm::Float64)
      #Unpack parameters in fuse_tank
      poiss = fuse_tank.inner_material.ν
      Eouter = fuse_tank.inner_material.E
      ρouter = fuse_tank.inner_material.ρ
      UTSouter = fuse_tank.inner_material.UTS
      ftankadd = fuse_tank.ftankadd
      wfb = fuse_tank.wfb
      nfweb = fuse_tank.nfweb
      ARtank = fuse_tank.ARtank
      θ_outer = fuse_tank.theta_outer

      θ1 = θ_outer[1]
      θ2 = θ_outer[2]

      Nmain = 2 #Tanks typically have two main support rings
      pc = 4 * pref #4*p_atm; Collapsing pressure, Eq. (7.11) in Barron (1985)
      s_a = UTSouter / 4

      #Calculate outer vessel geometry
      Rtank_outer = fuse_tank.Rfuse - fuse_tank.clearance_fuse
      Do = 2 * Rtank_outer #outside diameter

      Nstiff = Nmain + Ninterm #Total number of stiffeners
      L = l_cyl / (Nstiff - 1) #There are two stiffeners at the ends, so effective number of sections in skin is N - 1
      L_Do = L / Do

      #Find cylinder wall thickness. This applies to a short cylinder.
      pressure_res(t_D) = 2.42*Eouter*(t_D)^(5/2) / ( (1 - poiss^2)^(3/4) * (L_Do - 0.45*sqrt(t_D)) ) - pc
      t_Do = find_zero(pressure_res, 1e-3) #Find root with Roots.jl
      t_cyl = t_Do * Do

      #Find head wall thickness
      if ARtank == 2.0
            K1 = 0.90# See table 7.6 for D/D1=2.0 in Barron p. 367
      elseif ARtank == 1.0
            K1 = 0.50
      else  
            println("ARtank of heads not supported, see size_outer_tank()")
            K1=1.0
      end
      t_head = K1 * Do * sqrt(pc * sqrt(3*(1 - poiss^2))/ (0.5*Eouter))

      ## Areas
      wfblim = max( min( wfb , Rtank_outer) , 0.0 )
      thetafb = asin(wfblim / Rtank_outer)

      Shead = (2.0*π + 4.0*nfweb*thetafb)*(Rtank_outer)^2* (0.333 + 0.667*(1/ARtank)^1.6 )^0.625
      Scyl  = 2π*Rtank_outer*l_cyl  # Cross-sectional area

      Souter = Scyl + 2*Shead

      ## Volume and Weight
      Vcyl  = Scyl*t_cyl
      Vhead = Shead*t_head

      Wcyl  = Vcyl*ρouter*gee
      Whead =  Vhead*ρouter*gee

      Wtank_no_stiff = Wcyl + 2 * Whead

      # Size stiffeners
      tanktype = "outer"

      Wmainstiff = stiffener_weight(tanktype, Winnertank / Nmain, Rtank_outer, #Weight of one main stiffener, each one 
                                    s_a, ρouter, θ1, θ2, Nstiff, l_cyl, Eouter)  #carries half the inner vessel load
                                                                              
      Wintermstiff = stiffener_weight(tanktype, 0.0, Rtank_outer, 
                                    s_a, ρouter, θ1, θ2, Nstiff, l_cyl, Eouter) #Weight of one intermediate stiffener, which carries no load

      Wstiff = Nmain * Wmainstiff + Ninterm * Wintermstiff #Total stiffener weight

      Wtank = (Wtank_no_stiff + Wstiff) * (1 + ftankadd) #Find total tank weight, including additional mass factor

      l_outer = l_cyl + Do / ARtank + 2*t_head #Total length of outer vessel

      return Wtank, Wcyl, Whead, Wstiff, Souter, Shead, Scyl, t_cyl, t_head, l_outer
end

"""
    stiffener_weight(tanktype, W, Rtank, s_a, ρstiff, θ1, θ2 = 0.0, Nstiff = 2.0, l_cyl = 0, E = 0)
This function calculates the weight of a single stiffener in an inner or outer vessel for a given inner vessel weight.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `tanktype::String`: type of tank, options are "inner" or "outer".
    - `W::Float64`: load carried by a stiffener ring (N).
    - `Rtank::Float64`: tank radius (m).
    - `s_a::Float64`: maximum allowable stress in stiffener material (Pa).
    - `ρstiff::Float64`: stiffener density (kg/m^3).
    - `θ1::Float64`: angular position of bottom tank supports, measured from the bottom of the tank (rad).
    - `θ2::Float64`: angular position of top tank supports, measured from the bottom of the tank (rad). Only used with "outer" tank.
    - `Nstiff::Float64`: total number of stiffeners on outer vessel. Only used with "outer" tank.
    - `l_cyl::Float64`: length of cylindrical portion of tank (m). Only used with "outer" tank.
    - `E::Float64`: Young's modulus of stiffener material (Pa). Only used with "outer" tank.

    **Outputs:**
    - `Wstiff::Float64`: weight of a stiffener ring (N).
"""
function stiffener_weight(tanktype::String, W::Float64, Rtank::Float64, s_a::Float64, 
      ρstiff::Float64, θ1::Float64, θ2::Float64 = 0.0, Nstiff::Float64 = 2.0, l_cyl::Float64 = 0.0, E::Float64 = 0.0)
    
      if tanktype == "inner" 
            _, kmax = stiffeners_bendingM(θ1) #Find k = 2πM/(WR)
            Icollapse = 0 #inner vessel cannot collapse as it is pressurized

      elseif tanktype == "outer"
            _, kmax = stiffeners_bendingM_outer(θ1, θ2) #Find k = 2πM/(WR)
            pc = 4 * pref #Critical pressure is 4 times atmospheric pressure, Eq. (7.11) in Barron (1985)
            Do = 2 * Rtank #outer diameter

            L = l_cyl/ (Nstiff - 1) #Length of portion between supports
            Icollapse = pc * Do^3 * L / (24 * E) #Second moment of area needed to avoid collapse
      end

      Mmax = kmax * W * Rtank / (2π) #Maximum bending moment due to loads
      Z = Mmax / s_a #required section modulus to withstand bending loads

      #Assume sectional properties of a 100 x 100 I-beam
      W = 100e-3 #Flange width
      t_w = 7.1e-3 #Web thickness
      t_f = 8.8e-3 #Flange thickness

      #For an I-beam, I > W * H^2 * t_f / 2 + t_f^3 * W / 6
      #The required second moment of area is I = Icollapse + Z * (H/2 + t_f/2)
      #Find beam height by solving W * H^2 * t_f / 2 + t_f^3 * W / 6 = Icollapse + Z * (H/2 + t_f/2)

      a = t_f * W / 2 #Coefficients in quadratic equation
      b = -Z/2
      c = -1 * Icollapse - Z * t_f/2 + t_f^3 * W / 6

      H = (-b + sqrt(b^2 - 4 * a * c)) / (2 * a) #Solve quadratic eq.
      S = 2 * W * t_f + (H - t_f) * t_w #Beam cross-sectional area

      Wstiff = gee * ρstiff * S * 2π * Rtank #Weight of a single stiffener running along circumference
      return Wstiff
end

"""
    stiffeners_bendingM(θ)
This function can be used to calculate the bending moment distribution in a stiffener ring for an inner cryogenic tank.
It applies Eqs. (7.4) and (7.5) in Barron (1985) to find the bending moment distribution. The function returns the
maximum value of ``k = 2πM/(WR)`` on the ring's circumference.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `θ::Float64`: angular position of tank supports, measured from the bottom of the tank (rad).

    **Outputs:**
    - `ϕmax::Float64`: angular position of maximum bending moment on ring circumference (rad).
    - `kmax::Float64`: Maximum value of the ratio ``k = 2πM/(WR)`` on ring circumference.
"""
function stiffeners_bendingM(θ::Float64)
      ϕlist = LinRange(0.0, π, 361)
      k = zeros(length(ϕlist))

      for (i,ϕ) in enumerate(ϕlist)
            if 0 ≤ ϕ ≤ θ
                  k[i] = 0.5*cos(ϕ) + ϕ*sin(ϕ) - (π - θ)*sin(θ) + cos(θ) + cos(ϕ)*(sin(θ)^2)
            elseif θ ≤ ϕ ≤ π
                  k[i] = 0.5*cos(ϕ) - (π - ϕ)*sin(ϕ) + θ  + cos(θ) + cos(ϕ)*(sin(θ)^2)
                  #TODO this equation is Eq. (7.5) in Barron; however, this equation does not match the curves in Fig. 7.3
                  #Suspect error in Barron.
            end
      end
      kmax, imax = findmax(abs.(k))
      ϕmax = ϕlist[imax]
      return ϕmax, kmax
end

"""
    stiffeners_bendingM_outer(θ1,θ2)
This function can be used to calculate the bending moment distribution in a stiffener ring for an outer cryogenic tank.
It applies Eqs. (7.13)-(7.15) in Barron (1985) to find the bending moment distribution. The function returns the
maximum value of ``k = 2πM/(WR)`` on the ring's circumference.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `θ1::Float64`: angular position of bottom tank supports, measured from the bottom of the tank (rad).
    - `θ2::Float64`: angular position of top tank supports, measured from the bottom of the tank (rad).

    **Outputs:**
    - `ϕmax::Float64`: angular position of maximum bending moment on ring circumference (rad).
    - `kmax::Float64`: Maximum value of the ratio ``k = 2πM/(WR)`` on ring circumference.
"""
function stiffeners_bendingM_outer(θ1::Float64, θ2::Float64)
      ϕlist = LinRange(0.0, π, 361)
      k = zeros(length(ϕlist))

      for (i,ϕ) in enumerate(ϕlist)
            if 0 ≤ ϕ ≤ θ1
                  k[i] = cos(ϕ)*(sin(θ2)^2 - sin(θ1)^2 ) + (cos(θ2) - cos(θ1) ) +
                              - (π-θ2)*sin(θ2) +  (π-θ1)*sin(θ1)
            elseif θ1 ≤ ϕ ≤ θ2
                  k[i] = cos(ϕ)*(sin(θ2)^2 - sin(θ1)^2 ) + (cos(θ2) - cos(θ1) ) +
                              - (π-θ2)*sin(θ2) +  π*sin(ϕ) -θ1*sin(θ1)
            elseif θ2 ≤ ϕ ≤ π
                  k[i] = cos(ϕ)*(sin(θ2)^2 - sin(θ1)^2 ) + (cos(θ2) - cos(θ1) )  +
                              (θ2*sin(θ2) - θ1*sin(θ1)  )
            end
      end
      kmax, imax = findmax(abs.(k))
      ϕmax = ϕlist[imax]
      return ϕmax, kmax
end

"""
    optimize_outer_tank(fuse_tank, Winnertank, l_cyl, θ1, θ2)
This function optimizes the number of intermediate stiffener rings to minimize the weight of an outer vessel.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `fuse_tank::Struct`: structure with tank parameters.
    - `Winnertank::Float64`: weight of inner vessel and contents (N).
    - `l_cyl::Float64`: length of cylindrical portion of outer vessel (m).

    **Outputs:**
    - `Ninterm::Float64`: optimum number of intermediate stiffener rings.
"""
function optimize_outer_tank(fuse_tank, Winnertank::Float64, l_cyl::Float64)

      obj(x, grad) = size_outer_tank(fuse_tank, Winnertank, l_cyl, x[1])[1] #Minimize Wtank

      initial_x = [fuse_tank.Ninterm]

      #Set bounds
      lower = [0.0]
      upper = [50.0]

      #Use NLopt.jl to minimize function 
      opt = Opt(:LN_NELDERMEAD, length(initial_x))
      opt.lower_bounds = lower
      opt.upper_bounds = upper
      opt.ftol_rel = 1e-9
      opt.maxeval = 100  # Set the maximum number of function evaluations

      opt.min_objective = obj

      (minf,xopt,ret) = NLopt.optimize(opt, initial_x)
      Ninterm = xopt[1]
      return Ninterm
end

"""
      insulation_density_calc(material)

This function calculates the density of different insulation materials.
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `material::String`: material name.

      **Outputs:**
      - `ρ::Float64`: mass density (kg/m^3).
"""
function insulation_density_calc(material::String)
      if lowercase(material) == "rohacell41s"
            ρ = 35.0 #kg/m^3. From Brewer (1991)
      elseif lowercase(material) == "polyurethane27"
            ρ = 27.0 #kg/m^3
      elseif lowercase(material) == "polyurethane32"
            ρ = 32.0 #kg/m^3
      elseif lowercase(material) == "polyurethane35"
            ρ = 35.0 #kg/m^3
      elseif lowercase(material) == "mylar"
            ρ = 1390 #kg/m^3, https://www.matweb.com/search/datasheet_print.aspx?matguid=981d85aa72b0419bb4b26a3c06cb284d
      elseif lowercase(material) == "vacuum"
            ρ = 0 #kg/m^3
      else
            error("Insulation materials currently supported are
                  [rohacell41S, polyurethane27, polyurethane32, polyurethane35, vacuum],
                  but you supplied $material")
      end
      return ρ
end
