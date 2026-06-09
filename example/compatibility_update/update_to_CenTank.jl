"""
This script update the model saved from current branch to the CenTank branch readable format
"""

using TASOPT
using JLD2
include(__TASOPTindices__)

@kwdef mutable struct fuselage_tank_new
    """Fuel type name"""
    fueltype::String = ""
    """Fuel tank count"""
    tank_count::Int64 = 0
    """Fuel tank location"""
    placement::String = ""
    """Flag for insulation sizing"""
    sizes_insulation::Bool = false
    """Weight of fuel in one tank (N)"""
    Wfuelintank::Float64 = 0.0

    clearance_fuse::Float64 = 0.0

    """Vector with insulation layer thickness (m)"""
    t_insul::Vector{Float64} = Float64[]
    """Vector with insulation materials"""
    material_insul::Vector{ThermalInsulator} = Float64[]
    """Vector with insulation layer design indices"""
    iinsuldes::Vector{Int64} = Float64[]
    """Length of cylindrical portion of tank (m)"""
    l_cyl_inner::Float64 = 0.0
    """Length of inner tank (m)"""
    l_inner::Float64 = 0.0
    """Inner tank radius (m)"""
    Rinnertank::Float64 = 0.0
    """Vector with surface areas of insulation tank heads (m^2)"""
    Shead_insul::Vector{Float64} = Float64[]

    """Inner vessel material"""
    inner_material::StructuralAlloy = StructuralAlloy("Al-2219-T87")
    """Outer vessel material"""
    outer_material::StructuralAlloy = StructuralAlloy("Al-2219-T87")
    """Tank head aspect ratio"""
    ARtank::Float64 = 0.0
    """Angular location of inner vessel stiffeners"""
    theta_inner::Float64 = 0.0
    """Vector with angular location of outer vessel stiffeners"""
    theta_outer::Vector{Float64} = Float64[]
    """Number of intermediate stiffeners in outer vessel"""
    Ninterm::Float64 = 1.0
    
    """Venting pressure (Pa)"""
    pvent::Float64 = 0.0
    """Fill pressure (Pa)"""
    pinitial::Float64 = 0.0
    """Minimum allowable tank pressure (Pa)"""
    pmin::Float64 = 0.0
    """Departure hold time (s)"""
    t_hold_orig::Float64 = 0.0
    """Arrival hold time (s)"""
    t_hold_dest::Float64 = 0.0
    """Sea-level temperature for tank design (K)"""
    TSLtank::Vector{Float64} = []

    """Liquid fuel density (kg/m^3)"""
    rhofuel::Float64 = 0.0
    """Liquid fuel temperature in tank (K)"""
    Tfuel::Float64 = 0.0
    """Gas fuel density (kg/m^3)"""
    rhofuelgas::Float64 = 0.0
    """Fuel specific enthalpy of vaporization (J/kg)"""
    hvap::Float64 = 0.0
    """Percentage tank boiloff rate at start of cruise (%/h)"""
    boiloff_rate::Float64 = 0.0

    """Vessel additional mass fraction"""
    ftankadd::Float64 = 0.0
    """Vessel weld efficiency"""
    ew::Float64 = 0.0
    """Minimum ullage fraction"""
    ullage_frac::Float64 = 0.0
    """Heat leakage factor"""
    qfac::Float64 = 0.0
    """Pressure rise factor"""
    pfac::Float64 = 0.0

    """For addtional fuel tank"""
    ACT_eta_vol::Float64 = 0.88
    ACT_eta_wei::Float64 = 0.93
    ACT_A::Float64 = 0.0
    ACT_l::Float64 = 0.0
    ACT_W::Float64 = 0.0
    ACT_dx::Float64 = 0.0
    ACT_fuse_l_extend::Float64 = 0.0
end

@kwdef mutable struct Options_new
    #fuel options
    """Fuel type"""
    opt_fuel::TASOPT.FuelType.T
    """Fuel option index (non-driving; determined and used by gas calcs)"""
    ifuel::Int
    """Secondary fuel option index (Use only at certain phases of flight) """
    ifuel2nd::Int
    """Indicates presence of centerbox fuel tank, can only be true if has_wing_fuel is true"""
    has_centerbox_fuel::Bool
    """Indicates presence of wing fuel tanks """
    has_wing_fuel::Bool
    """Indicates presence of fuselage fuel tanks (non-driving; set by `fuse_tank` inputs)"""
    has_fuselage_fuel::Bool 
      #TODO: consider making ^ a driving parameter, rather than a reflection of fuse_tank parameters
      #Note: right now fuel can only be stored in the wings or the fuselage, not both
    """Indicates presense of additional center fuel tank (ACT)(Works with wing fuel but not fuselage fuel)"""
    has_ACT_fuel::Bool
    """Indicates the need to lengthen the airframe to compensate for any cargo space taken over by any ACT presence"""
    compensate_ACT::Bool
    
    #engine options
    """Engine location"""
    opt_engine_location::TASOPT.EngineLocation.T
    """Propulsion system architecture, performance and weight models set in ac.Engine"""
    opt_prop_sys_arch::TASOPT.PropSysArch.T
    """Calculate takeoff length and engine performance"""
    calculate_takeoff::Bool
    
    #fuselage/cabin options

    """Indicates if the aircraft has a double-decker fuselage configuration"""
    is_doubledecker::Bool

    #wing/stability options
    """Wing position strategy for longitudinal stability analysis: `Fixed` = static wing position, `FixedCLh` = move wing to achieve `CLh = CLhspec` in cruise, `MinStaticMargin` = move wing to achieve minimum static margin = `SMmin`"""
    opt_move_wing::TASOPT.WingMove.T

    #Trefftz plane options
    """Trefftz plane induced drag analysis configuration (discretization, k_tip, bunch, root_contraction)"""
    trefftz_config::TASOPT.aerodynamics.TrefftzPlaneConfig
end

function fuselage_tank_new(fuse_tank::fuselage_tank)
    return fuselage_tank_new(;
        fueltype = fuse_tank.fueltype,
        tank_count = fuse_tank.tank_count,
        placement = fuse_tank.placement,
        sizes_insulation = fuse_tank.sizes_insulation,
        Wfuelintank = fuse_tank.Wfuelintank,
        clearance_fuse = fuse_tank.clearance_fuse,
        t_insul = copy(fuse_tank.t_insul),
        material_insul = copy(fuse_tank.material_insul),
        iinsuldes = copy(fuse_tank.iinsuldes),
        l_cyl_inner = fuse_tank.l_cyl_inner,
        l_inner = fuse_tank.l_inner,
        Rinnertank = fuse_tank.Rinnertank,
        Shead_insul = copy(fuse_tank.Shead_insul),
        inner_material = fuse_tank.inner_material,
        outer_material = fuse_tank.outer_material,
        ARtank = fuse_tank.ARtank,
        theta_inner = fuse_tank.theta_inner,
        theta_outer = copy(fuse_tank.theta_outer),
        Ninterm = fuse_tank.Ninterm,
        pvent = fuse_tank.pvent,
        pinitial = fuse_tank.pinitial,
        pmin = fuse_tank.pmin,
        t_hold_orig = fuse_tank.t_hold_orig,
        t_hold_dest = fuse_tank.t_hold_dest,
        TSLtank = copy(fuse_tank.TSLtank),
        rhofuel = fuse_tank.rhofuel,
        Tfuel = fuse_tank.Tfuel,
        rhofuelgas = fuse_tank.rhofuelgas,
        hvap = fuse_tank.hvap,
        boiloff_rate = fuse_tank.boiloff_rate,
        ftankadd = fuse_tank.ftankadd,
        ew = fuse_tank.ew,
        ullage_frac = fuse_tank.ullage_frac,
        qfac = fuse_tank.qfac,
        pfac = fuse_tank.pfac
    )
end

getf(x, s::Symbol, default) = hasproperty(x, s) ? getproperty(x, s) : default
function Options_new(options)
    return Options_new(;
        opt_fuel = options.opt_fuel,
        ifuel = options.ifuel,
        ifuel2nd = options.ifuel2nd,
        has_centerbox_fuel = options.has_centerbox_fuel,
        has_wing_fuel = options.has_wing_fuel,
        has_fuselage_fuel = options.has_fuselage_fuel,
        has_ACT_fuel = false,
        compensate_ACT = false,
        opt_engine_location = options.opt_engine_location,
        opt_prop_sys_arch = options.opt_prop_sys_arch,
        calculate_takeoff = options.calculate_takeoff,
        is_doubledecker = options.is_doubledecker,
        opt_move_wing = options.opt_move_wing,
        trefftz_config = options.trefftz_config
    )
end

@kwdef mutable struct aircraft_tmp{WS}
    name::String
    description::String
    options::Union{TASOPT.Options,Options_new}

    parg::Vector{Float64}
    parm::Array{Float64,2}
    para::Array{Float64,3}
    pare::Array{Float64,3}

    is_sized::Vector{Bool}

    fuselage::TASOPT.Fuselage
    fuse_tank::Union{TASOPT.fuselage_tank, fuselage_tank_new}

    wing::TASOPT.Wing
    htail::TASOPT.Tail
    vtail::TASOPT.Tail
    engine::TASOPT.Engine
    landing_gear::TASOPT.LandingGear

    wake_system::WS
end

function aircraft_tmp(ac::aircraft)
    return aircraft_tmp(
        ac.name,
        ac.description,
        Options_new(ac.options),
        ac.parg,
        ac.parm,
        ac.para,
        ac.pare,
        ac.is_sized,
        ac.fuselage,
        fuselage_tank_new(ac.fuse_tank),  # <-- conversion happens here
        ac.wing,
        ac.htail,
        ac.vtail,
        ac.engine,
        ac.landing_gear,
        ac.wake_system
    )
end

function update_model!(ac::aircraft_tmp)
    parg = ac.parg #Vector{Float64}
    parm = ac.parm #Array{Float64, 2}
    pare = ac.pare #Array{Float64, 3}

    push!(parg,parg[igWfmax]/parg[igrhofuel]) #igVfmax
    push!(parg,parg[igWfuel]/parg[igrhofuel]) #igVfuel

    parm = vcat(parm, fill(parg[igWfuel]/parg[igrhofuel],size(parm,2))') #imVfuel
    pare = cat(pare, fill(parg[igrhofuel], 1, size(pare,2), size(pare,3)); dims=1)#ierhofuel_driven

    ac.parg = parg
    ac.parm = parm
    ac.pare = pare
end

#
loadpath = joinpath(@__DIR__,"../ModelSaved/acOptimized_BatOptEth/acOptimized_BatOptEth300.jld2")
savepath = joinpath(@__DIR__,"../ModelSaved/testSave.jld2")
@load loadpath ac
ac = aircraft_tmp(ac)
update_model!(ac)
@save savepath ac
