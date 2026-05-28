"""
This script compare the off-design performance between cases in terms of PFEI
and flight range
"""

using TASOPT
using DataFrames, CSV
include(__TASOPTindices__)
using Plots
include(joinpath(@__DIR__,"../Breguet_range_solve_offdes.jl"))
using .Breguet: Bre_off_des
using Statistics

# Other constants
epsR1R2 = 0.005

#### Helpers Functions
# Interpolation function for threshold range determination
function interp_range_at_threshold(range, y, i_peak, eps; side::Symbol)
    yth = y[i_peak] * (1 - eps)

    if side == :right
        # first point to the RIGHT of peak below threshold
        k2 = findnext(v -> v < yth, y, i_peak + 1)
        isnothing(k2) && return nothing
        k1 = k2 - 1
    elseif side == :left
        # first point to the LEFT of peak below threshold
        k1 = findprev(v -> v < yth, y, i_peak - 1)
        isnothing(k1) && return nothing
        k2 = k1 + 1
    else
        error("side must be :right or :left")
    end

    x1, x2 = range[k1], range[k2]
    y1, y2 = y[k1], y[k2]

    # linear interpolation: x at y=yth
    if y2 == y1
        return x1  # degenerate case
    end
    xth = x1 + (yth - y1) * (x2 - x1) / (y2 - y1)
    return xth
end

# Design parameters output
function init_results(case_keywords, des_ranges, offdes_ranges, fields_to_read)
    # case_keywords, des_ranges, offdes_ranges: Vector
    # fields_to_read: tuple of symbols
    noff = length(offdes_ranges)
    makevec() = fill!(Vector{Union{Missing,Float64}}(undef, noff), missing)
    Dict(
        kw => Dict(
            d => Dict(
                f => makevec() for f in fields_to_read
            ) for d in eachindex(des_ranges)
        ) for kw in eachindex(case_keywords)
    )
end

function init_results_3Layers(num_des_ranges, fields, num_offdes_ranges)
    """This function create a three layer dictionary"""
    # num_des_ranges: Int: number of design ranges
    # fields: tuple of symbols
    # num_offdes_ranges: Int: number of offdesign ranges
    makevec() = fill!(Vector{Union{Missing,Float64}}(undef, num_offdes_ranges), missing)
    output = Dict(d => Dict(
                    f => makevec() for f in fields
                  ) for d = 1:num_des_ranges
    )
    return output
end

function init_results_2Layers(num_des_ranges, fields)
    """This function create a two layer dictionary"""
    # num_des_ranges: Int: number of design ranges
    # fields: tuple of symbols
    makevec() = fill!(Vector{Union{Missing,Float64}}(undef, num_des_ranges), missing)
    output = Dict(f => makevec() for f in fields)
    return output
end

function find_frac_change(results,idx_base,idx_targ,numDesRan,xSymb,ySymb)
    #results: case*desRan*symb*offdes
    #idx_base,idx_targ: cases
    #numDesRan: number of design range
    #xSymb,ySymb: symbols to find difference, Ex. x: range, y: PFEI
    x_common = []
    y_change = []
    for i in 1:numDesRan
        #### Extract the base and target case
        results_base = results[idx_base][i]
        results_targ = results[idx_targ][i]
        x_base = results_base[xSymb]
        x_targ = results_targ[xSymb]
        y_base = results_base[ySymb]
        y_targ = results_targ[ySymb]
        #### Filter out the common cases using x
        minx = max(minimum(x_base),minimum(x_targ))
        maxx = min(maximum(x_base),maximum(x_targ))
        msk_base = (x_base .>= minx) .& (x_base .<= maxx)
        msk_targ = (x_targ .>= minx) .& (x_targ .<= maxx)
        x_base = x_base[msk_base]
        x_targ = x_targ[msk_targ]
        y_base = y_base[msk_base]
        y_targ = y_targ[msk_targ]
        @assert (round.(Int, x_base) == round.(Int, x_targ)) "Base$(x_base)/target$(x_targ) x do not match after filtering"
        #### Calculate the change
        push!(x_common, x_base)
        push!(y_change, (y_targ .- y_base) ./ y_base)
    end
    return x_common,y_change
end

function extract_acModel(ac_cur,idx_miss)
    """
    This function extract parameters from a aircraft model given a target mission
    Inputs:
        ac_cur: aircraft model
        idx_miss: index of the mission to extract data
    """
    #### overall performance data
    range_cur = ac_cur.parm[imRange, idx_miss] / 1852.0 #(nmi)
    PFEI_cur = ac_cur.parm[imPFEI, idx_miss] #(J/J)
    
    #### mass data
    massTO = ac_cur.parm[imWTO, idx_miss]/gee/1000.0 #Takeoff mass (Ton)
    massFuelTot = ac_cur.parm[imWfuel, idx_miss]/gee/1000.0 #Fuel mass (Ton)(Include reserved and burned)
    massPayload = ac_cur.parm[imWpay, idx_miss]/gee/1000.0 #(Ton)
    massEmpty = massTO - massFuelTot - massPayload #(Ton) empty weight
    
    #### fuel volume data
    rhoFuel = ac_cur.parg[igrhofuel] #kg/m3
    volFuel = massFuelTot * 1000.0 / rhoFuel #m3
    volFuelMax = ac_cur.parg[igWfmax] / gee / rhoFuel #Only work if design and off-design use the same fuel. Fuel volume maximum m3
    FuelVolumeFraction = volFuel/volFuelMax
    
    #### flight performance data
    LD_cruise = 0.5 * (ac_cur.para[iaCL, ipcruise1, idx_miss]/ac_cur.para[iaCD, ipcruise1, idx_miss] + 
                        ac_cur.para[iaCL, ipcruise2, idx_miss]/ac_cur.para[iaCD, ipcruise2, idx_miss]) #Averaged cruise lift-to-drag ratio
    LHV_cruise = 0.5 * (ac_cur.pare[iehfuel, ipcruise1, idx_miss] + ac_cur.pare[iehfuel, ipcruise2, idx_miss]) #Averaged cruise heating value (J/kg) (Include vaporization already)
    TSFC_cruise = 0.5 * (ac_cur.pare[ieTSFC, ipcruise1, idx_miss] + ac_cur.pare[ieTSFC, ipcruise2, idx_miss]) / gee #Averaged cruise thrust specfic heat consumption (kg/s/N)
    vel_cruise = 0.5 * (cos(ac_cur.para[iagamV, ipcruise1, idx_miss]) * ac_cur.pare[ieu0, ipcruise1, idx_miss] + 
                        cos(ac_cur.para[iagamV, ipcruise2, idx_miss]) * ac_cur.pare[ieu0, ipcruise2, idx_miss]) #Averaged cruise horizontal velocity (m/s)
    eta_total_cruise = (1.0/TSFC_cruise)*(vel_cruise/LHV_cruise) #total cruise engine efficiency

    #### derived data
    energy_flight = PFEI_cur * massPayload * (1000.0 * gee * 1852.0) * range_cur #(J) total flight energy 

    #### parameters for Breguet range
    frac_rese = ac_cur.parg[igfreserve] #W_reserveFuel / W_fuelburned
    PFEI_cru = (LHV_cruise*(ac_cur.para[iafracW, ipcruise1, idx_miss]-ac_cur.para[iafracW, ipcruise2, idx_miss])*ac_cur.parg[igWMTO]/gee)/
                (ac_cur.parm[imWpay, idx_miss] * ac_cur.parm[imRange, idx_miss]) #Cruise Only PFEI (J/J)
    range_cru = ac_cur.para[iaRange, ipcruise2, idx_miss] - ac_cur.para[iaRange, ipcruise1, idx_miss] #(m)
    @assert range_cru>0.0 "Find a negative cruise range, likely range to short"

    #### parameters for engine efficiency at cruise condition
    etaTherm_cru = Float64[]
    spePower_cru = Float64[]
    etaPropu_cru = Float64[]
    OPR_cru = Float64[]
    Tt_TurbIn_cru = Float64[]
    for phase in [ipcruise1,ipcruise2]
        ff_cru = ac_cur.pare[ieff, phase, idx_miss] #mdot_fuel / mdot_core
        BPR_cru = ac_cur.pare[ieBPR, phase, idx_miss] #mdot_BP / mdot_core
        mass_offtake_cru = ac_cur.pare[iemofft, phase, idx_miss] #kg/s single engine
        mass_core_cru = ac_cur.pare[iemcore, phase, idx_miss] #kg/s single engine
        u_coreExh_cru = ac_cur.pare[ieu6, phase, idx_miss] #m/s
        u_fanExh_cru = ac_cur.pare[ieu8, phase, idx_miss] #m/s
        u_inf_cru = ac_cur.pare[ieu0, phase, idx_miss] #m/s
        p_coreExh_cru = ac_cur.pare[iep6, phase, idx_miss] #Pa
        p_fanExh_cru = ac_cur.pare[iep8, phase, idx_miss] #Pa
        p_inf_cru = ac_cur.pare[iep0, phase, idx_miss] #Pa
        A_coreExh_cru = ac_cur.pare[ieA6, phase, idx_miss] #m2
        A_fanExh_cru = ac_cur.pare[ieA8, phase, idx_miss] #m2
        LHV_cru = ac_cur.pare[iehfuel, phase, idx_miss] #J/kg including vaporization heat
        Thrust_cru = ac_cur.pare[ieFsp, phase, idx_miss] * (u_inf_cru * mass_core_cru * (1.0 + BPR_cru)) #N
        Tt41_cru = ac_cur.pare[ieTt41, phase, idx_miss] #(K) turbine inlet temperature after cooling air
        P_Jet_cru = 0.5*(mass_core_cru*(1.0+ff_cru)-mass_offtake_cru)*u_coreExh_cru^2 +
                    0.5*mass_core_cru*BPR_cru*u_fanExh_cru^2 - 
                    0.5*mass_core_cru*(1.0 + BPR_cru)*u_inf_cru^2 + 
                    (p_coreExh_cru-p_inf_cru)*A_coreExh_cru*u_coreExh_cru +
                    (p_fanExh_cru-p_inf_cru)*A_fanExh_cru*u_fanExh_cru #Jet power (J/s)
        push!(OPR_cru, ac_cur.pare[iepid, phase, idx_miss]*ac_cur.pare[iepif, phase, idx_miss]*ac_cur.pare[iepilc, phase, idx_miss]*ac_cur.pare[iepihc, phase, idx_miss])
        push!(etaTherm_cru, P_Jet_cru/(mass_core_cru*ff_cru*LHV_cru)) #Thermal efficiency
        push!(spePower_cru, P_Jet_cru/(mass_core_cru*(1.0+ff_cru+BPR_cru)-mass_offtake_cru)) #J/kg
        push!(etaPropu_cru, (Thrust_cru*u_inf_cru)/P_Jet_cru)
        push!(Tt_TurbIn_cru, Tt41_cru)
    end
    spePower_cru = mean(spePower_cru) #J/kg
    etaTherm_cru = mean(etaTherm_cru)
    etaPropu_cru = mean(etaPropu_cru)            
    OPR_cru = mean(OPR_cru)
    Tt_TurbIn_cru = mean(Tt_TurbIn_cru) #K

    #### Return
    output = Dict(
        "range_nmi" => range_cur,
        "PFEI_JJ" => PFEI_cur,
        "massTO_Ton" => massTO,
        "massFuelTot_Ton" => massFuelTot,
        "massPayload_Ton" => massPayload,
        "massEmpty_Ton" => massEmpty,
        "rhoFuel_kgm3" => rhoFuel,
        "volFuel_m3" => volFuel,
        "LD_cru" => LD_cruise,
        "LHV_cru_Jkg" => LHV_cruise,
        "TSFC_cru_kgsN" => TSFC_cruise,
        "vel_cru_ms" => vel_cruise,
        "eta_total_cru" => eta_total_cruise,
        "ene_fli_J" => energy_flight,
        "frac_rese" => frac_rese,
        "PFEI_cru_JJ" => PFEI_cru,
        "range_cru_m" => range_cru,
        "eta_therm_cru" => etaTherm_cru,
        "spe_power_cru_Jkg" => spePower_cru,
        "eta_propu_cru" => etaPropu_cru,
        "OPR_cru" => OPR_cru,
        "Tt_turbin_cru_K" => Tt_TurbIn_cru,
        "FuelVolumeFraction"=> FuelVolumeFraction
    )
    return output
end

function findR1R2(range, massPay, volFuel, epsR1R2)
    """
    A function to find R1 and R2 range and their index
    Inputs:
        range, massPay, volFuel: Vector{Float64}, Ascending range
        epsR1R2: Float64: eps lower than
    """
    # R1
    massPay_Max_Val, massPay_Max_Idx = findmax(massPay) #assume your off-design mission reach maximum payload (if not, the first point would be picked)
    idx_R1_approx = findnext(massPay_ -> massPay_ < massPay_Max_Val*(1.0-epsR1R2), massPay, massPay_Max_Idx) - 1 #approximated R1 index for plot
    R1_interp = interp_range_at_threshold(range, massPay, massPay_Max_Idx, epsR1R2; side=:right) #inteprolated R1
    # R2
    volFuel_Max_Val, volFuel_Max_Idx = findmax(volFuel) #assume your off-design reach maximum fuel volume (threortically always for payload range diagram)
    idx_R2_approx = findprev(volFuel_ -> volFuel_ < volFuel_Max_Val*(1.0-epsR1R2), volFuel, volFuel_Max_Idx) + 1
    R2_interp = interp_range_at_threshold(range, volFuel, volFuel_Max_Idx, epsR1R2; side=:left) #inteprolated R2
    # Output
    output = Dict(
        "idx_R1_approx" => idx_R1_approx,
        "idx_R2_approx" => idx_R2_approx,
        "R1_interp" => R1_interp, #Same unit as the input range
        "R2_interp" => R2_interp,
    )
    return output
end

#### Setup IO
# Input case names - Retrofit
model_dir       = "../ModelSaved"
# Input case names - Ethanol
key_sized_Eth       = "acOptimized_BatOptEth" #180Pass_Opt: "acOptimized_BatOptEth", 3000nmiJet2EthRetroCase_but_optimized: "OptimizedJetToEth3000_"
des_range_sized_Eth = float.(collect(300:100:2900)) #Has to match with existings
# Input case names - Jet fuel
key_sized_Jet       = "acOptimized_BatOptJet" #180Pass_Opt: "CenterFuelTank_BatOptEth", 3000nmiJet2EthRetroCase_but_optimizedWithCenTank: "OptCenTankJetToEth3000_"
des_range_sized_Jet = float.(collect(300:100:3000)) #Has to match with existings
# Output directory
save_dir      = "../ModelProcessed"
save_name     = "Compare_at_design" #sub_folder will be created
# Fields to read out
const fields = (:range_nmi,   :PFEI_JJ, :massTO_Ton, :massFuelTot_Ton, :massPayload_Ton,
                :massEmpty_Ton, :rhoFuel_kgm3,  :volFuel_m3, :LD_cru, :LHV_cru_Jkg, :TSFC_cru_kgsN,
                :vel_cru_ms, :eta_total_cru, :ene_fli_J, :frac_rese, :PFEI_cru_JJ, :range_cru_m,
                :eta_therm_cru, :spe_power_cru_Jkg, :eta_propu_cru, :OPR_cru, :Tt_turbin_cru_K, :FuelVolumeFraction)

#### Create save directory
save_dir_sub  = joinpath(save_dir,save_name)
mkpath(save_dir_sub)

#### Initialization
results_sized_Eth = init_results_2Layers(length(des_range_sized_Eth), fields)
results_sized_Jet = init_results_2Layers(length(des_range_sized_Jet), fields)

#### Extract data for the sized missions
for (i, des_ranges_cur) in enumerate(des_range_sized_Eth)
    #### Read in the case
    ac_dir_sized = joinpath(model_dir,key_sized_Eth,key_sized_Eth*"$(round(Int,des_ranges_cur)).jld2")
    ac_sized = quickload_aircraft(ac_dir_sized)
    println("File, $(ac_dir_sized), read successfully")
    out_dict = extract_acModel(ac_sized,1)
    
    #### Output the data
    for f in fields
        results_sized_Eth[f][i] = out_dict[String(f)]
    end
end

for (i, des_ranges_cur) in enumerate(des_range_sized_Jet)
    #### Read in the case
    ac_dir_sized = joinpath(model_dir,key_sized_Jet,key_sized_Jet*"$(round(Int,des_ranges_cur)).jld2")
    ac_sized = quickload_aircraft(ac_dir_sized)
    println("File, $(ac_dir_sized), read successfully")
    out_dict = extract_acModel(ac_sized,1)
    
    #### Output the data
    for f in fields
        results_sized_Jet[f][i] = out_dict[String(f)]
    end
end

#### Plotting
# Create style
linestyles = repeat([:solid, :dash, :dot, :dashdot, :dashdotdot],1000)
linecolors = repeat([:blue, :red, :green, :orange, :purple, :brown, :pink, :gray, :black, :cyan,
                     :magenta, :teal, :navy, :maroon, :olive, :gold, :coral, :turquoise, :lime, :indigo], 1000) 
markers = repeat([:rect, :circle, :diamond, :utriangle, :dtriangle],1000)
# Plot PFEI
p1_1 = plot(xlabel="Design Ranges (nmi)", ylabel="PFEI (J/J)", dpi=800)
plot!(p1_1,ylims=(0.6,1.4))
global il = 1
plot!(p1_1, results_sized_Jet[:range_nmi], results_sized_Jet[:PFEI_JJ], marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="Optimized Jet Fuel Aircraft")
global il += 1
plot!(p1_1, results_sized_Eth[:range_nmi], results_sized_Eth[:PFEI_JJ], marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="Optimized Ethanol Aircraft")
savefig(p1_1, joinpath(save_dir_sub, "PFEI.png"))

# Plot Fuel Fraction
p1_2 = plot(xlabel="Design Ranges (nmi)", ylabel="Fuel Tank Used (%)", dpi=800)
global il = 1
plot!(p1_2, results_sized_Jet[:range_nmi], results_sized_Jet[:FuelVolumeFraction] .* 100.0, marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="Optimized Jet Fuel Aircraft")
global il += 1
plot!(p1_2, results_sized_Eth[:range_nmi], results_sized_Eth[:FuelVolumeFraction] .* 100.0, marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="Optimized Ethanol Aircraft")
savefig(p1_2, joinpath(save_dir_sub, "FuelTankUsed.png"))