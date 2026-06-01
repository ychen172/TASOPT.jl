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
# Design parameters output
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
    fuel_tank_frac = ac_cur.parm[imWfuel, idx_miss]/ac_cur.parg[igWfmax] #fraction of fuel tank used assuming the same fuel type

    #### fuel volume data
    rhoFuel = ac_cur.parg[igrhofuel] #kg/m3
    volFuel = massFuelTot * 1000.0 / rhoFuel #m3
    
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
        "fuel_tank_frac" => fuel_tank_frac
    )
    return output
end

#### Setup IO
# Input case names - Retrofit
model_dir       = "../ModelSaved"
key_design      = "acOptimized_BatOptJet"
key_offdes      = "acOptimized_BatOptJet_rerun"
des_ranges      = [500,1000,1500,2900,3000]
# Output directory
save_dir      = "../ModelProcessed"
save_name     = "Compare_Design_Offdesign" #sub_folder will be created
# Fields to read out
const fields = (:range_nmi,   :PFEI_JJ, :massTO_Ton, :massFuelTot_Ton, :massPayload_Ton,
                :massEmpty_Ton, :rhoFuel_kgm3,  :volFuel_m3, :LD_cru, :LHV_cru_Jkg, :TSFC_cru_kgsN,
                :vel_cru_ms, :eta_total_cru, :ene_fli_J, :frac_rese, :PFEI_cru_JJ, :range_cru_m,
                :eta_therm_cru, :spe_power_cru_Jkg, :eta_propu_cru, :OPR_cru, :Tt_turbin_cru_K, :fuel_tank_frac)

#### Create save directory
save_dir_sub  = joinpath(save_dir,save_name)
mkpath(save_dir_sub)

#### Initialization
results_design = init_results_2Layers(length(des_ranges), fields)
results_offdes = init_results_2Layers(length(des_ranges), fields)
  
#### Extract data for the design mission
for (i, des_range_cur) in enumerate(des_ranges)
    #### Read in the case
    ac_dir = joinpath(model_dir,key_design,key_design*"$(round(Int,des_range_cur)).jld2")
    ac = quickload_aircraft(ac_dir)
    out_dict = extract_acModel(ac,1)
    println("File, $(ac_dir), read successfully")
    #### Output the data
    for f in fields
        results_design[f][i] = out_dict[String(f)]
    end
    #### Read in the offdesign case
    ac_dir = joinpath(model_dir,key_offdes,key_offdes*"$(round(Int,des_range_cur))_$(round(Int,des_range_cur)).jld2")
    ac = quickload_aircraft(ac_dir)
    out_dict = extract_acModel(ac,2)
    println("File, $(ac_dir), read successfully")
    #### Output the data
    for f in fields
        results_offdes[f][i] = out_dict[String(f)]
    end
end

#### Plotting
# Create style
linestyles = repeat([:solid, :dash, :dot, :dashdot, :dashdotdot],1000)
linecolors = repeat([:blue, :red, :green, :orange, :purple, :brown, :pink, :gray, :black, :cyan,
                     :magenta, :teal, :navy, :maroon, :olive, :gold, :coral, :turquoise, :lime, :indigo], 1000) 
markers = repeat([:rect, :circle, :diamond, :utriangle, :dtriangle],1000)
# Plot PFEI
p1_1 = plot(xlabel="Ranges (nmi)", ylabel="PFEI (J/J)", dpi=800)
global il = 1
scatter!(p1_1, results_design[:range_nmi], results_design[:PFEI_JJ], marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="Design case")
global il += 1
scatter!(p1_1, results_offdes[:range_nmi], results_offdes[:PFEI_JJ], marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="Offdesign case")
savefig(p1_1, joinpath(save_dir_sub, "PFEI.png"))

p1_2 = plot(xlabel="Ranges (nmi)", ylabel="Takeoff Weight (Ton)", dpi=800)
global il = 1
scatter!(p1_2, results_design[:range_nmi], results_design[:massTO_Ton], marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="Design case")
global il += 1
scatter!(p1_2, results_offdes[:range_nmi], results_offdes[:massTO_Ton], marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="Offdesign case")
savefig(p1_2, joinpath(save_dir_sub, "massTO.png"))

p1_3 = plot(xlabel="Ranges (nmi)", ylabel="Fuel Mass (Ton)", dpi=800)
global il = 1
scatter!(p1_3, results_design[:range_nmi], results_design[:massFuelTot_Ton], marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="Design case")
global il += 1
scatter!(p1_3, results_offdes[:range_nmi], results_offdes[:massFuelTot_Ton], marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="Offdesign case")
savefig(p1_3, joinpath(save_dir_sub, "massFuel.png"))

p1_4 = plot(xlabel="Ranges (nmi)", ylabel="Payload Mass (Ton)", dpi=800)
global il = 1
scatter!(p1_4, results_design[:range_nmi], results_design[:massPayload_Ton], marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="Design case")
global il += 1
scatter!(p1_4, results_offdes[:range_nmi], results_offdes[:massPayload_Ton], marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="Offdesign case")
savefig(p1_4, joinpath(save_dir_sub, "massPayload.png"))

p1_5 = plot(xlabel="Ranges (nmi)", ylabel="Empty Mass (Ton)", dpi=800)
global il = 1
scatter!(p1_5, results_design[:range_nmi], results_design[:massEmpty_Ton], marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="Design case")
global il += 1
scatter!(p1_5, results_offdes[:range_nmi], results_offdes[:massEmpty_Ton], marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="Offdesign case")
savefig(p1_5, joinpath(save_dir_sub, "massEmpty.png"))

p1_6 = plot(xlabel="Ranges (nmi)", ylabel="Fractional Fuel Tank Used", dpi=800)
global il = 1
scatter!(p1_6, results_design[:range_nmi], results_design[:fuel_tank_frac], marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="Design case")
global il += 1
scatter!(p1_6, results_offdes[:range_nmi], results_offdes[:fuel_tank_frac], marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="Offdesign case")
savefig(p1_6, joinpath(save_dir_sub, "fraction_fuel_tank.png"))