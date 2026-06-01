"""
This script compare for off-design cases the change of R1 and R2 for two cases.
Here R1 and R2 are assumed to have been computed by bisection searching instead of linear interpolation.
"""

using TASOPT
using DataFrames, CSV
include(__TASOPTindices__)
using Plots
using Glob
using Statistics

#### For Plotting Style
linestyles = repeat([:solid, :dash, :dot, :dashdot, :dashdotdot],1000)
linecolors = repeat([:blue, :red, :green, :orange, :purple, :brown, :pink, :gray, :black, :cyan,
                     :magenta, :teal, :navy, :maroon, :olive, :gold, :coral, :turquoise, :lime, :indigo], 1000) 
markers = repeat([:rect, :circle, :diamond, :utriangle, :dtriangle],1000)

#### Helpers Functions
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

# Fields to read out
const fields = (:range_nmi,   :PFEI_JJ, :massTO_Ton, :massFuelTot_Ton, :massPayload_Ton,
                :massEmpty_Ton, :rhoFuel_kgm3,  :volFuel_m3, :LD_cru, :LHV_cru_Jkg, :TSFC_cru_kgsN,
                :vel_cru_ms, :eta_total_cru, :ene_fli_J, :frac_rese, :PFEI_cru_JJ, :range_cru_m,
                :eta_therm_cru, :spe_power_cru_Jkg, :eta_propu_cru, :OPR_cru, :Tt_turbin_cru_K, :fuel_tank_frac)

#### Setup IO
# Input case names - Retrofit
model_dir  = "../ModelSaved"
key1       = "acOptimized_BatOptJet_OffDes_Eth_" #Change of case 1 with respec to case 2
key2       = "acOptimized_BatOptJet_OffDes_Jet_"
nam1       = "Retrofit"
nam2       = "Baseline"
desRs      = collect(300:100:2900) #Has to be a common range for the comparison
# Output directory
save_dir      = "../ModelProcessed"
save_name     = "R1_R2_comparison" #sub_folder will be created

#### Create save directory
save_dir_sub  = joinpath(save_dir,save_name)
mkpath(save_dir_sub)

#### Initialization
case1_R1_res = init_results_2Layers(length(desRs), fields)
case1_R2_res = init_results_2Layers(length(desRs), fields)
case2_R1_res = init_results_2Layers(length(desRs), fields)
case2_R2_res = init_results_2Layers(length(desRs), fields)
  
#### Extract data for the design mission
for (i, desR) in enumerate(desRs)
    #### Read in the case
    case1_name_local = key1*"$(round(Int,desR))_"
    case2_name_local = key2*"$(round(Int,desR))_"
    case1_dir = joinpath(model_dir,key1,case1_name_local)
    case2_dir = joinpath(model_dir,key2,case2_name_local)
    #Search for R1 and R2 model
    case1R1_file = glob(case1_name_local*"R1_*.jld2", case1_dir)
    case1R2_file = glob(case1_name_local*"R2_*.jld2", case1_dir)
    case2R1_file = glob(case2_name_local*"R1_*.jld2", case2_dir)
    case2R2_file = glob(case2_name_local*"R2_*.jld2", case2_dir)
    @assert (length(case1R1_file) == 1) "Need one match for $(joinpath(case1_dir,case1_name_local*"R1_*.jld2"))"
    @assert (length(case1R2_file) == 1) "Need one match for $(joinpath(case1_dir,case1_name_local*"R2_*.jld2"))"
    @assert (length(case2R1_file) == 1) "Need one match for $(joinpath(case2_dir,case2_name_local*"R1_*.jld2"))"
    @assert (length(case2R2_file) == 1) "Need one match for $(joinpath(case2_dir,case2_name_local*"R2_*.jld2"))"
    
    #Readin the aircraft model
    ac_c1r1 = quickload_aircraft(case1R1_file[1])
    out_c1r1 = extract_acModel(ac_c1r1,2)
    println("File, $(case1R1_file[1]), read successfully")
    #
    ac_c1r2 = quickload_aircraft(case1R2_file[1])
    out_c1r2 = extract_acModel(ac_c1r2,2)
    println("File, $(case1R2_file[1]), read successfully")
    #
    ac_c2r1 = quickload_aircraft(case2R1_file[1])
    out_c2r1 = extract_acModel(ac_c2r1,2)
    println("File, $(case2R1_file[1]), read successfully")
    #
    ac_c2r2 = quickload_aircraft(case2R2_file[1])
    out_c2r2 = extract_acModel(ac_c2r2,2)
    println("File, $(case2R2_file[1]), read successfully")
    
    #### Output the data
    for f in fields
        case1_R1_res[f][i] = out_c1r1[String(f)]
        case1_R2_res[f][i] = out_c1r2[String(f)]
        case2_R1_res[f][i] = out_c2r1[String(f)]
        case2_R2_res[f][i] = out_c2r2[String(f)]
    end
end

#### Compute R1 and R2 change
R1Change = (case1_R1_res[:range_nmi] .- case2_R1_res[:range_nmi]) ./ case2_R1_res[:range_nmi]
R2Change = (case1_R2_res[:range_nmi] .- case2_R2_res[:range_nmi]) ./ case2_R2_res[:range_nmi]

#### Plotting
# Plot R1 R2
p1_1 = plot(xlabel="Design Range (nmi)", ylabel="Off-design Range (nmi)", dpi=800)
global il = 1
plot!(p1_1, desRs, case1_R1_res[:range_nmi], marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(nam1) R1")
global il += 1
plot!(p1_1, desRs, case1_R2_res[:range_nmi], marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(nam1) R2")
global il += 1
plot!(p1_1, desRs, case2_R1_res[:range_nmi], marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(nam2) R1")
global il += 1
plot!(p1_1, desRs, case2_R2_res[:range_nmi], marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(nam2) R2")
savefig(p1_1, joinpath(save_dir_sub, "R1_R2.png"))

# Plot R1 R2 Change
mskCha = desRs .>= 1500
p1_2 = plot(xlabel="Design Range (nmi)", ylabel="Ranges Reduction from Retrofitting (%)", dpi=800)
global il = 1
plot!(p1_2, desRs[mskCha], -R1Change[mskCha] .* 100.0, marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="R1")
global il += 1
plot!(p1_2, desRs[mskCha], -R2Change[mskCha] .* 100.0, marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="R2")
savefig(p1_2, joinpath(save_dir_sub, "R1_R2_Change.png"))