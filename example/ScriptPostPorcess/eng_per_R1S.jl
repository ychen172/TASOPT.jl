"""
This script compares and plots the engine performance between ethanol and jet fuel aircraft from R1 sizing process
    - Compare cruise phase engine efficiencies from the R1 sizing mission
"""

using TASOPT
include(__TASOPTindices__)
using Plots
using Glob
include(joinpath(@__DIR__,"../UtilityFunctions/Extract.jl"))
using .Extract: extract_acModel_compact!, init_results_2Layers, plot_cases_specified, extract_engine_performance

#### Setup IO
# Input case names - Retrofit
model_dir  = joinpath(__TASOPTroot__,"../example/ModelSaved")
save_dir   = joinpath(__TASOPTroot__,"../example/OutputAndPlots")
case_keys  = ["Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_","Opti_Eth_NoACT_V4_3_R1Sz_EtasEng_"]
case_Nams  = ["Jet Fuel"                         ,"Ethanol"                          ]
save_prefi = "R1_Size_Flight_Engine_Compare"

#### Create save directory if not exist
mkpath(save_dir)

#### Pre-allocate per-case result containers (one vector per case, filled over available ranges)
rans_per_case        = [Int[] for _ in case_keys]
eta_thermal_per_case = [Float64[] for _ in case_keys]
eta_propu_per_case   = [Float64[] for _ in case_keys]
eta_total_per_case   = [Float64[] for _ in case_keys]

#### Loop through each case to extract engine performance over cruise
for (j, key_j) in enumerate(case_keys)
    # Idenfity available ranges within the case
    case_dir = joinpath(model_dir,key_j)
    pattern = Regex("^" * key_j * "(\\d+)\\.jld2\$") #Match "<key><range>.jld2" only, exclude "..._optimization_history.jld2"
    rans_avail = Int[]
    for f in glob(key_j*"*.jld2", case_dir)
        m = match(pattern, basename(f))
        if m !== nothing
            push!(rans_avail, parse(Int, m.captures[1]))
        end
    end
    sort!(rans_avail)
    println("$(key_j): found $(length(rans_avail)) available design ranges")

    #### Loop through each range
    for (i, ran_i) in enumerate(rans_avail)
        # Read in aircraft model
        ac = quickload_aircraft(joinpath(model_dir,key_j,key_j*"$(round(Int,ran_i)).jld2"))
        println("$(key_j)*$(ran_i) loaded")
        # Extract engine performance at cruise 1 condition
        eta_thermal, eta_propulsion, eta_total = extract_engine_performance(
            mdot_core     = ac.pare[iemcore, ipcruise1, 1],
            mdotf_mdotc   = ac.pare[ieff, ipcruise1, 1],
            mdot_offtak   = ac.pare[iemofft, ipcruise1, 1],
            vel_core_exh  = ac.pare[ieu6, ipcruise1, 1],
            vel_fan_exh   = ac.pare[ieu8, ipcruise1, 1],
            vel_inlet     = ac.pare[ieu0, ipcruise1, 1],
            BPR           = ac.pare[ieBPR, ipcruise1, 1],
            area_core_exh = ac.pare[ieA6, ipcruise1, 1],
            area_fan_exh  = ac.pare[ieA8, ipcruise1, 1],
            P_core_exh    = ac.pare[iep6, ipcruise1, 1],
            P_fan_exh     = ac.pare[iep8, ipcruise1, 1],
            P_inlet       = ac.pare[iep0, ipcruise1, 1],
            LHV_fuel      = ac.pare[iehfuel, ipcruise1, 1],
            thrust        = ac.pare[ieFe, ipcruise1, 1])
        # Process and collect engine efficiencies at cruise 1 condition
        push!(rans_per_case[j], ran_i)
        push!(eta_thermal_per_case[j], eta_thermal)
        push!(eta_propu_per_case[j], eta_propulsion)
        push!(eta_total_per_case[j], eta_total)
    end

end

#### Plots
# Engine thermal, propulsive, and total effciency at cruise 1, each as its own plot, in percent
eta_thermal_pct_per_case = [100.0 .* v for v in eta_thermal_per_case]
eta_propu_pct_per_case   = [100.0 .* v for v in eta_propu_per_case]
eta_total_pct_per_case   = [100.0 .* v for v in eta_total_per_case]
plot_cases_specified("Design/R1 Range [nmi]", "Engine Thermal Efficiency at Cruise 1 [%]",
                     rans_per_case, eta_thermal_pct_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_ThermalEff.png"))
plot_cases_specified("Design/R1 Range [nmi]", "Engine Propulsive Efficiency at Cruise 1 [%]",
                     rans_per_case, eta_propu_pct_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_PropulsiveEff.png"))
plot_cases_specified("Design/R1 Range [nmi]", "Engine Total Efficiency at Cruise 1 [%]",
                     rans_per_case, eta_total_pct_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_TotalEff.png"))