"""
This script compares and plots the engine performance between ethanol and jet fuel aircraft from R1 sizing process
    - Compare cruise phase engine efficiencies from the R1 sizing mission
    - Compare rotate-phase engine efficiencies
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
save_prefi = "R1_Size_Fli_Eng_Perf"

#### Create save directory if not exist
mkpath(save_dir)

#### Helper: extract engine efficiencies at a given flight phase index
function eta_at_phase(ac, ip::Int)
    return extract_engine_performance(
        mdot_core     = ac.pare[iemcore, ip, 1],
        mdotf_mdotc   = ac.pare[ieff, ip, 1],
        mdot_offtak   = ac.pare[iemofft, ip, 1],
        vel_core_exh  = ac.pare[ieu6, ip, 1],
        vel_fan_exh   = ac.pare[ieu8, ip, 1],
        vel_inlet     = ac.pare[ieu0, ip, 1],
        BPR           = ac.pare[ieBPR, ip, 1],
        area_core_exh = ac.pare[ieA6, ip, 1],
        area_fan_exh  = ac.pare[ieA8, ip, 1],
        P_core_exh    = ac.pare[iep6, ip, 1],
        P_fan_exh     = ac.pare[iep8, ip, 1],
        P_inlet       = ac.pare[iep0, ip, 1],
        LHV_fuel      = ac.pare[iehfuel, ip, 1],
        thrust        = ac.pare[ieFe, ip, 1])
end

#### Helper: extract cycle diagnostics (overall pressure ratio, Tt4/Tt2 temperature ratio,
#### total turbine cooling air fraction, bypass ratio, combustion-gas specific heat) at a given flight phase index
function cycle_diag_at_phase(ac, ip::Int)
    OPR     = ac.pare[iepid, ip, 1]*ac.pare[iepif, ip, 1]*ac.pare[iepilc, ip, 1]*ac.pare[iepihc, ip, 1]
    Tt4_Tt2 = ac.pare[ieTt4, ip, 1]/ac.pare[ieTt2, ip, 1]
    cool_frac = ac.pare[ieepsc1, ip, 1] + ac.pare[ieepsc2, ip, 1] + ac.pare[ieepsc3, ip, 1] + ac.pare[ieepsc4, ip, 1]
    BPR = ac.pare[ieBPR, ip, 1]
    cpt4 = ac.pare[iecpt4, ip, 1]
    return OPR, Tt4_Tt2, cool_frac, BPR, cpt4
end

#### Pre-allocate per-case result containers (one vector per case, filled over available ranges)
rans_per_case           = [Int[] for _ in case_keys]
eta_thermal_per_case    = [Float64[] for _ in case_keys]
eta_propu_per_case      = [Float64[] for _ in case_keys]
eta_total_per_case      = [Float64[] for _ in case_keys]
eta_thermal_TO_per_case = [Float64[] for _ in case_keys]
eta_propu_TO_per_case   = [Float64[] for _ in case_keys]
eta_total_TO_per_case   = [Float64[] for _ in case_keys]
OPR_per_case            = [Float64[] for _ in case_keys]
Tt4Tt2_per_case         = [Float64[] for _ in case_keys]
coolfrac_per_case       = [Float64[] for _ in case_keys]
BPR_per_case            = [Float64[] for _ in case_keys]
OPR_TO_per_case         = [Float64[] for _ in case_keys]
Tt4Tt2_TO_per_case      = [Float64[] for _ in case_keys]
coolfrac_TO_per_case    = [Float64[] for _ in case_keys]
BPR_TO_per_case         = [Float64[] for _ in case_keys]
specpow_per_case        = [Float64[] for _ in case_keys] #Specific power: P_Jet per unit total mass flow, [kW/(kg/s)]
specpow_TO_per_case     = [Float64[] for _ in case_keys]
cpt4_per_case           = [Float64[] for _ in case_keys] #Combustion-gas specific heat at Tt4, [J/(kg*K)]
cpt4_TO_per_case        = [Float64[] for _ in case_keys]

#### Loop through each case to extract engine performance over cruise and takeoff
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
        eta_thermal, eta_propulsion, eta_total, specpow = eta_at_phase(ac, ipcruise1)
        OPR, Tt4Tt2, coolfrac, BPR, cpt4 = cycle_diag_at_phase(ac, ipcruise1)
        # Extract engine performance at rotate condition
        eta_thermal_TO, eta_propulsion_TO, eta_total_TO, specpow_TO = eta_at_phase(ac, iprotate)
        OPR_TO, Tt4Tt2_TO, coolfrac_TO, BPR_TO, cpt4_TO = cycle_diag_at_phase(ac, iprotate)
        # Process and collect engine efficiencies and cycle diagnostics
        push!(rans_per_case[j], ran_i)
        push!(eta_thermal_per_case[j], eta_thermal)
        push!(eta_propu_per_case[j], eta_propulsion)
        push!(eta_total_per_case[j], eta_total)
        push!(eta_thermal_TO_per_case[j], eta_thermal_TO)
        push!(eta_propu_TO_per_case[j], eta_propulsion_TO)
        push!(eta_total_TO_per_case[j], eta_total_TO)
        push!(OPR_per_case[j], OPR)
        push!(Tt4Tt2_per_case[j], Tt4Tt2)
        push!(coolfrac_per_case[j], coolfrac)
        push!(BPR_per_case[j], BPR)
        push!(OPR_TO_per_case[j], OPR_TO)
        push!(Tt4Tt2_TO_per_case[j], Tt4Tt2_TO)
        push!(coolfrac_TO_per_case[j], coolfrac_TO)
        push!(BPR_TO_per_case[j], BPR_TO)
        push!(specpow_per_case[j], specpow/1000.0) #[kW/(kg/s)]
        push!(specpow_TO_per_case[j], specpow_TO/1000.0)
        push!(cpt4_per_case[j], cpt4)
        push!(cpt4_TO_per_case[j], cpt4_TO)
    end

end

#### Plots
# Engine thermal, propulsive, and total effciency at cruise 1, each as its own plot, in percent
eta_thermal_pct_per_case = [100.0 .* v for v in eta_thermal_per_case]
eta_propu_pct_per_case   = [100.0 .* v for v in eta_propu_per_case]
eta_total_pct_per_case   = [100.0 .* v for v in eta_total_per_case]
plot_cases_specified("Design/R1 Range [nmi]", "Engine Thermal Efficiency at Cruise [%]",
                     rans_per_case, eta_thermal_pct_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_ThermalEff.png"))
plot_cases_specified("Design/R1 Range [nmi]", "Engine Propulsive Efficiency at Cruise [%]",
                     rans_per_case, eta_propu_pct_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_PropulsiveEff.png"))
plot_cases_specified("Design/R1 Range [nmi]", "Engine Total Efficiency at Cruise [%]",
                     rans_per_case, eta_total_pct_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_TotalEff.png"))

# Engine thermal, propulsive, and total efficiency at rotate, each as its own plot, in percent
eta_thermal_TO_pct_per_case = [100.0 .* v for v in eta_thermal_TO_per_case]
eta_propu_TO_pct_per_case   = [100.0 .* v for v in eta_propu_TO_per_case]
eta_total_TO_pct_per_case   = [100.0 .* v for v in eta_total_TO_per_case]
plot_cases_specified("Design/R1 Range [nmi]", "Engine Thermal Efficiency at Rotate [%]",
                     rans_per_case, eta_thermal_TO_pct_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_ThermalEff_Rotate.png"))
plot_cases_specified("Design/R1 Range [nmi]", "Engine Propulsive Efficiency at Rotate [%]",
                     rans_per_case, eta_propu_TO_pct_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_PropulsiveEff_Rotate.png"))
plot_cases_specified("Design/R1 Range [nmi]", "Engine Total Efficiency at Rotate [%]",
                     rans_per_case, eta_total_TO_pct_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_TotalEff_Rotate.png"))

# Overall pressure ratio, Tt4/Tt2 temperature ratio, and turbine cooling air fraction, at cruise and rotate
plot_cases_specified("Design/R1 Range [nmi]", "Overall Pressure Ratio at Cruise",
                     rans_per_case, OPR_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_OPR.png"))
plot_cases_specified("Design/R1 Range [nmi]", "Tt4/Tt2 at Cruise",
                     rans_per_case, Tt4Tt2_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_Tt4Tt2.png"))
plot_cases_specified("Design/R1 Range [nmi]", "Turbine Cooling Air Fraction at Cruise",
                     rans_per_case, coolfrac_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_CoolAirFrac.png"))
plot_cases_specified("Design/R1 Range [nmi]", "Overall Pressure Ratio at Rotate",
                     rans_per_case, OPR_TO_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_OPR_Rotate.png"))
plot_cases_specified("Design/R1 Range [nmi]", "Tt4/Tt2 at Rotate",
                     rans_per_case, Tt4Tt2_TO_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_Tt4Tt2_Rotate.png"))
plot_cases_specified("Design/R1 Range [nmi]", "Turbine Cooling Air Fraction at Rotate",
                     rans_per_case, coolfrac_TO_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_CoolAirFrac_Rotate.png"))

# Bypass ratio, at cruise and rotate
plot_cases_specified("Design/R1 Range [nmi]", "Bypass Ratio at Cruise",
                     rans_per_case, BPR_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_BPR.png"))
plot_cases_specified("Design/R1 Range [nmi]", "Bypass Ratio at Rotate",
                     rans_per_case, BPR_TO_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_BPR_Rotate.png"))

# Specific power, P_Jet/(mcore(1+BPR)), as a function of range, at cruise and rotate
plot_cases_specified("Design/R1 Range [nmi]", "Specific Power at Cruise [kW/(kg/s)]",
                     rans_per_case, specpow_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_SpecificPower.png"))
plot_cases_specified("Design/R1 Range [nmi]", "Specific Power at Rotate [kW/(kg/s)]",
                     rans_per_case, specpow_TO_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_SpecificPower_Rotate.png"))

# Combustion-gas specific heat (cpt4) against specific power, at cruise and rotate
plot_cases_specified("Specific Power at Cruise [kW/(kg/s)]", "cpt4 at Cruise [J/(kg*K)]",
                     specpow_per_case, cpt4_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_cpt4_vs_SpecificPower.png"))
plot_cases_specified("Specific Power at Rotate [kW/(kg/s)]", "cpt4 at Rotate [J/(kg*K)]",
                     specpow_TO_per_case, cpt4_TO_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_cpt4_vs_SpecificPower_Rotate.png"))