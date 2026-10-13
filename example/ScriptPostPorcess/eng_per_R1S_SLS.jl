"""
This script compares and plots the SLS engine performance between ethanol and jet fuel aircraft from R1 sizing process
    - Pick one common flight distance/sizing mission and compare the sea-level-static(SLS) engine performance with 0 offtake
    - Changing fuel type allows running a pre-sized engine with different fuel
"""

using TASOPT
include(__TASOPTindices__)
using Plots
include(joinpath(@__DIR__,"../UtilityFunctions/Extract.jl"))
using .Extract: plot_cases_specified, extract_engine_performance
include(joinpath(@__DIR__,"../UtilityFunctions/run_engine.jl"))
using .RunEngine: runOffDes

#### Default parameters for SLS
M0 = 0.0
P0 = 101320.0 #Pa
T0 = 288.2 #K
a0 = 340.2074661144284 #m/s
###############################

#### Input case names - Retrofit
model_dir  = joinpath(__TASOPTroot__,"../example/ModelSaved")
save_dir   = joinpath(__TASOPTroot__,"../example/OutputAndPlots")
case_keys  = ["Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_","Opti_Eth_NoACT_V4_3_R1Sz_EtasEng_"]
case_Nams  = ["Jet Fuel"                         ,"Ethanol"                          ]
fuel_index = [24                                 ,32                                 ]
rho_fuel   = [817.0                              ,789.0                              ]
hvap_fuel  = [358694.0                           ,918187.9                           ]
R1_range   = 2000 #Selected range for SLS comparison (exist in all cases)
n_Thr_ref  = 50 #Number of refinement between max and min thrust to sweep
save_prefi = "R1_Size_SLS_Engine_Compare_$(R1_range)nmi"

#### Create save directory if not exist
mkpath(save_dir)

#### Pre-allocate per-case result containers (one vector per case, filled over the thrust sweep)
thrust_per_case      = [Float64[] for _ in case_keys]
eta_thermal_per_case = [Float64[] for _ in case_keys]
OPR_per_case          = [Float64[] for _ in case_keys]

#### Loop through each case to extract engine performance for SLS
for (j, key_j) in enumerate(case_keys)
    #### Idenfity requested range is available in this case
    ac_path = joinpath(model_dir, key_j, key_j*"$(R1_range).jld2")
    isfile(ac_path) || error("Range $(R1_range) nmi not found for case $(key_j)")
    ac = quickload_aircraft(ac_path)
    println("$(key_j)*$(R1_range) loaded")

    #### Overwrite the fuel used for this engine's SLS sweep. fuel_index/hvap_fuel are indexed in
    #### parallel to case_keys (not derived from it), so e.g. setting fuel_index[1]=32 would run
    #### case_keys[1]'s (Jet's) own sized engine hardware on ethanol instead of its native fuel.
    ac.options.ifuel = fuel_index[j]
    ac.pare[iehvapcombustor,ipcruise1,1] = hvap_fuel[j]

    #### Extract maximum and minimum(non-zero) thrust from the R1 sizing mission
    Fe_all = ac.pare[ieFe, :, 1]
    Fe_max = maximum(Fe_all)
    Fe_min = minimum(Fe_all[Fe_all .> 0.0])

    #### Create thrust levels to test between min and max thrust (some might fail to converge at the lower thrust end)
    Fe_sweep = range(Fe_min, Fe_max, length=n_Thr_ref)

    #### Sweep through the SLS thrust point using runOffDes (no offtake, overwrite fuel type, skip failed to converged point)
    for Fe_target in Fe_sweep
        res = runOffDes(ac, M0, P0, T0, a0, Fe_target; zero_offtake=true)
        if !res.Lconv
            @warn "$(key_j): Fe=$(round(Fe_target/1000.0,digits=1))kN did not converge, skipping"
            continue
        end

        #### Post-process and collect engine efficiencies and cycle diagnostics (skip propulsive and total efficiency due to 0 velocity)
        eta_thermal, _, _ = extract_engine_performance(
            mdot_core     = res.mcore,
            mdotf_mdotc   = res.ff,
            mdot_offtak   = 0.0,
            vel_core_exh  = res.u6,
            vel_fan_exh   = res.u8,
            vel_inlet     = res.u0,
            BPR           = res.BPR,
            area_core_exh = res.A6,
            area_fan_exh  = res.A8,
            P_core_exh    = res.p6,
            P_fan_exh     = res.p8,
            P_inlet       = P0,
            LHV_fuel      = res.hfuel,
            thrust        = res.Fe)
        push!(thrust_per_case[j], res.Fe/1000.0) #[kN]
        push!(eta_thermal_per_case[j], eta_thermal)
        push!(OPR_per_case[j], res.OPR)
    end
end

#### Plots
# Engine thermal efficiencies as a function of thrust
eta_thermal_pct_per_case = [100.0 .* v for v in eta_thermal_per_case]
plot_cases_specified("Thrust [kN]", "Engine Thermal Efficiency at SLS [%]",
                     thrust_per_case, eta_thermal_pct_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_ThermalEff_vs_Thrust.png"))

# Engine thermal efficiencies as a function of pressure ratio
plot_cases_specified("Overall Pressure Ratio [-]", "Engine Thermal Efficiency at SLS [%]",
                     OPR_per_case, eta_thermal_pct_per_case, case_Nams,
                     joinpath(save_dir,"$(save_prefi)_ThermalEff_vs_OPR.png"))
