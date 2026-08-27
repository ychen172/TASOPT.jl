using TASOPT
using Plots
using CSV, DataFrames
include(__TASOPTindices__)
include(joinpath(__TASOPTroot__, "../example/utilities_engines/run_engine.jl"))
include(joinpath(__TASOPTroot__, "../example/utilities_for_optimization/objective_factory.jl"))
include(joinpath(__TASOPTroot__, "../example/utilities_engines/calibrate_engine.jl"))

## IO Directory
ac_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_/Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_2400.jld2")
miss_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/reference_engine_cycle_CFM5B/Target_Leap1B.csv")
save_key = "Leap1B"
save_dir = joinpath(__TASOPTroot__,"../example/ModelProcessed",save_key)
mkpath(save_dir)

## Load aircraft model and reference mission
ac = quickload_aircraft(ac_dir)
df = CSV.read(miss_dir, DataFrame)
Fn_N          = Float64.(df[:, "Thrust (kN)"] .* 1000.0) #[N]
BPR_Ref       = df[:, "BPR"]
OPR_Ref       = df[:, "OPR"]
Wfuel_Ref_kgs = Float64.(df[:, "Wf[lbm/s]"] .* 0.453592)

## Setup engine inlet conditions (EEDB reference condition)
M0 = 0.0
P0 = 101320.0 #Pa
T0 = 288.2 #K
a0 = 340.2074661144284 #m/s
Fan_Dia_Ref_m = 1.76 #m

## Candidates to compare, hardcoded directly from known TechRun results (not loaded from
## tech_opt result .jld2 files -- not all saved runs converged, and only the 2 x-vectors are
## needed here, not the rest of what's stored). numStageLC/numStageHC only apply when a
## candidate's x_eng has 4 elements (hard-ratio [BPR,pif,per_stage,Tt4]); ignored for 5-element
## [BPR,pif,pilc,pihc,Tt4] free-split candidates.
numStageLC = 3
numStageHC = 10
cases = [
    (label = "Free Split",
     x_tec = [0.9548557021511898, 0.927847659728918, 0.911476635195622, 0.9138021689093785, 0.9301471870102677, 0.9342863878824652, 0.9841208709069963, 1306.8933881534017],
     x_eng = [7.552368307247971, 1.7136459907314052, 1.2838292230496653, 32.92570978817826, 1585.5595451382987]),
    (label = "Fixed Ratio",
     x_tec = [0.972551729598183, 0.96, 0.944250570036917, 0.9411539029572416, 0.9296208234471574, 0.9551522914904217, 0.9870058969320323, 1303.5520348658742],
     x_eng = [7.483899219774775, 1.791880097483402, 1.3340410905093636, 1519.8172766551688]),
     (label = "Fixed Ratio 95% pib",
     x_tec = [0.9500099999999997, 0.9497451512602837, 0.9336044518105325, 0.9356643828707908, 0.9571704792856035, 0.9410155818926765, 0.9874838649583609, 1322.9524543776065],
     x_eng = [7.271241374728679, 1.80700247830106, 1.3303910509298469, 1518.0197174687125]),
]
case_legend = [c.label for c in cases]
numTry = length(cases)
numFn  = length(Fn_N)

Fn_N_Calc      = fill(NaN, numFn, numTry)
Wfuel_Calc_kgs = fill(NaN, numFn, numTry)
BPR_Calc       = fill(NaN, numFn, numTry)
OPR_Calc       = fill(NaN, numFn, numTry)
DFan_Calc_m    = fill(NaN, numFn, numTry)

for (idx, case) in enumerate(cases)
    ac_local = deepcopy(ac)

    ## Apply the 8 technology parameters (broadcast across every mission point)
    pib,epolf,epollc,epolhc,epolht,epollt,etab,Tmetal = case.x_tec
    ac_local.pare[iepib,   :, :] .= pib
    ac_local.pare[ieepolf, :, :] .= epolf
    ac_local.pare[ieepollc,:, :] .= epollc
    ac_local.pare[ieepolhc,:, :] .= epolhc
    ac_local.pare[ieepolht,:, :] .= epolht
    ac_local.pare[ieepollt,:, :] .= epollt
    ac_local.pare[ieetab,  :, :] .= etab
    ac_local.parg[igTmetal]       = Tmetal

    ## Apply the cycle design variables and resize -- auto-detect hard-ratio (4-elem) vs free-split (5-elem)
    try
        is_hard_ratio = length(case.x_eng) == 4
        if is_hard_ratio
            CaliEng.UpdAcEngMod!(ac_local, case.x_eng; iter_sizing=150,
                                 hard_lpc_hpc_stage_ratio=true, numStageLC=numStageLC, numStageHC=numStageHC)
        else
            CaliEng.UpdAcEngMod!(ac_local, case.x_eng; iter_sizing=150)
        end
    catch e
        @warn "Case \"$(case.label)\" failed to size, skipping" exception=e
        continue
    end

    ## Run the engine through the off-design missions -- zero_offtake=true to match EEDB's
    ## no-customer-offtake measurement convention
    for (idx_Fn,Fn_N_cur) in enumerate(Fn_N)
        try
            res = RunEngine.runOffDes(ac_local, M0, P0, T0, a0, Fn_N_cur; zero_offtake=true)
            Fn_N_Calc[idx_Fn, idx]      = res.Fe        #[N]
            Wfuel_Calc_kgs[idx_Fn, idx] = res.mcore * res.ff #[kg/s]
            BPR_Calc[idx_Fn, idx]       = res.BPR
            OPR_Calc[idx_Fn, idx]       = res.OPR
            DFan_Calc_m[idx_Fn, idx]    = ac_local.parg[igdfan] #[m]
        catch e
            @warn "Case \"$(case.label)\" failed at Fn_N=$(Fn_N_cur), leaving as NaN" exception=e
        end
    end
end

## Compare the engine performance
const MARKERS = [:rect, :circle, :diamond, :utriangle, :dtriangle, :cross, :xcross, :star4, :star5, :star6, :star7, :star8]

# Fuel flow rate
p = plot(xlabel="Thrust [kN]", ylabel="Fuel Flow [kg/s]", dpi=800, legend=:best)
for idx in 1:numTry
    plot!(p, Fn_N_Calc[:,idx] ./ 1000.0,  Wfuel_Calc_kgs[:,idx], marker=MARKERS[idx+1], label=case_legend[idx])
end
plot!(p, Fn_N ./ 1000.0,  Wfuel_Ref_kgs, marker=MARKERS[1], label="Reference")
savefig(p,joinpath(save_dir,"$(save_key)_Fuel_Flow.png"))

# Bypass ratio
p = plot(xlabel="Thrust [kN]", ylabel="Bypass Ratio", dpi=800, legend=:best)
for idx in 1:numTry
    plot!(p, Fn_N_Calc[:,idx] ./ 1000.0,  BPR_Calc[:,idx], marker=MARKERS[idx+1], label=case_legend[idx])
end
plot!(p, Fn_N ./ 1000.0,  BPR_Ref, marker=MARKERS[1], label="Reference")
savefig(p,joinpath(save_dir,"$(save_key)_Bypass_Ratio.png"))

# Overall Pressure Ratio
p = plot(xlabel="Thrust [kN]", ylabel="Overall Pressure Ratio", dpi=800, legend=:best)
for idx in 1:numTry
    plot!(p, Fn_N_Calc[:,idx] ./ 1000.0,  OPR_Calc[:,idx], marker=MARKERS[idx+1], label=case_legend[idx])
end
plot!(p, Fn_N ./ 1000.0,  OPR_Ref, marker=MARKERS[1], label="Reference")
savefig(p,joinpath(save_dir,"$(save_key)_Overall_Pressure_Ratio.png"))

# Fan Diameter
p = plot(xlabel="Thrust [kN]", ylabel="Fan Diameter [m]", dpi=800, legend=:best)
for idx in 1:numTry
    plot!(p, Fn_N_Calc[:,idx] ./ 1000.0,  DFan_Calc_m[:,idx], marker=MARKERS[idx+1], label=case_legend[idx])
end
plot!(p, Fn_N ./ 1000.0,  ones(length(Fn_N)) .* Fan_Dia_Ref_m, marker=MARKERS[1], label="Reference")
savefig(p,joinpath(save_dir,"$(save_key)_Fan_Diameter.png"))

println("DONE")
