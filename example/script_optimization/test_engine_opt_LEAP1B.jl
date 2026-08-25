using TASOPT
using CSV, DataFrames
include(__TASOPTindices__)
include(joinpath(__TASOPTroot__, "../example/utilities_engines/run_engine.jl"))
include(joinpath(__TASOPTroot__, "../example/utilities_for_optimization/objective_factory.jl"))
include(joinpath(__TASOPTroot__, "../example/utilities_engines/calibrate_engine.jl"))

#### Reference aircraft and EEDB data, matching example/utilities_engines/test.jl lines 9-23
ac_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_/Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_2400.jld2")
miss_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/reference_engine_cycle_CFM5B/Target_Leap1B.csv")

ac = quickload_aircraft(ac_dir)
df = CSV.read(miss_dir, DataFrame)
Fn_N          = Float64.(df[:, "Thrust (kN)"] .* 1000.0)      #[N]
BPR_ref       = df[:, "BPR"]                                  #Union{Missing,Float64} -- only the top 5 thrust points are reported
OPR_ref       = df[:, "OPR"]                                  #Union{Missing,Float64} -- only the top 5 thrust points are reported
WFuel_kgs_ref = Float64.(df[:, "Wf[lbm/s]"] .* 0.453592)      #[kg/s]

M0 = 0.0
P0 = 101320.0 #Pa
T0 = 288.2 #K
a0 = 340.2074661144284 #m/s
DFan_m_ref = 1.76 #m

#### Bounds for the 5 inner-loop cycle-design variables [BPR, pif, pilc, pihc, Tt4],
#### matching example/utilities_engines/test.jl lines 34-49 (Parameter(name, val, bon_up, bon_lo, d_val))
loBon_eng = [1.0,  1.25, 1.25, 1.25, 1000.0]
upBon_eng = [12.0, 4.0,  10.0, 60.0, 2000.0]

#### Bounds for the 8 outer-loop technology parameters [pib,epolf,epollc,epolhc,epolht,epollt,etab,Tmetal],
#### matching example/utilities_engines/test.jl lines 39-51
loBon_tec = [0.93, 0.87, 0.87, 0.87, 0.87, 0.87, 0.975, 1000.0]
upBon_tec = [0.98, 0.92, 0.92, 0.92, 0.92, 0.92, 0.999, 1330.0]

#### Initial guess: pull each parameter's actual cruise-point value from the reference aircraft, clamped to bounds
ini_eng_raw = [ac.pare[ieBPR,ipcruise1,1], ac.pare[iepif,ipcruise1,1], ac.pare[iepilc,ipcruise1,1],
               ac.pare[iepihc,ipcruise1,1], ac.pare[ieTt4,ipcruise1,1]]
ini_eng = clamp.(ini_eng_raw, loBon_eng, upBon_eng)

ini_tec_raw = [ac.pare[iepib,ipcruise1,1], ac.pare[ieepolf,ipcruise1,1], ac.pare[ieepollc,ipcruise1,1],
               ac.pare[ieepolhc,ipcruise1,1], ac.pare[ieepolht,ipcruise1,1], ac.pare[ieepollt,ipcruise1,1],
               ac.pare[ieetab,ipcruise1,1], ac.parg[igTmetal]]
ini_tec = clamp.(ini_tec_raw, loBon_tec, upBon_tec)

println("Initial guess [BPR, pif, pilc, pihc, Tt4]: ", ini_eng)
println("Initial guess [pib, epolf, epollc, epolhc, epolht, epollt, etab, Tmetal]: ", ini_tec)

#### Run the full double-loop technology calibration (outer: match LEAP-1B EEDB data; inner: minimum-TSFC cycle)
status, bestSol_tec, bestSol_eng, histTechPara, histEngPara, histPenl =
    CaliEng.tech_opt(ac; ini_tec=ini_tec, ini_eng=ini_eng,
                      upBon_tec=upBon_tec, upBon_eng=upBon_eng,
                      loBon_tec=loBon_tec, loBon_eng=loBon_eng,
                      Fn_N=Fn_N, WFuel_kgs_ref=WFuel_kgs_ref, OPR_ref=OPR_ref, BPR_ref=BPR_ref, DFan_m_ref=DFan_m_ref,
                      M0=M0, P0=P0, T0=T0, a0=a0,
                      printEvery=10, ftol_tec=1e-6, ftol_eng=1e-7, tol_coupling=1e-8, maxIter=50)

println()
println("Status: ", status)
println("Best Tech Parameters [pib, epolf, epollc, epolhc, epolht, epollt, etab, Tmetal]: ", bestSol_tec)
println("Best Engine Parameters [BPR, pif, pilc, pihc, Tt4]: ", bestSol_eng)
println("Number of fully-feasible evaluations recorded: ", length(histPenl))
