using TASOPT
using CSV, DataFrames
using JLD2
include(__TASOPTindices__)
include(joinpath(__TASOPTroot__, "../example/utilities_engines/run_engine.jl"))
include(joinpath(__TASOPTroot__, "../example/utilities_for_optimization/objective_factory.jl"))
include(joinpath(__TASOPTroot__, "../example/utilities_engines/calibrate_engine.jl"))
const min_fuse_radius = TASOPT.structures.find_minimum_radius_for_seats_per_row
const upd_fuse_pax! = TASOPT.structures.update_fuse_for_pax!

#### Reference aircraft and EEDB data, matching example/utilities_engines/test.jl lines 9-23
ac_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_/Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_2400.jld2")
miss_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/reference_engine_cycle_CFM5B/Target_Leap1B.csv")
save_nam = "tech_opt_LEAP1B_hard_ratio_result_Max8.jld2" #Save under current run folder
ac = quickload_aircraft(ac_dir)

#### Swap the structural material from the generic "TASOPT-Al" placeholder to a real, sourced alloy
#### (Al-2024-T4 -- the real baseline narrow-body material per Shukla et al. 2025, and what's
#### actually flying on the 737 family), uniformly across wing+fuselage (matches the paper's own
#### narrow-body treatment, and TASOPT's own structural model has no per-surface granularity to do
#### otherwise). Safety factors match read_input.jl's per-component values, not the generic default.
ac.wing.inboard.caps.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=1.5)
ac.wing.outboard.caps.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=1.5)
ac.wing.inboard.webs.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=1.5)
ac.wing.outboard.webs.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=1.5)
ac.fuselage.skin.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=3.0)
ac.fuselage.cone.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=2.0)
ac.fuselage.bendingmaterial_h.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=1.5)
ac.fuselage.bendingmaterial_v.material = ac.fuselage.bendingmaterial_h.material
ac.fuselage.floor.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=1.5)

#### Reset to a 210-pax (737MAX8 FAA exit-limit) payload/fuselage, same direct-overwrite pattern as
#### opt_from_multi_warm_starts_para_oag_Rerun.jl:182-193 / the R1Size MAX8 scripts.
num_pax = 210
wei_per_pass_N = 956.36773
num_pass_row = Int(6)
ac.parm[imWperpax,:] .= wei_per_pass_N
ac.parg[igWpaymax] = wei_per_pass_N*num_pax
ac.parm[imWpay,:] .= ac.parg[igWpaymax]*(1-1e-10)
ac.fuselage.cabin.exit_limit = num_pax
ac.fuselage.APU.W = 0.035*ac.parg[igWpaymax]
ac.fuselage.seat.W = 0.10*ac.parg[igWpaymax]
ac.fuselage.added_payload.W = 0.35*ac.parg[igWpaymax]
fuse_radius = min_fuse_radius(num_pass_row, ac)
ac.fuselage.layout.cross_section.radius = fuse_radius
ac.fuselage.cabin.front_seat_offset = 0.0
ac.fuselage.cabin.rear_seat_offset = 0.0
upd_fuse_pax!(ac)

#### R1 (max-payload design range) for 737MAX8, read off Boeing's own payload/range chart
#### (§3.2.3, "Payload/Range for Long Range Cruise: Model 737-8 / -8-200")
range_R1_nmi = 2621.0 #[nmi]
ac.parm[imRange,1] = range_R1_nmi*1852.0 #[m]
ac.para[iaMach, ipclimbn:ipdescent1, 1] .= 0.78

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

#### Hard-fix the LPC/HPC split at a physical 3-stage-LPC/10-stage-HPC ratio (equal pressure ratio
#### per stage), matching LEAP-1B's actual architecture. The inner-loop cycle search variable becomes
#### [BPR, pif, per_stage, Tt4] (4 elements) instead of [BPR, pif, pilc, pihc, Tt4] (5 elements) -- see
#### UpdAcEngMod!'s docstring in calibrate_engine.jl.
numStageLC = 3
numStageHC = 10

#### Bounds for the 7 inner-loop cycle+geometry design variables
#### [BPR, pif, per_stage, Tt4, CL_cruise, Sweep_wing, AR_wing].
#### per_stage bounds derived from the original pilc in [1.25,10] / pihc in [1.25,60] bounds'
#### OPR extremes: per_stage_lo = (1.25*1.25)^(1/13), per_stage_hi = (10*60)^(1/13).
#### CL_cruise/Sweep_wing/AR_wing bounds match the R1Size scripts' bound_glob (sensitivity-analysis-
#### justified -- see tasopt-engine-calibration-double-loop memory).
loBon_eng = [1.0,  1.25, 1.03, 1000.0, 0.3,  0.0,  5.0]
upBon_eng = [16.0, 4.0,  1.65, 2500.0, 1.00, 60.0, 20.0]

#### Bounds for the 8 outer-loop technology parameters [pib,epolf,epollc,epolhc,epolht,epollt,etab,Tmetal],
#### matching example/utilities_engines/test.jl lines 39-51
loBon_tec = [0.93, 0.86, 0.86, 0.86, 0.86, 0.86, 0.975, 1000.0]
upBon_tec = [0.98, 0.96, 0.96, 0.96, 0.96, 0.96, 0.999, 1500.0]

#### Design constraints for the inner engine-cycle search, checked against the sized aircraft
constraints = ObjectiveFactory.Constraint[
    ObjectiveFactory.Constraint(:(wing.layout.span);        lim_up=35.814,  pen_sca=2000.0),
    ObjectiveFactory.Constraint(:(parm[imlBF,1]);           lim_up=2400.0,  pen_sca=2000.0),
    ObjectiveFactory.Constraint(:(para[iagamV,ipclimbn,1]); lim_lo=0.015,   pen_sca=2000.0),
    ObjectiveFactory.Constraint(:(pare[ieTt3,:,1]);         lim_up=900.0,   pen_sca=2000.0),
    ObjectiveFactory.Constraint(:(pare[ieTmet1,:,1]);       lim_up=:(parg[igTmetal]), pen_sca=2000.0, eps_buff=1e-4),
    ObjectiveFactory.Constraint(:(parg[igdfan]);            lim_up=2.0,     pen_sca=2000.0),
    ObjectiveFactory.Constraint(:(parm[imWTO,1]);           lim_up=:(parg[igWMTO]),  pen_sca=2000.0, eps_buff=1e-4),
    ObjectiveFactory.Constraint(:(parm[imVfuel,1]);         lim_up=:(parg[igVfmax]), pen_sca=2000.0, eps_buff=1e-4),
]

#### Warm start: best known solution from tech_opt_LEAP1B_result_96MaxEff_CheckNegRateOfClimb.out
#### (TechRun#225, penalty=5.806405738287389, best-so-far unchanged since TechRun#180).
# ini_tec = [0.954855687966321, 0.9278477078626518, 0.9114766770119112, 0.9138021445004528,
#            0.9301470893763781, 0.9342864151766557, 0.984120884837075, 1306.8935808058604]
ini_tec = [0.9558329480883324, 0.9262420136229457, 0.907208300080772, 0.912959447456788,
           0.9429385982757191, 0.936171198797011, 0.9862773853225691, 1316.6276828021078]
#### That run's best [BPR,pif,pilc,pihc,Tt4] = [7.574856839182607, 1.7261951245581042,
#### 1.351869322008949, 31.289835401251533, 1583.707990298101]. Converted to [BPR,pif,per_stage,Tt4]
#### by taking OPR=pilc*pihc and redistributing across all 13 stages at equal pressure ratio per
#### stage: per_stage = OPR^(1/13) = 1.3338374976528422 (pilc^(1/3)==pihc^(1/10) both give this).
#### CL_cruise/Sweep_wing/AR_wing have no historical warm-start value (didn't exist as free variables
#### in that run) -- start them at the R1Size scripts' own bound_glob defaults instead.
# ini_eng = [7.574856839182607, 1.7261951245581042, 1.3338374976528422, 1583.707990298101,
#            0.6, 30.0, 10.0]
ini_eng = [8.165480728072914, 1.6358960754233371, 1.3241531220099012, 1524.7148624299641,
           0.5738509289814249, 28.44292299192375, 10.220177699402463]

ini_tec = clamp.(ini_tec, loBon_tec, upBon_tec)
ini_eng = clamp.(ini_eng, loBon_eng, upBon_eng)

println("Initial guess [BPR, pif, per_stage, Tt4, CL_cruise, Sweep_wing, AR_wing]: ", ini_eng)
println("Initial guess [pib, epolf, epollc, epolhc, epolht, epollt, etab, Tmetal]: ", ini_tec)

#### Reference values for the two new outer-loop matching terms (737MAX8, highest-weight-option
#### MTOW per Boeing's own convention -- see tasopt-r1size-reopt-vs-737max10 memory; wing AR shared
#### across the whole 737 MAX family)
MWTO_ref = 82.6*1000.0*gee #[N]
ARwing_ref = 10.16

#### Run the full double-loop technology calibration (outer: match LEAP-1B EEDB data + MWTO/wing-AR;
#### inner: minimum-PFEI cycle+geometry) with the 3-stage/10-stage LPC/HPC ratio hard-enforced.
status, bestSol_tec, bestSol_eng, histTechPara, histEngPara, histPenl =
    CaliEng.tech_opt(ac; ini_tec=ini_tec, ini_eng=ini_eng,
                     upBon_tec=upBon_tec, upBon_eng=upBon_eng,
                     loBon_tec=loBon_tec, loBon_eng=loBon_eng,
                     Fn_N=Fn_N, WFuel_kgs_ref=WFuel_kgs_ref, OPR_ref=OPR_ref, BPR_ref=BPR_ref, DFan_m_ref=DFan_m_ref,
                     MWTO_ref=MWTO_ref, ARwing_ref=ARwing_ref,
                     M0=M0, P0=P0, T0=T0, a0=a0,
                     printEvery=5, ftol_tec=1e-6, ftol_eng=1e-7, iter_sizing=150, maxIter=5000, optTyp=:LN_NELDERMEAD,
                     constraints=constraints,
                     hard_lpc_hpc_stage_ratio=true, numStageLC=numStageLC, numStageHC=numStageHC)

println()
println("Status: ", status)
println("Best Tech Parameters [pib, epolf, epollc, epolhc, epolht, epollt, etab, Tmetal]: ", bestSol_tec)
println("Best Engine Parameters [BPR, pif, per_stage, Tt4, CL_cruise, Sweep_wing, AR_wing]: ", bestSol_eng)
println("Number of fully-feasible evaluations recorded: ", length(histPenl))

#### Save results for now (quick JLD2 dump, not a formal output convention yet)
result_savepath = joinpath(__TASOPTroot__, "../example/script_optimization/$(save_nam)")
@save result_savepath status bestSol_tec bestSol_eng histTechPara histEngPara histPenl
println("Saved results to: ", result_savepath)