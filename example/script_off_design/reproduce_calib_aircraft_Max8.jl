"""
Reproduces a single sized aircraft model from engine calibration
"""

using TASOPT
using CSV, DataFrames
include(__TASOPTindices__)
include(joinpath(__TASOPTroot__, "../example/utilities_engines/run_engine.jl"))
include(joinpath(__TASOPTroot__, "../example/utilities_for_optimization/objective_factory.jl"))
include(joinpath(__TASOPTroot__, "../example/utilities_engines/calibrate_engine.jl"))
const min_fuse_radius = TASOPT.structures.find_minimum_radius_for_seats_per_row
const upd_fuse_pax! = TASOPT.structures.update_fuse_for_pax!

"""
    offdes_with_mbf_retry(ac, M0, P0, T0, a0, Fe_test;
                          retry_scales=[0.9999,0.999,0.998,0.995,0.99,0.98,0.95,0.9,1.3])

Wraps `RunEngine.runOffDes` with a fallback: if the cruise-seeded off-design solve fails to
converge, temporarily scale the aircraft's stored fan corrected-mass-flow initial guess
(`ac.pare[iembf,ipcruise1,1]`, what `runOffDes` reads its own initial guess from) by each of
`retry_scales` in turn and retry, restoring the original value afterward regardless of outcome.
For the 7791N idle-point EEDB convergence flake (2026-08-28), a fine-grained scan found the fix
is asymmetric: any *downward* nudge works, down to as little as -0.01% (`mbf*0.9999`), while every
upward nudge up to +5% fails (a larger jump, e.g. +30%, does eventually work upward too per an
earlier coarser sweep, kept here only as a last-resort fallback). Scales ordered smallest-to-
largest downward first so the retry stays as close to the original cruise-seeded guess as possible.
"""
function offdes_with_mbf_retry(ac, M0, P0, T0, a0, Fe_test;
                               retry_scales=[0.999999,0.99999,0.9999,0.999,0.998,0.995,0.99,0.98,0.95,0.9,1.3])
    res = RunEngine.runOffDes(ac, M0, P0, T0, a0, Fe_test; zero_offtake=true)
    res.Lconv && return res, 0
    mbf_orig = ac.pare[iembf, ipcruise1, 1]
    try
        for (i,s) in enumerate(retry_scales)
            ac.pare[iembf, ipcruise1, 1] = mbf_orig*s
            res = RunEngine.runOffDes(ac, M0, P0, T0, a0, Fe_test; zero_offtake=true)
            res.Lconv && return res, i
        end
    finally
        ac.pare[iembf, ipcruise1, 1] = mbf_orig
    end
    return res, -1 #never converged
end

#### Base warm-start airframe
ac_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_/Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_2400.jld2")
ac = quickload_aircraft(ac_dir)

#### Swap the structural material from the generic "TASOPT-Al" placeholder to a real, sourced alloy
#### (Al-2024-T4 -- the real baseline narrow-body material per Shukla et al. 2025)
ac.wing.inboard.caps.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=1.5)
ac.wing.outboard.caps.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=1.5)
ac.wing.inboard.webs.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=1.5)
ac.wing.outboard.webs.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=1.5)
ac.fuselage.skin.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=3.0)
ac.fuselage.cone.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=2.0)
ac.fuselage.bendingmaterial_h.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=1.5)
ac.fuselage.bendingmaterial_v.material = ac.fuselage.bendingmaterial_h.material
ac.fuselage.floor.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=1.5)

#### Reset to a 210-pax (737MAX8 FAA exit-limit) payload/fuselage, matching M78_210S
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

#### R1 (max-payload design range) for 737MAX8, and Mach 0.78 cruise (matching M78_210S)
range_R1_nmi = 2621.0 #[nmi]
ac.parm[imRange,1] = range_R1_nmi*1852.0 #[m]
ac.para[iaMach, ipclimbn:ipdescent1, 1] .= 0.78

#### Calibrated technology + engine-cycle/geometry parameters, from
#### test_engine_opt_LEAP1B_hard_ratio_M78_210S.jl's TechRun#200 (penalty 7.007138022178629).
#### hard_lpc_hpc_stage_ratio=true, numStageLC=3, numStageHC=10 (matching that calibration).
numStageLC = 3
numStageHC = 10
x_tec = [0.9571535292136226, 0.9331574934539477, 0.9136980491713441, 0.9184146215257285,
         0.9430533154646358, 0.9366603425373006, 0.9861806122272915, 1318.6185276152835]
x_eng = [7.695662307091953, 1.7210560910842305, 1.3265929460449386, 1521.6527309810594,
         0.6297346656329174, 24.797020200249143, 10.009251611881295]

#### Apply the technology level directly (matching UpdAcTecLvl!'s own unpacking)
pib,epolf,epollc,epolhc,epolht,epollt,etab,Tmetal = x_tec
ac.pare[iepib,   :, :] .= pib
ac.pare[ieepolf, :, :] .= epolf
ac.pare[ieepollc,:, :] .= epollc
ac.pare[ieepolhc,:, :] .= epolhc
ac.pare[ieepolht,:, :] .= epolht
ac.pare[ieepollt,:, :] .= epollt
ac.pare[ieetab,  :, :] .= etab
ac.parg[igTmetal]       = Tmetal

#### Apply the engine cycle+geometry and size the aircraft (UpdAcEngMod! calls size_aircraft! internally)
CaliEng.UpdAcEngMod!(ac, x_eng; iter_sizing=150, hard_lpc_hpc_stage_ratio=true,
                     numStageLC=numStageLC, numStageHC=numStageHC)

println("is_sized: ", ac.is_sized[1])
println("PFEI at design point: ", ac.parm[imPFEI,1])
println("WMTO (Ton): ", ac.parg[igWMTO]/gee/1000.0)
println("Wing AR: ", ac.wing.layout.AR)

#### EEDB off-design convergence check (log only, not part of any objective/penalty here)
miss_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/reference_engine_cycle_CFM5B/Target_Leap1B.csv")
df_eedb = CSV.read(miss_dir, DataFrame)
Fn_N_eedb = Float64.(df_eedb[:, "Thrust (kN)"] .* 1000.0)
M0 = 0.0; P0 = 101320.0; T0 = 288.2; a0 = 340.2074661144284

println()
println("--- EEDB off-design convergence check ---")
for Fn_N_cur in Fn_N_eedb
    res, retry_used = offdes_with_mbf_retry(ac, M0, P0, T0, a0, Fn_N_cur)
    status_str = retry_used == 0 ? "OK" : (retry_used == -1 ? "FAILED even after mbf retries" : "OK (needed mbf retry #$(retry_used))")
    println("Fn_N=$(Fn_N_cur/1000.0) kN -> Lconv=$(res.Lconv)  [$status_str]")
end

#### Save the resulting aircraft model
save_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/")
save_key = "Opti_Jet_NoACT_V4_CalibReproduced_Max8"
save_dir_actual = joinpath(save_dir,save_key*"_")
mkpath(save_dir_actual)
save_path = joinpath(save_dir_actual, "$(save_key)_$(round(Int,range_R1_nmi)).jld2")
quicksave_aircraft(ac, save_path)
println("Saved aircraft model to: ", save_path)
