using TASOPT
using CSV, DataFrames
include(__TASOPTindices__)
include(joinpath(__TASOPTroot__, "../example/utilities_engines/run_engine.jl"))
include(joinpath(__TASOPTroot__, "../example/utilities_for_optimization/objective_factory.jl"))
include(joinpath(__TASOPTroot__, "../example/utilities_engines/calibrate_engine.jl"))
include(joinpath(__TASOPTroot__, "utils/sensitivity.jl"))

ac_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Jet_NoACT_V3_14_R1Size_/Opti_Jet_NoACT_V3_14_R1Size_2400.jld2")
miss_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/reference_engine_cycle_CFM5B/Target_EEDB.csv")

ac = quickload_aircraft(ac_dir)
df = CSV.read(miss_dir, DataFrame)
Fn_N = df[:, "Thrust (kN)"] .* 1000.0 #[N]
BPR_Ref = df[:, "BPR"] #zeros(length(Fn_N))
OPR_Ref = df[:, "OPR"] #zeros(length(Fn_N))
Wfuel_Ref_kgs = df[:, "Wf[lbm/s]"] .* 0.453592 #zeros(length(Fn_N))

M0 = 0.0
P0 = 101320.0 #Pa
T0 = 288.2 #K
a0 = 340.2074661144284 #m/s

# for (i,Fn_N_cur) in enumerate(Fn_N)
#     res = RunEngine.runOffDes(ac,M0,P0,T0,a0,Fn_N_cur)
#     BPR_Ref[i] = res.BPR
#     OPR_Ref[i] = res.OPR
#     Wfuel_Ref_kgs[i] = res.mcore * res.ff #[kg/s]
# end

#### Design parameters
opt_parm = ObjectiveFactory.Parameter[]
push!(opt_parm,ObjectiveFactory.Parameter(:(pare[ieBPR,    ipclimb2:ipdescent4, 1]), 6.0, 12.0, 1.0, 0.5))
#
push!(opt_parm,ObjectiveFactory.Parameter(:(pare[iepif,    ipclimb2:ipdescent4, 1]), 2.0, 4.0, 1.25, 0.2))
push!(opt_parm,ObjectiveFactory.Parameter(:(pare[iepilc,   ipclimb2:ipdescent4, 1]), 3.0, 4.0,  2.0, 0.4))
push!(opt_parm,ObjectiveFactory.Parameter(:(pare[iepihc,   ipclimb2:ipdescent4, 1]), 10.0, 50.0, 1.25, 1.0))
push!(opt_parm,ObjectiveFactory.Parameter(:(pare[iepib,    ipclimb2:ipdescent4, 1]), 0.96, 0.98, 0.93, 0.005))
#
push!(opt_parm,ObjectiveFactory.Parameter(:(pare[ieepolf,  ipclimb2:ipdescent4, 1]), 0.90, 0.95, 0.80, 0.1))
push!(opt_parm,ObjectiveFactory.Parameter(:(pare[ieepollc, ipclimb2:ipdescent4, 1]), 0.90, 0.95, 0.80, 0.1))
push!(opt_parm,ObjectiveFactory.Parameter(:(pare[ieepolhc, ipclimb2:ipdescent4, 1]), 0.90, 0.95, 0.80, 0.1))
push!(opt_parm,ObjectiveFactory.Parameter(:(pare[ieepolht, ipclimb2:ipdescent4, 1]), 0.90, 0.95, 0.80, 0.1))
push!(opt_parm,ObjectiveFactory.Parameter(:(pare[ieepollt, ipclimb2:ipdescent4, 1]), 0.90, 0.95, 0.80, 0.1))
#
push!(opt_parm,ObjectiveFactory.Parameter(:(pare[ieetab,   ipclimb2:ipdescent4, 1]), 0.99, 0.999, 0.975, 0.002))
#
push!(opt_parm,ObjectiveFactory.Parameter(:(pare[ieTt4,    ipclimb2:ipdescent4, 1]), 1500.0, 2000.0, 1000.0, 100.0))

#### Extract initial guess
for (j,para_cur) in enumerate(opt_parm)
    # Get starting value from the warm start aircraft model
    try # Many of the variable (Need to check why cruise value were set to many phases)
        para_cur.val = getNestedProp_fromExpr(ac; parm_sym=para_cur.name)[5] #Hopefully get the cruise phase value if possible
    catch
        para_cur.val = getNestedProp_fromExpr(ac; parm_sym=para_cur.name)[1]
    end
    # Correct the local bounds be at least contained within global bounds
    para_cur.val = clamp(para_cur.val, para_cur.bon_lo, para_cur.bon_up)
end

# (; obj!, hist) = CaliEng.make_obj(ac, opt_parm, Fn_N,Wfuel_Ref_kgs, OPR_Ref, BPR_Ref;M0=M0, P0=P0, T0=T0, a0=a0)
# println(obj!([0.889,6.2379],nothing))

opt_parm_ori = deepcopy(opt_parm)

#### Optimization
status, hist, bestSol = CaliEng.optimize_match_EEDB!(ac,opt_parm,
                                                     Fn_N,Wfuel_Ref_kgs, 
                                                     OPR_Ref, BPR_Ref;
                                                     M0=M0, P0=P0, T0=T0, a0=a0,
                                                     max_iter_sizing=150,
                                                     max_iter_optim=1000000,
                                                     optimizer_type=:GN_CRS2_LM)

status, hist, bestSol = CaliEng.optimize_match_EEDB!(ac,opt_parm,
                                                     Fn_N,Wfuel_Ref_kgs, 
                                                     OPR_Ref, BPR_Ref;
                                                     M0=M0, P0=P0, T0=T0, a0=a0,
                                                     max_iter_sizing=150,
                                                     max_iter_optim=1000000,
                                                     optimizer_type=:LN_NELDERMEAD)

display(bestSol)
println(status)
for (j,para_cur) in enumerate(opt_parm_ori)
    println("Old values for $(para_cur.name): $(para_cur.val)")
end
for (j,para_cur) in enumerate(opt_parm)
    println("New values for $(para_cur.name): $(para_cur.val)")
end



