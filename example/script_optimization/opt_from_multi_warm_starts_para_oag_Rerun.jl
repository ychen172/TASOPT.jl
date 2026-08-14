"""
This script performs further optimization of the OAG seat-capacity sweep, warm-started from a
previous generation's own results (opt_from_multi_warm_starts_para_oag.jl or an earlier run of
this script) instead of the single fixed 230-passenger baseline.
Meant to be invoked by hand, once per generation: bump `this_save_key` (and point `prev_save_key`
at whatever you just produced) each time you rerun. The intended cycle is:
    1. Run this script (a "generation")
    2. screen_out_best_cases_oag.jl to pick the best case at each seat capacity so far
    3. compare_cases_oag.jl to sanity check the result
    4. Rerun this script again, chaining from the previous generation

`flag_cross_warm_start` toggles the warm-start source:
    - true  ("cross"): seed each seat capacity from an *adjacent* capacity's previous result
                       (alternating which side, by parity of the capacity's position), never itself.
                       Intended to avoid every case just re-converging to the same local optimum it
                       already had.
    - false ("self"):  seed each seat capacity from its *own* previous result, if it has one.
In both modes, if the preferred source doesn't have a saved model from the previous generation
(most notably: a seat capacity that has never converged at all, e.g. 340+), the warm-start source
falls back to the nearest seat capacity (by index distance, either direction, unlimited range) that
does have one. This is what lets each generation push the convergence boundary outward by one step
at a time, while every generation still attempts all seat capacities (so the case count never shrinks).
Cases with no realistic chance of converging (e.g. 340 while still Type-C wingspan constrained) are
expected to fail fast (NLopt's own ftol_rel stops a flat/failing local search after a handful of
evaluations - see the empirical bootstrap investigation), so re-attempting them every generation is cheap.
Use job array from orcd cluster to process multiple cases at the same time
Assume taskID starts from 1 instead of 0
Run by: sbatch ../../julia-jobArr-orcd.sh
"""

using Glob
using TASOPT
include(__TASOPTindices__)
include(joinpath(__TASOPTroot__,"utils/sensitivity.jl"))
include(joinpath(__TASOPTroot__, "../example/utilities_for_optimization/objective_factory.jl"))
using .ObjectiveFactory: Constraint, Parameter, OptHistory, Requirement, optimize_singlePt_PFEI!, optimizer_wrapper_global_local, size_aircraft_w_param!,
                         save_vec_struct_csv, load_csv_constraints, load_csv_parameters, load_csv_requirements, save_jld2, load_jld2, OffDesMission, ObjectiveVariable
include(joinpath(__TASOPTroot__,"../example/utilities_for_postprocessing/Extract.jl"))
const min_fuse_radius = TASOPT.structures.find_minimum_radius_for_seats_per_row
const upd_fuse_pax! = TASOPT.structures.update_fuse_for_pax!
const success_statuses = ObjectiveFactory.success_statuses

#### Optimization parameters
# Previous generation to warm-start from, and the key for this generation's output.
# For the very first rerun, prev_save_key should point at the original opt_from_multi_warm_starts_para_oag.jl output.
save_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/")
prev_save_key = "Opti_Jet_NoACT_OAG_6Seats_TypeC" # Bump these two by hand each invocation
this_save_key = "Opti_Jet_NoACT_OAG_6Seats_TypeC_V2"
flag_cross_warm_start = true # true: adjacent-neighbor cross warm-start. false: self warm-start (own previous result)
# Mission extraction directory
miss_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/OAG_Data_2024/OAG_Data_2024.csv")
# Optimization configuration parameters
flag_skip_global = true #Switch on to skip the global search if confident that warm start model can converge
max_iter_sizing = 150 # Maximum iterations for TASOPT sizing
optimizers = [:GN_CRS2_LM, :LN_NELDERMEAD] # Optimizer choice. [Global,Local]
max_num_iter_opt = [30000, 800, 10000] # Maximum number of optimization steps [Global, Local Coarse, Local Fine]
max_num_round_loc = [120, 30] #Maximum number of adaptive bounds refinements rounds for local search [Local Coarse, Local Fine]
span_glo_to_loc = 0.25 # span_local_search/span_global_search
# Misc
num_seats_per_row = 6
wei_per_pass_N = 956.36773 #Weight per passenger [N] (Assume a constant APU, seat, and added weight fractions)
pass_load_frac_off = 0.825 #Load factor of passengers
idx_fuel = 24 #Eth: 32, Jet: 24 #Assume off-design use the same fuel
rho_fuel_kgm3 = 817.0 #Eth: 789.0, Jet: 817.0 #kg/m3
hvap_fuel_Jkg = 358694.0 #Eth: 918187.9, Jet: 358694.0 #J/kg
objVar_des=ObjectiveVariable(:(parm[imWfuel,1]))
objVar_off=ObjectiveVariable(:(parm[imWfuel,2]))
pen_scale_PFEI = 2e5 #The scaling for penalty when using fuel burned per mission compared to the case that use PFEI
pen_scale_constraints = 1e4*pen_scale_PFEI
pen_scale_failed_sizing = 100.0*pen_scale_PFEI
# Setup an entry (value does not matter) for test ratio to be sweep later in optimization loop
mis_opt = Requirement[]
push!(mis_opt, Requirement(:(options.ifuel), Int(idx_fuel)))
push!(mis_opt, Requirement(:(parg[igrhofuel]), rho_fuel_kgm3))
push!(mis_opt, Requirement(:(pare[iehvap, :, 1]), hvap_fuel_Jkg))
push!(mis_opt, Requirement(:(pare[iehvapcombustor, :, 1]), hvap_fuel_Jkg))
push!(mis_opt, Requirement(:(options.has_ACT_fuel), false)) # Whether to allow additional center fuel tank (ACT)
push!(mis_opt, Requirement(:(options.compensate_ACT), false)) # Whether to increase the aircraft length to accmondate for the cargo space taken by ACT
push!(mis_opt, Requirement(:(fuse_tank.ACT_eta_vol), 1.00)) # Volumetric efficiency of ACT fuel tank
push!(mis_opt, Requirement(:(fuse_tank.ACT_eta_wei), 1.00)) # Gravimetric efficiency of ACT fuel tank
push!(mis_opt, Requirement(:(para[iaMach, ipclimbn:ipdescent1, 1]), 0.8)) # Cruise Mach number
push!(mis_opt, Requirement(:(vtail.opt_sizing), TailSizing.OEI)) # Vertical tail sizing mode (Engine out or fixed volume)
push!(mis_opt, Requirement(:(parg[igCLveout]), 0.47)) # Vertical tail CL at engine out condition
push!(mis_opt, Requirement(:(htail.opt_sizing), TailSizing.CLmaxFwdCG)) # Horizontal tail sizing mode (Forward CG or fixed volume)
push!(mis_opt, Requirement(:(htail.CL_max_fwd_CG), -0.7)) # Horizontal tail CL at forward CG condition
# Constraints for this optimization
con_opt = Constraint[]
push!(con_opt, Constraint(:(wing.layout.span); pen_sca=pen_scale_constraints, lim_up=35.814)) # [m] Type C wing span constraint
push!(con_opt, Constraint(:(parm[imlBF, 1]); pen_sca=pen_scale_constraints, lim_up=2400.0)) # [m] Maximum balanced field length for takeoff
push!(con_opt, Constraint(:(para[iagamV, ipclimbn, 1]); pen_sca=pen_scale_constraints, lim_lo=0.015)) # [rad] Minimum cruise climb angle at TOC
push!(con_opt, Constraint(:(pare[ieTt3, :, 1]); pen_sca=pen_scale_constraints, lim_up=900.0)) # [K] Maximum compressor outlet temperature
push!(con_opt, Constraint(:(pare[ieTmet1, :, 1]); pen_sca=pen_scale_constraints, lim_up=1333.33)) # [K] Maximum turbine metal temperature
push!(con_opt, Constraint(:(parg[igdfan]); pen_sca=pen_scale_constraints, lim_up=3.0)) # [m] Maximum fan diameter
push!(con_opt, Constraint(:(parm[imWTO, 1]); pen_sca=pen_scale_constraints, lim_up=:(parg[igWMTO]), eps_buff=1e-4)) # [N] Maximum takeoff weight
push!(con_opt, Constraint(:(parm[imVfuel, 1]); pen_sca=pen_scale_constraints, lim_up=:(parg[igVfmax]), eps_buff=1e-4)) # [m3] Maximum fuel volume
# Parameters to be optimized (Initial guess(Only for local search), Upper bound, Lower bound, Initial Step Size(Only for local search))
bound_glob = Parameter[] #Step size = 1/5 of local search span = 1/5 * 1/4 * global search span specified below
push!(bound_glob, Parameter(:(wing.layout.sweep), 30.0, 60.0, 0.0, 5.0)) # [deg] Wing sweep angle
push!(bound_glob, Parameter(:(wing.layout.AR), 10.0, 20.0, 5.0, 1.5)) # Wing aspect ratio
push!(bound_glob, Parameter(:(wing.inboard.cross_section.thickness_to_chord), 0.2, 0.6, 0.04, 0.05)) # Wing inner thickness to chord
push!(bound_glob, Parameter(:(wing.outboard.cross_section.thickness_to_chord), 0.2, 0.6, 0.04, 0.05))
push!(bound_glob, Parameter(:(wing.inboard.λ), 0.7, 1.0, 0.1, 0.1)) # Wing inner taper ratio
push!(bound_glob, Parameter(:(wing.outboard.λ), 0.3, 1.0, 0.1, 0.1))
push!(bound_glob, Parameter(:(para[iarcls, ipclimb2:ipdescent4, 1]), 1.0, 2.0, 0.4, 0.15)) # Wing cl ratio span break to root at cruise
push!(bound_glob, Parameter(:(para[iarclt, ipclimb2:ipdescent4, 1]), 1.0, 2.0, 0.4, 0.15)) # Wing cl ratio of tip to root at cruise
push!(bound_glob, Parameter(:(para[iaCL, ipclimb2:ipdescent4, 1]), 0.6, 1.00, 0.3, 0.07)) # Wing CL at cruise
push!(bound_glob, Parameter(:(para[iaalt, ipcruise1, 1]), 10000.0, 20000.0, 4000.0, 1600.0)) # [m] Cruise altitude
push!(bound_glob, Parameter(:(pare[iepif, ipclimb2:ipdescent4, 1]), 2.0, 4.0, 1.25, 0.2)) # Fan PR at cruise
push!(bound_glob, Parameter(:(pare[iepilc, ipclimb2:ipdescent4, 1]), 3.0, 10.0, 1.25, 0.4)) # LPC PR at cruise
push!(bound_glob, Parameter(:(pare[iepihc,ipclimb2:ipdescent4, 1]), 10.0, 50.0, 1.25, 1.0)) # HPC PR at cruise
push!(bound_glob, Parameter(:(pare[ieBPR, ipclimb2:ipdescent4, 1]), 8.0, 30.0, 1.0, 2.0)) # Fan BPR at cruise
push!(bound_glob, Parameter(:(pare[ieTt4, ipclimb2:ipdescent4, 1]), 1500.0, 2000.0, 1000.0, 100.0)) #[K] Turbine inlet temperature  at cruise
push!(bound_glob, Parameter(:(vtail.layout.AR), 2.0, 5.0, 1.0, 0.2)) #Vertical tail aspec ratio

#### Setup job-array for ORCD cluster run
@assert length(ARGS) >= 2 "Usage: julia xxx.jl task_id num_tasks"
task_id = parse(Int,ARGS[1]) #Current task id
num_tasks = parse(Int,ARGS[2]) #Total number of tasks
@assert 1 <= task_id <= num_tasks

#### Setup the previous-generation (read) and this-generation (write) save directories
prev_save_dir_actual = joinpath(save_dir,prev_save_key*"_")
save_dir_actual = joinpath(save_dir,this_save_key*"_")
mkpath(save_dir_actual)

#### Load the mission information
miss_off_des = Extract.read_oag(miss_dir)
seat_cap_keys = sort(collect(keys(miss_off_des)))
n_caps = length(seat_cap_keys)

#### Resolve which seat-capacity index to warm-start capacity index `i` from.
# prefer_cross=true : try the adjacent capacity first (alternating side by parity of i), never self.
# prefer_cross=false: try capacity i's own previous result first ("self").
# In both cases, if the preferred source has no saved model in the previous generation, fall back to
# the nearest seat capacity (by index distance, either direction) that does - unlimited range, erroring
# only if literally nothing in the previous generation converged.
function resolve_warm_start_idx(i::Int, n::Int, has_result::Function; prefer_cross::Bool=true)
    if !prefer_cross && has_result(i)
        return i
    end
    prefer_lower = iseven(i) # arbitrary but deterministic tie-break, reused at every search radius
    for d in 1:(n-1)
        cand_lower = i - d
        cand_upper = i + d
        order = prefer_lower ? (cand_lower, cand_upper) : (cand_upper, cand_lower)
        for cand in order
            (1 <= cand <= n) && has_result(cand) && return cand
        end
    end
    error("No seat capacity in the previous generation ($(prev_save_key)) has a saved model to warm-start capacity index $(i) from.")
end

#### Optimization
for i in task_id:num_tasks:n_caps
    sc = seat_cap_keys[i]

    #### Create inidividual log file for each task
    status_log = joinpath(save_dir_actual, "$(this_save_key)_$(round(Int,sc))_OptLog.txt")
    open(status_log, "w") do io
        println(io, "seat_capacity,warm_start_from_seat_capacity,Status")
    end

    #### Make copies for each task
    mis_opt_cur = deepcopy(mis_opt)
    con_opt_cur = deepcopy(con_opt)
    bound_glob_cur = deepcopy(bound_glob)

    #### Resolve and load the warm-start aircraft from the previous generation
    has_result = idx -> isfile(joinpath(prev_save_dir_actual, "$(prev_save_key)_$(round(Int,seat_cap_keys[idx])).jld2"))
    src_idx = resolve_warm_start_idx(i, n_caps, has_result; prefer_cross=flag_cross_warm_start)
    sc_src = seat_cap_keys[src_idx]
    println("Seat capacity $(sc): warm-starting from seat capacity $(sc_src) ($(src_idx==i ? "self" : "cross"))")
    ac = quickload_aircraft(joinpath(prev_save_dir_actual, "$(prev_save_key)_$(round(Int,sc_src)).jld2"))
    @assert ac.is_sized[1] #Make sure mission 1 is sized

    #### Load off-design missions
    ranges_off_nmi = miss_off_des[sc].ranges_nmi #Off-design ranges (nmi)
    weights_off = miss_off_des[sc].weights #Off-design missions' weightings (add to 1)

    #### Update the fueselage radius and layout
    ac.parm[imWperpax,:] .= wei_per_pass_N #Per-passenger weight update (N)
    ac.parg[igWpaymax] = wei_per_pass_N*sc #Maximum payload setup (N)
    ac.parm[imWpay,:] .= ac.parg[igWpaymax] #This will force the sizing mission to be R1 of the aircraft
    ac.fuselage.cabin.exit_limit = sc #Exit limit same as the maximum single-class seat layout (num pass)
    ac.fuselage.APU.W = 0.035*ac.parg[igWpaymax]
    ac.fuselage.seat.W = 0.10*ac.parg[igWpaymax]
    ac.fuselage.added_payload.W = 0.35*ac.parg[igWpaymax]
    fuse_radius = min_fuse_radius(round(Int,num_seats_per_row), ac) #Minimum fueslage radius given the number of seats per row
    ac.fuselage.layout.cross_section.radius = fuse_radius
    ac.fuselage.cabin.front_seat_offset = 0.0  # [m] Front seat offset from the front of cylinder
    ac.fuselage.cabin.rear_seat_offset = 0.0 # [m] Back seat offset from the back of cylinder
    upd_fuse_pax!(ac) #Update carbin geometry

    #### Use the fuselage radius to update the mission and bound requirements
    # Mission requirement
    push!(mis_opt_cur, Requirement(:(parg[igxeng]), ac.wing.layout.box_x - fuse_radius)) # [m] Engine axial location from the aircraft nose tip
    push!(mis_opt_cur, Requirement(:(parg[igdxeng2wbox]), fuse_radius)) # [m] Corresponding to above for wing box location scaling
    # Optimization parameters
    push!(bound_glob_cur, Parameter(:(parg[igyeng]), 3.0*fuse_radius, 10.0*fuse_radius, 1.75*fuse_radius, 0.2*fuse_radius)) # [m] Span wise engine location
    range_max_off_nmi = maximum(ranges_off_nmi)
    range_min_off_nmi = minimum(ranges_off_nmi)
    (range_max_off_nmi>10.0 && range_min_off_nmi>10.0 && range_max_off_nmi>range_min_off_nmi) || error("OAG minimum input range is required to be larger than 10 nmi")
    range_span_off_nmi = range_max_off_nmi-range_min_off_nmi
    push!(bound_glob_cur, Parameter(:(parm[imRange,1]), range_max_off_nmi*1852.0, (range_max_off_nmi+0.4*range_span_off_nmi)*1852.0, max(10.0,(range_min_off_nmi-0.2*range_span_off_nmi))*1852.0, 0.1*range_span_off_nmi*1852.0)) #[m]
    ac.parm[imRange,:] .= range_max_off_nmi*1852.0 # This is setup so that if use warm start, the sizing mission still stay with the sepecified starting point above

    #### If not running global search, then fetch the warm start variables from the initial aircraft model
    par_opt_cur = deepcopy(bound_glob_cur) #Use directly global bound if starting from global search
    if flag_skip_global
        for (j,para_cur) in enumerate(par_opt_cur)
            # Get starting value from the warm start aircraft model
            try # Many of the variable (Need to check why cruise value were set to many phases)
                para_cur.val = getNestedProp_fromExpr(ac; parm_sym=para_cur.name)[5] #Get the start-of-cruise phase value for those slices setup from climb2
            catch
                para_cur.val = getNestedProp_fromExpr(ac; parm_sym=para_cur.name)[1] #Works for both number and single element vector
            end
            # Correct the local bounds be at least contained within global bounds
            bon_up = para_cur.bon_up
            bon_lo = para_cur.bon_lo
            para_cur.val = clamp(para_cur.val, bon_lo, bon_up)
            (span_glo_to_loc<1.0) || error("Local bound span needs to be smaller than global")
            dx = (bon_up - bon_lo)*span_glo_to_loc
            para_cur.bon_up = para_cur.val+0.5*dx
            para_cur.bon_lo = para_cur.val-0.5*dx
            if para_cur.bon_up>bon_up
                para_cur.bon_up = bon_up
                para_cur.bon_lo = bon_up - dx
            elseif para_cur.bon_lo<bon_lo
                para_cur.bon_lo = bon_lo
                para_cur.bon_up = bon_lo + dx
            end
            para_cur.d_val = para_cur.d_val*span_glo_to_loc
        end
    end

    #### Save the mission requirement
    save_vec_struct_csv(joinpath(save_dir_actual, "$(this_save_key)_$(round(Int,sc))_mission_requirements.csv"), mis_opt_cur)
    save_vec_struct_csv(joinpath(save_dir_actual, "$(this_save_key)_$(round(Int,sc))_design_constraints.csv"), con_opt_cur)
    save_vec_struct_csv(joinpath(save_dir_actual, "$(this_save_key)_$(round(Int,sc))_global_bounds.csv"), bound_glob_cur)

    #### Prepare off-design missions
    num_miss_off = length(ranges_off_nmi)
    range_off_des_nmi = ranges_off_nmi #[nmi]
    wei_pay_off_des_N = fill(ac.parg[igWpaymax]*pass_load_frac_off, num_miss_off) #[N]
    idx_fuel_off_des = fill(Int(idx_fuel), num_miss_off)
    rho_fuel_off_des_kgm3 = fill(rho_fuel_kgm3, num_miss_off)
    hvap_fuel_off_des_Jkg = fill(hvap_fuel_Jkg, num_miss_off)
    off_des_miss = OffDesMission{Float64}(range_off_des_nmi, wei_pay_off_des_N, idx_fuel_off_des, rho_fuel_off_des_kgm3, hvap_fuel_off_des_Jkg)
    off_des_constraints = [:WPay,:MWTO,:VolFuel]
    PFEI_Weighting = [0.0;weights_off]

    #### Run the global local optimization process
    par_opt_found, status_found, hist_found =
    optimizer_wrapper_global_local(ac, par_opt_cur, bound_glob_cur;
                                   objVar_des=objVar_des, objVar_off=objVar_off, pen_failed_sizing=pen_scale_failed_sizing,
                                   miss_req=mis_opt_cur, constraints=con_opt_cur, max_iter_sizing=max_iter_sizing,
                                   off_des_miss=off_des_miss, off_des_constraints=off_des_constraints, PFEI_Weighting=PFEI_Weighting,
                                   optimizer_global=optimizers[1], max_iter_glo=max_num_iter_opt[1],   span_glo_to_loc=span_glo_to_loc, run_global=!flag_skip_global,
                                   optimizer_local =optimizers[2], max_iter_loc_C=max_num_iter_opt[2], max_round_loc_C=max_num_round_loc[1],
                                                                   max_iter_loc_F=max_num_iter_opt[3], max_round_loc_F=max_num_round_loc[2])

    ##### Post-process the optimization results
    if status_found in success_statuses
        # attempt to update an aircraft model to save
        ac_copy = deepcopy(ac)
        sized_succeeded = size_aircraft_w_param!(ac_copy; mission_req=mis_opt_cur, parameters=par_opt_found, max_iter_sizing=max_iter_sizing)
        # correct the status
        status_found = sized_succeeded ? status_found : (:FAILED_TO_REPRODUCE)
        # Only save the new model if successfuly sized
        if sized_succeeded
            quicksave_aircraft(ac_copy, joinpath(save_dir_actual, "$(this_save_key)_$(round(Int,sc)).jld2"))
        end
    end

    #### save the results
    open(status_log, "a") do io
        println(io, "$(round(Int,sc)),$(round(Int,sc_src)),$(string(status_found))")
    end
    save_vec_struct_csv(joinpath(save_dir_actual, "$(this_save_key)_$(round(Int,sc))_optimized_parameters.csv"), par_opt_found)
    save_jld2(joinpath(save_dir_actual, "$(this_save_key)_$(round(Int,sc))_optimization_history.jld2"), hist_found)
end
