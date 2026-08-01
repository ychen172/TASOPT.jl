"""
This script performs further optimization from multiple & adjacent cross warm-start model(s), but using parallel computing
Can create or remove optimization parameters. Will automatically extract initial guesses from warm start models.
multi-start means non-sequential optimization, case 1 result is not used as initial guess for case 2, because they are all started at the same time.
Use job array from orcd cluster to process multiple cases at the same time
Assume taskID starts from 1 instead of 0
"""

using TASOPT
include(__TASOPTindices__)
include(joinpath(__TASOPTroot__,"utils","sensitivity.jl"))
include(joinpath(__TASOPTroot__, "../example/utilities_for_optimization/objective_factory.jl"))
using .ObjectiveFactory: Constraint, Parameter, OptHistory, Requirement, optimize_singlePt_PFEI!, optimizer_wrapper_global_local, size_aircraft_w_param!,
                         save_vec_struct_csv, load_csv_constraints, load_csv_parameters, load_csv_requirements, save_jld2, load_jld2, OffDesMission
const min_fuse_radius = TASOPT.structures.find_minimum_radius_for_seats_per_row
const upd_fuse_pax! = TASOPT.structures.update_fuse_for_pax!
const success_statuses = ObjectiveFactory.success_statuses

#### Optimization parameters
# Prefix to read warm start data files
par_path_base_prefix = joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Jet_NoACT_V3_14_R1Size_/Opti_Jet_NoACT_V3_14_R1Size_") #Also the prefixed for aircraft model
# Path to save the models from optimization
save_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/")
save_key = "Opti_Jet_NoACT_V3_15P_R1Size"
# Optimization configuration parameters
flag_skip_global = true #Switch on to skip the global search if confident that warm start model can converge
max_iter_sizing = 150 # Maximum iterations for TASOPT sizing
optimizers = [:GN_CRS2_LM, :LN_NELDERMEAD] # Optimizer choice. [Global,Local]
max_num_iter_opt = [30000, 800, 10000] # Maximum number of optimization steps [Global, Local Coarse, Local Fine]
max_num_round_loc = [120, 30] #Maximum number of adaptive bounds refinements rounds for local search [Local Coarse, Local Fine]
span_glo_to_loc = 0.25 # span_local_search/span_global_search
# Setup optimization ranges
ranges_opti_nmi = collect(300.0:100.0:3000.0) #warm start should be available for these ranges
# Setup an entry (value does not matter) for test ratio to be sweep later in optimization loop
mis_opt = Requirement[]
push!(mis_opt, Requirement(:(parm[imRange,1]), 1e10)) #[m] different for each mission (overwritten)
push!(mis_opt, Requirement(:(parm[imWpay, 1]), 219964.5779*(1-1e-10))) #[N] Assume full load at R1
push!(mis_opt, Requirement(:(parg[igWpaymax]), 219964.5779)) #[N] Maximum payload weight (230 Pax with 956.36N per passenger)
push!(mis_opt, Requirement(:(options.ifuel), Int(24))) #Eth: 32       , Jet: 24            
push!(mis_opt, Requirement(:(parg[igrhofuel]), 817.0)) #Eth: 789.0    , Jet: 817.0 #kg/m3               
push!(mis_opt, Requirement(:(pare[iehvap, :, 1]), 358694.0)) #Eth: 918187.9 , Jet: 358694.0 #J/kg      
push!(mis_opt, Requirement(:(pare[iehvapcombustor, :, 1]), 358694.0)) #Eth: 918187.9 , Jet: 358694.0 #J/kg
push!(mis_opt, Requirement(:(options.has_ACT_fuel), false)) # Whether to allow additional center fuel tank (ACT)
push!(mis_opt, Requirement(:(options.compensate_ACT), false)) # Whether to increase the aircraft length to accmondate for the cargo space taken by ACT
push!(mis_opt, Requirement(:(fuse_tank.ACT_eta_vol), 1.00)) # Volumetric efficiency of ACT fuel tank
push!(mis_opt, Requirement(:(fuse_tank.ACT_eta_wei), 1.00)) # Gravimetric efficiency of ACT fuel tank
push!(mis_opt, Requirement(:(para[iaMach, ipclimbn:ipdescent1, 1]), 0.8)) # Cruise Mach number
push!(mis_opt, Requirement(:(vtail.opt_sizing), TailSizing.OEI)) # Vertical tail sizing mode (Engine out or fixed volume)
push!(mis_opt, Requirement(:(parg[igCLveout]), 0.47)) # Vertical tail CL at engine out condition
push!(mis_opt, Requirement(:(htail.opt_sizing), TailSizing.CLmaxFwdCG)) # Horizontal tail sizing mode (Forward CG or fixed volume)
push!(mis_opt, Requirement(:(htail.CL_max_fwd_CG), -0.7)) # Horizontal tail CL at forward CG condition
push!(mis_opt, Requirement(:(parg[igxeng]), 1e10)) # [m] Engine axial location (Though fixed value for each optimization, to be)
push!(mis_opt, Requirement(:(parg[igdxeng2wbox]), 1e10)) # [m] Same as above   (set by the default model before each optimization)
# Calibrated Fixed Design Parameters for Combustor
# push!(mis_opt, Requirement(:(pare[iepib, :, :]), 0.98)) #Combustor pressure ratio
# push!(mis_opt, Requirement(:(pare[ieepolf, :, :]), 0.92)) #Fan Poly Eff
# push!(mis_opt, Requirement(:(pare[ieepollc, :, :]), 0.92)) #LPC Poly Eff
# push!(mis_opt, Requirement(:(pare[ieepolhc, :, :]), 0.92)) #HPC Poly Eff
# push!(mis_opt, Requirement(:(pare[ieepolht, :, :]), 0.92)) #HPT Poly Eff
# push!(mis_opt, Requirement(:(pare[ieepollt, :, :]), 0.92)) #LPT Poly Eff
# push!(mis_opt, Requirement(:(pare[ieetab, :, :]), 0.9989)) #Combustion Eff
# push!(mis_opt, Requirement(:(parg[igTmetal]), 1280.0)) #Maximum metal temperature

# Constraints for this optimization
con_opt = Constraint[]
push!(con_opt, Constraint(:(wing.layout.span); pen_sca=1e4, lim_up=35.814)) # [m] Type C wing span constraint
push!(con_opt, Constraint(:(parm[imlBF, 1]); pen_sca=1e4, lim_up=2400.0)) # [m] Maximum balanced field length for takeoff
push!(con_opt, Constraint(:(para[iagamV, ipclimbn, 1]); pen_sca=1e4, lim_lo=0.015)) # [rad] Minimum cruise climb angle at TOC
push!(con_opt, Constraint(:(pare[ieTt3, :, 1]); pen_sca=1e4, lim_up=900.0)) # [K] Maximum compressor outlet temperature
push!(con_opt, Constraint(:(pare[ieTmet1, :, 1]); pen_sca=1e4, lim_up=1333.33)) # [K] Maximum turbine metal temperature
push!(con_opt, Constraint(:(parg[igdfan]); pen_sca=1e4, lim_up=2.0)) # [m] Maximum fan diameter
push!(con_opt, Constraint(:(parm[imWTO, 1]); pen_sca=1e4, lim_up=:(parg[igWMTO]), eps_buff=1e-4)) # [N] Maximum takeoff weight
push!(con_opt, Constraint(:(parm[imVfuel, 1]); pen_sca=1e4, lim_up=:(parg[igVfmax]), eps_buff=1e-4)) # [m3] Maximum fuel volume
# Parameters to be optimized (Initial guess(Only for local search), Upper bound, Lower bound, Initial Step Size(Only for local search))
bound_glob = Parameter[] #Step size = 1/5 of local search span = 1/5 * 1/4 * global search span specified below
push!(bound_glob, Parameter(:(wing.layout.sweep), 30.0, 60.0, 0.0, 5.0)) # [deg] Wing sweep angle 
push!(bound_glob, Parameter(:(wing.layout.AR), 10.0, 20.0, 5.0, 1.5)) # Wing aspect ratio
push!(bound_glob, Parameter(:(wing.inboard.cross_section.thickness_to_chord), 0.2, 0.6, 0.04, 0.05)) # Wing inner thickness to chord
push!(bound_glob, Parameter(:(wing.outboard.cross_section.thickness_to_chord), 0.2, 0.6, 0.04, 0.05))
push!(bound_glob, Parameter(:(wing.inboard.λ), 0.7, 1.0, 0.1, 0.1)) # Wing inner taper ratio
push!(bound_glob, Parameter(:(wing.outboard.λ), 0.3, 1.0, 0.1, 0.1))
push!(bound_glob, Parameter(:(parg[igyeng]), -1e10, -1e10, -1e10, 1e10)) # [m] Span wise engine location (Placeholder value here)
push!(bound_glob, Parameter(:(wing.layout.ηs), 0.3, 0.999, 0.1, 0.02)) # Wing span break location over half span
push!(bound_glob, Parameter(:(para[iarcls, ipclimb2:ipdescent4, 1]), 1.0, 2.0, 0.4, 0.15)) # Wing cl ratio span break to root at cruise
push!(bound_glob, Parameter(:(para[iarclt, ipclimb2:ipdescent4, 1]), 1.0, 2.0, 0.4, 0.15)) # Wing cl ratio of tip to root at cruise
push!(bound_glob, Parameter(:(para[iaCL, ipclimb2:ipdescent4, 1]), 0.6, 1.00, 0.3, 0.07)) # Wing CL at cruise
push!(bound_glob, Parameter(:(para[iaalt, ipcruise1, 1]), 10000.0, 20000.0, 4000.0, 1600.0)) # [m] Cruise altitude
push!(bound_glob, Parameter(:(pare[iepif, ipclimb2:ipdescent4, 1]), 2.0, 4.0, 1.25, 0.2)) # Fan PR at cruise 
push!(bound_glob, Parameter(:(pare[iepilc, ipclimb2:ipdescent4, 1]), 3.0, 4.0, 2.0, 0.4)) # LPC PR at cruise
push!(bound_glob, Parameter(:(pare[iepihc,ipclimb2:ipdescent4, 1]), 10.0, 50.0, 1.25, 1.0)) # HPC PR at cruise
push!(bound_glob, Parameter(:(pare[ieBPR, ipclimb2:ipdescent4, 1]), 8.0, 30.0, 1.0, 2.0)) # Fan BPR at cruise
push!(bound_glob, Parameter(:(pare[ieTt4, ipclimb2:ipdescent4, 1]), 1500.0, 2000.0, 1000.0, 100.0)) #[K] Turbine inlet temperature  at cruise
push!(bound_glob, Parameter(:(vtail.layout.AR), 2.0, 5.0, 1.0, 0.2)) #Vertical tail aspec ratio 

#### Setup job-array for ORCD cluster run
@assert length(ARGS) >= 2 "Usage: julia xxx.jl task_id num_tasks"
task_id = parse(Int,ARGS[1]) #Current task id
num_tasks = parse(Int,ARGS[2]) #Total number of tasks
@assert 1 <= task_id <= num_tasks

#### Setup the total save directory
save_dir_actual = joinpath(save_dir,save_key*"_") #Actual directory to save data
mkpath(save_dir_actual)

#### Optimization
for i in task_id:num_tasks:length(ranges_opti_nmi)
    #### Create inidividual log file for each task
    status_log = joinpath(save_dir_actual, "$(save_key)_$(round(Int,ranges_opti_nmi[i]))_OptLog.txt")
    open(status_log, "w") do io
        println(io, "Range_nmi,Status")
    end

    #### Make copies for each task
    mis_opt_cur = deepcopy(mis_opt)
    con_opt_cur = deepcopy(con_opt)
    bound_glob_cur = deepcopy(bound_glob)

    #### Load a sized default model regardless warm start or not
    if mod(i,2)==0 #Assume the total number of cases available for reading is even number
        ac = quickload_aircraft(par_path_base_prefix*"$(round(Int,ranges_opti_nmi[i-1])).jld2") #Even number iteration uses one case lower
    else
        ac = quickload_aircraft(par_path_base_prefix*"$(round(Int,ranges_opti_nmi[i+1])).jld2") #Odd number iteration uses one case higher
    end
    # ac = quickload_aircraft(par_path_base_prefix*"$(round(Int,ranges_opti_nmi[i])).jld2") #Multiple start
    # ac = quickload_aircraft(par_path_base_prefix) #Single start
    @assert ac.is_sized[1] #Make sure mission 1 is sized

    #### Update the fueselage radius and layout
    num_pass_row = Int(6) # Numebr of passengers per row
    fuse_radius = min_fuse_radius(num_pass_row, ac) #Minimum fueslage radius given the number of seats per row
    ac.fuselage.layout.cross_section.radius = fuse_radius
    ac.fuselage.cabin.front_seat_offset = 0.0  # [m] Front seat offset from the front of cylinder
    ac.fuselage.cabin.rear_seat_offset = 0.0 # [m] Back seat offset from the back of cylinder
    upd_fuse_pax!(ac) #Update carbin geometry
    
    #### Use the fuselage radius to update the mission and bound requirements
    mis_opt_cur[1].val = ranges_opti_nmi[i]*1852.0 #[m]
    mis_opt_cur[17].val = ac.wing.layout.box_x - fuse_radius
    mis_opt_cur[18].val = fuse_radius
    bound_glob_cur[7].val = fuse_radius * 3.0
    bound_glob_cur[7].bon_up = fuse_radius * 10.0
    bound_glob_cur[7].bon_lo = fuse_radius * 0.001 #unrealistic but should not have happened
    bound_glob_cur[7].d_val = fuse_radius * 0.2

    #### Setup special initial optimization parameters if choose to not run global search first but by using a warm start from a default model to directly do local search
    par_opt_cur = deepcopy(bound_glob_cur) #Use directly global bound if starting from global search
    if flag_skip_global
        for (j,para_cur) in enumerate(par_opt_cur)
            # Get starting value from the warm start aircraft model
            try # Many of the variable (Need to check why cruise value were set to many phases)
                para_cur.val = getNestedProp_fromExpr(ac; parm_sym=para_cur.name)[5] #Hopefully get the cruise phase value if possible
            catch
                para_cur.val = getNestedProp_fromExpr(ac; parm_sym=para_cur.name)[1]
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
        end
    end

    #### Save the mission requirement
    save_vec_struct_csv(joinpath(save_dir_actual, "$(save_key)_$(round(Int,ranges_opti_nmi[i]))_mission_requirements.csv"), mis_opt_cur)
    save_vec_struct_csv(joinpath(save_dir_actual, "$(save_key)_$(round(Int,ranges_opti_nmi[i]))_design_constraints.csv"), con_opt_cur)
    save_vec_struct_csv(joinpath(save_dir_actual, "$(save_key)_$(round(Int,ranges_opti_nmi[i]))_global_bounds.csv"), bound_glob_cur)

    #### Run the global local optimization process
    par_opt_found, status_found, hist_found = 
    optimizer_wrapper_global_local(ac, par_opt_cur, bound_glob_cur; miss_req=mis_opt_cur, constraints=con_opt_cur, max_iter_sizing=max_iter_sizing,
                                   optimizer_global=optimizers[1], max_iter_glo=max_num_iter_opt[1],   span_glo_to_loc=span_glo_to_loc, run_global=!flag_skip_global,
                                   optimizer_local =optimizers[2], max_iter_loc_C=max_num_iter_opt[2], max_round_loc_C=max_num_round_loc[1],
                                                                   max_iter_loc_F=max_num_iter_opt[3], max_round_loc_F=max_num_round_loc[2])
    
    #### Post-process the optimization results
    if status_found in success_statuses
        # attempt to update an aircraft model to save
        ac_copy = deepcopy(ac)
        sized_succeeded = size_aircraft_w_param!(ac_copy; mission_req=mis_opt_cur, parameters=par_opt_found, max_iter_sizing=max_iter_sizing)
        # correct the status
        status_found = sized_succeeded ? status_found : (:FAILED_TO_REPRODUCE)
        # Only save the new model if successfuly sized
        if sized_succeeded
            quicksave_aircraft(ac_copy, joinpath(save_dir_actual, "$(save_key)_$(round(Int,ranges_opti_nmi[i])).jld2"))
        end
    end
    
    #### save the results
    open(status_log, "a") do io
        println(io, "$(ranges_opti_nmi[i]),$(string(status_found))")
    end
    save_vec_struct_csv(joinpath(save_dir_actual, "$(save_key)_$(round(Int,ranges_opti_nmi[i]))_optimized_parameters.csv"), par_opt_found)
    save_jld2(joinpath(save_dir_actual, "$(save_key)_$(round(Int,ranges_opti_nmi[i]))_optimization_history.jld2"), hist_found)
end