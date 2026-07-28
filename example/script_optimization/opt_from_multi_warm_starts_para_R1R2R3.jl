"""
This script performs further optimization from multiple warm-start models
Can create or remove optimization parameters. Will automatically extract initial guesses from warm start models.
Extract the R1, R2, R3 missions from another model, and optimize for the PFEI of the R2 mission.
Also, optimize for the sizing mission range.
multi-start means non-sequential optimization, case 1 result is not used as initial guess for case 2.
Use job array from orcd cluster to process multiple cases at the same time
Assume taskID starts from 1 instead of 0
Run by: sbatch ../../julia-jobArr-orcd.sh
"""

using Glob
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
par_path_base_prefix = joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Eth_NoACT_V2_1_ULByEng_Aft3_/Opti_Eth_NoACT_V2_1_ULByEng_Aft3_") #Also the prefixed for aircraft model
# Path to save the models from optimization
save_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/")
save_key = "Opti_Eth_NoACT_V2_2_MatchR1R2R3"
# Mission extraction directory
miss_key = "Opti_Jet_NoACT_V2_1_ULByEng_Eth_PRD_"
miss_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/"*miss_key*"/")
# Optimization configuration parameters
flag_skip_global = true #Switch on to skip the global search if confident that warm start model can converge
max_iter_sizing = 150 # Maximum iterations for TASOPT sizing
optimizers = [:GN_CRS2_LM, :LN_NELDERMEAD] # Optimizer choice. [Global,Local]
max_num_iter_opt = [30000, 800, 10000] # Maximum number of optimization steps [Global, Local Coarse, Local Fine]
max_num_round_loc = [120, 30] #Maximum number of adaptive bounds refinements rounds for local search [Local Coarse, Local Fine]
span_glo_to_loc = 0.25 # span_local_search/span_global_search
# Setup optimization ranges
ranges_opti_nmi = [300]#collect(300.0:100.0:3000.0) #warm start reservoir to be selected from
# Setup an entry (value does not matter) for test ratio to be sweep later in optimization loop
mis_opt = Requirement[]
push!(mis_opt, Requirement(:(options.ifuel), Int(32)))                #Eth: 32       , Jet: 24
push!(mis_opt, Requirement(:(parg[igrhofuel]), 789.0))                #Eth: 789.0    , Jet: 817.0 #kg/m3
push!(mis_opt, Requirement(:(pare[iehvap, :, 1]), 918187.9))          #Eth: 918187.9 , Jet: 358694.0 #J/kg
push!(mis_opt, Requirement(:(pare[iehvapcombustor, :, 1]), 918187.9)) #Eth: 918187.9 , Jet: 358694.0 #J/kg
push!(mis_opt, Requirement(:(options.has_ACT_fuel), false)) # Whether to allow additional center fuel tank (ACT)
push!(mis_opt, Requirement(:(options.compensate_ACT), false)) # Whether to increase the aircraft length to accmondate for the cargo space taken by ACT
push!(mis_opt, Requirement(:(fuse_tank.ACT_eta_vol), 1.00)) # Volumetric efficiency of ACT fuel tank
push!(mis_opt, Requirement(:(fuse_tank.ACT_eta_wei), 1.00)) # Gravimetric efficiency of ACT fuel tank
push!(mis_opt, Requirement(:(para[iaMach, ipclimbn:ipdescent1, 1]), 0.8)) # Cruise Mach number
push!(mis_opt, Requirement(:(pare[iepilc, ipclimb2:ipdescent4, 1]), 3.0)) # Engine LPC pressure ratio at cruise
push!(mis_opt, Requirement(:(vtail.opt_sizing), TailSizing.OEI)) # Vertical tail sizing mode (Engine out or fixed volume)
push!(mis_opt, Requirement(:(parg[igCLveout]), 0.47)) # Vertical tail CL at engine out condition
push!(mis_opt, Requirement(:(htail.opt_sizing), TailSizing.CLmaxFwdCG)) # Horizontal tail sizing mode (Forward CG or fixed volume)
push!(mis_opt, Requirement(:(htail.CL_max_fwd_CG), -0.7)) # Horizontal tail CL at forward CG condition
push!(mis_opt, Requirement(:(parg[igxeng]), 1e10)) # [m] Engine axial location (Though fixed value for each optimization, to be)
push!(mis_opt, Requirement(:(parg[igdxeng2wbox]), 1e10)) # [m] Same as above   (set by the default model before each optimization)
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
push!(bound_glob, Parameter(:(pare[iepihc,ipclimb2:ipdescent4, 1]), 10.0, 50.0, 1.25, 1.0)) # HPC PR at cruise
push!(bound_glob, Parameter(:(pare[ieBPR, ipclimb2:ipdescent4, 1]), 8.0, 30.0, 1.0, 2.0)) # Fan BPR at cruise
push!(bound_glob, Parameter(:(pare[ieTt4, ipclimb2:ipdescent4, 1]), 1500.0, 2000.0, 1000.0, 100.0)) #[K] Turbine inlet temperature  at cruise
push!(bound_glob, Parameter(:(vtail.layout.AR), 2.0, 5.0, 1.0, 0.2)) #Vertical tail aspec ratio 
push!(bound_glob, Parameter(:(parm[imRange,1]), 1500.0*1852.0, 4000.0*1852.0, 300.0*1852.0, 100.0*1852.0)) #[m] To be optimized from the given R2 range
push!(bound_glob, Parameter(:(parm[imWpay, 1]), -1e10, -1e10, -1e10, 1e10)) #[N] To be optimized based on the aircraft's fixed maximum capacity

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

    #### Read in mission file
    miss_dir_cur = joinpath(miss_dir,miss_key*"$(round(Int,ranges_opti_nmi[i]))_")
    R1_key = miss_key*"$(round(Int,ranges_opti_nmi[i]))_R1_*.jld2"
    R1_key = only(glob(R1_key, miss_dir_cur)) #Find the only matched file
    ac_miss = quickload_aircraft(R1_key)
    Range_R1 = ac_miss.parm[imRange, 2] #[m]
    Wpay_R1 = ac_miss.parm[imWpay, 2] #[N]
    R2_key = miss_key*"$(round(Int,ranges_opti_nmi[i]))_R2_*.jld2"
    R2_key = only(glob(R2_key, miss_dir_cur))
    ac_miss = quickload_aircraft(R2_key)
    Range_R2 = ac_miss.parm[imRange, 2] #[m]
    Wpay_R2 = ac_miss.parm[imWpay, 2] #[N]
    R3_key = miss_key*"$(round(Int,ranges_opti_nmi[i]))_R3_*.jld2"
    R3_key = only(glob(R3_key, miss_dir_cur)) #Already contain the directory path
    ac_miss = quickload_aircraft(R3_key)
    Range_R3 = ac_miss.parm[imRange, 2] #[m]
    Wpay_R3 = ac_miss.parm[imWpay, 2] #[N]

    #### Load a sized default model regardless warm start or not
    i_cloest = argmin(abs.(ranges_opti_nmi .- Range_R2/1852.0)) #Find the cloest range to read the warm start model
    ac = quickload_aircraft(par_path_base_prefix*"$(round(Int,ranges_opti_nmi[i_cloest])).jld2")
    @assert ac.is_sized[1] #Make sure mission 1 is sized
    
    #### Update the fueselage radius and layout
    num_pass_row = Int(6) # Numebr of passengers per row
    fuse_radius = min_fuse_radius(num_pass_row, ac) #Minimum fueslage radius given the number of seats per row
    ac.fuselage.layout.cross_section.radius = fuse_radius
    ac.fuselage.cabin.front_seat_offset = 0.0  # [m] Front seat offset from the front of cylinder
    ac.fuselage.cabin.rear_seat_offset = 0.0 # [m] Back seat offset from the back of cylinder
    upd_fuse_pax!(ac) #Update carbin geometry
    
    #### Use the fuselage radius to update the mission and bound requirements
    mis_opt_cur[15].val = ac.wing.layout.box_x - fuse_radius
    mis_opt_cur[16].val = fuse_radius
    #
    bound_glob_cur[7].val = fuse_radius * 3.0
    bound_glob_cur[7].bon_up = fuse_radius * 10.0
    bound_glob_cur[7].bon_lo = fuse_radius * 0.001 #unrealistic but should not have happened
    bound_glob_cur[7].d_val = fuse_radius * 0.2
    #
    bound_glob_cur[19].val = ac.parg[igWpaymax]/2.0
    bound_glob_cur[19].bon_up = ac.parg[igWpaymax]
    bound_glob_cur[19].bon_lo = ac.parm[imWperpax,1]
    bound_glob_cur[19].d_val = ac.parg[igWpaymax]*0.2

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
            if j == 18
                para_cur.val = Range_R2
            elseif j == 19
                para_cur.val = Wpay_R2
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

    #### Prepare off-design mission for optimize at R2
    range_off_des_nmi = [Range_R1, Range_R2, Range_R3]./1852.0 #[nmi]
    wei_pay_off_des_N = [Wpay_R1, Wpay_R2, Wpay_R3] #[N]
    idx_fuel_off_des = fill(Int(32), 3)
    rho_fuel_off_des_kgm3 = fill(789.0, 3)
    hvap_fuel_off_des_Jkg = fill(918187.9, 3)
    off_des_miss = OffDesMission{Float64}(range_off_des_nmi, wei_pay_off_des_N, idx_fuel_off_des, rho_fuel_off_des_kgm3, hvap_fuel_off_des_Jkg)
    off_des_constraints = [:WPay,:MWTO,:VolFuel]
    PFEI_Weighting = [0.0,0.0,1.0,0.0] # Here we only use R2 PFEI for optimization

    #### Run the global local optimization process
    par_opt_found, status_found, hist_found = 
    optimizer_wrapper_global_local(ac, par_opt_cur, bound_glob_cur; miss_req=mis_opt_cur, constraints=con_opt_cur, max_iter_sizing=max_iter_sizing,
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