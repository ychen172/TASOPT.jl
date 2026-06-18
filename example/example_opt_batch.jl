using TASOPT
include(__TASOPTindices__)
include(joinpath(__TASOPTroot__, "../example/","utilities_for_optimization/","objective_factory.jl"))
using .ObjectiveFactory: Constraint, Parameter, OptHistory, Requirement, 
                         optimize_singlePt_PFEI!, optimizer_wrapper_global_local, 
                         save_vec_struct_csv, load_vec_struct_csv, save_jld2, load_jld2
const find_minimum_radius_for_seats_per_row = TASOPT.structures.find_minimum_radius_for_seats_per_row
const update_fuse_for_pax! = TASOPT.structures.update_fuse_for_pax!

#### Setup IO
# Initial data path
read_dir_ini = joinpath(__TASOPTroot__,"../example/","ModelSaved/","JetFuel_NoACT_Opt/","JetFuel_NoACT_Opt3000.jld2")
save_dir  = joinpath(__TASOPTroot__,"../example/","ModelSaved/")
# Save to
save_key = "Opt_Jet" #jld2

#### Setup mission requirement
# Some special parameters to fix
num_pass_per_row = Int(6)
front_seat_offset = 0.0 #3.048
rear_seat_offset = front_seat_offset/2.0
# Load and setup the aircraft (There are some overwriting cannot be fully automated by the solver)
ac = quickload_aircraft(read_dir_ini)
@assert ac.is_sized[1]
# Change the mode to size tails
ac.vtail.opt_sizing = TailSizing.OEI
ac.htail.opt_sizing = TailSizing.CLmaxFwdCG
# Calculate optimal fuselage layout
radius_fuselage = find_minimum_radius_for_seats_per_row(num_pass_per_row, ac) #Compute the set the current best cross-section
ac.fuselage.layout.cross_section.radius = radius_fuselage
ac.fuselage.cabin.front_seat_offset = front_seat_offset
ac.fuselage.cabin.rear_seat_offset = rear_seat_offset
update_fuse_for_pax!(ac) #Activate carbin sizing

# Additional parameter can go with the mission parameter
mis_opt = Requirement[]
push!(mis_opt, Requirement(:(parm[imRange,1]), 300*1852.0)) #[m] Being updated in loop later on collect(300:100:3000)
push!(mis_opt, Requirement(:(options.ifuel), Int(24))) #Fuel Index: Jet Fuel(24), Ethanol(32)
push!(mis_opt, Requirement(:(parg[igrhofuel]), 817.0)) #Fuel Density: Jet Fuel(817.0) (kg/m3), Ethanol(789.0) (kg/m3)
push!(mis_opt, Requirement(:(pare[iehvap, :, 1]), 358694.0)) #Heat of Vaporization: Jet Fuel(358694.0) (J/kg), Ethanol(918187.9) (J/kg)
push!(mis_opt, Requirement(:(pare[iehvapcombustor, :, 1]), 358694.0))
push!(mis_opt, Requirement(:(pare[iepilc, ipclimb2:ipdescent4, 1]), 3.0))
push!(mis_opt, Requirement(:(para[iaMach, ipclimbn:ipdescent1, 1]), 0.8))
push!(mis_opt, Requirement(:(para[iarcls, ipclimb2:ipdescent4, 1]), 1.08))
push!(mis_opt, Requirement(:(parg[igxeng]), ac.wing.layout.box_x - radius_fuselage))
push!(mis_opt, Requirement(:(parg[igCLveout]), 0.5))
push!(mis_opt, Requirement(:(htail.CL_max_fwd_CG), -0.7))

#### Setup constraints
con_opt = Constraint[]
push!(con_opt, Constraint(:(wing.layout.span); pen_sca=1e4, lim_up=35.814))
push!(con_opt, Constraint(:(parm[imlBF, 1]); pen_sca=1e4, lim_up=2400.0))
push!(con_opt, Constraint(:(para[iagamV, ipclimbn, 1]); pen_sca=1e4, lim_lo=0.015))
push!(con_opt, Constraint(:(pare[ieTt3, :, 1]); pen_sca=1e4, lim_up=900.0))
push!(con_opt, Constraint(:(pare[ieTmet1, :, 1]); pen_sca=1e4, lim_up=1333.33))
push!(con_opt, Constraint(:(parg[igdfan]); pen_sca=1e4, lim_up=2.0))
push!(con_opt, Constraint(:(parm[imWTO, 1]); pen_sca=1e4, lim_up=:(parg[igWMTO]), eps_buff=1e-4))
push!(con_opt, Constraint(:(parm[imVfuel, 1]); pen_sca=1e4, lim_up=:(parg[igVfmax]), eps_buff=1e-4))

#### Setup optimization parameters
par_opt = Parameter[]
push!(par_opt, Parameter(:(wing.layout.sweep), 30.0, 60.0, 0.0, 5.0))
push!(par_opt, Parameter(:(wing.layout.AR), 10.0, 20.0, 5.0, 1.5))
push!(par_opt, Parameter(:(wing.inboard.cross_section.thickness_to_chord), 0.2, 0.6, 0.04, 0.05))
push!(par_opt, Parameter(:(wing.outboard.cross_section.thickness_to_chord), 0.2, 0.6, 0.04, 0.05))
push!(par_opt, Parameter(:(wing.inboard.λ), 0.7, 1.0, 0.1, 0.1))
push!(par_opt, Parameter(:(wing.outboard.λ), 0.3, 1.0, 0.1, 0.1))
push!(par_opt, Parameter(:(parg[igyeng]), radius_fuselage*3.0, radius_fuselage*4.0, radius_fuselage*2.0, radius_fuselage*0.2))
push!(par_opt, Parameter(:(wing.layout.ηs), 0.3, 0.4, 0.2, 0.02))
push!(par_opt, Parameter(:(para[iarclt, ipclimb2:ipdescent4, 1]), 1.0, 2.0, 0.4, 0.15))

push!(par_opt, Parameter(:(para[iaCL, ipclimb2:ipdescent4, 1]), 0.6, 1.00, 0.3, 0.07))
push!(par_opt, Parameter(:(para[iaalt, ipcruise1, 1]), 10000.0, 20000.0, 4000.0, 1600.0))
push!(par_opt, Parameter(:(pare[iepif, ipclimb2:ipdescent4, 1]), 2.0, 4.0, 1.25, 0.2))
push!(par_opt, Parameter(:(pare[iepihc,ipclimb2:ipdescent4, 1]), 10.0, 50.0, 1.25, 1.0))
push!(par_opt, Parameter(:(pare[ieBPR, ipclimb2:ipdescent4, 1]), 8.0, 30.0, 1.0, 2))
push!(par_opt, Parameter(:(pare[ieTt4, ipclimb2:ipdescent4, 1]), 1500.0, 2000.0, 1000.0, 100.0))

push!(par_opt, Parameter(:(vtail.layout.AR), 2.0, 5.0, 1.0, 0.2))

#### Test single point optimization
mis_opt[1].val = 300*1852.0 #[m]
_,hist = optimize_singlePt_PFEI!(ac, par_opt; mission_req=mis_opt, constraints=con_opt, max_iter_optim=50000000, pen_failed_sizing=100.0, optimizer_type=:GN_CRS2_LM)

# #### Setup constraints
# constraints_opt = ConstraintsOpt()
# save_struct(constraints_opt, joinpath(save_dir, "$(save_name)_GlobalConstraints.csv"))

# #### Setup search ranges for optimized parameters
# bounds_opt_global                  = BoundsOpt() #This hard limits should be fixed across all cases
# bounds_opt_global.AR_lim           = (5.0,    20.0,    3.0) # Use the third entry to set the minimum half-local-search-range 
# bounds_opt_global.CL_lim           = (0.30,   1.00,    0.14)
# bounds_opt_global.sweep_deg_lim    = (0.0,    60.0,    5.0)
# bounds_opt_global.alt_cruise_m_lim = (4000.0, 20000.0, 2000.0)
# bounds_opt_global.taper_in_lim     = (0.025,  1.00,    0.2)
# bounds_opt_global.taper_out_lim    = (0.025,  1.00,    0.2)
# bounds_opt_global.tc_root_lim      = (0.04,   0.6,     0.1)
# bounds_opt_global.tc_span_lim      = (0.04,   0.6,     0.1)
# bounds_opt_global.rcls_lim         = (0.40,   2.0,     0.3)
# bounds_opt_global.rclt_lim         = (0.40,   2.0,     0.3)
# bounds_opt_global.Tt4_lim          = (1000.0, 2000.0,  200.0) 
# bounds_opt_global.PR_hpc_lim       = (1.25,   50.0,    5.0)
# bounds_opt_global.PR_fan_lim       = (1.25,   4.0,     0.2)
# bounds_opt_global.PR_lpc_lim       = (2.999,  3.001,   0.0)
# bounds_opt_global.BPR_lim          = (1.0,    30.0,    3.0)
# save_struct(bounds_opt_global, joinpath(save_dir, "$(save_name)_GlobalBounds.csv"))
# # Create a local boundary
# bounds_opt_local            = BoundsOpt()
# bounds_opt_local.Tt4_lim    = (bounds_opt_local.Tt4_lim[1], bounds_opt_local.Tt4_lim[2], 50.0) #Change step size to 50 K
# bounds_opt_local.PR_hpc_lim = (bounds_opt_local.PR_hpc_lim[1], bounds_opt_local.PR_hpc_lim[2], 0.2) # Change step size to 0.1 pressure ratio
# bounds_opt_local.PR_lpc_lim = (bounds_opt_global.PR_lpc_lim[1], bounds_opt_global.PR_lpc_lim[2], 0.00001) #Fixed Param: Consistent with global but fake search step
# # Clip the current local bounds by global bounds before the optimization
# # Given the large uncertaintin in the intial guess: Suppress the intial clipping to get a better starting solution. clip_loc_bound!(bounds_opt_local, bounds_opt_global)

# #### Initialize the log
# status_log = joinpath(save_dir, "$(save_name)_Log.txt")
# open(status_log, "w") do io
#     println(io, "range_nmi,status")
# end

# #### Optimization
# function main()
#     #### Load the starting model
#     ac = quickload_aircraft(joinpath(save_dir,"$(load_name).jld2"))
#     failsafe_name = load_name #something bad happen in the 1st step, reload this initial aircraft model for the 2nd step
#     bounds_local = deepcopy(bounds_opt_local)
#     failsafe_bounds_local = deepcopy(bounds_local) #Incase failsafe to old ac, also failsafe to old bound

#     #### Intialize the iterations
#     numRanges = length(range_lst)
#     idxRan = 1
#     flgReRun = false

#     #### Outer loop for the range
#     while idxRan <= numRanges || flgReRun
#         #### Check for re-run
#         if flgReRun
#             println("ReRun the previous range")
#             idxRan -= 1
#         end
#         ran_cur = range_lst[idxRan]
#         println("Runing range $(ran_cur)")

#         #### Update the range requirement
#         mission_req.range_des = (ran_cur * 1852.0)  #Design flight range (m)

#         #### Phase one optimization
#         # Initalization
#         countWhile = 0
#         flag_bound_change = true
#         while (flag_bound_change)
#             countWhile += 1
#             println("Range: $(ran_cur) on coarse run #$(countWhile)")
#             # Coarse optimization
#             status_cur, hist_optim_cur = optimize_rangefuel_fun!(ac; mission_req=mission_req, bounds_opt=bounds_local, constraints_opt=constraints_opt, iters_max_opt=iters_max_opt_coarse)

#             # Judge and adjust the bound
#             if (countWhile<1000) && (status_cur in (:SUCCESS, :STOPVAL_REACHED, :FTOL_REACHED, :XTOL_REACHED, :MAXEVAL_REACHED, :MAXTIME_REACHED))
#                 # Update the local bounds
#                 flag_bound_change = adjust_bounds!(ac, bounds_local, bounds_opt_global)
#                 if flag_bound_change
#                     clip_loc_bound!(bounds_local, bounds_opt_global)
#                 end
#             else
#                 ac = quickload_aircraft(joinpath(save_dir,"$(failsafe_name).jld2"))
#                 bounds_local = deepcopy(failsafe_bounds_local)
#                 status_cur = :Search_Range_Identification_Failed
#                 break
#             end
#         end

#         #### Phase two optimization
#         flgReRun = false # Should the current range be rerun
#         if (status_cur != :Search_Range_Identification_Failed)
#             println("Range: $(ran_cur) on fine run")
#             # Fine optimization
#             status_cur, hist_optim_cur = optimize_rangefuel_fun!(ac; mission_req=mission_req, bounds_opt=bounds_local, constraints_opt=constraints_opt, iters_max_opt=iters_max_opt_fine)

#             # Judge for bound adjustment
#             if status_cur in (:SUCCESS, :STOPVAL_REACHED, :FTOL_REACHED, :XTOL_REACHED, :MAXEVAL_REACHED, :MAXTIME_REACHED)
#                 # Check on bound
#                 flag_bound_change = adjust_bounds!(ac, bounds_local, bounds_opt_global)
#                 if flag_bound_change #If the flag bound change is true, need to re-run the current range(through phase 1 and phase 2 again)
#                     clip_loc_bound!(bounds_local, bounds_opt_global)
#                 end
#                 if !flag_bound_change #Nice, fine run succeed with perfect search range too, save the result
#                     # Save the optimized model
#                     quicksave_aircraft(ac, joinpath(save_dir, "$(save_name*string(round(Int,ran_cur))).jld2"))
#                     save_struct(bounds_local, joinpath(save_dir, "$(save_name*string(round(Int,ran_cur)))_BoundLocal.csv"))
#                     # Update the fail save
#                     failsafe_name = save_name*string(round(Int,ran_cur)) #Update the failsafe too
#                     failsafe_bounds_local = deepcopy(bounds_local)
#                     # Save the optimization history
#                     save_hist_compact!(joinpath(save_dir, "$(save_name*string(round(Int,ran_cur)))_History.jld2"), ran_cur, hist_optim_cur)
#                 else
#                     flgReRun = true
#                 end
#             else #Coarse search fine but fine search failed, fall back to old solution and skip the current iteration, keep the fine search error status
#                 ac = quickload_aircraft(joinpath(save_dir,"$(failsafe_name).jld2"))
#                 bounds_local = deepcopy(failsafe_bounds_local)
#             end
#         end

#         #### Store the log (Only if the current (phase 1 + phase 2 operations) are not being rerun)
#         if !flgReRun
#             open(status_log, "a") do io
#                 println(io, "$(round(Int,ran_cur)),$(string(status_cur))")
#             end
#         end

#         #### Update Index
#         idxRan += 1
#     end
# end
# main()