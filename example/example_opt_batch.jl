using TASOPT, NLopt
using JLD2
include(joinpath(@__DIR__, "optimize_rangefuel.jl"))
using .OptimizeRangeFuel: MissionReq, BoundsOpt, ConstraintsOpt, optimize_rangefuel_fun!, clip_loc_bound!, adjust_bounds!
include(__TASOPTindices__)
include(joinpath(@__DIR__, "postprocess.jl"))
using .PostProcess: save_hist_compact!, save_struct

#### Start from an initial aircraft model
# Load from
save_dir  = "ModelSaved"
load_name = "acOptimized_Bat" #jld2
# Save to
save_name = "acOptimized_BatOptJet" #jld2

#### Setup optimization parameters
iters_max_opt_coarse = 500 #Initial optimizations steps to find a good search bound
iters_max_opt_fine = 100000 #After the good search bound is estabilished

#### Setup the mission requirement
mission_req = MissionReq()
mission_req.idx_fuel = 24 #Fuel Index: Jet Fuel(24), Ethanol(32)
mission_req.rho_fuel = 817.0 #Fuel Density: Jet Fuel(817.0) (kg/m3), Ethanol(789.0) (kg/m3)
mission_req.hvap_fuel = 358694.0 #Heat of Vaporization: Jet Fuel(358694.0) (J/kg), Ethanol(918187.9) (J/kg)
# Range list to iterate through
range_lst = [301] #collect(300:100:3000) #[3000,2900] #Range in nmi (Need to be turned to meter for input) #collect(3000.0:-100:300) #28 cases

#### Setup constraints
constraints_opt = ConstraintsOpt()

#### Setup search ranges for optimized parameters
bounds_opt_global                  = BoundsOpt() #This hard limits should be fixed across all cases
bounds_opt_global.AR_lim           = (5.0,    20.0,    3.0) # Use the third entry to set the minimum half-local-search-range 
bounds_opt_global.CL_lim           = (0.30,   1.00,    0.14)
bounds_opt_global.sweep_deg_lim    = (0.0,    40.0,    5.0)
bounds_opt_global.alt_cruise_m_lim = (4000.0, 20000.0, 2000.0)
bounds_opt_global.taper_in_lim     = (0.025,  1.00,    0.2)
bounds_opt_global.taper_out_lim    = (0.025,  1.00,    0.2)
bounds_opt_global.tc_root_lim      = (0.04,   0.6,     0.1)
bounds_opt_global.tc_span_lim      = (0.04,   0.6,     0.1)
bounds_opt_global.rcls_lim         = (0.40,   2.0,     0.3)
bounds_opt_global.rclt_lim         = (0.40,   2.0,     0.3)
bounds_opt_global.Tt4_lim          = (1000.0, 2000.0,  200.0) 
bounds_opt_global.PR_hpc_lim       = (1.25,   50.0,    5.0)
bounds_opt_global.PR_fan_lim       = (1.25,   4.0,     0.2)
bounds_opt_global.PR_lpc_lim       = (2.999,  3.001,   0.0)
bounds_opt_global.BPR_lim          = (1.0,    30.0,    3.0)
save_struct(bounds_opt_global, joinpath(save_dir, "$(save_name)_GlobalBounds.csv"))
# Create a local boundary
bounds_opt_local            = BoundsOpt()
bounds_opt_local.Tt4_lim    = (bounds_opt_local.Tt4_lim[1], bounds_opt_local.Tt4_lim[2], 50.0) #Change step size to 50 K
bounds_opt_local.PR_hpc_lim = (bounds_opt_local.PR_hpc_lim[1], bounds_opt_local.PR_hpc_lim[2], 0.2) # Change step size to 0.1 pressure ratio
bounds_opt_local.PR_lpc_lim = (bounds_opt_global.PR_lpc_lim[1], bounds_opt_global.PR_lpc_lim[2], 0.00001) #Fixed Param: Consistent with global but fake search step
# Clip the current local bounds by global bounds before the optimization
clip_loc_bound!(bounds_opt_local, bounds_opt_global)

#### Initialize the log
status_log = joinpath(save_dir, "$(save_name)_Log.txt")
open(status_log, "w") do io
    println(io, "range_nmi,status")
end

#### Optimization
function main()
    #### Load the starting model
    ac = quickload_aircraft(joinpath(save_dir,"$(load_name).jld2"))
    failsafe_name = load_name #something bad happen in the 1st step, reload this initial aircraft model for the 2nd step
    bounds_local = deepcopy(bounds_opt_local)
    failsafe_bounds_local = deepcopy(bounds_local) #Incase failsafe to old ac, also failsafe to old bound

    #### Intialize the iterations
    numRanges = length(range_lst)
    idxRan = 1
    flgReRun = false

    #### Outer loop for the range
    while idxRan <= numRanges || flgReRun
        #### Check for re-run
        if flgReRun
            idxRan -= 1
        end
        ran_cur = range_lst[idxRan]

        #### Update the range requirement
        mission_req.range_des = (ran_cur * 1852.0)  #Design flight range (m)

        #### Phase one optimization
        # Initalization
        countWhile = 0
        flag_bound_change = true
        while (flag_bound_change)
            countWhile += 1
            # Coarse optimization
            status_cur, hist_optim_cur = optimize_rangefuel_fun!(ac; mission_req=mission_req, bounds_opt=bounds_local, constraints_opt=constraints_opt, iters_max_opt=iters_max_opt_coarse)

            # Judge and adjust the bound
            if (countWhile<1000) && (status_cur in (:SUCCESS, :STOPVAL_REACHED, :FTOL_REACHED, :XTOL_REACHED, :MAXEVAL_REACHED, :MAXTIME_REACHED))
                # Update the local bounds
                flag_bound_change = adjust_bounds!(ac, bounds_local, bounds_opt_global)
                if flag_bound_change
                    clip_loc_bound!(bounds_local, bounds_opt_global)
                end
            else
                ac = quickload_aircraft(joinpath(save_dir,"$(failsafe_name).jld2"))
                bounds_local = deepcopy(failsafe_bounds_local)
                status_cur = :Search_Range_Identification_Failed
                break
            end
        end

        #### Phase two optimization
        flgReRun = false # Should the current range be rerun
        if (status_cur != :Search_Range_Identification_Failed)
            # Fine optimization
            status_cur, hist_optim_cur = optimize_rangefuel_fun!(ac; mission_req=mission_req, bounds_opt=bounds_local, constraints_opt=constraints_opt, iters_max_opt=iters_max_opt_fine)

            # Judge for bound adjustment
            if status_cur in (:SUCCESS, :STOPVAL_REACHED, :FTOL_REACHED, :XTOL_REACHED, :MAXEVAL_REACHED, :MAXTIME_REACHED)
                # Check on bound
                flag_bound_change = adjust_bounds!(ac, bounds_local, bounds_opt_global)
                if flag_bound_change #If the flag bound change is true, need to re-run the current range(through phase 1 and phase 2 again)
                    clip_loc_bound!(bounds_local, bounds_opt_global)
                end
                if !flag_bound_change #Nice, fine run succeed with perfect search range too, save the result
                    # Save the optimized model
                    quicksave_aircraft(ac, joinpath(save_dir, "$(save_name*string(round(Int,ran_cur))).jld2"))
                    save_struct(bounds_local, joinpath(save_dir, "$(save_name*string(round(Int,ran_cur)))_BoundLocal.csv"))
                    # Update the fail save
                    failsafe_name = save_name*string(round(Int,ran_cur)) #Update the failsafe too
                    failsafe_bounds_local = deepcopy(bounds_local)
                    # Save the optimization history
                    save_hist_compact!(joinpath(save_dir, "$(save_name*string(round(Int,ran_cur)))_History.jld2"), ran_cur, hist_optim_cur)
                else
                    flgReRun = true
                end
            else #Coarse search fine but fine search failed, fall back to old solution and skip the current iteration, keep the fine search error status
                ac = quickload_aircraft(joinpath(save_dir,"$(failsafe_name).jld2"))
                bounds_local = deepcopy(failsafe_bounds_local)
            end
        end

        #### Store the log (Only if the current (phase 1 + phase 2 operations) are not being rerun)
        if !flgReRun
            open(status_log, "a") do io
                println(io, "$(round(Int,ran_cur)),$(string(status_cur))")
            end
        end

        #### Update Index
        idxRan += 1
    end
end
main()