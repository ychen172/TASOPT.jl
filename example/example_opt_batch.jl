using TASOPT, NLopt
using JLD2
include(joinpath(@__DIR__, "optimize_rangefuel.jl"))
using .OptimizeRangeFuel: MissionReq, BoundsOpt, ConstraintsOpt, optimize_rangefuel_fun!
include(__TASOPTindices__)
include(joinpath(@__DIR__, "postprocess.jl"))
using .PostProcess: save_hist_compact!

#### Start from an initial aircraft model
# Load from
save_dir  = "ModelSaved"
load_name = "acOptimized_Bat" #jld2
# Save to
save_name = "acOptimized_BatOptJet" #jld2

#### Setup optimization parameters
iters_max_opt = 1000 #1000 #Number of interations

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
bounds_opt_global.Tt4_lim          = (1000.0, 2000.0,  200.0) # Change to 50
bounds_opt_global.PR_hpc_lim       = (1.25,   50.0,    5.0)   # Change to 0.1
bounds_opt_global.PR_fan_lim       = (1.25,   4.0,     0.2)
bounds_opt_global.PR_lpc_lim       = (2.999,  3.001,   0.0)   # Change to 0.00001
bounds_opt_global.BPR_lim          = (1.0,    30.0,    3.0)
# Create a local boundary
bounds_opt_local                   = BoundsOpt()
bounds_opt_local.Tt4_lim = (bounds_opt_local.Tt4_lim[1], bounds_opt_local.Tt4_lim[2], 50.0)
bounds_opt_local.PR_hpc_lim = (bounds_opt_local.PR_hpc_lim[1], bounds_opt_local.PR_hpc_lim[2], 0.2)
bounds_opt_local.PR_lpc_lim = (bounds_opt_local.PR_lpc_lim[1], bounds_opt_local.PR_lpc_lim[2], 0.00001)

####Initialize the log
status_log = joinpath(save_dir, "$(save_name)_Log.txt")
open(status_log, "w") do io
    println(io, "range_nmi,status")
end

####Optimization function
function main()
    ####Load the starting model
    ac = quickload_aircraft(joinpath(save_dir,"$(load_name).jld2"))
    failsafe_name = load_name #something bad happen in the 1st step, reload this initial aircraft model for the 2nd step
    ####Optimization
    for (idx, ran_cur) in enumerate(range_lst)
        ####Update the mission requirement
        mission_req.range_des = (ran_cur * 1852.0)  #Design flight range (m)

        #### Run the optimization
        status_cur, hist_optim_cur = optimize_rangefuel_fun!(ac; mission_req=mission_req, bounds_opt=bounds_opt_local, constraints_opt=constraints_opt, iters_max_opt=iters_max_opt)

        #### Judging and save the result
        if status_cur in (:SUCCESS, :STOPVAL_REACHED, :FTOL_REACHED, :XTOL_REACHED, :MAXEVAL_REACHED, :MAXTIME_REACHED)
            quicksave_aircraft(ac, joinpath(save_dir, "$(save_name*string(round(Int,ran_cur))).jld2"))
            failsafe_name = save_name*string(round(Int,ran_cur))
        else
            ac = quickload_aircraft(joinpath(save_dir,"$(failsafe_name).jld2"))
        end

        ####Store the log
        open(status_log, "a") do io
            println(io, "$(round(Int,ran_cur)),$(string(status_cur))")
        end

        ####Store the optimization history
        save_hist_compact!(joinpath(save_dir, "$(save_name*string(round(Int,ran_cur)))_History.jld2"), ran_cur, hist_optim_cur)
    end
end
main()