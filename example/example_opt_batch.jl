using TASOPT, NLopt
using JLD2
include(joinpath(@__DIR__, "optimize_rangefuel.jl"))
using .OptimizeRangeFuel: MissionReq, BoundsOpt, ConstraintsOpt, optimize_rangefuel_fun!
include(__TASOPTindices__)

####History Saving function
function save_hist_compact!(jld2_path::AbstractString, ran_cur, hist)
    tag = "range_$(round(Int, ran_cur))"

    # violation "size" per evaluation (how many constraints violated)
    nviol = Int32.(length.(hist.violations))
    feasible = isempty.(hist.violations)

    # compact test_param into a matrix (cheaper than Vector{Vector})
    n = length(hist.test_param)
    X = if n == 0
        Matrix{Float32}(undef, 0, 0)
    else
        nvar = length(hist.test_param[1])
        M = Matrix{Float32}(undef, n, nvar)
        for i in 1:n
            @inbounds M[i, :] .= Float32.(hist.test_param[i])
        end
        M
    end

    jldopen(jld2_path, "a+") do f
        f["$tag/test_param"] = X
        f["$tag/penalty"] = Float32.(hist.penalty)
        f["$tag/PFEI"] = Float32.(hist.PFEI)
        f["$tag/nviol"] = nviol
        f["$tag/feasible"] = feasible
        f["$tag/n_eval"] = length(hist.penalty)
    end
end
#######################################Main Scirpt###########################################
####Start from an initial aircraft model
# Load from
save_dir  = "ModelSaved"
load_name = "acOptimized_Bat" #jld2
# Save to
save_name = "acOptimized_BatOpt" #jld2

####Setup the test conditions
#Fixed
idx_fuel = 24 #Fuel Index: Jet Fuel(24), Ethanol(32)
rho_fuel = 817.0 #Fuel Density: Jet Fuel(817.0) (kg/m3), Ethanol(789.0) (kg/m3)
hvap_fuel = 358694.0 #Heat of Vaporization: Jet Fuel(358694.0) (J/kg), Ethanol(918187.9) (J/kg)
iters_max_opt = 1000 #Number of interations
#Sweep
range_lst = [3000,2900] #Range in nmi (Need to be turned to meter for input) #collect(3000.0:-100:300) #28 cases

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
        ####Get the inputs to the optimizers
        # Intiailized inputs
        mission_req = MissionReq()
        bound_opt = BoundsOpt()
        constraints_opt = ConstraintsOpt()
        # Modify the mission requirements
        mission_req.idx_fuel = idx_fuel #Fuel Index: Jet Fuel(24), Ethanol(32)
        mission_req.rho_fuel = rho_fuel #Fuel Density: Jet Fuel(817.0) (kg/m3), Ethanol(789.0) (kg/m3)
        mission_req.hvap_fuel = hvap_fuel #Heat of Vaporization: Jet Fuel(358694.0) (J/kg), Ethanol(918187.9) (J/kg)
        # For sweep requirement
        mission_req.range_des = (ran_cur * 1852.0)  #Design flight range (m)

        #### Run the optimization
        status_cur, hist_optim_cur = optimize_rangefuel_fun!(ac; mission_req=mission_req, bounds_opt=bound_opt, constraints_opt=constraints_opt, iters_max_opt=iters_max_opt)

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