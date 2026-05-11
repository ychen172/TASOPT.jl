using TASOPT, NLopt
include(joinpath(@__DIR__, "objective_factory.jl"))
include(joinpath(@__DIR__, "optimize_rangefuel.jl"))
using .ObjectiveFactory: OptHistory, make_obj, best_feasible
using .OptimizeRangeFuel: MissionReq, BoundsOpt, ConstraintsOpt, optimize_rangefuel_fun!
include(__TASOPTindices__)

####Start an initial design case
save_dir  = "ModelSaved"
save_Name = "acOptimized_Bat" #jld2
load_dir  = "customized"
load_Name = "narrow_input" #toml
# IO
mkpath(save_dir)
ac = read_aircraft_model(joinpath(load_dir,load_Name*".toml"); templatefile = joinpath(load_dir,load_Name*".toml"))
size_aircraft!(ac)
# Get optmization input
mission_req = MissionReq()
bound_opt = BoundsOpt()
constraints_opt = ConstraintsOpt()
# Modify missions
mission_req.range_des = (3000.0 * 1852.0)  #Design flight range (m)
mission_req.idx_fuel = 24 #Fuel Index: Jet Fuel(24), Ethanol(32)
mission_req.rho_fuel = 817.0 #Fuel Density: Jet Fuel(817.0) (kg/m3), Ethanol(789.0) (kg/m3)
mission_req.hvap_fuel = 358694.0 #Heat of Vaporization: Jet Fuel(358694.0) (J/kg), Ethanol(918187.9) (J/kg)
# Run
status, hist_optim = optimize_rangefuel_fun!(ac; mission_req=mission_req, bounds_opt=bound_opt, constraints_opt=constraints_opt, maxeval=1000)
# Store result
# good_status = status in (:SUCCESS, :STOPVAL_REACHED, :FTOL_REACHED, :XTOL_REACHED, :MAXEVAL_REACHED, :MAXTIME_REACHED)
# bad_status  = status in (:FAILURE, :INVALID_ARGS, :OUT_OF_MEMORY, :ROUNDOFF_LIMITED, :FORCED_STOP)
if status in (:SUCCESS, :STOPVAL_REACHED, :FTOL_REACHED, :XTOL_REACHED, :MAXEVAL_REACHED, :MAXTIME_REACHED)
    quicksave_aircraft(ac,joinpath(save_dir, "$(save_Name).jld2"))
end