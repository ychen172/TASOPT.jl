using TASOPT, NLopt
include(joinpath(@__DIR__, "objective_factory.jl"))
using .ObjectiveFactory: OptHistory, make_obj
include(__TASOPTindices__)

# hist = OptHistory()
# obj  = make_obj(...)

####IO Prepare
# Make a folder for saving optimized aircraft model
save_dir = "ModelSaved"
mkpath(save_dir)

# Make a save name for the optimized model
save_Name = "acOptimized_Cus"

####Loading an Baseline Aircraft Model
ac = read_aircraft_model("./customized/narrow_input.toml")
size_aircraft!(ac)

####Modify the mission requirement for optimization (ToDo)

####Get the objective function for optimization

