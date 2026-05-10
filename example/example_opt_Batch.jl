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
ac = read_aircraft_model("./customized/narrow_input.toml"; templatefile = "./customized/narrow_input.toml")
size_aircraft!(ac)

####Modify the mission requirement for optimization (ToDo)

####Get the objective function for optimization
hist_optim = OptHistory() #Optimization history
## Setup constraints
max_span = 35.814 #Maximum Span(m)
max_lenField = 2438.4 #Maximum Balanced Field Length(m)
min_TOCGamma = 0.015 #Minimum Top of Climb Flight Angle (rad)
max_Tt3 = 900.0 #Maximum Combustor Inlet Temperature (K)
max_TMetal = 1333.33 #Maximum Metal Tempeature (K)
max_DiaFan = 2.0 #Maximum Fan Diameter (m)
constraints_optim = [max_span, max_lenField, min_TOCGamma, max_Tt3, max_TMetal, max_DiaFan]
## Get objective function
obj = make_obj(ac, constraints_optim, hist_optim)

