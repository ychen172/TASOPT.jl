"""
This script compare the design point performance between ethanol and jet
Fixed design bounds and Fixed constraints across different cases assumed
"""

using TASOPT
using DataFrames, CSV
include(__TASOPTindices__)
include(joinpath(@__DIR__, "postprocess.jl"))
using .PostProcess: ExtractDes
using .PostProcess.OptimizeRangeFuel: BoundsOpt, ConstraintsOpt

#### Setup IO
# I
model_dir   = "ModelSaved"
model_name  = "acOptimized_BatOptJet1500" #Frontal key name for the models(FuelRange)
bounds_name = "acOptimized_BatOptJet2800_BoundLocal" #.csv the bounds for the optimized parameters
# O
save_dir    = "ModelProcessed" #Outer Directory for saving models
save_name   = "BatOptJet1500" #The total save name for jet and ethanol comparison
mkpath(model_dir)
mkpath(save_dir)
save_dirSub = joinpath(save_dir,save_name) #Sub-directory to save comparison data
mkpath(save_dirSub)

#### Assumed default constraints used
constraints_opt = ConstraintsOpt()

#### Loading the bounds for optimization
bounds_opt = BoundsOpt()
df = CSV.read(joinpath(model_dir,"$(bounds_name).csv"), DataFrame)
fNames = fieldnames(BoundsOpt)
for fNameCur in fNames
    col = Symbol(fNameCur)
    vals = (Float64(df[1, col]), Float64(df[2, col]), Float64(df[3, col]))
    setfield!(bounds_opt, fNameCur, vals)
end

#### Setup parameters
ac = quickload_aircraft(joinpath(model_dir,"$(model_name).jld2"))
designParam = ExtractDes(ac, save_dirSub, save_name; flg_save=true, bounds_opt=bounds_opt, constraints_opt=constraints_opt)