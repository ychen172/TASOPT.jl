"""
This script extract parameters and print out for designed mission
"""

# 1) Load TASOPT
using TASOPT
using DataFrames, CSV
using Plots
include(__TASOPTindices__)
include(joinpath(@__DIR__, "postprocess.jl"))
using .PostProcess: ExtractDes
using .PostProcess.OptimizeRangeFuel: BoundsOpt, ConstraintsOpt

####Setup IO
load_dir = "ModelSaved"
load_name = "acOptimized_BatOptEth400"
save_dir = "ModelProcessed"
save_name = "BatOptEth400"
mkpath(load_dir)
mkpath(save_dir)

####Load model
ac = quickload_aircraft(joinpath(load_dir,"$(load_name).jld2"))

####Setup constraints and search bound (Need to be saved by optimization later)
bounds_opt = BoundsOpt()
constraints_opt = ConstraintsOpt() #Current setup use default values

####Extract and save the performance data
designParam = ExtractDes(ac, save_dir, save_name; flg_save=true, bounds_opt=bounds_opt, constraints_opt=constraints_opt)
println(designParam)