"""This script test the sensitivity for a group of parameters around optimized solution"""

using TASOPT
include(TASOPT.__TASOPTindices__)
include(joinpath(__TASOPTroot__,"utils","sensitivity.jl"))

modelPath = joinpath(@__DIR__,"../ModelSaved/acOptim_BatJet_CT/acOptim_BatJet_CT3000.jld2")
println(modelPath)
ac = quickload_aircraft(modelPath)
input_params = [:(ac.wing.layout.sweep)]
sensitivityVector = get_sensitivity(input_params; model_state=ac, eps=1e-3, optimizer=false, f_out_fn=nothing, diff_scheme=:backward, metric=:impact)
println(sensitivityVector)