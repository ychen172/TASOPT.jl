using TASOPT
using JLD2
include(__TASOPTindices__)
loadpath = joinpath(@__DIR__,"../ModelSaved/testSave.jld2")
savepath = joinpath(@__DIR__,"../ModelSaved/testSave2.jld2")
@load loadpath ac
# ac = aircraft_tmp(ac)
# update_model!(ac)
# @save savepath ac