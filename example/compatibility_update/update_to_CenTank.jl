"""
This script corresponds to the DualFuelUpdated verion of the function to do the conversion of model into currently TASOPT readable model
"""

using TASOPT
using JLD2
include(__TASOPTindices__)

function update_model!(ac)
    fuse_tank = TASOPT.fuselage_tank(;(f => getproperty(ac.fuse_tank, f) for f in fieldnames(TASOPT.fuselage_tank))...)
    options = TASOPT.Options(;(f => getproperty(ac.options, f) for f in fieldnames(TASOPT.Options))...)
    wake_system = TASOPT.aerodynamics.WakeSystem(options.trefftz_config)
    ac = TASOPT.aircraft{typeof(wake_system)}(ac.name, ac.description, options,
                                              ac.parg, ac.parm, ac.para, ac.pare, ac.is_sized,
                                              ac.fuselage, fuse_tank,
                                              ac.wing, ac.htail, ac.vtail,
                                              ac.engine, ac.landing_gear,
                                              wake_system)
    return ac
end

keyNameRead = "acOpt_BatJet_CT"
keyNameSave = "acOptim_BatJet_CT"
ranges = collect(300:100:3000)
loadpath = joinpath(@__DIR__,"../ModelSaved/$(keyNameRead)/$(keyNameRead)")
savepath = joinpath(@__DIR__,"../ModelSaved/$(keyNameSave)/$(keyNameSave)")
for range = ranges
    println("Range: $(range)")
    loadpath_cur = loadpath*"$(round(Int,range)).jld2"
    savepath_cur = savepath*"$(round(Int,range)).jld2"
    println("load: $(loadpath_cur)")
    println("save: $(savepath_cur)")
    @load loadpath_cur ac
    ac = update_model!(ac)
    mkpath(dirname(savepath_cur))
    @save savepath_cur ac
end