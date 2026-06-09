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

loadpath = joinpath(@__DIR__,"../ModelSaved/testSave.jld2")
savepath = joinpath(@__DIR__,"../ModelSaved/testSave2.jld2")
@load loadpath ac
ac = update_model!(ac)
@save savepath ac