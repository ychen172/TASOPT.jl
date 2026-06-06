"""
This script reset regression test with the default_sized and default_structures used for the test
Remember to update the reset_regression_test inside save_model.jl for additional output variables created
"""

using TASOPT
using Printf
include(TASOPT.__TASOPTindices__)
include(joinpath(TASOPT.__TASOPTroot__, "../src/IO/save_model.jl"))

#### Change for current reset
name_change_for_this_update = "change_rfmax_Jun5_26"

#### Load the run using the current model
ac = load_default_model()
@assert ac.fuselage.layout.radius == 1.9558 #Get the right model
size_aircraft!(ac; printiter=false)

#### Backup the old defaults
testdir = joinpath(TASOPT.__TASOPTroot__, "../test")
zipfile = joinpath(testdir,"backup_default_and_aero_BEFORE_$(name_change_for_this_update).zip")

run(`zip -j $zipfile \
    $(joinpath(testdir, "default_sized.jl")) \
    $(joinpath(testdir, "default_structures.jl")) \
    $(joinpath(testdir, "unit_test_aero.jl")) \
    $(joinpath(testdir, "weights.txt")) \
    $(joinpath(testdir, "aero.txt")) \
    $(joinpath(testdir, "geom.txt"))
    `)

#### Save the new defaults
reset_regression_test(ac)

#### Compute the offdesign
# Default input and offdesign
TASOPT.fly_mission!(ac, 2; printTO=false);
println("Use new `test_ac_off_design(ac, $(ac.parm[imPFEI, 2]), $(ac.parm[imWfuel, 2]),  $(ac.parm[imWTO, 2])`")
println("Use new PFEI_default: $(ac.parm[imPFEI])")
# Default wide
ac_W = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_wide.toml"))
@assert ac_W.fuselage.layout.radius ≈ 3.0988
size_aircraft!(ac_W; printiter=false);
println("Use new PFEI_wide: $(ac_W.parm[imPFEI])")
# Default regional
ac_R = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_regional.toml"))
@assert ac_R.fuselage.layout.radius ≈ 1.5113
size_aircraft!(ac_R; printiter=false);
println("Use new PFEI_region: $(ac_R.parm[imPFEI])")
# Default Hydrogen
ac_H = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/cryo_input.toml"))
@assert ac_H.fuselage.layout.radius ≈ 2.54
size_aircraft!(ac_H, iter=50; printiter=false);
println("Use new PFEI_hydrogen: $(ac_H.parm[imPFEI])")

#### Reset aeroperf_sweep
results = aeroperf_sweep(ac, [0.0, 0.8])
aerofields = (
   CLs      = nothing,
   CDs      = nothing,
   LDs      = nothing,
   CLhs     = nothing,
   CDis     = nothing,
   CDwings  = nothing,
   CDfuses  = nothing,
   CDhtails = nothing,
   CDvtails = nothing,
   CDothers = nothing,
   clpos    = nothing,
   clpss    = nothing,
   clpts    = nothing,
   cdfss    = nothing,
   cdpss    = nothing,
   cdwss    = nothing,
   cdss     = nothing,
)
for field in keys(aerofields)
   println(field, " = ", getproperty(results, field))
end

#### Replace print out Text
f = open(io->TASOPT.weight_buildup(ac,io=io), "weights.txt", "w")
f = open(io->TASOPT.aero(ac,io=io), "aero.txt", "w")
f = open(io->TASOPT.geometry(ac,io=io), "geom.txt", "w")