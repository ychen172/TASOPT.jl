"""
This script extract the combustor operating condtitions for R1S sizing point mission over aircraft with various design ranges.
"""

using TASOPT
include(__TASOPTindices__)
using Plots
include(joinpath(__TASOPTroot__,"../example/UtilityFunctions/objective_factory.jl"))
off_design_specified! = ObjectiveFactory.off_design_specified!
include(joinpath(__TASOPTroot__,"../example/UtilityFunctions/Extract.jl"))
using .Extract: read_oag,extract_combustion_inputs

#### Setup IO
# Input case names - OAG seat-capacity sweep
model_dir  = joinpath(__TASOPTroot__,"../example/ModelSaved")
caseKey   = "Opti_Eth_NoACT_V4_3_R1Sz_EtasEng_"
seat_ranges = collect(300:100:3000)
# Output folder name (within model_dir)
save_name     = "Combustor_$(caseKey)_FLI" #sub_folder will be created

#### Create save directory
save_dir_sub  = joinpath(model_dir,save_name)
mkpath(save_dir_sub)

#### Compute the oag missions engine and combustor performance and save
for (i, sr) in enumerate(seat_ranges)
    # Read in the case
    ac_dir = joinpath(model_dir,caseKey,caseKey*"$(sr).jld2")
    ac = quickload_aircraft(ac_dir)
    println("File, $(ac_dir), read successfully")
    # Extract the combustor operating condition
    extract_combustion_inputs(ac,1,save_name*"$(sr)",save_dir_sub)
end