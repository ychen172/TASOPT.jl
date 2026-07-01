"""
This script extract the input cycle parameters for Pycaso combustor emissions simulation
"""

using TASOPT
include(__TASOPTindices__)
include(joinpath(@__DIR__,"../utilities_for_postprocessing/Extract.jl"))
using .Extract: extract_combustion_inputs

#### Setup IO
model_dir = joinpath(__TASOPTroot__,"../example/ModelSaved")
read_key = "Opti_Eth_NoACT_V2_"
read_miss_idx = 1
ranges = collect(300:100:3000) #[nmi] must match with data exist
# Output
save_dir = joinpath(__TASOPTroot__,"../example/CombustorCycleSaved")
save_name = "Opti_Eth_NoACT_V2_CycIn"

#### Extract and save combustor data
save_sub_dir = joinpath(save_dir,save_name)
mkpath(save_sub_dir)
for ran_cur in ranges
    println("Extracting range $(ran_cur)")
    # Input aircraft model
    read_dir = joinpath(model_dir, read_key, read_key*"$(ran_cur).jld2")
    ac = quickload_aircraft(read_dir)
    # Read combustor inputs
    extract_combustion_inputs(ac, read_miss_idx, save_name*"_$(ran_cur)", save_sub_dir)
end
