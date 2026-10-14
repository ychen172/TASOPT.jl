"""
This script compares parameters across multiple cases
"""

using TASOPT
include(__TASOPTindices__)
using Plots
include(joinpath(@__DIR__,"../UtilityFunctions/Extract.jl"))
using .Extract: extract_acModel_compact!, init_results_2Layers, plot_cases_specified
include(joinpath(@__DIR__,"../UtilityFunctions/Breguet_range_solve_offdes.jl"))
using .Breguet: Bre_off_des

#### Setup IO
# Input case names - Retrofit
model_dir  = "../ModelSaved"
caseKeys   = ["Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_","Opti_Eth_NoACT_V4_3_R1Sz_EtasEng_"]
caseNames  = ["Jet Fuel"                         ,"Ethanol"                          ]
ranges     = collect(300:100:3000) 
# Output directory
save_dir      = "../ModelProcessed"
save_name     = "Aug_3_2026_Jet_Fuel_vs_Ethanol_Conventional" #sub_folder will be created
# Fields to read out
const fields = [:(parm[imRange,1]),:(parm[imVfuel,1])]

#### Create save directory
save_dir_sub  = joinpath(save_dir,save_name)
mkpath(save_dir_sub)

#### Initialization
dataset = [init_results_2Layers(length(ranges), fields) for _ in caseKeys] #[dataset[Expr][:],...]

#### Extract data for the sized missions
for (j, caseKey) in enumerate(caseKeys)
    for (i, ran) in enumerate(ranges)
        # Read in the case
        ac_dir = joinpath(model_dir,caseKey,caseKey*"$(round(Int,ran)).jld2")
        ac = quickload_aircraft(ac_dir)
        println("File, $(ac_dir), read successfully")
        # Extract data
        extract_acModel_compact!(ac, dataset[j], i)
        println("Data case $(j) at range $(ran) collected successfully")
    end
end

#### Plotting
plot_cases_specified("Design/R1 Range [nmi]", "Fuel Level [%]", 
                     [d[:(parm[imRange,1])] for d in dataset] ./ 1852.0, 
                     [d[:(parm[imVfuel,1])] .* 100.0 ./ d[:(parg[igVfmax])] for d in dataset],
                     caseNames, 
                     joinpath(save_dir_sub,"Fuel_Level.png"))
