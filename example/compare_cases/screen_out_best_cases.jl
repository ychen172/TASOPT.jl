"""
This script compares multiple version of optimization and select the best one for each mission and collect them all into one folder
"""

using TASOPT
include(__TASOPTindices__)
using Plots
include(joinpath(@__DIR__,"../utilities_for_postprocessing/Extract.jl"))
using .Extract: extract_acModel_compact!, init_results_2Layers, plot_cases_specified

#### Setup IO
# Input case names - Retrofit
model_dir  = "../ModelSaved"
caseKeys   = ["Opti_Jet_NoACT_V3_R1Size_","Opti_Jet_NoACT_V3_1_R1Size_","Opti_Jet_NoACT_V3_2_R1Size_","Opti_Jet_NoACT_V3_3_R1Size_","Opti_Jet_NoACT_V3_4_R1Size_","Opti_Jet_NoACT_V3_5_R1Size_","Opti_Jet_NoACT_V3_6_R1Size_","Opti_Jet_NoACT_V3_7_R1Size_"]
caseNames  = ["Jet Fuel 0",               "Jet Fuel 1",                 "Jet Fuel2",                   "Jet Fuel 3",                "Jet Fuel 4",                 "Jet Fuel 5",                 "Jet Fuel 6",                 "Jet Fuel 7"                 ]
ranges     = collect(300:100:3000) 
# Output directory
save_name  = "Opti_Jet_NoACT_V3_8_R1Size_" #Save the new model to
# Fields to read out
const fields = [:(parm[imPFEI,1])]

#### Create save directory
save_dir  = joinpath(model_dir,save_name)
mkpath(save_dir)

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

#### Screen out the best performance at each range
PFEI_collected = [d[:(parm[imPFEI,1])] for d in dataset] #nCases*nRange
PFEI_collected = reduce(hcat, PFEI_collected)
idx_col_best = map(argmin, eachrow(PFEI_collected))  

#### Save the copy the corresponding best model from each case into the new folder
for (idx,idx_best) in enumerate(idx_col_best)
    ran = round(Int, ranges[idx])
    #
    src = joinpath(model_dir,
                   caseKeys[idx_best],
                   caseKeys[idx_best] * "$(ran)_design_constraints.csv")
    dst = joinpath(save_dir,
                   save_name * "$(ran)_design_constraints.csv")
    cp(src, dst; force=true)
    #
    src = joinpath(model_dir,
                   caseKeys[idx_best],
                   caseKeys[idx_best] * "$(ran)_global_bounds.csv")
    dst = joinpath(save_dir,
                   save_name * "$(ran)_global_bounds.csv")
    cp(src, dst; force=true)
    #
    src = joinpath(model_dir,
                   caseKeys[idx_best],
                   caseKeys[idx_best] * "$(ran)_mission_requirements.csv")
    dst = joinpath(save_dir,
                   save_name * "$(ran)_mission_requirements.csv")
    cp(src, dst; force=true)
    #
    src = joinpath(model_dir,
                   caseKeys[idx_best],
                   caseKeys[idx_best] * "$(ran)_optimization_history.jld2")
    dst = joinpath(save_dir,
                   save_name * "$(ran)_optimization_history.jld2")
    cp(src, dst; force=true)
    #
    src = joinpath(model_dir,
                   caseKeys[idx_best],
                   caseKeys[idx_best] * "$(ran)_optimized_parameters.csv")
    dst = joinpath(save_dir,
                   save_name * "$(ran)_optimized_parameters.csv")
    cp(src, dst; force=true)
    #
    src = joinpath(model_dir,
                   caseKeys[idx_best],
                   caseKeys[idx_best] * "$(ran)_OptLog.txt")
    dst = joinpath(save_dir,
                   save_name * "$(ran)_OptLog.txt")
    cp(src, dst; force=true)
    #
    src = joinpath(model_dir,
                   caseKeys[idx_best],
                   caseKeys[idx_best] * "$(ran).jld2")
    dst = joinpath(save_dir,
                   save_name * "$(ran).jld2")
    cp(src, dst; force=true)
end

