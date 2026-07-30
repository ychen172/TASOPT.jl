"""
This script compares PRD performance between various cases at a single design range
"""

using Glob
using TASOPT
include(__TASOPTindices__)
using Plots
include(joinpath(@__DIR__,"../utilities_for_postprocessing/Extract.jl"))
using .Extract: extract_acModel_compact!, init_results_2Layers, plot_cases_specified

#### Setup IO
# Input case names - Retrofit
model_dir  = "../ModelSaved"
caseKeys   = ["Opti_Jet_NoACT_V2_1_ULByEng_Aft3_PRD_","Opti_Jet_NoACT_V2_1_ULByEng_Eth_PRD_", "Opti_Eth_NoACT_V2_2_MatchR1R2R3_PRD_"]
caseNames  = ["Baseline Jet Fuel"                    ,"Retrofitted Jet Fuel"                , "Optimized Ethanol"                   ]
desran     = 3000
offdesran  = collect(300:100:8000) 
# Output directory
save_dir      = "../ModelProcessed"
save_name     = "retrofitted_specially_sized_PRD" #sub_folder will be created
# Fields to read out
const fields = [:(parm[imRange,2]),:(parm[imWpay,2])]

#### Create save directory
save_dir_sub  = joinpath(save_dir,save_name)
mkpath(save_dir_sub)

#### Initialization
dataset = [init_results_2Layers(length(offdesran), fields) for _ in caseKeys] #[dataset[Expr][:],...]

#### Extract data for the sized missions
for (j, caseKey) in enumerate(caseKeys)
    for (i, ran) in enumerate(offdesran)
        # Read in the case
        ac_dir = joinpath(model_dir,caseKey,caseKey*"$(round(Int,desran))_",caseKey*"$(round(Int,desran))_$(round(Int,ran)).jld2")
        #
        try
            ac = quickload_aircraft(ac_dir)
            extract_acModel_compact!(ac, dataset[j], i)
            println("File, $(ac_dir), read successfully")
            println("Data case $(j) at range $(ran) collected successfully")
        catch
            #
        end
    end
end

#### Plotting
# PRD
plot_cases_specified("Off-design Ranges [nmi]", "Payload Weight [Ton]", 
                     [d[:(parm[imRange,2])] ./ 1852.0 for d in dataset],
                     [d[:(parm[imWpay,2])] ./ 9.81 ./ 1000.0 for d in dataset],
                     caseNames, 
                     joinpath(save_dir_sub,"PRD_$(desran).png"))