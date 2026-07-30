"""
This script compares R1, R2, R3 performance between retrofitted case and specially design case
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
caseKeys   = ["Opti_Jet_NoACT_V2_1_ULByEng_Eth_PRD_", "Opti_Eth_NoACT_V2_2_OffdesR1R2R3_"]
caseNames  = ["Retrofitted",                          "Optimized"                        ]
ranges     = collect(300:100:3000) 
# Output directory
save_dir      = "../ModelProcessed"
save_name     = "retrofitted_specially_sized" #sub_folder will be created
# Fields to read out
const fields = [:(parm[imRange,2]),:(parm[imWpay,2]),:(parm[imPFEI,2]),
                :(parg[igWpaymax]),:(parm[imVfuel,2]),:(parg[igVfmax]),
                :(parm[imWTO,2]),:(parg[igWMTO])]

#### Create save directory
save_dir_sub  = joinpath(save_dir,save_name)
mkpath(save_dir_sub)

#### Initialization
dataset_R1 = [init_results_2Layers(length(ranges), fields) for _ in caseKeys] #[dataset[Expr][:],...]
dataset_R2 = [init_results_2Layers(length(ranges), fields) for _ in caseKeys] #[dataset[Expr][:],...]
dataset_R3 = [init_results_2Layers(length(ranges), fields) for _ in caseKeys] #[dataset[Expr][:],...]

#### Extract data for the sized missions
for (j, caseKey) in enumerate(caseKeys)
    for (i, ran) in enumerate(ranges)
        # Read in the case
        ac_dir = joinpath(model_dir,caseKey,caseKey*"$(round(Int,ran))_")
        #
        R1_file = caseKey*"$(round(Int,ran))_R1_*.jld2"
        R1_file = only(glob(R1_file, ac_dir)) #Find the only matched file
        ac = quickload_aircraft(R1_file)
        extract_acModel_compact!(ac, dataset_R1[j], i)
        #
        R2_file = caseKey*"$(round(Int,ran))_R2_*.jld2"
        R2_file = only(glob(R2_file, ac_dir)) #Find the only matched file
        ac = quickload_aircraft(R2_file)
        extract_acModel_compact!(ac, dataset_R2[j], i)
        #
        R3_file = caseKey*"$(round(Int,ran))_R3_*.jld2"
        R3_file = only(glob(R3_file, ac_dir)) #Find the only matched file
        ac = quickload_aircraft(R3_file)
        extract_acModel_compact!(ac, dataset_R3[j], i)
        #
        println("File, $(ac_dir), read successfully")
        println("Data case $(j) at range $(ran) collected successfully")
    end
end

#### Plotting
# PFEI
plot_cases_specified("Sizing Range of the Baseline Jet Fuel Case [nmi]", "R1 PFEI [J/J]", 
                     [ranges for _ in dataset_R1],
                     [d[:(parm[imPFEI,2])] for d in dataset_R1],
                     caseNames, 
                     joinpath(save_dir_sub,"R1_PFEI.png"))
#
plot_cases_specified("Sizing Range of the Baseline Jet Fuel Case [nmi]", "R2 PFEI [J/J]", 
                     [ranges for _ in dataset_R2],
                     [d[:(parm[imPFEI,2])] for d in dataset_R2],
                     caseNames, 
                     joinpath(save_dir_sub,"R2_PFEI.png"))
#
plot_cases_specified("Sizing Range of the Baseline Jet Fuel Case [nmi]", "R3 PFEI [J/J]", 
                     [ranges for _ in dataset_R3],
                     [d[:(parm[imPFEI,2])] for d in dataset_R3],
                     caseNames, 
                     joinpath(save_dir_sub,"R3_PFEI.png"))
# Range
plot_cases_specified("Sizing Range of the Baseline Jet Fuel Case [nmi]", "R1 Range [nmi]", 
                     [ranges for _ in dataset_R1],
                     [d[:(parm[imRange,2])] for d in dataset_R1] ./ 1852,
                     caseNames, 
                     joinpath(save_dir_sub,"R1_Range.png"))
#
plot_cases_specified("Sizing Range of the Baseline Jet Fuel Case [nmi]", "R2 Range [nmi]", 
                     [ranges for _ in dataset_R2],
                     [d[:(parm[imRange,2])] for d in dataset_R2] ./ 1852,
                     caseNames, 
                     joinpath(save_dir_sub,"R2_Range.png"))
#
plot_cases_specified("Sizing Range of the Baseline Jet Fuel Case [nmi]", "R3 Range [nmi]", 
                     [ranges for _ in dataset_R3],
                     [d[:(parm[imRange,2])] for d in dataset_R3] ./ 1852,
                     caseNames, 
                     joinpath(save_dir_sub,"R3_Range.png"))
# Payload Weight
plot_cases_specified("Sizing Range of the Baseline Jet Fuel Case [nmi]", "R1 Payload Weight [Ton]", 
                     [ranges for _ in dataset_R1],
                     [d[:(parm[imWpay,2])] for d in dataset_R1] ./ 9.81 ./ 1000.0,
                     caseNames, 
                     joinpath(save_dir_sub,"R1_Wpay.png"))
#
plot_cases_specified("Sizing Range of the Baseline Jet Fuel Case [nmi]", "R2 Payload Weight [Ton]", 
                     [ranges for _ in dataset_R2],
                     [d[:(parm[imWpay,2])] for d in dataset_R2] ./ 9.81 ./ 1000.0,
                     caseNames, 
                     joinpath(save_dir_sub,"R2_Wpay.png"))
#
plot_cases_specified("Sizing Range of the Baseline Jet Fuel Case [nmi]", "R3 Payload Weight [Ton]", 
                     [ranges for _ in dataset_R3],
                     [d[:(parm[imWpay,2])] for d in dataset_R3] ./ 9.81 ./ 1000.0,
                     caseNames, 
                     joinpath(save_dir_sub,"R3_Wpay.png"))
# Payload Capacity
plot_cases_specified("Sizing Range of the Baseline Jet Fuel Case [nmi]", "R1 Payload Capacity (%)", 
                     [ranges for _ in dataset_R1],
                     [(100.0 .* d[:(parm[imWpay,2])] ./ d[:(parg[igWpaymax])]) for d in dataset_R1],
                     caseNames, 
                     joinpath(save_dir_sub,"R1_Wpay_Per.png"))
#
plot_cases_specified("Sizing Range of the Baseline Jet Fuel Case [nmi]", "R2 Payload Capacity (%)", 
                     [ranges for _ in dataset_R2],
                     [(100.0 .* d[:(parm[imWpay,2])] ./ d[:(parg[igWpaymax])]) for d in dataset_R2],
                     caseNames, 
                     joinpath(save_dir_sub,"R2_Wpay_Per.png"))
#
plot_cases_specified("Sizing Range of the Baseline Jet Fuel Case [nmi]", "R3 Payload Capacity (%)", 
                     [ranges for _ in dataset_R3],
                     [(100.0 .* d[:(parm[imWpay,2])] ./ d[:(parg[igWpaymax])]) for d in dataset_R3],
                     caseNames, 
                     joinpath(save_dir_sub,"R3_Wpay_Per.png"))
# Fuel Capacity
plot_cases_specified("Sizing Range of the Baseline Jet Fuel Case [nmi]", "R1 Fuel Volume Capacity (%)", 
                     [ranges for _ in dataset_R1],
                     [(100.0 .* d[:(parm[imVfuel,2])] ./ d[:(parg[igVfmax])]) for d in dataset_R1],
                     caseNames, 
                     joinpath(save_dir_sub,"R1_VFuel_Per.png"))
#
plot_cases_specified("Sizing Range of the Baseline Jet Fuel Case [nmi]", "R2 Fuel Volume Capacity (%)", 
                     [ranges for _ in dataset_R2],
                     [(100.0 .* d[:(parm[imVfuel,2])] ./ d[:(parg[igVfmax])]) for d in dataset_R2],
                     caseNames, 
                     joinpath(save_dir_sub,"R2_VFuel_Per.png"))
#
plot_cases_specified("Sizing Range of the Baseline Jet Fuel Case [nmi]", "R3 Fuel Volume Capacity (%)", 
                     [ranges for _ in dataset_R3],
                     [(100.0 .* d[:(parm[imVfuel,2])] ./ d[:(parg[igVfmax])]) for d in dataset_R3],
                     caseNames, 
                     joinpath(save_dir_sub,"R3_VFuel_Per.png"))
# WTO Capacity
plot_cases_specified("Sizing Range of the Baseline Jet Fuel Case [nmi]", "R1 Takeoff Weight Capacity (%)", 
                     [ranges for _ in dataset_R1],
                     [(100.0 .* d[:(parm[imWTO,2])] ./ d[:(parg[igWMTO])]) for d in dataset_R1],
                     caseNames, 
                     joinpath(save_dir_sub,"R1_WTO_Per.png"))
#
plot_cases_specified("Sizing Range of the Baseline Jet Fuel Case [nmi]", "R2 Takeoff Weight Capacity (%)", 
                     [ranges for _ in dataset_R2],
                     [(100.0 .* d[:(parm[imWTO,2])] ./ d[:(parg[igWMTO])]) for d in dataset_R2],
                     caseNames, 
                     joinpath(save_dir_sub,"R2_WTO_Per.png"))
#
plot_cases_specified("Sizing Range of the Baseline Jet Fuel Case [nmi]", "R3 Takeoff Weight Capacity (%)", 
                     [ranges for _ in dataset_R3],
                     [(100.0 .* d[:(parm[imWTO,2])] ./ d[:(parg[igWMTO])]) for d in dataset_R3],
                     caseNames, 
                     joinpath(save_dir_sub,"R3_WTO_Per.png"))