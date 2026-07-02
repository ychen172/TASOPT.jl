"""
This script compares parameters across multiple cases
"""

using TASOPT
include(__TASOPTindices__)
using Plots
include(joinpath(@__DIR__,"../utilities_for_postprocessing/Extract.jl"))
using .Extract: extract_acModel, init_results_2Layers, fill_results!, plot_cases
using Glob

#### Setup IO
# Input case names - Retrofit
model_dir  = "../ModelSaved"
caseKeys   = ["Opti_Jet_NoACT_to_Eth_MatR1R2R3_",  "Opti_Jet_NoACT_to_Eth_for_Eth_MatR1R2R3_", "Opti_Eth_NoACT_Mat_Jet_NoACT_to_Eth_MatR1R2R3_"]
caseNames  = ["Retrofitted Ethanol",               "Optimized by Retrofitting",                "Optimized Directly"]
idxMiss    = [2,2,2]
ranges     = collect(300:100:3000) 
# Output directory
save_dir      = "../ModelProcessed"
save_name     = "Retrofit_vs_Optimized" #sub_folder will be created
# Fields to read out
const fields = (:ran_des_nmi, :span_wing_m ,:sweep_wing_deg ,:AR_wing ,:ThiCho_wing_in ,:ThiCho_wing_out ,:taper_wing_in,
                :taper_wing_out ,:y_engine_m ,:x_engine_m ,:SpanBreak_wing ,:cls_wing ,:clt_wing ,:CL_cruise,
                :Alt_cruise_ft ,:PR_fan ,:PR_LPC ,:PR_HPC ,:BPR ,:Tt4 ,:AR_vtail)

#### Create save directory
save_dir_sub  = joinpath(save_dir,save_name)
mkpath(save_dir_sub)

#### Initialization
results = [init_results_2Layers(length(ranges), fields) for _ in caseKeys]

#### Extract data for the sized missions
for (i, caseKey) in enumerate(caseKeys)
    for (j, ran) in enumerate(ranges)
        # Read in the case
        ac_dir = joinpath(model_dir,caseKey,caseKey*"$(round(Int,ran))_")
        files = glob(caseKey*"$(round(Int,ran))_R2Mat_*.jld2", ac_dir)
        if length(files)!=1
            println("directory: $(ac_dir)")
            println("file name: "*caseKey*"$(round(Int,ran))_R2Mat_*.jld2")
            error("Expected exactly one matching file.")
        end
        ac_dir = files[1]
        ac = quickload_aircraft(ac_dir)
        println("File, $(ac_dir), read successfully")
        # Extract data
        out = extract_acModel(ac,idxMiss[i])
        # Store the data
        fill_results!(results[i],out,j)
    end
    results[i][:ran_des_nmi] = ranges
end

plot_cases("Reference Range (nmi)","Wing Span (m)",results,:ran_des_nmi,:span_wing_m,caseNames,joinpath(save_dir_sub, "WingSpan.png"))
plot_cases("Reference Range (nmi)","Wing Sweep (deg)",results,:ran_des_nmi,:sweep_wing_deg,caseNames,joinpath(save_dir_sub, "WingSweep.png"))
plot_cases("Reference Range (nmi)","Wing Aspect Ratio",results,:ran_des_nmi,:AR_wing,caseNames,joinpath(save_dir_sub, "WingAR.png"))
plot_cases("Reference Range (nmi)","Wing Inner Thickness to Chord",results,:ran_des_nmi,:ThiCho_wing_in,caseNames,joinpath(save_dir_sub, "WingTCIn.png"))
plot_cases("Reference Range (nmi)","Wing Outer Thickness to Chord",results,:ran_des_nmi,:ThiCho_wing_out,caseNames,joinpath(save_dir_sub, "WingTCOut.png"))
plot_cases("Reference Range (nmi)","Wing Inner Taper",results,:ran_des_nmi,:taper_wing_in,caseNames,joinpath(save_dir_sub, "WingTaperIn.png"))
plot_cases("Reference Range (nmi)","Wing Outer Taper",results,:ran_des_nmi,:taper_wing_out,caseNames,joinpath(save_dir_sub, "WingTaperOut.png"))
plot_cases("Reference Range (nmi)","Spanwise Engine Loc(m)",results,:ran_des_nmi,:y_engine_m,caseNames,joinpath(save_dir_sub, "yEngine.png"))
plot_cases("Reference Range (nmi)","Axial Engine Loc(m)",results,:ran_des_nmi,:x_engine_m,caseNames,joinpath(save_dir_sub, "xEngine.png"))
plot_cases("Reference Range (nmi)","Wing Span Break",results,:ran_des_nmi,:SpanBreak_wing,caseNames,joinpath(save_dir_sub, "WingSpanBreak.png"))
plot_cases("Reference Range (nmi)","Wing cl ratio Span to Root",results,:ran_des_nmi,:cls_wing,caseNames,joinpath(save_dir_sub, "Wingcls.png"))
plot_cases("Reference Range (nmi)","Wing cl ratio Tip to Root",results,:ran_des_nmi,:clt_wing,caseNames,joinpath(save_dir_sub, "Wingclt.png"))
plot_cases("Reference Range (nmi)","CL cruise",results,:ran_des_nmi,:CL_cruise,caseNames,joinpath(save_dir_sub, "CLcruise.png"))
plot_cases("Reference Range (nmi)","Alt cruise(ft)",results,:ran_des_nmi,:Alt_cruise_ft,caseNames,joinpath(save_dir_sub, "Altcruise.png"))
plot_cases("Reference Range (nmi)","Fan PR",results,:ran_des_nmi,:PR_fan,caseNames,joinpath(save_dir_sub, "PR_fan.png"))
plot_cases("Reference Range (nmi)","LPC PR",results,:ran_des_nmi,:PR_LPC,caseNames,joinpath(save_dir_sub, "PR_LPC.png"))
plot_cases("Reference Range (nmi)","HPC PR",results,:ran_des_nmi,:PR_HPC,caseNames,joinpath(save_dir_sub, "PR_HPC.png"))
plot_cases("Reference Range (nmi)","BPR",results,:ran_des_nmi,:BPR,caseNames,joinpath(save_dir_sub, "BPR.png"))
plot_cases("Reference Range (nmi)","Tt4(K)",results,:ran_des_nmi,:Tt4,caseNames,joinpath(save_dir_sub, "Tt4.png"))
plot_cases("Reference Range (nmi)","AR Vtail",results,:ran_des_nmi,:AR_vtail,caseNames,joinpath(save_dir_sub, "ARvtail.png"))