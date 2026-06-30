"""
This script compare the payload range envelop aross different design cases
"""

using TASOPT
include(__TASOPTindices__)
using Plots
include(joinpath(@__DIR__,"../utilities_for_postprocessing/Extract.jl"))
using .Extract: extract_acModel, init_results_2Layers, fill_results!, plot_cases
using Glob

#### Setup IO
read_dir = "../ModelSaved" 
ran_design = 1500
key_names = ["Opti_Jet_NoACT_to_Eth_OffDes_",      "Opti_Jet_NoACT_to_Eth_for_Eth_OffDes_",   "Opti_Eth_NoACT_Mat_Jet_NoACT_to_Eth_OffDes_"]
R1R2R3_names = ["Opti_Jet_NoACT_to_Eth_MatR1R2R3_","Opti_Jet_NoACT_to_Eth_for_Eth_MatR1R2R3_","Opti_Eth_NoACT_Mat_Jet_NoACT_to_Eth_MatR1R2R3_"]
lables    = ["Retrofit","Optimized through Retrofitting", "Optimized Directly"]
rans_offdes = collect(300:100:4000)

# Output directory
save_dir = "../ModelProcessed"
save_name = "Compare_Retrofit_Design_PRD"
# Fields to read out
const fields = (:range_nmi,:massPayload_Ton,:PFEI_JJ)

#### Create save directory
save_dir_sub  = joinpath(save_dir,save_name)
mkpath(save_dir_sub)

#### Initialization
results = [init_results_2Layers(length(rans_offdes), fields) for _ in key_names]
R1R2R3Restuls = [init_results_2Layers(length(rans_offdes), fields) for _ in [1,2,3]]
#### Extract data for the design mission
for (i, key_cur) in enumerate(key_names)
    j_last = 0
    for (j, ran_cur) in enumerate(rans_offdes)
        try
            # Read in aircraft model
            ac_dir = joinpath(read_dir,key_cur,key_cur*"$(ran_design)_",key_cur*"$(ran_design)_$(ran_cur).jld2")
            ac = quickload_aircraft(ac_dir)
            # Extract parameters
            results[i][:range_nmi][j] = ac.parm[imRange, 2] / 1852.0 #[nmi]
            results[i][:massPayload_Ton][j] = ac.parm[imWpay, 2] / gee / 1000.0 #[Ton]
            results[i][:PFEI_JJ][j] = ac.parm[imPFEI, 2]
            j_last = j
            # Read in aircraft models for R1 R2 R3 missions
            R1R2R3_dir = joinpath(read_dir,
                                  R1R2R3_names[i],
                                  R1R2R3_names[i]*"$(ran_design)_")
            files = glob(R1R2R3_names[i]*"$(ran_design)_R1Mat_*.jld2", R1R2R3_dir)
            @assert length(files) == 1 "Expected exactly one matching file."
            ac_dir = files[1]
            ac = quickload_aircraft(ac_dir)
            R1R2R3Restuls[i][:range_nmi][1] = ac.parm[imRange, 2] / 1852.0 #[nmi]
            R1R2R3Restuls[i][:massPayload_Ton][1] = ac.parm[imWpay, 2] / gee / 1000.0 #[Ton]
            R1R2R3Restuls[i][:PFEI_JJ][1] = ac.parm[imPFEI, 2]
            
            R1R2R3_dir = joinpath(read_dir,
                                  R1R2R3_names[i],
                                  R1R2R3_names[i]*"$(ran_design)_")
            files = glob(R1R2R3_names[i]*"$(ran_design)_R2Mat_*.jld2", R1R2R3_dir)
            @assert length(files) == 1 "Expected exactly one matching file."
            ac_dir = files[1]
            ac = quickload_aircraft(ac_dir)
            R1R2R3Restuls[i][:range_nmi][2] = ac.parm[imRange, 2] / 1852.0 #[nmi]
            R1R2R3Restuls[i][:massPayload_Ton][2] = ac.parm[imWpay, 2] / gee / 1000.0 #[Ton]
            R1R2R3Restuls[i][:PFEI_JJ][2] = ac.parm[imPFEI, 2]

            R1R2R3_dir = joinpath(read_dir,
                                  R1R2R3_names[i],
                                  R1R2R3_names[i]*"$(ran_design)_")
            files = glob(R1R2R3_names[i]*"$(ran_design)_R3Mat_*.jld2", R1R2R3_dir)
            @assert length(files) == 1 "Expected exactly one matching file."
            ac_dir = files[1]
            ac = quickload_aircraft(ac_dir)
            R1R2R3Restuls[i][:range_nmi][3] = ac.parm[imRange, 2] / 1852.0 #[nmi]
            R1R2R3Restuls[i][:massPayload_Ton][3] = ac.parm[imWpay, 2] / gee / 1000.0 #[Ton]
            R1R2R3Restuls[i][:PFEI_JJ][3] = ac.parm[imPFEI, 2]

            println("File, $(ac_dir), read successfully")
        catch e
            continue
        end
    end
    for field in keys(results[i])
        resize!(results[i][field], j_last)
    end
end

#### Compare the PRD
markers = [:circle, :square, :diamond, :utriangle, :dtriangle, :star5]
linestyles = [:solid, :dash, :dot, :dashdot]
p = plot(xlabel="Off-design Range (nmi)", ylabel="Payload Weight (Ton)", dpi=800)
for (i, key_cur) in enumerate(key_names)
    plot!(p, results[i][:range_nmi], results[i][:massPayload_Ton], marker=markers[i], linestyle=linestyles[i], lw=2, markerstrokewidth=0, label=lables[i])
    scatter!(p, R1R2R3Restuls[i][:range_nmi], R1R2R3Restuls[i][:massPayload_Ton], marker=markers[i], linestyle=linestyles[i], lw=2, markerstrokewidth=0, label=lables[i]*" R1,R2,R3")
end
savefig(p, joinpath(save_dir_sub, "PRD.png"))

markers = [:circle, :square, :diamond, :utriangle, :dtriangle, :star5]
linestyles = [:solid, :dash, :dot, :dashdot]
p = plot(xlabel="Off-design Range (nmi)", ylabel="PFEI (J/J)", dpi=800)
for (i, key_cur) in enumerate(key_names)
    plot!(p, results[i][:range_nmi], results[i][:PFEI_JJ], marker=markers[i], linestyle=linestyles[i], lw=2, markerstrokewidth=0, label=lables[i])
    # scatter!(p, R1R2R3Restuls[i][:range_nmi], R1R2R3Restuls[i][:PFEI_JJ], marker=markers[i], linestyle=linestyles[i], lw=2, markerstrokewidth=0, label=lables[i]*" R1,R2,R3")
end
savefig(p, joinpath(save_dir_sub, "PFEI.png"))