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
caseKeys   = ["Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_PRD_", "Opti_Eth_NoACT_V4_3_R1Sz_EtasEng_PRD_"]
caseNams   = ["Jet Fuel"                         , "Ethanol"                          ]
desRans    = [1000,2000,3000]
offDesRans = collect(300:100:8000) #Ensure that the first range is flyable across all the cases and all the design ranges
# Output directory
save_dir   = "../ModelProcessed"
save_name  = "Aug_3_2026_Jet_Fuel_vs_Ethanol_Conventional" #sub_folder will be created
# Fields to read out
const fields = [:(parm[imRange,2]),:(parm[imWpay,2]),:(parm[imPFEI,2])]

#### Create save directory
save_dir_sub  = joinpath(save_dir,save_name)
mkpath(save_dir_sub)

#### Initialization
dataset = [[init_results_2Layers(length(offDesRans), fields) for _ in desRans] for _ in caseKeys] #dataset[idxCase][idxRange][ExprForField][IdxOffDesRange]
datasetR1R2R3 = [[init_results_2Layers(3, fields) for _ in desRans] for _ in caseKeys] #dataset[idxCase][idxRange][ExprForField][index for R1 R2 R3]

#### Extract data for the sized missions
for (idxCase, caseKey) in enumerate(caseKeys)
    for (idxDesRan, desRange) in enumerate(desRans)
        model_dir_sub = joinpath(model_dir,caseKey,caseKey*"$(desRange)_")
        data = dataset[idxCase][idxDesRan]
        dataR1R2R3 = datasetR1R2R3[idxCase][idxDesRan]
        idxLast = 0
        for (idxOffDesRange, offRange) in enumerate(offDesRans)
            model_name = caseKey*"$(desRange)_"*"$(offRange).jld2"
            model_path = glob(model_name, model_dir_sub)
            if length(model_path)!=1
                continue
            end
            model_path = only(model_path)
            idxLast = idxOffDesRange
            ac = quickload_aircraft(model_path)
            extract_acModel_compact!(ac, data, idxOffDesRange)
        end
        for v in values(data)
            resize!(v,min(length(v),idxLast))
        end
        ####Extract R1 (Assume R1 R2 R3 all exist)
        model_name = caseKey*"$(desRange)_"*"R1_*.jld2"
        model_path = only(glob(model_name, model_dir_sub))
        ac = quickload_aircraft(model_path)
        extract_acModel_compact!(ac, dataR1R2R3, 1)
        #R2
        model_name = caseKey*"$(desRange)_"*"R2_*.jld2"
        model_path = only(glob(model_name, model_dir_sub))
        ac = quickload_aircraft(model_path)
        extract_acModel_compact!(ac, dataR1R2R3, 2)
        #R3
        model_name = caseKey*"$(desRange)_"*"R3_*.jld2"
        model_path = only(glob(model_name, model_dir_sub))
        ac = quickload_aircraft(model_path)
        extract_acModel_compact!(ac, dataR1R2R3, 3)
    end
end

#### Plotting
# Used for plotting
const LINESTYLES = [:solid, :dash, :dot, :dashdot, :dashdotdot]
const LINECOLORS = [:blue, :red, :green, :orange, :purple, :brown, :black, :magenta, :cyan, :olive]
const MARKERS    = [:rect, :circle, :diamond, :utriangle, :dtriangle]
# PRD
p = plot(xlabel="Off-design Ranges [nmi]", ylabel="Payload Weight [Ton]", dpi=800, legend=:best)
for (idxCase,dataCase) in enumerate(dataset)
    for (idxDesRan,dataDesign) in enumerate(dataCase)
        # Setup style
        i = (idxCase-1)*length(dataCase)+idxDesRan
        color = LINECOLORS[mod1(i, length(LINECOLORS))]
        marker = MARKERS[mod1(cld(i, length(LINECOLORS)), length(MARKERS))]
        linestyle = LINESTYLES[mod1(cld(i, length(LINECOLORS) * length(MARKERS)), length(LINESTYLES))]
        # Plot
        plot!(p, dataDesign[:(parm[imRange,2])] ./ 1852.0,
                 dataDesign[:(parm[imWpay,2])] ./ 9.81 ./ 1000.0, 
                 marker=marker, mc=color, msc=color, color=color, lw=2, 
                 linestyle=linestyle, label="$(caseNams[idxCase]) R1: $(desRans[idxDesRan]) nmi")
    end
end
for (idxCase,dataCase) in enumerate(dataset)
    for (idxDesRan,dataDesign) in enumerate(dataCase)
        # Plot R1 R2 R3
        scatter!(p, datasetR1R2R3[idxCase][idxDesRan][:(parm[imRange,2])] ./ 1852.0,
                    datasetR1R2R3[idxCase][idxDesRan][:(parm[imWpay,2])] ./ 9.81 ./ 1000.0, 
                    markershape=:circle, markersize=3, markercolor=:black, markerstrokecolor=:black, 
                    markerstrokewidth=0, label=nothing)
    end
end
xlims!(p, 0, 6000)
ylims!(p, 0, 25)
savefig(p,joinpath(save_dir_sub,"PRD.png"))

# PRD
p = plot(xlabel="Off-design Ranges [nmi]", ylabel="PFEI [J/J]", dpi=800, legend=:best)
for (idxCase,dataCase) in enumerate(dataset)
    for (idxDesRan,dataDesign) in enumerate(dataCase)
        # Setup style
        i = (idxCase-1)*length(dataCase)+idxDesRan
        color = LINECOLORS[mod1(i, length(LINECOLORS))]
        marker = MARKERS[mod1(cld(i, length(LINECOLORS)), length(MARKERS))]
        linestyle = LINESTYLES[mod1(cld(i, length(LINECOLORS) * length(MARKERS)), length(LINESTYLES))]
        # Plot
        plot!(p, dataDesign[:(parm[imRange,2])] ./ 1852.0,
                 dataDesign[:(parm[imPFEI,2])], 
                 marker=marker, mc=color, msc=color, color=color, lw=2, 
                 linestyle=linestyle, label="$(caseNams[idxCase]) R1: $(desRans[idxDesRan]) nmi")
    end
end
xlims!(p, 0, 6000)
ylims!(p, 0.6, 1)
savefig(p,joinpath(save_dir_sub,"PFEI.png"))

# for i in eachindex(xdata)
#         length(xdata[i])==length(ydata[i]) || error("xy data for case $(i) dont have the same number of data points $(length(xdata[i])),$(length(ydata[i]))")
#         color = LINECOLORS[mod1(i, length(LINECOLORS))]
#         marker = MARKERS[mod1(cld(i, length(LINECOLORS)), length(MARKERS))]
#         linestyle = LINESTYLES[mod1(cld(i, length(LINECOLORS) * length(MARKERS)), length(LINESTYLES))]
#         plot!(p, xdata[i], ydata[i], marker=marker, mc=color, msc=color, color=color, lw=lw, linestyle=linestyle, label=datalab[i])
#     end
    


# plot_cases_specified("Off-design Ranges [nmi]", "Payload Weight [Ton]",
#                     [dataDesign[:(parm[imRange,2])] ./ 1852.0 for dataCase in dataset for dataDesign in dataCase],
#                     [dataDesign[:(parm[imWpay,2])] ./ 9.81 ./ 1000.0 for dataCase in dataset for dataDesign in dataCase],
#                     [nothing for caseCur in caseKeys for desRanCur in desRans],
#                      joinpath(save_dir_sub,"PRD.png"))

#                      #"$(caseCur) & R1: $(desRanCur) nmi"