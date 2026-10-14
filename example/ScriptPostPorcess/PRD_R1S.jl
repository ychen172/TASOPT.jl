"""
This script compares payload-range performance of ethanol and jet fuel aircraft from R1 sizing process
    - For exemplary aircraft, plot out payload range performance comparison
"""

using Glob
using TASOPT
include(__TASOPTindices__)
using Plots

#### Setup IO
# Input case names - Retrofit
model_dir  = joinpath(__TASOPTroot__,"../example/ModelSaved")
save_dir   = joinpath(__TASOPTroot__,"../example/OutputAndPlots")
caseKeys   = ["Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_PRD_", "Opti_Eth_NoACT_V4_3_R1Sz_EtasEng_PRD_"]
caseNams   = ["Jet Fuel"                             , "Ethanol"                              ]
desRans    = [1000,1700,2400]
save_prefi = "R1_Size_PRD"

#### Create save directory if not exist
mkpath(save_dir)

#### Pre-allocate per-(case, design-range) result containers
Range_per_case = [[Float64[] for _ in desRans] for _ in caseKeys]
Wpay_per_case  = [[Float64[] for _ in desRans] for _ in caseKeys]

#### Extract data for the sized missions
for (idxCase, caseKey) in enumerate(caseKeys)
    for (idxDesRan, desRange) in enumerate(desRans)
        # PRD folder
        model_dir_sub = joinpath(model_dir,caseKey,caseKey*"$(desRange)_")

        # Directly discover available off-design ranges from the saved filenames
        # (e.g. "..._PRD_400_650.jld2" -> off-design range 650), rather than probing a wide
        # fixed candidate list. Anchored to purely-numeric filenames so R1/R2/R3 files are excluded.
        offDesRan_pattern = Regex("^" * caseKey * "$(desRange)_" * "(\\d+)\\.jld2\$")
        offDesRans_avail = Int[]
        for f in glob(caseKey*"$(desRange)_"*"*.jld2", model_dir_sub)
            m = match(offDesRan_pattern, basename(f))
            if m !== nothing
                push!(offDesRans_avail, parse(Int, m.captures[1]))
            end
        end
        sort!(offDesRans_avail)

        # Collect all off-design points (regular off-design ranges plus the R1/R2/R3 special
        # points) together, then sort by the actual flown range so R1/R2/R3 fold seamlessly
        # into the same payload-range curve instead of being plotted as a separate overlay.
        points = Tuple{Float64,Float64}[] # (Range [m], Wpay [N])
        for offRange in offDesRans_avail
            model_name = caseKey*"$(desRange)_"*"$(offRange).jld2"
            ac = quickload_aircraft(joinpath(model_dir_sub, model_name))
            push!(points, (ac.parm[imRange,2], ac.parm[imWpay,2]))
        end
        for Rtag in ("R1","R2","R3")
            model_name = caseKey*"$(desRange)_"*"$(Rtag)_*.jld2"
            ac = quickload_aircraft(only(glob(model_name, model_dir_sub)))
            push!(points, (ac.parm[imRange,2], ac.parm[imWpay,2]))
        end
        sort!(points, by=first)

        append!(Range_per_case[idxCase][idxDesRan], first.(points))
        append!(Wpay_per_case[idxCase][idxDesRan], last.(points))
    end
end

#### Plotting
# Used for plotting
const LINESTYLES = [:solid, :dash, :dot, :dashdot, :dashdotdot]
const LINECOLORS = [:blue, :red, :green, :orange, :purple, :brown, :black, :magenta, :cyan, :olive]
const MARKERS    = [:rect, :circle, :diamond, :utriangle, :dtriangle]

# Flatten the (case, design-range) grid into parallel lists, one entry per curve
case_combos = [(idxCase, idxDesRan) for idxCase in eachindex(caseKeys) for idxDesRan in eachindex(desRans)]
labels = ["$(caseNams[idxCase]) R1: $(desRans[idxDesRan]) nmi" for (idxCase, idxDesRan) in case_combos]

# PRD
p = plot(xlabel="Off-design Ranges [nmi]", ylabel="Payload Weight [Ton]", dpi=800, legend=:best)
for (i, (idxCase, idxDesRan)) in enumerate(case_combos)
    color = LINECOLORS[mod1(i, length(LINECOLORS))]
    marker = MARKERS[mod1(cld(i, length(LINECOLORS)), length(MARKERS))]
    linestyle = LINESTYLES[mod1(cld(i, length(LINECOLORS) * length(MARKERS)), length(LINESTYLES))]
    plot!(p, Range_per_case[idxCase][idxDesRan] ./ 1852.0,
             Wpay_per_case[idxCase][idxDesRan] ./ 9.81 ./ 1000.0,
             marker=marker, mc=color, msc=color, color=color, lw=2,
             linestyle=linestyle, label=labels[i])
end
xlims!(p, 0, 6000)
ylims!(p, 0, 25)
savefig(p,joinpath(save_dir,"$(save_prefi)_PRD.png"))