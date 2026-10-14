"""
This script compares and plots the empty weights between ethanol and jet fuel aircraft from R1 sizing process
    - Compare empty weight, fuel weight, passenger weight
    - Compare empty weight breakdown over engine weight, wing weight, fuselage weight, and others(hydraulic, landing gear, and tails)
"""

using TASOPT
include(__TASOPTindices__)
using Plots, StatsPlots
using Glob

#### Setup IO
# Input case names - Retrofit
model_dir  = joinpath(__TASOPTroot__,"../example/ModelSaved")
save_dir   = joinpath(__TASOPTroot__,"../example/OutputAndPlots")
case_keys  = ["Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_","Opti_Eth_NoACT_V4_3_R1Sz_EtasEng_"]
case_Nams  = ["Jet Fuel"                         ,"Ethanol"                          ]
save_prefi = "R1_Size_Fli_Weights"

#### Create save directory if not exist
mkpath(save_dir)

#### Pre-allocate per-case result containers (one vector per case, filled over available ranges), all in metric tons
rans_per_case    = [Int[] for _ in case_keys]
Wpay_per_case    = [Float64[] for _ in case_keys]
Wfuel_per_case   = [Float64[] for _ in case_keys]
Wempty_per_case  = [Float64[] for _ in case_keys]
Wfuse_per_case   = [Float64[] for _ in case_keys]
Wwing_per_case   = [Float64[] for _ in case_keys]
Weng_per_case    = [Float64[] for _ in case_keys]
Wothers_per_case = [Float64[] for _ in case_keys]
WhydLG_per_case  = [Float64[] for _ in case_keys] #Hydraulic + Nose Gear + Main Gear combined
Wtails_per_case  = [Float64[] for _ in case_keys] #Htail + Vtail combined
Minn_per_case    = [Float64[] for _ in case_keys] #Wing inboard root bending moment [kN*m]
WingS_per_case   = [Float64[] for _ in case_keys] #Wing planform area [m^2]
RootChord_per_case = [Float64[] for _ in case_keys] #Wing root chord [m]

#### Loop through each case to extract weights over ranges
for (j, key_j) in enumerate(case_keys)
    # Idenfity available ranges within each case
    case_dir = joinpath(model_dir,key_j)
    pattern = Regex("^" * key_j * "(\\d+)\\.jld2\$") #Match "<key><range>.jld2" only, exclude "..._optimization_history.jld2"
    rans_avail = Int[]
    for f in glob(key_j*"*.jld2", case_dir)
        m = match(pattern, basename(f))
        if m !== nothing
            push!(rans_avail, parse(Int, m.captures[1]))
        end
    end
    sort!(rans_avail)
    println("$(key_j): found $(length(rans_avail)) available design ranges")

    #### Loop through each range
    for (i, ran_i) in enumerate(rans_avail)
        # Read in aircraft model
        ac = quickload_aircraft(joinpath(model_dir,key_j,key_j*"$(round(Int,ran_i)).jld2"))
        println("$(key_j)*$(ran_i) loaded")
        # Extract aircraft total takeoff weight, initial fuel weight, passenger weight, and emtpy weights
        WTO   = ac.parm[imWTO,1]
        Wfuel = ac.parm[imWfuel,1]
        Wpay  = ac.parm[imWpay,1]
        Wempty = WTO - Wfuel - Wpay #Matches TASOPT.plot_details' own definition

        # Extract empty weight breakdown, engine weight, wing weight, fuselage weight, and others(hydraulic, landing gear, and tails)
        Wfuse   = ac.fuselage.weight
        Wwing   = ac.wing.weight
        Weng    = ac.parg[igWeng]
        Whpesys = ac.parg[igWMTO] * ac.fuselage.HPE_sys.W
        Wlgnose = ac.landing_gear.nose_gear.weight.W
        Wlgmain = ac.landing_gear.main_gear.weight.W
        Whtail  = ac.htail.weight
        Wvtail  = ac.vtail.weight
        Wftank  = ac.parg[igWftank]
        Wtesys  = ac.parg[igWtesys]
        Wothers = Whpesys + Wlgnose + Wlgmain + Whtail + Wvtail + Wftank + Wtesys #Everything in Wempty besides fuselage/wing/engine
        WhydLG  = Whpesys + Wlgnose + Wlgmain
        Wtails  = Whtail + Wvtail
        Minn    = ac.wing.inboard.moment
        WingS   = ac.wing.layout.S
        RootChord = ac.wing.layout.root_chord

        # Process and collect weights, converting from force [N] to mass [metric ton]
        push!(rans_per_case[j], ran_i)
        push!(Wpay_per_case[j], Wpay/TASOPT.gee/1000.0)
        push!(Wfuel_per_case[j], Wfuel/TASOPT.gee/1000.0)
        push!(Wempty_per_case[j], Wempty/TASOPT.gee/1000.0)
        push!(Wfuse_per_case[j], Wfuse/TASOPT.gee/1000.0)
        push!(Wwing_per_case[j], Wwing/TASOPT.gee/1000.0)
        push!(Weng_per_case[j], Weng/TASOPT.gee/1000.0)
        push!(Wothers_per_case[j], Wothers/TASOPT.gee/1000.0)
        push!(WhydLG_per_case[j], WhydLG/TASOPT.gee/1000.0)
        push!(Wtails_per_case[j], Wtails/TASOPT.gee/1000.0)
        push!(Minn_per_case[j], Minn/1000.0)
        push!(WingS_per_case[j], WingS)
        push!(RootChord_per_case[j], RootChord)
    end

end

#### Plots
#### Grouped + stacked + hatched bar chart: one bar per case at every range (Jet solid, Ethanol hatched),
#### each bar internally stacked by weight component. Each component gets its own hatch style (applied to
#### the Ethanol/hatched bars only) in addition to its own color, for easier segment-by-segment distinction.
const HATCH_STYLES = [:/, :x, :+, :-] #Distinct, GR-valid hatch styles; :., :*, :o, :0, and :\\ were tested and unreliable

#### Single shared color sequence, applied positionally bottom-up in every stacked bar/line plot below
#### (taken from the original Empty Weight breakdown bar, which used Fuselage=teal/Wing=royalblue/
#### Engine=crimson/Others=purple bottom-to-top), so color always means "stack position" consistently
#### across all four plots, not a per-plot arbitrary palette.
const BASE_COLORS = [:teal, :royalblue, :crimson, :purple]

#### Design-range onset of Ethanol's fan-diameter clipping (soft upper-bound Constraint in the R1Size
#### optimization: parg[igdfan] vs 2.0m fan-diameter limit). Jet never reaches this bound in the sweep
#### (maxes out at 91.2% of the limit at 3000nmi), so past this range the comparison is no longer
#### apples-to-apples: Ethanol's engine is pinned against the bound while Jet's keeps growing freely.
const DFAN_CLIP_RANGE = 2700.0 #nmi; Ethanol's fan diameter first reaches 2.0m here

#### Builds one stacked groupedbar call per side (Jet, Ethanol) using the full component matrix at once,
#### so StatsPlots handles the stacking math itself (proven reliable, unlike manual bottom-tracking).
function plot_stacked_comparison_bar(ranges::Vector{Float64}, jet_components, eth_components, comp_labels, colors,
                                      ylabel_str, title_str, save_path; bar_width=85.0, gap=10.0)
    x_jet = ranges .- (bar_width/2 + gap/2)
    x_eth = ranges .+ (bar_width/2 + gap/2)
    # groupedbar(...; bar_position=:stack) stacks the FIRST column on top and the LAST column at the
    # bottom, so the matrix/color/style columns are fed in reverse of comp_labels to render bottom-up
    # in the order comp_labels was given (e.g. Payload->Empty->Fuel bottom to top).
    mat_jet = hcat(reverse(jet_components)...)
    mat_eth = hcat(reverse(eth_components)...)
    ymax = max(maximum(sum(jet_components)), maximum(sum(eth_components)))
    colors_row = reshape(reverse(colors), 1, :)
    nk = length(comp_labels)
    styles = HATCH_STYLES[mod1.(1:nk, length(HATCH_STYLES))]
    styles_row = reshape(reverse(styles), 1, :)

    # The real stacking series carry no legend of their own (label=""); a 2-col legend fills row-major
    # in series-add order, so two full-matrix groupedbar calls can't produce "all Jet in col 1, all Eth
    # in col 2" directly. Instead, draw invisible NaN dummy series afterward, ordered Jet1,Eth1,Jet2,Eth2,...
    # so row-major fill puts every Jet entry in column 1 and every Eth entry in column 2, each with the
    # real color+hatch swatch matching its stacked segment.
    p = groupedbar(x_jet, mat_jet, bar_position=:stack, bar_width=bar_width, color=colors_row, linecolor=:black,
                    label="", xlabel="Design/R1 Range [nmi]", ylabel=ylabel_str, dpi=800,
                    legend=:topleft, legend_column=2, ylims=(0.0, 1.45*ymax), size=(1400, 900), titlefontsize=21, guidefontsize=24, tickfontsize=18, legendfontsize=18, left_margin=16Plots.mm, bottom_margin=12Plots.mm)
    groupedbar!(p, x_eth, mat_eth, bar_position=:stack, bar_width=bar_width, color=colors_row, linecolor=:black,
                fillstyle=styles_row, label="")

    for i in 1:nk
        bar!(p, [NaN], [NaN], color=colors[i], linecolor=:black, label=comp_labels[i]*" (Jet)")
        bar!(p, [NaN], [NaN], color=colors[i], linecolor=:black, fillstyle=styles[i], label=comp_labels[i]*" (Eth)")
    end

    # Mark where the fan-diameter constraint starts clipping, to flag where the Jet-vs-Ethanol
    # comparison is no longer apples-to-apples (see DFAN_CLIP_RANGE definition above)
    vline!(p, [DFAN_CLIP_RANGE], line=(:dash, 2, :gray40), label="Max DiaFan Lim (Eth)")

    savefig(p, save_path)
    return p
end

jet_idx = findfirst(==("Jet Fuel"), case_Nams)
eth_idx = findfirst(==("Ethanol"), case_Nams)

#### Only compare ranges that actually converged for BOTH fuels, then keep every 2nd common range so bars can be wider
ranges_common_all = sort(intersect(rans_per_case[jet_idx], rans_per_case[eth_idx]))
ranges_common = Float64.(ranges_common_all[1:2:end])
println("Common ranges available: $(length(ranges_common_all)); plotting every 2nd => $(length(ranges_common)) bars")

#### Index lookup per case: for each common range, find its position within that case's own (possibly
#### differently-ordered/sized) result arrays, rather than assuming the two cases line up index-for-index
jet_pos = [findfirst(==(round(Int,r)), rans_per_case[jet_idx]) for r in ranges_common]
eth_pos = [findfirst(==(round(Int,r)), rans_per_case[eth_idx]) for r in ranges_common]

plot_stacked_comparison_bar(ranges_common,
    [Wpay_per_case[jet_idx][jet_pos], Wempty_per_case[jet_idx][jet_pos], Wfuel_per_case[jet_idx][jet_pos]],
    [Wpay_per_case[eth_idx][eth_pos], Wempty_per_case[eth_idx][eth_pos], Wfuel_per_case[eth_idx][eth_pos]],
    ["Payload", "Empty", "Fuel"], BASE_COLORS[1:3],
    "Weight [metric ton]", "Total Weight Breakdown: Jet (solid) vs Ethanol (hatched)",
    joinpath(save_dir, "$(save_prefi)_TotalWeight_BarCompare.png"))

plot_stacked_comparison_bar(ranges_common,
    [Wfuse_per_case[jet_idx][jet_pos], Wwing_per_case[jet_idx][jet_pos], Weng_per_case[jet_idx][jet_pos], Wothers_per_case[jet_idx][jet_pos]],
    [Wfuse_per_case[eth_idx][eth_pos], Wwing_per_case[eth_idx][eth_pos], Weng_per_case[eth_idx][eth_pos], Wothers_per_case[eth_idx][eth_pos]],
    ["Fuselage", "Wing", "Engine", "Others"], BASE_COLORS,
    "Weight [metric ton]", "Empty Weight Breakdown: Jet (solid) vs Ethanol (hatched)",
    joinpath(save_dir, "$(save_prefi)_EmptyWeight_BarCompare.png"))

#### Breakdown of "Others" into 2 groups: Hydraulics & Landing Gear combined (all fixed WMTO mass
#### fractions, so they move together and aren't independently meaningful - see discussion above;
#### Hydraulic alone is too small to warrant its own line/segment), and Tails (Htail+Vtail combined),
#### the genuinely independent structural driver. (Wftank/Wtesys are always 0 for these wing-fuel
#### cases, so they're omitted; both groups are pre-combined in the extraction loop above.)
plot_stacked_comparison_bar(ranges_common,
    [WhydLG_per_case[jet_idx][jet_pos], Wtails_per_case[jet_idx][jet_pos]],
    [WhydLG_per_case[eth_idx][eth_pos], Wtails_per_case[eth_idx][eth_pos]],
    ["Hydraulics & Landing Gear", "Vertical & Horizontal Tails"], BASE_COLORS[1:2],
    "Weight [metric ton]", "Others Breakdown: Jet (solid) vs Ethanol (hatched)",
    joinpath(save_dir, "$(save_prefi)_Others_BarCompare.png"))

#### Relative weight increase of Ethanol vs Jet, (Eth-Jet)/Jet*100, for each of the 4 empty-weight
#### components, across the FULL common-range sweep (not just every 2nd point, since this is a line
#### plot, not bars). As design range grows, required thrust grows too, and Engine/Others pull away
#### from Fuselage/Wing in relative terms, showing them as the main drivers of Ethanol's empty-weight penalty.
jet_pos_all = [findfirst(==(r), rans_per_case[jet_idx]) for r in ranges_common_all]
eth_pos_all = [findfirst(==(r), rans_per_case[eth_idx]) for r in ranges_common_all]

pct_diff(jet_vec, eth_vec) = [(eth_vec[eth_pos_all[k]] - jet_vec[jet_pos_all[k]]) / jet_vec[jet_pos_all[k]] * 100.0
                              for k in eachindex(ranges_common_all)]

pct_fuse   = pct_diff(Wfuse_per_case[jet_idx], Wfuse_per_case[eth_idx])
pct_wing   = pct_diff(Wwing_per_case[jet_idx], Wwing_per_case[eth_idx])
pct_eng    = pct_diff(Weng_per_case[jet_idx], Weng_per_case[eth_idx])
pct_others = pct_diff(Wothers_per_case[jet_idx], Wothers_per_case[eth_idx])

#### Colors match BASE_COLORS positionally - these 4 components are exactly EmptyWeight_BarCompare's
#### breakdown, so the line plot's colors are drawn straight from that "corresponding breakup bar plot".
p = plot(ranges_common_all, pct_fuse, label="Fuselage", color=BASE_COLORS[1], lw=2.5, marker=:circle, markersize=3,
         xlabel="Design/R1 Range [nmi]", ylabel="Ethanol vs Jet Weight Difference [%]", dpi=800,
         legend=:topleft, size=(1400, 900), titlefontsize=21, guidefontsize=24, tickfontsize=18, legendfontsize=18, left_margin=16Plots.mm, bottom_margin=12Plots.mm)
plot!(p, ranges_common_all, pct_wing, label="Wing", color=BASE_COLORS[2], lw=2.5, marker=:circle, markersize=3)
plot!(p, ranges_common_all, pct_eng, label="Engine", color=BASE_COLORS[3], lw=2.5, marker=:circle, markersize=3)
plot!(p, ranges_common_all, pct_others, label="Others", color=BASE_COLORS[4], lw=2.5, marker=:circle, markersize=3)
hline!(p, [0.0], line=(:dash, 1, :black), label="")
vline!(p, [DFAN_CLIP_RANGE], line=(:dash, 2, :gray40), label="Max DiaFan Lim (Eth)")
savefig(p, joinpath(save_dir, "$(save_prefi)_ComponentPctDiff_vs_Range.png"))

#### Same relative-increase line plot, but for the "Others" breakdown (Hydraulics & Landing Gear, Vertical & Horizontal Tails)
#### - colors match BASE_COLORS[1:2] positionally, drawn from Others_BarCompare's own palette.
pct_hydLG = pct_diff(WhydLG_per_case[jet_idx], WhydLG_per_case[eth_idx])
pct_tails = pct_diff(Wtails_per_case[jet_idx], Wtails_per_case[eth_idx])

p2 = plot(ranges_common_all, pct_hydLG, label="Hydraulics & Landing Gear", color=BASE_COLORS[1], lw=2.5, marker=:circle, markersize=3,
          xlabel="Design/R1 Range [nmi]", ylabel="Ethanol vs Jet Weight Difference [%]", dpi=800,
          legend=:topleft, size=(1400, 900), titlefontsize=21, guidefontsize=24, tickfontsize=18, legendfontsize=18, left_margin=16Plots.mm, bottom_margin=12Plots.mm)
plot!(p2, ranges_common_all, pct_tails, label="Vertical & Horizontal Tails", color=BASE_COLORS[2], lw=2.5, marker=:circle, markersize=3)
hline!(p2, [0.0], line=(:dash, 1, :black), label="")
vline!(p2, [DFAN_CLIP_RANGE], line=(:dash, 2, :gray40), label="Max DiaFan Lim (Eth)")
savefig(p2, joinpath(save_dir, "$(save_prefi)_Others_ComponentPctDiff_vs_Range.png"))

#### Wing inboard (root) bending moment (the load actually used to size the inboard spar, see
#### wing_weights.jl) across the range sweep, Jet vs Ethanol. Outboard omitted - inboard alone is
#### enough to show the trend, and the outboard curve tracked it closely anyway.
Minn_jet_all = [Minn_per_case[jet_idx][jet_pos_all[k]] for k in eachindex(ranges_common_all)]
Minn_eth_all = [Minn_per_case[eth_idx][eth_pos_all[k]] for k in eachindex(ranges_common_all)]

p3 = plot(ranges_common_all, Minn_jet_all, label="Jet", color=BASE_COLORS[1], lw=2.5, marker=:circle, markersize=3,
          xlabel="Design/R1 Range [nmi]", ylabel="Bending Moment [kN*m]", dpi=800,
          legend=:topleft, size=(1400, 900), titlefontsize=21, guidefontsize=24, tickfontsize=18, legendfontsize=18, left_margin=16Plots.mm, bottom_margin=12Plots.mm)
plot!(p3, ranges_common_all, Minn_eth_all, label="Ethanol", color=BASE_COLORS[2], lw=2.5, marker=:circle, markersize=3)
vline!(p3, [DFAN_CLIP_RANGE], line=(:dash, 2, :gray40), label="Max DiaFan Lim (Eth)")
savefig(p3, joinpath(save_dir, "$(save_prefi)_WingBendingMoment_vs_Range.png"))

#### Wing root chord, Jet vs Ethanol, across the range sweep - the geometric quantity that drives the
#### c^3 structural-efficiency effect discussed above (root chord feeds the inboard spar-cap sizing).
RootChord_jet_all = [RootChord_per_case[jet_idx][jet_pos_all[k]] for k in eachindex(ranges_common_all)]
RootChord_eth_all = [RootChord_per_case[eth_idx][eth_pos_all[k]] for k in eachindex(ranges_common_all)]

p3b = plot(ranges_common_all, RootChord_jet_all, label="Jet", color=BASE_COLORS[1], lw=2.5, marker=:circle, markersize=3,
           xlabel="Design/R1 Range [nmi]", ylabel="Root Chord [m]", dpi=800,
           legend=:topleft, size=(1400, 900), titlefontsize=21, guidefontsize=24, tickfontsize=18, legendfontsize=18, left_margin=16Plots.mm, bottom_margin=12Plots.mm)
plot!(p3b, ranges_common_all, RootChord_eth_all, label="Ethanol", color=BASE_COLORS[2], lw=2.5, marker=:circle, markersize=3)
vline!(p3b, [DFAN_CLIP_RANGE], line=(:dash, 2, :gray40), label="Max DiaFan Lim (Eth)")
savefig(p3b, joinpath(save_dir, "$(save_prefi)_WingRootChord_vs_Range.png"))

#### Wing "density" = wing structural weight / wing planform area, Jet vs Ethanol, across the range sweep
WingDensity_jet = [1000.0 * Wwing_per_case[jet_idx][jet_pos_all[k]] / WingS_per_case[jet_idx][jet_pos_all[k]]
                   for k in eachindex(ranges_common_all)] #[kg/m^2]
WingDensity_eth = [1000.0 * Wwing_per_case[eth_idx][eth_pos_all[k]] / WingS_per_case[eth_idx][eth_pos_all[k]]
                   for k in eachindex(ranges_common_all)] #[kg/m^2]

p4 = plot(ranges_common_all, WingDensity_jet, label="Jet", color=BASE_COLORS[1], lw=2.5, marker=:circle, markersize=3,
          xlabel="Design/R1 Range [nmi]", ylabel="Wing Mass Density [kg/m^2]", dpi=800,
          legend=:bottomleft, size=(1400, 900), titlefontsize=21, guidefontsize=24, tickfontsize=18, legendfontsize=18, left_margin=16Plots.mm, bottom_margin=12Plots.mm)
plot!(p4, ranges_common_all, WingDensity_eth, label="Ethanol", color=BASE_COLORS[2], lw=2.5, marker=:circle, markersize=3)
vline!(p4, [DFAN_CLIP_RANGE], line=(:dash, 2, :gray40), label="Max DiaFan Lim (Eth)")
savefig(p4, joinpath(save_dir, "$(save_prefi)_WingDensity_vs_Range.png"))
