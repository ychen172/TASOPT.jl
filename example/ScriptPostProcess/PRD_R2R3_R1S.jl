"""
This script compares the payload-range diagram's R2 (max-fuel range) and R3 (ferry range) points
between ethanol and jet fuel aircraft, across the full R1 (design range) sweep
    - For every common R1 sizing point, extract each fuel's R2 and R3 flown range
    - Plot Ethanol's percentage reduction in R2 and R3 (vs Jet, same R1) as a function of R1
"""

using Glob
using TASOPT
include(__TASOPTindices__)
using Plots

#### Setup IO
model_dir  = joinpath(__TASOPTroot__,"../example/ModelSaved")
save_dir   = joinpath(__TASOPTroot__,"../example/OutputAndPlots")
caseKeys   = ["Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_PRD_", "Opti_Eth_NoACT_V4_3_R1Sz_EtasEng_PRD_"]
caseNams   = ["Jet Fuel"                             , "Ethanol"                              ]
save_prefi = "R1_Size_PRD_R2R3"

#### Create save directory if not exist
mkpath(save_dir)

#### Helper: discover available R1 (design range) subfolders for a case, directly from disk
function available_R1_ranges(caseKey)
    case_dir = joinpath(model_dir, caseKey)
    pattern = Regex("^" * caseKey * "(\\d+)_\$")
    rans = Int[]
    for d in readdir(case_dir)
        isdir(joinpath(case_dir, d)) || continue
        m = match(pattern, d)
        if m !== nothing
            push!(rans, parse(Int, m.captures[1]))
        end
    end
    sort!(rans)
    return rans
end

#### Helper: load the R2 or R3 special point for a given case/R1 and return its flown range [nmi]
function load_PRD_point_range_nmi(caseKey, R1_range, Rtag)
    sub_dir = joinpath(model_dir, caseKey, caseKey*"$(R1_range)_")
    pattern = caseKey*"$(R1_range)_"*"$(Rtag)_*.jld2"
    matches = glob(pattern, sub_dir)
    length(matches) == 1 || error("Expected exactly 1 match for $(pattern) in $(sub_dir), found $(length(matches))")
    ac = quickload_aircraft(matches[1])
    return ac.parm[imRange,2] / 1852.0 #[nmi], off-design (2nd) mission slot
end

#### Discover R1 ranges available for both fuels, keep only the common ones
R1_avail = [available_R1_ranges(k) for k in caseKeys]
println("Available R1 ranges per case:")
for (nam,avail) in zip(caseNams,R1_avail)
    println("  $(nam): $(length(avail)) ranges")
end
R1_common = sort(intersect(R1_avail...))
println("Common R1 ranges: $(length(R1_common))")

#### Extract R2 and R3 for both fuels at every common R1 point
jet_idx = findfirst(==("Jet Fuel"), caseNams)
eth_idx = findfirst(==("Ethanol"), caseNams)

R2_jet = Float64[]; R2_eth = Float64[]
R3_jet = Float64[]; R3_eth = Float64[]
for R1 in R1_common
    push!(R2_jet, load_PRD_point_range_nmi(caseKeys[jet_idx], R1, "R2"))
    push!(R2_eth, load_PRD_point_range_nmi(caseKeys[eth_idx], R1, "R2"))
    push!(R3_jet, load_PRD_point_range_nmi(caseKeys[jet_idx], R1, "R3"))
    push!(R3_eth, load_PRD_point_range_nmi(caseKeys[eth_idx], R1, "R3"))
    println("R1=$(R1) nmi: Jet R2=$(round(R2_jet[end],digits=0)), R3=$(round(R3_jet[end],digits=0))  ||  Eth R2=$(round(R2_eth[end],digits=0)), R3=$(round(R3_eth[end],digits=0))")
end

#### Percentage reduction of Ethanol's R2/R3 vs Jet's, at the same R1: (Jet-Eth)/Jet*100
pct_R2_reduction = 100.0 .* (R2_jet .- R2_eth) ./ R2_jet
pct_R3_reduction = 100.0 .* (R3_jet .- R3_eth) ./ R3_jet

#### Design-range onset of Ethanol's fan-diameter clipping (soft upper-bound constraint in the R1Size
#### optimization: parg[igdfan] vs 2.0m fan-diameter limit) - matches DFAN_CLIP_RANGE in the other R1S scripts
# const DFAN_CLIP_RANGE = 2700.0 #nmi

p = plot(R1_common, pct_R2_reduction, label="R2 (Fuel Tank Full)", color=:royalblue, lw=2.5, marker=:circle, markersize=3,
         xlabel="Design/R1 Range [nmi]", ylabel="Range Reduction [%]", dpi=800,
         legend=:topright, size=(1400, 900), titlefontsize=21, guidefontsize=24, tickfontsize=18, legendfontsize=18,
         left_margin=16Plots.mm, bottom_margin=12Plots.mm)
plot!(p, R1_common, pct_R3_reduction, label="R3 (One Passenger)", color=:crimson, lw=2.5, marker=:circle, markersize=3)
hline!(p, [0.0], line=(:dash, 1, :black), label="")
# vline!(p, [DFAN_CLIP_RANGE], line=(:dash, 2, :gray40), label="Max DiaFan Lim (Eth)")
savefig(p, joinpath(save_dir, "$(save_prefi)_PctReduction_vs_R1.png"))
