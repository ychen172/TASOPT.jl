"""
This script compares payload-range performance of ethanol and jet fuel aircraft from OAG sizing
process, for a handful of example seat capacities
    - For exemplary aircraft seat capacities, plot out payload range performance comparison
"""

using Glob
using TASOPT
include(__TASOPTindices__)
using Plots
include(joinpath(@__DIR__,"../UtilityFunctions/Extract.jl"))
using .Extract: read_oag

#### Setup IO
model_dir  = joinpath(__TASOPTroot__,"../example/ModelSaved")
save_dir   = joinpath(__TASOPTroot__,"../example/OutputAndPlots")
# PRD top-level folders
case_folders    = ["Opti_Jet_NoACT_OAG_MI_Tail_V3_PRD", "Opti_Eth_NoACT_OAG_MI_Tail_V3_PRD"]
case_name_front = ["Opti_Jet_NoACT_OAG_Ml_*"          , "Opti_Eth_NoACT_OAG_Ml_*"          ] #"Use * for outside of pattern matching"
case_name_back  = ["_PRD","_PRD"] #For subfolder name, after _Seat 
caseNams        = ["Jet Fuel"                         , "Ethanol"                          ]
seat_caps       = [100, 150, 200]
save_prefi      = "OAG_PRD"
# Required OAG mission set data
mission_file    = "OffDesignMissions_50_300_300_Tail.csv"
pass_load_frac_off = 0.825 #Load factor of passengers (lower 4 missions) (fraction of ac.parg[igWpaymax])
pass_load_frac_tail = 0.850 #The farthest mission payload fraction (tail mission)
# Plotting
plt_marker      = [:diamond, :circle, :rect, :utriangle, :dtriangle, :star5, :hexagon, :pentagon, :xcross, :cross] #Apply to seats
plt_linestyle   = [:solid  , :dash] #Apply to cases
plt_linecolor   = [:red    , :purple] #Apply to cases

#### Create save directory if not exist
mkpath(save_dir)

#### Load the required OAG mission set (4 bulk + 1 tail range per seat capacity)
miss_off_des = read_oag(joinpath(model_dir, mission_file))

#### Pre-allocate per-(case, seat-capacity) result containers
Range_per_case = [[Float64[] for _ in seat_caps] for _ in case_folders]
Wpay_per_case  = [[Float64[] for _ in seat_caps] for _ in case_folders]
# Required-mission scatter points (4 bulk + 1 tail), one set per (case, seat capacity)
MissionRange_per_case = [[Float64[] for _ in seat_caps] for _ in case_folders]
MissionWpay_per_case  = [[Float64[] for _ in seat_caps] for _ in case_folders]

#### Extract data for the sized missions
for (idxCase, case_folder) in enumerate(case_folders)
    for (idxSeat, sc) in enumerate(seat_caps)
        # Find this seat capacity's PRD subfolder
        case_dir = joinpath(model_dir, case_folder)
        sub_pattern = case_name_front[idxCase]*"_$(sc)"*case_name_back[idxCase]
        sub_matches = glob(sub_pattern, case_dir)
        length(sub_matches)==1 || error("Expected exactly 1 PRD subfolder match for $(sub_pattern) in $(case_dir), found $(length(sub_matches)): $(sub_matches)")
        model_dir_sub = sub_matches[1]
        prefix = basename(model_dir_sub) #Folder name is the prefix for the file names

        # Directly discover available off-design ranges from the saved filenames
        offDesRan_pattern = Regex("^" * prefix * "(\\d+)\\.jld2\$")
        offDesRans_avail = Int[]
        for f in glob(prefix*"*.jld2", model_dir_sub)
            m = match(offDesRan_pattern, basename(f))
            if m !== nothing
                push!(offDesRans_avail, parse(Int, m.captures[1]))
            end
        end
        sort!(offDesRans_avail)

        # Collect all off-design points
        points = Tuple{Float64,Float64}[] # (Range [m], Wpay [N])
        for offRange in offDesRans_avail
            model_name = prefix*"$(offRange).jld2"
            ac = quickload_aircraft(joinpath(model_dir_sub, model_name))
            push!(points, (ac.parm[imRange,2], ac.parm[imWpay,2]))
        end
        for Rtag in ("R1","R2","R3")
            model_name = prefix*"$(Rtag)_*.jld2"
            ac = quickload_aircraft(only(glob(model_name, model_dir_sub)))
            push!(points, (ac.parm[imRange,2], ac.parm[imWpay,2]))
        end
        sort!(points, by=first)

        append!(Range_per_case[idxCase][idxSeat], first.(points))
        append!(Wpay_per_case[idxCase][idxSeat], last.(points))

        # Required OAG missions for this seat capacity
        Wpaymax = quickload_aircraft(only(glob(prefix*"R1_*.jld2", model_dir_sub))).parg[igWpaymax]
        miss_ranges_nmi = miss_off_des[sc].ranges_nmi
        miss_wpay_N = fill(Wpaymax*pass_load_frac_off, length(miss_ranges_nmi))
        miss_wpay_N[end] = Wpaymax*pass_load_frac_tail
        append!(MissionRange_per_case[idxCase][idxSeat], miss_ranges_nmi)
        append!(MissionWpay_per_case[idxCase][idxSeat], miss_wpay_N)
    end
end

#### Extract the R1 mission (max payload, max range) for every seat capacity with a PRD subfolder
R1_seats_per_case = [Int[]     for _ in case_folders]
R1_range_per_case = [Float64[] for _ in case_folders] # [m]
R1_wpay_per_case  = [Float64[] for _ in case_folders] # [N]
for (idxCase, case_folder) in enumerate(case_folders)
    case_dir = joinpath(model_dir, case_folder)
    seat_pattern = Regex("_(\\d+)" * case_name_back[idxCase] * "\$") #Seat capacity right before the subfolder suffix
    for model_dir_sub in glob(case_name_front[idxCase]*case_name_back[idxCase], case_dir)
        prefix = basename(model_dir_sub)
        m = match(seat_pattern, prefix)
        m === nothing && continue
        ac = quickload_aircraft(only(glob(prefix*"R1_*.jld2", model_dir_sub)))
        push!(R1_seats_per_case[idxCase], parse(Int, m.captures[1]))
        push!(R1_range_per_case[idxCase], ac.parm[imRange,2])
        push!(R1_wpay_per_case[idxCase], ac.parm[imWpay,2])
    end
    # Order by seat capacity so the curve runs from small to large aircraft
    idx_sort = sortperm(R1_seats_per_case[idxCase])
    permute!(R1_seats_per_case[idxCase], idx_sort)
    permute!(R1_range_per_case[idxCase], idx_sort)
    permute!(R1_wpay_per_case[idxCase], idx_sort)
    println("$(caseNams[idxCase]) R1 seat capacities: $(R1_seats_per_case[idxCase])")
end

#### Plotting
# Flatten the (case, seat-capacity) grid into parallel lists, one entry per curve
case_combos = [(idxCase, idxSeat) for idxCase in eachindex(case_folders) for idxSeat in eachindex(seat_caps)]
labels = ["$(caseNams[idxCase]), $(seat_caps[idxSeat]) Seats" for (idxCase, idxSeat) in case_combos]

# PRD
p = plot(xlabel="Off-design Ranges [nmi]", ylabel="Payload Weight [Ton]", dpi=800, legend=:best)
for (i, (idxCase, idxSeat)) in enumerate(case_combos)
    color = plt_linecolor[idxCase]
    marker = plt_marker[idxSeat]
    linestyle = plt_linestyle[idxCase]
    plot!(p, Range_per_case[idxCase][idxSeat] ./ 1852.0,
             Wpay_per_case[idxCase][idxSeat] ./ 9.81 ./ 1000.0,
             marker=marker, mc=color, msc=color, color=color, lw=2,
             linestyle=linestyle, label=labels[i])
end

# Required OAG mission set (4 bulk + 1 tail)
for (i, (idxCase, idxSeat)) in enumerate(case_combos)
    lbl = i==1 ? "Target Missions" : "" #Only the first oag mission scatter get a legend
    scatter!(p, MissionRange_per_case[idxCase][idxSeat],
                MissionWpay_per_case[idxCase][idxSeat] ./ 9.81 ./ 1000.0,
                marker=plt_marker[idxSeat], markercolor=:black, markerstrokecolor=:black,
                label=lbl)
end
# xlims!(p, 0, 6000)
ylims!(p, 0, 30)
savefig(p,joinpath(save_dir,"$(save_prefi)_PRD.png"))

# R1 payload vs range across all seat capacities, one curve per case
p_R1 = plot(xlabel="R1 Range [nmi]", ylabel="R1 Payload Weight [Ton]", dpi=800, legend=:best)
for idxCase in eachindex(case_folders)
    color = plt_linecolor[idxCase]
    plot!(p_R1, R1_range_per_case[idxCase] ./ 1852.0,
                R1_wpay_per_case[idxCase] ./ 9.81 ./ 1000.0,
                marker=:circle, mc=color, msc=color, color=color, lw=2,
                linestyle=plt_linestyle[idxCase], label=caseNams[idxCase])
end
savefig(p_R1,joinpath(save_dir,"$(save_prefi)_R1.png"))
