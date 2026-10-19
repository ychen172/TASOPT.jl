"""
This script compares and plots the aerodynamic/geometric parameters between ethanol and jet fuel aircraft from R1 sizing process
    - Compare wing span, wing aspect ratio, wing surface area, and lift-to-drag ratio at start of cruise
"""

using TASOPT
include(__TASOPTindices__)
using Plots
using Glob

#### Setup IO
# Input case names - Retrofit
model_dir  = joinpath(__TASOPTroot__,"../example/ModelSaved")
save_dir   = joinpath(__TASOPTroot__,"../example/OutputAndPlots")
case_keys  = ["Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_","Opti_Eth_NoACT_V4_3_R1Sz_EtasEng_"]
case_Nams  = ["Jet Fuel"                         ,"Ethanol"                          ]
save_prefi = "R1_Size_Fli_Aero"

#### Create save directory if not exist
mkpath(save_dir)

#### Pre-allocate per-case result containers (one vector per case, filled over available ranges)
rans_per_case = [Int[] for _ in case_keys]
span_per_case = [Float64[] for _ in case_keys]
AR_per_case   = [Float64[] for _ in case_keys]
LD_per_case   = [Float64[] for _ in case_keys]
S_per_case    = [Float64[] for _ in case_keys]
frac_nace_fuse_drag = [Float64[] for _ in case_keys]

#### Loop through each case to extract aero/geometry parameters over ranges
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

        # Extract wing geometry and start-of-cruise aero performance
        span = ac.wing.layout.span
        AR   = ac.wing.layout.AR
        S    = ac.wing.layout.S
        CL   = ac.para[iaCL, ipcruise1, 1]
        CD   = ac.para[iaCD, ipcruise1, 1]

        CDfuse  = ac.para[iaCDfuse, ipcruise1, 1]
        CDi     = ac.para[iaCDi, ipcruise1, 1]
        CDwing  = ac.para[iaCDwing, ipcruise1, 1]
        CDhtail = ac.para[iaCDhtail, ipcruise1, 1]
        CDvtail = ac.para[iaCDvtail, ipcruise1, 1]
        CDnace  = ac.para[iaCDnace, ipcruise1, 1]

        push!(rans_per_case[j], ran_i)
        push!(span_per_case[j], span)
        push!(AR_per_case[j], AR)
        push!(S_per_case[j], S)
        push!(LD_per_case[j], CL/CD)
        push!(frac_nace_fuse_drag[j], (CDi+CDwing)/(CDi+CDwing+CDfuse+CDhtail+CDvtail+CDnace))
    end

end

#### Plots
#### Single shared color sequence, same as empty_weight_R1S.jl's BASE_COLORS, applied positionally
#### (Jet first, Ethanol second) across all plots in this script.
const BASE_COLORS = [:teal, :royalblue]

#### Design-range onset of each fuel's own wing-span clipping (soft upper-bound Constraint in the
#### R1Size optimization: wing.layout.span vs 35.814m Type-C gate limit). Jet clips later (~800nmi)
#### than Ethanol (~500nmi) since Ethanol's heavier fuel load pushes it against the span bound sooner.
const SPAN_CLIP_JET = 800.0 #nmi
const SPAN_CLIP_ETH = 500.0 #nmi

#### Builds one 2-curve comparison line plot (Jet vs Ethanol), styled consistently with empty_weight_R1S.jl.
#### show_span_clip adds two vertical dashed lines marking where each fuel's wing span hits the gate limit.
function plot_two_curves(ranges_per_case, y_per_case, case_Nams, ylabel_str, save_path;
                          legend_pos=:topleft, show_span_clip=false)
    p = plot(ranges_per_case[1], y_per_case[1], label=case_Nams[1], color=BASE_COLORS[1], lw=2.5, marker=:circle, markersize=3,
             xlabel="Design/R1 Range [nmi]", ylabel=ylabel_str, dpi=800,
             legend=legend_pos, size=(1400, 900), titlefontsize=21, guidefontsize=24, tickfontsize=18, legendfontsize=18,
             left_margin=16Plots.mm, bottom_margin=12Plots.mm)
    plot!(p, ranges_per_case[2], y_per_case[2], label=case_Nams[2], color=BASE_COLORS[2], lw=2.5, marker=:circle, markersize=3)
    if show_span_clip
        vline!(p, [SPAN_CLIP_JET], line=(:dash, 2, BASE_COLORS[1]), label="Max Wingspan Lim (Jet)")
        vline!(p, [SPAN_CLIP_ETH], line=(:dash, 2, BASE_COLORS[2]), label="Max Wingspan Lim (Eth)")
    end
    savefig(p, save_path)
    return p
end

plot_two_curves(rans_per_case, AR_per_case, case_Nams, "Wing Aspect Ratio",
                joinpath(save_dir, "$(save_prefi)_AspectRatio.png");
                legend_pos=:topright, show_span_clip=true)

plot_two_curves(rans_per_case, LD_per_case, case_Nams, "Lift-to-Drag Ratio at Cruise",
                joinpath(save_dir, "$(save_prefi)_LD_Cruise.png");
                legend_pos=:bottom, show_span_clip=true)

plot_two_curves(rans_per_case, S_per_case, case_Nams, "Wing Surface Area [m^2]",
                joinpath(save_dir, "$(save_prefi)_WingSurfaceArea.png");
                show_span_clip=true)

plot_two_curves(rans_per_case, frac_nace_fuse_drag.*100, case_Nams, "Wing Drag Fraction [%]",
                joinpath(save_dir, "$(save_prefi)_WingDragFraction.png");
                show_span_clip=true)
