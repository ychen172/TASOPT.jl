# Overlay the top-view wireframe of several saved aircraft models on one plot, each aircraft
# drawn in a single flat color (not per-part), with every aircraft's wing box x-position shifted
# to x=0 so they're directly comparable regardless of their own absolute fuselage-station
# convention.
using TASOPT
include(TASOPT.__TASOPTindices__)
using Plots
include(joinpath(@__DIR__, "../UtilityFunctions/Plot_Aircraft/AircraftGeometryFrame.jl"))
include(joinpath(@__DIR__, "../UtilityFunctions/Plot_Aircraft/DepthRendering.jl"))

save_dir = joinpath(__TASOPTroot__, "../example/OutputAndPlots")
mkpath(save_dir)

"""Shift every x-coordinate in a `Dict{Symbol,GeomLines}` by `dx` (in place on new arrays,
original `components` left untouched)."""
function shift_frame_x(components::Dict{Symbol,GeomLines}, dx::Float64)
    shifted = Dict{Symbol,GeomLines}()
    for (part, frame) in components
        shifted[part] = [(xs .+ dx, ys, zs) for (xs, ys, zs) in frame]
    end
    return shifted
end

#### Models to overlay, each with a single flat color for ALL of its components.
MODELS = [
    (file="Opti_Jet_NoACT_OAG_Ml_6Se_TyC_24Bf_TailCorr_V4_210.jld2", color=:red),
    (file="Opti_Jet_NoACT_OAG_Ml_4Se_TyC_24Bf_TailCorr_V4_100.jld2", color=:blue),
    (file="Opti_Jet_NoACT_OAG_Ml_4Se_TyC_24Bf_TailCorr_V4_50.jld2",  color=:black),
]

p = nothing
for model in MODELS
    path = joinpath(__TASOPTroot__, "../example/ModelSaved/$(model.file)")
    ac = quickload_aircraft(path)
    components = aircraft_frame(ac; upswept_tailcone=true, nose_droop_frac=0.2)

    dx = -ac.wing.layout.box_x #shift so the wing box x-position lands at x=0
    components = shift_frame_x(components, dx)

    flat_colors = Dict(part => model.color for part in keys(components))
    global p = plot_view_by_part(components; xlabel="x [m]", ylabel="y [m]",
        title="Top View Overlay (wing-box-aligned)", VIEW_ANGLES[:top]..., part_colors=flat_colors, p=p)
end

path = joinpath(save_dir, "TopView_Overlay.png")
savefig(p, path)
println("Saved top-view overlay plot to $path")
