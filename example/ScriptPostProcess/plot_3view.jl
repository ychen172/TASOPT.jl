# Plot a simple 3-view (top / side / front) wireframe drawing of an aircraft, built from
# aircraft_frame in UtilityFunctions/AircraftGeometryFrame.jl. Single-bubble fuselage only;
# no interior/windows.
using TASOPT
include(TASOPT.__TASOPTindices__)
using Plots
include(joinpath(@__DIR__, "../UtilityFunctions/Plot_Aircraft/AircraftGeometryFrame.jl"))
include(joinpath(@__DIR__, "../UtilityFunctions/Plot_Aircraft/DepthRendering.jl"))

save_dir = joinpath(__TASOPTroot__, "../example/OutputAndPlots")
mkpath(save_dir)

"""
    plot_three_view(ac; save_dir, basename)

Build and save both the part-colored and flat-closest depth-shaded 4-panel (top/side/front/iso)
3-view wireframes for aircraft `ac`, as "<basename>_ThreeView.png" and
"<basename>_ThreeView_DepthFlatClosest.png" inside `save_dir`.
"""
function plot_three_view(ac; save_dir::String, basename::String)
    components = aircraft_frame(ac; upswept_tailcone=true, nose_droop_frac=0.2)

    p_top   = plot_view_by_part(components; xlabel="x [m]", ylabel="y [m]", title="Top View",   VIEW_ANGLES[:top]...)
    p_side  = plot_view_by_part(components; xlabel="x [m]", ylabel="z [m]", title="Side View",  VIEW_ANGLES[:side]...)
    p_front = plot_view_by_part(components; xlabel="y [m]", ylabel="z [m]", title="Front View", VIEW_ANGLES[:front]...)
    p_iso   = plot_view_by_part(components; xlabel="", ylabel="", title="3D View (45°az, 45°el)", VIEW_ANGLES[:iso]...)
    p = plot(p_top, p_side, p_front, p_iso, layout=(2, 2), size=(1600, 1600))
    path = joinpath(save_dir, "$(basename)_ThreeView.png")
    savefig(p, path)
    println("Saved 3-view plot to $path")

    p_top_fc   = plot_view_depth_flatclosest(components; xlabel="x [m]", ylabel="y [m]", title="Top View (flat-closest gray)",   VIEW_ANGLES[:top]...)
    p_side_fc  = plot_view_depth_flatclosest(components; xlabel="x [m]", ylabel="z [m]", title="Side View (flat-closest gray)",  VIEW_ANGLES[:side]...)
    p_front_fc = plot_view_depth_flatclosest(components; xlabel="y [m]", ylabel="z [m]", title="Front View (flat-closest gray)", VIEW_ANGLES[:front]...)
    p_iso_fc   = plot_view_depth_flatclosest(components; xlabel="", ylabel="", title="3D View (flat-closest gray)", VIEW_ANGLES[:iso]...)
    p_flatclosest = plot(p_top_fc, p_side_fc, p_front_fc, p_iso_fc, layout=(2, 2), size=(1600, 1600))
    path_fc = joinpath(save_dir, "$(basename)_ThreeView_DepthFlatClosest.png")
    savefig(p_flatclosest, path_fc)
    println("Saved flat-closest-gray 3-view plot to $path_fc")
    return nothing
end

#### Aircraft models to process: each a .jld2 saved via quicksave_aircraft, loaded AS-IS (no
#### size_aircraft! call -- these are already-sized saved models). `basename` sets the output
#### filename prefix for that model's two plots. Extend this list (or load+size a fresh aircraft
#### and wrap it in a similar NamedTuple) to add more models.
MODELS = [
    (file=joinpath(__TASOPTroot__, "../example/ModelSaved/Opti_Eth_NoACT_OAG_Ml_4Se_TyC_24Bf_TailCorr_V4_50.jld2"),  basename="Opti_Eth_4Se_50"),
    (file=joinpath(__TASOPTroot__, "../example/ModelSaved/Opti_Eth_NoACT_OAG_Ml_6Se_TyC_24Bf_TailCorr_V4_210.jld2"), basename="Opti_Eth_6Se_210"),
    (file=joinpath(__TASOPTroot__, "../example/ModelSaved/Opti_Jet_NoACT_OAG_Ml_4Se_TyC_24Bf_TailCorr_V4_50.jld2"),  basename="Opti_Jet_4Se_50"),
    (file=joinpath(__TASOPTroot__, "../example/ModelSaved/Opti_Jet_NoACT_OAG_Ml_6Se_TyC_24Bf_TailCorr_V4_210.jld2"), basename="Opti_Jet_6Se_210"),
]

for model in MODELS
    ac = quickload_aircraft(model.file)
    plot_three_view(ac; save_dir=save_dir, basename=model.basename)
end
