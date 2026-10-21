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

#### Load and size the default aircraft model
ac = load_default_model()
size_aircraft!(ac)

components = aircraft_frame(ac; upswept_tailcone=true, nose_droop_frac=0.2)

#### Part-colored 3-view
p_top   = plot_view_by_part(components; xlabel="x [m]", ylabel="y [m]", title="Top View",   VIEW_ANGLES[:top]...)
p_side  = plot_view_by_part(components; xlabel="x [m]", ylabel="z [m]", title="Side View",  VIEW_ANGLES[:side]...)
p_front = plot_view_by_part(components; xlabel="y [m]", ylabel="z [m]", title="Front View", VIEW_ANGLES[:front]...)
p_iso   = plot_view_by_part(components; xlabel="", ylabel="", title="3D View (45°az, 45°el)", VIEW_ANGLES[:iso]...)

p = plot(p_top, p_side, p_front, p_iso, layout=(2, 2), size=(1600, 1600))
savefig(p, joinpath(save_dir, "ThreeView.png"))
println("Saved 3-view plot to $(joinpath(save_dir, "ThreeView.png"))")

#### Flat-closest depth-shaded 3-view
p_top_fc   = plot_view_depth_flatclosest(components; xlabel="x [m]", ylabel="y [m]", title="Top View (flat-closest gray)",   VIEW_ANGLES[:top]...)
p_side_fc  = plot_view_depth_flatclosest(components; xlabel="x [m]", ylabel="z [m]", title="Side View (flat-closest gray)",  VIEW_ANGLES[:side]...)
p_front_fc = plot_view_depth_flatclosest(components; xlabel="y [m]", ylabel="z [m]", title="Front View (flat-closest gray)", VIEW_ANGLES[:front]...)
p_iso_fc   = plot_view_depth_flatclosest(components; xlabel="", ylabel="", title="3D View (flat-closest gray)", VIEW_ANGLES[:iso]...)

p_flatclosest = plot(p_top_fc, p_side_fc, p_front_fc, p_iso_fc, layout=(2, 2), size=(1600, 1600))
savefig(p_flatclosest, joinpath(save_dir, "ThreeView_DepthFlatClosest.png"))
println("Saved flat-closest-gray 3-view plot to $(joinpath(save_dir, "ThreeView_DepthFlatClosest.png"))")
