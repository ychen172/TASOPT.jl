# Plot a simple 3-view (top / side / front) wireframe drawing of an aircraft, built from
# aircraft_frame in UtilityFunctions/AircraftGeometryFrame.jl. Single-bubble fuselage only;
# no interior/windows.
using TASOPT
include(TASOPT.__TASOPTindices__)
using Plots
include(joinpath(@__DIR__, "../UtilityFunctions/Plot_Aircraft/AircraftGeometryFrame.jl"))

save_dir = joinpath(__TASOPTroot__, "../example/OutputAndPlots")
mkpath(save_dir)

#### Load and size the default aircraft model
ac = load_default_model()
size_aircraft!(ac)

components = aircraft_frame(ac; upswept_tailcone=true, nose_droop_frac=0.2)

const PART_COLORS = Dict(:fuselage => :black, :wing => :royalblue, :htail => :seagreen,
    :vtail => :purple, :engines => :crimson, :landing_gear => :gray40)

function plot_view_by_part(components::Dict{Symbol,GeomLines}; view::Symbol, xlabel::String, ylabel::String, title::String)
    p = plot(xlabel=xlabel, ylabel=ylabel, title=title, aspect_ratio=:equal, legend=false, dpi=300)
    for (part, frame) in components
        color = get(PART_COLORS, part, :black)
        for (xs, ys, zs) in frame
            u, v = view == :top ? (xs, ys) : view == :side ? (xs, zs) : (ys, zs)
            plot!(p, u, v, color=color, lw=1)
        end
    end
    return p
end

p_top   = plot_view_by_part(components; view=:top,   xlabel="x [m]", ylabel="y [m]", title="Top View")
p_side  = plot_view_by_part(components; view=:side,  xlabel="x [m]", ylabel="z [m]", title="Side View")
p_front = plot_view_by_part(components; view=:front, xlabel="y [m]", ylabel="z [m]", title="Front View")

p = plot(p_top, p_side, p_front, layout=(3, 1), size=(1000, 1600))
savefig(p, joinpath(save_dir, "ThreeView.png"))
println("Saved 3-view plot to $(joinpath(save_dir, "ThreeView.png"))")
