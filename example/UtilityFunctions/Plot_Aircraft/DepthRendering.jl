"""
Plots.jl-dependent rendering helpers for an `aircraft_frame(ac)` wireframe: turning its
`Dict{Symbol,GeomLines}` into an actual top/side/front/iso 2D plot, with either per-part colors
or depth-based grayscale shading. Kept separate from `AircraftGeometryFrame.jl`, which is pure
geometry/math with no `Plots.jl` dependency and no rendering-backend assumptions -- this file is
where that choice gets made. Assumes `AircraftGeometryFrame.jl` (for `GeomLines`,
`axonometric_xy`, `axonometric_depth`) has already been included by the caller.
"""

const PART_COLORS = Dict(:fuselage => :black, :wing => :royalblue, :htail => :seagreen,
    :vtail => :purple, :engines => :crimson, :landing_gear => :gray40)

"""
    VIEW_ANGLES

The classic orthographic top/side/front views are just specific `(azimuth_deg, elevation_deg)`
pairs of the same general axonometric projection (`axonometric_xy`/`axonometric_depth`) used for
the oblique 3D view -- confirmed by both direct derivation and visual check:
- `:top   => (0, 90)`: looking straight down. `(u,v)` and depth both match the plain `(x,y)`/`-z`
  formulas exactly.
- `:side  => (0, 0)`: looking from the side, level. `(u,v)` and depth both match the plain
  `(x,z)`/`y` formulas exactly. (Depth convention: with the nose drawn at the left and z up, the
  out-of-page/near direction works out to `-y`, the LEFT (port) wing -- see `point_depth`.)
- `:front => (90, 0)`: nose-on. Depth matches the plain `x` formula exactly (nose nearest). The
  `(u,v)` position comes out mirrored (`u=-y` instead of `y`) relative to the plain `(y,z)`
  formula, but since every TASOPT aircraft is left-right symmetric, that mirroring is invisible
  in the rendered shape -- only the depth/shading convention needs to match, and it does.
So there is exactly one formula for position and one for depth, shared by every view; `:top`/
`:side`/`:front` are just named entries in this table, not separately hand-coded cases.
"""
const VIEW_ANGLES = Dict(
    :top   => (azimuth_deg=0.0, elevation_deg=90.0),
    :side  => (azimuth_deg=0.0, elevation_deg=0.0),
    :front => (azimuth_deg=90.0, elevation_deg=0.0),
    :iso   => (azimuth_deg=45.0, elevation_deg=45.0),
)

"""
Per-point "depth into the screen" along one polyline: the shared `axonometric_depth` rotation,
at whatever `(azimuth_deg, elevation_deg)` the caller passes (see `VIEW_ANGLES`). Convention:
smaller returned value = nearer the viewer = should be drawn darker.
"""
point_depth(xs::Vector{Float64}, ys::Vector{Float64}, zs::Vector{Float64};
    azimuth_deg::Float64, elevation_deg::Float64) =
    axonometric_depth(xs, ys, zs; azimuth_deg=azimuth_deg, elevation_deg=elevation_deg)

"""
    depth_sorted_segments(entries)

Painter's-algorithm sort, at 2-point *segment* granularity rather than whole-line granularity.
A whole-line average depth is wrong whenever a single polyline itself spans a wide depth range --
e.g. a fuselage meridian (`axisym_body_frame_x`) runs the *entire* nose-to-tail length at one
circumferential angle, so in the front view (depth = x) its nose end is the nearest point in the
whole plot while its own tail end is among the farthest: no single order for that whole line can
be right. Splitting into 2-point segments and sorting each independently (farthest first, nearest
last) lets the nose end of that very meridian draw on top of far structure while its tail end
still draws underneath near structure, same as every other line.

`entries` is a `Vector` of `(u, v, d, payload)` -- equal-length coordinate/depth vectors for one
polyline, plus whatever per-line payload the caller needs carried through (a color, `nothing`,
...). Returns the flattened, sorted list of `(u2, v2, d2, payload)` 2-point segments.
"""
function depth_sorted_segments(entries)
    segs = Tuple{Vector{Float64},Vector{Float64},Vector{Float64},Any}[]
    for (u, v, d, payload) in entries
        for i in 1:length(u)-1
            push!(segs, (u[i:i+1], v[i:i+1], d[i:i+1], payload))
        end
    end
    sort!(segs, by=s -> (s[3][1] + s[3][2]) / 2, rev=true)
    return segs
end

#### Both plotting functions below draw segments back-to-front (farthest depth first, nearest
#### last), painter's-algorithm style, so nearer structure visually overlaps/occludes farther
#### structure instead of draw order being arbitrary dict/frame iteration order. Both take a
#### plain (azimuth_deg, elevation_deg) pair -- look up VIEW_ANGLES[:top/:side/:front/:iso] at
#### the call site -- rather than a named `view` symbol, since every view is just a point in the
#### same (az, el) space.

function plot_view_by_part(components::Dict{Symbol,GeomLines}; xlabel::String, ylabel::String, title::String,
    azimuth_deg::Float64, elevation_deg::Float64)

    entries = Tuple{Vector{Float64},Vector{Float64},Vector{Float64},Any}[]
    for (part, frame) in components
        color = get(PART_COLORS, part, :black)
        for (xs, ys, zs) in frame
            u, v = axonometric_xy(xs, ys, zs; azimuth_deg=azimuth_deg, elevation_deg=elevation_deg)
            d = point_depth(xs, ys, zs; azimuth_deg=azimuth_deg, elevation_deg=elevation_deg)
            push!(entries, (u, v, d, color))
        end
    end
    segs = depth_sorted_segments(entries)

    p = plot(xlabel=xlabel, ylabel=ylabel, title=title, aspect_ratio=:equal, legend=false, dpi=800)
    for (u, v, _, color) in segs
        plot!(p, u, v, color=color, lw=1)
    end
    return p
end

"""Flat-per-segment depth-shaded variant: give each 2-point segment ONE flat color taken from
its closer (darker/smaller-depth) endpoint -- so a segment never looks lighter than its own
nearest point actually is. Preferred over interpolating each segment's color along its length
(an earlier variant, dropped): the flat-closest coloring makes components like the horizontal
tail read more clearly against the fuselage in the side view."""
function plot_view_depth_flatclosest(components::Dict{Symbol,GeomLines}; xlabel::String, ylabel::String, title::String,
    azimuth_deg::Float64, elevation_deg::Float64, gray_near::Float64=0.0, gray_far::Float64=0.82)

    entries = Tuple{Vector{Float64},Vector{Float64},Vector{Float64},Nothing}[]
    for frame in values(components), (xs, ys, zs) in frame
        u, v = axonometric_xy(xs, ys, zs; azimuth_deg=azimuth_deg, elevation_deg=elevation_deg)
        d = point_depth(xs, ys, zs; azimuth_deg=azimuth_deg, elevation_deg=elevation_deg)
        push!(entries, (u, v, d, nothing))
    end

    dmin = minimum(minimum(e[3]) for e in entries)
    dmax = maximum(maximum(e[3]) for e in entries)
    span = dmax - dmin

    segs = depth_sorted_segments(entries)

    p = plot(xlabel=xlabel, ylabel=ylabel, title=title, aspect_ratio=:equal, legend=false, dpi=800)
    for (u, v, d, _) in segs
        dmin_seg = min(d[1], d[2]) #closer (darker) of the two endpoints
        t = span > 0 ? (dmin_seg - dmin) / span : 0.0
        g = gray_near + t * (gray_far - gray_near)
        plot!(p, u, v, color=Plots.Gray(g), lw=1.5)
    end
    return p
end
