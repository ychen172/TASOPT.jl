"""
Build a simple 3D wireframe ("frame-net") representation of an aircraft's external geometry,
component by component, for single-bubble-fuselage aircraft. No interior/windows/seats.

Each component function returns a `Vector{NTuple{3,Vector{Float64}}}`: a list of polylines,
each a (xs, ys, zs) tuple of equal-length coordinate vectors. `aircraft_frame(ac)` collects
every component into one `Dict{Symbol,Any}`, and `project_view` flattens that into 2D line
segments for a top/side/front view.

Axes follow TASOPT's own convention: x = longitudinal (nose to tail), y = lateral (spanwise),
z = vertical (z_wing/z_Htail sign convention, positive up).
"""

const GeomLines = Vector{NTuple{3,Vector{Float64}}}

include(joinpath(@__DIR__, "AirfoilShape.jl"))

#### Generic axisymmetric-body helpers ---------------------------------------

"""Circle of radius `r` in the y-z plane at station `x`, centered at (x, y0, z0)."""
function _circle_yz(x::Float64, r::Float64, y0::Float64, z0::Float64; n_theta::Int=16)
    θ = range(0, 2π, length=n_theta + 1)
    xs = fill(x, n_theta + 1)
    ys = y0 .+ r .* sin.(θ)
    zs = z0 .+ r .* cos.(θ)
    return (xs, ys, zs)
end

"""Circle of radius `r` in the x-y plane at height `z`, centered at (x0, y0, z)."""
function _circle_xy(z::Float64, r::Float64, x0::Float64, y0::Float64; n_theta::Int=16)
    θ = range(0, 2π, length=n_theta + 1)
    xs = x0 .+ r .* cos.(θ)
    ys = y0 .+ r .* sin.(θ)
    zs = fill(z, n_theta + 1)
    return (xs, ys, zs)
end

"""Circle of radius `r` in the x-z plane at station `y`, centered at (x0, y, z0)."""
function _circle_xz(y::Float64, r::Float64, x0::Float64, z0::Float64; n_theta::Int=16)
    θ = range(0, 2π, length=n_theta + 1)
    xs = x0 .+ r .* sin.(θ)
    ys = fill(y, n_theta + 1)
    zs = z0 .+ r .* cos.(θ)
    return (xs, ys, zs)
end

"""
    axisym_body_frame_x(xs, rs; y0, z0, n_theta, n_rings)

Wireframe of a body of revolution about a line parallel to the x-axis (through (y0,z0)),
with radius `rs[i]` at station `xs[i]`. Used for the fuselage and the engine nacelles.
"""
function axisym_body_frame_x(xs::Vector{Float64}, rs::Vector{Float64};
    y0::Float64=0.0, z0::Union{Float64,Vector{Float64}}=0.0, n_theta::Int=8, n_rings::Int=8)

    @assert length(xs) == length(rs)
    z0v = z0 isa Vector ? z0 : fill(z0, length(xs))
    @assert length(z0v) == length(xs) "z0 vector must match xs in length (one center-height per station)"
    frame = GeomLines()

    n_rings = min(n_rings, length(xs))
    ring_idx = unique(round.(Int, range(1, length(xs), length=n_rings)))
    for i in ring_idx
        push!(frame, _circle_yz(xs[i], rs[i], y0, z0v[i]; n_theta=n_theta))
    end

    θs = range(0, 2π, length=n_theta + 1)[1:end-1]
    for θ in θs
        ys = y0 .+ rs .* sin(θ)
        zs = z0v .+ rs .* cos(θ)
        push!(frame, (copy(xs), ys, zs))
    end
    return frame
end

"""
    axisym_body_frame_z(zs, rs; x0, y0, n_theta, n_rings)

Wireframe of a body of revolution about a line parallel to the z-axis (through (x0,y0)),
with radius `rs[i]` at height `zs[i]`. Used for the landing gear struts.
"""
function axisym_body_frame_z(zs::Vector{Float64}, rs::Vector{Float64};
    x0::Float64=0.0, y0::Float64=0.0, n_theta::Int=8, n_rings::Int=8)

    @assert length(zs) == length(rs)
    frame = GeomLines()

    n_rings = min(n_rings, length(zs))
    ring_idx = unique(round.(Int, range(1, length(zs), length=n_rings)))
    for i in ring_idx
        push!(frame, _circle_xy(zs[i], rs[i], x0, y0; n_theta=n_theta))
    end

    θs = range(0, 2π, length=n_theta + 1)[1:end-1]
    for θ in θs
        xs = x0 .+ rs .* cos(θ)
        ys = y0 .+ rs .* sin(θ)
        push!(frame, (xs, ys, copy(zs)))
    end
    return frame
end

"""
    axisym_body_frame_y(ys, rs; x0, z0, n_theta, n_rings)

Wireframe of a body of revolution about a line parallel to the y-axis (through (x0,z0)), with
radius `rs[i]` at station `ys[i]`. Used for the landing-gear wheels: a real wheel's axle points
along y, so its circular rolling profile is visible in the x-z (side-view) plane, not as a flat
circle in the ground (x-y) plane.
"""
function axisym_body_frame_y(ys::Vector{Float64}, rs::Vector{Float64};
    x0::Float64=0.0, z0::Float64=0.0, n_theta::Int=8, n_rings::Int=8)

    @assert length(ys) == length(rs)
    frame = GeomLines()

    n_rings = min(n_rings, length(ys))
    ring_idx = unique(round.(Int, range(1, length(ys), length=n_rings)))
    for i in ring_idx
        push!(frame, _circle_xz(ys[i], rs[i], x0, z0; n_theta=n_theta))
    end

    θs = range(0, 2π, length=n_theta + 1)[1:end-1]
    for θ in θs
        xs = x0 .+ rs .* sin(θ)
        zs = z0 .+ rs .* cos(θ)
        push!(frame, (xs, copy(ys), zs))
    end
    return frame
end

#### Fuselage ------------------------------------------------------------------

"""
    fuselage_radius(fuse, x)

TASOPT's own axisymmetric profile-drag radius function R(x) (see
`docs/src/aero/theory_fuse_profile_drag.md`): nose/cylinder/tail-cone radius at longitudinal
station `x`, the same curve already used for the top-view outline in `stickfig`. Only valid for
`fuse.layout` with zero bubble offset (single-bubble).
"""
function fuselage_radius(fuse, x::Float64)
    Rfuse = fuse.layout.radius
    anose = fuse.layout.nose_radius
    btail = fuse.layout.tail_radius

    xnose   = fuse.layout.x_nose
    xblend1 = fuse.layout.x_start_cylinder
    xblend2 = fuse.layout.x_end_cylinder
    xend    = fuse.layout.x_end

    tapers_to_point = TASOPT.compare_strings(fuse.layout.opt_tapers_to, "point")
    dytail_frac = tapers_to_point ? 1.0 : 0.2 #Fraction of Rfuse the tail tapers down by

    if x <= xnose
        return 0.0
    elseif x < xblend1
        fracx = (xblend1 - x) / (xblend1 - xnose)
        return Rfuse * (1 - fracx^anose)^(1 / anose)
    elseif x <= xblend2
        return Rfuse
    elseif x < xend
        fracx = (x - xblend2) / (xend - xblend2)
        return Rfuse * (1 - dytail_frac * fracx^btail)
    else
        return Rfuse * (1 - dytail_frac)
    end
end

"""
    fuselage_skirt_height(fuse, x)

TASOPT's real single-bubble cross-section isn't a plain circle even at zero side-lobes: its
structural shell area/perimeter (`area`/`get_perimeter` in `fuselage_geometry.jl`) are
`π*R^2 + 2*R*ΔR` / `2π*R + 2*ΔR`, i.e. a circle of radius `R` with a flat-sided rectangular
"skirt" of height `ΔR = bubble_lower_downward_shift` hanging below it (upper semicircle, two
straight sides, lower semicircle shifted down by `ΔR`) -- not the double-bubble side-lobe
geometry (`n_webs`/`wfb`), which TASOPT's own `MultiBubble` cross-section adds separately.
This taper law (same `fracx` shape as `fuselage_radius`, just applied to `ΔR` instead of `R`)
carries that skirt height from its full cylinder value down to `0` at the nose tip and tail end,
so the belly bulge fades out smoothly rather than cutting off abruptly at the cylinder boundary.
"""
function fuselage_skirt_height(fuse, x::Float64)
    ΔRmax = fuse.layout.bubble_lower_downward_shift
    ΔRmax == 0.0 && return 0.0

    anose = fuse.layout.nose_radius
    btail = fuse.layout.tail_radius
    xnose   = fuse.layout.x_nose
    xblend1 = fuse.layout.x_start_cylinder
    xblend2 = fuse.layout.x_end_cylinder
    xend    = fuse.layout.x_end
    tapers_to_point = TASOPT.compare_strings(fuse.layout.opt_tapers_to, "point")
    dytail_frac = tapers_to_point ? 1.0 : 0.2

    if x <= xnose
        return 0.0
    elseif x < xblend1
        fracx = (xblend1 - x) / (xblend1 - xnose)
        return ΔRmax * (1 - fracx^anose)^(1 / anose)
    elseif x <= xblend2
        return ΔRmax
    elseif x < xend
        fracx = (x - xblend2) / (xend - xblend2)
        return ΔRmax * (1 - dytail_frac * fracx^btail)
    else
        return ΔRmax * (1 - dytail_frac)
    end
end

"""
    fuselage_cross_section_outline(R, ΔR; n_arc=8) -> (ys, zs)

Closed-loop "circle + skirt" cross-section outline (see `fuselage_skirt_height`): upper
semicircle of radius `R` centered at `z=0` (unchanged from a plain circle), straight sides of
length `ΔR`, lower semicircle of radius `R` centered at `z=-ΔR`. Degenerates exactly to a plain
circle when `ΔR=0` (no discontinuity). `n_arc` points per semicircle; the returned loop repeats
its first point at the end (closed).
"""
function fuselage_cross_section_outline(R::Float64, ΔR::Float64; n_arc::Int=8)
    θ_upper = range(0, π, length=n_arc + 1)
    θ_lower = range(π, 2π, length=n_arc + 1)

    ys = Float64[]; zs = Float64[]
    for θ in θ_upper
        push!(ys, R * cos(θ)); push!(zs, R * sin(θ))
    end
    for θ in θ_lower
        push!(ys, R * cos(θ)); push!(zs, R * sin(θ) - ΔR)
    end
    push!(ys, ys[1]); push!(zs, zs[1]) #close the loop
    return ys, zs
end

"""
    fuselage_tailcone_z_bend(fuse, x)

Artificial, drawing-only vertical offset of the fuselage's local center axis within the tail
cone, chosen so the *top* surface stays flat at `Rfuse` (`R(x)+z(x) = Rfuse`) while the bottom
sweeps upward -- the classic real-aircraft "upswept tail" look. **Not** a TASOPT-modeled
quantity. Zero for `x <= x_end_cylinder` (nose and cylinder untouched, continuous with the
tail-cone formula at that boundary since `R(x_end_cylinder) = Rfuse` there).
"""
function fuselage_tailcone_z_bend(fuse, x::Float64)
    xblend2 = fuse.layout.x_end_cylinder
    x <= xblend2 && return 0.0
    return fuse.layout.radius - fuselage_radius(fuse, x)
end

"""
    fuselage_nose_z_droop(fuse, x, droop_frac)

Artificial, drawing-only vertical offset of the fuselage's local center axis within the nose
taper, a **plain linear** droop (not following `R(x)`'s own nonlinear shape the way the tail
bend does): `0` at the cylinder start (`x_start_cylinder`), down to `-droop_frac*Rfuse` at the
nose tip (`x_nose`), linear in between -- **clamped** so the *bottom* surface (`z(x)-R(x)`) never
dips below the cylinder's own belly line (`z_fuse-Rfuse = -Rfuse`); a straight-line center droop
combined with the nose radius shrinking toward the tip could otherwise push the bottom skin
below that reference partway along the taper, even though both ends (`z=0` at the cylinder,
`z=-droop_frac*Rfuse` at the point tip where `R=0`) are individually fine. `droop_frac=0` (the
default when this option is off) means no droop at all. **Not** a TASOPT-modeled quantity.
"""
function fuselage_nose_z_droop(fuse, x::Float64, droop_frac::Float64)
    droop_frac == 0.0 && return 0.0
    Rfuse   = fuse.layout.radius
    xnose   = fuse.layout.x_nose
    xblend1 = fuse.layout.x_start_cylinder
    x >= xblend1 && return 0.0

    fracx = clamp((xblend1 - x) / (xblend1 - xnose), 0.0, 1.0) #0 at xblend1, 1 at/before xnose
    z_raw = -droop_frac * Rfuse * fracx
    z_floor = fuselage_radius(fuse, x) - Rfuse #z(x) - R(x) >= -Rfuse  <=>  z(x) >= R(x) - Rfuse
    return max(z_raw, z_floor)
end

"""
    fuselage_frame(fuse; n_theta=24, n_rings=8, n_profile=160, upswept_tailcone=false, nose_droop_frac=0.0)

Wireframe of the real "circle + skirt" single-bubble fuselage shape (`fuselage_radius` +
`fuselage_skirt_height`, see `fuselage_cross_section_outline`) -- not a plain circle (that was an
earlier simplification in this file; TASOPT's own structural shell area/perimeter formulas
include the skirt even at zero side-lobes). With `upswept_tailcone=true`, the tail cone's center
axis is bent upward per `fuselage_tailcone_z_bend`; with `nose_droop_frac != 0`, the nose's
center axis droops down (linearly) per `fuselage_nose_z_droop` -- both visual-only, not TASOPT
quantities, and independent of each other (different x-ranges, so they never interact). Defaults
(`false`/`0.0`) still draw the plain circle+skirt shape, just without either bend.

`n_theta` here means the same thing it always has (longitudinal-line/ring angular resolution),
but is now split in half internally (`n_arc = n_theta÷2`) between the upper and lower semicircles
of `fuselage_cross_section_outline`, so the total point count per ring stays comparable to the
old plain-circle version.
"""
function fuselage_frame(fuse; n_theta::Int=24, n_rings::Int=8, n_profile::Int=160,
    upswept_tailcone::Bool=false, nose_droop_frac::Float64=0.0)

    xnose = fuse.layout.x_nose
    xend  = fuse.layout.x_end

    xs = collect(range(xnose, xend, length=n_profile))
    rs = fuselage_radius.(Ref(fuse), xs)
    drs = fuselage_skirt_height.(Ref(fuse), xs)
    z_tail = upswept_tailcone ? fuselage_tailcone_z_bend.(Ref(fuse), xs) : zeros(length(xs))
    z_nose = fuselage_nose_z_droop.(Ref(fuse), xs, nose_droop_frac)
    z0 = z_tail .+ z_nose

    n_arc = max(n_theta ÷ 2, 1)
    outlines = [fuselage_cross_section_outline(rs[i], drs[i]; n_arc=n_arc) for i in eachindex(xs)]
    n_pts = length(outlines[1][1]) #same for every station: 2*(n_arc+1)+1 (closed loop)

    frame = GeomLines()

    #Cross-section rings at a handful of x stations
    n_rings_eff = min(n_rings, length(xs))
    ring_idx = unique(round.(Int, range(1, length(xs), length=n_rings_eff)))
    for i in ring_idx
        ys, zs = outlines[i]
        push!(frame, (fill(xs[i], n_pts), ys, zs .+ z0[i]))
    end

    #Longitudinal meridian lines: same outline-point index, swept across every x station
    #(last point duplicates the first to close each ring's loop -- skip it here, it's not a
    #distinct longitudinal line)
    for k in 1:n_pts-1
        ys_k = [outlines[i][1][k] for i in eachindex(xs)]
        zs_k = [outlines[i][2][k] + z0[i] for i in eachindex(xs)]
        push!(frame, (copy(xs), ys_k, zs_k))
    end

    return frame
end

#### Wing / horizontal tail / vertical tail -----------------------------------

"""
    wing_planform_outline(layout, λs, λt; xcLE=-0.4, xcTE=0.6, mirrored=true)

Planform outline (root LE -> break LE -> tip LE -> tip TE -> break TE -> root TE), same geometry
`stickfig` draws. Returns `(xw, yw)`, both length-6: `xw` the longitudinal position of each
vertex, `yw` its distance from the centerline (always >= 0).

`layout.span`/`root_span` are interpreted differently depending on `mirrored`:
- `mirrored=true` (wing, htail): `span` is the *full* span across both sides, so each drawn
  side only extends out to `span/2` -- this is what `stickfig` assumes.
- `mirrored=false` (vtail): after `size_aircraft!`, `vtail.layout.span` is already the real,
  single-fin height (see `size_aircraft.jl:560-588`, where an internal doubling trick used to
  reuse the wing span-loading code is explicitly undone before the final value is stored) --
  so the full, un-halved value is used directly.
"""
function wing_planform_outline(layout, λs::Float64, λt::Float64; xcLE::Float64=-0.4, xcTE::Float64=0.6,
    mirrored::Bool=true)

    co = layout.root_chord
    cs = co * λs
    ct = co * λt

    scale = mirrored ? 0.5 : 1.0
    bo = layout.root_span
    b  = layout.span
    bs = layout.ηs * b

    dx = layout.box_x
    tanL = tan(layout.sweep * π / 180.0)
    xoff_s = tanL * (bs - bo) * scale
    xoff_t = tanL * (b  - bo) * scale

    xw = [co * xcLE + dx, xoff_s + cs * xcLE + dx, xoff_t + ct * xcLE + dx,
          xoff_t + ct * xcTE + dx, xoff_s + cs * xcTE + dx, co * xcTE + dx]
    yw = [bo, bs, b, b, bs, bo] .* scale
    return xw, yw
end

"""
    wing_local_chord_toc(wing, span_coord; mirrored=true)

Local chord, thickness-to-chord, and leading-edge x position at spanwise position `span_coord`
(an absolute distance from the centerline/root axis -- for `mirrored=true`, wing/htail, one
side's y; for `mirrored=false`, vtail, height above the root). Uses the same two-panel
root->break->tip taper/interpolation law `wing_drag.jl`'s profile-drag integration does: the
inboard panel (root->break) linearly interpolates thickness ratio from
`inboard.thickness_to_chord` to `outboard.thickness_to_chord`; the outboard panel (break->tip)
holds that ratio flat at `outboard.thickness_to_chord` (chord still tapers normally via
`outboard.λ` -- only the *ratio* is held flat, see `wing_drag.jl:105,138`). For htail/vtail,
where the break degenerates to the root (`ηs=ηo`, `inboard.λ=1`), this collapses correctly to a
single outboard-only taper using only `outboard.thickness_to_chord` the whole way -- confirmed
against the real sized default aircraft, where `ηs==ηo` and `inboard.λ==1.0` for both tails.

The inboard/outboard split uses `η < ηs` (strict), not `<=`: at the exact root/break boundary
both branches agree for a real wing (continuous), but for htail/vtail's degenerate zero-width
inboard panel (`ηs==ηo`), only `<` reliably lands in the outboard branch at that single
coincident point -- `<=` would pick `inboard.thickness_to_chord` there instead, silently wrong
whenever it differs from `outboard.thickness_to_chord` (masked on this aircraft only because the
two happen to be numerically equal).
"""
function wing_local_chord_toc(wing, span_coord::Float64; mirrored::Bool=true)
    layout = wing.layout
    co = layout.root_chord
    λs = wing.inboard.λ
    λt = wing.outboard.λ
    bo = layout.root_span
    b  = layout.span
    ηs = layout.ηs
    ηo = bo / b
    b_half = mirrored ? b / 2 : b
    η = span_coord / b_half

    if η < ηs
        frac  = (η - ηo) / max(ηs - ηo, eps())
        chord = co * (1 + (λs - 1) * frac)
        toc   = wing.inboard.cross_section.thickness_to_chord * (1 - frac) +
                wing.outboard.cross_section.thickness_to_chord * frac
    else
        frac  = (η - ηs) / (1 - ηs)
        chord = co * (λs + (λt - λs) * frac)
        toc   = wing.outboard.cross_section.thickness_to_chord
    end

    tanL = tan(layout.sweep * π / 180.0)
    xoff_scale = mirrored ? 0.5 : 1.0
    xLE = tanL * (span_coord - bo * xoff_scale) - 0.4 * chord + layout.box_x
    return chord, toc, xLE
end

"""
    wing_frame(wing; dihedral_angle=0.0, z_root_override=nothing, mirror=true, y_stations=nothing, symmetric=false, dat_path=DEFAULT_AIRFOIL_DAT)

Complete wireframe for the wing: flat LE/TE planform outline (extended inward to the centerline
as a flat "centerbox", see below), real airfoil cross-sections (`AirfoilShape.jl`) at the
centerbox and at a handful of spanwise stations (default: root, break, tip), and two continuous
spanwise lines (upper crest, lower trough) threading through those cross-sections' thickness
envelope. The leading/trailing edges are *not* duplicated as separate spanwise lines here --
confirmed they'd trace the exact same points as the flat outline already does (same formula).

`dihedral_angle` is the same borrowed landing-gear assumption `engine_frame`/`landing_gear_frame`
already use (TASOPT has no real wing dihedral anywhere else) -- pass `0.0` (the default) for a
flat wing. `z_root_override`, when given, replaces `wing.layout.z` as the root reference
entirely. `wing.layout.z` is otherwise undocumented anywhere in TASOPT (used only in the
landing-gear clearance formula and this drawing code), so `aircraft_frame` defines it instead as
the wing's leading-edge height, with the dihedral extension's root value (`y=0`) set flush with
the fuselage's own true belly line -- a drawing-only convention choice, not a recovered fact.

The inner "wingbox" section (`y=0` to `y=bo/2`, inside the fuselage) is drawn **flat/horizontal**
(dihedral doesn't apply there -- dihedral is conventionally taken to start at the fuselage side,
not the centerline) but at the dihedral-*extrapolated* height, i.e. whatever z the normal
`z_root + y*tan(dihedral_angle)` formula gives at `y=bo/2` -- not `z_root` itself. That's what
makes the inner box meet the outboard panel with zero kink: the outboard panel's own z at
`y=bo/2` is exactly that same extrapolated value. The box's cross-sections (`y=0` and `y=bo/2`)
use one **fixed** airfoil: the root's own chord/thickness-to-chord
(`wing_local_chord_toc(wing, bo/2)`), not a continuation of the outboard taper law, which isn't
even defined for `y<bo/2`. New centerline points reuse the root's own LE/TE x position unswept,
consistent with a simple rectangular carry-through box.
"""
function wing_frame(wing; dihedral_angle::Float64=0.0, z_root_override::Union{Nothing,Float64}=nothing,
    mirror::Bool=true, y_stations::Union{Nothing,Vector{Float64}}=nothing,
    symmetric::Bool=false, dat_path::String=DEFAULT_AIRFOIL_DAT)

    layout = wing.layout
    λs = wing.inboard.λ
    λt = wing.outboard.λ
    bo, b, ηs = layout.root_span, layout.span, layout.ηs
    z_root = z_root_override !== nothing ? z_root_override : layout.z

    #### Flat LE/TE planform outline, extended inward to the centerline (flat centerbox)
    xw, yw = wing_planform_outline(layout, λs, λt; mirrored=true)
    zw = z_root .+ yw .* tan(dihedral_angle)
    z_box = zw[1] #dihedral-extrapolated height at y=bo/2 (yw[1]==yw[end]==bo/2), flat inboard of it
    xw_ext = [xw[1]; xw; xw[end]]
    yw_ext = [0.0; yw; 0.0]
    zw_ext = [z_box; zw; z_box]

    frame = GeomLines()
    push!(frame, (xw_ext, yw_ext, zw_ext))
    if mirror
        push!(frame, (xw_ext, -yw_ext, zw_ext))
    end

    #### Real airfoil cross-sections + crest/trough spanwise lines, centerbox through tip
    chord0, toc0, xLE0 = wing_local_chord_toc(wing, bo / 2) #Fixed root cross-section for the centerbox
    stations = y_stations !== nothing ? y_stations : [bo / 2, ηs * b / 2, b / 2]
    station_data = [(0.0, chord0, toc0, xLE0, z_box)]
    for y in stations
        chord, toc, xLE = wing_local_chord_toc(wing, y)
        push!(station_data, (y, chord, toc, xLE, z_root + y * tan(dihedral_angle)))
    end

    crests = NTuple{3,Float64}[]; troughs = NTuple{3,Float64}[]
    for (y, chord, toc, xLE, z) in station_data
        xu, yu, xl, yl = airfoil_section(chord, toc; dat_path=dat_path, xLE=xLE, yref=z, symmetric=symmetric)
        push!(frame, (xu, fill(y, length(xu)), yu))
        push!(frame, (xl, fill(y, length(xl)), yl))
        if mirror
            push!(frame, (xu, fill(-y, length(xu)), yu))
            push!(frame, (xl, fill(-y, length(xl)), yl))
        end
        _, _, crest, trough = airfoil_key_points(xu, yu, xl, yl)
        push!(crests, (crest[1], y, crest[2]))
        push!(troughs, (trough[1], y, trough[2]))
    end

    for pts in (crests, troughs)
        allpts = mirror ? vcat([(x, -y, z) for (x, y, z) in reverse(pts)], pts) : pts
        xs = [p[1] for p in allpts]; ys = [p[2] for p in allpts]; zs = [p[3] for p in allpts]
        push!(frame, (xs, ys, zs))
    end

    return frame
end

"""
    htail_frame(htail; dihedral_angle=0.0, extra_z=0.0, z_root_override=nothing, mirror=true, y_stations=nothing, symmetric=false, dat_path=DEFAULT_AIRFOIL_DAT)

Complete wireframe for the horizontal tail: flat LE/TE planform outline (from `y=bo/2` outward
only -- unlike the wing, the two htail halves are **not** connected by a carry-through box; they
stay visually separated at the centerline), real airfoil cross-sections at a handful of spanwise
stations (default: root, break, tip), and the two continuous spanwise crest/trough lines through
them (LE/TE spanwise lines are skipped, same reasoning as `wing_frame`: the flat outline already
traces them exactly).

`extra_z` is a uniform rigid vertical shift, meant for `aircraft_frame`'s `upswept_tailcone`
option (the htail sits on the bent tail-cone centerline, not the flat `z=0` axis) -- applied as
one constant for the whole tail rather than per-point, since evaluating the bend separately at
each chord point would tilt the chord (an artificial, un-modeled incidence angle), not just
translate it. `dihedral_angle`/`z_root_override` work the same as `wing_frame`, except htail
keeps its own `layout.z` convention by default (not overridden by `aircraft_frame`).
"""
function htail_frame(htail; dihedral_angle::Float64=0.0, extra_z::Float64=0.0,
    z_root_override::Union{Nothing,Float64}=nothing, mirror::Bool=true,
    y_stations::Union{Nothing,Vector{Float64}}=nothing, symmetric::Bool=false,
    dat_path::String=DEFAULT_AIRFOIL_DAT)

    layout = htail.layout
    λs = htail.inboard.λ
    λt = htail.outboard.λ
    bo, b, ηs = layout.root_span, layout.span, layout.ηs
    z_root = z_root_override !== nothing ? z_root_override : layout.z
    z0 = z_root + extra_z

    xw, yw = wing_planform_outline(layout, λs, λt; mirrored=true)
    zw = z0 .+ yw .* tan(dihedral_angle)

    frame = GeomLines()
    push!(frame, (xw, yw, zw))
    if mirror
        push!(frame, (xw, -yw, zw))
    end

    stations = y_stations !== nothing ? y_stations : [bo / 2, ηs * b / 2, b / 2]
    crests = NTuple{3,Float64}[]; troughs = NTuple{3,Float64}[]
    for y in stations
        chord, toc, xLE = wing_local_chord_toc(htail, y)
        z = z0 + y * tan(dihedral_angle)
        xu, yu, xl, yl = airfoil_section(chord, toc; dat_path=dat_path, xLE=xLE, yref=z, symmetric=symmetric)
        push!(frame, (xu, fill(y, length(xu)), yu))
        push!(frame, (xl, fill(y, length(xl)), yl))
        if mirror
            push!(frame, (xu, fill(-y, length(xu)), yu))
            push!(frame, (xl, fill(-y, length(xl)), yl))
        end
        _, _, crest, trough = airfoil_key_points(xu, yu, xl, yl)
        push!(crests, (crest[1], y, crest[2]))
        push!(troughs, (trough[1], y, trough[2]))
    end

    for pts in (crests, troughs)
        allpts = mirror ? vcat([(x, -y, z) for (x, y, z) in reverse(pts)], pts) : pts
        xs = [p[1] for p in allpts]; ys = [p[2] for p in allpts]; zs = [p[3] for p in allpts]
        push!(frame, (xs, ys, zs))
    end
    return frame
end

"""
    vtail_frame(vtail; root_offset=0.0, z_stations=nothing, symmetric=true, dat_path=DEFAULT_AIRFOIL_DAT)

Complete wireframe for the vertical tail: flat LE/TE planform outline (un-mirrored, `layout`'s
spanwise axis mapped to height `z` starting from `root_offset`, no dihedral), real airfoil
cross-sections at a handful of height stations (default: root, break, tip), and two continuous
spanwise (height-wise) crest/trough lines through them -- un-mirrored, single root->tip sweep,
unlike the wing/htail, and the "thickness" direction here is lateral `y`, not vertical `z` (LE/TE
spanwise lines skipped for the same reason as `wing_frame`/`htail_frame`: the flat outline
already traces them exactly). `symmetric` defaults to `true` -- a real vertical tail is almost
always a symmetric section, unlike the cambered wing/htail.
"""
function vtail_frame(vtail; root_offset::Float64=0.0, z_stations::Union{Nothing,Vector{Float64}}=nothing,
    symmetric::Bool=true, dat_path::String=DEFAULT_AIRFOIL_DAT)

    layout = vtail.layout
    λs = vtail.inboard.λ
    λt = vtail.outboard.λ
    bo, b, ηs = layout.root_span, layout.span, layout.ηs

    xw, yw = wing_planform_outline(layout, λs, λt; mirrored=false)
    frame = GeomLines()
    push!(frame, (xw, zeros(length(xw)), yw .+ root_offset))

    stations = z_stations !== nothing ? z_stations : [bo, ηs * b, b]
    crests = NTuple{3,Float64}[]; troughs = NTuple{3,Float64}[]
    for h in stations
        chord, toc, xLE = wing_local_chord_toc(vtail, h; mirrored=false)
        z = root_offset + h
        xu, yu, xl, yl = airfoil_section(chord, toc; dat_path=dat_path, xLE=xLE, yref=0.0, symmetric=symmetric)
        push!(frame, (xu, yu, fill(z, length(xu))))
        push!(frame, (xl, yl, fill(z, length(xl))))
        _, _, crest, trough = airfoil_key_points(xu, yu, xl, yl)
        push!(crests, (crest[1], crest[2], z))
        push!(troughs, (trough[1], trough[2], z))
    end

    for pts in (crests, troughs)
        xs = [p[1] for p in pts]; ys = [p[2] for p in pts]; zs = [p[3] for p in pts]
        push!(frame, (xs, ys, zs))
    end
    return frame
end

#### Engines -------------------------------------------------------------------

"""
    engine_frame(parg, wing; dihedral_angle=0.0, n_theta=10, z_root_override=nothing)

Wireframe of the wing-mounted engine nacelles (aft BLI fan not modeled yet). Reuses the exact
spanwise/chordwise engine-placement math from `stickfig`. The nacelle is drawn with its top
flush against the wing's dihedral-extended `z` reference at that station (no pylon-drop gap
modeled, matching TASOPT's own landing-gear clearance formula, which makes the same
flush-against-`z` assumption -- see `size_landing_gear.jl`). `z_root_override`, same convention
as `wing_frame`, replaces `wing.layout.z` as the root reference (see that docstring for why
`aircraft_frame` now passes `-Rfuse` here, defining `z` as the wing's leading-edge height,
flush with the fuselage belly at the root).
"""
function engine_frame(parg, wing; dihedral_angle::Float64=0.0, n_theta::Int=10,
    z_root_override::Union{Nothing,Float64}=nothing)
    D = parg[igdfan]
    neng = parg[igneng]
    lnace = parg[iglnace]

    co = wing.layout.root_chord
    bo = wing.layout.root_span
    b  = wing.layout.span
    ηs = wing.layout.ηs
    λs = wing.inboard.λ
    λt = wing.outboard.λ
    dy = 2 * D

    if neng == 2
        yi = [ηs * b / 2]
    else
        yi = collect(range(bo / 2 + dy, b / 2 * 3 / 4, length=Int(neng / 2)))
    end

    ηi = yi ./ (b / 2)
    ηo = bo / b
    ci = similar(yi)
    for (i, η) in enumerate(ηi)
        ci[i] = η <= ηs ? co * (1 + (λs - 1) * (η - ηo) / (ηs - ηo)) :
                           co * (λs + (λt - λs) * (η - ηs) / (1 - ηs))
    end

    tanL = tan(wing.layout.sweep * π / 180.0)
    xi = @. tanL * (yi - bo / 2) - 0.4 * ci + wing.layout.box_x - 1.0

    z_root = z_root_override !== nothing ? z_root_override : wing.layout.z

    frame = GeomLines()
    for (y_i, x_i) in zip(yi, xi)
        z_wing_local = z_root + y_i * tan(dihedral_angle)
        chord_i, toc_i, _ = wing_local_chord_toc(wing, y_i)
        #z is the wing's leading-edge height (see aircraft_frame); the engine actually hangs
        #~40% chord aft of the LE, so drop an extra half local thickness to approximate the
        #true lower-surface height there (drawing-only, not a TASOPT quantity).
        z_center = z_wing_local - chord_i * toc_i / 2 - D / 2

        xs_nace = [x_i, x_i + lnace]
        rs_nace = [D / 2, D / 2]

        append!(frame, axisym_body_frame_x(xs_nace, rs_nace; y0=y_i, z0=z_center, n_theta=n_theta, n_rings=2))
        append!(frame, axisym_body_frame_x(xs_nace, rs_nace; y0=-y_i, z0=z_center, n_theta=n_theta, n_rings=2))
    end
    return frame
end

#### Landing gear ----------------------------------------------------------------

#Real 737 MAX tire sizes/spacing (Boeing D6-38A004 Rev K, July 2025, Sec. 7.2.1,
#literature/737MAX_LandingGearSizes.pdf): nose "27x7.75R15" = 27in diameter x 7.75in width, dual
#wheels at 1ft4in track; main "H44.5x16.5R21" = 44.5in diameter x 16.5in width, dual wheels per
#strut at 2ft10in track. Both wheel pairs sit side-by-side in y (the footprint diagram's two
#dimensions are both measured against the same fuselage-centerline axis, same orientation for
#nose and main -- not fore-aft tandem, a bogie layout the 737 doesn't actually use). Six wheels
#total (2 nose + 2x2 main) -- real reference hardware used purely for a plausible drawing, same
#spirit as the b737b airfoil, not a TASOPT-derived value.
const NOSE_WHEEL_DIAMETER = 27 * in_to_m
const NOSE_WHEEL_WIDTH = 7.75 * in_to_m
const NOSE_WHEEL_TRACK = (1 * 12 + 4) * in_to_m
const MAIN_WHEEL_DIAMETER = 44.5 * in_to_m
const MAIN_WHEEL_WIDTH = 16.5 * in_to_m
const MAIN_WHEEL_TRACK = (2 * 12 + 10) * in_to_m

#The strut itself is approximated as a cylinder whose diameter equals the gap between the two
#wheels it connects (track minus one wheel's width) -- the strut/axle housing that extends from
#between the wheels up into the fuselage, not the thin generic placeholder radius used before
#real wheel data was available.
const NOSE_STRUT_DIAMETER = NOSE_WHEEL_TRACK - NOSE_WHEEL_WIDTH
const MAIN_STRUT_DIAMETER = MAIN_WHEEL_TRACK - MAIN_WHEEL_WIDTH

"""
    landing_gear_frame(landing_gear, xend, Rfuse, Deng, yeng, span, wing_z, dihedral_angle; nose_extra_z=0.0, kwargs...)

Wireframe of 9 cylinders: 3 landing-gear struts (`axisym_body_frame_z`, axis vertical along z)
and 6 real-sized wheels (`axisym_body_frame_y`, axis along y -- a wheel's axle really does point
spanwise, so its rolling/diameter profile is visible in the side-view x-z plane, not as a flat
circle lying in the ground plane). The nose strut ends in 2 side-by-side (y) wheels, each main
strut (2 of them) ends in another 2 side-by-side (y) wheels -- not fore-aft tandem, a bogie
layout the 737 doesn't actually use (the footprint diagram measures both wheel-pairs' spacing
against the same fuselage-centerline axis). Each wheel cylinder's center sits at
`ground_z + radius` so its bottom, not its center, is flush with the ground. Each strut stops at
that same `ground_z + radius` height (the wheel's own axle height) rather than continuing down
to `ground_z` itself, so it doesn't visually rub through the bottom of the tire.

Each strut's own diameter is approximated as the gap between the two wheels it connects (track
minus one wheel's width -- `NOSE_STRUT_DIAMETER`/`MAIN_STRUT_DIAMETER`), representing the
axle/strut housing that extends from between the wheels up into the fuselage, rather than a thin
generic placeholder.

Diameters/widths/spacing default to real 737 MAX values (`NOSE_WHEEL_DIAMETER` etc., sourced
from Boeing D6-38A004 -- see the comment above those constants) -- real reference hardware used
purely for a plausible drawing, not a TASOPT-derived value (TASOPT doesn't model wheel
count/size/geometry at all; `size_landing_gear.jl` only tracks `number_struts`/`wheels_per_strut`
*counts* for the historical-correlations weight formula). If
`landing_gear.model == "historical_correlations"`, the real sized strut lengths are used.
Otherwise (e.g. "mass_fractions", which only gives weight, not geometry) the exact same
clearance/tailstrike formula from `size_landing_gear.jl` is reapplied here, purely for drawing,
using the TOML defaults (`tailstrike_angle=10°`, `engine_ground_clearance=23in`,
`y_main_landing_gear_halfspan_fraction=0.2`) as fixed plotting-only assumptions. Nothing is
written back to the aircraft.

`dihedral_angle` (radians) is **not** re-derived here -- it must be passed in already resolved,
the same value given to `wing_frame`/`engine_frame`, so all three stay consistent by
construction (see `aircraft_frame`'s single `dihedral_angle_override` hook).

Both struts' lengths are measured from the dihedral-raised wing height at their own spanwise
attachment point (`wing_z` for the nose strut, `wing_z + y_main*tan(dihedral)` for the main
struts) -- that's what makes both struts reach the *same* flat ground, and makes the engine's
drawn clearance above that ground come out equal to `ground_clearance` exactly, consistent with
how `size_landing_gear.jl`'s formula was derived.

`nose_extra_z` is a uniform shift applied only to the nose strut's top (not the main struts,
which sit nowhere near the nose), meant for `aircraft_frame`'s `nose_droop_frac` option --
unlike htail's rigid-shift caveat, there's no tilt risk here since a strut is a single vertical
line, not a chord.

`wing_thickness_at_yeng` (the local absolute wing thickness -- chord*t/c -- at the engine's
spanwise station) is added as a half-thickness margin onto `l_clearance`'s fallback formula, to
account for `z` now meaning the wing's leading-edge height rather than its lower surface (see
`aircraft_frame`'s docstring) -- the engine actually hangs an extra half-thickness below that
reference. Leave at `0.0` to reproduce the plain formula.
"""
function landing_gear_frame(landing_gear, xend::Float64, Rfuse::Float64,
    Deng::Float64, yeng::Float64, span::Float64, wing_z::Float64, dihedral_angle::Float64;
    nose_extra_z::Float64=0.0, wing_thickness_at_yeng::Float64=0.0,
    fallback_tailstrike_deg::Float64=10.0,
    fallback_clearance_m::Float64=23 * in_to_m, fallback_y_frac::Float64=0.2,
    n_theta::Int=8,
    nose_wheel_diameter::Float64=NOSE_WHEEL_DIAMETER, nose_wheel_width::Float64=NOSE_WHEEL_WIDTH,
    nose_wheel_track::Float64=NOSE_WHEEL_TRACK, nose_strut_diameter::Float64=NOSE_STRUT_DIAMETER,
    main_wheel_diameter::Float64=MAIN_WHEEL_DIAMETER, main_wheel_width::Float64=MAIN_WHEEL_WIDTH,
    main_wheel_track::Float64=MAIN_WHEEL_TRACK, main_strut_diameter::Float64=MAIN_STRUT_DIAMETER)

    model_has_geometry = TASOPT.compare_strings(landing_gear.model, "historical_correlations")

    tailstrike_angle   = model_has_geometry ? landing_gear.tailstrike_angle       : fallback_tailstrike_deg * π / 180
    ground_clearance   = model_has_geometry ? landing_gear.engine_ground_clearance : fallback_clearance_m
    y_frac             = model_has_geometry ? landing_gear.main_gear.y_offset_halfspan_fraction : fallback_y_frac

    x_main = landing_gear.main_gear.weight.x
    x_nose = landing_gear.nose_gear.weight.x

    if model_has_geometry && landing_gear.nose_gear.length > 0.0
        lgnose_length = landing_gear.nose_gear.length
        lgmain_length = landing_gear.main_gear.length
        y_main        = landing_gear.main_gear.weight.y
    else
        l_tailstrike  = (xend - x_main) * tan(tailstrike_angle) - 2 * Rfuse
        #+wing_thickness_at_yeng/2: z is now the wing's *leading-edge* height (see aircraft_frame),
        #not the lower surface, so the engine (mounted ~40% chord aft of the LE) actually hangs an
        #extra half-thickness below that reference -- drawing-only, not a change to TASOPT's own formula.
        l_clearance   = ground_clearance + Deng - yeng * tan(dihedral_angle) + wing_thickness_at_yeng / 2
        lgnose_length = max(l_clearance, l_tailstrike)
        y_main        = span / 2 * y_frac
        lgmain_length = lgnose_length + y_main * tan(dihedral_angle)
    end

    #Ground is a single fixed reference, independent of nose_extra_z: both struts' *bottoms* must
    #stay flush on the same real ground regardless of the (cosmetic-only) nose-droop top shift --
    #only each strut's top follows its own attachment point; the strut is drawn shorter/longer to
    #compensate, never translated as a rigid body (same reasoning as the gear-length derivation
    #itself: ground_z = wing_z - lgnose_length is exactly what both branches already agree on).
    ground_z = wing_z - lgnose_length
    z_nose_top = wing_z + nose_extra_z
    z_main_top = wing_z + y_main * tan(dihedral_angle)

    frame = GeomLines()

    #Nose gear: single strut on centerline (diameter = the gap between the two wheels it
    #connects, extending up into the fuselage), dual wheels side-by-side at the bottom. Wheels
    #are true cylinders with axle along y (axisym_body_frame_y), tire bottom flush with the
    #ground (center at ground_z + radius, not ground_z itself -- that would bury half the tire).
    #The strut itself stops at the wheel's axle height (ground_z + radius), not the ground, so it
    #doesn't visually rub through the bottom of the tire.
    r_nose_wheel = nose_wheel_diameter / 2
    append!(frame, axisym_body_frame_z([ground_z + r_nose_wheel, z_nose_top], [nose_strut_diameter / 2, nose_strut_diameter / 2];
        x0=x_nose, y0=0.0, n_theta=n_theta, n_rings=2))
    for sgn in (1.0, -1.0)
        y_c = sgn * nose_wheel_track / 2
        ys_wheel = [y_c - nose_wheel_width / 2, y_c + nose_wheel_width / 2]
        append!(frame, axisym_body_frame_y(ys_wheel, [r_nose_wheel, r_nose_wheel];
            x0=x_nose, z0=ground_z + r_nose_wheel, n_theta=n_theta, n_rings=2))
    end

    #Main gear: two struts (symmetric about centerline), each with dual wheels side-by-side in y
    #(not fore-aft tandem -- the 737's footprint diagram measures both wheel-pair spacings
    #against the same fuselage-centerline axis, same orientation as the nose)
    r_main_wheel = main_wheel_diameter / 2
    for sgn in (1.0, -1.0)
        y_m = sgn * y_main
        append!(frame, axisym_body_frame_z([ground_z + r_main_wheel, z_main_top], [main_strut_diameter / 2, main_strut_diameter / 2];
            x0=x_main, y0=y_m, n_theta=n_theta, n_rings=2))
        for xsgn in (1.0, -1.0)
            y_c = y_m + xsgn * main_wheel_track / 2
            ys_wheel = [y_c - main_wheel_width / 2, y_c + main_wheel_width / 2]
            append!(frame, axisym_body_frame_y(ys_wheel, [r_main_wheel, r_main_wheel];
                x0=x_main, z0=ground_z + r_main_wheel, n_theta=n_theta, n_rings=2))
        end
    end
    return frame
end

#### Top-level orchestration and view projection ------------------------------

"""
    aircraft_frame(ac; dihedral_angle_override=nothing, upswept_tailcone=false, nose_droop_frac=0.0)

Collect every component's wireframe into one `Dict{Symbol,GeomLines}`. Single-bubble fuselage
aircraft only.

The wing dihedral angle used for drawing the wing itself, the htail, the engine height, and the
landing gear length/position is resolved **once**, here, and passed into
`wing_frame`/`engine_frame`/`landing_gear_frame` unchanged -- they never re-derive it
themselves, so there's a single source of truth and no risk of them drifting apart. By default
that's `landing_gear.wing_dihedral_angle` when `landing_gear.model == "historical_correlations"`,
or the TOML's own `6°` default otherwise (TASOPT has no real wing dihedral anywhere else, and no
htail dihedral at all -- see `landing_gear.jl`'s own `#TODO this should be a wing param`
comment; applying the same borrowed value to htail is a drawing-only choice, not a TASOPT
model). Pass `dihedral_angle_override` (radians) to try a different value everywhere at once.

`upswept_tailcone` bends the fuselage's tail-cone center axis upward
(`fuselage_tailcone_z_bend`, visual-only) and, consistently, rigidly shifts the whole htail up by
that same bend evaluated once at `htail.layout.box_x` (a uniform shift, not a per-point one --
see `htail_frame`'s docstring for why). The vtail needs no equivalent shift: it's anchored to the
tail cone's *top* surface, which stays flat at `Rfuse` by construction.

`nose_droop_frac` linearly droops the fuselage's nose center axis down by up to
`nose_droop_frac*Rfuse` at the very tip (`fuselage_nose_z_droop`, also visual-only), and
consistently shifts the nose landing gear strut's top down by that same bend evaluated once at
`nose_gear.weight.x` -- a uniform shift, same idea as the htail's, but with no tilt concern since
a strut is a single vertical line, not a chord.

`wing.layout.z` is undocumented anywhere in real TASOPT (used only in the landing-gear clearance
formula and this drawing code -- confirmed by searching the whole source tree). For the wing
specifically (not htail/vtail, which keep their own `layout.z` conventions), it's defined here as
the leading-edge height, with the dihedral extension's root (`y=0`) value set flush with the
fuselage's own true belly line at the wing's `box_x` station -- `-(Rfuse + fuselage_skirt_height(...))`,
not just `-Rfuse`, now that the fuselage is drawn with its real circle+skirt cross-section rather
than a plain circle. This is a drawing-only convention choice, not a recovered TASOPT fact, and
it replaces `wing.layout.z` entirely for the wing's planform outline, its airfoil sections, the
engine nacelles, and both landing-gear struts, all driven from this one `wing_z_root` value for
consistency.
"""
function aircraft_frame(ac; dihedral_angle_override::Union{Nothing,Float64}=nothing,
    upswept_tailcone::Bool=false, nose_droop_frac::Float64=0.0)
    fuse   = ac.fuselage
    wing   = ac.wing
    htail  = ac.htail
    vtail  = ac.vtail
    parg   = ac.parg
    landing_gear = ac.landing_gear

    Rfuse = fuse.layout.radius
    model_has_geometry = TASOPT.compare_strings(landing_gear.model, "historical_correlations")
    dihedral_angle = if dihedral_angle_override !== nothing
        dihedral_angle_override
    else
        model_has_geometry ? landing_gear.wing_dihedral_angle : 6.0 * π / 180
    end

    #Vertical tail root starts at the fuselage skin (constant Rfuse, single-bubble cylinder section)
    vtail_root_offset = Rfuse + vtail.layout.z
    htail_extra_z = upswept_tailcone ? fuselage_tailcone_z_bend(fuse, htail.layout.box_x) : 0.0
    nose_gear_extra_z = nose_droop_frac != 0.0 ? fuselage_nose_z_droop(fuse, landing_gear.nose_gear.weight.x, nose_droop_frac) : 0.0
    #wing.layout.z is otherwise undocumented anywhere in TASOPT; for the wing specifically we
    #define it here as the leading-edge height, with the dihedral extension's root (y=0) value
    #flush with the fuselage's own true belly (circle + skirt, not just -Rfuse) at the wing's own
    #x station -- a drawing-only convention choice, not a recovered fact.
    wing_z_root = -(Rfuse + fuselage_skirt_height(fuse, wing.layout.box_x))
    yeng = parg[igyeng]
    chord_eng, toc_eng, _ = wing_local_chord_toc(wing, yeng)
    wing_thickness_at_yeng = chord_eng * toc_eng

    components = Dict{Symbol,GeomLines}()
    components[:fuselage] = fuselage_frame(fuse; upswept_tailcone=upswept_tailcone, nose_droop_frac=nose_droop_frac)
    components[:wing]     = wing_frame(wing; dihedral_angle=dihedral_angle, z_root_override=wing_z_root)
    components[:htail]    = htail_frame(htail; dihedral_angle=dihedral_angle, extra_z=htail_extra_z) #htail keeps its own layout.z convention, unchanged
    components[:vtail]    = vtail_frame(vtail; root_offset=vtail_root_offset)
    components[:engines]  = engine_frame(parg, wing; dihedral_angle=dihedral_angle, z_root_override=wing_z_root)
    components[:landing_gear] = landing_gear_frame(landing_gear, fuse.layout.x_end,
        Rfuse, parg[igdfan], yeng, wing.layout.span, wing_z_root, dihedral_angle;
        nose_extra_z=nose_gear_extra_z, wing_thickness_at_yeng=wing_thickness_at_yeng)
    return components
end

"""
    project_view(components; view=:top)

Flatten an `aircraft_frame(ac)` dict into a list of 2D `(u, v)` polylines for `view`
in `(:top, :side, :front)`.
"""
function project_view(components::Dict{Symbol,GeomLines}; view::Symbol=:top)
    lines2d = Tuple{Vector{Float64},Vector{Float64}}[]
    for frame in values(components), (xs, ys, zs) in frame
        if view == :top
            push!(lines2d, (xs, ys))
        elseif view == :side
            push!(lines2d, (xs, zs))
        elseif view == :front
            push!(lines2d, (ys, zs))
        else
            error("Unknown view $view, expected :top, :side, or :front")
        end
    end
    return lines2d
end
