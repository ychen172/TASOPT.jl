# Scale a real, named reference airfoil shape to a target chord and thickness-to-chord ratio,
# for cross-sectional visualization only. This is NOT the airfoil TASOPT itself uses for drag
# (that's an undocumented bespoke MSES family designed by Drela, see src/airfoil_data/README) --
# it's a real published airfoil used purely as a plausible shape for plotting a wing section.
#
# Default reference: Boeing 737 midspan airfoil ("b737b-il"), UIUC Airfoil Coordinates Database
# (https://m-selig.ae.illinois.edu/ads/coord/b737b.dat), Lednicer format.
# Max thickness 12.5% at 29.7% chord. Max camber 0.8% at 10% chord.

const DEFAULT_AIRFOIL_DAT = joinpath(@__DIR__, "airfoil_data", "b737b.dat")
const REF_AIRFOIL_MAX_TOC = 0.125 #Reference airfoil's own documented max thickness/chord

"""
    read_airfoil_dat(path=DEFAULT_AIRFOIL_DAT) -> (xu, yu, xl, yl)

Parse a Lednicer-format .dat airfoil coordinate file (name line, point-count line, blank line,
upper surface LE->TE, blank line, lower surface LE->TE). Coordinates are normalized to chord=1,
re-zeroed so the leading edge (the shared x=0 point, identical on both surfaces in the raw data)
sits exactly at y=0 -- the raw `b737b.dat` file's own LE is at y=0.0088, not 0, which would
otherwise put the drawn LE slightly off whatever height it's told to sit at.
"""
function read_airfoil_dat(path::String=DEFAULT_AIRFOIL_DAT)
    lines = readlines(path)
    nu, nl = (round(Int, parse(Float64, s)) for s in split(lines[2]))

    coords = Vector{Tuple{Float64,Float64}}[]
    cur = Tuple{Float64,Float64}[]
    for line in lines[3:end]
        s = strip(line)
        if isempty(s)
            if !isempty(cur)
                push!(coords, cur)
                cur = Tuple{Float64,Float64}[]
            end
            continue
        end
        x, y = parse.(Float64, split(s))
        push!(cur, (x, y))
    end
    !isempty(cur) && push!(coords, cur)
    @assert length(coords) == 2 "Expected 2 coordinate blocks (upper, lower surfaces), got $(length(coords))"
    upper, lower = coords
    @assert length(upper) == nu && length(lower) == nl "Point counts don't match header"

    xu = [p[1] for p in upper]; yu = [p[2] for p in upper]
    xl = [p[1] for p in lower]; yl = [p[2] for p in lower]

    @assert xu[1] == 0.0 && xl[1] == 0.0 && yu[1] == yl[1] "Expected both surfaces' first point to be the shared leading edge at x=0"
    y_LE = yu[1]
    yu .-= y_LE
    yl .-= y_LE
    return xu, yu, xl, yl
end

"""
    symmetrize_airfoil(xu, yu, xl, yl) -> (xu, yu, xl, yl)

Strip the camber, keeping only the thickness envelope, so the result is a symmetric airfoil
with the same thickness distribution (half the envelope above y=0, half below). Used for the
vertical tail, which (unlike a cambered lifting wing/htail section) is almost always a symmetric
section in real aircraft. The lower surface is linearly interpolated onto the upper surface's
own x-grid so both outputs share one x-grid.
"""
function symmetrize_airfoil(xu::Vector{Float64}, yu::Vector{Float64}, xl::Vector{Float64}, yl::Vector{Float64})
    yl_on_xu = map(xu) do x
        j = clamp(searchsortedlast(xl, x), 1, length(xl) - 1)
        x0, x1 = xl[j], xl[j+1]
        t = x1 > x0 ? (x - x0) / (x1 - x0) : 0.0
        yl[j] + t * (yl[j+1] - yl[j])
    end
    halft = (yu .- yl_on_xu) ./ 2
    return xu, halft, xu, -halft
end

"""
    airfoil_section(chord, thickness_to_chord; dat_path=DEFAULT_AIRFOIL_DAT, xLE=0.0, yref=0.0, symmetric=false)

Scale the reference airfoil shape to a target `chord` [m] and maximum `thickness_to_chord`.
The thickness is matched by uniformly scaling the reference shape's y-coordinates by
`thickness_to_chord / REF_AIRFOIL_MAX_TOC` -- the chordwise (x) point distribution is kept as
the reference airfoil's. This also scales the apparent camber proportionally, which isn't
physically exact (real camber doesn't generally scale with thickness), but is a reasonable
simplification for a schematic cross-section. With `symmetric=true`, the camber is first
stripped via `symmetrize_airfoil` (same thickness envelope, zero camber) -- use this for the
vertical tail. Returns `(xu, yu, xl, yl)` in real units, with the leading edge at `(xLE, yref)`.
"""
function airfoil_section(chord::Float64, thickness_to_chord::Float64;
    dat_path::String=DEFAULT_AIRFOIL_DAT, xLE::Float64=0.0, yref::Float64=0.0, symmetric::Bool=false)

    xu, yu, xl, yl = read_airfoil_dat(dat_path)
    symmetric && ((xu, yu, xl, yl) = symmetrize_airfoil(xu, yu, xl, yl))
    scale = thickness_to_chord / REF_AIRFOIL_MAX_TOC

    xu_s = xLE .+ xu .* chord
    yu_s = yref .+ yu .* scale .* chord
    xl_s = xLE .+ xl .* chord
    yl_s = yref .+ yl .* scale .* chord
    return xu_s, yu_s, xl_s, yl_s
end

"""
    airfoil_key_points(xu, yu, xl, yl) -> (le, te, crest, trough)

Four named `(x, y)` points on an already-scaled/placed airfoil section (as returned by
`airfoil_section`): the leading edge (`xu[1]`, shared with `xl[1]`), the trailing edge (`xu[end]`,
with `y` averaged between the upper and lower TE -- they aren't identical in the reference
airfoil, which has a small nonzero TE thickness), the upper surface's highest point ("crest",
`argmax(yu)`), and the lower surface's lowest point ("trough", `argmin(yl)`). Used to build
continuous spanwise wireframe lines through a sequence of cross-sections, not just isolated
per-station outlines.
"""
function airfoil_key_points(xu::Vector{Float64}, yu::Vector{Float64}, xl::Vector{Float64}, yl::Vector{Float64})
    le = (xu[1], yu[1])
    te = (xu[end], (yu[end] + yl[end]) / 2)
    i_top = argmax(yu)
    crest = (xu[i_top], yu[i_top])
    i_bot = argmin(yl)
    trough = (xl[i_bot], yl[i_bot])
    return le, te, crest, trough
end

"""
    plot_airfoil_section(chord, thickness_to_chord; kwargs...)

Plot the scaled airfoil section (upper and lower surfaces). `kwargs` are passed through to
`airfoil_section`.
"""
function plot_airfoil_section(chord::Float64, thickness_to_chord::Float64; kwargs...)
    xu, yu, xl, yl = airfoil_section(chord, thickness_to_chord; kwargs...)
    p = plot(xu, yu, label="", color=:black, lw=1.5, aspect_ratio=:equal,
        xlabel="x [m]", ylabel="y [m]",
        title="Airfoil Section (c=$(round(chord,digits=2))m, t/c=$(round(thickness_to_chord*100,digits=1))%)")
    plot!(p, xl, yl, label="", color=:black, lw=1.5)
    return p
end
