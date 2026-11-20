# Diagnostics of the hole-filling experiments, d = 1: the largest empty disk, the fill of a hole,
# the discrepancy of the marker density, and the errors of target weights by velocity band.

"""
    EllipticHole(x, v, a, b; L = 2π)

The elliptic hole `((x − x_h)/a)² + ((v − v_h)/b)² < 1` around `(x_h, v_h) = (x, v)` in
`[0, L) × ℝ`, with positions taken modulo `L`.
"""
struct EllipticHole{T <: Real}
    x::T
    v::T
    a::T
    b::T
    L::T
end

function EllipticHole(x::Real, v::Real, a::Real, b::Real; L::Real = 2π)
    EllipticHole(promote(
        float(x), v, a, b, L)...)
end

"""
    inside_hole(h::EllipticHole, x, v)

Whether the point `(x, v)` lies in the hole.
"""
inside_hole(h::EllipticHole, x, v) = ((mod(x, h.L) - h.x) / h.a)^2 + ((v - h.v) / h.b)^2 < 1

"""
    hole_outline(h::EllipticHole; n = 200) -> (x, v)

`n` points of the boundary of the hole.
"""
function hole_outline(h::EllipticHole; n = 200)
    t = range(0, 2π; length = n)
    h.x .+ h.a .* cos.(t), h.v .+ h.b .* sin.(t)
end

"""
    hole_fill(h::EllipticHole, x, v; volume)

The number of markers in the hole relative to an equidistribution of the `length(x)` markers over
the phase-space `volume`: `#{k : (x_k, v_k) ∈ hole} / (π a b N / volume)`.
"""
function hole_fill(h::EllipticHole, x, v; volume)
    count(k -> inside_hole(h, x[k], v[k]), eachindex(x)) /
    (π * h.a * h.b * length(x) / volume)
end

"""
    largest_empty_disk(x, v; L, dx, dv, xprobes = range(0, L; length = 241)[1:(end - 1)],
        vprobes = range(-3.5, 3.5; length = 201))

The largest distance from a probe point to its nearest marker, in the units `dx`, `dv` of the
lattice spacing (positions modulo `L`, minimum-image distance in `x`): the radius of the largest
empty disk among the probes. A perfect lattice gives `√2/2 ≈ 0.71`.
"""
function largest_empty_disk(
        x, v; L, dx, dv, xprobes = range(0, L; length = 241)[1:(end - 1)],
        vprobes = range(-3.5, 3.5; length = 201))
    xs = mod.(x, L)
    worst = zero(float(eltype(xs)))
    for px in xprobes, pv in vprobes

        dmin = convert(typeof(worst), Inf)
        for k in eachindex(xs)
            ddx = abs(px - xs[k])
            ddx = min(ddx, L - ddx) / dx
            ddv = (pv - v[k]) / dv
            dmin = min(dmin, ddx^2 + ddv^2)
        end
        worst = max(worst, dmin)
    end
    sqrt(worst)
end

"""
    local_empty_disk(h::EllipticHole, x, v; dx, dv, half_widths = (1.0, 0.8), n = (81, 65))

[`largest_empty_disk`](@ref) on `n` probes of the box `x_h ± half_widths[1]`,
`v_h ± half_widths[2]` around the hole.
"""
function local_empty_disk(
        h::EllipticHole, x, v; dx, dv, half_widths = (1.0, 0.8), n = (81, 65))
    L = h.L
    xs = mod.(x, L)
    worst = zero(float(eltype(xs)))
    for px in range(h.x - half_widths[1], h.x + half_widths[1]; length = n[1]),
        pv in range(h.v - half_widths[2], h.v + half_widths[2]; length = n[2])

        d = minimum(
            k -> (min(abs(px - xs[k]), L - abs(px - xs[k])) / dx)^2 + ((pv - v[k]) / dv)^2,
            eachindex(xs))
        worst = max(worst, d)
    end
    sqrt(worst)
end

"""
    density_discrepancy(hp::HatProjection, x, v; L, V)

The root mean square deviation of the lumped marker density `g_h` ([`lumped_projection`](@ref))
from its mean `N/2` in scaled variables, relative to the mean, on the nodes away from the
velocity boundary (two node rows dropped at each end).
"""
function density_discrepancy(hp::HatProjection, x, v; L, V)
    g = lumped_projection(hp, x ./ L, v ./ V, ones(length(x)))
    gbar = length(x) / 2
    inner = g[:, 3:(end - 2)]
    sqrt(sum(abs2, inner .- gbar) / length(inner)) / gbar
end

"""
    weight_errors_by_band(v, w, wref; bands = ((0, 1), (1, 2), (2, 3), (3, 4)))

The root mean square and the largest absolute value of the relative error `w_k / wref_k − 1` of
the markers with `|v_k|` in each band `[a, b)`, as a tuple of `(; rms, max)`.
"""
function weight_errors_by_band(v, w, wref; bands = ((0, 1), (1, 2), (2, 3), (3, 4)))
    rel = w ./ wref .- 1
    map(bands) do band
        idx = findall(k -> band[1] <= abs(v[k]) < band[2], eachindex(v))
        (rms = sqrt(sum(abs2, rel[idx]) / length(idx)), max = maximum(abs, rel[idx]))
    end
end
