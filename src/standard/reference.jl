# The reference of the standard test set: a resample of the markers onto a lattice of the box,
# weighted by a density. It keeps the total weight only; it is not symplectic.

using .Projection: hat_interpolate, lumped_projection

"""
    GridResample(; L, V, projection = HatProjection(16, 17), density = nothing)

The reference remap of the standard test set on the box `[−L/2, L/2) × [−V, V]`.
`remap(method, mk)` places the `n × n` cell centres `(−L/2 + (i − 0.5) dx, −V + (j − 0.5) dv)`,
`dx = L/n`, `dv = 2V/n`, `i` outer and `j` inner, for `N = n²` markers, and gives each the
weight `c f(x_i, v_j) dx dv`, with `c = Σ w / Σ f dx dv` so that the total weight is kept, and
the mass of `mk`. With `density = nothing`, `f` is the lumped projection of `mk` on the hats
`projection` of the box (a marker with `|v| > V` lands on the edge nodes, by the clamp of the
hats); otherwise `f = density(x, v)`. A lattice point where `f` vanishes keeps its marker with
the weight zero.

The `info` holds the `variant` (`:lumped` or `:exact`), `raw_total = Σ f dx dv`, the
`scale = c` and the number `zero_weights` of zero weights. An `ArgumentError` refuses `d × N`
markers and an `N` that is no square `n²` with `n ≥ 1`; a `DomainError` refuses a density that
is not positive at any lattice point.
"""
struct GridResample{T <: Real, D} <: AbstractRemap
    L::T
    V::T
    projection::HatProjection{T}
    density::D
end

function GridResample(;
        L::Real, V::Real, projection::HatProjection = HatProjection(16, 17),
        density = nothing)
    GridResample(float(L), float(V), projection, density)
end

"""
    REFERENCE_METHODS

The keys of the two variants of the reference [`GridResample`](@ref) for
[`case_method`](@ref): `:grid_lumped` (the lumped projection of the markers) and `:grid_exact`
(the exact density of the case).
"""
const REFERENCE_METHODS = (:grid_lumped, :grid_exact)

function remap(method::GridResample, mk::Markers)
    mk.x isa AbstractVector || throw(ArgumentError("GridResample acts on 1D-1V markers"))
    N = length(mk)
    n = _lattice_side(N)
    n === nothing && throw(ArgumentError(
        "GridResample places an n × n lattice, so N must be a square n^2 >= 1, got N = $N"))
    L, V, hp = method.L, method.V, method.projection
    dx, dv = L / n, 2V / n
    x, v = _cell_centres(L, V, n)
    f = if method.density === nothing
        g = lumped_projection(hp, mk.x ./ L, mk.v ./ V, mk.w)
        [hat_interpolate(hp, g, x[k] / L, v[k] / V) / (L * V) for k in eachindex(x)]
    else
        method.density.(x, v)
    end
    raw = f .* (dx * dv)
    raw_total = sum(raw)
    raw_total > 0 || throw(DomainError(raw_total,
        "GridResample needs a density that is positive at some lattice point"))
    scale = sum(mk.w) / raw_total
    w = raw .* scale
    RemapResult(Markers(x, v, w; m = mk.m),
        (; variant = method.density === nothing ? :lumped : :exact, raw_total, scale,
            zero_weights = count(iszero, w)))
end

export GridResample, REFERENCE_METHODS
