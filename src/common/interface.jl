"""
    AbstractRemap

A redistribution of weighted markers. Each subtype implements [`remap`](@ref).
"""
abstract type AbstractRemap end

"""
    RemapResult(markers, info)

The markers after a remap and a `NamedTuple` of diagnostics (iterations, residuals, multipliers,
the homotopy parameter reached, …), whose fields depend on the method.
"""
struct RemapResult{M <: Markers, I <: NamedTuple}
    markers::M
    info::I
end

"""
    remap(method::AbstractRemap, markers::Markers; kwargs...) -> RemapResult

Apply the redistribution `method` to `markers`. The methods are listed in the documentation of
each construction: the gauge rezoning and the current-quantile matching of
[`Remapping`](@ref), the kick–drift and non-separable maps and the reweightings of
[`Projection`](@ref), the spline-projection redistribution of [`Redistribution`](@ref) and the
sliced lifts of [`SlicedRedistribution`](@ref).
"""
function remap end

"""
    conservation_report(before::Markers, after::Markers; grid = nothing)

The changes of the total weight, the total momentum and the kinetic energy between two marker
sets, and, for a [`PeriodicGrid`](@ref) `grid`, the largest changes of the grid moments
`M^(0)`, `M^(1)`, `M^(2)` relative to the corresponding moments of `|v|`.
"""
function conservation_report(before::Markers, after::Markers; grid = nothing)
    r = (weight = total_weight(after) - total_weight(before),
        momentum = norm(total_momentum(after) .- total_momentum(before)),
        kinetic_energy = kinetic_energy(after) - kinetic_energy(before))
    grid === nothing && return r
    rel(k) = maximum(abs.(grid_moment(grid, after, k) .- grid_moment(grid, before, k)) ./
                     grid_moment(grid, before.x, abs.(before.v), before.w, k))
    merge(r, (charge = rel(0), current = rel(1), energy = rel(2)))
end

export AbstractRemap, RemapResult, remap, conservation_report
