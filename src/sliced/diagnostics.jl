# Diagnostics of the equidistribution of markers in a phase-space box [−Lx, Lx) × [−Vmax, Vmax].

"""
    periodic_wrap(x, Lx)

The position `x` mapped into `[−Lx, Lx)`: `mod(x + Lx, 2Lx) − Lx`.
"""
periodic_wrap(x::Real, Lx::Real) = mod(x + Lx, 2Lx) - Lx

"""
    phase_quadrant_counts(x, v) -> (n₁, n₂, n₃, n₄)

The numbers of markers in the four quadrants of the `(x, v)` plane, counterclockwise from
`x ≥ 0, v ≥ 0`: `Q1 = {x ≥ 0, v ≥ 0}`, `Q2 = {x < 0, v ≥ 0}`, `Q3 = {x < 0, v < 0}`,
`Q4 = {x ≥ 0, v < 0}`.
"""
function phase_quadrant_counts(x::AbstractVector, v::AbstractVector)
    (sum((x .>= 0) .& (v .>= 0)), sum((x .< 0) .& (v .>= 0)), sum((x .< 0) .& (v .< 0)),
        sum((x .>= 0) .& (v .< 0)))
end

"""
    phase_cell_statistics(x, v; gx = 10, gv = 10, Lx = π, Vmax = 2.25)
        -> (; empty, min_count, max_count, rstd, counts)

The marker counts of the `gx × gv` equal cells of `[−Lx, Lx) × [−Vmax, Vmax]` (positions wrapped
periodically; markers with `|v| > Vmax` are not counted), the number of empty cells, the smallest
and largest count, and the relative standard deviation `rstd` of the counts (standard deviation
over mean, with the mean over all cells).
"""
function phase_cell_statistics(x::AbstractVector, v::AbstractVector; gx::Integer = 10,
        gv::Integer = 10, Lx::Real = π, Vmax::Real = 2.25)
    counts = zeros(Int, gx, gv)
    for a in eachindex(x, v)
        ix = clamp(floor(Int, (periodic_wrap(x[a], Lx) + Lx) / (2Lx) * gx) + 1, 1, gx)
        if abs(v[a]) <= Vmax
            iv = clamp(floor(Int, (v[a] + Vmax) / (2Vmax) * gv) + 1, 1, gv)
            counts[ix, iv] += 1
        end
    end
    mean_count = sum(counts) / length(counts)
    rstd = sqrt(sum((counts .- mean_count) .^ 2) / length(counts)) / mean_count
    (; empty = count(==(0), counts),
        min_count = minimum(counts), max_count = maximum(counts),
        rstd, counts)
end

"""
    cell_count_history(snapshots, times; gx = 10, gv = 10, Lx = π, Vmax = 2.25)
        -> (; times, empty, rstd)

The number of empty cells and the relative standard deviation of the cell counts of
[`phase_cell_statistics`](@ref) for each marker set of `snapshots` (1D-1V `Markers`).
"""
function cell_count_history(snapshots::AbstractVector{<:Markers}, times::AbstractVector;
        gx::Integer = 10, gv::Integer = 10, Lx::Real = π, Vmax::Real = 2.25)
    stats = [phase_cell_statistics(mk.x, mk.v; gx = gx, gv = gv, Lx = Lx, Vmax = Vmax)
             for mk in snapshots]
    (; times = collect(times), empty = [s.empty for s in stats],
        rstd = [s.rstd for s in stats])
end

"""
    smoothed_density_difference(x₀, v₀, w₀, x₁, v₁, w₁; hx = 0.45, hv = 0.35, Lx = π,
                                Vmax = 2.25, ngrid = 30, margin = 0.1)

The relative `L²` difference `‖f₁ − f₀‖ / ‖f₀‖` of the kernel-smoothed densities
`f_k(x, v) = Σ_a w_a K(x − x_a, v − v_a)` of two weighted marker sets, on the `ngrid × ngrid`
points of `[−Lx + margin, Lx − margin] × [−Vmax, Vmax]`. The kernel is that of
[`kernel_density_ratio`](@ref), periodic in `x`.
"""
function smoothed_density_difference(x₀::AbstractVector, v₀::AbstractVector,
        w₀::AbstractVector, x₁::AbstractVector, v₁::AbstractVector, w₁::AbstractVector;
        hx::Real = 0.45, hv::Real = 0.35, Lx::Real = π, Vmax::Real = 2.25,
        ngrid::Integer = 30, margin::Real = 0.1)
    s = Lx / π
    inv_hx2 = 1.0 / hx^2
    inv_hv2 = 1.0 / hv^2
    T = promote_type(eltype(w₀), eltype(w₁), Float64)
    f_err = zero(T)
    f_nrm = zero(T)
    for xg in range(-Lx + margin, Lx - margin; length = ngrid),
        vg in range(-Vmax, Vmax; length = ngrid)

        fb = zero(T)
        fa = zero(T)
        for c in eachindex(x₀, v₀, w₀, x₁, v₁, w₁)
            dx1 = _periodic_distance2(xg - x₀[c], s)
            dx2 = _periodic_distance2(xg - x₁[c], s)
            fb += w₀[c] * exp(-(dx1 * inv_hx2 + (vg - v₀[c])^2 * inv_hv2) / 2)
            fa += w₁[c] * exp(-(dx2 * inv_hx2 + (vg - v₁[c])^2 * inv_hv2) / 2)
        end
        f_err += (fa - fb)^2
        f_nrm += fb^2
    end
    sqrt(f_err / f_nrm)
end
