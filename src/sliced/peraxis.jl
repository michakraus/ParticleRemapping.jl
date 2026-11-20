# The sliced lifts in d > 1 dimensions, applied along each canonical axis.

"""
    per_axis_sliced_lift(x, p, bins, taus; kwargs...) -> (x̃, p̃)

The [`split_quantile_ot_lift`](@ref) of each axis `ℓ` of `d × N` positions and momenta,
`(x̃[ℓ, :], p̃[ℓ, :]) = Ψ_ℓ(x[ℓ, :], p[ℓ, :])`, with the slices `bins[ℓ]` and the targets `taus`
(the same for all axes). The pairs `(x_ℓ, p_ℓ)` are canonically conjugate and the lifts of
different axes act on different pairs, so the product is symplectic for `Σ_ℓ dx_ℓ ∧ dp_ℓ` and
keeps each component of `Σ p`. The keyword arguments go to every lift.
"""
function per_axis_sliced_lift(x::AbstractMatrix, p::AbstractMatrix, bins::AbstractVector,
        taus; kwargs...)
    size(x) == size(p) || throw(DimensionMismatch("x is $(size(x)), p is $(size(p))"))
    length(bins) == size(x, 1) ||
        throw(DimensionMismatch("$(size(x, 1)) axes but $(length(bins)) slicings"))
    rows = [split_quantile_ot_lift(x[ℓ, :], p[ℓ, :], bins[ℓ], taus; kwargs...)
            for ℓ in axes(x, 1)]
    permutedims(reduce(hcat, first.(rows))), permutedims(reduce(hcat, last.(rows)))
end

"""
    remap(method::SlicedQuantileRemap, mk::Markers{T, <:AbstractMatrix}) -> RemapResult

The sliced redistribution of `d × N` markers: for each axis `ℓ` the lifts `Ψ_q`, `Ψ_p` of the
pair `(x_ℓ, p_ℓ)`, sliced by `v_ℓ` and by `q̃_ℓ` (one [`SlicedLiftPlan`](@ref) per axis); then
one set of multipliers `μ` from the `2d`-dimensional [`kernel_density_ratio`](@ref), the
[`canonical_weight_squeeze`](@ref) of each axis with this `μ`, and `w̄ = w_ref μ`,
`v̄ = p̄ / (m w̄)`. Each lift sees the two-dimensional marginal `(x_ℓ, v_ℓ)` only: correlations
between different axes are not equalised. `info` holds the vector of plans and the intermediate
`d × N` states.
"""
function remap(method::SlicedQuantileRemap, mk::Markers{T, <:AbstractMatrix}) where {T}
    m = mk.m
    w_ref = sum(mk.w) / length(mk)
    p = momentum(mk)
    d = size(mk.x, 1)
    plans = [_plan(method, mk.x[ℓ, :], mk.v[ℓ, :], p[ℓ, :], m, w_ref) for ℓ in 1:d]
    lifted = [apply_sliced_lift(plans[ℓ], mk.x[ℓ, :], p[ℓ, :]) for ℓ in 1:d]
    by_axis(s) = permutedims(reduce(hcat, [getfield(l, s) for l in lifted]))
    q̃, p̃, q̂, p̂ = by_axis(:q̃), by_axis(:p̃), by_axis(:q̂), by_axis(:p̂)
    q_anc = permutedims(reduce(hcat, [plan.q_anc for plan in plans]))
    μ = squeeze_multipliers(q̂, p̂, q_anc, mk.x, mk.v, mk.w; m = m, w_ref = w_ref,
        hx = method.hx, hv = method.hv, ratio_floor = method.ratio_floor,
        iterations = method.iterations, Lx = method.Lx)
    q̄, p̄ = _squeeze_axes(q̂, p̂, μ, q_anc)
    w̄ = w_ref .* μ
    v̄ = p̄ ./ (m .* reshape(w̄, 1, :))
    RemapResult(Markers(q̄, v̄, w̄; m = m), (; q̃, p̃, q̂, p̂, q̄, p̄, μ, w_ref, plan = plans))
end
