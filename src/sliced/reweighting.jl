# The density ratio that sets the slot multipliers μ of the weight squeeze and the final weights.
#
# The kernel is Gaussian in v and periodic in x with period 2Lx:
#     K(Δx, Δv) = exp(−(Σ_ℓ d(Δx_ℓ)² / hx² + Σ_ℓ Δv_ℓ² / hv²) / 2),
#     d(Δx)² = 2 s² (1 − cos(Δx / s)),  s = Lx / π,
# which equals Δx² to leading order and is exactly periodic. Positions and velocities are vectors
# (1D-1V) or d × N matrices.

# The squared periodic distance 2 s² (1 − cos(Δ / s)) with s = Lx / π.
_periodic_distance2(Δ, s) = s^2 * (2.0 * (1.0 - cos(Δ / s)))

# A vector of N values as a 1 × N matrix; a d × N matrix as it is.
_axes_by_markers(x::AbstractVector) = reshape(x, 1, :)
_axes_by_markers(x::AbstractMatrix) = x

# One value per marker, broadcastable against x: itself for a vector x, a row for a matrix x.
_per_marker(μ::AbstractVector, ::AbstractVector) = μ
_per_marker(μ::AbstractVector, ::AbstractMatrix) = reshape(μ, 1, :)

# The exponent of the kernel between marker a of (X, V) and marker c of (Y, U).
function _kernel_exponent(X, V, a, Y, U, c, s, inv_hx2, inv_hv2, ::Type{T}) where {T}
    e = zero(T)
    for ℓ in axes(X, 1)
        e += _periodic_distance2(X[ℓ, a] - Y[ℓ, c], s) * inv_hx2
    end
    for ℓ in axes(V, 1)
        e += (V[ℓ, a] - U[ℓ, c])^2 * inv_hv2
    end
    e
end

"""
    kernel_density_ratio(xe, ve, x, v, w; hx = 0.45, hv = 0.35, Lx = π)

The ratio `r_a = Σ_c w_c K(xe_a − x_c, ve_a − v_c) / Σ_c K(xe_a − xe_c, ve_a − ve_c)` at the
evaluation points `(xe, ve)`: a kernel estimate of the density `f` of the weighted markers
`(x, v, w)` divided by a kernel estimate of the density of the evaluation points themselves. It is
the weight that a marker at `(xe_a, ve_a)` carries when the evaluation points sample `f`. The
kernel is Gaussian with widths `hx`, `hv` and periodic in each position component with period
`2Lx`. Positions and velocities are vectors (1D-1V) or `d × N` matrices.
"""
function kernel_density_ratio(
        xe::AbstractVecOrMat, ve::AbstractVecOrMat, x::AbstractVecOrMat,
        v::AbstractVecOrMat, w::AbstractVector; hx::Real = 0.45, hv::Real = 0.35,
        Lx::Real = π)
    Xe, Ve, X, V = _axes_by_markers.((xe, ve, x, v))
    s = Lx / π
    inv_hx2 = 1.0 / hx^2
    inv_hv2 = 1.0 / hv^2
    T = promote_type(eltype(xe), eltype(ve), eltype(x), eltype(v), eltype(w), Float64)
    r = Vector{T}(undef, size(Xe, 2))
    for a in axes(Xe, 2)
        fn = zero(T)
        rd = zero(T)
        for c in axes(X, 2)
            e = _kernel_exponent(Xe, Ve, a, X, V, c, s, inv_hx2, inv_hv2, T)
            fn += w[c] * exp(-e / 2)
        end
        for c in axes(Xe, 2)
            e = _kernel_exponent(Xe, Ve, a, Xe, Ve, c, s, inv_hx2, inv_hv2, T)
            rd += exp(-e / 2)
        end
        r[a] = fn / rd
    end
    r
end

# The weight squeeze along each axis of d × N positions and momenta, with one common μ.
function _squeeze_axes(q̂::AbstractVector, p̂::AbstractVector, μ, q_anc::AbstractVector)
    canonical_weight_squeeze(q̂, p̂, μ, q_anc)
end

function _squeeze_axes(q̂::AbstractMatrix, p̂::AbstractMatrix, μ, q_anc::AbstractMatrix)
    rows = [canonical_weight_squeeze(q̂[ℓ, :], p̂[ℓ, :], μ, q_anc[ℓ, :]) for ℓ in axes(q̂, 1)]
    permutedims(reduce(hcat, first.(rows))), permutedims(reduce(hcat, last.(rows)))
end

"""
    squeeze_multipliers(q̂, p̂, q_anc, x, v, w; m, w_ref = sum(w) / length(w), hx = 0.45,
                        hv = 0.35, ratio_floor = 0.08, iterations = 3, Lx = π)

The slot multipliers `μ` (`Σ μ = N`) of the [`canonical_weight_squeeze`](@ref) from the density
ratio of the input markers `(x, v, w)`, by a predictor–corrector loop. Start from the points
`(q̂, p̂ / (m w_ref))` and repeat `iterations` times:

1. `μ_a ∝ r_a + ratio_floor · w_ref` with `r` the [`kernel_density_ratio`](@ref) at the current
   points, normalised to `Σ μ = N`;
2. squeeze `(q̂, p̂)` with these `μ` (along each axis for `d × N` arrays) and evaluate at the
   squeezed points `(q̄, p̄ / (m w_ref μ))`.

The result depends on the state, so a squeeze with it is canonical for frozen `μ` only.
"""
function squeeze_multipliers(q̂::AbstractVecOrMat, p̂::AbstractVecOrMat,
        q_anc::AbstractVecOrMat, x::AbstractVecOrMat, v::AbstractVecOrMat, w::AbstractVector;
        m::Real, w_ref::Real = sum(w) / length(w), hx::Real = 0.45, hv::Real = 0.35,
        ratio_floor::Real = 0.08, iterations::Integer = 3, Lx::Real = π)
    N = length(w)
    μ = ones(N)
    q_eval = copy(q̂)
    v_eval = p̂ ./ (m .* w_ref)
    for _ in 1:iterations
        μ_raw = kernel_density_ratio(q_eval, v_eval, x, v, w; hx = hx, hv = hv, Lx = Lx) .+
                ratio_floor * w_ref
        μ = (N / sum(μ_raw)) .* μ_raw
        q_eval, p_eval = _squeeze_axes(q̂, p̂, μ, q_anc)
        v_eval = p_eval ./ (m .* (w_ref .* _per_marker(μ, p_eval)))
    end
    μ
end
