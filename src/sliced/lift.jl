# The sliced quantile optimal-transport lift: the gradient of a strictly convex potential S(u) and
# its cotangent lift, and the dual lift that acts on the momenta.

# The attenuation factor κ of a slice with targets τ and mollifier width ε: the positive root of
# κ² + ρ κ − 1 = 0 with ρ = ε / (half the target span). It tends to 1 as ε → 0.
function _attenuation(taus::AbstractVector, ε::Real)
    ρ = ε / max((taus[end] - taus[1]) / 2, 1e-12)
    sqrt(4.0 + ρ^2) / 2.0 - ρ / 2.0
end

"""
    sliced_quantile_potential(u, bins, taus; β_mean = 0.60, β_intra = 0.72, ε = 0.4)

The strictly convex potential of the sliced quantile lift,

```
S(u) = (1/2N) (Σ_c u_c)² + (1 − β_mean)/2 Σ_b n_b (ū_b − ū)² + (1 − β_intra)/2 Σ_b Σ_{a∈b} (u_a − ū_b)²
     + Σ_b β_intra / (2 n_b κ_b) Σ_{r,s} |τ_r − τ_s| √((u_{a_r} − u_{a_s})² + ε²),
```

with `ū` the mean of all `u`, `ū_b` the mean over slice `b` (the indices `bins[b]`, `n_b` of
them, in the order of the targets `taus[b]`) and `κ_b` the attenuation factor of the slice. Its
gradient is the map of [`split_quantile_ot_lift`](@ref).
"""
function sliced_quantile_potential(u::AbstractVector{T}, bins, taus; β_mean::Real = 0.60,
        β_intra::Real = 0.72, ε::Real = 0.4) where {T}
    N = length(u)
    ū = sum(u) / N
    S = sum(u)^2 / (2N)
    for (b, idx) in pairs(bins)
        nb = length(idx)
        ūb = sum(u[idx]) / nb
        S += (1 - β_mean) / 2 * nb * (ūb - ū)^2
        S += (1 - β_intra) / 2 * sum(abs2, u[idx] .- ūb)
        τ = taus[b]
        c = β_intra / (2nb * _attenuation(τ, ε))
        for r in 1:nb, s in 1:nb

            S += c * abs(τ[r] - τ[s]) * sqrt((u[idx[r]] - u[idx[s]])^2 + ε^2)
        end
    end
    S
end

"""
    split_quantile_ot_lift(u, v, bins, taus; anchors = nothing, β_mean = 0.60, β_intra = 0.72,
                           ε = 0.4, with_cotangent = true) -> (ũ, ṽ)

The anchored cotangent lift of the gradient map `ũ = ∇S(u)` of the
[`sliced_quantile_potential`](@ref) `S`:

```
ũ_a = β_mean ū + (β_intra − β_mean) ū_b + (1 − β_intra) u_a
      + β_intra / (n_b κ_b) Σ_{s} |τ_r − τ_s| (u_a − u_{a_s}) / √((u_a − u_{a_s})² + ε²),
ṽ   = v_anc + (∇²S(u))⁻¹ (v − v_anc),
```

for marker `a = a_r` of slice `b`. Within each slice the ordered markers relax towards the
targets `taus[b]`; the slice means relax towards the mean by the factor `β_mean`, the deviations
within a slice by `β_intra`. The Hessian `∇²S` is symmetric with `∇²S ⪰ (1 − max(β_mean, β_intra)) I`
and `∇²S 1 = 1`, so the lift is a symplectic diffeomorphism of `(u, v)` that keeps `Σ u` and
`Σ v`. The anchors `v_anc` (zero without `anchors`) are constants.

`(u, v)` is `(q, p)` for the lift in the positions and `(p, q)` for the [`dual_quantile_ot_lift`](@ref).
With `with_cotangent = false` the conjugate variable is returned unchanged: the map is then not
symplectic (a control).
"""
function split_quantile_ot_lift(u::AbstractVector, v_conj::AbstractVector,
        bins::AbstractVector{<:AbstractVector{<:Integer}},
        taus::AbstractVector{<:AbstractVector{<:Real}};
        anchors::Union{Nothing, AbstractVector{<:Real}} = nothing, β_mean::Real = 0.60,
        β_intra::Real = 0.72, ε::Real = 0.4, with_cotangent::Bool = true)
    T = promote_type(eltype(u), eltype(v_conj))
    unew, H = _gradient_hessian(T, u, bins, taus, β_mean, β_intra, ε)
    with_cotangent || return unew, Vector{T}(v_conj)
    anchors === nothing && return unew, H \ Vector{T}(v_conj)
    unew, (H \ (v_conj .- T.(anchors))) .+ T.(anchors)
end

"""
    sliced_quantile_hessian(u, bins, taus; β_mean = 0.60, β_intra = 0.72, ε = 0.4)

The Hessian `∇²S(u)` of the [`sliced_quantile_potential`](@ref): the Jacobian of the map
`u ↦ ∇S(u)` of [`split_quantile_ot_lift`](@ref), a symmetric matrix with
`∇²S ⪰ (1 − max(β_mean, β_intra)) I` and `∇²S 1 = 1`.
"""
function sliced_quantile_hessian(u::AbstractVector, bins, taus; β_mean::Real = 0.60,
        β_intra::Real = 0.72, ε::Real = 0.4)
    last(_gradient_hessian(eltype(u), u, bins, taus, β_mean, β_intra, ε))
end

# The gradient ∇S(u) and the Hessian ∇²S(u) of the sliced quantile potential, in element type T.
function _gradient_hessian(::Type{T}, u, bins, taus, β_mean, β_intra, ε) where {T}
    N = length(u)
    ubar = sum(u) / N
    unew = Vector{T}(undef, N)
    H = fill(T(β_mean) / N, N, N)
    eps2 = T(ε)^2
    for (b, idx) in pairs(bins)
        nb = length(idx)
        ubar_b = sum(u[idx]) / nb
        bm_diff = T(β_intra - β_mean) / nb
        for ia in 1:nb, ic in 1:nb

            H[idx[ia], idx[ic]] += bm_diff
        end
        τ = taus[b]
        coef = T(β_intra / _attenuation(τ, ε)) / nb
        for r in 1:nb
            a = idx[r]
            ua = u[a]
            rep_a = zero(T)
            diag_h = zero(T)
            for s in 1:nb
                r == s && continue
                c = idx[s]
                wt = T(abs(τ[r] - τ[s]))
                diff_u = ua - u[c]
                denom = sqrt(diff_u^2 + eps2)
                rep_a += wt * diff_u / denom
                hac = -coef * wt * eps2 / (denom * denom * denom)
                H[a, c] += hac
                diag_h -= hac
            end
            unew[a] = T(β_mean) * ubar + T(β_intra - β_mean) * ubar_b +
                      (one(T) - T(β_intra)) * ua + coef * rep_a
            H[a, a] += (one(T) - T(β_intra)) + diag_h
        end
    end
    unew, H
end

"""
    dual_quantile_ot_lift(q, p, bins, taus; anchors = nothing, kwargs...) -> (q̂, p̂)

The dual lift: [`split_quantile_ot_lift`](@ref) with the roles of positions and momenta
exchanged, `p̂ = ∇S(p)`, `q̂ = q_anc + (∇²S(p))⁻¹ (q − q_anc)`. It equalises the momenta within
each slice and is symplectic for the same reasons, since `(q, p) ↦ (p, −q)` is.
"""
function dual_quantile_ot_lift(q::AbstractVector, p::AbstractVector, bins, taus; kwargs...)
    p̂, q̂ = split_quantile_ot_lift(p, q, bins, taus; kwargs...)
    q̂, p̂
end
