# The canonical weight squeeze T_μ: the linear cotangent lift (q, p) ↦ (q_anc + M⁻ᵀ (q − q_anc), M p)
# with M = diag(μ) (I + (1/N) 1 (1 − μ)ᵀ), in O(N) operations by the Sherman–Morrison formula.

"""
    canonical_weight_squeeze(q, p, μ, q_anc) -> (q̄, p̄)

The momentum-conserving canonical weight squeeze with positive slot multipliers `μ`,
`Σ_a μ_a = N`:

```
p̄_a = μ_a (p_a + (1/N) Σ_c (1 − μ_c) p_c),
q̄_a = q_anc,a + ((q_a − q_anc,a) − (1 − μ_a) (1/N) Σ_c (q_c − q_anc,c)) / μ_a .
```

It is the linear cotangent lift `p̄ = M p`, `q̄ − q_anc = M⁻ᵀ (q − q_anc)` of
`M = diag(μ) (I + (1/N) 1 (1 − μ)ᵀ)`, whose inverse transpose the Sherman–Morrison formula gives
in closed form because `Σ μ = N`. It keeps `Σ dq ∧ dp` and `Σ p`. The anchors `q_anc` are
constants. Throws an `ArgumentError` unless `Σ μ = N` to round-off.
"""
function canonical_weight_squeeze(q::AbstractVector, p::AbstractVector, μ::AbstractVector,
        q_anc::AbstractVector{<:Real})
    N = length(q)
    abs(sum(μ) - N) <= sqrt(eps(Float64)) * N ||
        throw(ArgumentError("the slot multipliers must sum to N = $N, not $(sum(μ))"))
    TR = promote_type(eltype(q), eltype(p), eltype(μ))
    p_corr = zero(TR)
    q_corr = zero(TR)
    for c in 1:N
        p_corr += (one(TR) - TR(μ[c])) * TR(p[c])
        q_corr += TR(q[c]) - TR(q_anc[c])
    end
    p_corr /= N
    q_corr /= N
    pbar = Vector{TR}(undef, N)
    qbar = Vector{TR}(undef, N)
    for a in 1:N
        ma = TR(μ[a])
        pbar[a] = ma * (TR(p[a]) + p_corr)
        qbar[a] = TR(q_anc[a]) + ((TR(q[a]) - TR(q_anc[a])) - (one(TR) - ma) * q_corr) / ma
    end
    qbar, pbar
end
