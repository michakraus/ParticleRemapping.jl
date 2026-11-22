# The parametrisation of the remap Hamiltonian by tensor-product cubic B-splines in (ξ, v, ln w),
# d = 1: the kick–drift map and the non-separable implicit-midpoint map, with their parameter
# vectors θ, stage by stage.

"""
    AbstractRemapParametrisation

A family of symplectic maps `R_θ` of 1D-1V markers with a parameter vector `θ`:
[`KickDriftParametrisation`](@ref) and [`NonSeparableParametrisation`](@ref). Each implements
[`remap_stages`](@ref) and [`remap_parameter_count`](@ref).
"""
abstract type AbstractRemapParametrisation end

# The `n` cubic splines on `[min ln w, max ln w]`; for equal weights on an interval of unit length
# from `ln w`, so that the stages are the same for every marker.
function _log_weight_basis(w, n)
    lw = log.(w)
    lo, hi = extrema(lw)
    UniformBasis(lo, (hi > lo ? hi - lo : one(hi)) / (n - 1), n, false)
end

"""
    KickDriftParametrisation(w; variant = :weight_dependent, L = 2π, V = 4, nξ = 32, nv = 40,
        nw = 16, stages = 2)

The kick–drift map `R_θ = D_S ∘ K_S ∘ … ∘ D_1 ∘ K_1` of `stages` stages for markers with the
weights `w`, with `nξ` periodic cubic splines `B_i` on `[0, L)`, `nv` cubic splines `B̂_j` on
`[−V, V]` and `nw` cubic splines `C_l` on `[min ln w, max ln w]` (for equal weights on
`[ln w, ln w + 1]`). The `variant` selects the
generating functions:

- `:weight_dependent`: `κ = w Σ λ_il B_i(ξ) C_l(ln w)`, `η = w Σ μ_jl B̂_j(v) C_l(ln w)`;
- `:control_A`: `κ = w̄ Σ λ_i B_i(ξ)`, `η = w̄ Σ μ_j B̂_j(p/(m w̄))` with the mean weight `w̄`, no
  dependence on the weight in `(x, p)`;
- `:control_B`: `κ = w Σ λ_i B_i(ξ)`, `η = w Σ μ_j B̂_j(v)`, one area-preserving map of `(x, v)`.

The parameters of stage `s` are `(vec(λ^s), vec(μ^s))`, stage after stage.
"""
struct KickDriftParametrisation{T <: Real} <: AbstractRemapParametrisation
    variant::Symbol
    Bξ::UniformBasis{T}
    Bv::UniformBasis{T}
    Bl::UniformBasis{T}
    wbar::T
    stages::Int
end

function KickDriftParametrisation(w::AbstractVector; variant::Symbol = :weight_dependent,
        L = 2π, V = 4.0, nξ = 32, nv = 40, nw = 16, stages = 2)
    variant in (:weight_dependent, :control_A, :control_B) ||
        throw(ArgumentError("unknown kick–drift variant $variant"))
    Bl = _log_weight_basis(w, nw)
    KickDriftParametrisation(variant, UniformBasis(0.0, L / nξ, nξ, true),
        UniformBasis(-V, 2V / (nv - 1), nv, false), Bl, float(sum(w) / length(w)), Int(stages))
end

_nw(par::KickDriftParametrisation) = par.variant === :weight_dependent ? length(par.Bl) : 1

"""
    remap_parameter_count(par::AbstractRemapParametrisation)

The length of the parameter vector `θ` of the map `par`.
"""
function remap_parameter_count(par::KickDriftParametrisation)
    par.stages * (length(par.Bξ) + length(par.Bv)) * _nw(par)
end

# Σ c[i, l] B_i'(s) C_l(ℓ) over the nonzero basis functions; Bl === nothing for one constant C.
function _dsum(c, Bs, s, Bl, ℓ)
    acc = zero(eltype(c)) * zero(s)
    eachnonzero(Bs, s) do i, _, dB
        if Bl === nothing
            acc += c[i, 1] * dB
        else
            eachnonzero(Bl, ℓ) do l, C, _
                acc += c[i, l] * dB * C
            end
        end
    end
    acc
end

# The kick and drift gradient functions gκ(ξ, w), gη(v, w) of one stage with coefficients λ, μ.
function _stage_functions(par::KickDriftParametrisation, λ, μ)
    Bξ, Bv, Bl, wbar = par.Bξ, par.Bv, par.Bl, par.wbar
    if par.variant === :weight_dependent
        gκ = (ξ, w) -> [w * _dsum(λ, Bξ, ξ[1], Bl, log(w))]
        gη = (v, w) -> [w * _dsum(μ, Bv, v[1], Bl, log(w))]
    elseif par.variant === :control_A
        gκ = (ξ, w) -> [wbar * _dsum(λ, Bξ, ξ[1], nothing, nothing)]
        gη = (v, w) -> [w * _dsum(μ, Bv, v[1] * w / wbar, nothing, nothing)]
    else
        gκ = (ξ, w) -> [w * _dsum(λ, Bξ, ξ[1], nothing, nothing)]
        gη = (v, w) -> [w * _dsum(μ, Bv, v[1], nothing, nothing)]
    end
    gκ, gη
end

"""
    kick_drift_map(par::KickDriftParametrisation, θ, x, p, w; m = 1, history = nothing)
        -> (X, P)

The image `R_θ(x, p)` of the markers with lifted positions `x`, momenta `p` and weights `w`
(vectors) under the kick–drift map `par` ([`centroid_kick`](@ref) and [`velocity_drift`](@ref)).
A vector `history` receives `(X, P)` after each kick and each drift.
"""
function kick_drift_map(par::KickDriftParametrisation, θ, x, p, w; m = 1, history = nothing)
    X = reshape(x, 1, :)
    P = reshape(p, 1, :) .+ zero(eltype(θ))
    nk, nd = length(par.Bξ) * _nw(par), length(par.Bv) * _nw(par)
    off = 0
    for _ in 1:(par.stages)
        λ = reshape(view(θ, (off + 1):(off + nk)), length(par.Bξ), :)
        μ = reshape(view(θ, (off + nk + 1):(off + nk + nd)), length(par.Bv), :)
        off += nk + nd
        gκ, gη = _stage_functions(par, λ, μ)
        P = centroid_kick(X, P, w, gκ)
        history === nothing || push!(history, (vec(copy(X)), vec(copy(P))))
        X = velocity_drift(X, P, w, m, gη)
        history === nothing || push!(history, (vec(copy(X)), vec(copy(P))))
    end
    vec(X), vec(P)
end

"""
    NonSeparableParametrisation(w; L = 2π, V = 4, nξ = 16, nv = 17, nℓ = 6, stages = 2,
        weighted = true, iteration_limit = 50)

The map of `stages` non-separable implicit-midpoint stages for markers with the weights `w`, with
`φ_s(ξ, v, ℓ) = Σ λ^s_ijl B_i(ξ) B̂_j(v) C_l(ℓ) / m`, `ℓ = ln w`, so that `ψ_s = m w φ_s`; `nξ`
periodic cubic splines on `[0, L)`, `nv` on `[−V, V]` and `nℓ` on `[min ln w, max ln w]` (for
equal weights on `[ln w, ln w + 1]`). With
`weighted = false` (control C) every marker evaluates `φ_s` at the mean of `ln w`, which makes the
map one area-preserving map of `(x, v)`. A stage whose fixed-point iteration needs more than
`iteration_limit` iterations makes the map return `NaN` positions, so that a fit rejects such
parameters: the limit keeps the contraction rate of the iteration below about 0.55. The parameters
of stage `s` are `vec(λ^s)` with the linear index of `(i, j, l)`.
"""
struct NonSeparableParametrisation{T <: Real} <: AbstractRemapParametrisation
    Bξ::UniformBasis{T}
    Bv::UniformBasis{T}
    Bℓ::UniformBasis{T}
    stages::Int
    weighted::Bool
    iteration_limit::Int
end

function NonSeparableParametrisation(w::AbstractVector; L = 2π, V = 4.0, nξ = 16, nv = 17,
        nℓ = 6, stages = 2, weighted::Bool = true, iteration_limit = 50)
    Bℓ = _log_weight_basis(w, nℓ)
    NonSeparableParametrisation(UniformBasis(0.0, L / nξ, nξ, true),
        UniformBasis(-V, 2V / (nv - 1), nv, false), Bℓ, Int(stages), weighted,
        Int(iteration_limit))
end

function _stage_count(par::NonSeparableParametrisation)
    length(par.Bξ) * length(par.Bv) * length(par.Bℓ)
end
remap_parameter_count(par::NonSeparableParametrisation) = par.stages * _stage_count(par)

# The terms (linear index of λ_ijl, ∂_v, ∂_ξ of B_i(ξ) B̂_j(v) C_l(ℓ) / m) that do not vanish.
function _φterms(par::NonSeparableParametrisation, ξ, v, ℓ, m)
    T = promote_type(typeof(ξ), typeof(v))
    terms = Tuple{Int, T, T}[]
    lin = LinearIndices((length(par.Bξ), length(par.Bv), length(par.Bℓ)))
    eachnonzero(par.Bξ, ξ) do i, Bi, dBi
        eachnonzero(par.Bv, v) do j, Bj, dBj
            eachnonzero(par.Bℓ, ℓ) do l, Cl, _
                push!(terms, (lin[i, j, l], Bi * dBj * Cl / m, dBi * Bj * Cl / m))
            end
        end
    end
    terms
end

# (∂_v φ, ∂_ξ φ) for the coefficients λ of one stage.
function _φslopes(par::NonSeparableParametrisation, λ, ξ, v, ℓ, m)
    a = b = zero(promote_type(eltype(λ), typeof(ξ), typeof(v)))
    for (q, ∂a, ∂b) in _φterms(par, ξ, v, ℓ, m)
        a += λ[q] * ∂a
        b += λ[q] * ∂b
    end
    a, b
end

"""
    nonseparable_map(par::NonSeparableParametrisation, θ, x, p, w; m = 1,
        sensitivities = false, centroid_terms = true, maxit = par.iteration_limit,
        history = nothing) -> (; X, P, iterations, sens)

The image of the markers with lifted positions `x`, momenta `p` and weights `w` (vectors) under
the stages of `par`, each solved by [`nonseparable_stage`](@ref) with at most `maxit` iterations;
`iterations` is the largest number that a stage needed. Without convergence `X` and `P` are `NaN`.

With `sensitivities = true`, `sens = (; dξ, dv, dXc)` holds the derivatives of the shifted
positions `ξ_k`, of the velocities `v_k` and of the centroid `X_c` with respect to `θ` (`dξ`, `dv`
are `length(θ) × N`), from the linearised step equations: one `2 × 2` system per marker and one
`2 × 2` system for the shift `c` and the boost `b`. With `centroid_terms = false` the terms of `c`
and `b` are dropped, the case that must fail. A vector `history` receives `(X, P)` after each
stage.
"""
function nonseparable_map(par::NonSeparableParametrisation, θ, x, p, w; m = 1,
        sensitivities = false, centroid_terms = true, maxit = par.iteration_limit,
        history = nothing)
    Nm, nθ = length(x), length(θ)
    nθs = _stage_count(par)
    α = w ./ sum(w)
    ℓfixed = sum(log, w) / Nm
    ℓof(ω) = par.weighted ? log(ω) : ℓfixed
    ℓ = ℓof.(w)
    X, P = reshape(x, 1, :), reshape(p, 1, :)
    dξ, dv, dXc = zeros(nθ, Nm), zeros(nθ, Nm), zeros(nθ)
    Sg = Diagonal([1.0, -1.0])
    iterations = 0
    for s in 1:(par.stages)
        cols = ((s - 1) * nθs + 1):(s * nθs)
        λ = view(θ, cols)
        gφ = (ξ, v, w) -> begin
            a, b = _φslopes(par, λ, ξ[1], v[1], ℓof(w), m)
            ([b], [a])
        end
        X1, P1, its, ok = nonseparable_stage(X, P, w, m, gφ; maxit)
        ok || return (;
            X = fill(NaN, Nm), P = fill(NaN, Nm), iterations = its, sens = nothing)
        iterations = max(iterations, its)
        if sensitivities
            ξ0, ξ1 = vec(X) .- dot(α, X), vec(X1) .- dot(α, X1)
            v0, v1 = vec(P) ./ (m .* w), vec(P1) ./ (m .* w)
            dξn, dvn = similar(dξ), similar(dv)
            rhsξ, rhsv = zeros(nθ), zeros(nθ)
            T = zeros(2, 2)
            AiS = Vector{Matrix{Float64}}(undef, Nm)
            for k in 1:Nm
                ξb, vb = (ξ0[k] + ξ1[k]) / 2, (v0[k] + v1[k]) / 2
                H = ForwardDiff.jacobian(
                    y -> collect(_φslopes(par, λ, y[1], y[2], ℓ[k], m)),
                    [ξb, vb])
                Ai = inv(I - Sg * H / 2)
                C = Ai * (I + Sg * H / 2)
                AiS[k] = Ai * Sg
                a, b = view(dξ, :, k), view(dv, :, k)
                dξn[:, k] .= C[1, 1] .* a .+ C[1, 2] .* b
                dvn[:, k] .= C[2, 1] .* a .+ C[2, 2] .* b
                for (q, ∂a, ∂b) in _φterms(par, ξb, vb, ℓ[k], m)
                    dξn[cols[q], k] += AiS[k][1, 1] * ∂a + AiS[k][1, 2] * ∂b
                    dvn[cols[q], k] += AiS[k][2, 1] * ∂a + AiS[k][2, 2] * ∂b
                    rhsξ[cols[q]] += α[k] * ∂a
                    rhsv[cols[q]] += α[k] * ∂b
                end
                ma, mb = (a .+ view(dξn, :, k)) ./ 2, (b .+ view(dvn, :, k)) ./ 2
                rhsξ .+= α[k] .* (H[1, 1] .* ma .+ H[1, 2] .* mb)
                rhsv .+= α[k] .* (H[2, 1] .* ma .+ H[2, 2] .* mb)
                T .+= α[k] .* (H * AiS[k]) ./ 2
            end
            dg = centroid_terms ? (I + T) \ permutedims(hcat(rhsξ, rhsv)) : zeros(2, nθ)
            for k in 1:Nm
                dξn[:, k] .-= AiS[k][1, 1] .* dg[1, :] .+ AiS[k][1, 2] .* dg[2, :]
                dvn[:, k] .-= AiS[k][2, 1] .* dg[1, :] .+ AiS[k][2, 2] .* dg[2, :]
            end
            dξ, dv = dξn, dvn
            dXc .+= dg[1, :]
        end
        X, P = X1, P1
        history === nothing || push!(history, (vec(copy(X)), vec(copy(P))))
    end
    (; X = vec(X), P = vec(P), iterations, sens = sensitivities ? (; dξ, dv, dXc) : nothing)
end

"""
    remap_stages(par::AbstractRemapParametrisation, θ, x, p, w; m = 1, history = nothing)
        -> (X, P)

The image `R_θ(x, p)`: [`kick_drift_map`](@ref) or the positions and momenta of
[`nonseparable_map`](@ref).
"""
function remap_stages(par::KickDriftParametrisation, θ, x, p, w; m = 1, history = nothing)
    kick_drift_map(par, θ, x, p, w; m, history)
end

function remap_stages(
        par::NonSeparableParametrisation, θ, x, p, w; m = 1, history = nothing)
    r = nonseparable_map(par, θ, x, p, w; m, history)
    r.X, r.P
end

"""
    nonseparable_sensitivity_error(par, θ, x, p, w; m = 1, centroid_terms = true, ndir = 3,
        maxit = 200, rng = Random.default_rng())

The largest relative deviation of the sensitivities of [`nonseparable_map`](@ref) from forward-mode
derivatives along `ndir` random directions of `θ`. The partials of the fixed-point iteration
converge more slowly than its values, so the check allows `maxit` iterations.
"""
function nonseparable_sensitivity_error(
        par::NonSeparableParametrisation, θ, x, p, w; m = 1,
        centroid_terms = true, ndir = 3, maxit = 200, rng::AbstractRNG = Random.default_rng())
    N = length(x)
    s = nonseparable_map(par, θ, x, p, w; m, sensitivities = true, centroid_terms, maxit).sens
    err = zero(float(eltype(x)))
    for _ in 1:ndir
        δ = randn(rng, length(θ))
        d = ForwardDiff.derivative(
            t -> begin
                r = nonseparable_map(par, θ .+ t .* δ, x, p, w; m, maxit)
                vcat(r.X, r.P)
            end, 0.0)
        dX, dV = d[1:N], d[(N + 1):end] ./ (m .* w)
        err = max(
            err, maximum(abs, dX .- permutedims(s.dξ .+ s.dXc) * δ) / maximum(abs, dX),
            maximum(abs, dV .- permutedims(s.dv) * δ) / maximum(abs, dV))
    end
    err
end
