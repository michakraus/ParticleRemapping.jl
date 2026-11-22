# The charge-conserving reweighting, d = 1: the weights after the map that minimise the relative
# entropy to the target weights under the constraint of the grid charge.

"""
    spline_grid_charge(b::UniformBasis, x, w)

The grid charge `ρ_i = Σ_k w_k Λ_i(x_k)` of markers at the positions `x` with the weights `w` on
the periodic cubic B-splines `b`, for the charge `q = 1` of a marker of unit weight.
"""
function spline_grid_charge(b::UniformBasis, x, w)
    ρ = zeros(eltype(w), b.n)
    for k in eachindex(x, w)
        eachnonzero((i, Λ, _) -> (ρ[i] += w[k] * Λ), b, x[k])
    end
    ρ
end

"""
    spline_charge_matrix(b::UniformBasis, x)

The matrix `A` with `A[i, k] = Λ_i(x_k)`, so that the grid charge is `A w`.
"""
function spline_charge_matrix(b::UniformBasis, x)
    A = zeros(b.n, length(x))
    for k in eachindex(x)
        eachnonzero((i, Λ, _) -> (A[i, k] += Λ), b, x[k])
    end
    A
end

"""
    charge_reweight(b::UniformBasis, x, wstar, ρ; tol = 1e-13, maxit = 50)
        -> (; w, λ, iterations, converged)

The weights `w_k = w*_k exp(Σ_i λ_i Λ_i(x_k))` with the grid charge `A w = ρ`
([`spline_charge_matrix`](@ref)). The multipliers `λ` minimise the convex dual
`Φ(λ) = Σ_k w_k(λ) − λᵀρ`; Newton's method with a backtracking (Armijo) line search on `Φ` solves
for them. The line search accepts a step whose change of `Φ` lies within the rounding error of
`Φ`, so that the full Newton step is taken near the solution. The iteration stops when `|A w − ρ| < tol max(ρ)`, and reports `converged = false` when
the Hessian `A diag(w) Aᵀ` is not positive definite (a shape function without markers) or the
iteration does not converge in `maxit` iterations.
"""
function charge_reweight(b::UniformBasis, x, wstar, ρ; tol = 1e-13, maxit = 50)
    A = spline_charge_matrix(b, x)
    λ = zeros(b.n)
    weights(λ) = wstar .* exp.(A' * λ)
    Φ(λ) = sum(weights(λ)) - dot(λ, ρ)
    for it in 1:maxit
        w = weights(λ)
        F = A * w - ρ
        maximum(abs, F) < tol * maximum(ρ) &&
            return (; w, λ, iterations = it - 1, converged = true)
        H = Symmetric(A * Diagonal(w) * A')
        isposdef(H) || return (; w, λ, iterations = it - 1, converged = false)
        δ = -(H \ F)
        s, Φ0 = 1.0, Φ(λ)
        # the rounding error of Φ: near the solution the decrease of a Newton step falls below it
        ε = 8 * eps(sum(w) + abs(dot(λ, ρ)))
        while Φ(λ .+ s .* δ) > Φ0 + 1e-4 * s * dot(F, δ) + ε && s > 1e-12
            s /= 2
        end
        λ = λ .+ s .* δ
    end
    w = weights(λ)
    (; w, λ, iterations = maxit, converged = maximum(abs, A * w - ρ) < tol * maximum(ρ))
end

"""
    additive_charge_correction(b::UniformBasis, x, wstar, ρ)

The least-norm additive correction `w = w* + Aᵀ(AAᵀ)⁻¹(ρ − A w*)`. It reproduces the grid charge
`ρ` but can give negative weights where small target weights meet an excess of charge: the case
that must fail against [`charge_reweight`](@ref).
"""
function additive_charge_correction(b::UniformBasis, x, wstar, ρ)
    A = spline_charge_matrix(b, x)
    wstar .+ A' * ((A * A') \ (ρ .- A * wstar))
end

"""
    ChargeConservingReweighting(basis::UniformBasis; tol = 1e-13, maxit = 50)
    ChargeConservingReweighting(L, n; kwargs...)

The charge-conserving reweighting on the periodic cubic B-splines `basis` (or `n` splines on
`[0, L)`): the weights after the map that minimise the relative entropy to the target weights
under the constraint of the grid charge, by [`charge_reweight`](@ref).

    remap(method::ChargeConservingReweighting, markers::Markers; reference_charge)

The markers carry the positions after the map, the target weights `w*` and the velocities
`p / (m w*)`; the result keeps the positions and the momenta and replaces the weights by `w′` and
the velocities by `p / (m w′)`. `info` holds the multipliers `λ`, the Newton `iterations`,
`converged`, the relative grid-charge residuals of `w*` and of `w′` against `reference_charge`,
and `max_log_change = max_k |ln(w′_k / w*_k)|`. Without a solution (a shape function with charge
and no marker) it throws an `ErrorException`.
"""
struct ChargeConservingReweighting{T <: Real} <: AbstractRemap
    basis::UniformBasis{T}
    tol::T
    maxit::Int
end

function ChargeConservingReweighting(basis::UniformBasis; tol = 1e-13, maxit = 50)
    ChargeConservingReweighting(basis, oftype(basis.h, tol), Int(maxit))
end

function ChargeConservingReweighting(L::Real, n::Integer; kwargs...)
    ChargeConservingReweighting(UniformBasis(0.0, L / n, n, true); kwargs...)
end

function remap(method::ChargeConservingReweighting, mk::Markers; reference_charge)
    mk.x isa AbstractVector ||
        throw(ArgumentError("the charge-conserving reweighting acts on 1D-1V markers"))
    ρ = reference_charge
    b = method.basis
    p = momentum(mk)
    r = charge_reweight(b, mk.x, mk.w, ρ; tol = method.tol, maxit = method.maxit)
    r.converged || error("no positive weights reproduce the reference grid charge")
    residual(w) = maximum(abs, spline_grid_charge(b, mk.x, w) .- ρ) / maximum(ρ)
    info = (; λ = r.λ, iterations = r.iterations, converged = r.converged,
        charge_residual_target = residual(mk.w), charge_residual = residual(r.w),
        max_log_change = maximum(abs, log.(r.w ./ mk.w)))
    RemapResult(Markers(copy(mk.x), velocity(p, r.w, mk.m), r.w; m = mk.m), info)
end
