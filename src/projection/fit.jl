# The fit of the parameters of the remap map to the target points: Levenberg–Marquardt on the
# regularised least-squares objective.

# The Jacobian of the residual: forward-mode for the explicit kick–drift stages, the linearised step
# equations for the non-separable stages.
function _residual_jacobian(
        par::KickDriftParametrisation, residual, θ, x, p, w, wt, m, ℓx, ℓv)
    ForwardDiff.jacobian(residual, θ)
end

function _residual_jacobian(
        par::NonSeparableParametrisation, residual, θ, x, p, w, wt, m, ℓx,
        ℓv)
    s = nonseparable_map(par, θ, x, p, w; m, sensitivities = true).sens
    vcat(permutedims(s.dξ .+ s.dXc) ./ ℓx, permutedims(s.dv) .* (w ./ (wt .* ℓv)))
end

"""
    fit_remap_parameters(par, x, p, w, xt, vt, wt; m = 1, ℓx, ℓv, β = 1e-3, maxit = 40,
        trace = false) -> (; θ, F0, F, iterations)

The parameters `θ` of the map `par` that minimise

    J(θ) = ½ Σ_k (|x^S_k − x*_k|² / ℓx² + |p^S_k/(m w*_k) − v*_k|² / ℓv²) + (β/2) |θ|²,

with `(x^S, p^S) = R_θ(x, p)` ([`remap_stages`](@ref)) and the targets `(xt, vt, wt)`, by the
Levenberg–Marquardt method from `θ = 0`: the damped Gauss–Newton step
`(A + μ diag(A)) δ = −g` with `A = JᵀJ + β I`, the damping `μ` divided by 3 after an accepted step
and multiplied by 4 after a rejected one. The iteration stops after `maxit` iterations, when no step
decreases `J` within 20 trials, or when the relative decrease falls below `1e-6`. Returns the
parameters, the objective at `θ = 0` and at the end, and the number of iterations.
"""
function fit_remap_parameters(
        par::AbstractRemapParametrisation, x, p, w, xt, vt, wt; m = 1,
        ℓx, ℓv, β = 1e-3, maxit = 40, trace = false)
    function residual(θ)
        X, P = remap_stages(par, θ, x, p, w; m)
        vcat((X .- xt) ./ ℓx, (P ./ (m .* wt) .- vt) ./ ℓv)
    end
    θ = zeros(remap_parameter_count(par))
    obj(θ, r) = 0.5 * dot(r, r) + 0.5β * dot(θ, θ)
    r = residual(θ)
    F = obj(θ, r)
    F0 = F
    μlm = 1e-2
    for it in 1:maxit
        Jm = _residual_jacobian(par, residual, θ, x, p, w, wt, m, ℓx, ℓv)
        g = Jm' * r + β * θ
        A = Symmetric(Jm' * Jm + β * I)
        accepted = false
        for _ in 1:20
            δ = -((A + μlm * Diagonal(diag(A) .+ 1e-12)) \ g)
            θn = θ + δ
            rn = residual(θn)
            Fn = obj(θn, rn)
            if Fn < F
                accepted = true
                rel = (F - Fn) / F
                trace &&
                    @printf("    it %2d: F/F0 = %.4f, |g| = %.2e, μ = %.1e\n", it, Fn / F0,
                        norm(g), μlm)
                θ, r, F = θn, rn, Fn
                μlm = max(μlm / 3, 1e-8)
                rel < 1e-6 && return (; θ, F0, F, iterations = it)
                break
            end
            μlm *= 4
        end
        accepted || return (; θ, F0, F, iterations = it)
    end
    (; θ, F0, F, iterations = maxit)
end
