# The fit of the parameters of the remap map to the target points: GeometricSolvers' Gauss–Newton
# with a dogleg step, on the regularised least-squares objective.

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

# The solver of the fit; any other symbol is refused before any evaluation.
function _check_solver(solver::Symbol)
    solver === :gauss_newton ||
        throw(ArgumentError("unknown solver :$solver; use :gauss_newton"))
    solver
end

"""
    fit_remap_parameters(par, x, p, w, xt, vt, wt; m = 1, ℓx, ℓv, β = 1e-3, maxit = 40,
        trace = false, solver = :gauss_newton, θ0 = nothing)
        -> (; θ, F0, F, iterations, status)

The parameters `θ` of the map `par` that minimise

    J(θ) = ½ Σ_k (|x^S_k − x*_k|² / ℓx² + |p^S_k/(m w*_k) − v*_k|² / ℓv²) + (β/2) |θ|²,

with `(x^S, p^S) = R_θ(x, p)` ([`remap_stages`](@ref)) and the targets `(xt, vt, wt)`, from the
start `θ0` (`θ = 0` for `nothing`; a finite vector of [`remap_parameter_count`](@ref) entries,
converted to `Float64`), with at most `maxit ≥ 0` iterations, an integer.

The solver, `solver = :gauss_newton`, is GeometricSolvers' `GaussNewton(; globalization =
DogLeg())` on the least-squares problem of the residual `[r; √β θ]` with the Jacobian
`[∂r/∂θ; √β I]`, with `max_iterations = maxit` and GeometricSolvers' default tolerances; `status`
is its `SolverStatus`, whose `code` says why the solve ended (`SUCCESS` on a convergence test,
with `stationary` set for the stationarity test; `MAXITERS` after the budget). A trial `θ` with a
non-finite residual is left to the step rule of the solver. With `trace = true` it prints the
return code, the iterations and `F/F0` after the solve.

Returns the parameters, the objective `F0` at `θ = 0` (whatever the start) and `F` at the returned
`θ`, the number of iterations and the status. A `solver` other than `:gauss_newton`, a negative
`maxit`, a `maxit` above `typemax(Int32)` or a `β` that is negative, infinite or `NaN` throws an
`ArgumentError`, a non-integer `maxit` a `TypeError`, a `θ0` of the wrong length a
`DimensionMismatch`, and a non-finite `θ0` an `ArgumentError`.
"""
function fit_remap_parameters(
        par::AbstractRemapParametrisation, x, p, w, xt, vt, wt; m = 1,
        ℓx, ℓv, β = 1e-3, maxit::Integer = 40, trace = false,
        solver::Symbol = :gauss_newton,
        θ0::Union{Nothing, AbstractVector{<:Real}} = nothing)
    _check_solver(solver)
    maxit >= 0 || throw(ArgumentError("maxit must be non-negative, got $maxit"))
    maxit > typemax(Int32) &&
        throw(ArgumentError(
            "maxit must be at most $(typemax(Int32)), got $maxit"))
    isfinite(β) && β >= 0 ||
        throw(ArgumentError("β must be finite and non-negative, got $β"))
    n = remap_parameter_count(par)
    if θ0 !== nothing
        length(θ0) == n || throw(DimensionMismatch(
            "θ0 has length $(length(θ0)), the parametrisation has $n parameters"))
        all(isfinite, θ0) || throw(ArgumentError("θ0 must be finite"))
    end
    function residual(θ)
        X, P = remap_stages(par, θ, x, p, w; m)
        vcat((X .- xt) ./ ℓx, (P ./ (m .* wt) .- vt) ./ ℓv)
    end
    θ = zeros(n)
    obj(θ, r) = 0.5 * dot(r, r) + 0.5β * dot(θ, θ)
    r = residual(θ)
    F0 = obj(θ, r)
    if θ0 !== nothing
        θ = collect(Float64, θ0)
    end
    # GaussNewton with DogLeg on the augmented residual [r; √β θ], with `nr` rows of `r`
    nr, sβ = length(r), sqrt(β)
    function F!(y, θ, _)
        y[1:nr] .= residual(θ)
        y[(nr + 1):end] .= sβ .* θ
        y
    end
    function J!(Jy, θ, _)
        fill!(Jy, 0)
        Jy[1:nr, :] .= _residual_jacobian(par, residual, θ, x, p, w, wt, m, ℓx, ℓv)
        for i in 1:n
            Jy[nr + i, i] = sβ
        end
        Jy
    end
    prob = LeastSquaresProblem(F!, J!; output_length = nr + n)
    θ, status = solve(prob, θ, GaussNewton(; globalization = DogLeg()); max_iterations = maxit)
    F = obj(θ, residual(θ))
    trace && @printf("    %s after %d iterations: F/F0 = %.4f\n", status.code,
        status.iterations, F / F0)
    (; θ, F0, F, iterations = Int(status.iterations), status)
end
