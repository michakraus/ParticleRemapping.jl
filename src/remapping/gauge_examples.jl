# The worked examples of the gauge rezoning: one charge constraint in closed form, with the fold
# of its solution branch, and eight markers on two nodes.

"""
    lambert_w0(z)

The principal branch `W₀` of the Lambert function, the inverse of `u ↦ u eᵘ` on `[−1, ∞)`, for
`z ≥ −1/e`. A Newton iteration safeguarded by bisection converges also at the branch point
`z = −1/e`, where `W₀ = −1`.
"""
function lambert_w0(z::Real)
    z ≥ -exp(-one(z)) || throw(DomainError(z, "outside the principal branch, z < -1/e"))
    lo, hi = -one(float(z)), z < ℯ ? one(float(z)) : log(z)
    w = (lo + hi) / 2
    for _ in 1:200
        f = w * exp(w) - z
        f > 0 ? (hi = w) : (lo = w)
        wn = w - f / (exp(w) * (w + 1))
        wn = lo < wn < hi ? wn : (lo + hi) / 2
        abs(wn - w) ≤ 4eps(typeof(w)) * (1 + abs(w)) && return wn
        w = wn
    end
    w
end

"""
    charge_constraint_example(; v = [0.5, 0.6, 1.4, 3.0], w = [0.1, 0.4, 0.4, 0.1], s = 1)

The worked example of a single charge constraint: `n` markers on the support of one node with
`S = 1` for all of them, `K = {0}`. Returns `(; v, w, s2, Vstar, μstar)` with the equispaced
targets `V*` of the group and `μ* = ln(V*/v)`.
"""
function charge_constraint_example(; v::AbstractVector = [0.5, 0.6, 1.4, 3.0],
        w::AbstractVector = [0.1, 0.4, 0.4, 0.1], s::Real = 1.0)
    issorted(v) || throw(ArgumentError("the speeds v must be sorted"))
    n = length(v)
    Vstar = [v[1] + (r - 1) * (v[end] - v[1]) / (n - 1) for r in 1:n]
    (; v = collect(float.(v)), w = collect(float.(w)), s2 = float(s)^2, Vstar,
        μstar = log.(Vstar ./ v))
end

"""
    charge_stretches(ex, A, β = 1)

The stationary log-stretches of [`charge_constraint_example`](@ref) for the multiplier `A` and the
targets `β μ*`: `μ_a = β μ*_a + W₀(s² A e^{−β μ*_a})` ([`lambert_w0`](@ref)).
"""
function charge_stretches(ex, A::Real, β::Real = 1)
    β .* ex.μstar .+ [lambert_w0(ex.s2 * A * exp(-β * μs)) for μs in ex.μstar]
end

"""
    charge_constraint_residual(ex, A, β = 1)

The charge constraint `c⁽⁰⁾(A) = Σ_a w_a (e^{−μ_a(A)} − 1)` at the stationary stretches
[`charge_stretches`](@ref); it decreases strictly in `A`.
"""
function charge_constraint_residual(ex, A::Real, β::Real = 1)
    sum(ex.w .* (exp.(-charge_stretches(ex, A,
        β)) .- 1))
end

"""
    charge_multiplier_limit(ex, β = 1)

The smallest multiplier for which every marker of [`charge_constraint_example`](@ref) stays on the
principal branch of the Lambert function, `A_lim = −e^{β min_a μ*_a} / (s² e)`.
"""
charge_multiplier_limit(ex, β::Real = 1) = -exp(β * minimum(ex.μstar)) / (ex.s2 * exp(1))

"""
    charge_stationary_point(ex, β = 1)

The minimiser of [`charge_constraint_example`](@ref) for the targets `β μ*`, or `nothing` beyond
the fold. A minimiser exists if and only if `c⁽⁰⁾(A_lim) ≥ 0`; then the multiplier `A` is the
root of [`charge_constraint_residual`](@ref), found by bisection. Returns
`(; A, μ, λ, V, W, Hw)` with the curvatures `Hw = H_a/w_a = 1/s² + A e^{−μ_a}`.
"""
function charge_stationary_point(ex, β::Real = 1)
    lo, hi = charge_multiplier_limit(ex, β), 10.0
    charge_constraint_residual(ex, lo, β) > 0 || return nothing
    charge_constraint_residual(ex, hi, β) < 0 || error("no upper bracket of the multiplier")
    for _ in 1:200
        A = (lo + hi) / 2
        charge_constraint_residual(ex, A, β) > 0 ? (lo = A) : (hi = A)
        hi - lo < 1e-16 * (1 + abs(A)) && break
    end
    A = (lo + hi) / 2
    μ = charge_stretches(ex, A, β)
    (; A, μ, λ = exp.(μ), V = exp.(μ) .* ex.v,
        W = ex.w .* exp.(-μ), Hw = 1 / ex.s2 .+
                                   A .* exp.(-μ))
end

"""
    charge_fold(ex; βmax = 10)

The fold `β_f` of [`charge_constraint_example`](@ref): the largest target scale with a minimiser,
by bisection between `1` and `βmax`.
"""
function charge_fold(ex; βmax::Real = 10.0)
    lo, hi = 1.0, float(βmax)
    charge_stationary_point(ex, lo) === nothing && error("fold below β = 1")
    charge_stationary_point(ex, hi) === nothing || error("no fold below β = $βmax")
    for _ in 1:100
        β = (lo + hi) / 2
        charge_stationary_point(ex, β) === nothing ? (hi = β) : (lo = β)
    end
    lo
end

"""
    eight_marker_example(; m = 1)

The worked example of the gauge rezoning step by step: the periodic interval `[0, 2)` with two
nodes, eight markers of weight `0.25`, four in each cell, with positive velocities, `K = {0, 2}`,
`s = 1`. Returns `(; markers, grid, method)`, where `method` is the [`GaugeRezoning`](@ref) of
the example: the Schur-complement solver, `Δτ₀ = 1/4`, and the tolerances `10⁻¹⁴` on the scaled
constraints and `10⁻¹²` on the step.
"""
function eight_marker_example(; m::Real = 1.0)
    x = [0.2, 0.4, 0.6, 0.8, 1.2, 1.4, 1.6, 1.8]
    v = [1.00, 1.05, 1.10, 2.50, 0.50, 1.80, 1.85, 1.90]
    grid = PeriodicGrid(2.0, 2)
    method = GaugeRezoning(
        grid; K = (0, 2), s = 1.0, dτ0 = 0.25, dτmin = 1e-3, tol = 1e-14,
        steptol = 1e-12, solver = :schur)
    (; markers = Markers(x, v, fill(0.25, 8); m), grid, method)
end

"""
    eight_marker_homotopy(ex = eight_marker_example(); npath = 81)

Everything the figures of [`eight_marker_example`](@ref) show: the `result` of the
[`GaugeRezoning`](@ref) with its homotopy history (both solvers compared at every step), the
Newton attempt `direct` at `τ = 1` from `μ = 0` without homotopy, and a fine continuation `path`
over `npath` equispaced `τ`, each warm-started from the previous one, with the stretches `μ`, the
multipliers `α` and the smallest curvature `min_curvature` at each `τ`.
"""
function eight_marker_homotopy(ex = eight_marker_example(); npath::Integer = 81)
    method = ex.method
    result = remap(method, ex.markers; compare_solvers = true)
    prob = result.info.problem
    nc = length(result.info.α)
    newton = (; tol = method.tol, steptol = method.steptol, maxit = method.maxit,
        solver = method.solver)
    direct = gauge_newton(prob, 1.0, zeros(8), zeros(nc); newton..., compare_solvers = true)
    τs = range(0, 1; length = npath)
    μs, αs = [zeros(8)], [zeros(nc)]
    for τ in τs[2:end]
        r = gauge_newton(prob, τ, μs[end], αs[end]; newton...)
        r.converged || error("the continuation fails at τ = $τ")
        push!(μs, r.μ)
        push!(αs, r.α)
    end
    minH = [minimum(gauge_hessian(prob, μ, α) ./ prob.w) for (μ, α) in zip(μs, αs)]
    (; example = ex, result, direct,
        path = (; τ = collect(τs), μ = μs, α = αs,
            min_curvature = minH))
end
