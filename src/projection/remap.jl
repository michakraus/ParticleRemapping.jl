# The full algorithm of the Hamiltonian particle projection for 1D-1V markers on [0, L) × [−V, V]:
# projection, target points, lift, fit, map, reweighting, and optionally the charge-conserving
# reweighting.

"""
    KickDriftRemap(; variant = :weight_dependent, L = 2π, V = 4, nξ = 32, nv = 40, nw = 16,
        stages = 2, projection = HatProjection(32, 33), moser_steps = 40, β = 1e-3, maxit = 40,
        density = nothing, charge = nothing, solver = :levenberg_marquardt, seed = nothing)

The Hamiltonian particle projection with the kick–drift map of [`KickDriftParametrisation`](@ref)
(`variant`, `nξ`, `nv`, `nw`, `stages`). The target points come from the density-interpolation
flow on the hats `projection` with `moser_steps` Runge–Kutta steps, the target weights from the
lumped projection `f_h` (`density = nothing`) or from `density(xs, vs)`; the fit is
[`fit_remap_parameters`](@ref) with the regularisation `β`, at most `maxit` iterations and the
`solver` (`:levenberg_marquardt` or `:gauss_newton`). It starts at `θ = 0` for `seed = nothing`,
and at `θ0 = 1e-3 .* randn(Xoshiro(seed), n)`, with `n` the number of parameters, for a
non-negative integer `seed`, stored as an `Int`. An unknown `solver`, or a `seed` that is
negative or above `typemax(Int)`, throws an `ArgumentError`. With `charge` a
[`ChargeConservingReweighting`](@ref), its weights replace the target weights.

    remap(method, markers::Markers; targets = nothing, reference_charge = nothing)

runs the algorithm on 1D-1V markers. Precomputed `targets` (a `NamedTuple` `(; x, v, w)` of
[`projection_targets`](@ref)) skip the first two steps. The charge-conserving reweighting keeps
`reference_charge`, by default the grid charge of `markers`. The result holds the markers
`(x′ mod L, p′/(m w′), w′)`; `info` holds the parameters `θ`, the objective `F0` at `θ = 0` and
`F` after the fit, `objective_ratio = F/F0`, the fit `iterations`, its `status` (the
`SolverStatus` of `:gauss_newton`, `nothing` for `:levenberg_marquardt`), the lifted positions
`X` and momenta `P` after the map, the `targets`, the parametrisation, the largest number of
fixed-point iterations of a stage (`stage_iterations`, `nothing` for explicit stages), the times
of the fit and of one map, the relative changes of the total momentum and of the total weight,
and, with `charge`, the charge-conserving weights `wc` and the diagnostics of
[`ChargeConservingReweighting`](@ref). Target weights that are not positive throw a
`DomainError`: the lumped projection vanishes at targets between markers when its hats are finer
than the marker spacing.
"""
struct KickDriftRemap{T <: Real, D, C, S <: Union{Nothing, Int}} <: AbstractRemap
    variant::Symbol
    L::T
    V::T
    nξ::Int
    nv::Int
    nw::Int
    stages::Int
    projection::HatProjection{T}
    moser_steps::Int
    β::T
    maxit::Int
    density::D
    charge::C
    solver::Symbol
    seed::S
end

# The seed of the start of the fit: `nothing`, or a non-negative integer stored as an `Int`.
_remap_seed(::Nothing) = nothing
function _remap_seed(seed::Integer)
    seed >= 0 || throw(ArgumentError("the seed must be non-negative, got $seed"))
    seed <= typemax(Int) ||
        throw(ArgumentError("the seed must be at most typemax(Int), got $seed"))
    Int(seed)
end

function KickDriftRemap(; variant::Symbol = :weight_dependent, L = 2π, V = 4.0, nξ = 32,
        nv = 40, nw = 16, stages = 2, projection = HatProjection(32, 33), moser_steps = 40,
        β = 1e-3, maxit = 40, density = nothing, charge = nothing,
        solver::Symbol = :levenberg_marquardt, seed::Union{Nothing, Integer} = nothing)
    KickDriftRemap(
        variant, float(L), float(V), nξ, nv, nw, stages, projection, moser_steps,
        float(β), maxit, density, charge, _check_solver(solver), _remap_seed(seed))
end

"""
    NonSeparableRemap(; L = 2π, V = 4, nξ = 16, nv = 17, nℓ = 6, stages = 2, weighted = true,
        iteration_limit = 50, projection = HatProjection(32, 33), moser_steps = 40, β = 1e-3,
        maxit = 40, density = nothing, charge = nothing, solver = :levenberg_marquardt,
        seed = nothing)

The Hamiltonian particle projection with the non-separable implicit-midpoint stages of
[`NonSeparableParametrisation`](@ref) (`nξ`, `nv`, `nℓ`, `stages`, `weighted`,
`iteration_limit`); otherwise as [`KickDriftRemap`](@ref), with the same `remap` method.
"""
struct NonSeparableRemap{T <: Real, D, C, S <: Union{Nothing, Int}} <: AbstractRemap
    L::T
    V::T
    nξ::Int
    nv::Int
    nℓ::Int
    stages::Int
    weighted::Bool
    iteration_limit::Int
    projection::HatProjection{T}
    moser_steps::Int
    β::T
    maxit::Int
    density::D
    charge::C
    solver::Symbol
    seed::S
end

function NonSeparableRemap(; L = 2π, V = 4.0, nξ = 16, nv = 17, nℓ = 6, stages = 2,
        weighted::Bool = true, iteration_limit = 50, projection = HatProjection(32, 33),
        moser_steps = 40, β = 1e-3, maxit = 40, density = nothing, charge = nothing,
        solver::Symbol = :levenberg_marquardt, seed::Union{Nothing, Integer} = nothing)
    NonSeparableRemap(float(L), float(V), nξ, nv, nℓ, stages, weighted, iteration_limit,
        projection, moser_steps, float(β), maxit, density, charge, _check_solver(solver),
        _remap_seed(seed))
end

"""
    remap_parametrisation(method, w)

The parametrisation of the map of `method` ([`KickDriftRemap`](@ref) or
[`NonSeparableRemap`](@ref)) for markers with the weights `w`.
"""
function remap_parametrisation(method::KickDriftRemap, w)
    KickDriftParametrisation(w; method.variant, method.L, method.V, method.nξ, method.nv,
        method.nw, method.stages)
end

function remap_parametrisation(method::NonSeparableRemap, w)
    NonSeparableParametrisation(w; method.L, method.V, method.nξ, method.nv, method.nℓ,
        method.stages, method.weighted, method.iteration_limit)
end

# The map with the fitted θ, and the largest number of fixed-point iterations of a stage.
function _apply_map(par::KickDriftParametrisation, θ, x, p, w, m)
    (kick_drift_map(par, θ, x, p, w; m)...,
        nothing)
end
function _apply_map(par::NonSeparableParametrisation, θ, x, p, w, m)
    r = nonseparable_map(par, θ, x, p, w; m)
    r.X, r.P, r.iterations
end

function remap(
        method::Union{KickDriftRemap, NonSeparableRemap}, mk::Markers; targets = nothing,
        reference_charge = nothing)
    mk.x isa AbstractVector ||
        throw(ArgumentError("the Hamiltonian particle projection acts on 1D-1V markers"))
    L, V, hp = method.L, method.V, method.projection
    x, v, w, m = mk.x, mk.v, mk.w, mk.m
    p = momentum(mk)
    tg = targets === nothing ?
         projection_targets(hp, x, v, w; L, V, steps = method.moser_steps,
        density = method.density) : targets
    all(>(0), tg.w) || throw(DomainError(minimum(tg.w),
        "the target weights must be positive; a lumped projection on hats finer than the marker spacing vanishes at some targets"))
    par = remap_parametrisation(method, w)
    θ0 = method.seed === nothing ? nothing :
         1e-3 .* randn(Xoshiro(method.seed), remap_parameter_count(par))
    fit_time = @elapsed fr = fit_remap_parameters(par, x, p, w, tg.x, tg.v, tg.w; m,
        ℓx = L / hp.nX, ℓv = 2V / (hp.nV - 1), method.β, method.maxit, method.solver, θ0)
    map_time = @elapsed X, P, stage_iterations = _apply_map(par, fr.θ, x, p, w, m)
    info = (; fr.θ, fr.F0, fr.F, objective_ratio = fr.F / fr.F0, fr.iterations, fr.status,
        X, P, targets = tg, parametrisation = par, stage_iterations, fit_time, map_time,
        momentum_change = abs(sum(P) - sum(p)) / sum(abs, p),
        weight_change = abs(sum(tg.w) - sum(w)) / sum(w))
    method.charge === nothing &&
        return RemapResult(Markers(mod.(X, L), P ./ (m .* tg.w), tg.w; m), info)
    b = method.charge.basis
    ρ = reference_charge === nothing ? spline_grid_charge(b, x, w) : reference_charge
    r = charge_reweight(b, X, tg.w, ρ; tol = method.charge.tol, maxit = method.charge.maxit)
    residual(u) = maximum(abs, spline_grid_charge(b, X, u) .- ρ) / maximum(ρ)
    info = (; info..., wc = r.w, λ = r.λ, charge_iterations = r.iterations,
        charge_converged = r.converged, charge_residual_target = residual(tg.w),
        charge_residual = residual(r.w), max_log_change = maximum(abs, log.(r.w ./ tg.w)))
    RemapResult(Markers(mod.(X, L), P ./ (m .* r.w), r.w; m), info)
end
