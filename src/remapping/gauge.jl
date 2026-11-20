# Construction 1: velocity rezoning by the gauge step (1D-1V).
#
# Each marker keeps x_a and p_a = m w_a v_a; its velocity is stretched, V_a = e^{μ_a} v_a, and its
# weight shrinks, W_a = w_a e^{-μ_a}. The log-stretches μ minimise a weighted distance to targets
# μ* subject to grid-moment constraints c_i^(k)(μ) = 0, k ∈ K, solved by Newton's method on the
# KKT system with a homotopy in the targets.

"""
    velocity_groups(grid::PeriodicGrid, x, v)

The groups of markers of the gauge rezoning: one vector of marker indices for each pair of a
cell `[x_c, x_{c+1})` and a sign of the velocity, as `(c, s) => indices` pairs sorted by `(c, s)`.
"""
function velocity_groups(grid::PeriodicGrid, x::AbstractVector, v::AbstractVector)
    groups = Dict{Tuple{Int, Int}, Vector{Int}}()
    for a in eachindex(x, v)
        push!(get!(groups, (cellof(grid, x[a]), Int(sign(v[a]))), Int[]), a)
    end
    sort!(collect(groups); by = first)
end

"""
    gauge_targets(grid::PeriodicGrid, x, v)

The target log-stretches `μ*` of the gauge rezoning. In each group of
[`velocity_groups`](@ref) with `n ≥ 3` markers the speeds are sorted and mapped onto the `n`
equispaced speeds between the smallest and the largest speed of the group, keeping their order;
`μ*_a = ln(|V*_a| / |v_a|)`. Markers in smaller groups get `μ*_a = 0`.
"""
function gauge_targets(grid::PeriodicGrid, x::AbstractVector, v::AbstractVector)
    μstar = zeros(float(eltype(v)), length(v))
    for (_, idx) in velocity_groups(grid, x, v)
        length(idx) < 3 && continue
        o = idx[sortperm(abs.(v[idx]))]
        lo, hi = abs(v[o[1]]), abs(v[o[end]])
        n = length(o)
        for (r, a) in enumerate(o)
            μstar[a] = log((lo + (r - 1) * (hi - lo) / (n - 1)) / abs(v[a]))
        end
    end
    μstar
end

"""
    speed_gap_ratio(u)

The ratio of the largest to the mean gap between the sorted values `|u|`: `1` for equispaced
speeds.
"""
function speed_gap_ratio(u::AbstractVector)
    d = diff(sort(abs.(u)))
    maximum(d) / (sum(d) / length(d))
end

"""
    mean_gap_ratio(grid::PeriodicGrid, x, v)

The measure of evenness of the gauge rezoning: the mean of [`speed_gap_ratio`](@ref) over the
groups of [`velocity_groups`](@ref) with at least three markers.
"""
function mean_gap_ratio(grid::PeriodicGrid, x::AbstractVector, v::AbstractVector)
    rs = [speed_gap_ratio(v[idx])
          for (_, idx) in velocity_groups(grid, x, v) if length(idx) ≥ 3]
    sum(rs) / length(rs)
end

"""
    GaugeRezoningProblem(grid, markers; K = (0, 2), s = 1, targets = gauge_targets(…))

The data of the constrained problem of the gauge rezoning of 1D-1V `markers` on a
[`PeriodicGrid`](@ref): the two nodes and the shape-function values of each marker, the
constrained moments `K`, the squared scale `s²`, the targets `μ*` and the scales
`M_i^(k)[w, |v|]` of the constraints. Constraint `(k, i)` is row `(q − 1) N_g + i` for the `q`-th
entry `k` of `K`.
"""
struct GaugeRezoningProblem{T <: Real, G <: PeriodicGrid, KT <: Tuple}
    grid::G
    x::Vector{T}
    v::Vector{T}
    w::Vector{T}
    m::T
    nodes::Matrix{Int}
    S::Matrix{T}
    K::KT
    s2::T
    μstar::Vector{T}
    scale::Vector{T}
end

function GaugeRezoningProblem(
        grid::PeriodicGrid, mk::Markers; K::Tuple = (0, 2), s::Real = 1,
        targets::AbstractVector = gauge_targets(grid, mk.x, mk.v))
    mk.x isa AbstractVector ||
        throw(ArgumentError("the gauge rezoning acts on 1D-1V markers"))
    1 in K &&
        throw(ArgumentError("the first moment is kept by every gauge step; drop k = 1 from K"))
    N = length(mk)
    nodes = Matrix{Int}(undef, N, 2)
    S = Matrix{eltype(mk.v)}(undef, N, 2)
    for a in 1:N
        i, j, Si, Sj = shape(grid, mk.x[a])
        nodes[a, 1], nodes[a, 2], S[a, 1], S[a, 2] = i, j, Si, Sj
    end
    scale = reduce(vcat, [grid_moment(grid, mk.x, abs.(mk.v), mk.w, k) for k in K])
    T = eltype(mk.v)
    GaugeRezoningProblem(
        grid, collect(mk.x), collect(mk.v), collect(mk.w), mk.m, nodes, S, K,
        T(s)^2, collect(T, targets), scale)
end

_nconstraints(prob::GaugeRezoningProblem) = prob.grid.Ng * length(prob.K)

# the multiplier α^(k) interpolated to marker a, A_a^(k) = Σ_i α_i^(k) S_i(x_a)
function _interpolate(prob::GaugeRezoningProblem, α, q, a)
    off = (q - 1) * prob.grid.Ng
    α[off + prob.nodes[a, 1]] * prob.S[a, 1] + α[off + prob.nodes[a, 2]] * prob.S[a, 2]
end

"""
    gauge_constraints(prob::GaugeRezoningProblem, μ) -> (c, G)

The constraints `c_i^(k)(μ) = Σ_a w_a S_i(x_a) v_a^k (e^{(k−1)μ_a} − 1)` for `k ∈ K` and their
sparse Jacobian `G_{(k,i),a} = (k − 1) w_a S_i(x_a) v_a^k e^{(k−1)μ_a}`, which has at most
`2|K|` non-zero entries per column.
"""
function gauge_constraints(prob::GaugeRezoningProblem{T}, μ::AbstractVector) where {T}
    Ng = prob.grid.Ng
    Tc = promote_type(T, eltype(μ))
    c = zeros(Tc, _nconstraints(prob))
    rows, cols, vals = Int[], Int[], Tc[]
    for (q, k) in enumerate(prob.K), a in eachindex(μ), j in 1:2
        row = (q - 1) * Ng + prob.nodes[a, j]
        base = prob.w[a] * prob.S[a, j] * prob.v[a]^k
        c[row] += base * (exp((k - 1) * μ[a]) - 1)
        push!(rows, row)
        push!(cols, a)
        push!(vals, base * (k - 1) * exp((k - 1) * μ[a]))
    end
    c, dropzeros!(sparse(rows, cols, vals, length(c), length(μ)))
end

"""
    gauge_objective(prob::GaugeRezoningProblem, μ, τ = 1)

The objective `J(μ) = Σ_a w_a (μ_a − τ μ*_a)² / (2 s²)` with the targets scaled by the homotopy
parameter `τ`.
"""
function gauge_objective(prob::GaugeRezoningProblem, μ::AbstractVector, τ::Real = 1)
    sum(prob.w .* (μ .- τ .* prob.μstar) .^ 2) / (2prob.s2)
end

"""
    gauge_objective_gradient(prob::GaugeRezoningProblem, μ, τ = 1)

The gradient `∇J(μ)_a = w_a (μ_a − τ μ*_a) / s²` of [`gauge_objective`](@ref).
"""
function gauge_objective_gradient(prob::GaugeRezoningProblem, μ::AbstractVector, τ::Real = 1)
    prob.w .* (μ .- τ .* prob.μstar) ./ prob.s2
end

"""
    gauge_stationarity(prob::GaugeRezoningProblem, μ, α, τ = 1)

The gradient of the Lagrangian `J(μ) + Σ α_i^(k) c_i^(k)(μ)` with respect to `μ`, which
decouples marker by marker:
`w_a [(μ_a − τ μ*_a)/s² + Σ_k (k − 1) A_a^(k) v_a^k e^{(k−1)μ_a}]` with the interpolated
multipliers `A_a^(k) = Σ_i α_i^(k) S_i(x_a)`.
"""
function gauge_stationarity(
        prob::GaugeRezoningProblem, μ::AbstractVector, α::AbstractVector,
        τ::Real = 1)
    g = gauge_objective_gradient(prob, μ, τ)
    for (q, k) in enumerate(prob.K), a in eachindex(μ)

        A = _interpolate(prob, α, q, a)
        g[a] += (k - 1) * prob.w[a] * A * prob.v[a]^k * exp((k - 1) * μ[a])
    end
    g
end

"""
    gauge_hessian(prob::GaugeRezoningProblem, μ, α)

The diagonal of the Hessian of the Lagrangian with respect to `μ`,
`H_a = w_a [1/s² + Σ_k (k − 1)² A_a^(k) v_a^k e^{(k−1)μ_a}]`; the off-diagonal entries vanish.
"""
function gauge_hessian(prob::GaugeRezoningProblem, μ::AbstractVector, α::AbstractVector)
    H = prob.w ./ prob.s2
    for (q, k) in enumerate(prob.K), a in eachindex(μ)

        A = _interpolate(prob, α, q, a)
        H[a] += A * prob.w[a] * prob.v[a]^k * (k - 1)^2 * exp((k - 1) * μ[a])
    end
    H
end

"""
    gauge_kkt_matrix(H, G)

The sparse KKT matrix `[diag(H) Gᵀ; G 0]` of the Newton step of the gauge rezoning.
"""
function gauge_kkt_matrix(H::AbstractVector, G::AbstractMatrix)
    nc = size(G, 1)
    Gs = sparse(G)
    [spdiagm(0 => H) Gs'; Gs spzeros(eltype(H), nc, nc)]
end

"""
    gauge_schur_matrix(H, G)

The Schur complement `G diag(H)⁻¹ Gᵀ` of the KKT matrix, of size `|K| N_g`: symmetric, periodic
block-tridiagonal with blocks of size `|K|`, and positive definite when every `H_a > 0` and `G`
has full row rank.
"""
gauge_schur_matrix(H::AbstractVector, G::AbstractMatrix) = G * Diagonal(1 ./ H) * G'

"""
    gauge_newton_step(H, G, ∇J, c; solver = :schur) -> (Δμ, α⁺)

The Newton step of the gauge rezoning, the solution of the KKT system
`[diag(H) Gᵀ; G 0] [Δμ; α⁺] = −[∇J; c]`. With `solver = :schur` the multipliers solve the
Schur-complement system `G diag(H)⁻¹ Gᵀ α⁺ = c − G diag(H)⁻¹ ∇J` and the stretches follow
marker by marker, `Δμ = −diag(H)⁻¹ (∇J + Gᵀ α⁺)`; with `solver = :kkt` the full sparse KKT
system is solved.
"""
function gauge_newton_step(H::AbstractVector, G::AbstractMatrix, gJ::AbstractVector,
        c::AbstractVector; solver::Symbol = :schur)
    if solver === :schur
        αn = Matrix(gauge_schur_matrix(H, G)) \ (c .- G * (gJ ./ H))
        Δμ = -(gJ .+ G' * αn) ./ H
        return Δμ, αn
    elseif solver === :kkt
        N = length(H)
        sol = gauge_kkt_matrix(H, G) \ vcat(-gJ, -c)
        return sol[1:N], sol[(N + 1):end]
    end
    throw(ArgumentError("solver must be :schur or :kkt, got :$solver"))
end

"""
    gauge_newton(prob, τ, μ, α; tol = 1e-13, steptol = 1e-10, maxit = 30, solver = :schur,
                 compare_solvers = false)

Newton's method on the optimality conditions of the gauge rezoning at the homotopy parameter
`τ` (step 2 of the algorithm *gauge rezoning*), from the stretches `μ` and multipliers `α`. Each
iteration evaluates `c`, `G`, `∇J` and `H`; rejects if some `H_a ≤ 0`; solves for the Newton
step ([`gauge_newton_step`](@ref)); stops when the scaled constraint residual is below `tol` and
`max |Δμ_a| < steptol`; and otherwise takes the first step length `t ∈ {1, 1/2, …}` that does
not increase the `ℓ¹` merit function ([`gauge_merit`](@ref)).

Returns `(; converged, μ, α, iterations)`, where `iterations` holds one `(; it, residual, step,
min_curvature, t, solver_difference)` per iteration (`step = NaN` on a rejection, `t = 0` on
convergence). With `compare_solvers = true` each step is solved both ways and
`solver_difference` is the largest difference between the Schur-complement and the full KKT
solution.
"""
function gauge_newton(prob::GaugeRezoningProblem{T}, τ::Real, μ0::AbstractVector,
        α0::AbstractVector; tol::Real = 1e-13, steptol::Real = 1e-10, maxit::Integer = 30,
        solver::Symbol = :schur, compare_solvers::Bool = false) where {T}
    μ, α = collect(T, μ0), collect(T, α0)
    ν = zero(T)
    iterations = _NewtonRecord[]
    record!(it, res, step, minHw, t, δ) = push!(iterations,
        (; it, residual = res, step, min_curvature = minHw, t, solver_difference = δ))
    for it in 0:(maxit - 1)
        c, G = gauge_constraints(prob, μ)
        gJ = gauge_objective_gradient(prob, μ, τ)
        H = gauge_hessian(prob, μ, α)
        res = maximum(abs.(c) ./ prob.scale)
        minHw = minimum(H ./ prob.w)
        # A non-positive H_a gives a direction of negative curvature in the kernel of G (one
        # marker against the many others of its nodes): the stationary point is not a minimum.
        if !(minHw > 0)
            record!(it, res, NaN, minHw, NaN, NaN)
            return (; converged = false, μ, α, iterations)
        end
        Δμ, αn = gauge_newton_step(H, G, gJ, c; solver)
        δsolve = NaN
        if compare_solvers
            Δμ2, αn2 = gauge_newton_step(
                H, G, gJ, c; solver = solver === :schur ? :kkt : :schur)
            δsolve = max(maximum(abs.(Δμ .- Δμ2)), maximum(abs.(αn .- αn2)))
        end
        step = maximum(abs.(Δμ))
        if res < tol && step < steptol
            record!(it, res, step, minHw, 0.0, δsolve)
            return (; converged = true, μ, α, iterations)
        end
        ν = max(ν, 2maximum(abs.(αn .* prob.scale)))
        m0 = gauge_merit(prob, μ, τ, ν)
        t = one(T)
        while gauge_merit(prob, μ .+ t .* Δμ, τ, ν) > m0 && t > 1e-10
            t /= 2
        end
        if t ≤ 1e-10
            res < 1e-10 || return (; converged = false, μ, α, iterations)
            t = one(T)
        end
        record!(it, res, step, minHw, t, δsolve)
        μ .+= t .* Δμ
        α = αn
    end
    (; converged = false, μ, α, iterations)
end

const _NewtonRecord = NamedTuple{
    (:it, :residual, :step, :min_curvature, :t,
        :solver_difference),
    Tuple{Int, Float64, Float64, Float64, Float64, Float64}}

"""
    gauge_merit(prob::GaugeRezoningProblem, μ, τ, ν)

The `ℓ¹` merit function `Ψ_ν(μ) = J(μ) + ν Σ_{i,k} |c_i^(k)(μ)| / M_i^(k)[w, |v|]` of the line
search, with each constraint scaled by the corresponding moment of `|v|`.
"""
function gauge_merit(prob::GaugeRezoningProblem, μ::AbstractVector, τ::Real, ν::Real)
    gauge_objective(prob, μ, τ) + ν * sum(abs.(gauge_constraints(prob, μ)[1]) ./ prob.scale)
end

"""
    gauge_homotopy(prob; dτ0 = 0.05, dτmin = 1e-3, newton...) -> (; μ, α, τ, history)

The homotopy of the algorithm *gauge rezoning*: from `μ = 0`, `α = 0`, `τ = 0`, try
`τ' = min(1, τ + Δτ)` with [`gauge_newton`](@ref) (keywords `newton...`) started at the last
accepted point. On success accept `τ'` and set `Δτ ← min(2Δτ, Δτ₀)`; on failure halve `Δτ`.
Stop at `τ = 1` or when `Δτ < Δτ_min`, a fold of the solution branch.

`history` lists every attempt as `(; τ, accepted, μ, α, iterations)`, starting with the point
`τ = 0`; `iterations` is the record of [`gauge_newton`](@ref).
"""
function gauge_homotopy(
        prob::GaugeRezoningProblem{T}; dτ0::Real = 0.05, dτmin::Real = 1e-3,
        newton...) where {T}
    N, nc = length(prob.w), _nconstraints(prob)
    μ, α = zeros(T, N), zeros(T, nc)
    τ, dτ = zero(T), T(dτ0)
    history = [(;
        τ, accepted = true, μ = copy(μ), α = copy(α), iterations = _NewtonRecord[])]
    while τ < 1 && dτ ≥ dτmin
        τn = min(one(T), τ + dτ)
        r = gauge_newton(prob, τn, μ, α; newton...)
        push!(history, (; τ = τn, accepted = r.converged, μ = r.μ, α = r.α,
            iterations = r.iterations))
        if r.converged
            τ, μ, α = τn, r.μ, r.α
            dτ = min(2dτ, T(dτ0))
        else
            dτ /= 2
        end
    end
    (; μ, α, τ, history)
end

"""
    GaugeRezoning(grid; K = (0, 2), s = 1, dτ0 = 0.05, dτmin = 1e-3, tol = 1e-13,
                  steptol = 1e-10, maxit = 30, solver = :schur)

The velocity rezoning by the gauge step on the [`PeriodicGrid`](@ref) `grid`: every marker keeps
its position and its momentum, its velocity is stretched by `λ_a = e^{μ_a}` and its weight divided
by `λ_a`. The stretches move the speeds of each group of cell and sign towards an equispaced set
([`gauge_targets`](@ref)), subject to the conservation of the grid moments `k ∈ K` (`0`: charge,
`2`: kinetic energy, `3`: heat flux). Newton's method ([`gauge_newton`](@ref), with tolerance
`tol` on the scaled constraints, `steptol` on the step, at most `maxit` iterations, and the
Schur-complement or full KKT `solver`) follows the minimiser by a homotopy in the targets
([`gauge_homotopy`](@ref)) with the initial step `dτ0` and the smallest step `dτmin`; the scale
`s` only scales the multipliers.
"""
struct GaugeRezoning{G <: PeriodicGrid, KT <: Tuple, T <: Real} <: AbstractRemap
    grid::G
    K::KT
    s::T
    dτ0::T
    dτmin::T
    tol::T
    steptol::T
    maxit::Int
    solver::Symbol

    function GaugeRezoning(grid::G, K::KT, s::T, dτ0::T, dτmin::T, tol::T, steptol::T,
            maxit::Int, solver::Symbol) where {G <: PeriodicGrid, KT <: Tuple, T <: Real}
        solver in (:schur, :kkt) ||
            throw(ArgumentError("solver must be :schur or :kkt, got :$solver"))
        new{G, KT, T}(grid, K, s, dτ0, dτmin, tol, steptol, maxit, solver)
    end
end

function GaugeRezoning(
        grid::PeriodicGrid; K::Tuple = (0, 2), s::Real = 1, dτ0::Real = 0.05,
        dτmin::Real = 1e-3, tol::Real = 1e-13, steptol::Real = 1e-10, maxit::Integer = 30,
        solver::Symbol = :schur)
    T = float(promote_type(
        typeof(s), typeof(dτ0), typeof(dτmin), typeof(tol), typeof(steptol)))
    GaugeRezoning(grid, K, T(s), T(dτ0), T(dτmin), T(tol), T(steptol), Int(maxit), solver)
end

"""
    GaugeRezoningProblem(method::GaugeRezoning, markers; targets = gauge_targets(…))

The problem of the gauge rezoning `method` for `markers`.
"""
function GaugeRezoningProblem(method::GaugeRezoning, mk::Markers;
        targets::AbstractVector = gauge_targets(method.grid, mk.x, mk.v))
    GaugeRezoningProblem(method.grid, mk; K = method.K, s = method.s, targets)
end

"""
    remap(method::GaugeRezoning, markers::Markers; targets = gauge_targets(…),
          compare_solvers = false) -> RemapResult

Run the algorithm *gauge rezoning*: the homotopy of [`gauge_homotopy`](@ref), then the gauge step
`W_a = w_a e^{−μ_a}`, `V_a = p_a / (m W_a)` with the last accepted `μ`, so that `P_a = p_a` holds
to round-off. The `info` holds the log-stretches `μ`, the targets `μstar`, the multipliers `α`,
the homotopy parameter `τ` reached, the smallest curvature `min_curvature = min_a H_a/w_a` there,
the `problem` and the `history` of the homotopy.
"""
function remap(method::GaugeRezoning, mk::Markers;
        targets::AbstractVector = gauge_targets(method.grid, mk.x, mk.v),
        compare_solvers::Bool = false)
    prob = GaugeRezoningProblem(method, mk; targets)
    h = gauge_homotopy(prob; dτ0 = method.dτ0, dτmin = method.dτmin, tol = method.tol,
        steptol = method.steptol, maxit = method.maxit, solver = method.solver,
        compare_solvers)
    W = mk.w .* exp.(-h.μ)
    after = from_canonical(copy(mk.x), momentum(mk), W; m = mk.m)
    minHw = minimum(gauge_hessian(prob, h.μ, h.α) ./ prob.w)
    RemapResult(after,
        (; μ = h.μ, μstar = prob.μstar, α = h.α, τ = h.τ, min_curvature = minHw,
            problem = prob, history = h.history))
end
