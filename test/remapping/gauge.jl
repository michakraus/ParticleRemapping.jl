using ParticleRemapping
using ForwardDiff
using LinearAlgebra
using Random
using Test

Random.seed!(31)

# A small marker set on four nodes, with weights, velocities and signs of their own.
function small_markers(; N = 40, L = 3.0, m = 1.3)
    x = L .* rand(N)
    v = (0.3 .+ 2 .* rand(N)) .* rand((-1, 1), N)
    w = 0.2 .+ rand(N)
    Markers(x, v, w; m), PeriodicGrid(L, 4)
end

@testset "groups, targets and gap ratios" begin
    mk, grid = small_markers()
    groups = velocity_groups(grid, mk.x, mk.v)
    @test issorted(first.(groups))
    @test sort(reduce(vcat, last.(groups))) == 1:length(mk)
    for ((c, s), idx) in groups
        @test all(a -> cellof(grid, mk.x[a]) == c && sign(mk.v[a]) == s, idx)
    end
    μstar = gauge_targets(grid, mk.x, mk.v)
    V = mk.v .* exp.(μstar)
    for (_, idx) in groups
        u = sort(abs.(mk.v[idx]))
        if length(idx) ≥ 3
            # equispaced between the extreme speeds, which keep their values, in the old order
            @test sort(abs.(V[idx])) ≈ range(u[1], u[end]; length = length(idx)) rtol=8eps()
            @test sortperm(abs.(V[idx])) == sortperm(abs.(mk.v[idx]))
            @test speed_gap_ratio(V[idx]) ≈ 1 rtol=1e-12
        else
            @test all(iszero, μstar[idx])
        end
    end
    @test speed_gap_ratio([0.0, 1.0, 4.0]) ≈ 3 / 2
    @test speed_gap_ratio([-1.0, 0.5, 1.5]) ≈ 1         # speeds 0.5, 1, 1.5
    @test speed_gap_ratio([-1.0, 0.5, 2.0]) ≈ 4 / 3     # speeds 0.5, 1, 2
    @test mean_gap_ratio(grid, mk.x, V) ≈ 1 rtol=1e-12
    @test mean_gap_ratio(grid, mk.x, mk.v) > 1
end

@testset "the gauge step: exact properties for any μ" begin
    mk, grid = small_markers()
    μ = 0.4 .* randn(length(mk))
    W = mk.w .* exp.(-μ)
    after = from_canonical(copy(mk.x), momentum(mk), W; m = mk.m)
    @test after.v ≈ exp.(μ) .* mk.v rtol=8eps()
    @test momentum(after) ≈ momentum(mk) rtol=8eps()
    # the grid current is unchanged, the other moments change by c^(k)(μ)
    @test grid_moment(grid, after, 1) ≈ grid_moment(grid, mk, 1) rtol=1e-14
    prob = GaugeRezoningProblem(grid, mk; K = (0, 2, 3))
    c, _ = gauge_constraints(prob, μ)
    for (q, k) in enumerate((0, 2, 3))
        rows = ((q - 1) * grid.Ng + 1):(q * grid.Ng)
        @test c[rows] ≈ grid_moment(grid, after, k) .- grid_moment(grid, mk, k) rtol=1e-12
    end
    # the weak identity for a test function φ(x, v)
    φ(x, v) = sin(x) * v^2 + exp(-v^2) * x
    lhs = sum(after.w .* φ.(after.x, after.v)) - sum(mk.w .* φ.(mk.x, mk.v))
    rhs = sum(mk.w .* (exp.(-μ) .* φ.(mk.x, exp.(μ) .* mk.v) .- φ.(mk.x, mk.v)))
    @test lhs ≈ rhs rtol=1e-12
    @test_throws r"^ArgumentError: the first moment is kept" GaugeRezoningProblem(grid, mk;
        K = (0, 1))
end

@testset "derivatives of the constrained problem" begin
    mk, grid = small_markers()
    prob = GaugeRezoningProblem(grid, mk; K = (0, 2, 3), s = 0.7)
    μ = 0.3 .* randn(length(mk))
    α = randn(3grid.Ng)
    τ = 0.6
    # the Jacobian G against automatic differentiation of the constraints
    _, G = gauge_constraints(prob, μ)
    @test Matrix(G) ≈ ForwardDiff.jacobian(μ -> gauge_constraints(prob, μ)[1], μ) rtol=1e-13
    @test maximum(count(!iszero, G; dims = 1)) == 2 * 3
    # the Lagrangian J + αᵀc: its gradient is the stationarity, its Hessian is diagonal
    lag(μ) = gauge_objective(prob, μ, τ) + dot(α, gauge_constraints(prob, μ)[1])
    @test gauge_objective_gradient(prob, μ, τ) ≈
          ForwardDiff.gradient(μ -> gauge_objective(prob, μ, τ), μ) rtol=1e-13
    @test gauge_stationarity(prob, μ, α, τ) ≈ ForwardDiff.gradient(lag, μ) rtol=1e-12
    Hlag = ForwardDiff.hessian(lag, μ)
    @test gauge_hessian(prob, μ, α) ≈ diag(Hlag) rtol=1e-12
    @test norm(Hlag - Diagonal(diag(Hlag))) == 0
    # the scale s only scales the multipliers: the stationarity for (s, α) is that for (1, s²α)/s²
    prob1 = GaugeRezoningProblem(grid, mk; K = (0, 2, 3), s = 1.0)
    @test gauge_stationarity(prob, μ, α, τ) ≈
          gauge_stationarity(prob1, μ, 0.49 .* α, τ) ./ 0.49 rtol=1e-12
end

@testset "the Newton step: Schur complement and full KKT" begin
    mk, grid = small_markers()
    prob = GaugeRezoningProblem(grid, mk)
    μ = 0.1 .* randn(length(mk))
    α = 0.01 .* randn(2grid.Ng)
    c, G = gauge_constraints(prob, μ)
    H = gauge_hessian(prob, μ, α)
    @test all(>(0), H)
    gJ = gauge_objective_gradient(prob, μ, 0.5)
    Δμ, αn = gauge_newton_step(H, G, gJ, c; solver = :schur)
    Δμk, αk = gauge_newton_step(H, G, gJ, c; solver = :kkt)
    @test Δμ ≈ Δμk rtol=1e-10
    @test αn ≈ αk rtol=1e-10
    # the step solves both block rows of the KKT system
    @test H .* Δμ .+ G' * αn ≈ -gJ atol=1e-13
    @test G * Δμ ≈ -c atol=1e-13
    S = gauge_schur_matrix(H, G)
    @test S ≈ S' rtol=8eps()
    @test isposdef(Symmetric(Matrix(S)))
    # periodic block-tridiagonal: nodes i and j couple only when they share a cell
    Ng = grid.Ng
    for (q, r) in Iterators.product(1:2, 1:2), i in 1:Ng, j in 1:Ng
        d = min(mod(i - j, Ng), mod(j - i, Ng))
        d > 1 && @test S[(q - 1) * Ng + i, (r - 1) * Ng + j] == 0
    end
    @test size(gauge_kkt_matrix(H, G)) == (length(mk) + 2Ng, length(mk) + 2Ng)
    @test_throws r"^ArgumentError: solver must be :schur or :kkt" gauge_newton_step(
        H, G, gJ, c;
        solver = :lu)
    @test_throws r"^ArgumentError: solver must be :schur or :kkt" GaugeRezoning(grid;
        solver = :lu)
end

@testset "remap: homotopy, constraints and the minimiser" begin
    mk, grid = small_markers(; N = 120)
    for solver in (:schur, :kkt), K in ((0, 2), (0,), (0, 2, 3))

        method = GaugeRezoning(grid; K, solver, dτ0 = 0.25)
        r = remap(method, mk)
        after, info = r.markers, r.info
        @test info.τ > 0
        @test after.v != mk.v
        prob = info.problem
        c, G = gauge_constraints(prob, info.μ)
        @test maximum(abs.(c) ./ prob.scale) < method.tol
        # first-order conditions and the second-order sufficient condition
        @test maximum(abs, gauge_stationarity(prob, info.μ, info.α, info.τ)) < 1e-9
        @test info.min_curvature > 0
        @test momentum(after) ≈ momentum(mk) rtol=8eps()
        @test after.x == mk.x
        @test all(sign.(after.v) .== sign.(mk.v))
        rep = conservation_report(mk, after; grid)
        @test rep.current < 1e-13
        0 in K && @test rep.charge < method.tol
        2 in K && @test rep.energy < method.tol
        # the homotopy history starts at τ = 0 and ends at the τ reached
        @test info.history[1].τ == 0
        @test last(filter(h -> h.accepted, info.history)).τ == info.τ
    end
end
