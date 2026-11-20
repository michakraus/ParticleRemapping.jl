# The dynamical low-rank solver for ∂_t f + v ∂_x f + E(x) ∂_v f = 0 on [0, 2π) × [−6, 6], d = 1:
# free streaming keeps the numerical rank 3 and converges with the step; the full grid matches the
# characteristics; the BUG solution converges to the full grid with the step (second order: the BUG
# integrator is not exact at full rank); the truncation error stays below n_steps ϑ; the
# evaluation off the grid. Cases that must fail: the start at rank 1, and the fixed rank 1.

using ParticleRemapping
using Printf
using Random
using Test

const L, Vd = 2π, 6.0
const nx, nv = 64, 256
const T, dt = 3.0, 0.01
f0(x, v) = (1 + 0.3 * cos(x)) * exp(-v^2 / 2) / sqrt(2π)
E(x) = -0.5 * sin(x)

l2(g, F) = sqrt(g.hx * g.hv * sum(abs2, F))

Random.seed!(20261122)

@testset "Dynamical low-rank solver" begin
    @testset "free streaming, rank 3" begin
        g = VlasovGrid(nx, nv, L, Vd, x -> 0.0)
        F0 = f0.(g.x, g.v')
        Ff = vlasov_full_solve(g, F0, T, dt)
        Fe = [(1 + 0.3 * cos(x - v * T)) * exp(-v^2 / 2) / sqrt(2π) for x in g.x, v in g.v]
        f, ranks = dlra_solve(g, F0, T, dt; tol = 1e-10, rmin = 4)
        f2, _ = dlra_solve(g, F0, T, dt / 2; tol = 1e-10, rmin = 4)
        Ff2 = vlasov_full_solve(g, F0, T, dt / 2)
        f1, ranks1 = dlra_solve(g, F0, T, dt; tol = 1e-10, rmin = 1)
        err = l2(g, lowrank_full(f) .- Ff) / l2(g, Ff)
        err2 = l2(g, lowrank_full(f2) .- Ff2) / l2(g, Ff2)
        err_full = l2(g, Ff .- Fe) / l2(g, Fe)
        err1 = l2(g, lowrank_full(f1) .- Ff) / l2(g, Ff)
        @printf("  free streaming: largest rank %d; error %.1e, %.1e at half the step; full grid against the closed form %.1e\n",
            maximum(ranks), err, err2, err_full)
        @printf("  free streaming from rank 1: largest rank %d, error %.2f\n",
            maximum(ranks1),
            err1)
        @test maximum(ranks) <= 3
        @test err < 1e-5 && err / err2 > 3        # second order in the step: 4 in the limit
        @test err_full < 1e-5
        # must fail: from rank 1 the projected dynamics stays at rank 1 (the note quotes 0.28)
        @test maximum(ranks1) == 1 && err1 > 0.25
    end

    g = VlasovGrid(nx, nv, L, Vd, E)
    F0 = f0.(g.x, g.v')
    Ff = vlasov_full_solve(g, F0, T, dt)
    function exact(g)
        x, v = vlasov_characteristics(E, repeat(g.x, length(g.v)),
            repeat(g.v; inner = length(g.x)), -T)
        reshape(f0.(x, v), length(g.x), length(g.v))
    end

    @testset "full grid against characteristics" begin
        err = l2(g, Ff .- exact(g)) / l2(g, Ff)
        gc = VlasovGrid(nx, nv ÷ 2, L, Vd, E)
        Fc = vlasov_full_solve(gc, f0.(gc.x, gc.v'), T, dt)
        errc = l2(gc, Fc .- exact(gc)) / l2(gc, Fc)
        @printf("  full grid against characteristics: %.1e (%.1e on half the velocity points)\n",
            err, errc)
        @test err < 1e-4                          # the note quotes 2.7e-5
        @test errc / err > 8                      # fourth order in v: 16 in the limit
    end

    @testset "BUG with negligible truncation" begin
        f, ranks = dlra_solve(g, F0, T, dt; tol = 1e-13, rmin = 4)
        f2, _ = dlra_solve(g, F0, T, dt / 2; tol = 1e-13, rmin = 4)
        Ff2 = vlasov_full_solve(g, F0, T, dt / 2)
        err = l2(g, lowrank_full(f) .- Ff) / l2(g, Ff)
        err2 = l2(g, lowrank_full(f2) .- Ff2) / l2(g, Ff2)
        @printf("  ϑ = 1e-13: largest rank %d; error %.1e, %.1e at half the step\n",
            maximum(ranks), err, err2)
        @test err < 1e-6 && err / err2 > 3        # second order in the step
    end

    @testset "BUG with truncation, ϑ = $ϑ" for ϑ in (1e-4, 1e-6)
        f, ranks = dlra_solve(g, F0, T, dt; tol = ϑ, rmin = 4)
        err = l2(g, lowrank_full(f) .- Ff) / l2(g, F0)
        mass0 = g.hx * g.hv * sum(F0)
        dm = abs(lowrank_mass(g, f) - mass0) / mass0
        @printf("  ϑ = %.0e: largest rank %d; error %.1e (bound %.1e); mass error %.1e\n",
            ϑ,
            maximum(ranks), err, length(ranks) * ϑ, dm)
        @test err < length(ranks) * ϑ
        @test maximum(ranks) < nx ÷ 2
        @test dm < 10ϑ                            # the note quotes 6.8e-5 and 3.2e-7
    end

    @testset "evaluation at points" begin
        f = lowrank_factorise(g, F0; tol = 1e-14)
        x, v = L .* rand(500), 8 .* rand(500) .- 4
        @test maximum(abs, lowrank_evaluate(g, f, x, v) .- f0.(x, v)) / maximum(F0) < 1e-4
    end

    @testset "numerical rank" begin
        σ = [3.0, 2.0, 1.0, 0.1, 0.01]
        @test numerical_rank(σ, 0.0) == 5
        @test numerical_rank(σ, 0.02) == 4
        @test numerical_rank(σ, sqrt(0.1^2 + 0.01^2)) == 3
        @test numerical_rank(σ, 10.0) == 1
    end

    @testset "must fail: fixed rank 1" begin
        f, _ = dlra_solve(g, F0, T, dt; tol = 0.0, rmin = 1, rmax = 1)
        err = l2(g, lowrank_full(f) .- Ff) / l2(g, Ff)
        @printf("  fixed rank 1: error %.2f\n", err)
        @test lowrank_rank(f) == 1
        @test err > 0.15                          # the note quotes 0.16
    end
end
