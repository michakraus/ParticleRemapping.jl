# The lumped projection onto hat functions and the density-interpolation targets.

using ParticleRemapping
using LinearAlgebra
using Random
using Test

Random.seed!(20261123)

@testset "Lumped projection and targets" begin
    hp = HatProjection(32, 33)
    @test hp.hs == 1 / 32 && hp.hr == 2 / 32

    @testset "partition of unity and normalisation" begin
        s, r = rand(500), 2 .* rand(500) .- 1
        q = rand(500)
        g = lumped_projection(hp, s, r, q)
        @test all(>=(0), g)
        # ∫ g_h = Σ q_k: the hat of node (i, j) has the integral hs hr (half at the clamped ends)
        vol = [hp.hs * hp.hr * ((j == 1 || j == hp.nV) ? 0.5 : 1.0)
               for i in 1:(hp.nX),
        j in 1:(hp.nV)]
        @test sum(g .* vol) ≈ sum(q) rtol = 8eps() * 500
        # interpolation reproduces the nodal values and a bilinear function
        f(s, r) = 0.3 + 0.2 * r
        G = [f((i - 1) * hp.hs, -1 + (j - 1) * hp.hr) for i in 1:(hp.nX), j in 1:(hp.nV)]
        @test hat_interpolate(hp, G, 0.2, 0.37) ≈ f(0.2, 0.37) rtol = 8eps()
        @test hat_interpolate(hp, G, 1.2, 0.37) ≈ f(0.2, 0.37) rtol = 8eps()   # periodic in s
    end

    @testset "Laplacian" begin
        A = moser_laplacian(hp)
        # constants span the null space: periodic in s, Neumann in r
        @test norm(A * ones(size(A, 1))) < 1e-9 * norm(A)
        @test rank(A) == size(A, 1) - 1
    end

    @testset "targets equidistribute the markers" begin
        prob = hole_filling_problem(; hole_x = π / 2)
        tg = projection_targets(hp, prob.x0, prob.v0, prob.w0; prob.L, prob.V,
            density = (xs, vs) -> prob.density.(xs, vs))
        @test sum(tg.w) ≈ sum(prob.w0) rtol = 8eps() * prob.N   # the total weight is exact
        @test all(abs.(tg.v) .<= prob.V)
        # the density deviation drops, and the hole of the markers is closed by the targets
        @test density_discrepancy(hp, tg.x, tg.v; prob.L, prob.V) <
              density_discrepancy(hp, prob.x0, prob.v0; prob.L, prob.V)
        @test hole_fill(prob.hole, tg.x, tg.v; prob.volume) > 0.8
        @test hole_fill(prob.hole, prob.x0, prob.v0; prob.volume) == 0
        # the target weights are c f at the targets with c ≈ |Ω| / N
        c = tg.w ./ prob.density.(tg.x, tg.v)
        @test maximum(c) ≈ minimum(c) rtol = 8eps()
        @test c[1] ≈ prob.volume / prob.N rtol = 0.02
        # the flow history starts at the markers and ends at the targets
        hist = moser_flow_history(prob; projection = hp)
        @test length(hist) == 41
        @test hist[1][1] ≈ prob.x0 rtol = 8eps()
        @test hist[1][2] ≈ prob.v0 rtol = 8eps()
        @test hist[end][1] ≈ tg.x rtol = 8eps()
        @test hist[end][2] ≈ tg.v rtol = 8eps()
        # without a density the weights come from the lumped projection f_h
        wp = projection_target_weights(hp, prob.x0, prob.v0, prob.w0, tg.x, tg.v; prob.L,
            prob.V)
        @test sum(wp) ≈ sum(prob.w0) rtol = 8eps() * prob.N
        @test maximum(abs, wp ./ tg.w .- 1) > 0.1   # the lumped projection has a large error
    end
end
