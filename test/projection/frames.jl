# The data of the animations and of the pictures of the algorithm steps, on small runs.

using ParticleRemapping
using Random
using Test

Random.seed!(20261126)

prob = hole_filling_problem(; nL = 16, hole_a = 0.9, hole_b = 0.8)
e = hole_filling_experiment(; problem = prob, projection = HatProjection(12, 13),
    variants = (:weight_dependent, :non_separable), maxit = 2)

@testset "Frames" begin
    @testset "remap animation" begin
        for v in (:weight_dependent, :non_separable)
            r = e.results[v]
            frames = remap_animation_frames(e, v; n = 4, nreweight = 3)
            @test length(frames) == 1 + 2 * 4 + 2 * 3
            @test frames[1] == (prob.x0, prob.v0, "before")
            # the end of the stages is the map with the fitted θ
            @test frames[9][1] ≈ r.X rtol = 8eps()
            @test frames[9][2] ≈ r.P ./ (prob.m .* prob.w0) rtol = 8eps()
            # the reweightings end at v1 and vc
            @test frames[12][2] ≈ r.v1 rtol = 8eps()
            @test frames[end][2] ≈ r.vc rtol = 8eps()
        end
    end

    @testset "remap animation of a remap result" begin
        # without a charge-conserving reweighting: the frames end at the remapped markers
        r = remap(KickDriftRemap(; projection = HatProjection(12, 13), maxit = 2), prob.markers)
        frames = remap_animation_frames(prob, r; n = 4, nreweight = 3)
        @test length(frames) == 1 + 2 * 4 + 3
        @test frames[1] == (prob.x0, prob.v0, "before")
        @test [t for (_, _, t) in frames[[2, 6]]] == ["stage 1", "stage 2"]
        @test mod.(frames[9][1], prob.L) ≈ r.markers.x rtol = 8eps()
        @test frames[end][2] ≈ r.markers.v rtol = 8eps()
        @test frames[end][3] == "reweighting to the target weights"
        # with it: three more frames, which end at the charge-conserving markers
        rc = remap(
            KickDriftRemap(; projection = HatProjection(12, 13), maxit = 2,
                charge = ChargeConservingReweighting(prob.L, 16)),
            prob.markers)
        frames = remap_animation_frames(prob, rc; n = 4, nreweight = 3)
        @test length(frames) == 1 + 2 * 4 + 2 * 3
        @test frames[end][2] ≈ rc.markers.v rtol = 8eps()
        @test frames[end][3] == "charge-conserving reweighting"
    end

    @testset "stage states" begin
        s = remap_stage_states(e, :weight_dependent)
        @test [t for (_, _, t) in s] ==
              ["before", "kick 1", "drift 1", "kick 2", "drift 2", "reweighting",
            "charge-conserving"]
        @test s[2][1] == prob.x0                     # a kick keeps the positions
        @test s[5][1] ≈ e.results[:weight_dependent].X rtol = 8eps()
        s = remap_stage_states(e, :non_separable)
        @test [t for (_, _, t) in s][2:3] == ["stage 1", "stage 2"]
    end

    @testset "DLRA push" begin
        d = dlra_weights_experiment(; problem = prob, fit = false, tolerances = (1e-4,))
        frames = dlra_push_frames(d; nsub = 50)
        @test length(frames) == 1 + round(Int, d.T / (50 * d.dt))
        @test frames[end].t ≈ d.T
        @test frames[end].x ≈ d.xT rtol = 1e-12
        @test all(f.rank >= 1 for f in frames)
        @test size(frames[end].F) == (length(d.grid.x), length(d.grid.v))
    end
end
