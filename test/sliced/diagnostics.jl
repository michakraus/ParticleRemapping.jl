using ParticleRemapping
using Random: MersenneTwister
using Test

@testset "periodic wrap" begin
    @test periodic_wrap(0.25, 1.0) == 0.25
    @test periodic_wrap(1.25, 1.0) == -0.75
    @test periodic_wrap(-1.0, 1.0) == -1.0
    @test periodic_wrap(1.0, 1.0) == -1.0
    @test isapprox(periodic_wrap(7.0, π), 7.0 - 2π; rtol = 8 * eps())
end

@testset "quadrant counts" begin
    x = [0.5, -0.5, -0.5, 0.5, 0.0, -0.1]
    v = [0.5, 0.5, -0.5, -0.5, 0.0, 0.2]
    @test phase_quadrant_counts(x, v) == (2, 2, 1, 1)
end

@testset "cell statistics" begin
    # one marker in the centre of each cell of a 4 × 2 grid on [−2, 2) × [−1, 1]
    xc = [-1.5, -0.5, 0.5, 1.5]
    vc = [-0.5, 0.5]
    x = [xi for xi in xc for _ in vc]
    v = [vi for _ in xc for vi in vc]
    s = phase_cell_statistics(x, v; gx = 4, gv = 2, Lx = 2.0, Vmax = 1.0)
    @test s.empty == 0
    @test s.rstd == 0
    @test s.counts == ones(Int, 4, 2)
    # move one marker into another cell, wrap one periodically, drop one with |v| > Vmax
    x2, v2 = copy(x), copy(v)
    x2[1] = 0.6                   # into cell (3, 1)
    x2[3] = -0.5 + 4.0            # the same cell after wrapping
    v2[8] = 1.5                   # not counted
    s2 = phase_cell_statistics(x2, v2; gx = 4, gv = 2, Lx = 2.0, Vmax = 1.0)
    @test s2.counts == [0 1; 1 1; 2 1; 1 0]
    @test (s2.empty, s2.min_count, s2.max_count) == (2, 0, 2)
    mean_count = 7 / 8
    @test isapprox(s2.rstd, sqrt(sum((s2.counts .- mean_count) .^ 2) / 8) / mean_count;
        rtol = 8 * eps())
    h = cell_count_history([Markers(x, v, ones(8)), Markers(x2, v2, ones(8))], [0.0, 1.0];
        gx = 4, gv = 2, Lx = 2.0, Vmax = 1.0)
    @test h.empty == [0, 2]
    @test h.rstd == [0.0, s2.rstd]
end

@testset "smoothed density difference" begin
    rng = MersenneTwister(3)
    x, v, w = 2π .* rand(rng, 50) .- π, randn(rng, 50), rand(rng, 50)
    @test smoothed_density_difference(x, v, w, x, v, w) == 0
    # doubling the weights doubles the density: relative difference 1
    @test isapprox(smoothed_density_difference(x, v, w, x, v, 2w), 1; rtol = 8 * eps() * 50)
    # periodic: shifting the positions by 2π changes nothing beyond round-off
    @test smoothed_density_difference(x, v, w, x .+ 2π, v, w) < 8 * eps() * 50 * 2π
end
