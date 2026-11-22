using ParticleRemapping
using Test

@testset "densities of the standard experiments" begin
    @test checkerboard_density(0.5, 0.5) == 0.5
    @test checkerboard_density(-0.5, 0.5) == 0.0
    @test checkerboard_density(-0.5, 0.5; variant = :positive) == 0.25
    @test checkerboard_density(1.5, 0.5; variant = :positive) == 0.0
    @test_throws ArgumentError checkerboard_density(0.0, 0.0; variant = :uniform)
    # the density of the long run at t = 0 is f₀ of the lattice, normalised like the weights
    mk = electric_wave_markers()
    f = [electric_wave_density(x, v; steps = 0) for (x, v) in zip(mk.x, mk.v)]
    @test isapprox(f .* (2π / 30) .* (2.8 / 30), mk.w; rtol = 1e-12)
    # the inverse of the discrete steps: f is constant along the computed characteristics
    T = standard_wave_markers(; steps = 40)
    f0 = [electric_wave_density(x, v; steps = 0) for (x, v) in zip(mk.x, mk.v)]
    fT = [electric_wave_density(x, v; steps = 40) for (x, v) in zip(T.x, T.v)]
    @test isapprox(fT, f0; rtol = 1e-9)
    # the represented density vanishes beyond the lattice, the positive one does not
    @test electric_wave_density(0.0, 3.0; steps = 0) == 0
    @test electric_wave_density(0.0, 3.0; steps = 0, variant = :positive) > 0
end

@testset "standard checkerboard" begin
    for key in (:gauge, :transport), variant in (:represented, :positive)

        r = standard_checkerboard(key; variant)
        @test r.error === nothing
        # positions and signs of v are kept: the empty quadrants stay empty
        @test r.quadrants_after == r.quadrants_before == (200, 0, 200, 0)
        @test r.max_shift == 0
        @test abs(r.report.momentum) < 1e-14
    end
    @test abs(standard_checkerboard(:gauge).report.weight) < 1e-14
    sliced = standard_checkerboard(:sliced)
    @test sliced.quadrants_after[2] + sliced.quadrants_after[4] > 100
    @test 0.1 < sliced.hole_weight < 0.3
    @test abs(sliced.report.weight) < 1e-14 && abs(sliced.report.momentum) < 1e-14
    # the projection needs a positive density on its box
    for key in (:projection, :nonseparable)
        r = standard_checkerboard(key)
        @test r.after === nothing
        @test occursin("target weights must be positive", r.error)
    end
    # the spline-projection redistribution: its relaxation leaves the box at full size
    @test occursin("Infs or NaNs", standard_checkerboard(:spline).error)
    @test_throws ArgumentError standard_checkerboard(:gauge; variant = :uniform)
    @test_throws ArgumentError standard_method(:unknown, :checkerboard)
end

@testset "standard long run" begin
    sliced = standard_long_time(:sliced)
    @test sliced.stats_before.empty == 14
    @test sliced.stats_after.empty == 0
    @test sliced.stats_after.rstd < 0.65 * sliced.stats_before.rstd
    @test sliced.l2_difference < 0.12
    gauge = standard_long_time(:gauge)
    @test gauge.after.x ≈ gauge.before.x
    @test abs(gauge.report.weight) < 1e-13 && abs(gauge.report.kinetic_energy) < 1e-13
    transport = standard_long_time(:transport; variant = :positive)
    @test transport.stats_after.empty == 0
    # the tails of the positive density are not sampled: the total weight falls
    @test -0.3 < transport.report.weight < -0.1
    @test occursin("target weights must be positive", standard_long_time(:projection).error)
    @test occursin("Infs or NaNs", standard_long_time(:spline).error)
end

@testset "standard remap history" begin
    h = standard_remap_history(:sliced)
    @test h.error === nothing
    @test h.remapped.remap_times ≈ [5, 10, 15, 20]
    # no cell is empty right after a remap; between remaps a few cells empty again
    after_remap = [findfirst(≈(t), h.h_remap.times) for t in h.remapped.remap_times]
    @test all(==(0), h.h_remap.empty[after_remap])
    @test maximum(h.h_remap.empty[h.h_remap.times .>= 5]) <= 3
    @test h.rstd_remap < 0.7 * h.rstd_plain
end
