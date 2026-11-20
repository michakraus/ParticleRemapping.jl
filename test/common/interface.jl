using ParticleRemapping
using Random
using Test

Random.seed!(14)

@testset "conservation report" begin
    g = PeriodicGrid(2π, 8)
    x, v, w = 2π .* rand(100), randn(100), 0.5 .+ rand(100)
    before = Markers(x, v, w)
    after = gauge_step(before, 0.8 .+ 0.4 .* rand(100))
    r = conservation_report(before, after; grid = g)
    # the gauge step keeps the momenta and so the grid current
    @test abs(r.momentum) < 1e-13
    @test r.current < 1e-13
    # it changes the charge and the kinetic energy unless constrained
    @test r.charge > 1e-3
    @test r.energy > 1e-3
    @test r.weight ≈ total_weight(after) - total_weight(before)
    @test keys(conservation_report(before, after)) == (:weight, :momentum, :kinetic_energy)
end
