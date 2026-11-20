using ParticleRemapping
using Random
using Test

Random.seed!(11)

@testset "construction" begin
    mk = Markers([0.1, 0.2], [1.0, -2.0], [0.5, 0.25]; m = 2)
    @test mk.m === 2.0
    @test length(mk) == 2
    @test spatial_dimension(mk) == 1
    @test_throws r"^DimensionMismatch: x has size" Markers([0.1, 0.2], [1.0], [1.0, 1.0])
    @test_throws r"^DimensionMismatch: 2 markers but 3 weights" Markers(
        [0.1, 0.2], [1.0, 2.0],
        [1.0, 1.0, 1.0])
    mk2 = Markers(rand(3, 5), rand(3, 5), rand(5); m = 1.5)
    @test spatial_dimension(mk2) == 3
end

@testset "momentum and velocity" begin
    x, v, w, m = rand(6), randn(6), 0.5 .+ rand(6), 1.7
    mk = Markers(x, v, w; m = m)
    @test momentum(mk) ≈ [m * w[a] * v[a] for a in 1:6] rtol=8eps()
    @test total_momentum(mk) ≈ sum(m * w[a] * v[a] for a in 1:6) rtol=8eps()
    @test kinetic_energy(mk) ≈ m / 2 * sum(w[a] * v[a]^2 for a in 1:6) rtol=8eps()
    @test total_weight(mk) ≈ sum(w)
    # d = 2: one column per marker, the weight multiplies every component
    X, V = rand(2, 4), randn(2, 4)
    W = [0.3, 0.7, 1.1, 2.0]
    mk2 = Markers(X, V, W; m = m)
    @test momentum(mk2) ≈ [m * W[a] * V[i, a] for i in 1:2, a in 1:4] rtol=8eps()
    @test total_momentum(mk2) ≈ [sum(m * W[a] * V[i, a] for a in 1:4) for i in 1:2] rtol=8eps()
    # from_canonical inverts canonical
    y, p = canonical(mk2)
    back = from_canonical(y, p, W; m = m)
    @test back.v ≈ V rtol=8eps()
    @test back.x == X
end

@testset "gauge step" begin
    x, v, w = rand(5), randn(5), 0.5 .+ rand(5)
    mk = Markers(x, v, w; m = 1.3)
    λ = 0.5 .+ rand(5)
    g = gauge_step(mk, λ)
    @test g.x == mk.x
    @test g.v ≈ λ .* v rtol=8eps()
    @test g.w ≈ w ./ λ rtol=8eps()
    # the canonical momenta, and so the total momentum, do not change
    @test momentum(g) ≈ momentum(mk) rtol=8eps()
    # the weight does change: the gauge step is not the identity on (x, v, w)
    @test abs(total_weight(g) - total_weight(mk)) > 1e-3
    @test gauge_step(mk, 2.0).v ≈ 2 .* v
    @test_throws r"^DomainError.*the gauge stretch must be positive"s gauge_step(mk, [
        1.0, -1.0, 1.0,
        1.0, 1.0])
end
