using ParticleRemapping
using Random
using Test

Random.seed!(13)

g = PeriodicGrid(2π, 16)

@testset "cells and shape functions" begin
    @test g.Δx ≈ 2π / 16
    @test cellof(g, 0.0) == 1
    @test cellof(g, g.Δx * 3.5) == 4
    @test cellof(g, -0.1) == 16
    @test cellof(g, 2π) == 1
    i, j, Si, Sj = shape(g, g.Δx * 15.25)
    @test (i, j) == (16, 1)
    @test Si ≈ 0.75
    @test Sj ≈ 0.25
    x = 2π .* rand(50)
    S = shape_matrix(g, x)
    @test all(sum(S; dims = 2) .≈ 1)
    @test all(count(!iszero, S; dims = 2) .<= 2)
end

@testset "grid moments" begin
    x, v, w = 2π .* rand(200), randn(200), 0.5 .+ rand(200)
    for k in 0:3
        M = grid_moment(g, x, v, w, k)
        # partition of unity: the moments sum to Σ w v^k
        @test sum(M) ≈ sum(w .* v .^ k) rtol=1e-13
        @test M ≈ shape_matrix(g, x)' * (w .* v .^ k) rtol=1e-13
    end
    mk = Markers(x, v, w)
    @test grid_moment(g, mk, 1) == grid_moment(g, x, v, w, 1)
end
