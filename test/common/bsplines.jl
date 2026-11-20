using ParticleRemapping
using ForwardDiff
using Test

@testset "cubic B-spline" begin
    @test bspline3(0.0) ≈ 2 / 3
    @test bspline3(1.0) ≈ 1 / 6
    @test bspline3(2.0) == 0
    @test bspline3(-1.5) ≈ bspline3(1.5)
    # integer translates sum to one
    for t in range(0, 1; length = 7)
        @test sum(bspline3(t - j) for j in -3:3) ≈ 1
    end
    for t in (-1.7, -0.4, 0.3, 1.2)
        @test dbspline3(t) ≈ ForwardDiff.derivative(bspline3, t) rtol=1e-12
    end
end

@testset "uniform basis" begin
    b = UniformBasis(0.0, 2π / 12, 12, true)
    @test length(b) == 12
    for s in (0.0, 0.7, 3.1, 6.2, -0.5, 7.0)
        @test sum(basis_values(b, s)) ≈ 1
        n = Ref(0)
        eachnonzero((j, B, dB) -> (n[] += 1; @test(1 <= j <= 12)), b, s)
        @test n[] == 4
    end
    # derivatives agree with forward-mode differentiation of the values
    s = 1.3
    dB = zeros(12)
    eachnonzero((j, B, d) -> (dB[j] += d), b, s)
    @test dB ≈ ForwardDiff.derivative(t -> basis_values(b, t), s) rtol=1e-12
    # a non-periodic basis drops the splines beyond its ends
    bn = UniformBasis(-1.0, 0.5, 5, false)
    @test sum(basis_values(bn, 0.0)) ≈ 1
    @test sum(basis_values(bn, -1.0)) < 1
end
