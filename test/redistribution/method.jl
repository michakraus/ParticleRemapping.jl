# The whole redistribution through `remap`: projection, targets, parameter solve, symplectic
# step and weights, on the verification set-up.

using ParticleRemapping
using ForwardDiff
using LinearAlgebra
using Test

ex = redistribution_example()
mk = ex.markers

@testset "Set-up: lattice markers, leapfrog clustering" begin
    @test length(mk) == 16
    @test mk.m == 1.5
    @test all(0 .<= mk.x .< ex.Lx)
    @test ex.initial.x[1] ≈ π / 4 rtol = 8 * eps()
    @test ex.initial.v[1] ≈ -1.8 rtol = 8 * eps()
    @test ex.initial.w == mk.w                       # leapfrog keeps the weights
    @test norm(mk.x - ex.initial.x) > 1                # and moves the markers
end

@testset "remap: implicit midpoint in one step" begin
    res = remap(SplineRedistribution(), mk)
    after, info = res.markers, res.info
    @test info.target_error < 1e-13
    @test only(info.residual) < 1e-13
    @test abs(total_weight(after) - total_weight(mk)) < 1e-14
    @test abs(total_momentum(after) - total_momentum(mk)) < 1e-15
    @test all(>(0), after.w)
    @test after.w == info.target_weights
    @test after.v ≈ info.targets_v rtol = 1e-12
    @test mod.(info.targets_x, ex.Lx) ≈ after.x rtol = 1e-12
    @test info.cv_after < info.cv_before
    @test only(info.contraction) > 1
end

@testset "remap: symplectic Euler, one step and two substeps" begin
    one = remap(SplineRedistribution(; scheme = SymplecticEulerStep()), mk)
    @test one.info.target_error > 0.1                   # the caustic control
    two = remap(SplineRedistribution(; scheme = SymplecticEulerStep(), substeps = 2), mk)
    @test two.info.target_error < 1e-13
    @test all(<(1), two.info.contraction)
    @test abs(total_momentum(two.markers) - total_momentum(mk)) < 1e-15
end

@testset "remap: the whole map is canonical in (q, p) at fixed parameters" begin
    res = remap(SplineRedistribution(), mk)
    H = res.info.hamiltonian
    θ = only(res.info.θ)
    q, p = canonical(mk)
    step = forward_redistribution_step(H, ImplicitMidpointStep(), q, p, θ)
    J = canonical_J(length(q))
    @test norm(transpose(step.jacobian) * J * step.jacobian - J) < 1e-13
    # the end state in (x, v, w) has the momenta of the step
    @test momentum(res.markers) ≈ step.z[(length(q) + 1):end] rtol = 1e-12
end

@testset "remap rejects markers in more than one dimension" begin
    mk2 = Markers(rand(2, 4), rand(2, 4), ones(4))
    @test_throws r"^ArgumentError: SplineRedistribution redistributes 1D-1V markers" remap(
        SplineRedistribution(), mk2)
end
