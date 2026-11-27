# The full algorithm through the method types on a small lattice with a hole: the kick–drift
# remap, the non-separable remap, with and without the charge-conserving reweighting.

using ParticleRemapping
using Random
using Test

Random.seed!(20261125)

prob = hole_filling_problem(; nL = 16, hole_a = 0.9, hole_b = 0.8)
mk = prob.markers
hp = HatProjection(12, 13)

@testset "Remap methods" begin
    # the four methods with the default solver, Gauss–Newton
    methods = (KickDriftRemap(; nξ = 8, nv = 10, nw = 4, projection = hp, maxit = 5),
        KickDriftRemap(; variant = :control_B, nξ = 8, nv = 10, projection = hp, maxit = 5),
        NonSeparableRemap(; nξ = 6, nv = 7, nℓ = 3, projection = hp, maxit = 5),
        NonSeparableRemap(; nξ = 6, nv = 7, nℓ = 3, projection = hp, maxit = 5,
            weighted = false))
    for method in methods
        variant = method isa KickDriftRemap ? method.variant :
                  "weighted = $(method.weighted)"
        @testset "$(nameof(typeof(method))) $variant" begin
            r = remap(method, mk)
            i = r.info
            @test i.F < i.F0 && i.objective_ratio ≈ i.F / i.F0
            @test method.solver === :gauss_newton
            @test i.status !== nothing
            @test nameof(typeof(i.status)) === :SolverStatus
            @test i.iterations == i.status.iterations <= 5
            @test i.momentum_change < 1e-14
            @test i.weight_change < 1e-14
            @test total_momentum(r.markers) ≈ total_momentum(mk) atol = 1e-14 *
                                                                        sum(abs, prob.p0)
            @test total_weight(r.markers) ≈ total_weight(mk) rtol = 1e-14
            # the result: positions modulo L, the target weights, momenta P
            @test all(0 .<= r.markers.x .< prob.L)
            @test r.markers.w == i.targets.w
            @test momentum(r.markers) ≈ i.P rtol = 8eps()
            @test r.markers.x ≈ mod.(i.X, prob.L) rtol = 8eps()
            @test length(i.θ) == remap_parameter_count(i.parametrisation)
            @test (i.stage_iterations === nothing) == (method isa KickDriftRemap)
        end
    end

    @testset "charge-conserving reweighting" begin
        charge = ChargeConservingReweighting(prob.L, 8)
        r = remap(
            NonSeparableRemap(; nξ = 6, nv = 7, nℓ = 3, projection = hp, maxit = 3,
                charge), mk)
        i = r.info
        @test i.charge_converged
        @test i.charge_residual < 1e-12
        @test i.charge_residual_target > 1e-3       # must fail: the target weights alone
        ρ0 = spline_grid_charge(charge.basis, mk.x, mk.w)
        @test spline_grid_charge(charge.basis, r.markers.x, r.markers.w) ≈ ρ0 rtol = 1e-12
        @test r.markers.w == i.wc
        @test momentum(r.markers) ≈ i.P rtol = 8eps()
        @test i.max_log_change ≈ maximum(abs, log.(i.wc ./ i.targets.w)) rtol = 8eps()
    end

    @testset "projected target weights" begin
        # without a density the target weights come from the lumped projection
        tg = projection_targets(hp, prob.x0, prob.v0, prob.w0; prob.L, prob.V)
        r = remap(
            KickDriftRemap(; variant = :control_B, nξ = 8, nv = 10, projection = hp,
                maxit = 2), mk)
        @test r.info.targets.w ≈ tg.w rtol = 8eps()
    end

    @testset "the seed of the start" begin
        # with no iteration the fit returns its start: the seeded draw, or θ = 0 without a seed
        seeded(seed) = remap(
            KickDriftRemap(; nξ = 8, nv = 10, nw = 4, projection = hp, maxit = 0,
                solver = :gauss_newton, seed), mk).info.θ
        θ3 = seeded(3)
        @test θ3 == 1e-3 .* randn(Xoshiro(3), length(θ3))
        @test seeded(Int32(3)) == θ3
        @test seeded(nothing) == zeros(length(θ3))
    end

    @testset "nonpositive target weights" begin
        # hats finer than the marker spacing: f_h vanishes at some targets
        @test_throws r"^DomainError(?s:.)*the target weights must be positive" remap(
            KickDriftRemap(; projection = HatProjection(32, 33)), mk)
    end

    @testset "methods of the experiment" begin
        @test hole_remap_method(:control_A).variant === :control_A
        @test hole_remap_method(:non_separable).weighted
        @test !hole_remap_method(:control_C).weighted
        @test_throws r"^ArgumentError: unknown variant" hole_remap_method(:other)
        mw = hole_remap_method(:weight_dependent; solver = :gauss_newton, seed = 5)
        @test mw.solver === :gauss_newton && mw.seed === 5
        mn = hole_remap_method(:non_separable; solver = :gauss_newton)
        @test mn.solver === :gauss_newton && mn.seed === nothing
        @test KickDriftRemap(; seed = Int32(3)).seed === 3
        @test NonSeparableRemap(; seed = Int32(3)).seed === 3
        @test KickDriftRemap(; seed = 0).seed === 0 &&
              NonSeparableRemap(; seed = 0).seed === 0
        @test KickDriftRemap().solver === :gauss_newton
        @test NonSeparableRemap().solver === :gauss_newton
        for M in (KickDriftRemap, NonSeparableRemap)
            @test_throws r"^ArgumentError: unknown solver :levenberg_marquardt; use :gauss_newton$" M(;
                solver = :levenberg_marquardt)
            @test_throws r"^ArgumentError: unknown solver :newton; use :gauss_newton$" M(;
                solver = :newton)
            @test_throws r"^ArgumentError: the seed must be non-negative, got -1$" M(;
                seed = -1)
            @test_throws r"^ArgumentError: the seed must be at most typemax\(Int\), got 18446744073709551615$" M(;
                seed = typemax(UInt64))
            @test M(; seed = UInt64(typemax(Int))).seed === typemax(Int)
        end
        @test_throws r"^ArgumentError: the Hamiltonian particle projection acts on 1D-1V" remap(
            KickDriftRemap(), Markers(rand(2, 3), rand(2, 3), rand(3)))
    end

    # The remap examples of the manual, as the pages run them: each fit ends after its budget,
    # with the return code `MAXITERS` that the page states.
    @testset "the examples of the manual" begin
        maxiters(r, budget) = string(r.info.status.code) == "MAXITERS" &&
                              !r.info.status.stationary && r.info.iterations == budget
        # the home page
        home = remap(
            NonSeparableRemap(; nξ = 6, nv = 7, nℓ = 3, projection = hp, maxit = 5), mk)
        @test maxiters(home, 5)
        # the algorithm page and the overview of the projection
        Random.seed!(20261119)
        small = hole_filling_problem(; hole_x = π / 2, nL = 24, hole_a = 0.9, hole_b = 0.8)
        algorithm = remap(KickDriftRemap(; projection = hp, maxit = 10), small.markers)
        @test maxiters(algorithm, 10)
        overview = remap(
            KickDriftRemap(; projection = hp, maxit = 20,
                density = (xs, vs) -> small.density.(xs, vs)),
            small.markers)
        @test maxiters(overview, 20)
    end
end
