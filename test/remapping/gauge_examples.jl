using ParticleRemapping
using Test

# A number quoted with d decimals in the note agrees to half a unit of its last digit.
agrees(x, quoted, d) = abs(x - quoted) ≤ 0.5 * 10.0^(-d) + 8eps(abs(quoted))

@testset "Lambert function" begin
    @test lambert_w0(0.0) ≈ 0 atol=8eps()
    @test lambert_w0(exp(1)) ≈ 1 rtol=8eps()
    @test lambert_w0(-exp(-1)) ≈ -1 atol=1e-7          # the branch point, a double root
    for z in (-0.3, -0.01, 0.2, 3.0, 50.0)
        w = lambert_w0(z)
        @test w * exp(w) ≈ z rtol=8eps()
        @test w ≥ -1
    end
    @test_throws r"^DomainError with -0\.4:\noutside the principal branch" lambert_w0(-0.4)
end

@testset "Example 1: one charge constraint in closed form" begin
    ex = charge_constraint_example()
    @test all(agrees.(ex.Vstar, [0.5, 1.3333, 2.1667, 3.0], 4))
    @test all(agrees.(ex.μstar, [0, 0.7985, 0.4367, 0], 4))
    s = charge_stationary_point(ex)
    # (1a) the Lambert-W formula solves the stationarity condition (μ − μ*)/s² = A e^{−μ}
    @test maximum(abs.((s.μ .- ex.μstar) ./ ex.s2 .- s.A .* exp.(-s.μ))) < 8eps()
    # (1a) the charge constraint holds
    @test sum(s.W) ≈ sum(ex.w) rtol=8eps()
    # control: the targets violate it, Σ w e^{−μ*} = 0.6385 instead of 1
    @test agrees(sum(ex.w .* exp.(-ex.μstar)), 0.6385, 4)
    # the curvature on the principal branch, H_a/w_a = (1 + W₀(z_a))/s²
    z = ex.s2 * s.A .* exp.(-ex.μstar)
    @test s.Hw ≈ (1 .+ lambert_w0.(z)) ./ ex.s2 rtol=1e-14
    # the numbers of the note
    @test agrees(s.A, -0.35537, 5)
    @test agrees(charge_multiplier_limit(ex), -0.36788, 5)
    @test charge_multiplier_limit(ex) ≈ -exp(-1) rtol=8eps()
    @test all(agrees.(s.μ, [-0.7595, 0.6043, 0.1222, -0.7595], 4))
    @test all(agrees.(s.λ, [0.4679, 1.8300, 1.1300, 0.4679], 4))
    @test all(agrees.(s.V, [0.2339, 1.0980, 1.5820, 1.4037], 4))
    @test all(agrees.(s.W, [0.2137, 0.2186, 0.3540, 0.2137], 4))
    @test all(agrees.(s.Hw, [0.2405, 0.8058, 0.6855, 0.2405], 4))
    @test agrees(sum(ex.w .* ex.v .^ 2), 1.8530, 4)
    @test agrees(sum(s.W .* s.V .^ 2), 1.5823, 4)
    @test s.V[4] < s.V[3]                              # the order of the speeds is lost
    # (1b) the fold: H_a/w_a ∝ (β_f − β)^{1/2} → 0, with β_f found to round-off by bisection
    βf = charge_fold(ex)
    @test agrees(βf, 1.3205, 4)
    sf = charge_stationary_point(ex, βf)
    @test minimum(sf.Hw) < 100sqrt(eps())
    @test sf.μ[[1, 4]] ≈ [-1, -1] atol=1e-7
    @test sf.A ≈ charge_multiplier_limit(ex, βf) atol=1e-12
    # (1b) no stationary point beyond the fold
    @test charge_stationary_point(ex, 1.001βf) === nothing
    # without a fold: heavier end markers (the research log)
    @test_throws r"^no fold below β = 10" charge_fold(charge_constraint_example(;
        w = [0.3, 0.2, 0.4, 0.1]))
end

@testset "Example 2: eight markers on two nodes" begin
    data = eight_marker_homotopy()
    ex, result = data.example, data.result
    mk, info = ex.markers, result.info
    prob = info.problem
    @test all(agrees.(
        prob.S[:, 1] .* (prob.nodes[:, 1] .== 1) .+
        prob.S[:, 2] .* (prob.nodes[:, 2] .== 1),
        [0.8, 0.6, 0.4, 0.2, 0.2, 0.4, 0.6, 0.8], 12))
    @test all(agrees.(mk.v .* exp.(info.μstar),
        [1.0, 1.5, 2.0, 2.5, 0.5, 0.9667, 1.4333, 1.9], 4))
    @test all(agrees.(info.μstar, [0, 0.3567, 0.5978, 0, 0, -0.6217, -0.2552, 0], 4))
    @test prob.scale[1:2] ≈ [1, 1] rtol=8eps()
    # control: the targets violate the constraints, c(μ*) = (0.0398, 0.0609, −0.0958, −0.1063)
    cstar, _ = gauge_constraints(prob, info.μstar)
    @test all(agrees.(cstar, [0.0398, 0.0609, -0.0958, -0.1063], 4))
    @test maximum(abs.(cstar) ./ prob.scale) > 1e-2

    # control: Newton at τ = 1 from μ = 0 without homotopy is rejected (min H/w = −0.065)
    direct = data.direct
    @test !direct.converged
    @test agrees(direct.iterations[1].step, 0.56, 2)
    @test isnan(direct.iterations[end].step)
    @test agrees(direct.iterations[end].min_curvature, -0.065, 3)

    # (2a) the Schur-complement step equals the full KKT step, at every step
    δ = [l.solver_difference for h in info.history for l in h.iterations]
    push!(δ, direct.iterations[1].solver_difference)
    @test maximum(δ) < 8eps()

    # (2b) the homotopy reaches τ = 1 in four accepted steps, without a rejection
    steps = info.history[2:end]
    @test [h.τ for h in steps] == [0.25, 0.5, 0.75, 1.0]
    @test all(h -> h.accepted, steps)
    @test all(h -> all(l -> l.t ∈ (0, 1), h.iterations), steps)    # full Newton steps
    last_step = steps[end].iterations
    @test length(last_step) == 5
    # the note rounds 1.95e-8 to 2.0e-8
    @test all(agrees.([l.residual for l in last_step[2:4]] ./ [1e-3, 1e-5, 1e-8],
        [3.50, 5.80, 1.95], 2))
    @test all(agrees.([l.step for l in last_step[1:4]] ./ [1e-1, 1e-2, 1e-4, 1e-7],
        [1.65, 1.81, 3.42, 1.19], 2))
    @test all(agrees.(
        [l.min_curvature for l in last_step], [
            0.740, 0.791, 0.824, 0.826, 0.826], 3))
    @test last_step[end].residual < ex.method.tol
    @test last_step[end].step < ex.method.steptol
    @test info.τ == 1

    # (2b) stationarity (to the last Newton step) and the constraints, and every H_a > 0
    c, _ = gauge_constraints(prob, info.μ)
    @test maximum(abs, gauge_stationarity(prob, info.μ, info.α, 1.0)) < ex.method.steptol
    @test maximum(abs.(c) ./ prob.scale) < ex.method.tol
    @test info.min_curvature > 0
    Hw = gauge_hessian(prob, info.μ, info.α) ./ prob.w
    @test all(agrees.(Hw, [0.908, 0.971, 1.024, 1.185, 1.104, 1.087, 0.951, 0.826], 3))

    # (2c) the grid current and every momentum are unchanged
    after = result.markers
    @test grid_moment(ex.grid, after, 1) ≈ grid_moment(ex.grid, mk, 1) rtol=8eps()
    @test momentum(after) ≈ momentum(mk) rtol=8eps()

    # the result of the note
    @test all(agrees.(info.α, [-0.1184, 0.1688, -0.0446, 0.0260], 4))
    @test all(agrees.(info.μ,
        [-0.0336, 0.3806, 0.6315, 0.0312, 0.0977, -0.5262, -0.2145, 0.0591], 4))
    @test all(agrees.(exp.(info.μ),
        [0.9669, 1.4632, 1.8805, 1.0317, 1.1026, 0.5908, 0.8070, 1.0609], 4))
    @test all(agrees.(after.v,
        [0.9669, 1.5363, 2.0685, 2.5792, 0.5513, 1.0635, 1.4929, 2.0157], 4))
    @test all(agrees.(after.w,
        [0.2585, 0.1709, 0.1329, 0.2423, 0.2267, 0.4231, 0.3098, 0.2356], 4))
    @test agrees(speed_gap_ratio(mk.v[1:4]), 2.80, 2)
    @test agrees(speed_gap_ratio(after.v[1:4]), 1.06, 2)
    @test agrees(speed_gap_ratio(mk.v[5:8]), 2.79, 2)
    @test agrees(speed_gap_ratio(after.v[5:8]), 1.07, 2)
    # the end markers, whose targets are μ* = 0, move by up to 0.1
    @test 0.05 < maximum(abs.(info.μ[[1, 4, 5, 8]])) < 0.1

    # the fine continuation of the figures ends at the same minimiser
    @test data.path.τ[end] == 1
    @test data.path.μ[end] ≈ info.μ atol=1e-12
    @test minimum(data.path.min_curvature) > 0
end
