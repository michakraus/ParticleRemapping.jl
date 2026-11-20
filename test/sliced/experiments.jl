# Suites 1–5a of the verification of the sliced redistribution, with every control, and the probes
# of the constructions that do not work. The thresholds are those of the verification script; the
# reference numbers come from a cold run of the script (Julia 1.13.4, aarch64), which uses its own
# Dual numbers. The defects are round-off; the controls agree with the script to the round-off of
# the Jacobians (different pivoting of the solves, cond(∇²S) ≲ 10), hence rtol = 1e-13.

using ParticleRemapping
using Test

@testset "1. factorisation and pull-back on (x, v, w)" begin
    r = factorisation_pullback_experiment()
    @test r.err_full < 1e-13
    @test r.err_no_recoil > 1e-1
    @test isapprox(r.err_no_recoil, 0.5073615519131363; rtol = 1e-13)
    @test r.err_on_leaf < 1e-14
    @test r.err_across_dw > 1e-1
    @test isapprox(r.err_across_dw, 0.9740596015484777; rtol = 1e-13)
end

@testset "2. pseudo-time flow with and without recoil" begin
    r = pseudo_time_flow_experiment()
    @test r.err_recoil < 1e-13
    @test r.err_no_recoil > 1e-2
    @test isapprox(r.err_no_recoil, 0.030573267254231595; rtol = 1e-13)
    # with fewer fixed-point sweeps the step is not converged and not symplectic
    @test pseudo_time_flow_experiment(; iterations = 3).err_recoil > 1e-6
end

@testset "3. autonomous map in d = 1 and d = 2" begin
    r = autonomous_map_experiment()
    @test r.err_symp_1d < 1e-13
    @test r.err_ptot_1d < 1e-15
    @test r.err_nocot > 1e-1
    @test isapprox(r.err_nocot, 0.7600914624090467; rtol = 1e-13)
    @test r.err_badmu > 1e-2
    @test isapprox(r.err_badmu, 10.786509290207738; rtol = 1e-13)
    @test r.err_symp_2d < 1e-13
    @test r.err_pvec_2d < 1e-15
end

@testset "4. composition with a later flow step" begin
    r = composition_defect_experiment()
    @test r.defect > 1e-2
    @test isapprox(r.defect, 0.03572970047023824; rtol = 1e-13)
    @test r.err_formula < 1e-14
    @test r.err_const_w < 1e-15
    # the extended flow: θ̇ = −∂H/∂w̄ with the form Σ dq ∧ dp + Σ dw̄ ∧ dθ
    @test r.err_extended < 1e-14
    @test r.err_extended_wrong_sign > 1e-1
end

@testset "5a. hole filling on the checkerboard" begin
    r = checkerboard_experiment()
    @test r.quadrants_before == (200, 0, 200, 0)
    @test r.quadrants_marginal[2] + r.quadrants_marginal[4] < 20
    @test r.quadrants_marginal == (200, 0, 200, 0)
    @test all(c -> 65 <= c <= 135, r.quadrants_sliced)
    @test r.quadrants_sliced == (98, 101, 112, 89)
    # the lifts keep the total momentum; the weights do not change
    @test isapprox(sum(momentum(r.sliced)), sum(momentum(r.before));
        atol = 8 * eps() * 400 * maximum(abs, momentum(r.before)))
    @test r.sliced.w == r.before.w
end

@testset "what does not work" begin
    # Gaussian repulsion with a quartic wall: DΦ = I − α ∇²U loses definiteness, then the run
    # leaves the finite numbers
    r = repulsion_wall_probe()
    @test r.first_fold.λ < 0
    @test (r.first_fold.pass, r.first_fold.frame) == (3, 1)
    @test r.first_nonfinite == (4, 3)
    # with a weaker wall the map stays a diffeomorphism over the same passes
    @test repulsion_wall_probe(; k_w = 0.0).first_fold === nothing
    # one convex quantile lift of all markers in rotated frames: two markers reach Q2
    @test rotated_frame_lift_probe() ==
          (quadrants_before = (200, 0, 200, 0), quadrants_after = (198, 2, 200, 0))
    # reweighting the filled checkerboard at fixed momenta: small weights where the input
    # density vanishes stretch the velocities out of [−1, 1]
    c = checkerboard_reweighting_control()
    @test c.weight_ratio[1] < 0.2
    @test c.max_speed > 5
    @test maximum(abs, checkerboard_experiment().sliced.v) <= 1.05
end
