# Suite 5b of the verification of the sliced redistribution: N = 900 markers in a two-wave field
# to T = 25, then one sliced remap; with the controls that skip the squeeze or shrink the
# mollifier. The thresholds are those of the verification script; the reference numbers come from
# a cold run of the script (Julia 1.13.4, aarch64), which the port reproduces to round-off.

using ParticleRemapping
using Test

@testset "5b. long-time run with one remap" begin
    r = long_time_remap_experiment()
    sb, sa = r.stats_before, r.stats_after
    @test sb.empty >= 10
    @test sa.empty == 0
    @test sa.rstd < 0.65 * sb.rstd
    @test r.l2_difference < 0.12
    @test r.momentum_error < 1e-15
    @test r.weight_error < 1e-15
    @test (sb.empty, sb.min_count, sb.max_count) == (14, 0, 20)
    @test (sa.empty, sa.min_count, sa.max_count) == (0, 1, 15)
    @test isapprox(sb.rstd, 0.5976086723909207; rtol = 8 * eps())
    @test isapprox(sa.rstd, 0.367197786605922; rtol = 8 * eps() * 900)
    @test isapprox(r.l2_difference, 0.10414395551641199; rtol = 8 * eps() * 900)
    # the Hessians of the two lifts are bounded below by 1 − max β
    @test r.min_hessian_q >= 1 - 0.78 - 8 * eps() * 900
    @test r.min_hessian_p >= 1 - 0.80 - 8 * eps() * 900
end

@testset "controls of the long-time run" begin
    base = long_time_remap_experiment()
    # weights from the density ratio set at fixed momenta, without the squeeze: velocities
    # stretch and cells empty again
    c = long_time_remap_experiment(; squeeze = false)
    @test c.stats_after.empty > 0
    @test c.stats_after.rstd > base.stats_before.rstd
    @test c.l2_difference > 3 * base.l2_difference
    @test c.momentum_error < 1e-15
    # β_mean = β_intra = 0.88: the slice-mean modes have the eigenvalue 1 − 0.88 of ∇²S
    # exactly, whatever ε, and amplify the dual variable by 1 / 0.12
    method = SlicedQuantileRemap(; β_mean_q = 0.88, β_intra_q = 0.88, β_mean_p = 0.88,
        β_intra_p = 0.88)
    s = long_time_remap_experiment(; method)
    @test isapprox(s.min_hessian_q, 0.12; atol = 8 * eps() * 900)
    @test isapprox(s.min_hessian_p, 0.12; atol = 8 * eps() * 900)
    # a mollifier of 0.1 target spacings instead of 2 moves the markers further from their
    # anchors in the dual lift
    excursion(ex) = maximum(abs, ex.info.q̂ .- ex.info.plan.q_anc)
    narrow = long_time_remap_experiment(;
        method = SlicedQuantileRemap(; ε_q_spacings = 0.1, ε_p_spacings = 0.1))
    @test excursion(narrow) > 1.4 * excursion(base)
end
