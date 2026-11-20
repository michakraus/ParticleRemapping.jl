# Verification of the sliced symplectic redistribution: suites 1–5 of the source verification
# script, with every control, and the probes of the constructions that do not work. Prints every
# number and checks it against the threshold of the source script.
#
#     julia --startup-file=no --project=scripts scripts/sliced/verify_sliced_redistribution.jl

using ParticleRemapping
using Printf: @printf
using Test

show_number(name, x) = (@printf("  %-34s %.16g\n", name, x); flush(stdout))
show_value(name, x) = (println("  ", rpad(name, 34), " ", x); flush(stdout))

@testset "sliced redistribution" begin
    @testset "1. factorisation and pull-back on (x, v, w)" begin
        println("1. factorisation and pull-back on (x, v, w)")
        r = factorisation_pullback_experiment()
        foreach(k -> show_number(string(k), r[k]), keys(r))
        @test r.err_full < 1e-13
        @test r.err_no_recoil > 1e-1
        @test r.err_on_leaf < 1e-14
        @test r.err_across_dw > 1e-1
    end

    @testset "2. pseudo-time flow" begin
        println("2. pseudo-time flow with and without recoil")
        r = pseudo_time_flow_experiment()
        foreach(k -> show_number(string(k), r[k]), keys(r))
        @test r.err_recoil < 1e-13
        @test r.err_no_recoil > 1e-2
        r3 = pseudo_time_flow_experiment(; iterations = 3)
        show_number("err_recoil, 3 sweeps", r3.err_recoil)
        @test r3.err_recoil > 1e-6
    end

    @testset "3. autonomous map" begin
        println("3. autonomous map in d = 1 and d = 2")
        r = autonomous_map_experiment()
        for k in (:err_symp_1d, :err_ptot_1d, :err_nocot, :err_badmu, :err_symp_2d, :err_pvec_2d)
            show_number(string(k), r[k])
        end
        @test r.err_symp_1d < 1e-13
        @test r.err_ptot_1d < 1e-15
        @test r.err_nocot > 1e-1
        @test r.err_badmu > 1e-2
        @test r.err_symp_2d < 1e-13
        @test r.err_pvec_2d < 1e-15
    end

    @testset "4. composition" begin
        println("4. composition with a later flow step")
        r = composition_defect_experiment()
        foreach(k -> show_number(string(k), r[k]), keys(r))
        @test r.defect > 1e-2
        @test r.err_formula < 1e-14
        @test r.err_const_w < 1e-15
        @test r.err_extended < 1e-14
        @test r.err_extended_wrong_sign > 1e-1
    end

    @testset "5a. checkerboard" begin
        println("5a. hole filling on the checkerboard (N = 400)")
        r = checkerboard_experiment()
        show_value("quadrants before", r.quadrants_before)
        show_value("quadrants sliced (B_v = 10, B_x = 20)", r.quadrants_sliced)
        show_value("quadrants marginal (one slice)", r.quadrants_marginal)
        show_number("max |v| after the sliced lifts", maximum(abs, r.sliced.v))
        show_number("max |x| after the sliced lifts", maximum(abs, r.sliced.x))
        show_number("max |x| after the marginal control", maximum(abs, r.marginal.x))
        show_number("max |v| after the marginal control", maximum(abs, r.marginal.v))
        @test r.quadrants_marginal[2] + r.quadrants_marginal[4] < 20
        @test all(c -> 65 <= c <= 135, r.quadrants_sliced)
    end

    @testset "5b. long-time run" begin
        println("5b. long-time run (N = 900, T = 25), one remap")
        t = @elapsed r = long_time_remap_experiment()
        for (label, s) in (("before", r.stats_before), ("after", r.stats_after))
            show_value("cells $label: empty, min, max", (s.empty, s.min_count, s.max_count))
            show_number("cells $label: relative std", s.rstd)
        end
        show_number("smoothed L2 density difference", r.l2_difference)
        show_number("total momentum error", r.momentum_error)
        show_number("total weight error", r.weight_error)
        show_number("smallest eigenvalue of ∇²S, Ψ_q", r.min_hessian_q)
        show_number("smallest eigenvalue of ∇²S, Ψ_p", r.min_hessian_p)
        show_number("run time (s)", t)
        @test r.stats_before.empty >= 10
        @test r.stats_after.empty == 0
        @test r.stats_after.rstd < 0.65 * r.stats_before.rstd
        @test r.l2_difference < 0.12
        @test r.momentum_error < 1e-15
        @test r.weight_error < 1e-15
    end

    @testset "what does not work" begin
        println("probes of constructions that do not work")
        r = repulsion_wall_probe()
        show_value("repulsion + wall: first fold", r.first_fold)
        show_value("repulsion + wall: first non-finite", r.first_nonfinite)
        show_value("repulsion, no wall: first fold", repulsion_wall_probe(; k_w = 0.0).first_fold)
        @test r.first_fold.λ < 0
        q = rotated_frame_lift_probe()
        show_value("one-slice lift in rotated frames", q)
        @test q.quadrants_after[2] + q.quadrants_after[4] <= 2
        c = checkerboard_reweighting_control()
        show_value("reweighting at fixed p: w̄/w range", c.weight_ratio)
        show_number("reweighting at fixed p: max |v̄|", c.max_speed)
        @test c.max_speed > 5
    end
end
