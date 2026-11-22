# The numbers of the standard experiments: the tables of the pages "Numerical experiments" of every
# chapter and of the comparison.
#
#     julia --startup-file=no --project=scripts scripts/standard/standard_experiments.jl
#
# Runs every construction on the checkerboard and on the long run, with both densities, and the
# runs with a remap after every 100 steps, about 15 minutes. The spline-projection redistribution
# also runs its reduced sizes.

using ParticleRemapping
using Printf

short(e) = first(replace(e, "\n" => " "), 100)

println("== the checkerboard (before: quadrants (200, 0, 200, 0), 32 of 64 cells empty)")
for variant in (:represented, :positive), key in STANDARD_METHODS

    r = standard_checkerboard(key; variant)
    if r.after === nothing
        println(rpad("$key, $variant", 28), "fails: ", short(r.error))
    else
        @printf("%-28s quadrants %s, weight in Q2 ∪ Q4 %.3f, empty %d, rstd %.3f, max shift %.2f, ΔW %.1e, ΔP %.1e\n",
            "$key, $variant", r.quadrants_after, r.hole_weight, r.stats_after.empty,
            r.stats_after.rstd, r.max_shift, r.report.weight, r.report.momentum)
    end
end

println("\n== the long run (before: 14 empty cells, rstd 0.598)")
for variant in (:represented, :positive), key in STANDARD_METHODS

    r = standard_long_time(key; variant)
    if r.after === nothing
        println(rpad("$key, $variant", 28), "fails: ", short(r.error))
    else
        s = r.stats_after
        @printf("%-28s empty %d, min %d, max %d, rstd %.3f, L2 %.3f, ΔW %.1e, ΔK %.3g, max |v| %.3g\n",
            "$key, $variant", s.empty, s.min_count, s.max_count, s.rstd, r.l2_difference,
            r.report.weight, r.report.kinetic_energy, r.max_v)
    end
end

println("\n== the spline-projection redistribution at reduced sizes")
for N in (16, 36, 64)
    r = standard_checkerboard(:spline; N)
    @printf("checkerboard N = %d: quadrants %s → %s, target error %.3g\n", N,
        r.quadrants_before, r.quadrants_after, r.result.info.target_error)
end
for n in (4, 6, 8)
    r = standard_long_time(:spline; Nx = n, Nv = n)
    @printf("long run %d × %d: empty %d → %d, target error %.3g\n", n, n,
        r.stats_before.empty, r.stats_after.empty, r.result.info.target_error)
end

println("\n== remaps after every 100 steps (mean rstd for t ≥ 5; empty cells at t = 0, 12.5, 25)")
for (key, variant) in ((:projection, :positive), (:nonseparable, :positive),
    (:gauge, :represented), (:transport, :represented), (:transport, :positive),
    (:sliced, :represented))
    h = standard_remap_history(key; variant)
    if h.error === nothing
        @printf("%-28s %.3f → %.3f, empty %s\n", "$key, $variant", h.rstd_plain,
            h.rstd_remap, h.h_remap.empty[[1, 26, 51]])
    else
        println(rpad("$key, $variant", 28), "fails: ", short(h.error))
    end
end
