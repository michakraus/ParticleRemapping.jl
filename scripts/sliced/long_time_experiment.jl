# The long-time experiment of the sliced redistribution: N = 900 markers in a two-wave electric
# field to T = 25, one remap at the end (suite 5b) and its variants, and a run with a remap after
# every 100 steps. Prints every number and checks the claims of the documentation.
#
#     julia --startup-file=no --project=scripts scripts/sliced/long_time_experiment.jl

using ParticleRemapping
using Printf: @printf
using Test

function row(label, ex)
    s = ex.stats_after
    @printf("%-28s %3d %3d %3d  %.4f  %.4f  %.4f  %.4f  %.3f  %.3f  %.3f\n",
        label, s.empty,
        s.min_count, s.max_count, s.rstd, ex.l2_difference, ex.min_hessian_q,
        ex.min_hessian_p, maximum(abs, ex.info.q̂ .- ex.info.plan.q_anc),
        maximum(abs, ex.info.q̄ .- ex.info.plan.q_anc), maximum(abs, ex.after.v))
    flush(stdout)
end

println("one remap at T = 25: variants")
println("variant                      emp min max  rstd    L2      λ_q     λ_p     " *
        "|q̂-a|  |q̄-a|  max|v̄|")
base = long_time_remap_experiment()
sb = base.stats_before
@printf("%-28s %3d %3d %3d  %.4f\n",
    "before the remap", sb.empty, sb.min_count, sb.max_count,
    sb.rstd)
row("sliced map (defaults)", base)
nosq = long_time_remap_experiment(; squeeze = false)
row("no squeeze (control)", nosq)
variants = [
    "β_mean = β_intra = 0.88" => SlicedQuantileRemap(; β_mean_q = 0.88, β_intra_q = 0.88,
        β_mean_p = 0.88, β_intra_p = 0.88),
    "ε = 0.1 spacings" => SlicedQuantileRemap(; ε_q_spacings = 0.1, ε_p_spacings = 0.1),
    "ε = 0.5 spacings" => SlicedQuantileRemap(; ε_q_spacings = 0.5, ε_p_spacings = 0.5),
    "ε = 1 spacing" => SlicedQuantileRemap(; ε_q_spacings = 1.0, ε_p_spacings = 1.0)]
results = Dict(label => long_time_remap_experiment(; method)
for (label, method) in variants)
for (label, _) in variants
    row(label, results[label])
end

println("\nruns to T = 25, cell statistics every 10 steps")
initial = electric_wave_markers()
t_plain = @elapsed plain = electric_wave_simulation(initial; record_every = 10)
t_remap = @elapsed remapped = electric_wave_simulation(initial; record_every = 10,
    remap_method = SlicedQuantileRemap(), remap_every = 100)
h_plain = cell_count_history(plain.snapshots, plain.times)
h_remap = cell_count_history(remapped.snapshots, remapped.times)
println("remap times: ", remapped.remap_times)
@printf("%6s  %12s  %12s\n", "t", "no remap", "remap/100")
for k in eachindex(h_plain.times)
    (k - 1) % 5 == 0 || continue
    @printf("%6.1f  %3d  %.4f  %3d  %.4f\n", h_plain.times[k], h_plain.empty[k],
        h_plain.rstd[k], h_remap.empty[k], h_remap.rstd[k])
end
@printf("run times: %.2f s without remaps, %.2f s with %d remaps\n", t_plain, t_remap,
    length(remapped.remap_times))
mean_after(h, from) = sum(h.rstd[from:end]) / length(h.rstd[from:end])
@printf("mean relative std for t ≥ 5: %.4f without remaps, %.4f with remaps\n",
    mean_after(h_plain, 11), mean_after(h_remap, 11))
@printf("total weight at T: %.16g with remaps, %.16g without\n", sum(remapped.markers.w),
    sum(plain.markers.w))

@testset "long-time claims" begin
    @test sb.empty == 14 && base.stats_after.empty == 0
    @test nosq.stats_after.empty > 0
    @test nosq.stats_after.rstd > sb.rstd
    @test isapprox(results["β_mean = β_intra = 0.88"].min_hessian_q, 0.12; atol = 1e-12)
    @test maximum(abs, results["ε = 0.1 spacings"].info.q̂ .- base.info.plan.q_anc) >
          maximum(abs, base.info.q̂ .- base.info.plan.q_anc)
    @test plain.markers.w == initial.w
    @test isapprox(sum(remapped.markers.w), sum(initial.w); rtol = 1e-13)
end
