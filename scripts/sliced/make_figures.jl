# The figures (PDF) and animations (MP4) of the sliced redistribution.
#
#     julia --startup-file=no --project=scripts scripts/sliced/make_figures.jl [--results-dir=DIR]
#
# Without `--results-dir` the files go to results/sliced/ in the package root.

using ParticleRemapping
using CairoMakie

function results_directory(args)
    for a in args
        startswith(a, "--results-dir=") && return a[(length("--results-dir=") + 1):end]
    end
    joinpath(@__DIR__, "..", "..", "results", "sliced")
end

function main(args)
    dir = mkpath(results_directory(args))
    written = String[]
    save_pdf(name, fig) = (
        path = joinpath(dir, name); save(path, fig); push!(written, path))

    stages = sliced_stages_example()
    save_pdf("sliced-stages.pdf", plot_sliced_stages(stages.before, stages.result))
    push!(written,
        animate_sliced_stages(joinpath(dir, "sliced-stages.mp4"), stages.before, stages.result))

    save_pdf("checkerboard.pdf", plot_checkerboard_filling(checkerboard_experiment()))
    save_pdf("wave-remap.pdf", plot_wave_remap(long_time_remap_experiment()))

    initial = electric_wave_markers()
    plain = electric_wave_simulation(initial; record_every = 10)
    remapped = electric_wave_simulation(initial; record_every = 10,
        remap_method = SlicedQuantileRemap(), remap_every = 100)
    h_plain = cell_count_history(plain.snapshots, plain.times)
    h_remap = cell_count_history(remapped.snapshots, remapped.times)
    save_pdf("cell-count-history.pdf",
        plot_cell_count_history([h_plain, h_remap];
            labels = ["no remap", "remap every 100 steps"],
            remap_times = remapped.remap_times))
    push!(written, animate_long_time_run(joinpath(dir, "long-time-run.mp4"), remapped, h_remap))
    push!(written,
        animate_long_time_run(joinpath(dir, "long-time-run-no-remap.mp4"), plain, h_plain))
    foreach(println, written)
end

main(ARGS)
