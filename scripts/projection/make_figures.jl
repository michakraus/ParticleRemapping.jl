# Figures and animations of the Hamiltonian particle projection, to results/projection/:
#
# - the standard pictures of the standard example of the chapter: projection-remap.pdf (S1),
#   projection-stages.pdf (S2) and projection-stages.mp4 (A1); the changes of the conserved
#   quantities (T1) are printed;
# - the structural statements: projection-{mean-field, fold, pm-v}.pdf;
# - the hole-filling experiment: projection-{hole-pi2, hole-pi, charge}.pdf, the steps of the
#   algorithm projection-{projection, target-flow, stages-kick-drift-pi2,
#   stages-non-separable-pi}.pdf and the animations
#   projection-{remap-non-separable-pi, remap-kick-drift-pi2, moser-targets}.mp4;
# - the target weights after a push: projection-dlra.pdf and projection-dlra.mp4.
#
#   julia --startup-file=no --project=scripts scripts/projection/make_figures.jl [--results-dir=DIR]
#
# The output goes to DIR, by default results/projection/ in the package root. Runtime: about ten
# minutes (two hole-filling runs with three maps each, and the low-rank experiment).

using ParticleRemapping
using CairoMakie
using Random

function results_dir(args)
    i = findfirst(startswith("--results-dir="), args)
    i === nothing ? joinpath(@__DIR__, "..", "..", "results", "projection") :
    split(args[i], "="; limit = 2)[2]
end

hole_decoration(hole) = (ax, _) -> plot_hole_outline!(ax, hole)

function main(args = ARGS)
    set_theme!(remap_theme())
    dir = mkpath(results_dir(args))
    savefig(name, fig) = (path = joinpath(dir, "projection-$name.pdf"); save(path, fig);
        println("wrote ", path))
    function animate(name, frames, hole)
        path = joinpath(dir, "projection-$name.mp4")
        animate_stages(path, frames; period = (0, hole.L), decorate! = hole_decoration(hole))
        println("wrote ", path)
    end

    # the standard example: the kick–drift remap of 24 × 24 markers with a hole at x_h = π/2
    Random.seed!(20261119)
    small = hole_filling_problem(; hole_x = π / 2, nL = 24, hole_a = 0.9, hole_b = 0.8)
    method = KickDriftRemap(; projection = HatProjection(12, 13), maxit = 20,
        density = (xs, vs) -> small.density.(xs, vs))
    result = remap(method, small.markers)
    tg = result.info.targets
    println("standard example: objective ratio ", result.info.objective_ratio, ", ",
        result.info.iterations, " iterations")
    println("T1: ",
        conservation_report(small.markers, result.markers; grid = PeriodicGrid(2π, 12)))
    savefig("remap",
        plot_remap(
            small.markers, result.markers; titles = ("(a) before", "(b) after the remap"),
            targets = (tg.x, tg.v), period = (0, 2π)))
    states = remap_states(method, small.markers, result)
    savefig("stages",
        plot_stages(states; columns = 3, period = (0, 2π), window = (0, 2π, -5, 5),
            decorate! = hole_decoration(small.hole), targets = (tg.x, tg.v),
            target_panels = [length(states)]))
    animate("stages", remap_frames(method, small.markers, result), small.hole)

    # the structural statements
    Random.seed!(20261119)
    savefig("mean-field", plot_mean_field_step(mean_field_step_example(mean_field_example())))
    savefig("fold", plot_implicit_weight_fold())
    savefig("pm-v", plot_pm_v_coupling(pm_v_coupling()))

    # hole filling at x_h = π/2 and π
    variants = (:weight_dependent, :non_separable, :control_C)
    Random.seed!(20261119)
    e2 = hole_filling_experiment(; hole_x = π / 2, variants)
    e1 = hole_filling_experiment(; hole_x = π, variants)
    savefig("hole-pi2", plot_hole_filling(e2))
    savefig("hole-pi", plot_hole_filling(e1))
    savefig("charge", plot_charge_correction(e2, e1; labels = ["x_h = π/2", "x_h = π"]))

    # the steps of the algorithm
    savefig("projection", plot_lumped_projection(e1.problem))
    savefig("target-flow",
        plot_target_flow(e1.problem, e1.targets; window = (π - 1.2, π + 1.2, 0, 2.4)))
    savefig("stages-kick-drift-pi2",
        plot_hole_stages(remap_states(experiment_remap(e2, :weight_dependent)...),
            e2.problem.hole))
    savefig("stages-non-separable-pi",
        plot_hole_stages(remap_states(experiment_remap(e1, :non_separable)...),
            e1.problem.hole; columns = 5))
    animate(
        "remap-non-separable-pi", remap_frames(experiment_remap(e1, :non_separable)...),
        e1.problem.hole)
    animate(
        "remap-kick-drift-pi2", remap_frames(experiment_remap(e2, :weight_dependent)...),
        e2.problem.hole)
    history = moser_flow_history(e1.problem)
    animate("moser-targets", [fill(history[1], 10); history; fill(history[end], 10)],
        e1.problem.hole)

    # target weights after a push
    Random.seed!(20261119)
    d = dlra_weights_experiment()
    savefig("dlra", plot_dlra_weights(d))
    path = joinpath(dir, "projection-dlra.mp4")
    animate_dlra(path, dlra_push_frames(d), d.grid)
    println("wrote ", path)
end

main()
