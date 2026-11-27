# Figures and animations of the Hamiltonian particle projection, to results/projection/:
#
# - the standard pictures of the standard example of the chapter (`lattice_example()`):
#   projection-remap.pdf (S1), projection-stages.pdf (S2) and projection-stages.mp4 (A1); the
#   changes of the conserved quantities (T1) are printed;
# - the construction steps on `lattice_example()`: projection-{projection, target-flow}.pdf and
#   projection-moser-targets.mp4;
# - the structural statements: projection-{mean-field, fold, pm-v}.pdf;
# - the hole-filling experiment: projection-{hole-pi2, hole-pi, charge, stages-kick-drift-pi2,
#   stages-non-separable-pi}.pdf and the animations
#   projection-{remap-non-separable-pi, remap-kick-drift-pi2}.mp4;
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
    function animate(name, frames; period = (0, 2π), hole = nothing)
        path = joinpath(dir, "projection-$name.mp4")
        if hole === nothing
            animate_stages(path, frames; period)
        else
            animate_stages(
                path, frames; period = (0, hole.L), decorate! = hole_decoration(hole))
        end
        println("wrote ", path)
    end

    # the standard example: the kick–drift remap of the mildly deformed 24 × 24 lattice
    Random.seed!(20261119)
    ex = lattice_example()
    hp = HatProjection(12, 13)
    method = KickDriftRemap(;
        L = ex.L, V = ex.V, projection = hp, nξ = 16, nv = 20, nw = 8,
        maxit = 20, density = (xs, vs) -> ex.density.(xs, vs))
    result = remap(method, ex.markers)
    tg = result.info.targets
    println("standard example: objective ratio ", result.info.objective_ratio, ", ",
        result.info.iterations, " iterations, ", result.info.status.code)
    println("T1: ",
        conservation_report(ex.markers, result.markers; grid = PeriodicGrid(2π, 12)))
    savefig("remap",
        plot_remap(
            ex.markers, result.markers; titles = ("(a) before", "(b) after the remap"),
            targets = (tg.x, tg.v), period = (0, ex.L)))
    states = remap_states(method, ex.markers, result)
    savefig("stages",
        plot_stages(states; columns = 3, period = (0, ex.L), window = (0, ex.L, -6, 6),
            targets = (tg.x, tg.v), target_panels = [length(states)]))
    animate("stages", remap_frames(method, ex.markers, result); period = (0, ex.L))

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

    # the steps of the algorithm: construction on the deformed lattice, hole stages on e1/e2
    savefig("projection", plot_lumped_projection(ex; projection = hp))
    savefig("target-flow", plot_target_flow(ex, tg))
    savefig("stages-kick-drift-pi2",
        plot_hole_stages(remap_states(experiment_remap(e2, :weight_dependent)...),
            e2.problem.hole))
    savefig("stages-non-separable-pi",
        plot_hole_stages(remap_states(experiment_remap(e1, :non_separable)...),
            e1.problem.hole; columns = 5))
    animate(
        "remap-non-separable-pi", remap_frames(experiment_remap(e1, :non_separable)...);
        hole = e1.problem.hole)
    animate(
        "remap-kick-drift-pi2", remap_frames(experiment_remap(e2, :weight_dependent)...);
        hole = e2.problem.hole)
    history = moser_flow_history(ex; projection = hp)
    animate("moser-targets", [fill(history[1], 10); history; fill(history[end], 10)];
        period = (0, ex.L))

    # target weights after a push
    Random.seed!(20261119)
    d = dlra_weights_experiment()
    savefig("dlra", plot_dlra_weights(d))
    path = joinpath(dir, "projection-dlra.mp4")
    animate_dlra(path, dlra_push_frames(d), d.grid)
    println("wrote ", path)
end

main()
