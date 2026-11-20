# Figures and animations of the Hamiltonian particle projection: the seven figures
# projection-{hole-pi2, hole-pi, charge, pm-v, mean-field, fold, dlra}.pdf, the four animations
# projection-{remap-non-separable-pi, remap-kick-drift-pi2, moser-targets, dlra}.mp4, and the
# pictures of the algorithm steps projection-{projection, target-flow, stages-kick-drift-pi2,
# stages-non-separable-pi}.pdf.
#
#   julia --startup-file=no --project=scripts scripts/projection/make_figures.jl [--results-dir=DIR]
#
# The output goes to DIR, by default results/projection/ in the package root. Runtime: about seven
# minutes (two hole-filling runs with three maps each, and the low-rank experiment).

using ParticleRemapping
using CairoMakie
using Random

function results_dir(args)
    i = findfirst(startswith("--results-dir="), args)
    i === nothing ? joinpath(@__DIR__, "..", "..", "results", "projection") :
    split(args[i], "="; limit = 2)[2]
end

function main(args = ARGS)
    dir = mkpath(results_dir(args))
    savefig(name, fig) = (path = joinpath(dir, "projection-$name.pdf"); save(path, fig);
        println("wrote ", path))
    anim(name) = joinpath(dir, "projection-$name.mp4")

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
    savefig("target-flow", plot_target_flow(e1.problem, e1.targets))
    savefig("stages-kick-drift-pi2",
        plot_remap_stages(remap_stage_states(e2, :weight_dependent), e2.problem.hole))
    savefig("stages-non-separable-pi",
        plot_remap_stages(remap_stage_states(e1, :non_separable), e1.problem.hole))

    # target weights after a push
    Random.seed!(20261119)
    d = dlra_weights_experiment()
    savefig("dlra", plot_dlra_weights(d))

    # animations
    for (e, v, name) in ((e1, :non_separable, "remap-non-separable-pi"),
        (e2, :weight_dependent, "remap-kick-drift-pi2"))
        animate_remap(anim(name), remap_animation_frames(e, v), e.problem.hole)
        println("wrote ", anim(name))
    end
    animate_moser_targets(anim("moser-targets"), moser_flow_history(e1.problem),
        e1.problem.hole)
    println("wrote ", anim("moser-targets"))
    animate_dlra(anim("dlra"), dlra_push_frames(d), d.grid)
    println("wrote ", anim("dlra"))
end

main()
