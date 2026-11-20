# Figures and animations of the gauge rezoning and the canonical transport.
#
# Writes into results/remapping/:
# - the standard pictures of the two remaps on the examples of the overview page:
#   remapping-{gauge,transport}-remap.pdf (S1), remapping-{gauge,transport}-stages.pdf (S2) and
#   remapping-{gauge,transport}-stages.mp4 (A1);
# - the figures of the worked and numerical examples, remapping-{ex1,ex2,gauge,ex3,ex4}.pdf, and
#   their animations remapping-{ex1-fold,ex2-homotopy,gauge-homotopy,ex3-lift,ex4-matching}.mp4;
# - the pictures of the algorithms and of the canonical structure,
#   remapping-{kkt,gauge-orbits,sheet}.pdf.
#
# Run, from the package root:
#   julia --startup-file=no --project=scripts scripts/remapping/make_figures.jl

using ParticleRemapping
using CairoMakie

set_theme!(remap_theme())

const OUT = joinpath(@__DIR__, "..", "..", "results", "remapping")
mkpath(OUT)

function savefig(name, f)
    path = joinpath(OUT, "remapping-$name.pdf")
    save(path, f)
    println("wrote ", relpath(path))
end
function animation(name, writer)
    path = writer(joinpath(OUT, "remapping-$name.mp4"))
    println("wrote ", relpath(path))
end

# the standard example of the gauge rezoning: S1, S2, A1
gauge = gauge_rezoning_setup()
mk, grid = gauge.markers, gauge.grid
r = remap(gauge.method, mk)
savefig("gauge-remap",
    plot_remap(mk, r.markers;
        titles = ("(a) before", "(b) after, τ = $(round(r.info.τ; digits = 3))")))
savefig("gauge-stages",
    plot_stages(remap_states(gauge.method, mk, r); columns = 2, panel_size = (400, 280)))
animation("gauge-stages",
    path -> animate_stages(path, remap_frames(gauge.method, mk, r)))
println("    gauge rezoning: ", conservation_report(mk, r.markers; grid))

# the standard example of the canonical transport: S1, S2, A1
transport = quantile_matching_setup(; N = 10_000)
tm = transport.markers
rt = remap(transport.admissible, tm)
period = (0, transport.grid.L)
savefig("transport-remap",
    plot_remap(tm, rt.markers; titles = ("(a) before", "(b) after"), period))
savefig("transport-stages",
    plot_stages(remap_states(transport.admissible, tm, rt); period, panel_size = (
        300, 260)))
animation("transport-stages",
    path -> animate_stages(path, remap_frames(transport.admissible, tm, rt); period))
println(
    "    canonical transport: ", conservation_report(tm, rt.markers; grid = transport.grid),
    ", largest |R| = ", rt.info.max_residual,
    ", W/w in ", extrema(rt.markers.w ./ tm.w))

# the charge-constraint example and its fold
ex1 = charge_constraint_example()
savefig("ex1", plot_charge_example(ex1))
animation("ex1-fold", path -> animate_charge_fold(path, ex1))

# the eight-marker example
ex2 = eight_marker_homotopy()
savefig("ex2", plot_eight_marker_example(ex2))
animation("ex2-homotopy", path -> animate_eight_marker_homotopy(path, ex2))
savefig("kkt", plot_kkt_structure(ex2.result.info.problem, ex2.result.info.μ, ex2.result.info.α))

# the numerical example of the gauge rezoning, and its homotopy in finer steps, in the first two
# cells, for the animation
savefig("gauge", plot_gauge_experiment(mk, r))
fine = remap(gauge_rezoning_setup(; dτ0 = 0.01).method, mk)
println("    animation homotopy: ", count(h -> h.accepted, fine.info.history) - 1,
    " accepted steps up to τ = ", round(fine.info.τ; digits = 4))
sel = findall(<(2grid.Δx), mk.x)
zoom = [RemapStage(Markers(f.markers.x[sel], f.markers.v[sel], f.markers.w[sel]; m = mk.m),
            "cells 1 and 2, " * f.label)
        for f in remap_frames(gauge.method, mk, fine; frames_per_step = 2)]
vl = 1.05maximum(abs, fine.markers.v[sel])
animation("gauge-homotopy",
    path -> animate_stages(path, zoom; window = (0, 2grid.Δx, -vl, vl)))

# the cotangent lift: closed form and numerical example
lift = quantile_matching_setup()
savefig("ex3", plot_lift_example(lift, quantile_matching_experiment(lift)))
frames = lift_matching_frames(lift)
animation("ex3-lift", path -> animate_lift_matching(path, frames, lift))

# the five-marker matching
ex4 = quantile_matching_example()
savefig("ex4", plot_matching_example(ex4))
animation("ex4-matching", path -> animate_matching_example(path, ex4))

# the canonical structure
savefig("gauge-orbits", plot_gauge_orbits())
savefig("sheet", plot_velocity_momentum_sheet())
