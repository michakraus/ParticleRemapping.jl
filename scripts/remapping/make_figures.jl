# Figures and animations of the Hamiltonian remapping.
#
# Writes the figures remapping-{ex1,ex2,gauge,ex3,ex4}.pdf of the note, the animations
# remapping-{ex1-fold,ex2-homotopy,gauge-homotopy,ex3-lift,ex4-matching}.mp4, and the pictures
# of the algorithms remapping-{kkt,gauge-orbits,sheet}.pdf into results/remapping/.
#
# Run, from the package root:
#   julia --startup-file=no --project=scripts scripts/remapping/make_figures.jl

using ParticleRemapping
using CairoMakie

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

# the charge-constraint example and its fold
ex1 = charge_constraint_example()
savefig("ex1", plot_charge_example(ex1))
animation("ex1-fold", path -> animate_charge_fold(path, ex1))

# the eight-marker example
ex2 = eight_marker_homotopy()
savefig("ex2", plot_eight_marker_example(ex2))
animation("ex2-homotopy", path -> animate_eight_marker_homotopy(path, ex2))
savefig("kkt", plot_kkt_structure(ex2.result.info.problem, ex2.result.info.μ, ex2.result.info.α))

# the numerical example of the gauge rezoning, and its homotopy in finer steps for the animation
gauge = gauge_rezoning_setup()
savefig("gauge", plot_gauge_experiment(gauge.markers, remap(gauge.method, gauge.markers)))
fine = remap(gauge_rezoning_setup(; dτ0 = 0.01).method, gauge.markers)
println("    animation homotopy: ", count(h -> h.accepted, fine.info.history) - 1,
    " accepted steps up to τ = ", round(fine.info.τ; digits = 4))
animation("gauge-homotopy", path -> animate_gauge_homotopy(path, gauge.markers, fine))

# the cotangent lift: closed form and numerical example
transport = quantile_matching_setup()
savefig("ex3", plot_lift_example(transport, quantile_matching_experiment(transport)))
frames = lift_matching_frames(transport)
animation("ex3-lift", path -> animate_lift_matching(path, frames, transport))

# the five-marker matching
ex4 = quantile_matching_example()
savefig("ex4", plot_matching_example(ex4))
animation("ex4-matching", path -> animate_matching_example(path, ex4))

# the canonical structure
savefig("gauge-orbits", plot_gauge_orbits())
savefig("sheet", plot_velocity_momentum_sheet())
