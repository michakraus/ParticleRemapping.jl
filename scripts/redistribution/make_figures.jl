# Figures and animations of the spline-projection redistribution.
#
#     julia --startup-file=no --project=scripts scripts/redistribution/make_figures.jl [--results-dir=DIR]
#
# Writes PDFs and MP4s, all named redistribution-…, into results/redistribution/ (or DIR).

using ParticleRemapping
using CairoMakie
using LinearAlgebra
using Printf
using Random

set_theme!(remap_theme())

const RESULTS = let arg = findfirst(startswith("--results-dir="), ARGS)
    arg === nothing ? joinpath(@__DIR__, "..", "..", "results", "redistribution") :
    split(ARGS[arg], "="; limit = 2)[2]
end
mkpath(RESULTS)
out(name) = joinpath(RESULTS, name)

# the verification set-up, the standard example of the chapter
ex = redistribution_example()
mk = ex.markers
method = SplineRedistribution()
res = remap(method, mk)
info = res.info
x̄, v̄ = info.targets_x, info.targets_v
period = (0.0, ex.Lx)

# S1, S2, A1: the standard pictures of the remap
save(out("redistribution-remap.pdf"),
    plot_remap(mk, res.markers; targets = (x̄, v̄), period = period,
        titles = ("(a) before: clustered markers", "(b) after: redistributed markers")))
save(out("redistribution-stages.pdf"),
    plot_stages(remap_states(method, mk, res); targets = (x̄, v̄), period = period))
animate_stages(out("redistribution-stages.mp4"), remap_frames(method, mk, res);
    targets = (x̄, v̄), period = period)

# the velocity–momentum fold of a Maxwellian
vth, m, ΔΩ = 1.3, 1.5, 0.2
f(v) = isotropic_maxwellian(v; vth = vth)
vs = collect(range(-4.5, 4.5; length = 601))
save(out("redistribution-fold.pdf"),
    plot_velocity_momentum_fold(vs, velocity_fold_momentum.(f, vs; m = m, ΔΩ = ΔΩ),
        velocity_fold_determinant.(f, vs; m = m, ΔΩ = ΔΩ); vth = vth,
        pmax = maxwellian_fold_momentum_bound(; vth = vth, m = m, ΔΩ = ΔΩ)))

# the L² and the KKT projections
q, v, w = mk.x, mk.v, mk.w
p = momentum(mk)
ps = PhaseSpaceSplines()
save(out("redistribution-projections.pdf"),
    plot_spline_projections(phase_space_spline_grid(ps, l2_spline_projection(ps, q, v, w)),
        phase_space_spline_grid(ps, info.coefficients); markers = mk, targets = (x̄, v̄)))

# the relaxation of the targets: the ten steps of the method, and forty longer steps for the
# animation
cv(h) = kernel_density_cv(mod.(h[1], ex.Lx), h[2]; Lx = ex.Lx, σx = ex.σx, σv = ex.σv)
save(out("redistribution-relaxation.pdf"),
    plot_target_relaxation(info.relaxation; Lx = ex.Lx, cv = cv.(info.relaxation)))
long = declustering_relaxation(q, v; Lx = ex.Lx, σx = ex.σx, σv = ex.σv, step = 0.02,
    iterations = 40)
animate_target_relaxation(out("redistribution-relaxation.mp4"), long.history; Lx = ex.Lx)

# outside the verification set-up (KNOWN_ISSUES K2): 900 markers, the default relaxation
rng = MersenneTwister(1)
x9, v9 = 2π .* rand(rng, 900), randn(rng, 900)
r9 = declustering_relaxation(x9, v9; Lx = method.Lx, σx = method.σx, σv = method.σv)
fig9 = plot_target_relaxation(r9.history; Lx = method.Lx)
hlines!(content(fig9[1, 1]), [-method.Vmax, method.Vmax]; color = :gray25, linestyle = :dash)
save(out("redistribution-relaxation-900.pdf"), fig9)

# the forward Newton solves: midpoint (one step), symplectic Euler (one step, which ends on a
# spurious root, and two substeps)
H = info.hamiltonian
p̄ = info.target_weights .* mk.m .* v̄
se1 = redistribution_substeps(H, SymplecticEulerStep(), q, p, x̄, p̄)
se2 = redistribution_substeps(H, SymplecticEulerStep(), q, p, x̄, p̄; substeps = 2)
save(out("redistribution-newton.pdf"),
    plot_newton_convergence(only(info.history), only(se1.history), vcat(se2.history...);
        labels = [
            "midpoint, 1 step", "symplectic Euler, 1 step", "symplectic Euler, 2 substeps"]))

# the closed-loop defect
cl_ex = closed_loop_example()
cl = target_matching_closed_loop(cl_ex.H, cl_ex.z)
lhs, rhs = closed_loop_defect(cl.jacobians...)
vl = variational_closed_loop(cl_ex.z)
J4 = canonical_J(4)
save(out("redistribution-closed-loop.pdf"),
    plot_closed_loop_defect(lhs, rhs;
        variational = transpose(vl.jacobians[1]) * J4 * vl.jacobians[1] - J4))

for file in sort(filter(startswith("redistribution-"), readdir(RESULTS)))
    @printf("%-36s %8d bytes\n", file, filesize(joinpath(RESULTS, file)))
end
flush(stdout)
