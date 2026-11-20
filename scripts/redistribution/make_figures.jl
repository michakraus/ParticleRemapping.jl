# Figures and animations of the symplectic particle redistribution.
#
#     julia --startup-file=no --project=scripts scripts/redistribution/make_figures.jl [--results-dir=DIR]
#
# Writes PDFs and MP4s into results/redistribution/ (or DIR).

using ParticleRemapping
using CairoMakie
using LinearAlgebra
using Printf

const RESULTS = let arg = findfirst(startswith("--results-dir="), ARGS)
    arg === nothing ? joinpath(@__DIR__, "..", "..", "results", "redistribution") :
    split(ARGS[arg], "="; limit = 2)[2]
end
mkpath(RESULTS)
out(name) = joinpath(RESULTS, name)

# the velocity–momentum fold of a Maxwellian
vth, m, ΔΩ = 1.3, 1.5, 0.2
f(v) = isotropic_maxwellian(v; vth = vth)
vs = collect(range(-4.5, 4.5; length = 601))
save(out("velocity_momentum_fold.pdf"),
    plot_velocity_momentum_fold(vs, velocity_fold_momentum.(f, vs; m = m, ΔΩ = ΔΩ),
        velocity_fold_determinant.(f, vs; m = m, ΔΩ = ΔΩ); vth = vth,
        pmax = maxwellian_fold_momentum_bound(; vth = vth, m = m, ΔΩ = ΔΩ)))

# the verification set-up and the steps of the redistribution
ex = redistribution_example()
mk = ex.markers
q, v, w = mk.x, mk.v, mk.w
p = momentum(mk)
ps = PhaseSpaceSplines()
res = remap(SplineRedistribution(), mk)
info = res.info
x̄, v̄ = info.targets_x, info.targets_v

save(out("spline_projections.pdf"),
    plot_spline_projections(phase_space_spline_grid(ps, l2_spline_projection(ps, q, v, w)),
        phase_space_spline_grid(ps, info.coefficients); markers = mk, targets = (x̄, v̄)))

cv(h) = kernel_density_cv(mod.(h[1], ex.Lx), h[2]; Lx = ex.Lx, σx = ex.σx, σv = ex.σv)
save(out("target_relaxation.pdf"),
    plot_target_relaxation(info.relaxation; Lx = ex.Lx, cv = cv.(info.relaxation)))

save(out("redistribution_markers.pdf"),
    plot_redistribution_markers(mk, res.markers; targets = (x̄, v̄), Lx = ex.Lx))

# the forward Newton solves: midpoint (one step), symplectic Euler (one step, which ends on a
# spurious root, and two substeps)
H = info.hamiltonian
p̄ = info.target_weights .* mk.m .* v̄
se1 = redistribution_substeps(H, SymplecticEulerStep(), q, p, x̄, p̄)
se2 = redistribution_substeps(H, SymplecticEulerStep(), q, p, x̄, p̄; substeps = 2)
save(out("newton_convergence.pdf"),
    plot_newton_convergence(only(info.history), only(se1.history), vcat(se2.history...);
        labels = [
            "midpoint, 1 step", "symplectic Euler, 1 step", "symplectic Euler, 2 substeps"]))

# the closed-loop defect
cl_ex = closed_loop_example()
cl = target_matching_closed_loop(cl_ex.H, cl_ex.z)
lhs, rhs = closed_loop_defect(cl.jacobians...)
vl = variational_closed_loop(cl_ex.z)
J4 = canonical_J(4)
save(out("closed_loop_defect.pdf"),
    plot_closed_loop_defect(lhs, rhs;
        variational = transpose(vl.jacobians[1]) * J4 * vl.jacobians[1] - J4))

# animations: the symplectic step in pseudo-time and the relaxation of the targets
path = pseudotime_path(H, ImplicitMidpointStep(), q, p, only(info.θ);
    s = range(0, 1; length = 49))
animate_redistribution(out("redistribution_step.mp4"), path, vcat(x̄, p̄); p_scale = info.p_scale)
animate_target_relaxation(out("target_relaxation.mp4"), info.relaxation; Lx = ex.Lx)

for file in sort(readdir(RESULTS))
    @printf("%-32s %8d bytes\n", file, filesize(joinpath(RESULTS, file)))
end
flush(stdout)
