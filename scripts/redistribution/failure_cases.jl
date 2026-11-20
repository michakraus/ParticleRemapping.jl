# Where the spline-projection redistribution works and where it fails: the numbers of the
# experiments page, of the limitations page and of KNOWN_ISSUES.md K2 and K3, and the basis sizes
# of the section on higher dimensions.
#
#     julia --startup-file=no --project=scripts scripts/redistribution/failure_cases.jl
#
# Runs about a minute; the hole-filling case (K3) takes most of it.

using ParticleRemapping
using LinearAlgebra
using Printf
using Random

println("Julia ", VERSION, ", check_bounds = ", Base.JLOptions().check_bounds)

# 1. The verification set-up: 16 markers, one implicit midpoint step (the standard example).
ex = redistribution_example()
mk = ex.markers
t = @elapsed res = remap(SplineRedistribution(), mk)
t = @elapsed res = remap(SplineRedistribution(), mk)
rep = conservation_report(mk, res.markers)
info = res.info
@printf("\n1. verification set-up, N = %d, one midpoint step (%.2f s warm)\n", length(mk),
    t)
@printf("   target error %.2e, contraction number %.3f, residual %.1e\n",
    info.target_error, only(info.contraction), only(info.residual))
@printf("   change of weight %.1e, momentum %.1e, kinetic energy %.4f (of %.4f)\n",
    rep.weight, rep.momentum, rep.kinetic_energy, kinetic_energy(mk))
@printf("   kernel-density CV %.4f -> %.4f, weights %.4f .. %.4f -> %.4f .. %.4f\n",
    info.cv_before, info.cv_after, extrema(mk.w)..., extrema(res.markers.w)...)
for (name, method) in (
    ("midpoint, 3 substeps", SplineRedistribution(; substeps = 3)),
    ("symplectic Euler, 1 step", SplineRedistribution(; scheme = SymplecticEulerStep())),
    ("symplectic Euler, 2 substeps",
    SplineRedistribution(; scheme = SymplecticEulerStep(), substeps = 2)))
    r = remap(method, mk)
    rr = conservation_report(mk, r.markers)
    @printf("   %-30s target error %.2e, contraction %s, momentum change %.1e\n", name,
        r.info.target_error, join((@sprintf("%.3f", c) for c in r.info.contraction), ", "),
        rr.momentum)
end
flush(stdout)

# 2. K2: 900 markers with the default relaxation step.
rng = MersenneTwister(1)
N2 = 900
x2 = 2π .* rand(rng, N2)
v2 = randn(rng, N2)
mk2 = Markers(x2, v2, fill(1 / N2, N2); m = 1.0)
d = SplineRedistribution()
rel = declustering_relaxation(x2, v2; Lx = d.Lx, σx = d.σx, σv = d.σv, step = d.relax_step,
    iterations = d.relax_iterations, v_wall = d.v_wall, k_wall = d.k_wall)
@printf("\n2. K2: %d markers, default relaxation (step %.3f, %d steps)\n", N2,
    d.relax_step, d.relax_iterations)
@printf("   max |v| before %.2f, max |v̄| after %.2f (Vmax = %.1f), largest displacement %.2f\n",
    maximum(abs, v2), maximum(abs, rel.v), d.Vmax,
    maximum(hypot.(rel.x .- x2, rel.v .- v2)))
@printf("   targets outside the velocity box: %d\n", count(>(d.Vmax), abs.(rel.v)))
message = try
    remap(d, mk2)
    "no error"
catch e
    first(split(sprint(showerror, e), '\n'))
end
println("   remap: ", message)
flush(stdout)

# 3. K3: the 244 markers of the hole-filling problem, relaxation step scaled by 16/256.
prob = hole_filling_problem(; nL = 16, hole_a = 0.9, hole_b = 0.8)
m3 = SplineRedistribution(; Vmax = 4.0, relax_step = 0.006 * 16 / 256)
t3 = @elapsed r3 = remap(m3, prob.markers)
rep3 = conservation_report(prob.markers, r3.markers)
@printf("\n3. K3: hole-filling problem, N = %d, relax_step = %.6f (%.0f s)\n", prob.N,
    m3.relax_step, t3)
@printf("   target error %.3g, contraction number %.3g, Newton residual %.3g\n",
    r3.info.target_error, only(r3.info.contraction), only(r3.info.residual))
@printf("   change of weight %.1e, momentum %.3g, kinetic energy %.4g (of %.4g)\n",
    rep3.weight, rep3.momentum, rep3.kinetic_energy, kinetic_energy(prob.markers))
flush(stdout)

# 4. Basis sizes in d dimensions (d = 1 is the package; d = 3 is 3d3v), for the defaults
#    of SplineRedistribution: nx = 8, nv = 10 cubic cells of the phase-space spline,
#    nx_cells = 10, np_cells = 8 quartic cells of the pairwise Hamiltonian.
H1 = info.hamiltonian
@printf("\n4. basis sizes (package, d = 1: N_θ = %d, M_x = %d, M_p = %d)\n",
    hamiltonian_parameter_count(H1), H1.Mx, H1.Mp)
@printf("   phase-space spline of the package: |I| = %d\n", length(info.coefficients))
for dim in 1:3
    I_size = 8^dim * (10 + 3)^dim                # periodic in x, clamped in v, degree 3
    Mx, Mp = 10^dim, (8 + 4)^dim                  # periodic quartic, clamped quartic
    Nθ = Mp + Mx * Mp^2
    @printf("   d = %d: |I| = %d, (k+1)^(2d) = %d per row of E; M_x = %d, M_p = %d, N_θ = %.4g; ",
        dim, I_size, 4^(2dim), Mx, Mp, Nθ)
    @printf("nonzero (m, l, r) per pair (k+1)^(3d) = %d; 2Nd - d at N = 16^(2d): %d\n",
        5^(3dim), 2 * 16^(2dim) * dim - dim)
end
flush(stdout)
