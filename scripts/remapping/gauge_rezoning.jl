# The numerical example of the gauge rezoning: 6400 clustered markers on 16 nodes, with the
# constrained moments K = {0, 2} and K = {0, 2, 3}, and the unconstrained targets as a control.
# Prints the numbers of the table of the gauge rezoning and checks them.
#
# Run, from the package root:
#   julia --startup-file=no --project=scripts scripts/remapping/gauge_rezoning.jl

using ParticleRemapping
using Printf

results = Pair{String, Bool}[]
function check(name, ok, detail = "")
    push!(results, name => ok)
    println(rpad(name, 58), detail, ok ? "  ok" : "  UNEXPECTED")
end

setup = gauge_rezoning_setup()
mk, grid = setup.markers, setup.grid
function relres(after, k)
    maximum(abs.(grid_moment(grid, after, k) .- grid_moment(grid, mk, k)) ./
            grid_moment(grid, mk.x, abs.(mk.v), mk.w, k))
end
gap(after) = mean_gap_ratio(grid, mk.x, after.v)

t = @elapsed r = remap(setup.method, mk)
after = r.markers
@printf("K = {0,2}: τ = %.4f, min H_a/w_a = %.2e, %.1f s\n", r.info.τ, r.info.min_curvature,
    t)
check("(1) charge constraint, K = {0,2}", relres(after, 0) < 1e-13,
    @sprintf("%.1e", relres(after, 0)))
check("(1) energy constraint, K = {0,2}", relres(after, 2) < 1e-13,
    @sprintf("%.1e", relres(after, 2)))
check("(2) current, not constrained", relres(after, 1) < 1e-14,
    @sprintf("%.1e", relres(after, 1)))
δp = maximum(abs.(momentum(after) .- momentum(mk)) ./ abs.(momentum(mk)))
check("(2) p_a = m W_a V_a for every marker", δp < 8eps(), @sprintf("%.1e", δp))
δE = abs(kinetic_energy(after) / kinetic_energy(mk) - 1)
check("(2) global kinetic energy", δE < 1e-12, @sprintf("%.1e", δE))
check("(3) velocity gaps more even", gap(after) < gap(mk),
    @sprintf("%.2f -> %.2f", gap(mk), gap(after)))
check("    signs of v kept", all(sign.(after.v) .== sign.(mk.v)))
@printf("    stretch range: λ ∈ [%.3f, %.3f]\n", extrema(exp.(r.info.μ))...)
check("Control: K = {0,2} changes the third moment", relres(after, 3) > 1e-3,
    @sprintf("%.2f", relres(after, 3)))

unconstrained = from_canonical(
    mk.x, momentum(mk), mk.w .* exp.(-gauge_targets(grid, mk.x, mk.v));
    m = mk.m)
check("Control: unconstrained targets violate the charge", relres(unconstrained, 0) > 1e-3,
    @sprintf("%.2f", relres(unconstrained, 0)))
@printf("    unconstrained targets: largest/mean gap %.2f\n", gap(unconstrained))

r3 = remap(gauge_rezoning_setup(; K = (0, 2, 3)).method, mk)
res3 = maximum(relres(r3.markers, k) for k in (0, 2, 3))
@printf("K = {0,2,3}: τ = %.4f, largest/mean gap %.2f\n", r3.info.τ, gap(r3.markers))
check("(4) K = {0,2,3}: all three constraints", res3 < 1e-13, @sprintf("%.1e", res3))

println("\nTable (gauge rezoning):")
@printf("  τ reached                 %.3f   %.3f   1\n", r.info.τ, r3.info.τ)
@printf("  largest/mean gap (%.1f)   %.1f    %.1f    %.2f\n", gap(mk), gap(after),
    gap(r3.markers), gap(unconstrained))
@printf("  charge residual           %.1e  %.1e  %.2f\n", relres(after, 0),
    relres(r3.markers, 0), relres(unconstrained, 0))
@printf("  kinetic energy residual   %.1e  %.1e\n", relres(after, 2), relres(r3.markers, 2))
@printf("  current residual          %.1e\n", relres(after, 1))
@printf("  change of M^(3)           %.2f     %.1e\n", relres(after, 3),
    relres(r3.markers, 3))
@printf("  range of λ                [%.2f, %.1f]\n", extrema(exp.(r.info.μ))...)

nfail = count(!last, results)
println("\n", length(results) - nfail, "/", length(results), " checks as expected")
nfail == 0 || error("$nfail check(s) not as expected")
