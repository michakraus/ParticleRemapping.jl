# The numerical example of the canonical transport: 2e5 markers sampled from
# f = (1 + 0.4 cos x) M(v), the admissible cotangent lift of a = J^{-1/2} and the inadmissible
# lift of a = 1 + 0.5 cos x at t = 1.5, each followed by the current-quantile matching. Prints the
# numbers of the table of the canonical transport and checks them.
#
# Run, from the package root:
#   julia --startup-file=no --project=scripts scripts/remapping/current_quantile_matching.jl

using ParticleRemapping
using Printf

results = Pair{String, Bool}[]
function check(name, ok, detail = "")
    push!(results, name => ok)
    println(rpad(name, 62), detail, ok ? "  ok" : "  UNEXPECTED")
end

setup = quantile_matching_setup()
t = @elapsed ex = quantile_matching_experiment(setup)
@printf("%d markers, four configurations in %.1f s\n", length(ex.original), t)
err(mk) = binned_weight_error(setup.target, setup.grid, mk.x, mk.v, mk.w)
e_orig, e_id, e_c = err(ex.original), err(ex.identity.markers), err(ex.admissible.markers)
e_kept, e_b = err(ex.kept), err(ex.inadmissible.markers)
r_id, r_c, r_b = (getproperty(ex, k).info.max_residual
for k in (:identity, :admissible, :inadmissible))
d = maximum(abs.(mod.(ex.admissible.info.X .- ex.original.x .+ π, 2π) .- π))

@printf("    admissible lift moves markers by up to %.2f in x\n", d)
check("(c) ψ = id + matching: error below sampling noise", e_id < e_orig,
    @sprintf("%.4f < %.4f", e_id, e_orig))
check("(c) admissible ψ + matching: error below sampling noise", e_c < e_orig,
    @sprintf("%.4f < %.4f", e_c, e_orig))
check("Control: admissible ψ, no matching (W = w): error large", e_kept > 3e_orig,
    @sprintf("%.4f", e_kept))
check("Control: inadmissible ψ + matching: residual large", r_b > 10r_c,
    @sprintf("%.1f", r_b))
mk, P = ex.admissible.markers, ex.admissible.info.P
δ = maximum(abs.(mk.m .* mk.w .* mk.v .- P) ./ abs.(P))
check("    gauge identity m W V = P (admissible case)", δ < 8eps(), @sprintf("%.1e", δ))

println("\nTable (canonical transport with current-quantile matching):")
@printf("  original markers (sampling noise)             %.4f   --\n", e_orig)
@printf("  ψ = id, then matching                         %.4f   %.3f\n", e_id, r_id)
@printf("  admissible ψ, then matching                   %.4f   %.3f\n", e_c, r_c)
@printf("  admissible ψ, old weights kept (W = w)        %.4f   --\n", e_kept)
@printf("  inadmissible ψ, then matching                 %.3f    %.1f\n", e_b, r_b)

nfail = count(!last, results)
println("\n", length(results) - nfail, "/", length(results), " checks as expected")
nfail == 0 || error("$nfail check(s) not as expected")
