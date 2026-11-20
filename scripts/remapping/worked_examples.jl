# The four worked examples of the Hamiltonian remapping: one charge constraint in closed form and
# its fold, eight markers on two nodes step by step, the admissible cotangent lift in closed form,
# and the current-quantile matching of five markers. Prints every number of their tables and
# checks the claims, each with a control that must fail.
#
# Run, from the package root:
#   julia --startup-file=no --project=scripts scripts/remapping/worked_examples.jl

using ParticleRemapping
using Printf

results = Pair{String, Bool}[]
function check(name, ok, detail = "")
    push!(results, name => ok)
    println(rpad(name, 62), detail, ok ? "  ok" : "  UNEXPECTED")
end
row(name, x; d = 4) = println("    ", rpad(name, 10), " = ", round.(x; digits = d))

println("Example 1: charge constraint only, one node")
ex1 = charge_constraint_example()
s = charge_stationary_point(ex1)
row("v", ex1.v)
row("V*", ex1.Vstar)
row("μ*", ex1.μstar)
@printf("    A = %.5f, A_lim = %.5f\n", s.A, charge_multiplier_limit(ex1))
foreach(k -> row(string(k), getproperty(s, k)), (:μ, :λ, :V, :W, :Hw))
@printf("    Σ w e^{-μ*} = %.4f, Σ w v² = %.4f, Σ W V² = %.4f\n",
    sum(ex1.w .* exp.(-ex1.μstar)), sum(ex1.w .* ex1.v .^ 2), sum(s.W .* s.V .^ 2))
stat = maximum(abs.((s.μ .- ex1.μstar) ./ ex1.s2 .- s.A .* exp.(-s.μ)))
check("(1a) Lambert-W stationarity residual", stat < 8eps(), @sprintf("%.1e", stat))
check("(1a) charge constraint Σ W = Σ w", abs(sum(s.W) - sum(ex1.w)) < 8eps(),
    @sprintf("%.1e", abs(sum(s.W) - sum(ex1.w))))
ct = abs(sum(ex1.w .* exp.(-ex1.μstar)) - sum(ex1.w))
check("Control: targets violate the charge constraint", ct > 1e-2, @sprintf("%.3f", ct))
βf = charge_fold(ex1)
sf = charge_stationary_point(ex1, βf)
@printf("    fold at β_f = %.4f: A = %.6f, min H/w = %.1e, μ = %s\n",
    βf, sf.A, minimum(sf.Hw),
    string(round.(sf.μ; digits = 4)))
check("(1b) min H_a/w_a → 0 at the fold", minimum(sf.Hw) < 1e-6, @sprintf("%.1e",
    minimum(sf.Hw)))
check("(1b) no stationary point beyond the fold", charge_stationary_point(ex1, 1.001βf) ===
                                                  nothing)

println("\nExample 2: eight markers, two nodes, K = {0, 2}")
data = eight_marker_homotopy()
info, mk = data.result.info, data.example.markers
prob = info.problem
row("μ*", info.μstar)
row("V*", mk.v .* exp.(info.μstar))
cstar = gauge_constraints(prob, info.μstar)[1]
row("c(μ*)", cstar)
println("    full targets from μ = 0, without homotopy:")
for l in data.direct.iterations
    @printf("      it %d: max|c|/scale = %.2e, max|Δμ| = %.2e, min H/w = %.4f\n", l.it,
        l.residual, l.step, l.min_curvature)
end
check("Control: Newton at τ = 1 without homotopy is rejected", !data.direct.converged)
δ = maximum(l.solver_difference for h in info.history for l in h.iterations; init = 0.0)
δ = max(δ, data.direct.iterations[1].solver_difference)
check("(2a) Schur-complement step equals the KKT step", δ < 8eps(), @sprintf("%.1e", δ))
for h in info.history[2:end]
    @printf("    τ = %.4f%s\n", h.τ, h.accepted ? "" : " (rejected)")
    for l in h.iterations
        @printf("      it %d: max|c|/scale = %.2e, max|Δμ| = %.2e, min H/w = %.4f, t = %.3g\n",
            l.it, l.residual, l.step, l.min_curvature, l.t)
    end
end
check("(2b) homotopy reaches τ = 1", info.τ == 1,
    @sprintf("in %d steps", count(h -> h.accepted, info.history) - 1))
c = gauge_constraints(prob, info.μ)[1]
st = maximum(abs, gauge_stationarity(prob, info.μ, info.α, 1.0))
check("(2b) stationarity and constraints",
    st < 1e-12 && maximum(abs.(c) ./ prob.scale) < 1e-14,
    @sprintf("stat %.1e, constr %.1e", st, maximum(abs.(c) ./ prob.scale)))
check("(2b) every H_a > 0", info.min_curvature > 0, @sprintf("min H/w = %.4f",
    info.min_curvature))
after = data.result.markers
row("α", info.α)
row("μ", info.μ)
row("λ", exp.(info.μ))
row("V", after.v)
row("W", after.w)
row("H/w", gauge_hessian(prob, info.μ, info.α) ./ prob.w; d = 3)
cur = maximum(abs.(grid_moment(data.example.grid, after, 1) .-
                   grid_moment(data.example.grid, mk, 1)))
check("(2c) grid current unchanged", cur < 8eps(), @sprintf("%.1e", cur))
check("(2c) every momentum m W V = p", maximum(abs.(momentum(after) .- momentum(mk))) <
                                       8eps())
check(
    "Control: targets violate the constraints", maximum(abs.(cstar) ./ prob.scale) > 1e-2,
    @sprintf("%.3f", maximum(abs.(cstar) ./ prob.scale)))
@printf("    gap ratio cell 1: %.2f -> %.2f, cell 2: %.2f -> %.2f\n",
    speed_gap_ratio(mk.v[1:4]), speed_gap_ratio(after.v[1:4]),
    speed_gap_ratio(mk.v[5:8]), speed_gap_ratio(after.v[5:8]))

println("\nExample 3: cotangent lift for f = J(x) M(v), equal weights")
J(x) = cosine_profile(x)
xs = [0.0, π / 2, π, 4.0]
devs = Float64[]
for x in xs
    X, dφ = admissible_lift(x, 1.5)
    push!(devs, abs(J(X) * dφ^2 / J(x) - 1))
    @printf("    x = %.4f -> X = %.4f, φ' = %.4f, W/w = 1/φ' = %.4f\n", x, X, dφ, 1 / dφ)
end
check("(3a) J(φ_t) φ_t'² = J for the closed-form flow", maximum(devs) < 8eps(),
    @sprintf("%.1e", maximum(devs)))
Vdev = maximum(abs(admissible_lift_matching(x, v, 1.5, 0.01).V - v)
for x in xs,
v in (-1.3, 0.4, 2.2))
check("(3a) matched velocity V = v", Vdev < 1e-14, @sprintf("%.1e", Vdev))
th = lift_theta(2π) / 2
xg = range(0, 2π; length = 401)
dens = [sqrt(J(x) * J(admissible_lift(x, th)[1])) for x in xg]
@printf("    t = Θ(2π)/2 = %.4f: max/min of the marker density %.4f (t = 0: %.4f)\n", th,
    maximum(dens) / minimum(dens), maximum(J.(xg)) / minimum(J.(xg)))
check("(3b) half-period lift: marker density max/min < 1.2", maximum(dens) / minimum(dens) <
                                                             1.2)
check("Control: at t = 0 the marker density ∝ J is uneven", maximum(J.(xg)) /
                                                            minimum(J.(xg)) > 2)
R = [lift_compatibility_residual(J, x,
         cotangent_lift_flow(x -> 1 + 0.5cos(x), x -> -0.5sin(x), x, 1.5)...) for x in xs]
row("R(x), a = 1 + 0.5 cos x", R; d = 2)
check("Control: a = 1 + 0.5 cos x is incompatible", maximum(abs.(R)) > 1e-2,
    @sprintf("%.2f", maximum(abs.(R))))

println("\nExample 4: current-quantile matching in one cell")
e = quantile_matching_example()
row("|P|", e.aP)
row("Π", e.Π; d = 5)
@printf("    Σ|P|/(mΔx) = %.5f, F̂(∞) = %.5f, R = %.5f\n", sum(e.aP) / (e.m * e.Δx),
    e.Jbar * maxwellian_current(Inf), e.R)
row("Π/(1+R)", e.Πc; d = 5)
row("V", e.V)
row("W", e.W)
@printf("    Σ W = %.3f, Σ W V² = %.3f; Σ w = %.3f, Σ w v² = %.3f\n", sum(e.W),
    sum(e.W .* e.V .^ 2), 5e.w0, sum(e.w0 .* e.v .^ 2))
check("(4a) m W V = P", maximum(abs.(e.m .* e.W .* e.V .- e.aP)) < 8eps())
dm = maximum(abs.(e.Jbar .* maxwellian_current.(e.V) .- e.Πc))
check("(4a) F̂(V_r) = Π_r/(1+R)", dm < 8eps(), @sprintf("%.1e", dm))
do_ = maximum(abs.(e.Jbar .* maxwellian_current.(e.v) .- e.Πc))
check("Control: the old velocities miss the quantiles", do_ > 1e-2, @sprintf("%.3f", do_))

nfail = count(!last, results)
println("\n", length(results) - nfail, "/", length(results), " checks as expected")
nfail == 0 || error("$nfail check(s) not as expected")
