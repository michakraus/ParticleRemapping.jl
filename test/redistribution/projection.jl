# R2: the L² projection (continuous moments, Control D), the positive weighted KKT projection
# (positive weights, exact discrete moments) and the raw evaluation (Control C).

using ParticleRemapping
using LinearAlgebra
using Test

ex = redistribution_example()
mk = ex.markers
x, v, w, m = mk.x, mk.v, mk.w, mk.m
p = w .* m .* v
ps = PhaseSpaceSplines()

@testset "L² projection: continuous moments, negative overshoot (Control D)" begin
    l2 = spline_moments(ps, l2_spline_projection(ps, x, v, w); m = m)
    @test abs(l2.mass - sum(w)) < 1e-14
    @test abs(l2.momentum - sum(p)) < 1e-14
    # Control D: the L² spline is negative in the sparse regions
    @test l2.minimum < -0.1
    # the lumped prior is positive everywhere
    @test minimum(lumped_spline_projection(ps, x, v, w)) > 0
end

@testset "Weighted KKT projection: positive weights, exact discrete moments" begin
    relaxed = declustering_relaxation(x, v; Lx = ex.Lx, σx = ex.σx, σv = ex.σv)
    x̄, v̄ = relaxed.x, relaxed.v
    ΔΩ = number_density_volumes(ps, x̄, v̄)
    @test all(>(0), ΔΩ)
    c_pos = lumped_spline_projection(ps, x, v, w)
    E = spline_evaluation_matrix(ps, x̄, v̄)
    A = kkt_moment_matrix(E, ΔΩ, v̄; m = m)
    @test rank(A) == 2
    kkt = kkt_redistribution_weights(ps, c_pos, x̄, v̄, ΔΩ, sum(w), sum(p); m = m)
    w̄ = kkt.weights
    @test A * kkt.coefficients ≈ [sum(w), sum(p)] rtol = 1e-14
    @test minimum(w̄) > 0
    @test abs(sum(w̄) - sum(w)) < 1e-14
    @test abs(sum(w̄ .* m .* v̄) - sum(p)) < 1e-14
    # Control C: the raw evaluation of the prior misses the discrete momentum
    @test abs(sum(kkt.raw_weights .* m .* v̄) - sum(p)) > 1e-4
    # stationarity: the gradient W (c − c_pos) of the objective is orthogonal to every feasible
    # direction, the null space of A; the difference c − c_pos carries the round-off
    # eps |c_I|, which W = diag(c_pos)⁻² amplifies to eps / c_pos_I
    Δc = kkt.coefficients .- c_pos
    Z = nullspace(A)
    g = Δc ./ c_pos .^ 2
    @test norm(transpose(Z) * g) <= 16 * eps() * norm(1 ./ c_pos)
    # the branch recovery of the coarse-basis regime: from a perturbed start, Newton's method
    # on m ΔΩ f_h(q̄, v) v = p̄ returns the targets on their branch
    p̄ = w̄ .* m .* v̄
    rec = recover_velocity_branch(ps, kkt.coefficients, x̄, p̄, v̄ .+ 0.02, ΔΩ; m = m)
    @test rec.v ≈ v̄ rtol = 1e-12
    @test rec.w ≈ w̄ rtol = 1e-12
    # the slope is ∂_v (m ΔΩ f_h v): compare with a central difference (step 1e-5, truncation
    # O(h²) relative)
    h = 1e-5
    P(vs) = m .* ΔΩ .* evaluate_phase_space_spline(ps, kkt.coefficients, x̄, vs) .* vs
    @test rec.slope ≈ (P(v̄ .+ h) .- P(v̄ .- h)) ./ 2h rtol = 1e-6
    # a feasible competitor has a larger objective
    d = Z * normalize(ones(size(Z, 2)))
    objective(c) = sum(abs2, (c .- c_pos) ./ c_pos) / 2
    @test objective(kkt.coefficients .+ 1e-3 .* d) > objective(kkt.coefficients)
end
