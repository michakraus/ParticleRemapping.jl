# R4: the forward Newton solves of the implicit midpoint and the Type-2 symplectic Euler step
# with the solved parameters: convergence to the targets, symplecticity (by the implicit
# function theorem and by finite differences), momentum conservation, and the caustic control.

using ParticleRemapping
using LinearAlgebra
using Test

ex = redistribution_example()
mk = ex.markers
q, v, w, m = mk.x, mk.v, mk.w, mk.m
p = w .* m .* v
N = length(q)
ps = PhaseSpaceSplines()
relaxed = declustering_relaxation(q, v; Lx = ex.Lx, σx = ex.σx, σv = ex.σv)
q̄, v̄ = relaxed.x, relaxed.v
ΔΩ = number_density_volumes(ps, q̄, v̄)
kkt = kkt_redistribution_weights(
    ps, lumped_spline_projection(ps, q, v, w), q̄, v̄, ΔΩ, sum(w),
    sum(p); m = m)
p̄ = kkt.weights .* m .* v̄
p_scale = 2.5 * maximum(max.(abs.(p), abs.(p̄)))
H = PairwiseSplineHamiltonian(ex.Lx, p_scale; nx_cells = 10, np_cells = 8)
S = redistribution_scaling(N, p_scale)
J = canonical_J(N)
z̄ = vcat(q̄, p̄)

θ_mid = solve_redistribution_parameters(H, ImplicitMidpointStep(), q, p, q̄, p̄)
θ_se = solve_redistribution_parameters(H, SymplecticEulerStep(), q, p, q̄, p̄)
mid = forward_redistribution_step(H, ImplicitMidpointStep(), q, p, θ_mid)
se = forward_redistribution_step(H, SymplecticEulerStep(), q, p, θ_se)

# fourth-order differences with steps 1e-6 (positions) and 1e-6 p_s (momenta)
fd_steps = vcat(fill(1e-6, N), fill(1e-6 * p_scale, N))
function fd4(F, z)
    D = zeros(length(F(z)), length(z))
    for j in eachindex(z)
        e = zeros(length(z))
        e[j] = fd_steps[j]
        D[:, j] = (-F(z .+ 2e) .+ 8.0 .* F(z .+ e) .- 8.0 .* F(z .- e) .+ F(z .- 2e)) ./
                  (12.0 * e[j])
    end
    D
end

@testset "Implicit midpoint: one step reaches the targets" begin
    @test mid.residual < 1e-13
    @test norm(S .* (mid.z - z̄)) < 1e-13
    @test norm(mid.z - vcat(q, p)) > 0.1                      # the step moves the markers
    @test abs(sum(mid.z[(N + 1):(2N)]) - sum(p)) < 1e-15
    @test norm(transpose(mid.jacobian) * J * mid.jacobian - J) < 1e-13
    D_fd = fd4(
        z -> forward_redistribution_step(H, ImplicitMidpointStep(), z[1:N],
            z[(N + 1):(2N)], θ_mid).z,
        vcat(q, p))
    @test norm(mid.jacobian - D_fd) / norm(mid.jacobian) < 1e-9
    @test norm(transpose(D_fd) * J * D_fd - J) < 1e-7
    # the contraction bound is sufficient, not necessary: it is violated here, and the
    # homotopy still reaches the target branch
    @test redistribution_contraction_number(H, ImplicitMidpointStep(), q, p, q̄, p̄, θ_mid) >
          1
end

@testset "Symplectic Euler: symplectic step, caustic at one step (control)" begin
    @test se.residual < 1e-13
    @test abs(sum(se.z[(N + 1):(2N)]) - sum(p)) < 1e-15
    @test norm(transpose(se.jacobian) * J * se.jacobian - J) < 1e-13
    D_fd = fd4(
        z -> forward_redistribution_step(H, SymplecticEulerStep(), z[1:N],
            z[(N + 1):(2N)], θ_se).z,
        vcat(q, p))
    @test norm(se.jacobian - D_fd) / norm(se.jacobian) < 1e-9
    @test norm(transpose(D_fd) * J * D_fd - J) < 1e-7
    # ‖τ ∇_q ∇_p̄ H‖ > 1: the identity-connected branch ends on another root, not on z̄
    @test redistribution_contraction_number(H, SymplecticEulerStep(), q, p, q̄, p̄, θ_se) > 1
    @test norm(S .* (se.z - z̄)) > 0.1
    # z̄ is a root of the same equations: the step equations hold there to round-off
    Hq, Hp, _ = pairwise_hamiltonian_derivatives(H, q, p̄, θ_se)
    @test norm(p̄ .- p .+ Hq) / p_scale < 1e-13
    @test norm(q .+ Hp .- q̄) < 1e-13
end

@testset "Substeps restore the contraction bound" begin
    for scheme in (SymplecticEulerStep(), ImplicitMidpointStep())
        K = scheme isa SymplecticEulerStep ? 2 : 3
        sub = redistribution_substeps(H, scheme, q, p, q̄, p̄; substeps = K)
        @test length(sub.θ) == K
        @test all(<(1), sub.contraction)
        @test all(<(1e-13), sub.residual)
        @test norm(S .* (sub.z - z̄)) < 1e-13
        @test abs(sum(sub.z[(N + 1):(2N)]) - sum(p)) < 1e-15
    end
    # the single step: the bound fails for both schemes
    one = redistribution_substeps(H, SymplecticEulerStep(), q, p, q̄, p̄)
    @test only(one.contraction) > 1
end

@testset "Pseudo-time path" begin
    path = pseudotime_path(
        H, ImplicitMidpointStep(), q, p, θ_mid; s = range(0, 1; length = 5))
    @test path[1] == vcat(q, p)
    @test norm(S .* (path[end] - z̄)) < 1e-13
    @test all(abs(sum(z[(N + 1):(2N)]) - sum(p)) < 1e-15 for z in path)
end
