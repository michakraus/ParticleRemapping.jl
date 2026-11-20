# R5: the closed-loop identity Φ*ω − ω = −τ Σ d g_k ∧ dΘ_k for a target-matching parameter map
# (Control E), the variational closure that restores Φ*ω = ω, and the fixed-grid no-go.

using ParticleRemapping
using ForwardDiff
using LinearAlgebra
using Test

ex = closed_loop_example()
z, H = ex.z, ex.H
n = length(z) ÷ 2
J = canonical_J(n)

# fourth-order differences with the steps of the note's check: 1e-5 and 1e-5 p_s
fd_steps = vcat(fill(1e-5, n), fill(1e-5 * H.p_scale, n))
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

@testset "Control E: target matching, implicit midpoint" begin
    cl = target_matching_closed_loop(H, z)
    lhs, rhs = closed_loop_defect(cl.jacobians...)
    @test norm(lhs) > 1e-3
    # exact Jacobians: the identity holds to round-off, amplified by the linear solves
    @test norm(lhs - rhs) / norm(lhs) < 1e3 * eps()
    # the note's check with finite differences
    DΦ = fd4(y -> target_matching_closed_loop(H, y; jacobians = false).z, z)
    DΘ = fd4(y -> target_matching_parameters(H, y), z)
    Dg = fd4(y -> target_matching_closed_loop(H, y; jacobians = false).g, z)
    lhs_fd, rhs_fd = closed_loop_defect(DΦ, DΘ, Dg)
    @test norm(lhs_fd) > 1e-3
    @test norm(lhs_fd - rhs_fd) / norm(lhs_fd) < 1e-8
    @test norm(lhs_fd - lhs) / norm(lhs) < 1e-8
    # at fixed parameters the same step is symplectic
    step = forward_redistribution_step(H, ImplicitMidpointStep(), z[1:n], z[(n + 1):(2n)],
        cl.θ)
    @test norm(transpose(step.jacobian) * J * step.jacobian - J) < 1e-13
    @test step.z ≈ cl.z rtol = 1e-14
end

@testset "Control E: target matching, symplectic Euler" begin
    cl = target_matching_closed_loop(H, z; scheme = SymplecticEulerStep())
    lhs, rhs = closed_loop_defect(cl.jacobians...)
    @test norm(lhs) > 1e-3
    @test norm(lhs - rhs) / norm(lhs) < 1e3 * eps()
end

@testset "Variational closure: closed-loop symplecticity" begin
    vl = variational_closed_loop(z)
    DΦ, DΘ = vl.jacobians
    # the parameters move with the state, and still Φ*ω = ω
    @test norm(DΘ) > 1e-2
    @test norm(transpose(DΦ) * J * DΦ - J) < 1e-12
    @test abs(sum(vl.z[(n + 1):(2n)]) - sum(z[(n + 1):(2n)])) < 1e-15
    @test norm(vl.z - z) > 1e-2                                  # the step moves the markers
    # finite differences reach the same Jacobian to their truncation and round-off
    DΦ_fd = fd4(y -> variational_closed_loop(y; jacobians = false).z, z)
    @test norm(DΦ_fd - DΦ) / norm(DΦ) < 1e-8
    @test norm(transpose(DΦ_fd) * J * DΦ_fd - J) < 1e-8
    # stationarity: θ = θ₀ + Ψ(z_{1/2}) / κ
    zm = 0.5 .* (z .+ vl.z)
    Ψ, _, _ = cotangent_lift_modes(zm[1:n], zm[(n + 1):(2n)], 3)
    @test vl.θ ≈ [0.06, 0.03, 0.015] .+ Ψ ./ 2 rtol = 1e-14
end

@testset "Corollary (i): positions reset onto a fixed grid" begin
    q_grid = collect(range(0, 2π; length = n + 1)[1:n])
    D = [zeros(n, 2n); zeros(n, n) I(n)]
    @test D ≈ ForwardDiff.jacobian(y -> fixed_grid_closed_loop(y, q_grid), z)
    @test transpose(D) * J * D == zeros(2n, 2n)
end
