# R3: the symmetric ensemble Hamiltonians (translation and permutation invariance, the Fourier
# moment and centroid-frame forms) and Controls A, B and C.

using ParticleRemapping
using ForwardDiff
using LinearAlgebra
using Random
using Test

Random.seed!(1234)

# the verification set-up up to the targets (as in R2)
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
p̄_raw = kkt.raw_weights .* m .* v̄
p_scale = 2.5 * maximum(max.(abs.(p), abs.(p̄)))
H = PairwiseSplineHamiltonian(ex.Lx, p_scale; nx_cells = 10, np_cells = 8)
S = redistribution_scaling(N, p_scale)

q_mid, p_mid = 0.5 .* (q .+ q̄), 0.5 .* (p .+ p̄)
θ = solve_redistribution_parameters(H, ImplicitMidpointStep(), q, p, q̄, p̄)

@testset "Pairwise spline Hamiltonian: symmetry" begin
    @test hamiltonian_parameter_count(H) == 12 + 10 * 12^2
    Hq, Hp, _ = pairwise_hamiltonian_derivatives(H, q_mid, p_mid, θ)
    # Noether: Σ_a ∇_{q_a} H = 0
    @test abs(sum(Hq)) < 1e-15
    # translation invariance
    Hq_s, Hp_s, _ = pairwise_hamiltonian_derivatives(H, q_mid .+ 1.2345, p_mid, θ)
    @test norm(Hq_s - Hq) < 1e-14
    @test norm(Hp_s - Hp) < 1e-14
    # permutation equivariance
    perm = sortperm([sin(17.3 * k) for k in 1:N])
    Hq_p, Hp_p, _ = pairwise_hamiltonian_derivatives(H, q_mid[perm], p_mid[perm], θ)
    @test norm(Hq_p - Hq[perm]) < 1e-14
    @test norm(Hp_p - Hp[perm]) < 1e-14
    @test ensemble_hamiltonian(H, q_mid[perm], p_mid[perm], θ) ≈
          ensemble_hamiltonian(H, q_mid, p_mid, θ) rtol = 1e-13
end

@testset "Pairwise spline Hamiltonian: derivatives against automatic differentiation" begin
    z = vcat(q_mid, p_mid)
    Hfun(y) = ensemble_hamiltonian(H, y[1:N], y[(N + 1):(2N)], θ)
    Hq, Hp, Hzz = pairwise_hamiltonian_derivatives(H, q_mid, p_mid, θ)
    g = ForwardDiff.gradient(Hfun, z)
    @test vcat(Hq, Hp) ≈ g rtol = 1e-12
    @test Hzz ≈ ForwardDiff.hessian(Hfun, z) rtol = 1e-12
    @test Hzz ≈ transpose(Hzz)
    dq, dp = basis_hamiltonian_gradients(H, q_mid, p_mid)
    @test dq * θ ≈ Hq rtol = 1e-12
    @test dp * θ ≈ Hp rtol = 1e-12
    @test maximum(abs, sum(dq; dims = 1)) < 1e-15
end

@testset "Parameter system: rank 2Nd − d, Control C" begin
    for scheme in (ImplicitMidpointStep(), SymplecticEulerStep())
        G = redistribution_parameter_matrix(H, scheme, q, p, q̄, p̄)
        @test rank(Diagonal(S) * G) == 2N - 1
        # every column lies in M₀: no net momentum
        @test maximum(abs, sum(G[(N + 1):(2N), :]; dims = 1)) < 1e-15
        θs = solve_redistribution_parameters(H, scheme, q, p, q̄, p̄)
        @test norm(S .* (G * θs - vcat(q̄ .- q, p̄ .- p))) < 1e-13
    end
    # Control C: the raw momenta leave the compatible subspace, so G θ = Δz has no solution
    G_raw = redistribution_parameter_matrix(H, ImplicitMidpointStep(), q, p, q̄, p̄_raw)
    Δz_raw = vcat(q̄ .- q, p̄_raw .- p)
    θ_raw = G_raw \ Δz_raw
    @test norm(G_raw * θ_raw - Δz_raw) > 1e-4
    # the least-squares residual is the net momentum spread over N markers
    @test norm(G_raw * θ_raw - Δz_raw) ≈ abs(sum(p̄_raw) - sum(p)) / sqrt(N) rtol = 1e-8
end

@testset "Control A: uncoupled Hamiltonians change the total momentum" begin
    z_unc = uncoupled_control_step(vcat(q, p), p_scale)
    @test abs(sum(z_unc[(N + 1):(2N)]) - sum(p)) > 1e-2
end

@testset "Control B: mean-force subtraction breaks symplecticity" begin
    z0 = vcat(q, p)
    z_adhoc = mean_force_control_step(z0, p_scale)
    @test abs(sum(z_adhoc[(N + 1):(2N)]) - sum(p)) < 1e-15
    D = ForwardDiff.jacobian(z -> mean_force_control_step(z, p_scale), z0)
    J = canonical_J(N)
    @test norm(transpose(D) * J * D - J) > 1e-2
    # the vector field of Lemma (ii): ∂q̇_b/∂q_a = 0 but ∂ṗ_a/∂p_b ≠ 0, so J⁻¹ DX is not symmetric
    h_q(qa, pa) = -sin(qa) * (1 + pa / 2)
    h_p(qa, pa) = cos(qa) / 2
    X(z) = vcat(h_p.(z[1:N], z[(N + 1):(2N)]),
        -h_q.(z[1:N], z[(N + 1):(2N)]) .+ sum(h_q.(z[1:N], z[(N + 1):(2N)])) / N)
    DX = ForwardDiff.jacobian(X, z0)
    @test sum(DX[(N + 1):(2N), :]; dims = 1) ≈ zeros(1, 2N) atol = 1e-15
    @test norm(J * DX - transpose(J * DX)) > 1e-2
end

@testset "Fourier-moment form: O(N) factorisation" begin
    F = FourierMomentHamiltonian([1, 2, 3], p_scale; np_cells = 4)
    θF = randn(hamiltonian_parameter_count(F))
    @test hamiltonian_parameter_count(F) == 8 + 2 * 3 * 64
    E = fourier_pairwise_hamiltonian(F, q_mid, p_mid, θF)
    @test ensemble_hamiltonian(F, q_mid, p_mid, θF) ≈ E rtol = 1e-12
    Hq, Hp = ensemble_hamiltonian_gradient(F, q_mid, p_mid, θF)
    g = ForwardDiff.gradient(
        y -> fourier_pairwise_hamiltonian(F, y[1:N], y[(N + 1):(2N)], θF),
        vcat(q_mid, p_mid))
    @test vcat(Hq, Hp) ≈ g rtol = 1e-12
    @test abs(sum(Hq)) < 1e-12 * norm(Hq)
    perm = randperm(N)
    @test ensemble_hamiltonian(F, q_mid[perm] .+ 0.77, p_mid[perm], θF) ≈ E rtol = 1e-12
end

@testset "Centroid-frame form on ℝ: recoil-coupled gradients" begin
    C = CentroidFrameHamiltonian(4.0, p_scale; nq_cells = 6, np_cells = 4)
    θC = randn(hamiltonian_parameter_count(C))
    qc = q_mid .- 3.0                      # positions on ℝ, centroid away from 0
    E = ensemble_hamiltonian(C, qc, p_mid, θC)
    Hq, Hp = ensemble_hamiltonian_gradient(C, qc, p_mid, θC)
    g = ForwardDiff.gradient(y -> ensemble_hamiltonian(C, y[1:N], y[(N + 1):(2N)], θC),
        vcat(qc, p_mid))
    @test vcat(Hq, Hp) ≈ g rtol = 1e-12
    @test abs(sum(Hq)) < 1e-12 * norm(Hq)
    perm = randperm(N)
    @test ensemble_hamiltonian(C, qc[perm] .+ 1.9, p_mid[perm], θC) ≈ E rtol = 1e-12
    # the field J ∇H is Hamiltonian: the Jacobian of the gradient formula is symmetric
    D = ForwardDiff.jacobian(
        y -> vcat(ensemble_hamiltonian_gradient(C, y[1:N], y[(N + 1):(2N)],
            θC)...),
        vcat(qc, p_mid))
    @test D ≈ transpose(D) rtol = 1e-12
end
