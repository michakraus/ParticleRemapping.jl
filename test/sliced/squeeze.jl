using ParticleRemapping
using ForwardDiff: ForwardDiff
using LinearAlgebra: Diagonal, I, cond, det, opnorm
using Random: MersenneTwister
using Test

rng = MersenneTwister(11)
N = 10
q, p = randn(rng, N), randn(rng, N)
q_anc = uniform_quantile_targets(2.0, N)
raw = 0.3 .+ rand(rng, N)
μ = (N / sum(raw)) .* raw
M = Diagonal(μ) * (I + (1 / N) * ones(N) * transpose(1 .- μ))

@testset "Sherman–Morrison form of the squeeze" begin
    q̄, p̄ = canonical_weight_squeeze(q, p, μ, q_anc)
    tol = 8 * eps() * N * cond(M) * (maximum(abs, p) + maximum(abs, q .- q_anc))
    @test isapprox(p̄, M * p; atol = tol)
    @test isapprox(q̄ .- q_anc, transpose(inv(M)) * (q .- q_anc); atol = tol)
    # det(I + (1/N) 1 (1 − μ)ᵀ) = 1 + (N − Σμ)/N = 1, so det M = Π μ
    @test isapprox(det(M), prod(μ); rtol = 8 * eps() * N * cond(M))
    @test isapprox(sum(p̄), sum(p); atol = 8 * eps() * N * maximum(abs, p̄))
    # the squeeze multiplies p by μ_a around the common shift (1/N) Σ (1 − μ_c) p_c
    @test isapprox(p̄ ./ μ .- p, fill(sum((1 .- μ) .* p) / N, N); atol = tol)
end

@testset "symplecticity" begin
    sq(z) = vcat(canonical_weight_squeeze(z[1:N], z[(N + 1):(2N)], μ, q_anc)...)
    J = ForwardDiff.jacobian(sq, vcat(q, p))
    @test maximum(abs, transpose(J) * canonical_J(N) * J - canonical_J(N)) <
          8 * eps() * N * cond(J) * opnorm(J)^2
    # J is the block matrix of the linear cotangent lift
    @test isapprox(J, [transpose(inv(M)) zeros(N, N); zeros(N, N) M];
        atol = 8 * eps() * N * cond(M) * opnorm(M))
end

@testset "Σμ = N is needed" begin
    μ_bad = 1.1 .* μ
    @test_throws r"^ArgumentError: the slot multipliers must sum to N = 10" canonical_weight_squeeze(
        q, p, μ_bad, q_anc)
    # for Σμ = N (1 + δ) the matrix M changes the total momentum by δ Σ (1 − μ_c) p_c
    M_bad = Diagonal(μ_bad) * (I + (1 / N) * ones(N) * transpose(1 .- μ_bad))
    δ = sum(μ_bad) / N - 1
    @test isapprox(sum(M_bad * p) - sum(p), δ * sum((1 .- μ_bad) .* p);
        atol = 8 * eps() * N * opnorm(M_bad) * maximum(abs, p))
    @test abs(sum(M_bad * p) - sum(p)) > 1e-3
end
