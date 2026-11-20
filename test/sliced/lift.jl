using ParticleRemapping
using ForwardDiff: ForwardDiff
using LinearAlgebra: Symmetric, cond, eigmin, opnorm
using Random: MersenneTwister
using Test

# Three slices of four markers with targets of different spans, and anchors that differ per slice.
rng = MersenneTwister(7)
N = 12
u = 0.8 .* randn(rng, N)
v = randn(rng, N)
bins = sorted_slices(randn(rng, N), u, 3)
taus = [uniform_quantile_targets(L, 4) for L in (1.3, 0.7, 2.1)]
anchors = slice_anchors(bins, [-0.4, 0.1, 0.9])
scale = N * (maximum(abs, u) + 2.1)

function lifted(z; kw...)
    vcat(split_quantile_ot_lift(z[1:N], z[(N + 1):(2N)], bins, taus; kw...)...)
end
function dual_lifted(z; kw...)
    vcat(dual_quantile_ot_lift(z[1:N], z[(N + 1):(2N)], bins, taus; kw...)...)
end
defect(J) = maximum(abs, transpose(J) * canonical_J(N) * J - canonical_J(N))

@testset "potential, gradient and Hessian (β_mean = $βm, β_intra = $βi, ε = $ε)" for (βm, βi, ε) in Iterators.product(
    (0.55, 0.8), (0.6, 0.7), (0.3, 1.1))
    kw = (β_mean = βm, β_intra = βi, ε = ε)
    ũ, ṽ = split_quantile_ot_lift(u, v, bins, taus; kw...)
    S(x) = sliced_quantile_potential(x, bins, taus; kw...)
    # the lift is the gradient of the potential, its Jacobian the Hessian
    @test isapprox(ũ, ForwardDiff.gradient(S, u); atol = 8 * eps() * scale)
    H = sliced_quantile_hessian(u, bins, taus; kw...)
    @test isapprox(H, ForwardDiff.hessian(S, u); atol = 8 * eps() * scale)
    @test H == transpose(H)
    # ∇²S 1 = 1 and ∇²S ⪰ (1 − max β) I
    @test isapprox(H * ones(N), ones(N); atol = 8 * eps() * N)
    @test eigmin(Symmetric(H)) >= 1 - max(βm, βi) - 8 * eps() * N
    # a slice-constant vector with zero mean is an eigenvector with eigenvalue 1 − β_mean
    y = zeros(N)
    for (b, c) in zip(bins, (2.0, -1.5, -0.5))
        y[b] .= c
    end
    @test isapprox(H * y, (1 - βm) .* y; atol = 8 * eps() * N)
    # the cotangent part solves ∇²S ṽ = v, and the sums of u and v are kept
    @test isapprox(H * ṽ, v; atol = 8 * eps() * N * cond(H))
    @test isapprox(sum(ũ), sum(u); atol = 8 * eps() * scale)
    @test isapprox(sum(ṽ), sum(v); atol = 8 * eps() * N * cond(H))
    # translation equivariance: u + c ↦ ũ + c, ṽ unchanged
    ũc, ṽc = split_quantile_ot_lift(u .+ 0.37, v, bins, taus; kw...)
    @test isapprox(ũc, ũ .+ 0.37; atol = 8 * eps() * scale)
    @test isapprox(ṽc, ṽ; atol = 8 * eps() * N * cond(H))
end

@testset "symplecticity of the lifts" begin
    z = vcat(u, v)
    for kw in ((;), (; anchors = anchors), (;
        anchors = anchors, β_mean = 0.8, β_intra = 0.6))
        J = ForwardDiff.jacobian(x -> lifted(x; kw...), z)
        tol = 8 * eps() * N * cond(J) * opnorm(J)^2
        @test defect(J) < tol
        Jd = ForwardDiff.jacobian(x -> dual_lifted(x; kw...), z)
        @test defect(Jd) < 8 * eps() * N * cond(Jd) * opnorm(Jd)^2
    end
    # anchored lift: ṽ − v_anc = (∇²S)⁻¹ (v − v_anc), so Σ ṽ = Σ v for any anchors
    ũ, ṽ = split_quantile_ot_lift(u, v, bins, taus; anchors = anchors)
    H = sliced_quantile_hessian(u, bins, taus)
    @test isapprox(H * (ṽ .- anchors), v .- anchors; atol = 8 * eps() * N * cond(H))
    @test isapprox(sum(ṽ), sum(v); atol = 8 * eps() * N * cond(H))
    # control: without the cotangent part the map is not symplectic
    Jc = ForwardDiff.jacobian(x -> lifted(x; with_cotangent = false), z)
    @test defect(Jc) > 1e-1
    @test lifted(z; with_cotangent = false)[(N + 1):(2N)] == v
end

@testset "quantile limit" begin
    # β_mean = β_intra = 1 and ε → 0: the markers of a slice, in the order of the targets, go to
    # the mean plus the targets, ũ_a = ū + τ_r + O(ε/L), since κ → 1 and the mollified sign → ±1.
    ε = 1e-9
    ũ, _ = split_quantile_ot_lift(u, v, bins, taus; β_mean = 1.0, β_intra = 1.0, ε = ε,
        with_cotangent = false)
    expected = zeros(N)
    for (b, idx) in pairs(bins), (r, a) in pairs(idx)

        expected[a] = sum(u) / N + taus[b][r]
    end
    @test isapprox(ũ, expected; atol = 10 * ε / 0.35)
    # the mutant that keeps each slice mean instead of the global mean misses it
    mutant = similar(expected)
    for (b, idx) in pairs(bins), (r, a) in pairs(idx)

        mutant[a] = sum(u[idx]) / length(idx) + taus[b][r]
    end
    @test maximum(abs, ũ - mutant) > 1e-3
end
