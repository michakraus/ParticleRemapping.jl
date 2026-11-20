using ParticleRemapping
using ForwardDiff: ForwardDiff
using LinearAlgebra: cond, opnorm
using Random: MersenneTwister
using Test

rng = MersenneTwister(31)
d, N = 2, 36
x = 2π .* rand(rng, d, N) .- π
v = 0.8 .* randn(rng, d, N)
w = 0.5 .+ rand(rng, N)
mk = Markers(x, v, w; m = 1.1)
method = SlicedQuantileRemap(; B_v = 3, B_x = 6)
result = remap(method, mk)
info = result.info

@testset "per-axis lift" begin
    bins = [sorted_slices(v[ℓ, :], x[ℓ, :], 3) for ℓ in 1:d]
    taus = [uniform_quantile_targets(π, 12) for _ in 1:3]
    p = momentum(mk)
    x̃, p̃ = per_axis_sliced_lift(x, p, bins, taus; β_mean = 0.7)
    for ℓ in 1:d
        @test (x̃[ℓ, :], p̃[ℓ, :]) ==
              split_quantile_ot_lift(x[ℓ, :], p[ℓ, :], bins[ℓ], taus; β_mean = 0.7)
    end
    @test_throws DimensionMismatch per_axis_sliced_lift(x, p, bins[1:1], taus)
end

@testset "remap of 2D-2V markers" begin
    after = result.markers
    @test size(after.x) == (d, N)
    @test isapprox(sum(after.w), sum(w); rtol = 8 * eps() * N)
    @test isapprox(total_momentum(after), total_momentum(mk);
        atol = 8 * eps() * N * maximum(abs, momentum(mk)))
    @test all(>(0), after.w)
    @test length(info.plan) == d
    # each axis is sliced by its own velocity component
    @test info.plan[1].bins_v == sorted_slices(v[1, :], x[1, :], 3)
    @test info.plan[2].bins_v == sorted_slices(v[2, :], x[2, :], 3)
end

@testset "the frozen per-axis map is symplectic" begin
    q_anc = permutedims(reduce(hcat, [pl.q_anc for pl in info.plan]))
    function frozen(z)
        X = reshape(z[1:(d * N)], d, N)
        P = reshape(z[(d * N + 1):end], d, N)
        out = [begin
                   st = apply_sliced_lift(info.plan[ℓ], X[ℓ, :], P[ℓ, :])
                   canonical_weight_squeeze(st.q̂, st.p̂, info.μ, q_anc[ℓ, :])
               end
               for ℓ in 1:d]
        Q̄ = permutedims(reduce(hcat, first.(out)))
        P̄ = permutedims(reduce(hcat, last.(out)))
        vcat(vec(Q̄), vec(P̄))
    end
    z = vcat(vec(x), vec(momentum(mk)))
    @test isapprox(frozen(z), vcat(vec(info.q̄), vec(info.p̄)); rtol = 8 * eps())
    J = ForwardDiff.jacobian(frozen, z)
    Ω = canonical_J(d * N)
    @test maximum(abs, transpose(J) * Ω * J - Ω) < 8 * eps() * N * cond(J) * opnorm(J)^2
end
