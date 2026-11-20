using ParticleRemapping
using Random: MersenneTwister
using Test

@testset "uniform quantile targets" begin
    τ = uniform_quantile_targets(2.0, 4)
    @test τ == [-1.5, -0.5, 0.5, 1.5]
    # midpoints of n equal cells: symmetric about 0, spacing 2L/n, half a spacing from the ends
    τ = uniform_quantile_targets(π, 7)
    @test isapprox(sum(τ), 0; atol = 8 * eps() * 7π)
    @test all(isapprox.(diff(τ), 2π / 7; atol = 8 * eps() * π))
    @test isapprox(τ[1] + π, π / 7; rtol = 8 * eps())
end

@testset "sorted slices and anchors" begin
    rng = MersenneTwister(1)
    N, B = 24, 4
    key, order = randn(rng, N), randn(rng, N)
    bins = sorted_slices(key, order, B)
    @test length(bins) == B
    @test all(length.(bins) .== N ÷ B)
    @test sort(reduce(vcat, bins)) == 1:N
    # every key of slice b is below every key of slice b + 1
    @test all(maximum(key[bins[b]]) < minimum(key[bins[b + 1]]) for b in 1:(B - 1))
    # each slice is ordered by the second key
    @test all(issorted(order[idx]) for idx in bins)
    @test_throws r"^ArgumentError: 24 markers do not split into 5 equal slices" sorted_slices(
        key, order, 5)

    centres = uniform_quantile_targets(1.0, B)
    anc = slice_anchors(bins, centres)
    @test all(anc[a] == centres[b] for b in 1:B for a in bins[b])
end
