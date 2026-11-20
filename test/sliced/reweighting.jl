using ParticleRemapping
using Random: MersenneTwister
using Test

rng = MersenneTwister(21)
N = 30
x, v = 2.0 .* rand(rng, N) .- 1.0, randn(rng, N)
w = 0.5 .+ rand(rng, N)

@testset "kernel density ratio" begin
    # one source and one evaluation point: the ratio is w K(Δx, Δv) with the periodic kernel
    for (Lx, hx, hv) in Iterators.product((1.3, 4.0), (0.3, 0.6), (0.25, 0.5))
        s = Lx / π
        Δx, Δv = 0.7, -0.4
        K = exp(-(2s^2 * (1 - cos(Δx / s)) / hx^2 + Δv^2 / hv^2) / 2)
        r = kernel_density_ratio([0.2 + Δx], [0.1 + Δv], [0.2], [0.1], [1.7]; hx, hv, Lx)
        @test isapprox(r[1], 1.7 * K; rtol = 8 * eps())
    end
    # evaluation points equal to equal-weight sources: the ratio is the weight
    r = kernel_density_ratio(x, v, x, v, fill(0.3, N); Lx = 1.0)
    @test isapprox(r, fill(0.3, N); rtol = 8 * eps() * N)
    # periodic in x with period 2Lx
    r1 = kernel_density_ratio(x, v, x .+ 0.1, v, w; Lx = 1.0)
    r2 = kernel_density_ratio(x .+ 2.0, v, x .+ 0.1, v, w; Lx = 1.0)
    @test isapprox(r1, r2; rtol = 8 * eps() * N)
    # a 1 × N matrix gives the vector result; in 2D the kernel is the product over the axes
    @test kernel_density_ratio(reshape(x, 1, :), reshape(v, 1, :), reshape(x, 1, :),
        reshape(v, 1, :), w) == kernel_density_ratio(x, v, x, v, w)
    X, V = [0.3; -0.2;;], [0.1; 0.4;;]
    r2d = kernel_density_ratio(X .+ [0.5, 0.2], V .+ [-0.3, 0.1], X, V, [2.0])
    k1 = kernel_density_ratio([0.8], [-0.2], [0.3], [0.1], [1.0])[1]
    k2 = kernel_density_ratio([0.0], [0.5], [-0.2], [0.4], [1.0])[1]
    @test isapprox(r2d[1], 2.0 * k1 * k2; rtol = 8 * eps())
end

@testset "squeeze multipliers" begin
    m, w_ref = 1.3, sum(w) / N
    q̂ = x .+ 0.05 .* randn(rng, N)
    p̂ = m .* w_ref .* v
    q_anc = round.(q̂; digits = 1)
    for iterations in (1, 3)
        μ = squeeze_multipliers(q̂, p̂, q_anc, x, v, w; m, iterations, Lx = 1.0)
        @test isapprox(sum(μ), N; rtol = 8 * eps() * N)
        @test all(>(0), μ)
    end
    # one pass: μ ∝ r(q̂, p̂ / (m w_ref)) + ratio_floor w_ref
    for ratio_floor in (0.0, 0.08)
        μ1 = squeeze_multipliers(q̂, p̂, q_anc, x, v, w; m, iterations = 1, Lx = 1.0,
            ratio_floor)
        raw = kernel_density_ratio(q̂, p̂ ./ (m * w_ref), x, v, w; Lx = 1.0) .+
              ratio_floor * w_ref
        @test isapprox(μ1, N .* raw ./ sum(raw); rtol = 8 * eps() * N)
    end
    # the corrector passes change the multipliers
    @test squeeze_multipliers(q̂, p̂, q_anc, x, v, w; m, iterations = 3, Lx = 1.0) !=
          squeeze_multipliers(q̂, p̂, q_anc, x, v, w; m, iterations = 1, Lx = 1.0)
end
