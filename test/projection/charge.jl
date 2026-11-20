# The charge-conserving reweighting on a random ensemble in d = 1 after a kick–drift map: the grid
# charge, positivity, the total weight and the momentum; the minimal relative entropy, checked on
# the null space of the charge matrix independently of the dual problem; the identity. Cases that
# must fail: the target weights alone, the additive least-norm correction, and positions after
# which a shape function with charge has no marker.

using ParticleRemapping
using LinearAlgebra
using Printf
using Random
using Test

const L, m = 2π, 1.0
const N = 400
const grid = UniformBasis(0.0, L / 16, 16, true)

Random.seed!(20261121)

x = L .* rand(1, N)
v = randn(1, N)
w = exp.(0.5 .* randn(N))
w .*= L / sum(w)
p = m .* w' .* v

# a kick–drift map
gκ(ξ, w) = [0.4 * w * cos(ξ[1]) + 0.2 * sin(2ξ[1])]
gη(v, w) = [0.3 * w * v[1] + 0.2 * sin(v[1])]
ps = centroid_kick(x, p, w, gκ)
xs = vec(velocity_drift(x, ps, w, m, gη))

ρ0 = spline_grid_charge(grid, vec(x), w)
wstar = w .* exp.(0.3 .* randn(N))
wstar .*= sum(w) / sum(wstar)

relative_entropy(u, ustar) = sum(u .* log.(u ./ ustar) .- u .+ ustar)

@testset "Charge-conserving reweighting" begin
    w1, λ, its, ok = charge_reweight(grid, xs, wstar, ρ0)
    A = spline_charge_matrix(grid, xs)
    res = maximum(abs, A * w1 .- ρ0) / maximum(ρ0)
    res_star = maximum(abs, A * wstar .- ρ0) / maximum(ρ0)
    @printf("  Newton: %d iterations; charge residual %.1e (target weights alone %.2f)\n",
        its,
        res, res_star)

    @testset "grid charge, positivity, total weight, momentum" begin
        @test ok
        @test its == 5                             # the note quotes 5 Newton iterations
        @test res < 1e-13                          # the note quotes 7.1e-16
        @test A * w1 ≈ spline_grid_charge(grid, xs, w1) rtol = 8eps()
        @test abs(sum(w1) - sum(w)) < 1e-13 * sum(w)
        @test minimum(w1) > 0
        @test abs(sum(ps) - sum(p)) < 1e-13 * sum(abs, p)
        @test w1 ≈ wstar .* exp.(A' * λ) rtol = 8eps()
    end

    @testset "minimal relative entropy" begin
        Z = nullspace(A)
        g = log.(w1 ./ wstar)                      # gradient of the relative entropy at w1
        @printf("  stationarity on the null space of A: %.1e\n", norm(Z' * g) / norm(g))
        @test norm(Z' * g) / norm(g) < 1e-10       # the note quotes 9.1e-16
        E0 = relative_entropy(w1, wstar)
        worse = minimum(1:50) do _
            δ = Z * randn(size(Z, 2))
            δ .*= 0.05 * minimum(w1) / maximum(abs, δ)
            min(relative_entropy(w1 .+ δ, wstar), relative_entropy(w1 .- δ, wstar)) - E0
        end
        @test worse > 0
    end

    @testset "identity" begin
        w2, λ2, its2, ok2 = charge_reweight(grid, vec(x), w, ρ0)
        @test ok2 && its2 == 0 && w2 == w
    end

    @testset "must fail: target weights alone" begin
        @test res_star > 0.28                      # the note quotes 0.29
    end

    @testset "must fail: additive correction" begin
        # an excess of charge on [0, 1) where half of the target weights are small
        wbad = copy(w)
        idx = findall(k -> xs[k] < 1, eachindex(xs))
        wbad[idx[1:2:end]] .*= 3
        wbad[idx[2:2:end]] .*= 1e-3
        wadd = additive_charge_correction(grid, xs, wbad, ρ0)
        wexp, _, _, okexp = charge_reweight(grid, xs, wbad, ρ0)
        @printf("  additive: min w = %.3f; exponential: min w = %.1e\n", minimum(wadd),
            minimum(wexp))
        @test maximum(abs, A * wadd .- ρ0) < 1e-12 * maximum(ρ0)
        @test minimum(wadd) < -0.01                # the note quotes −0.015
        @test okexp && minimum(wexp) > 0
    end

    @testset "must fail: a shape function without markers" begin
        xhalf = mod.(xs, L) ./ 2                   # every marker in [0, L/2)
        @test !charge_reweight(grid, xhalf, wstar, ρ0).converged
        mk = Markers(xhalf, ps[1, :] ./ (m .* wstar), wstar; m)
        @test_throws r"no positive weights reproduce" remap(
            ChargeConservingReweighting(grid), mk;
            reference_charge = ρ0)
    end

    @testset "the method type" begin
        mk = Markers(xs, ps[1, :] ./ (m .* wstar), wstar; m)
        r = remap(ChargeConservingReweighting(L, 16), mk; reference_charge = ρ0)
        @test r.markers.w ≈ w1 rtol = 1e-12
        @test r.markers.x == xs
        @test momentum(r.markers) ≈ ps[1, :] rtol = 8eps()
        @test r.info.converged && r.info.iterations == its
        @test r.info.charge_residual < 1e-13
        @test r.info.charge_residual_target ≈ res_star rtol = 1e-12
        @test r.info.max_log_change ≈ maximum(abs, log.(w1 ./ wstar)) rtol = 1e-12
    end
end
