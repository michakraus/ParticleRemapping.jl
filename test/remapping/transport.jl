using ForwardDiff
using LinearAlgebra
using ParticleRemapping
using Random
using Test

Random.seed!(7)

@testset "Simpson's rule" begin
    # exact for cubics, fourth order otherwise
    @test simpson_rule(x -> x^3 - 2x, -1.0, 2.0; n = 2) ≈ 3.75 - 3 rtol=8eps()
    e1 = abs(simpson_rule(exp, 0.0, 1.0; n = 10) - (exp(1) - 1))
    e2 = abs(simpson_rule(exp, 0.0, 1.0; n = 20) - (exp(1) - 1))
    @test e1 / e2 ≈ 16 rtol=0.01
    @test_throws r"^ArgumentError: Simpson's rule needs an even" simpson_rule(exp, 0, 1; n = 3)
end

@testset "(a) the identity behind the matching proposition" begin
    m = 1.3
    ñ(X, P) = (1 + 0.3sin(X)) * exp(-P^2 / 2) * (1 + 0.2P^2)
    W(X, V) = 1 + 0.3cos(X) + 0.2tanh(V)
    dWV(X, V) = W(X, V) + 0.2V * sech(V)^2
    f′ = represented_density(ñ, W, dWV; m)
    ffrozen = represented_density(ñ, W, dWV; m, frozen = true)
    errs = [cumulative_current_error(f′, ñ, W, X, V; m)
            for X in (0.4, 2.1), V in (-2.5, -0.7, 0.9, 3.0)]
    errs_frozen = [cumulative_current_error(ffrozen, ñ, W, X, V; m)
                   for X in (0.4, 2.1), V in (-2.5, -0.7, 0.9, 3.0)]
    # Simpson with 4000 subintervals of width ≤ 1e-3: error h⁴ |g⁗|/180 ≪ 1e-12 (note: 1.8e-14)
    @test maximum(errs) < 1e-12
    # control: the frozen-weight density fails by 5e-2
    @test minimum(errs_frozen) ≈ 0.05 atol=0.005
    # dWV is the derivative of W V
    @test dWV(0.4, 0.9) ≈
          (W(0.4, 0.9 + 1e-6) * (0.9 + 1e-6) - W(0.4, 0.9 - 1e-6) * (0.9 - 1e-6)) / 2e-6 rtol=1e-9
end

@testset "(b) the cotangent lift keeps the current profile iff a²J is constant" begin
    α = 0.4
    J(x) = cosine_profile(x, α)
    a_ok(x) = J(x)^(-1 / 2)
    da_ok(x) = 0.5α * sin(x) * J(x)^(-3 / 2)
    a_bad(x) = 1 + 0.5cos(x)
    da_bad(x) = -0.5sin(x)
    xs = range(0, 2π; length = 9)[1:(end - 1)]
    dev(a, da) = maximum(xs) do x
        φ, ξ = cotangent_lift_flow(a, da, x, 1.5)
        abs(J(φ) * ξ^2 / J(x) - 1)
    end
    # classical Runge–Kutta, dt = 1e-3: global error of order dt⁴ = 1e-12 (note: 6.4e-15)
    @test dev(a_ok, da_ok) < 1e-10
    # control: a = 1 + 0.5 cos x violates it by 5.7
    @test dev(a_bad, da_bad) ≈ 5.7 atol=0.05
    # R = J(x)/(φ'² J(X)) − 1 is the same deviation, inverted
    φ, ξ = cotangent_lift_flow(a_bad, da_bad, 1.0, 1.5)
    @test lift_compatibility_residual(J, 1.0, φ, ξ) ≈ 1 / (J(φ) * ξ^2 / J(1.0)) - 1 rtol=8eps()
end

@testset "cotangent lift: canonical and periodic" begin
    a(x) = 1 + 0.3sin(x)
    da(x) = 0.3cos(x)
    x, p = 2π .* rand(5), randn(5)
    X, P = cotangent_lift(a, da, x, p, 0.7; dt = 1e-3)
    @test all(0 .≤ X .< 2π)
    # the flow of ẋ = a(x) against a fine Runge–Kutta reference and P φ' = p
    for i in eachindex(x)
        φ, ξ = cotangent_lift_flow(a, da, x[i], 0.7; dt = 1e-4)
        @test mod(φ, 2π) ≈ X[i] atol=1e-10
        @test P[i] * ξ ≈ p[i] rtol=1e-10
    end
    # canonical for a fixed a: the map of one marker (x, p) has det D = 1
    lift1(z) = collect(cotangent_lift(a, da, [z[1]], [z[2]], 0.7; L = 1e3, dt = 1e-3))
    @test symplecticity_defect(z -> vcat(lift1(z)...), [0.4, 1.3]; jacobian = fd_jacobian)[1] <
          1e-8
end

@testset "current-quantile matching on a grid" begin
    grid = PeriodicGrid(2π, 4)
    target = SeparableMaxwellianTarget(grid, 0.4)
    # the cell means of 1 + 0.4 cos x
    for c in 1:4
        @test target.Jbar[c] ≈
              simpson_rule(x -> 1 + 0.4cos(x), (c - 1) * π / 2, c * π / 2) / (π / 2) rtol=1e-12
        @test cell_cumulative_current(target, c, Inf) ≈ cell_current_total(target, c, 1) rtol=8eps()
        y = 0.3cell_current_total(target, c, -1)
        @test cell_cumulative_current(target, c, cell_current_inverse(target, c, -1, y)) ≈ y rtol=1e-14
    end
    @test maxwellian_current(maxwellian_current_inverse(0.2)) ≈ 0.2 rtol=1e-14
    @test maxwellian_current(-1.2) == maxwellian_current(1.2)
    X = 2π .* rand(400)
    P = 0.01 .* randn(400)
    V, W, R = current_quantile_match(grid, X, P, target)
    @test W .* V ≈ P rtol=8eps()
    @test all(>(0), W)
    @test all(sign.(V) .== sign.(P))
    @test size(R) == (4, 2)
    for c in 1:4, (j, s) in enumerate((1, -1))

        idx = findall(a -> cellof(grid, X[a]) == c && sign(P[a]) == s, eachindex(X))
        @test R[c, j] ≈ sum(abs.(P[idx])) / grid.Δx / cell_current_total(target, c, s) - 1 rtol=1e-12
        # the order of |P| is the order of |V|
        @test sortperm(abs.(P[idx])) == sortperm(abs.(V[idx]))
    end
    @test_throws r"^ArgumentError: a marker with P = 0" current_quantile_match(
        grid, X, [0.0;
                  P[2:end]], target)
    @test_throws r"^ArgumentError: the momenta \|P\| must be sorted" branch_quantile_matching(
        [0.2, 0.1], 1.0, identity)
    # the method: ψ = id keeps positions and momenta
    mk = Markers(X, P ./ 0.01, fill(0.01, 400))
    r = remap(CurrentQuantileMatching(grid, (x, p) -> (x, p), target), mk)
    @test r.markers.x == X
    @test momentum(r.markers) ≈ momentum(mk) rtol=8eps()
    @test r.info.max_residual == maximum(abs, R)
    # binned weights: the total of the markers inside the bins
    H = binned_weights(grid, X, V, W, range(-4, 4; length = 9))
    @test sum(H) ≈ sum(W[-4 .≤ V .< 4]) rtol=1e-12
    @test binned_weight_error(target, grid, X, V, W) ≥ 0
end

@testset "the matching along rays in two velocity dimensions" begin
    # In d velocity dimensions the gauge step keeps the direction, P = m W V, and the identity of
    # the matching proposition holds along each ray V = r θ:
    #     ∫_0^r f′(r′θ) r′^d dr′ = (1/m) ∫_0^{ρ(r)} ñ(ρ′θ) ρ′^d dρ′,  ρ(r) = m W(rθ) r,
    # with f′ = ñ(m W V) |det ∂(m W V)/∂V| W. The determinant comes from ForwardDiff here, not
    # from the closed form m^d W^{d−1} ∂_r(W r) of the derivation.
    m, d = 1.3, 2
    ñ(P) = exp(-(P[1]^2 + 2P[2]^2) / 2) * (1 + 0.2P[1] + 0.1P[1] * P[2])
    W(V) = 1 + 0.2tanh(V[1]) + 0.1V[2]^2 / (1 + V[2]^2)
    Pof(V) = m * W(V) .* V
    detP(V) = abs(LinearAlgebra.det(ForwardDiff.jacobian(Pof, V)))
    f′(V) = ñ(Pof(V)) * detP(V) * W(V)
    ffrozen(V) = ñ(Pof(V)) * m^d * W(V)^d * W(V)
    function ray_error(f, θ, r)
        lhs = simpson_rule(s -> f(s .* θ) * s^d, 0.0, r)
        ρ = m * W(r .* θ) * r
        rhs = simpson_rule(s -> ñ(s .* θ) * s^d, 0.0, ρ) / m
        abs(lhs - rhs) / abs(rhs)
    end
    θs = [[cos(t), sin(t)] for t in (0.3, 1.9, 4.0)]
    errs = [ray_error(f′, θ, r) for θ in θs, r in (0.7, 2.5)]
    errs_frozen = [ray_error(ffrozen, θ, r) for θ in θs, r in (0.7, 2.5)]
    # Simpson with 4000 subintervals: error far below 1e-12
    @test maximum(errs) < 1e-12
    # control: the frozen-weight determinant m^d W^d misses the radial derivative of W
    @test minimum(errs_frozen) > 1e-3
    # the closed form of the determinant, m^d W^{d−1} ∂_r(W r)
    V, h = [0.8, -0.6], 1e-6
    θ = V ./ hypot(V...)
    dWr = (W(V .+ h .* θ) * (1 + h / hypot(V...)) -
           W(V .- h .* θ) * (1 - h / hypot(V...))) /
          (2h / hypot(V...))
    @test detP(V) ≈ m^d * W(V)^(d - 1) * dWr rtol = 1e-8
end

@testset "DensityTarget against the separable target" begin
    grid = PeriodicGrid(2π, 8)
    sep = SeparableMaxwellianTarget(grid, 0.4)
    # the midpoint rule on 64 points per cell and the trapezoidal rule with h = 0.002 in v are
    # accurate to about 1e-6 here; the tolerances sit above that and below a wrong branch or cell
    tab = DensityTarget(grid, (x, v) -> (1 + 0.4cos(x)) * exp(-v^2 / 2) / sqrt(2π);
        vmax = 8.0, nx = 64, nv = 4000)
    for c in (1, 4, 8), s in (1, -1)

        @test isapprox(cell_current_total(tab, c, s), cell_current_total(sep, c, s);
            rtol = 1e-4)
        y = 0.6 * cell_current_total(sep, c, s)
        @test isapprox(
            cell_current_inverse(tab, c, s, y), cell_current_inverse(sep, c, s, y);
            rtol = 1e-3)
    end
    @test isapprox(
        cell_cumulative_current(tab, 3, -1.2), cell_cumulative_current(sep, 3, 1.2);
        rtol = 1e-4)
    # an asymmetric density: the two branches differ
    skew = DensityTarget(grid, (x, v) -> v > 0 ? 2exp(-v^2 / 2) : exp(-v^2 / 2); vmax = 8.0)
    @test isapprox(cell_current_total(skew, 2, 1), 2cell_current_total(skew, 2, -1);
        rtol = 1e-6)
end
