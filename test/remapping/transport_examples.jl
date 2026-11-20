using ParticleRemapping
using Test

# A number quoted with d decimals in the note agrees to half a unit of its last digit.
agrees(x, quoted, d) = abs(x - quoted) ≤ 0.5 * 10.0^(-d) + 8eps(abs(quoted))

@testset "Example 3: the admissible cotangent lift in closed form" begin
    α, t = 0.4, 1.5
    J(x) = cosine_profile(x, α)
    xs = [0.0, π / 2, π, 4.0]
    # Θ' = J^{1/2}, against the central difference (truncation h² Θ'''/6 ≈ 1e-11)
    h = 1e-5
    for x in (0.7, 2.9)
        @test (lift_theta(x + h) - lift_theta(x - h)) / 2h ≈ sqrt(J(x)) atol=1e-9
        @test lift_theta_inverse(lift_theta(x)) ≈ x atol=1e-13
    end
    maps = [admissible_lift(x, t) for x in xs]
    # (3a) J(φ_t) φ_t'² = J
    @test maximum(abs(J(X) * dφ^2 / J(x) - 1) for (x, (X, dφ)) in zip(xs, maps)) < 8eps()
    # the closed form is the flow of ẋ = J^{−1/2}: against the Runge–Kutta flow (dt⁴ ≈ 1e-12)
    a(x) = J(x)^(-1 / 2)
    da(x) = 0.5α * sin(x) * J(x)^(-3 / 2)
    for (x, (X, dφ)) in zip(xs, maps)
        φ, ξ = cotangent_lift_flow(a, da, x, t)
        @test φ ≈ X atol=1e-9
        @test ξ ≈ dφ atol=1e-9
    end
    # the table of the note at t = 1.5
    @test all(agrees.(first.(maps), [1.3197, 3.3328, 4.8586, 5.4954], 4))
    @test all(agrees.(last.(maps), [1.1285, 1.2832, 0.7530, 0.7590], 4))
    @test all(agrees.(1 ./ last.(maps), [0.8862, 0.7793, 1.3281, 1.3176], 4))
    # (3a) the matched velocity is v and the weight w0/φ'
    w0 = 0.01
    for x in xs, v in (-1.3, 0.4, 2.2)

        r = admissible_lift_matching(x, v, t, w0)
        @test r.V ≈ v atol=1e-14
        @test r.W ≈ w0 / r.dφ rtol=1e-13
        @test r.P ≈ r.W * r.V rtol=8eps()
    end
    # (3b) at half a period of Θ the marker density (J(x) J(X))^{1/2} has max/min 1.18
    th = lift_theta(2π) / 2
    @test agrees(th, 3.1089, 4)
    xg = range(0, 2π; length = 401)
    dens = [sqrt(J(x) * J(admissible_lift(x, th)[1])) for x in xg]
    @test maximum(dens) / minimum(dens) < 1.2
    @test agrees(maximum(dens) / minimum(dens), 1.18, 2)
    # control: at t = 0 it is J, with max/min 7/3
    @test maximum(J.(xg)) / minimum(J.(xg)) ≈ 7 / 3 rtol=8eps()
    # control: the lift of a = 1 + 0.5 cos x is incompatible, R = 3.70, 4.03, −0.53, −0.85
    R = map(xs) do x
        φ, ξ = cotangent_lift_flow(x -> 1 + 0.5cos(x), x -> -0.5sin(x), x, t)
        lift_compatibility_residual(J, x, φ, ξ)
    end
    @test all(agrees.(R, [3.70, 4.03, -0.53, -0.85], 2))
end

@testset "Example 4: five markers in one cell" begin
    e = quantile_matching_example()
    @test e.aP ≈ [0.02, 0.05, 0.075, 0.11, 0.145] rtol=8eps()
    @test e.Π ≈ [0.01, 0.045, 0.1075, 0.2, 0.3275] rtol=8eps()
    @test sum(e.aP) / (e.m * e.Δx) ≈ 0.4 rtol=8eps()
    @test agrees(e.Jbar * maxwellian_current(Inf), 0.39894, 5)
    @test agrees(e.R, 0.00265, 5)
    @test all(agrees.(e.Πc, [0.00997, 0.04488, 0.10722, 0.19947, 0.32663], 5))
    @test agrees(sqrt(-2log(1 - sqrt(2π) * e.Πc[5])), 1.8482, 4)
    @test all(agrees.(e.V, [0.2250, 0.4886, 0.7912, 1.1774, 1.8482], 4))
    @test all(agrees.(e.W, [0.0889, 0.1023, 0.0948, 0.0934, 0.0785], 4))
    # (4a) m W V = P and the markers sit at the closed current quantiles
    @test e.m .* e.W .* e.V ≈ e.aP rtol=8eps()
    @test maximum(abs.(e.Jbar .* maxwellian_current.(e.V) .- e.Πc)) < 8eps(maximum(e.Πc))
    # control: the old velocities miss them by up to 0.067
    @test agrees(maximum(abs.(e.Jbar .* maxwellian_current.(e.v) .- e.Πc)), 0.067, 3)
    # charge and kinetic energy only to quantile accuracy
    @test agrees(sum(e.W), 0.458, 3)
    @test agrees(sum(e.W .* e.V .^ 2), 0.486, 3)
    @test sum(fill(e.w0, 5)) ≈ 0.5 rtol=8eps()
    @test agrees(sum(e.w0 .* e.v .^ 2), 0.417, 3)
    @test_throws r"^ArgumentError: the example lies on the positive branch" quantile_matching_example(;
        v = [-0.2, 0.5])
end
