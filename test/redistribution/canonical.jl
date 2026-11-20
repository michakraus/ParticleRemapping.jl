# R1: the pull-back of ω to (x, v, w), its kernel, the factorisation freedom, and the
# velocity–momentum fold.

using ParticleRemapping
using ForwardDiff
using LinearAlgebra
using Test

@testset "Pull-back of ω to (x, v, w)" begin
    m = 1.7
    # two markers, (x₁, v₁, w₁, x₂, v₂, w₂)
    y0 = [0.4, -0.8, 0.25, 1.9, 1.3, 0.40]
    Ω = presymplectic_matrix(y0; m = m)

    # automatic differentiation: exact up to round-off
    Ω_ad = kinematic_pullback(y0; m = m)
    @test norm(Ω_ad - Ω) <= 8 * eps() * norm(Ω)

    # fourth-order central differences: π is quadratic, so the truncation error vanishes and
    # the round-off of the step, eps/h, sets the error
    h = 1e-5
    Ω_fd = kinematic_pullback(y0; m = m, jacobian = (F, y) -> fd_jacobian(F, y; h = h, order = 4))
    @test norm(Ω_fd - Ω) < 100 * eps() / h

    # the kernel: N Euler scaling fields, rank 2N
    K = euler_scaling_fields(y0)
    @test norm(Ω * K) <= 8 * eps() * norm(Ω) * norm(K)
    @test rank(Ω) == 4
    @test rank(K) == 2
end

@testset "Factorisation freedom: a canonical map and any positive split" begin
    m = 1.3
    y0 = [0.4, -0.8, 0.25, 1.9, 1.3, 0.40, 3.1, 0.2, 0.33]
    # a symplectic map of (q, p) (a kick and a drift), then a state-dependent weight
    function canonical_then_split(y; keep_momentum = true)
        x, v, w = y[1:3:end], y[2:3:end], y[3:3:end]
        p = m .* w .* v
        p̄ = p .- 0.3 .* sin.(x)
        x̄ = x .+ 0.5 .* p̄
        w̄ = 1 .+ 0.2 .* cos.(x) .* w .^ 2 .+ 0.1 .* sum(abs2, v)
        v̄ = keep_momentum ? p̄ ./ (m .* w̄) : p̄ ./ (m .* w)
        collect(Iterators.flatten(zip(x̄, v̄, w̄)))
    end
    @test pullback_residual(canonical_then_split, y0; m = m, jacobian = ForwardDiff.jacobian) <
          64 * eps()
    # control: the same weights with the velocity kept from the old weights change p̄
    @test pullback_residual(y -> canonical_then_split(y; keep_momentum = false), y0; m = m,
        jacobian = ForwardDiff.jacobian) > 1e-2
end

@testset "Velocity–momentum fold" begin
    vth = 1.3
    f(v) = isotropic_maxwellian(v; vth = vth)
    # 1D: the determinant f + v f' = (1 − v²/v_th²) f vanishes on the thermal sphere and is
    # negative beyond it
    @test abs(velocity_fold_determinant(f, vth)) < 8 * eps() * f(vth)
    @test velocity_fold_determinant(f, 0.5vth) > 0
    @test velocity_fold_determinant(f, 1.5vth) < 0
    for v in (-2.2, -0.4, 0.3, 1.1, 2.5)
        @test velocity_fold_determinant(f, v) ≈ (1 - v^2 / vth^2) * f(v) rtol = 8 * eps()
    end

    # the bound p_max = m ΔΩ v_th C e^{−1/2}, reached on |v| = v_th
    m, ΔΩ = 1.5, 0.2
    pmax = maxwellian_fold_momentum_bound(; vth = vth, m = m, ΔΩ = ΔΩ)
    vs = range(-6, 6; length = 2001)
    @test maximum(abs.(velocity_fold_momentum.(f, vs; m = m, ΔΩ = ΔΩ))) <= pmax
    @test abs(velocity_fold_momentum(f, vth; m = m, ΔΩ = ΔΩ)) ≈ pmax rtol = 8 * eps()
    # both branches |v| < v_th and |v| > v_th reach the same momentum
    P = v -> velocity_fold_momentum(f, v; m = m, ΔΩ = ΔΩ)
    target = P(0.5vth)
    lo, hi = vth, 10vth                      # P decreases on (v_th, ∞): bisect for the root
    for _ in 1:200
        mid = (lo + hi) / 2
        P(mid) > target ? (lo = mid) : (hi = mid)
    end
    @test lo > 1.5vth
    @test P(lo) ≈ target rtol = 1e-12

    # 2D: the determinant identity (Sylvester) against the Jacobian of P, anisotropic density
    g(v) = exp(-v[1]^2 / 2 - v[2]^2 / 3 + 0.3 * v[1] * v[2] / 4)
    for v in ([0.3, -0.7], [1.4, 0.9], [-2.0, 1.5])
        DP = ForwardDiff.jacobian(u -> velocity_fold_momentum(g, u; m = m, ΔΩ = ΔΩ), v)
        @test det(DP) ≈ velocity_fold_determinant(g, v; m = m, ΔΩ = ΔΩ) rtol = 1e-13
    end
    # 2D Maxwellian: zero on the thermal circle, negative outside
    h2(v) = isotropic_maxwellian(v; vth = vth)
    @test abs(velocity_fold_determinant(h2, vth .* [cos(0.7), sin(0.7)])) < 8 * eps()
    @test velocity_fold_determinant(h2, 1.2vth .* [cos(0.7), sin(0.7)]) < 0
end
