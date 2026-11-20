using ParticleRemapping
using ForwardDiff: ForwardDiff
using LinearAlgebra: cond, opnorm
using Test

ex = sliced_stages_example()
before, result = ex.before, ex.result
info = result.info
plan = info.plan
N = length(before)
m = before.m
w_ref = sum(before.w) / N

@testset "method defaults" begin
    method = SlicedQuantileRemap()
    @test (method.Lx, method.Vmax, method.B_v, method.B_x) == (Float64(π), 2.25, 15, 90)
    @test (method.β_mean_q, method.β_intra_q, method.β_mean_p, method.β_intra_p) ==
          (0.75, 0.78, 0.78, 0.80)
    @test (method.ε_q_spacings, method.ε_p_spacings, method.hx, method.hv) ==
          (2.0, 2.0, 0.45, 0.35)
    @test (method.ratio_floor, method.iterations) == (0.08, 3)
end

@testset "plan" begin
    @test length(plan.bins_v) == 5 && all(length.(plan.bins_v) .== 24)
    @test length(plan.bins_x) == 12 && all(length.(plan.bins_x) .== 10)
    # velocity slices ordered by position, position slices ordered by the lifted momentum
    @test all(issorted(before.x[idx]) for idx in plan.bins_v)
    @test all(maximum(before.v[plan.bins_v[b]]) <= minimum(before.v[plan.bins_v[b + 1]])
    for b in 1:4)
    @test all(issorted(info.p̃[idx]) for idx in plan.bins_x)
    @test all(maximum(info.q̃[plan.bins_x[b]]) <= minimum(info.q̃[plan.bins_x[b + 1]])
    for b in 1:11)
    # anchors at the slice centres, targets at the uniform quantiles
    vb = uniform_quantile_targets(2.25, 5)
    @test all(plan.p_anc[a] == m * w_ref * vb[b] for b in 1:5 for a in plan.bins_v[b])
    xb = uniform_quantile_targets(Float64(π), 12)
    @test all(plan.q_anc[a] == xb[b] for b in 1:12 for a in plan.bins_x[b])
    @test all(τ == uniform_quantile_targets(Float64(π), 24) for τ in plan.taus_x)
    vj = uniform_quantile_targets(2.25, 10)
    @test all(τ == m * w_ref .* vj for τ in plan.taus_p)
    @test all(plan.slot_v[idx] == vj for idx in plan.bins_x)
    # mollifier widths: two target spacings
    @test plan.ε_q == 2.0 * (2π / 24)
    @test plan.ε_p == 2.0 * (2 * maximum(abs, m * w_ref .* vj) / 10)
end

@testset "conservation and the final stage" begin
    after = result.markers
    @test isapprox(sum(after.w), sum(before.w); rtol = 8 * eps() * N)
    @test isapprox(sum(momentum(after)), sum(momentum(before));
        atol = 8 * eps() * N * maximum(abs, momentum(before)))
    @test isapprox(sum(info.p̃), sum(momentum(before)); atol = 8 * eps() * N^2 * w_ref)
    @test all(>(0), after.w)
    @test after.w == w_ref .* info.μ
    @test after.x == info.q̄
    @test isapprox(momentum(after), info.p̄; rtol = 8 * eps())
    q̄, p̄ = canonical_weight_squeeze(info.q̂, info.p̂, info.μ, plan.q_anc)
    @test q̄ == info.q̄ && p̄ == info.p̄
end

@testset "the frozen map is symplectic" begin
    function frozen(z; with_cotangent = true)
        st = apply_sliced_lift(plan, z[1:N], z[(N + 1):(2N)]; with_cotangent)
        vcat(canonical_weight_squeeze(st.q̂, st.p̂, info.μ, plan.q_anc)...)
    end
    z = vcat(before.x, momentum(before))
    @test isapprox(frozen(z), vcat(info.q̄, info.p̄); rtol = 8 * eps())
    J = ForwardDiff.jacobian(frozen, z)
    Ω = canonical_J(N)
    @test maximum(abs, transpose(J) * Ω * J - Ω) < 8 * eps() * N * cond(J) * opnorm(J)^2
    # control: without the cotangent parts
    Jc = ForwardDiff.jacobian(x -> frozen(x; with_cotangent = false), z)
    @test maximum(abs, transpose(Jc) * Ω * Jc - Ω) > 1e-1
end
