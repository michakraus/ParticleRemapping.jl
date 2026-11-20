using ParticleRemapping
using Test

# Central differences with h = 1e-5: truncation error h²/6 times third derivatives of order 10,
# round-off eps/h times entries of order 10, so a symplectic map leaves residuals below 1e-8.
const FD_TOL = 1e-8

@testset "symplecticity condition, the seven checks" begin
    m = 1.7
    samples = [[0.3, -1.2, 0.8], [1.1, 0.4, 2.5], [-0.7, 2.0, 0.15], [2.2, -0.3, 1.3]]
    checks = symplectic_condition_checks(; m, samples)
    @test length(checks) == 7
    r = Dict(c.name => c.residual for c in checks)
    for c in checks
        c.holds && @test c.residual < FD_TOL
    end
    @test count(c -> c.holds, checks) == 4

    vmax = maximum(abs(z[2]) for z in samples)
    # control A: X = x + v gives DΦᵀΩDΦ − Ω with the entries ±m v
    @test r["control: X = x + v"] ≈ m * vmax atol=FD_TOL
    # free streaming: Φ_t*ω − ω = m t v dv ∧ dw with t = 0.9
    @test r["free streaming on (x, v, w)"] ≈ m * 0.9 * vmax atol=FD_TOL
    # control B: det ∂(X, P)/∂(x, v) = m w q(v) with q the derivative of (v + 0.4v³)/(1 + 1.2v²)
    # times (1 + 1.2v²), so the residual is m w |q(v) − 1|
    q(v) = ((1 + 1.2v^2)^2 - 2.4v * (v + 0.4v^3)) / (1 + 1.2v^2)^2
    @test r["control: V = v + 0.4 v³, W = w/det DT"] ≈
          maximum(m * z[3] * abs(q(z[2]) - 1) for z in samples) atol=FD_TOL
    # the margins the note quotes: 3.4 and 1.2
    @test r["control: X = x + v"] ≈ 3.4 atol=0.05
    @test r["control: V = v + 0.4 v³, W = w/det DT"] ≈ 1.2 atol=0.05
end

@testset "frozen-weight residual" begin
    # a drift x ↦ x + v is canonical for a frozen weight, a stretch v ↦ 2v is not (det = 2 m w)
    @test frozen_weight_residual((x, v) -> (x + v, v), (x, v, w) -> w, 0.3, 0.7, 1.5;
        m = 2) < FD_TOL
    @test frozen_weight_residual((x, v) -> (x, 2v), (x, v, w) -> w, 0.3, 0.7, 1.5; m = 2) ≈
          2 * 1.5 atol=FD_TOL
end
