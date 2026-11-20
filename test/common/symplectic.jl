using ParticleRemapping
using LinearAlgebra
using Random
using Test

Random.seed!(12)

@testset "pack and the canonical matrix" begin
    x, p = rand(2, 3), rand(2, 3)
    z = pack(x, p)
    @test unpack(z, 2, 3) == (x, p)
    J = canonical_J(3)
    @test J' == -J
    @test J * J == -I
end

@testset "symplecticity defect" begin
    # a kick p ↦ p + ∇χ(q) with a coupled potential is symplectic
    χgrad(q) = [sin(q[1] - q[2]) + q[1]^2, -sin(q[1] - q[2])]
    kickmap(z) = vcat(z[1:2], z[3:4] .+ χgrad(z[1:2]))
    z = randn(4)
    δ, nrm = symplecticity_defect(kickmap, z)
    @test δ < 8eps() * nrm^2
    # a scaling of one position by 2 is not: DRᵀ J DR − J has the entries ±(2 − 1)
    scaled(z) = [2z[1], z[2], z[3], z[4]]
    @test symplecticity_defect(scaled, z)[1] ≈ 1 rtol=8eps()
    # the same with finite differences
    @test symplecticity_defect(kickmap, z; jacobian = fd_jacobian)[1] < 1e-9
end

@testset "finite-difference Jacobian" begin
    F(z) = [sin(z[1]) * z[2], exp(z[2])]
    z = [0.3, -0.4]
    Jex = [cos(z[1])*z[2] sin(z[1]); 0 exp(z[2])]
    # truncation error h² f'''/6 for order 2 and h⁴ f⁽⁵⁾/30 for order 4
    @test maximum(abs, fd_jacobian(F, z; h = 1e-3) - Jex) < 1e-6
    @test maximum(abs, fd_jacobian(F, z; h = 1e-3) - Jex) > 1e-9
    @test maximum(abs, fd_jacobian(F, z; h = 1e-3, order = 4) - Jex) < 1e-11
    @test_throws r"^ArgumentError: order must be 2 or 4" fd_jacobian(F, z; order = 3)
end

@testset "presymplectic form on (x, v, w)" begin
    m = 1.7
    z = [0.3, 1.2, 0.8, -0.5, -0.7, 1.4]
    Ω = presymplectic_matrix(z; m = m)
    @test Ω' == -Ω
    # the kernel w ∂_w − v ∂_v of each marker
    for a in 0:1
        ξ = zeros(6)
        ξ[3a + 2], ξ[3a + 3] = -z[3a + 2], z[3a + 3]
        @test norm(Ω * ξ) < 8eps()
    end
    @test rank(Ω) == 4
    # a state-dependent gauge step preserves the form; streaming x ↦ x + v does not
    λ(z) = 1 + 0.3 * sin(z[1] + z[5]) + 0.2 * z[3]^2
    gauge(z) = [z[1], λ(z) * z[2], z[3] / λ(z), z[4], λ(z) * z[5], z[6] / λ(z)]
    @test pullback_residual(gauge, z; m = m) < 1e-8
    stream(z) = [z[1] + z[2], z[2], z[3], z[4] + z[5], z[5], z[6]]
    @test pullback_residual(stream, z; m = m) > 0.1
end
