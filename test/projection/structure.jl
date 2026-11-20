# The structural facts of the Hamiltonian particle projection: the compressibility of a
# weight-dependent Hamiltonian, the fold of the implicit weight, and one map for all markers under a
# mean-field Hamiltonian, with the cases that must fail.

using ParticleRemapping
using ForwardDiff
using LinearAlgebra
using Printf
using Random
using Test

# ---- compressibility {W, ∂_wΨ} -----------------------------------------------------------------

# A random smooth Ψ(x, p, w) as a sum of modes that couple x, p and ln w.
function random_Ψ(d; modes = 5)
    K = [randn(d) for _ in 1:modes]
    L = [randn(d) for _ in 1:modes]
    c = randn(modes)
    a = randn(modes)
    φ = 2π * rand(modes)
    (x, p, w) -> sum(a[j] * sin(dot(K[j], x) + dot(L[j], p) + c[j] * log(w) + φ[j])
    for j in 1:modes) + 0.3 * sum(abs2, p) * log(w)
end

# A random positive weight field W(x, p).
function random_W(d)
    q = randn(d)
    r = randn(d)
    (x, p) -> 1.5 + 0.5 * sin(dot(q, x) + dot(r, p)) + 0.2 * cos(dot(r, x))
end

function compressibility_case(d)
    Ψ = random_Ψ(d)
    W = random_W(d)
    errs, errs_sign, errs_frozen = Float64[], Float64[], Float64[]
    Wζ(ζ) = W(ζ[1:d], ζ[(d + 1):(2d)])
    for _ in 1:20
        ζ0 = randn(2d)
        w0 = Wζ(ζ0)
        Ψw(ζ) = ForwardDiff.derivative(w -> Ψ(ζ[1:d], ζ[(d + 1):(2d)], w), w0)
        lhs = effective_divergence(Ψ, W, ζ0, d)
        rhs = compressibility_bracket(Ψ, W, ζ0, d)
        rhs_total = compressibility_bracket(Ψ, W, ζ0, d; total = true)
        push!(errs, max(abs(lhs - rhs), abs(lhs - rhs_total)) / max(1, abs(lhs)))
        # must fail: the opposite sign {∂_wΨ, W}, and the divergence at frozen weight
        push!(errs_sign, abs(lhs - canonical_bracket(Ψw, Wζ, ζ0, d)) / max(1, abs(lhs)))
        push!(errs_frozen, abs(lhs - frozen_weight_divergence(Ψ, W, ζ0, d)) /
                           max(1, abs(lhs)))
    end
    maximum(errs), minimum(errs_sign), minimum(errs_frozen)
end

# ---- the fold of the implicit weight -----------------------------------------------------------

# A random positive density in v: a mixture of three Gaussians.
function random_f(d)
    μ = [0.5 * randn(d) for _ in 1:3]
    σ = 0.7 .+ rand(3)
    a = 0.5 .+ rand(3)
    v -> sum(a[j] * exp(-sum(abs2, v - μ[j]) / (2σ[j]^2)) for j in 1:3)
end

function fold_case(d; m = 1.3, c = 0.7)
    f = random_f(d)
    errs, errs_control = Float64[], Float64[]
    for _ in 1:20
        v = randn(d)
        detAD = det(ForwardDiff.jacobian(u -> implicit_weight_map(f, u; m, c), v))
        formula = fold_determinant(f, v; m, c)
        # must fail for d ≥ 2: f^d in place of f^(d−1)
        control = formula * f(v)
        push!(errs, abs(detAD - formula) / max(abs(detAD), 1e-300))
        push!(errs_control, abs(detAD - control) / max(abs(detAD), 1e-300))
    end
    maximum(errs), minimum(errs_control)
end

Random.seed!(20261119)

@testset "Compressibility and fold" begin
    # Each testset starts from the same state of the random number generator. The note quotes
    # |div u − {W, ∂_wΨ}| ≤ 1.1e-15 (round-off of second derivatives of O(1) functions), the
    # controls ≥ 7.8e-2 and ≥ 3.9e-2.
    @testset "compressibility, d = $d" for d in 1:3
        err, err_sign, err_frozen = compressibility_case(d)
        @printf("  d = %d: |div u − {W, ∂_wΨ}| ≤ %.1e; wrong sign ≥ %.1e, frozen W ≥ %.1e\n",
            d, err, err_sign, err_frozen)
        @test err < 16eps()
        @test err_sign > 7.8e-2        # must fail
        @test err_frozen > 3.9e-2      # must fail
    end

    # The note quotes a relative error ≤ 1.7e-15 and the control f^d ≥ 2.4e-3.
    @testset "fold determinant, d = $d" for d in 1:3
        err, err_control = fold_case(d)
        @printf("  d = %d: fold determinant relative error ≤ %.1e; control f^d ≥ %.1e\n",
            d,
            err, err_control)
        @test err < 16eps()
        d >= 2 && @test err_control > 2.4e-3   # must fail
    end

    @testset "fold of a Maxwellian" begin
        vth = 1.0
        fv(v) = exp(-sum(abs2, v) / (2vth^2)) / (2π * vth^2)
        detΦ(v) = det(ForwardDiff.jacobian(u -> implicit_weight_map(fv, u), v))
        e = normalize(randn(2))
        @test detΦ(0.99vth * e) > 0
        @test detΦ(1.01vth * e) < 0
        @test abs(detΦ(vth * e)) < 1e-14
        # p_max = v_th f_M(v_th) = e^{-1/2}/√(2π) for m = c = v_th = 1
        @test maxwellian_fold_peak() ≈ exp(-0.5) / sqrt(2π) rtol = 8eps()
        pmax = maxwellian_fold_peak()
        @test maxwellian_fold_roots(0.5pmax) == 2
        @test maxwellian_fold_roots(1.5pmax) == 0
    end
end

# ---- one map for all markers under a mean-field Hamiltonian ------------------------------------

Random.seed!(20261119)
ex = mean_field_example()
H = mean_field_hamiltonian(ex.h, ex.φ)

@testset "Mean-field Hamiltonian: one map for all markers" begin
    @test length(ex.Z0) == 274

    @testset "∇_{z_k} H = ∇Ψ_μ(z_k)" begin
        G = ensemble_gradient(H, ex.Z0)
        Ψ = mean_field_single_hamiltonian(ex.h, ex.φ, ex.Z0)
        err = maximum(norm(G[k] - ForwardDiff.gradient(Ψ, ex.Z0[k]))
        for k in eachindex(ex.Z0))
        @printf("  %d markers: max |∇_{z_k}H − ∇Ψ(z_k)| = %.1e\n", length(ex.Z0), err)
        @test err < 16eps()
    end

    for (name, ensemble, single) in (("symplectic Euler", ensemble_symplectic_euler,
        particle_symplectic_euler),
        ("implicit midpoint", ensemble_implicit_midpoint, particle_implicit_midpoint))
        @testset "one map, $name" begin
            Z1, Zhat = ensemble(H, ex.Z0, ex.τ)
            Ψhat = mean_field_single_hamiltonian(ex.h, ex.φ, Zhat)
            ψ(ζ) = single(Ψhat, ζ, ex.τ)
            # the note quotes |ψ(z_k) − z'_k| ≤ 6.8e-16 while markers move up to 0.39
            err = maximum(norm(ψ(ex.Z0[k]) - Z1[k]) for k in eachindex(ex.Z0))
            moved = maximum(norm(Z1[k] - ex.Z0[k]) for k in eachindex(ex.Z0))
            dets = [det(ForwardDiff.jacobian(ψ, randn(2))) for _ in 1:20]
            @test err < 16eps()
            @test moved > 0.37
            @test maximum(abs.(dets .- 1)) < 8eps()

            # holes persist: the image of the hole keeps its area and holds no marker
            @test count(ζ -> inside_polygon(ex.circle, ζ), ex.Z0) == 0
            img = ψ.(ex.circle)
            A0, A1 = shoelace_area(ex.circle), shoelace_area(img)
            @printf("  %s: |ψ(z_k) − z'_k| ≤ %.1e, moved up to %.2f, |det Dψ − 1| ≤ %.1e, hole area %.6f → %.6f\n",
                name, err, moved, maximum(abs.(dets .- 1)), A0, A1)
            @test A0 ≈ 0.384845 atol = 5e-7      # the polygon area quoted in the note
            @test abs(A1 - A0) / A0 < 1e-5
            @test count(ζ -> inside_polygon(img, ζ), Z1) == 0
        end
    end

    @testset "must fail: a weight-dependent Hamiltonian has no single map" begin
        # two markers at the same point with the weights 0.5 and 2
        Z = vcat(ex.Z0, [[0.5, -0.4], [0.5, -0.4]])
        w = vcat(ones(length(ex.Z0)), [0.5, 2.0])
        Hw(Z) = (
            Λ = sum(ex.φ, Z) / length(Z); sum(k -> ex.h_w(Z[k], w[k], Λ), eachindex(Z)))
        Zn, _ = ensemble_implicit_midpoint(Hw, Z, ex.τ)
        gap = norm(Zn[end] - Zn[end - 1])
        @printf("  control: images of two coincident markers differ by %.2f\n", gap)
        @test gap > 0.1                    # the note quotes 0.14
    end
end
