# The stages of the remap map on small ensembles, d = 1 and 2: kicks with the centroid shift and
# drifts (symplectic, momentum-conserving, one symplectic Euler step of H^K + H^D), the stage
# Jacobians in (x, v), the reweighting at fixed (x, p), and the non-separable implicit-midpoint
# stages in their velocity form. Every case that must fail runs: the naive kick, parameters that
# depend on the state, the trapezoidal rule, the stage without the centroid shift, and the wrong
# forms of the stage Jacobian.

using ParticleRemapping
using ForwardDiff
using LinearAlgebra
using Printf
using Random
using Test

const m = 1.3
const L = 2π

# Random smooth stage functions, periodic in ξ with period L and depending on ln w. `amp` scales
# the kick potential; it is the parameter of the check of state-dependent parameters.
function random_stage(d)
    K = [rand(-2:2, d) for _ in 1:4]
    b = randn(4)
    a = 0.3 * randn(4)
    φ = 2π * rand(4)
    Q = [randn(d) for _ in 1:4]
    e = randn(4)
    c = 0.3 * randn(4)
    ψ = 2π * rand(4)
    κ(ξ, w, amp) = amp * sum(a[j] * sin(dot(K[j], ξ) * 2π / L + b[j] * log(w) + φ[j])
    for j in 1:4)
    η(v, w) = sum(c[j] * cos(dot(Q[j], v) + e[j] * log(w) + ψ[j]) for j in 1:4) +
              0.05 * log(w) * sum(abs2, v)
    gκ(amp) = (ξ, w) -> ForwardDiff.gradient(s -> κ(s, w, amp), ξ)
    gη = (v, w) -> ForwardDiff.gradient(u -> η(u, w), v)
    (; κ, η, gκ, gη)
end

function R_stages(z, d, N, w, stages; amp = 1.0, kickfun = centroid_kick)
    x, p = unpack(z, d, N)
    for st in stages
        p = kickfun(x, p, w, st.gκ(amp))
        x = velocity_drift(x, p, w, m, st.gη)
    end
    pack(x, p)
end

# The tolerance of a symplecticity defect: round-off in a product of two Jacobians of norm ‖DR‖.
defect_tol(nDR) = 100 * eps() * max(1, nDR)^2

Random.seed!(20261119)

@testset "Kick and drift" begin
    for d in (1, 2)
        N = 6
        w = 0.1 .+ 2 * rand(N)
        x0 = L * rand(d, N)
        p0 = m * w' .* randn(d, N)
        z0 = pack(x0, p0)
        stages = [random_stage(d), random_stage(d)]
        R(z) = R_stages(z, d, N, w, stages)

        # the note quotes a defect ≤ 2.9e-15 and |ΔP| ≤ 8.9e-16
        @testset "d = $d, symplectic and momentum-conserving" begin
            defect, nDR = symplecticity_defect(R, z0)
            x1, p1 = unpack(R(z0), d, N)
            ΔP = maximum(abs, sum(p1; dims = 2) - sum(p0; dims = 2))
            @printf("  d = %d: kick–drift symplecticity defect %.1e (‖DR‖ = %.1f), |ΔP| = %.1e\n",
                d, defect, nDR, ΔP)
            @test defect < defect_tol(nDR)
            @test ΔP < 100 * eps() * max(1, sum(abs, p0) + sum(abs, p1))
            @test maximum(abs, p1 - p0) > 1e-2
        end

        @testset "d = $d, the momentum correction is a Galilean boost" begin
            st = stages[1]
            α = w ./ sum(w)
            Xc = x0 * α
            G = reduce(hcat, [st.gκ(1.0)(x0[:, k] .- Xc, w[k]) for k in 1:N])
            S = sum(G; dims = 2)
            boost = (S .* α') ./ (m .* w')
            # equal to S / (m w_tot) for every marker, a closed form independent of k
            @printf("  d = %d: boost spread over the markers %.1e\n", d,
                maximum(abs, boost .- boost[:, 1]))
            @test boost ≈ repeat(S ./ (m * sum(w)), 1, N) rtol = 8eps()
            @test maximum(abs, boost .- boost[:, 1]) <
                  100 * eps() * max(1, maximum(abs, boost))
            @test maximum(abs, boost) > 1e-3
        end

        @testset "d = $d, a kick and a drift are symplectic Euler of H^K + H^D" begin
            st = stages[1]
            HK(xv) = (x = reshape(xv, d, N); α = w ./ sum(w); Xc = x * α;
                sum(k -> st.κ(x[:, k] .- Xc, w[k], 1.0), 1:N))
            HD(pv) = (
                p = reshape(pv, d, N); sum(k -> st.η(p[:, k] ./ (m * w[k]), w[k]), 1:N))
            p_se = p0 .- reshape(ForwardDiff.gradient(HK, vec(x0)), d, N)
            x_se = x0 .+ reshape(ForwardDiff.gradient(HD, vec(p_se)), d, N)
            p_kd = centroid_kick(x0, p0, w, st.gκ(1.0))
            x_kd = velocity_drift(x0, p_kd, w, m, st.gη)
            err = max(maximum(abs, p_se - p_kd), maximum(abs, x_se - x_kd))
            @printf("  d = %d: |kick–drift − symplectic Euler of H^K + H^D| = %.1e\n", d,
                err)
            @test err < 1e-13                     # the note quotes 3.1e-15
        end

        @testset "d = $d, must fail: the naive kick" begin
            Rn(z) = R_stages(z, d, N, w, stages; kickfun = naive_kick)
            defect, _ = symplecticity_defect(Rn, z0)
            x1, p1 = unpack(Rn(z0), d, N)
            ΔP = maximum(abs, sum(p1; dims = 2) - sum(p0; dims = 2))
            @printf("  d = %d: naive kick: defect %.2f, |ΔP| = %.1e\n", d, defect, ΔP)
            @test defect > (d == 1 ? 0.85 : 1.5)   # the note quotes 0.86 and 1.6
            @test ΔP < 1e-12
        end

        @testset "d = $d, must fail: parameters that depend on the state" begin
            θ(z) = 1 + 0.5 * sum(sin, unpack(z, d, N)[1]) / N
            Rstate(z) = R_stages(z, d, N, w, stages; amp = θ(z))
            θ0 = θ(z0)
            Rfrozen(z) = R_stages(z, d, N, w, stages; amp = θ0)
            defect_state, _ = symplecticity_defect(Rstate, z0)
            defect_frozen, nDR = symplecticity_defect(Rfrozen, z0)
            @printf("  d = %d: θ(z): defect %.2f; θ frozen: defect %.1e\n", d, defect_state,
                defect_frozen)
            @test defect_state > 0.15             # the note quotes 0.17 and 0.16
            @test defect_frozen < defect_tol(nDR)
        end

        if d == 1
            @testset "d = 1, stage Jacobians in (x, v)" begin
                W(x, v) = (1 + 0.3 * cos(x)) * exp(-v^2 / 2)
                Δv(x, ℓ) = 0.2 * sin(x + 0.7) * (1 + 0.4ℓ + 0.1ℓ^2)
                Δx(v, ℓ) = 0.3 * cos(1.3v) * sin(0.5ℓ)
                kick_xv(z) = [z[1], z[2] + Δv(z[1], log(W(z[1], z[2])))]
                drift_xv(z) = [z[1] + Δx(z[2], log(W(z[1], z[2]))), z[2]]
                err = 0.0
                for _ in 1:20
                    z = [L * rand(), 3 * randn()]
                    dK, dD = kick_drift_jacobian_determinants(Δv, Δx, W, z)
                    err = max(err, abs(det(ForwardDiff.jacobian(kick_xv, z)) - dK),
                        abs(det(ForwardDiff.jacobian(drift_xv, z)) - dD))
                end
                @test err < 8eps()                # the note quotes 2.2e-16
                # control: without a dependence on ℓ the kick preserves the area exactly
                kick_c(z) = [z[1], z[2] + Δv(z[1], 0.0)]
                @test det(ForwardDiff.jacobian(kick_c, [1.0, 0.5])) == 1
                # the kick couples v and −v: equal velocity change, opposite compression
                zp, zm = [π / 2, 1.2], [π / 2, -1.2]
                dKp = det(ForwardDiff.jacobian(kick_xv, zp)) - 1
                dKm = det(ForwardDiff.jacobian(kick_xv, zm)) - 1
                @printf("  stage Jacobians in (x, v): error %.1e; kick at ±v: det − 1 = %+.3f, %+.3f\n",
                    err, dKp, dKm)
                @test kick_xv(zp)[2] - zp[2] ≈ kick_xv(zm)[2] - zm[2] rtol = 1e-14
                @test dKp ≈ -dKm rtol = 1e-12
                @test abs(dKp) ≈ 0.047 atol = 5e-4  # the note quotes ∓0.047
            end
        end

        @testset "d = $d, reweighting at fixed (x, p)" begin
            w_new = 0.1 .+ 2 * rand(N)
            v_old = p0 ./ (m .* w')
            v_new = p0 ./ (m .* w_new')
            @test maximum(abs, v_new - (w ./ w_new)' .* v_old) < 1e-14
            @test sum(m .* w_new' .* v_new; dims = 2) ≈ sum(p0; dims = 2) rtol = 1e-14
        end
    end
end

# ---- non-separable stages ----------------------------------------------------------------------

# A random smooth φ(ξ, v, w), periodic in ξ with period L, with terms that couple ξ, v and ln w.
function random_phi(d; amp = 0.12)
    K = [rand(-2:2, d) for _ in 1:4]
    Q = [randn(d) for _ in 1:4]
    b = randn(4)
    a = amp * randn(4)
    ϕ = 2π * rand(4)
    (ξ, v, w) -> sum(a[j] *
                     sin(dot(K[j], ξ) * 2π / L + dot(Q[j], v) + b[j] * log(w) + ϕ[j])
    for j in 1:4) + 0.02 * log(w) * sum(abs2, v)
end

function gradients(φ)
    (ξ, v, w) -> (
        ForwardDiff.gradient(s -> φ(s, v, w), ξ), ForwardDiff.gradient(u -> φ(ξ, u, w), v))
end

function stage_map(φs, d, N, w; centroid = true)
    z -> begin
        x, p = unpack(z, d, N)
        for φ in φs
            if centroid
                x, p, _, ok = nonseparable_stage(x, p, w, m, gradients(φ))
                ok || error("nonseparable_stage did not converge")
            else
                z1 = hamiltonian_one_step(
                    nonseparable_hamiltonian(φ, d, N, w; m,
                        centroid = false), pack(x, p))
                x, p = unpack(z1, d, N)
            end
        end
        pack(x, p)
    end
end

Random.seed!(20261120)

@testset "Non-separable stages" begin
    for d in (1, 2)
        N = 7
        w = 0.1 .+ 2 * rand(N)
        x0 = L * rand(d, N)
        p0 = m * w' .* randn(d, N)
        z0 = pack(x0, p0)
        φs = [random_phi(d), random_phi(d)]

        @testset "d = $d, the velocity form is the implicit midpoint rule of H" begin
            x1, p1, its, ok = nonseparable_stage(x0, p0, w, m, gradients(φs[1]))
            H = nonseparable_hamiltonian(φs[1], d, N, w; m)
            z1 = hamiltonian_one_step(H, z0)
            J = canonical_J(d * N)
            z1v = pack(x1, p1)
            res = maximum(abs, z1v - z0 - J * ForwardDiff.gradient(H, (z0 + z1v) / 2))
            @printf("  d = %d: %d iterations; |velocity form − implicit midpoint| %.1e, residual %.1e\n",
                d, its, maximum(abs, z1v - z1), res)
            @test ok
            @test maximum(abs, z1v - z1) < 1e-13   # the note quotes 4.4e-16
            @test res < 1e-13
            @test maximum(abs, z1v - z0) > 1e-2
        end

        @testset "d = $d, symplectic and momentum-conserving" begin
            R = stage_map(φs, d, N, w)
            defect, nDR = symplecticity_defect(R, z0)
            x1, p1 = unpack(R(z0), d, N)
            ΔP = maximum(abs, sum(p1; dims = 2) - sum(p0; dims = 2))
            @printf("  d = %d: two non-separable stages: defect %.1e (‖DR‖ = %.1f), |ΔP| = %.1e\n",
                d, defect, nDR, ΔP)
            @test defect < defect_tol(nDR)          # the note quotes 1.1e-15
            @test ΔP < 100 * eps() * max(1, sum(abs, p0) + sum(abs, p1))
        end

        @testset "d = $d, cases that must fail" begin
            # the trapezoidal rule keeps the linear invariant P and loses symplecticity
            H = nonseparable_hamiltonian(φs[1], d, N, w; m)
            T(z) = hamiltonian_one_step(H, z; rule = :trapezoidal)
            defect_T, _ = symplecticity_defect(T, z0)
            x1, p1 = unpack(T(z0), d, N)
            ΔP_T = maximum(abs, sum(p1; dims = 2) - sum(p0; dims = 2))
            # without the centroid shift the map is symplectic and loses P
            Rn = stage_map(φs, d, N, w; centroid = false)
            defect_n, nDRn = symplecticity_defect(Rn, z0)
            x1, p1 = unpack(Rn(z0), d, N)
            ΔP_n = maximum(abs, sum(p1; dims = 2) - sum(p0; dims = 2))
            @printf("  d = %d: trapezoidal defect %.1e with |ΔP| %.1e; no centroid shift: defect %.1e, |ΔP| %.2f\n",
                d, defect_T, ΔP_T, defect_n, ΔP_n)
            @test defect_T > (d == 1 ? 7e-4 : 3.5e-4)   # the note quotes 7.5e-4 and 3.6e-4
            @test ΔP_T < 1e-13
            @test defect_n < defect_tol(nDRn)
            @test ΔP_n > (d == 1 ? 0.28 : 1.6)        # the note quotes 0.29 and 1.7
        end

        @testset "d = $d, stage Jacobian and divergence in (x, v)" begin
            φ = random_phi(d; amp = 0.3)
            A = randn(2d)
            W(y) = exp(0.4 * sin(dot(A[1:d], y[1:d]) + 0.7) -
                       0.3 * sum(abs2, y[(d + 1):end]) +
                       0.2 * dot(A[(d + 1):end], y[(d + 1):end]))
            err_det = err_first = err_wrong = err_div = 0.0
            for _ in 1:20
                y = vcat(L * rand(d), 1.5 * randn(d))
                y1 = marker_midpoint_step(φ, W, y)
                detD = det(ForwardDiff.jacobian(s -> marker_midpoint_step(φ, W, s), y))
                f = nonseparable_jacobian_determinant(φ, W, y, y1)
                err_det = max(err_det, abs(detD - f.exact))
                err_first = max(err_first, abs(detD - f.first))
                err_wrong = max(err_wrong, abs(detD - f.wrong))
                div, br = nonseparable_divergence(φ, W, y)
                err_div = max(err_div, abs(div - br))
            end
            @printf("  d = %d: |det D − formula| %.1e; first order ≥ %.3f; wrong sign ≥ %.3f; |div u − {W, ∂_wφ}| %.1e\n",
                d, err_det, err_first, err_wrong, err_div)
            @test err_det < 1e-13                   # the note quotes 4.4e-15
            @test err_div < 1e-13                   # the note quotes 8.9e-16
            @test err_first > 0.3                   # must fail; measured ≥ 0.315
            @test err_wrong > 0.45                  # must fail; measured ≥ 0.475
        end
    end

    @testset "d = 1, the ±v coupling at ∂_x W = 0" begin
        c = pm_v_coupling(; vs = [1.2, -1.2])
        @printf("  x = π: separable det D − 1 = %+.3f, %+.3f; non-separable %+.3f, %+.3f at v = ±1.2\n",
            c.separable..., c.nonseparable...)
        # a separable φ: opposite signs at ±v (the note quotes ∓0.120)
        @test c.separable[1] * c.separable[2] < 0
        @test c.separable[1] ≈ -0.120 atol = 5e-4
        @test c.separable[2] ≈ 0.120 atol = 5e-4
        # the non-separable φ: the same sign (the note quotes −0.144 at both)
        @test c.nonseparable[1] * c.nonseparable[2] > 0
        @test c.nonseparable[1] ≈ -0.144 atol = 5e-4
        @test c.nonseparable[2] ≈ -0.144 atol = 5e-4
    end
end

@testset "Guards of the stage solvers" begin
    H(z) = 0.5 * sum(abs2, z) + 3 * sum(z .^ 4)
    @test_throws r"^ArgumentError: rule must be :midpoint or :trapezoidal" hamiltonian_one_step(
        H, [0.1, 0.2]; rule = :euler)
    # a step too large for the fixed-point iteration
    @test_throws r"did not converge" hamiltonian_one_step(z -> 40 * H(z), [1.0, 1.0];
        maxit = 20)
    @test_throws r"did not converge" marker_midpoint_step(
        (ξ, v, w) -> 40 * (ξ[1]^4 + v[1]^4),
        y -> 1.0, [1.0, 1.0]; maxit = 20)
    # without convergence the non-separable stage returns its input and `false`
    gφ(ξ, v, w) = ([40 * ξ[1]^3], [40 * v[1]^3])
    x, p = [1.0 -1.0], [1.0 2.0]
    x1, p1, its, ok = nonseparable_stage(x, p, [1.0, 2.0], 1.0, gφ; maxit = 5)
    @test !ok && its == 5 && x1 == x && p1 == p
end
