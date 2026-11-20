# Verification run of the symplectic particle redistribution (testsets R1–R5 of the note):
# prints every number that the note and the research log quote, next to the quoted value, and
# checks the statements and the controls.
#
#     julia --startup-file=no --project=scripts scripts/redistribution/verify_symplectic_redistribution.jl
#
# Two set-ups of the relaxation are run: the note's script (step 0.006, ten steps) and the one
# that produced the numbers of the research log (step 0.02, eight steps).

using ParticleRemapping
using ForwardDiff
using LinearAlgebra
using Printf
using Test

const J16 = canonical_J(16)

row(name, value, quoted = "") = @printf("  %-52s %12.4e   %s\n", name, value, quoted)

# fourth-order (Richardson) central differences, in the order of terms of the note's script
function fd4(F, z, steps)
    D = zeros(length(F(z)), length(z))
    for j in eachindex(z)
        e = zeros(length(z))
        e[j] = steps[j]
        D[:, j] = (-F(z .+ 2e) .+ 8.0 .* F(z .+ e) .- 8.0 .* F(z .- e) .+ F(z .- 2e)) ./
                  (12.0 * e[j])
    end
    D
end

# R1
function r1()
    m = 1.7
    y0 = [0.4, -0.8, 0.25, 1.9, 1.3, 0.40]
    Ω = presymplectic_matrix(y0; m = m)
    fd = norm(kinematic_pullback(y0; m = m, jacobian = (F, y) -> fd_jacobian(F, y; order = 4)) -
              Ω)
    ad = norm(kinematic_pullback(y0; m = m) - Ω)
    vth = 1.3
    f(v) = isotropic_maxwellian(v; vth = vth)
    D = velocity_fold_determinant(f, vth)
    println("R1  pull-back of ω and velocity–momentum fold")
    row("‖π*ω − ω_analytic‖, fourth-order differences h = 1e-5", fd, "note: < 1e-13")
    row("‖π*ω − ω_analytic‖, automatic differentiation", ad)
    row("det D_v P at |v| = v_th", D)
    row("det D_v P at |v| = 0.5 v_th", velocity_fold_determinant(f, 0.5vth), "> 0")
    row("det D_v P at |v| = 1.5 v_th", velocity_fold_determinant(f, 1.5vth), "< 0")
    @test ad <= 8 * eps() * norm(Ω)
    @test fd < 100 * eps() / 1e-5
    @test abs(D) < 1e-15
    (; fd, ad)
end

# R2–R4 on one relaxation set-up
function r2_r4(; relax_step = 0.006, relax_iterations = 10, nx0 = 4, nv0 = 4)
    ex = redistribution_example(; nx0 = nx0, nv0 = nv0)
    mk = ex.markers
    q, v, w, m = mk.x, mk.v, mk.w, mk.m
    p = w .* m .* v
    N = length(q)
    ps = PhaseSpaceSplines()
    l2 = spline_moments(ps, l2_spline_projection(ps, q, v, w); m = m)
    relaxed = declustering_relaxation(q, v; Lx = ex.Lx, σx = ex.σx, σv = ex.σv,
        step = relax_step, iterations = relax_iterations)
    q̄, v̄ = relaxed.x, relaxed.v
    cv(xs, vs) = kernel_density_cv(mod.(xs, ex.Lx), vs; Lx = ex.Lx, σx = ex.σx, σv = ex.σv)
    ΔΩ = number_density_volumes(ps, q̄, v̄)
    kkt = kkt_redistribution_weights(ps, lumped_spline_projection(ps, q, v, w), q̄, v̄, ΔΩ,
        sum(w), sum(p); m = m)
    w̄ = kkt.weights
    p̄ = w̄ .* m .* v̄
    p̄_raw = kkt.raw_weights .* m .* v̄
    out = (; N, cv_before = cv(q, v), cv_after = cv(q̄, v̄), min_l2 = l2.minimum,
        mass_l2 = abs(l2.mass - sum(w)), mom_l2 = abs(l2.momentum - sum(p)),
        min_w̄ = minimum(w̄), mass_kkt = abs(sum(w̄) - sum(w)), mom_kkt = abs(sum(p̄) - sum(p)),
        raw_defect = abs(sum(p̄_raw) - sum(p)))
    N == 16 || return out

    p_scale = 2.5 * maximum(max.(abs.(p), abs.(p̄)))
    H = PairwiseSplineHamiltonian(ex.Lx, p_scale; nx_cells = 10, np_cells = 8)
    S = redistribution_scaling(N, p_scale)
    G_raw = redistribution_parameter_matrix(H, ImplicitMidpointStep(), q, p, q̄, p̄_raw)
    Δz_raw = vcat(q̄ .- q, p̄_raw .- p)
    inconsistent = norm(G_raw * (G_raw \ Δz_raw) - Δz_raw)
    G = redistribution_parameter_matrix(H, ImplicitMidpointStep(), q, p, q̄, p̄)
    θ_mid = solve_redistribution_parameters(H, ImplicitMidpointStep(), q, p, q̄, p̄)
    θ_se = solve_redistribution_parameters(H, SymplecticEulerStep(), q, p, q̄, p̄)
    q_mid, p_mid = 0.5 .* (q .+ q̄), 0.5 .* (p .+ p̄)
    Hq, Hp, _ = pairwise_hamiltonian_derivatives(H, q_mid, p_mid, θ_mid)
    Hq_s, Hp_s, _ = pairwise_hamiltonian_derivatives(H, q_mid .+ 1.2345, p_mid, θ_mid)
    perm = sortperm([sin(17.3 * k) for k in 1:N])
    Hq_p, Hp_p, _ = pairwise_hamiltonian_derivatives(H, q_mid[perm], p_mid[perm], θ_mid)
    z_unc = uncoupled_control_step(vcat(q, p), p_scale)
    D_adhoc = ForwardDiff.jacobian(z -> mean_force_control_step(z, p_scale), vcat(q, p))

    z̄ = vcat(q̄, p̄)
    mid = forward_redistribution_step(H, ImplicitMidpointStep(), q, p, θ_mid)
    se = forward_redistribution_step(H, SymplecticEulerStep(), q, p, θ_se)
    se_fine = forward_redistribution_step(H, SymplecticEulerStep(), q, p, θ_se;
        homotopy = (1:1000) ./ 1000)
    se2 = redistribution_substeps(H, SymplecticEulerStep(), q, p, q̄, p̄; substeps = 2)
    mid4 = redistribution_substeps(H, ImplicitMidpointStep(), q, p, q̄, p̄; substeps = 4)
    steps = vcat(fill(1e-6, N), fill(1e-6 * p_scale, N))
    D_mid = fd4(
        z -> forward_redistribution_step(H, ImplicitMidpointStep(), z[1:N],
            z[(N + 1):(2N)], θ_mid).z,
        vcat(q, p),
        steps)
    D_se = fd4(
        z -> forward_redistribution_step(H, SymplecticEulerStep(), z[1:N],
            z[(N + 1):(2N)], θ_se).z,
        vcat(q, p),
        steps)
    symp(D) = norm(transpose(D) * J16 * D - J16)
    merge(out,
        (; p_scale, rank_G = rank(Diagonal(S) * G), inconsistent,
            noether = abs(sum(Hq)), shift = max(norm(Hq_s - Hq), norm(Hp_s - Hp)),
            perm = max(norm(Hq_p - Hq[perm]), norm(Hp_p - Hp[perm])),
            control_a = abs(sum(z_unc[(N + 1):(2N)]) - sum(p)),
            control_b = symp(D_adhoc),
            res_mid = mid.residual, res_se = se.residual,
            hit_mid = norm(S .* (mid.z - z̄)), hit_se = norm(S .* (se.z - z̄)),
            hit_se2 = norm(S .* (se2.z - z̄)), hit_mid4 = norm(S .* (mid4.z - z̄)),
            hit_se_fine = norm(S .* (se_fine.z - z̄)),
            se_target_residual = norm(p̄ .- p .+
                                      pairwise_hamiltonian_derivatives(H, q, p̄, θ_se)[1]) /
                                 p_scale,
            κ_mid = redistribution_contraction_number(
                H, ImplicitMidpointStep(), q, p, q̄, p̄,
                θ_mid),
            κ_se = redistribution_contraction_number(H, SymplecticEulerStep(), q, p, q̄, p̄,
                θ_se),
            κ_se2 = maximum(se2.contraction), κ_mid4 = maximum(mid4.contraction),
            ift_fd_mid = norm(mid.jacobian - D_mid) / norm(mid.jacobian),
            ift_fd_se = norm(se.jacobian - D_se) / norm(se.jacobian),
            symp_mid_ift = symp(mid.jacobian), symp_se_ift = symp(se.jacobian),
            symp_mid_fd = symp(D_mid), symp_se_fd = symp(D_se),
            mom_mid = abs(sum(mid.z[(N + 1):(2N)]) - sum(p)),
            mom_se = abs(sum(se.z[(N + 1):(2N)]) - sum(p))))
end

function print_r2_r4(r, label; log = false)
    L(s) = log ? "log: " * s : ""
    println("R2  ($label)")
    row("kernel-density CV before relaxation", r.cv_before, L("0.220"))
    row("kernel-density CV after relaxation", r.cv_after, L("0.158"))
    row("|∫ f_h^{L²} − Σ w|", r.mass_l2, "note: < 1e-14")
    row("|∫ m v f_h^{L²} − Σ p|", r.mom_l2, "note: < 1e-14")
    row("Control D: min f_h^{L²}", r.min_l2, "< 0")
    row("min w̄ (KKT)", r.min_w̄, L("0.0227"))
    row("|Σ w̄ − Σ w| (KKT)", r.mass_kkt, L("1.11e-16"))
    row("|Σ p̄ − Σ p| (KKT)", r.mom_kkt, L("0.0"))
    row("Control C: |Σ p̄_raw − Σ p|", r.raw_defect, L("0.215") * "  note: > 1e-4")
    println("R3  ($label)")
    row("rank S G (2 N d − d = 31)", r.rank_G)
    row("Control C: ‖G θ − Δz_raw‖ (least squares)", r.inconsistent, "note: > 1e-4")
    row("|Σ_a ∇_{q_a} H_θ|", r.noether, L("4.86e-17") * "  note: < 1e-15")
    row("translation: ‖∇H(q + ξ) − ∇H(q)‖", r.shift, "note: < 1e-14")
    row("permutation: ‖∇H(σz) − σ∇H(z)‖", r.perm, L("2.10e-16") * "  note: < 1e-14")
    row("Control A: |Σ p̄ − Σ p|, uncoupled", r.control_a, "note: > 1e-2")
    row("Control B: ‖DΨᵀ J DΨ − J‖, mean-force subtraction", r.control_b, "note: > 1e-2")
    println("R4  ($label)")
    row("midpoint: Newton residual", r.res_mid, "note: < 1e-13")
    row("midpoint: ‖S (z_end − z̄)‖", r.hit_mid, L("8.04e-15") * "  note: < 1e-13")
    row("midpoint: (τ/2)‖S J ∇²H S⁻¹‖ at z_1/2", r.κ_mid, "bound: < 1")
    row("midpoint, 4 substeps: ‖S (z_end − z̄)‖", r.hit_mid4)
    row("midpoint, 4 substeps: largest contraction number", r.κ_mid4)
    row("midpoint: IFT vs FD Jacobian", r.ift_fd_mid, L("6.48e-10") * "  note: < 1e-9")
    row("midpoint: ‖DΨᵀJDΨ − J‖ (IFT)", r.symp_mid_ift, "note: < 1e-13")
    row("midpoint: ‖DΨᵀJDΨ − J‖ (FD)", r.symp_mid_fd, L("2.73e-8") * "  note: < 1e-7")
    row("midpoint: |Σ p_end − Σ p|", r.mom_mid, L("1.67e-16") * "  note: < 1e-15")
    row("SE: Newton residual", r.res_se, "note: < 1e-13")
    row("SE: ‖S (z_end − z̄)‖", r.hit_se, L("0.150 (caustic)") * "  note: < 1e-13")
    row("SE: ‖τ ∇_q ∇_p̄ H‖ at (q, p̄)", r.κ_se, "bound: < 1")
    row("SE, 1000 homotopy stages: ‖S (z_end − z̄)‖", r.hit_se_fine)
    row("SE: residual of the step equations at z̄", r.se_target_residual)
    row("SE, 2 substeps: ‖S (z_end − z̄)‖", r.hit_se2)
    row("SE, 2 substeps: largest contraction number", r.κ_se2)
    row("SE: IFT vs FD Jacobian", r.ift_fd_se, L("2.47e-10") * "  note: < 1e-9")
    row("SE: ‖DΨᵀJDΨ − J‖ (IFT)", r.symp_se_ift, "note: < 1e-13")
    row("SE: ‖DΨᵀJDΨ − J‖ (FD)", r.symp_se_fd, L("2.90e-8") * "  note: < 1e-7")
    row("SE: |Σ p_end − Σ p|", r.mom_se, "note: < 1e-15")
end

function check_r2_r4(r)
    @test r.mass_l2 < 1e-14
    @test r.mom_l2 < 1e-14
    @test r.min_l2 < 0
    @test r.raw_defect > 1e-4
    @test r.min_w̄ > 0
    @test r.mass_kkt < 1e-14
    @test r.mom_kkt < 1e-14
    @test r.cv_after < r.cv_before
    @test r.rank_G == 31
    @test r.inconsistent > 1e-4
    @test r.noether < 1e-15
    @test r.shift < 1e-14
    @test r.perm < 1e-14
    @test r.control_a > 1e-2
    @test r.control_b > 1e-2
    @test r.res_mid < 1e-13
    @test r.res_se < 1e-13
    @test r.hit_mid < 1e-13
    @test r.hit_se > 0.1 && r.κ_se > 1          # the caustic control at one step
    @test r.hit_se_fine > 0.1                    # not a jump of the homotopy
    @test r.se_target_residual < 1e-13           # z̄ is a root the branch does not reach
    @test r.hit_se2 < 1e-13 && r.κ_se2 < 1      # two substeps restore the bound
    @test r.hit_mid4 < 1e-13 && r.κ_mid4 < 1
    @test r.ift_fd_mid < 1e-9
    @test r.ift_fd_se < 1e-9
    @test r.symp_mid_ift < 1e-13 && r.symp_mid_fd < 1e-7
    @test r.symp_se_ift < 1e-13 && r.symp_se_fd < 1e-7
    @test r.mom_mid < 1e-15
    @test r.mom_se < 1e-15
end

# R5
function r5()
    ex = closed_loop_example()
    z, H = ex.z, ex.H
    n = length(z) ÷ 2
    J4 = canonical_J(n)
    steps = vcat(fill(1e-5, n), fill(1e-5 * H.p_scale, n))
    cl = target_matching_closed_loop(H, z)
    lhs, rhs = closed_loop_defect(cl.jacobians...)
    DΦ = fd4(y -> target_matching_closed_loop(H, y; jacobians = false).z, z, steps)
    DΘ = fd4(y -> target_matching_parameters(H, y), z, steps)
    Dg = fd4(y -> target_matching_closed_loop(H, y; jacobians = false).g, z, steps)
    lhs_fd, rhs_fd = closed_loop_defect(DΦ, DΘ, Dg)
    cl_se = target_matching_closed_loop(H, z; scheme = SymplecticEulerStep())
    lhs_se, rhs_se = closed_loop_defect(cl_se.jacobians...)
    vl = variational_closed_loop(z)
    DΦv, DΘv = vl.jacobians
    DΦv_fd = fd4(y -> variational_closed_loop(y; jacobians = false).z, z, steps)
    r = (; defect = norm(lhs), relerr = norm(lhs - rhs) / norm(lhs),
        defect_fd = norm(lhs_fd), relerr_fd = norm(lhs_fd - rhs_fd) / norm(lhs_fd),
        defect_se = norm(lhs_se), relerr_se = norm(lhs_se - rhs_se) / norm(lhs_se),
        var = norm(transpose(DΦv) * J4 * DΦv - J4),
        var_fd = norm(transpose(DΦv_fd) * J4 * DΦv_fd - J4), DΘ = norm(DΘv),
        mom = abs(sum(vl.z[(n + 1):(2n)]) - sum(z[(n + 1):(2n)])))
    println("R5  closed loop (4 markers, periodic quartic basis with 5 cells)")
    row("Control E: ‖DΦᵀJDΦ − J‖ (midpoint, exact)", r.defect, "note: > 1e-3")
    row("Control E: relative error of the identity (exact)", r.relerr)
    row("Control E: ‖DΦᵀJDΦ − J‖ (midpoint, FD)", r.defect_fd)
    row("Control E: relative error of the identity (FD)", r.relerr_fd, "note: < 1e-8")
    row("Control E: ‖DΦᵀJDΦ − J‖ (symplectic Euler, exact)", r.defect_se)
    row("Control E: relative error (symplectic Euler, exact)", r.relerr_se)
    row("variational closure: ‖DΦᵀJDΦ − J‖ (exact)", r.var, "note: < 1e-12")
    row("variational closure: ‖DΦᵀJDΦ − J‖ (FD)", r.var_fd, "note: < 1e-12")
    row("variational closure: ‖DΘ‖", r.DΘ, "note: > 1e-2")
    row("variational closure: |Σ p̄ − Σ p|", r.mom, "note: < 1e-15")
    @test r.defect > 1e-3
    @test r.relerr < 1e3 * eps()
    @test r.relerr_fd < 1e-8
    @test r.defect_se > 1e-3
    @test r.relerr_se < 1e3 * eps()
    @test r.var < 1e-12
    @test r.var_fd < 1e-8
    @test r.DΘ > 1e-2
    @test r.mom < 1e-15
    r
end

t = @elapsed @testset "Symplectic particle redistribution" begin
    r1()
    println()
    a = r2_r4()
    print_r2_r4(a, "script: relaxation step 0.006, 10 steps")
    check_r2_r4(a)
    println()
    b = r2_r4(; relax_step = 0.02, relax_iterations = 8)
    print_r2_r4(b, "research log: relaxation step 0.02, 8 steps"; log = true)
    check_r2_r4(b)
    println()
    c = r2_r4(; nx0 = 6, nv0 = 6)
    println("R2  (36 markers, 6 × 6 lattice)")
    row("Control D: min f_h^{L²}", c.min_l2, "log: −0.140")
    @test c.min_l2 < 0
    println()
    r5()
end
@printf("\nruntime %.1f s\n", t)
flush(stdout)
