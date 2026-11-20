using ParticleRemapping
using Test

# A number quoted with d decimals in the note agrees to half a unit of its last digit.
agrees(x, quoted, d) = abs(x - quoted) ≤ 0.5 * 10.0^(-d) + 8eps(abs(quoted))

@testset "numerical example of the gauge rezoning (6400 markers, 16 nodes)" begin
    setup = gauge_rezoning_setup()
    mk, grid, method = setup.markers, setup.grid, setup.method
    @test length(mk) == 6400
    @test method.solver === :kkt
    # max over the nodes of |M^(k)[W, V] − M^(k)[w, v]| / M^(k)[w, |v|]
    relres(after, k) = maximum(abs.(grid_moment(grid, after, k) .-
                                    grid_moment(grid, mk, k)) ./
                               grid_moment(grid, mk.x, abs.(mk.v), mk.w, k))

    r = remap(method, mk)
    after = r.markers
    @test agrees(r.info.τ, 0.594, 3)
    @test agrees(r.info.min_curvature / 1e-4, 8.5, 1)
    # (1) charge and energy constraints, to the Newton tolerance (note: 2.6e-15, 1.9e-15)
    @test relres(after, 0) < method.tol
    @test relres(after, 2) < method.tol
    # (2) the current, not constrained, to round-off (note: 1.2e-16)
    @test relres(after, 1) < 1e-14
    # (2) every momentum p_a = m W_a V_a (note: 2e-16)
    @test maximum(abs.(momentum(after) .- momentum(mk)) ./ abs.(momentum(mk))) < 8eps()
    # (2) the global kinetic energy, the sum of the local ones
    @test abs(kinetic_energy(after) / kinetic_energy(mk) - 1) < grid.Ng * method.tol
    # (3) the velocity gaps more even, 77.5 → 27.5, and every sign kept
    gb, ga = mean_gap_ratio(grid, mk.x, mk.v), mean_gap_ratio(grid, mk.x, after.v)
    @test ga < gb
    @test agrees(gb, 77.5, 1)
    @test agrees(ga, 27.5, 1)
    @test all(sign.(after.v) .== sign.(mk.v))
    λ = exp.(r.info.μ)
    @test agrees(minimum(λ), 0.29, 2)
    @test agrees(maximum(λ), 8.1, 1)
    # control: K = {0, 2} changes the third moment by 22 %
    @test agrees(relres(after, 3), 0.22, 2)
    # control: the unconstrained targets violate the charge by 45 % and spread perfectly
    Wt = mk.w .* exp.(-gauge_targets(grid, mk.x, mk.v))
    unconstrained = from_canonical(mk.x, momentum(mk), Wt; m = mk.m)
    @test agrees(relres(unconstrained, 0), 0.45, 2)
    @test agrees(mean_gap_ratio(grid, mk.x, unconstrained.v), 1.00, 2)
    # (4) K = {0, 2, 3}: all three constraints; the homotopy stops earlier
    m3 = gauge_rezoning_setup(; K = (0, 2, 3)).method
    r3 = remap(m3, mk)
    @test maximum(relres(r3.markers, k) for k in (0, 2, 3)) < m3.tol
    @test agrees(r3.info.τ, 0.402, 3)
    @test agrees(mean_gap_ratio(grid, mk.x, r3.markers.v), 39.2, 1)
    # the Schur-complement solver follows the same branch to the same fold
    rs = remap(gauge_rezoning_setup(; solver = :schur).method, mk)
    @test rs.info.τ == r.info.τ
    @test rs.info.μ ≈ r.info.μ atol=1e-10
end

@testset "numerical example of the canonical transport (2e5 markers)" begin
    setup = quantile_matching_setup()
    ex = quantile_matching_experiment(setup)
    @test length(ex.original) == 200_000
    err(mk) = binned_weight_error(setup.target, setup.grid, mk.x, mk.v, mk.w)
    e_orig = err(ex.original)
    e_id, e_c = err(ex.identity.markers), err(ex.admissible.markers)
    e_kept, e_b = err(ex.kept), err(ex.inadmissible.markers)
    rc, rb = ex.admissible.info.max_residual, ex.inadmissible.info.max_residual
    # (c) after ψ = id and after the admissible lift the matching is below the sampling noise
    @test e_id < e_orig
    @test e_c < e_orig
    # control: the admissible lift without matching (W = w) has an error ten times the noise
    @test e_kept > 3e_orig
    # control: after the inadmissible lift the compatibility residual is large
    @test rb > 10rc
    # the gauge identity m W V = P after the admissible lift
    mk = ex.admissible.markers
    @test maximum(abs.(mk.m .* mk.w .* mk.v .- ex.admissible.info.P) ./
                  abs.(ex.admissible.info.P)) < 8eps()
    # the table of the note
    @test agrees(e_orig, 0.0273, 4)
    @test agrees(e_id, 0.0113, 4)
    @test agrees(ex.identity.info.max_residual, 0.038, 3)
    @test agrees(e_c, 0.0117, 4)
    @test agrees(rc, 0.035, 3)
    @test agrees(e_kept, 0.2605, 4)
    @test agrees(e_b, 1.277, 3)
    @test agrees(rb, 5.4, 1)
    # the admissible lift moves markers by up to 1.85 in x
    d = maximum(abs.(mod.(ex.admissible.info.X .- ex.original.x .+ π, 2π) .- π))
    @test agrees(d, 1.85, 2)
    # the cell counts and charges follow (J(x) J(X))^{1/2} and J(X) of the closed form
    grid = setup.grid
    cells = cellof.(Ref(grid), mk.x)
    for c in 1:(grid.Ng)
        @test sum(mk.w[cells .== c]) / grid.Δx ≈ setup.target.Jbar[c] rtol=0.05
    end
end
