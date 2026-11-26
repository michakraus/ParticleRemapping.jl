# The parametrised maps on a small lattice: θ = 0 is the identity, every θ gives a symplectic map
# that conserves the total momentum, the controls B and C move every marker by one map of (x, v)
# (two coincident markers of different weights stay coincident, which the weight-dependent maps
# break; that this map preserves the area is the stage-Jacobian test of `stages.jl`), and the
# sensitivities of the non-separable map agree with forward-mode derivatives (and miss them
# without the terms of the centroid shift and the boost).

using ParticleRemapping
using ForwardDiff
using LinearAlgebra
using Random
using Test

Random.seed!(20261124)

prob = hole_filling_problem(; nL = 12, hole_a = 0.9, hole_b = 0.8)
x, p, w, m = prob.x0, prob.p0, prob.w0, prob.m
N = prob.N

function map_of(par, θ)
    z -> begin
        X, P = remap_stages(par, θ, z[1:N], z[(N + 1):end], w; m)
        vcat(X, P)
    end
end

@testset "Parametrised maps" begin
    kd = Dict(v => KickDriftParametrisation(w; variant = v, nξ = 8, nv = 10, nw = 4)
    for v in (:weight_dependent, :control_A, :control_B))
    ns = NonSeparableParametrisation(w; nξ = 6, nv = 7, nℓ = 3)
    nsC = NonSeparableParametrisation(w; nξ = 6, nv = 7, nℓ = 3, weighted = false)

    @testset "parameter counts" begin
        @test remap_parameter_count(kd[:weight_dependent]) == 2 * (8 + 10) * 4
        @test remap_parameter_count(kd[:control_A]) == 2 * (8 + 10)
        @test remap_parameter_count(ns) == 2 * 6 * 7 * 3
        @test_throws r"^ArgumentError: unknown kick–drift variant" KickDriftParametrisation(
            w;
            variant = :other)
    end

    for par in (values(kd)..., ns, nsC)
        n = remap_parameter_count(par)
        @testset "$(nameof(typeof(par))) $(par isa KickDriftParametrisation ? par.variant : par.weighted)" begin
            X0, P0 = remap_stages(par, zeros(n), x, p, w; m)
            @test X0 ≈ x rtol = 8eps()
            @test P0 ≈ p atol = 8eps() * maximum(abs, p)
            θ = 0.02 .* randn(n)
            X, P = remap_stages(par, θ, x, p, w; m)
            @test maximum(abs, X .- x) > 1e-3        # the map moves the markers
            @test abs(sum(P) - sum(p)) < 1e-14 * sum(abs, p)
            defect, nDR = symplecticity_defect(map_of(par, θ), vcat(x, p))
            @test defect < 100 * eps() * max(1, nDR)^2
        end
    end

    @testset "one map of (x, v) for the controls B and C" begin
        # two more markers at one point of (x, v), with two different weights of the lattice
        xq, vq = 2.0, 0.7
        wa, wb = w[findfirst(>(0.01), w)], w[findfirst(<(1e-3), w)]
        x2, v2, w2 = [x; xq; xq], [prob.v0; vq; vq], [w; wa; wb]
        p2 = m .* w2 .* v2
        gap(par) = begin
            θ = 0.02 .* randn(remap_parameter_count(par))
            X, P = remap_stages(par, θ, x2, p2, w2; m)
            V = P ./ (m .* w2)
            max(abs(X[end] - X[end - 1]), abs(V[end] - V[end - 1]))
        end
        kd2(v) = KickDriftParametrisation(w2; variant = v, nξ = 8, nv = 10, nw = 4)
        ns2(weighted) = NonSeparableParametrisation(w2; nξ = 6, nv = 7, nℓ = 3, weighted)
        @test gap(kd2(:control_B)) < 1e-13
        @test gap(ns2(false)) < 1e-13
        # the maps that depend on the weight separate the two markers: must fail
        @test gap(kd2(:weight_dependent)) > 1e-4
        @test gap(ns2(true)) > 1e-4
    end

    @testset "kick–drift history" begin
        par = kd[:weight_dependent]
        θ = 0.02 .* randn(remap_parameter_count(par))
        history = Tuple{Vector{Float64}, Vector{Float64}}[]
        X, P = kick_drift_map(par, θ, x, p, w; m, history)
        @test length(history) == 4
        @test history[1][1] == x                     # a kick keeps the positions
        @test history[2][2] == history[1][2]         # a drift keeps the momenta
        @test history[end] == (X, P)
    end

    @testset "non-separable sensitivities" begin
        θ = 0.02 .* randn(remap_parameter_count(ns))
        r = nonseparable_map(ns, θ, x, p, w; m, sensitivities = true)
        @test r.iterations < ns.iteration_limit
        err = nonseparable_sensitivity_error(ns, θ, x, p, w; m)
        err_nc = nonseparable_sensitivity_error(ns, θ, x, p, w; m, centroid_terms = false)
        @test err < 1e-10
        @test err_nc > 1e-4                          # must fail
        # a stage beyond the iteration limit returns NaN, which the fit rejects
        big = nonseparable_map(ns, 50 .* θ, x, p, w; m, maxit = 5)
        @test all(isnan, big.X) && big.sens === nothing
    end

    @testset "fit" begin
        tg = projection_targets(HatProjection(12, 13), x, prob.v0, w; prob.L, prob.V,
            density = (xs, vs) -> prob.density.(xs, vs))
        for par in (kd[:weight_dependent], ns)
            fr = fit_remap_parameters(par, x, p, w, tg.x, tg.v, tg.w; m, ℓx = prob.L / 12,
                ℓv = 2prob.V / 12, maxit = 5)
            @test fr.F < fr.F0
            @test 1 <= fr.iterations <= 5
            X, P = remap_stages(par, fr.θ, x, p, w; m)
            r = vcat((X .- tg.x) ./ (prob.L / 12), (P ./ (m .* tg.w) .- tg.v) ./
                                                   (2prob.V / 12))
            @test fr.F ≈ 0.5 * dot(r, r) + 0.5e-3 * dot(fr.θ, fr.θ) rtol = 1e-12
        end
    end
end

@testset "equal weights" begin
    # equal weights give the ln w basis a unit span, not a zero spacing, which would make every
    # evaluation NaN; the map is then the same area-preserving map for every marker
    we = fill(0.01, N)
    pe = p .* (we ./ w)
    for par in (KickDriftParametrisation(we; nξ = 8, nv = 10, nw = 4),
        NonSeparableParametrisation(we; nξ = 6, nv = 7, nℓ = 3))
        θp = 0.05 .* randn(remap_parameter_count(par))
        X, P = remap_stages(par, θp, x, pe, we; m)
        @test all(isfinite, X) && all(isfinite, P)
        @test X != x
        Φ = z -> vcat(remap_stages(par, θp, z[1:N], z[(N + 1):end], we; m)...)
        @test ParticleRemapping.symplecticity_defect(Φ, vcat(x, pe))[1] < 1e-10
    end
end

@testset "fit with Gauss–Newton" begin
    # the small problem of `remap.jl` with the targets of the exact density on hats 12 × 13
    sp = hole_filling_problem(; nL = 16, hole_a = 0.9, hole_b = 0.8)
    xs, ps, ws, ms = sp.x0, sp.p0, sp.w0, sp.m
    tgs = projection_targets(HatProjection(12, 13), xs, sp.v0, ws; sp.L, sp.V,
        density = (a, b) -> sp.density.(a, b))
    ℓx, ℓv = sp.L / 12, 2sp.V / 12
    kds = KickDriftParametrisation(ws; sp.L, sp.V, nξ = 8, nv = 10, nw = 4)
    nss = NonSeparableParametrisation(ws; sp.L, sp.V, nξ = 6, nv = 7, nℓ = 3)
    n = remap_parameter_count(kds)
    fit(par; kwargs...) = fit_remap_parameters(
        par, xs, ps, ws, tgs.x, tgs.v, tgs.w; m = ms,
        ℓx, ℓv, kwargs...)
    # the objective ½|r|² + ½β|θ|² with the residual recomputed from the stages
    function objective(par, θ, β)
        X, P = remap_stages(par, θ, xs, ps, ws; m = ms)
        r = vcat((X .- tgs.x) ./ ℓx, (P ./ (ms .* tgs.w) .- tgs.v) ./ ℓv)
        0.5 * dot(r, r) + 0.5β * dot(θ, θ)
    end
    # the gradient Jmᵀr + βθ of the objective, with the Jacobian Jm of the package's residual
    function gradient(par, θ, β)
        residual(θ) = begin
            X, P = remap_stages(par, θ, xs, ps, ws; m = ms)
            vcat((X .- tgs.x) ./ ℓx, (P ./ (ms .* tgs.w) .- tgs.v) ./ ℓv)
        end
        Jm = ParticleRemapping.Projection._residual_jacobian(
            par, residual, θ, xs, ps, ws, tgs.w, ms, ℓx, ℓv)
        Jm' * residual(θ) + β * θ
    end
    g0 = norm(gradient(kds, zeros(n), 0))
    code(fr) = string(fr.status.code)

    @testset "the default is Gauss–Newton" begin
        for par in (kds, nss)
            a = fit(par; maxit = 5)
            b = fit(par; maxit = 5, solver = :gauss_newton)
            @test a.θ == b.θ && a.F == b.F && a.iterations == b.iterations
            @test nameof(typeof(a.status)) === :SolverStatus && a.status !== nothing
            @test a.F < a.F0
            @test 1 <= a.iterations == a.status.iterations <= 5
            @test a.F ≈ objective(par, a.θ, 1e-3) rtol = 1e-12
        end
    end

    # the default solver, with the budget of the fit tests
    fgn = fit(kds; maxit = 20)
    @testset "Gauss–Newton with DogLeg" begin
        @test code(fgn) == "SUCCESS"
        @test fgn.status.stationary
        @test fgn.iterations isa Int
        @test fgn.iterations == fgn.status.iterations <= 20
        @test fgn.F < fgn.F0
        @test fgn.F ≈ objective(kds, fgn.θ, 1e-3) rtol = 1e-12
        @test fgn.F0 ≈ objective(kds, zeros(n), 1e-3) rtol = 1e-12
        # θ is stationary for the objective, not only for the solver's model of it
        @test norm(gradient(kds, fgn.θ, 1e-3)) < 1e-8 * g0
        # the bound of the budget is inclusive: the maxit guard runs before the length guard, so
        # typemax(Int32) passes it and the length guard answers, and no solve runs
        @test_throws DimensionMismatch fit(
            kds; solver = :gauss_newton, maxit = typemax(Int32), θ0 = zeros(n + 1))
    end

    @testset "the regularisation rows" begin
        f10 = fit(kds; solver = :gauss_newton, maxit = 20, β = 10)
        @test norm(f10.θ) < 0.8 * norm(fgn.θ)
        @test f10.F ≈ objective(kds, f10.θ, 10) rtol = 1e-12
        @test norm(gradient(kds, f10.θ, 10)) < 1e-8 * g0
    end

    @testset "the budget of Gauss–Newton" begin
        f3 = fit(nss; solver = :gauss_newton, maxit = 3)
        @test code(f3) == "MAXITERS"
        @test f3.iterations == 3
    end

    @testset "the start and the zero budget" begin
        θ0 = 1e-3 .* randn(Xoshiro(7), n)
        f0 = fit(kds; maxit = 0, θ0)
        fz = fit(kds; maxit = 0)
        @test f0.θ == θ0
        @test f0.θ !== θ0
        @test fz.θ == zeros(n)
        @test f0.iterations == 0 && fz.iterations == 0
        @test f0.F0 == fz.F0
        @test fz.F == fz.F0
        @test f0.F ≈ objective(kds, θ0, 1e-3) rtol = 1e-12
        @test f0.F != f0.F0
        @test code(f0) == "MAXITERS" && f0.status.iterations == 0
        # a Float32 start is converted to Float64
        f32 = fit(kds; maxit = 0, θ0 = Float32.(θ0))
        @test f32.θ isa Vector{Float64} && f32.θ == Float32.(θ0)
        # iterations is returned as an Int for every integer type of maxit
        for cap in (Int32(0), false, Int128(0))
            @test fit(kds; maxit = cap).iterations === 0
        end
        for cap in (Int32(1), true, Int128(1))
            @test fit(kds; maxit = cap).iterations === 1
        end
    end

    @testset "non-finite start objective" begin
        nss1 = NonSeparableParametrisation(ws; sp.L, sp.V, nξ = 6, nv = 7, nℓ = 3,
            iteration_limit = 1)
        θbad = randn(Xoshiro(11), remap_parameter_count(nss1))
        fbad = fit(nss1; maxit = 5, θ0 = θbad)
        @test fbad.θ == θbad && fbad.θ !== θbad
        @test isfinite(fbad.F0)
        @test isnan(fbad.F)
        @test fbad.iterations === 0
        @test code(fbad) == "NONFINITE"
    end

    @testset "trace output" begin
        out = mktemp() do _, io
            redirect_stdout(io) do
                fit(kds; solver = :gauss_newton, maxit = 0, trace = true)
                flush(stdout)
            end
            seekstart(io)
            read(io, String)
        end
        @test occursin("MAXITERS after 0 iterations: F/F0 = 1.0000", out)
    end

    @testset "the non-separable map with Gauss–Newton" begin
        fns = fit(nss; solver = :gauss_newton, maxit = 20)
        @test all(isfinite, fns.θ)
        @test isfinite(fns.F)
        @test fns.F < fns.F0
        @test fns.iterations <= 20
        @test code(fns) in ("SUCCESS", "MAXITERS", "LINESEARCH_FAILED")
        @test fns.F ≈ objective(nss, fns.θ, 1e-3) rtol = 1e-12
    end

    @testset "invalid arguments" begin
        # Gauss–Newton is the one solver; every other symbol is refused before any evaluation
        for solver in (:levenberg_marquardt, :newton)
            @test_throws Regex("^ArgumentError: unknown solver :$solver; use :gauss_newton\$") fit(
                kds; solver)
        end
        @test_throws Regex("^DimensionMismatch: θ0 has length $(n + 1), the " *
                           "parametrisation has $n parameters\$") fit(
            kds; θ0 = zeros(n + 1))
        @test_throws r"^ArgumentError: θ0 must be finite$" fit(
            kds; θ0 = [NaN; zeros(n - 1)])
        @test_throws r"^ArgumentError: θ0 must be finite$" fit(
            kds; θ0 = [zeros(n - 1); -Inf])
        @test_throws r"^ArgumentError: maxit must be non-negative, got -1$" fit(
            kds; maxit = -1)
        @test_throws r"^ArgumentError: β must be finite and non-negative, got -1.0$" fit(
            kds; β = -1.0)
        @test_throws r"^ArgumentError: β must be finite and non-negative, got NaN$" fit(
            kds; β = NaN)
        @test_throws r"^ArgumentError: β must be finite and non-negative, got Inf$" fit(
            kds; β = Inf)
        @test_throws r"^ArgumentError: β must be finite and non-negative, got -Inf$" fit(
            kds; β = -Inf)
        # β = 0 is allowed: the regularisation rows are zero rows
        @test fit(kds; maxit = 0, β = 0).iterations == 0
        @test_throws r"^ArgumentError: maxit must be at most 2147483647 for :gauss_newton, got 2147483648$" fit(
            kds; maxit = Int(typemax(Int32)) + 1)
        # the maxit guards run before the length guard
        @test_throws r"^ArgumentError: maxit must be at most 2147483647" fit(
            kds; maxit = Int(typemax(Int32)) + 1, θ0 = zeros(n + 1))
        @test_throws TypeError fit(kds; maxit = 3.0)
    end
end
