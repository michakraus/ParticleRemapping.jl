# The hole-filling experiment at both hole positions: the whole construction (projection,
# targets, fit, map, reweighting, charge-conserving reweighting) for the weight-dependent kicks and
# drifts, the non-separable stages and the controls A, B, C, each fitted by Gauss–Newton with
# DogLeg. Each number of the tables of the documentation is checked with the digits it is quoted
# with, a number quoted as a range within that range widened by its width, and the return code of
# each fit by its name; the values of each fitted row are printed. At x_h = π the fit of the kicks
# and drifts from θ = 0 ends on a mirror-symmetric saddle, and a fit from a seeded start fills the
# hole. Cases that must fail: the controls B and C (one area-preserving map of (x, v)) fitted from
# θ = 0 do not shrink the hole beyond the reweighting alone; the kicks and drifts at the saddle do
# not fill the hole at x_h = π; the target weights alone miss the grid charge; the sensitivities
# without the terms of the shift and the boost miss the derivatives.

using ForwardDiff
using LinearAlgebra
using ParticleRemapping
using Random
using Test

# One unit of the last digit of the number `quoted`.
function quoted_unit(quoted::AbstractString)
    m = match(r"^-?(\d+)(?:\.(\d+))?(?:e(-?\d+))?$", quoted)
    decimals = m.captures[2] === nothing ? 0 : length(m.captures[2])
    exponent = m.captures[3] === nothing ? 0 : parse(Int, m.captures[3])
    10.0^(exponent - decimals)
end

# `x` agrees with the number `quoted` to the digits that `quoted` shows: within half a unit of its
# last digit. A range `lo–hi` holds `x` once it is widened on each side by its width, and by at
# least one unit of the last digit of `hi`.
function agrees(x, quoted::AbstractString)
    if occursin('–', quoted)
        lo, hi = split(quoted, '–')
        d = max(parse(Float64, hi) - parse(Float64, lo), quoted_unit(hi))
        return parse(Float64, lo) - d <= x <= parse(Float64, hi) + d
    end
    abs(x - parse(Float64, quoted)) <= 0.5 * quoted_unit(quoted) * (1 + 1e-9)
end

# The stationarity ratio ‖JᵀF‖ / (‖J‖_F ‖F‖) of the fit of the variant `v` of the experiment `e` at
# its θ: `F = [r; √β θ]` is the residual of the fit with the regularisation, and
# `J = [∂r/∂θ; √β I]` its Jacobian, by forward mode.
function stationarity_ratio(e, v; β = 1e-3)
    prob, tg, hp, r = e.problem, e.targets, e.projection, e.results[v]
    ℓx, ℓv = prob.L / hp.nX, 2prob.V / (hp.nV - 1)
    function residual(θ)
        X, P = remap_stages(r.parametrisation, θ, prob.x0, prob.p0, prob.w0; prob.m)
        vcat((X .- tg.x) ./ ℓx, (P ./ (prob.m .* tg.w) .- tg.v) ./ ℓv)
    end
    Jr = ForwardDiff.jacobian(residual, r.θ)
    res = residual(r.θ)
    g = Jr' * res .+ β .* r.θ
    norm(g) /
    (sqrt(sum(abs2, Jr) + β * length(r.θ)) * sqrt(dot(res, res) + β * dot(r.θ, r.θ)))
end

# The return code of a fit by its name, with "stationary" where the stationarity test ended it.
stop(status) = string(status.code) * (status.stationary ? ", stationary" : "")

# Whether the test asserts the return code of the row of the variant `v` at `hole_x` on the Julia
# `version`. The fit of the kicks and drifts from θ = 0 at x_h = π ends on the saddle with
# `SUCCESS, stationary` on Julia 1.13.4 and 1.14.0-beta1; on a Julia with a prerelease tag the
# test asserts its numbers only, not its return code or its stationarity (KNOWN_ISSUES.md, K17).
function code_asserted(hole_x, v; version = VERSION)
    isempty(version.prerelease) || !(hole_x ≈ π && v === :weight_dependent)
end

# The experiment with `verbose = true`, and what it prints.
function experiment_and_output(; kwargs...)
    mktemp() do path, io
        e = redirect_stdout(() -> hole_filling_experiment(; verbose = true, kwargs...), io)
        close(io)
        e, read(path, String)
    end
end

# The line of the output `out` of the experiment for the variant `v`.
printed(out, v) = only(filter(l -> startswith(strip(l), "$v "), split(out, '\n')))

const VARIANTS = (:weight_dependent, :control_A, :control_B, :non_separable, :control_C)

# The tables of the documentation: largest empty disk near the hole, fill, objective; then with the
# charge-conserving weights the charge residual of w*, the empty disk, the fill and max |ln(w′/w*)|;
# and the return code of the fit.
const TABLES = Dict(
    π / 2 => (before = ("2.75", "0.00"), targets = ("1.58", "0.86"),
        rows = Dict(:θ0 => ("1.84", "0.51", "", "9.1e-3", "1.84", "0.51", "0.014", ""),
            :weight_dependent => ("1.43", "0.90", "0.258", "8.9e-2", "1.48", "0.78",
                "0.175–0.176", "MAXITERS"),
            :control_A =>
                ("1.67", "0.66", "0.834–0.835", "2.7e-2", "1.66", "0.70", "0.057",
                    "MAXITERS"),
            :control_B =>
                ("1.84", "0.62", "0.699", "2.8e-2", "1.83", "0.47", "0.061–0.062",
                    "MAXITERS"),
            :non_separable =>
                ("1.56", "0.90", "0.438", "9.1e-2", "1.57", "0.70", "0.166",
                    "MAXITERS"),
            :control_C => ("1.88", "0.54–0.55", "0.667–0.668", "3.2e-2", "1.85–1.86",
                "0.54–0.55", "0.071", "LINESEARCH_FAILED"))),
    Float64(π) => (before = ("2.75", "0.00"), targets = ("1.58", "0.86"),
        rows = Dict(:θ0 => ("1.88", "0.47", "", "4.2e-3", "1.91", "0.47", "0.008", ""),
            :weight_dependent => ("1.81", "0.47", "0.927", "4.2e-3", "1.84", "0.47",
                "0.008", "SUCCESS, stationary"),
            :control_A => ("1.88", "0.47", "0.999", "4.2e-3", "1.91", "0.47", "0.008",
                "SUCCESS, stationary"),
            :control_B => ("1.88", "0.55", "0.972", "4.2e-3", "1.86", "0.47", "0.008",
                "SUCCESS, stationary"),
            :non_separable =>
                ("1.51–1.52", "0.93–0.94", "0.520",
                    "6.4e-2", "1.55–1.56", "0.70", "0.169",
                    "LINESEARCH_FAILED"),
            :control_C => ("2.04", "0.54–0.55", "0.713–0.714", "1.8e-2", "2.05–2.06",
                "0.54–0.55", "0.053", "LINESEARCH_FAILED"))))

# The rows "seeded" of control A and control B at x_h = π: the range of each number over the
# seeds 101:105, and the return codes that the fits from these seeds end on.
const SEEDED_CONTROLS = Dict(
    :control_A => (
        ("1.49–1.61", "0.81–0.90", "0.654–0.711", "4.5e-2–5.2e-2", "1.50–1.74",
            "0.50–0.55", "0.102–0.116"),
        ("MAXITERS",)),
    :control_B => (
        ("1.51–1.98", "0.62–0.78", "0.749–0.937", "1.6e-2–2.1e-2", "1.45–1.96",
            "0.46–0.78", "0.047–0.071"),
        ("MAXITERS", "SUCCESS, stationary")))

@testset "the return codes that the test asserts" begin
    @test code_asserted(Float64(π), :weight_dependent; version = v"1.13.4")
    @test !code_asserted(Float64(π), :weight_dependent; version = v"1.15.0-DEV.1003")
    @test code_asserted(Float64(π), :control_A; version = v"1.15.0-DEV.1003")
    @test code_asserted(π / 2, :weight_dependent; version = v"1.15.0-DEV.1003")
end

@testset "Hole filling, x_h = $(round(hole_x; digits = 4))" for hole_x in (π / 2, Float64(π))
    Random.seed!(20261119)
    e, output = experiment_and_output(; hole_x)
    res = e.results
    near(v) = res[v].diagnostics.empty
    fill(v) = res[v].diagnostics.fill
    table = TABLES[hole_x]

    @testset "set-up" begin
        @test e.problem.N == 2280
        @test e.problem.nL^2 - e.problem.N == 24
        @test length(res[:weight_dependent].θ) == 2304
        @test length(res[:non_separable].θ) == 3264
    end

    @testset "momentum and total weight" begin
        for v in VARIANTS
            # relative; the source note quotes ≤ 1.1e-17, the manual below 1e-15
            @test res[v].ΔP < 1e-15
            @test res[v].Δw < 1e-12
        end
    end

    @testset "charge-conserving reweighting" begin
        for v in (:θ0, VARIANTS...)
            r = res[v]
            @test r.charge_ok && r.charge_res_c < 1e-12 && minimum(r.wc) > 0
            @test r.charge_its <= 16          # the note quotes at most 16 Newton iterations
            @test r.charge_res > 1e-3         # must fail: the target weights alone
        end
    end

    @testset "tables" begin
        @test agrees(e.before.empty, table.before[1]) &&
              agrees(e.before.fill, table.before[2])
        @test agrees(e.target_diagnostics.empty, table.targets[1])
        @test agrees(e.target_diagnostics.fill, table.targets[2])
        for (v, row) in table.rows
            r = res[v]
            v === :θ0 || println("hole table, x_h = $(round(hole_x; digits = 4)), $v: ",
                join(
                    string.((
                        r.diagnostics.empty, r.diagnostics.fill, r.ratio, r.charge_res,
                        r.diagnostics_charge.empty, r.diagnostics_charge.fill, r.dlnw)),
                    " "),
                ", ", stop(r.status))
            @test agrees(r.diagnostics.empty, row[1])
            @test agrees(r.diagnostics.fill, row[2])
            v === :θ0 || @test agrees(r.ratio, row[3])
            @test agrees(r.charge_res, row[4])
            @test agrees(r.diagnostics_charge.empty, row[5])
            @test agrees(r.diagnostics_charge.fill, row[6])
            @test agrees(r.dlnw, row[7])
            if v !== :θ0 && code_asserted(hole_x, v)
                @test stop(r.status) == row[8]
                # the line that `verbose = true` prints for the variant names its return code
                @test occursin(row[8], printed(output, v))
            end
        end
    end

    @testset "sensitivities of the non-separable map" begin
        Random.seed!(20261119)
        prob = e.problem
        r = res[:non_separable]
        @test r.stage_iterations == 50        # the iteration limit is active at the fitted θ
        err = nonseparable_sensitivity_error(
            r.parametrisation, r.θ, prob.x0, prob.p0, prob.w0)
        err_nc = nonseparable_sensitivity_error(r.parametrisation, r.θ, prob.x0, prob.p0,
            prob.w0; centroid_terms = false)
        # the source note quotes 2.8e-13 and ≥ 2.4e-2, the manual 2.3e-13 and ≥ 2.9e-2
        @test err < 1e-10
        @test err_nc > 2e-2                   # must fail
    end

    @testset "fill" begin
        # one area-preserving map of (x, v), fitted from θ = 0, does not shrink the hole beyond
        # the reweighting
        @test near(:control_B) > 0.95 * near(:θ0)
        @test near(:control_C) > 0.95 * near(:θ0)
        if hole_x ≈ π / 2
            @test near(:weight_dependent) < 0.9 * near(:θ0)
            @test fill(:weight_dependent) > fill(:control_B) + 0.2
            @test near(:non_separable) < 0.9 * near(:θ0)
            @test fill(:non_separable) > fill(:control_C) + 0.2
        else
            # the fit of the kicks and drifts from θ = 0 ends on the mirror-symmetric saddle of
            # the objective, where the ±v coupling keeps the hole open
            r = res[:weight_dependent]
            if code_asserted(hole_x, :weight_dependent)
                @test string(r.status.code) == "SUCCESS" && r.status.stationary
                @test stationarity_ratio(e, :weight_dependent) < 1e-8
            end
            @test near(:weight_dependent) > 0.9 * near(:θ0)
            # the non-separable stages do not have the ±v coupling
            @test near(:non_separable) < 0.9 * near(:θ0)
            @test fill(:non_separable) > fill(:control_C) + 0.2
            @test fill(:non_separable) > fill(:weight_dependent) + 0.2
        end
    end

    # A seeded start breaks the mirror symmetry, and the kicks and drifts fill the hole at
    # x_h = π. The bounds hold the range of the seeds 101:105 of the script with a margin.
    hole_x ≈ π && @testset "the seeded weight-dependent fit" begin
        es = hole_filling_experiment(;
            hole_x, variants = (:weight_dependent,), seed = 101, targets = e.targets)
        r = es.results[:weight_dependent]
        prob = es.problem
        @test es.targets === e.targets
        @test string(r.status.code) == "MAXITERS" && r.iterations == 40
        @test r.diagnostics.empty < 1.55
        @test r.diagnostics.empty < 0.9 * es.results[:θ0].diagnostics.empty
        @test r.diagnostics.fill > 0.85
        mirror = EllipticHole(hole_x, -1.2, 0.45, 0.40; prob.L)
        @test local_empty_disk(mirror, r.X, r.v1; prob.dx, prob.dv) < 1.0
        @test 0.25 < r.ratio < 0.40
    end

    # Control A and control B from the seeded start of the seed 101 at x_h = π: each number
    # within its quoted range, widened as `agrees` widens a range.
    hole_x ≈ π && @testset "the seeded controls" begin
        es = hole_filling_experiment(;
            hole_x, variants = (:control_A, :control_B), seed = 101, targets = e.targets)
        for (v, (quoted, stops)) in SEEDED_CONTROLS
            r = es.results[v]
            @test stop(r.status) in stops
            values = (r.diagnostics.empty, r.diagnostics.fill, r.ratio, r.charge_res,
                r.diagnostics_charge.empty, r.diagnostics_charge.fill, r.dlnw)
            for (x, q) in zip(values, quoted)
                @test agrees(x, q)
            end
        end
    end

    @testset "the remap result" begin
        prob = e.problem
        method = NonSeparableRemap(; charge = ChargeConservingReweighting(prob.L, 32),
            density = (xs, vs) -> prob.density.(xs, vs), solver = :gauss_newton)
        out = remap(method, prob.markers; targets = e.targets)
        r = res[:non_separable]
        @test out.info.X == r.X && out.info.P == r.P
        @test out.markers.x == mod.(r.X, prob.L)
        @test out.markers.w == r.wc
        @test out.markers.v == r.vc
        @test all(0 .<= out.markers.x .< prob.L)
        @test abs(total_momentum(out.markers) - total_momentum(prob.markers)) <
              1e-12 * sum(abs, prob.p0)
        @test total_weight(out.markers) ≈ total_weight(prob.markers) rtol = 1e-12
    end
end

# The note quotes the errors at x_h = π; at π/2 the bands [1, 2) and [3, 4) differ in the second
# digit.
@testset "Relative error of the projected target weights, x_h = $(round(hole_x; digits = 4))" for (hole_x,
    quoted) in ((Float64(π), ("0.10", "0.17", "0.13", "0.21")),
    (π / 2, ("0.10", "0.18", "0.13", "0.22")))
    errs = projected_weight_errors(hole_filling_problem(; hole_x))
    for (e, q) in zip(errs, quoted)
        @test agrees(e.rms, q)
    end
    @test agrees(maximum(e.max for e in errs), "0.96")
end
