# The hole-filling experiment at both hole positions: the whole construction (projection,
# targets, fit, map, reweighting, charge-conserving reweighting) for the weight-dependent kicks and
# drifts, the non-separable stages and the controls A, B, C. Each number of the tables of the
# documentation is checked with the digits it is quoted with. Cases that must fail: the controls
# B and C (one area-preserving map of (x, v)) do not shrink the hole beyond the reweighting alone;
# the kicks and drifts do not fill the hole at x_h = π; the target weights alone miss the grid
# charge; the sensitivities without the terms of the shift and the boost miss the derivatives.

using ParticleRemapping
using Random
using Test

# `x` agrees with the number `quoted` to the digits that `quoted` shows: within half a unit of its
# last digit.
function agrees(x, quoted::AbstractString)
    m = match(r"^-?(\d+)(?:\.(\d+))?(?:e(-?\d+))?$", quoted)
    decimals = m.captures[2] === nothing ? 0 : length(m.captures[2])
    exponent = m.captures[3] === nothing ? 0 : parse(Int, m.captures[3])
    abs(x - parse(Float64, quoted)) <= 0.5 * 10.0^(exponent - decimals) * (1 + 1e-9)
end

const VARIANTS = (:weight_dependent, :control_A, :control_B, :non_separable, :control_C)

# The tables of the documentation: largest empty disk near the hole, fill, objective; then with the
# charge-conserving weights the charge residual of w*, the empty disk, the fill and max |ln(w′/w*)|.
const TABLES = Dict(
    π / 2 => (before = ("2.75", "0.00"), targets = ("1.58", "0.86"),
        rows = Dict(:θ0 => ("1.84", "0.51", "", "9.1e-3", "1.84", "0.51", "0.014"),
            :weight_dependent =>
                ("1.58", "0.94", "0.264", "9.0e-2", "1.56", "0.70", "0.173"),
            :control_A => ("1.65", "0.82", "0.795", "5.4e-2", "1.90", "0.78", "0.142"),
            :control_B => ("1.84", "0.62", "0.699", "2.8e-2", "1.83", "0.47", "0.061"),
            :non_separable =>
                ("1.54", "0.94", "0.439", "9.2e-2", "1.56", "0.70", "0.165"),
            :control_C => ("1.80", "0.55", "0.663", "3.1e-2", "1.83", "0.55", "0.069"))),
    Float64(π) => (before = ("2.75", "0.00"), targets = ("1.58", "0.86"),
        rows = Dict(:θ0 => ("1.88", "0.47", "", "4.2e-3", "1.91", "0.47", "0.008"),
            :weight_dependent =>
                ("1.81", "0.47", "0.927", "4.2e-3", "1.84", "0.47", "0.008"),
            :control_A => ("1.88", "0.47", "0.999", "4.2e-3", "1.91", "0.47", "0.008"),
            :control_B => ("1.88", "0.55", "0.972", "4.2e-3", "1.86", "0.47", "0.008"),
            :non_separable =>
                ("1.50", "0.94", "0.521", "6.3e-2", "1.53", "0.70", "0.165"),
            :control_C => ("2.00", "0.55", "0.716", "1.6e-2", "2.03", "0.55", "0.050"))))

@testset "Hole filling, x_h = $(round(hole_x; digits = 4))" for hole_x in (π / 2, Float64(π))
    Random.seed!(20261119)
    e = hole_filling_experiment(; hole_x)
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
            @test res[v].ΔP < 1e-12            # the note quotes ≤ 1.1e-17 relative
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
            @test agrees(r.diagnostics.empty, row[1])
            @test agrees(r.diagnostics.fill, row[2])
            v === :θ0 || @test agrees(r.ratio, row[3])
            @test agrees(r.charge_res, row[4])
            @test agrees(r.diagnostics_charge.empty, row[5])
            @test agrees(r.diagnostics_charge.fill, row[6])
            @test agrees(r.dlnw, row[7])
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
        @test err < 1e-10                     # the note quotes 2.8e-13
        @test err_nc > 2e-2                   # must fail; the note quotes ≥ 2.4e-2
    end

    @testset "fill" begin
        # one area-preserving map of (x, v) does not shrink the hole beyond the reweighting
        @test near(:control_B) > 0.95 * near(:θ0)
        @test near(:control_C) > 0.95 * near(:θ0)
        if hole_x ≈ π / 2
            @test near(:weight_dependent) < 0.9 * near(:θ0)
            @test fill(:weight_dependent) > fill(:control_B) + 0.2
            @test near(:non_separable) < 0.9 * near(:θ0)
            @test fill(:non_separable) > fill(:control_C) + 0.2
        else
            # the limit of the separable map where ∂_x W = 0
            @test near(:weight_dependent) > 0.9 * near(:θ0)
            # the non-separable stages do not have the ±v coupling
            @test near(:non_separable) < 0.9 * near(:θ0)
            @test fill(:non_separable) > fill(:control_C) + 0.2
            @test fill(:non_separable) > fill(:weight_dependent) + 0.2
        end
    end

    @testset "the remap result" begin
        prob = e.problem
        method = NonSeparableRemap(; charge = ChargeConservingReweighting(prob.L, 32),
            density = (xs, vs) -> prob.density.(xs, vs))
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
