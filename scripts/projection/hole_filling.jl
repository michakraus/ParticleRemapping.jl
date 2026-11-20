# The hole-filling experiment of the Hamiltonian particle projection at both hole positions,
# x_h = π/2 (∂_x W ≠ 0) and x_h = π (∂_x W = 0): prints the two tables of the documentation (the
# remap with the target weights, and the charge-conserving reweighting), the cost of each map, the
# sensitivities of the non-separable map, and the relative error of the projected target weights,
# and checks every number against the value the note quotes.
#
#   julia --startup-file=no --project=scripts scripts/projection/hole_filling.jl
#
# Runtime: about 3 minutes.

using ParticleRemapping
using Printf
using Random
using Test

# `x` agrees with `quoted` to the digits that `quoted` shows.
function agrees(x, quoted::AbstractString)
    m = match(r"^-?(\d+)(?:\.(\d+))?(?:e(-?\d+))?$", quoted)
    decimals = m.captures[2] === nothing ? 0 : length(m.captures[2])
    exponent = m.captures[3] === nothing ? 0 : parse(Int, m.captures[3])
    abs(x - parse(Float64, quoted)) <= 0.5 * 10.0^(exponent - decimals) * (1 + 1e-9)
end

const ROWS = (:θ0 => "θ = 0, reweighting only", :weight_dependent => "weight-dependent",
    :control_A => "control A", :control_B => "control B", :non_separable => "non-separable",
    :control_C => "control C")

# The numbers of the note, row by row: empty disk, fill, objective; charge, disk, fill, max |ln|.
const NOTE = Dict(
    π / 2 => Dict(:θ0 => ("1.84", "0.51", "", "9.1e-3", "1.84", "0.51", "0.014"),
        :weight_dependent =>
            ("1.58", "0.94", "0.264", "9.0e-2", "1.56", "0.70", "0.173"),
        :control_A => ("1.65", "0.82", "0.795", "5.4e-2", "1.90", "0.78", "0.142"),
        :control_B => ("1.84", "0.62", "0.699", "2.8e-2", "1.83", "0.47", "0.061"),
        :non_separable => ("1.54", "0.94", "0.439", "9.2e-2", "1.56", "0.70", "0.165"),
        :control_C => ("1.80", "0.55", "0.663", "3.1e-2", "1.83", "0.55", "0.069")),
    Float64(π) => Dict(:θ0 => ("1.88", "0.47", "", "4.2e-3", "1.91", "0.47", "0.008"),
        :weight_dependent =>
            ("1.81", "0.47", "0.927", "4.2e-3", "1.84", "0.47", "0.008"),
        :control_A => ("1.88", "0.47", "0.999", "4.2e-3", "1.91", "0.47", "0.008"),
        :control_B => ("1.88", "0.55", "0.972", "4.2e-3", "1.86", "0.47", "0.008"),
        :non_separable => ("1.50", "0.94", "0.521", "6.3e-2", "1.53", "0.70", "0.165"),
        :control_C => ("2.00", "0.55", "0.716", "1.6e-2", "2.03", "0.55", "0.050")))

function report(hole_x)
    Random.seed!(20261119)
    e = hole_filling_experiment(; hole_x)
    prob, res = e.problem, e.results
    @printf("\nhole at x_h = %.4f: N = %d (%d markers removed); a perfect lattice has the empty disk %.2f\n",
        hole_x, prob.N, prob.nL^2 - prob.N, sqrt(0.5))
    @printf("%-26s %10s %6s %9s | %8s %6s %6s %8s\n",
        "", "empty disk", "fill", "objective",
        "charge", "disk", "fill", "max|ln|")
    @printf("%-26s %10.2f %6.2f\n", "before", e.before.empty, e.before.fill)
    @printf("%-26s %10.2f %6.2f\n", "targets", e.target_diagnostics.empty,
        e.target_diagnostics.fill)
    for (v, name) in ROWS
        r = res[v]
        @printf("%-26s %10.2f %6.2f %9s | %8.1e %6.2f %6.2f %8.3f\n", name,
            r.diagnostics.empty, r.diagnostics.fill,
            v === :θ0 ? "1" : @sprintf("%.3f", r.ratio), r.charge_res,
            r.diagnostics_charge.empty, r.diagnostics_charge.fill, r.dlnw)
    end
    for (v, name) in ROWS[2:end]
        r = res[v]
        @printf("  %-17s %4d parameters, fit %2d iterations in %4.1f s, one map %.3f s%s; |ΔP|/Σ|p| %.1e, Newton %d iterations, charge residual %.1e\n",
            name, length(r.θ), r.iterations, r.fit_time, r.map_time,
            r.stage_iterations === nothing ? "" :
            @sprintf(", %d fixed-point iterations", r.stage_iterations),
            r.ΔP, r.charge_its, r.charge_res_c)
    end
    Random.seed!(20261119)
    r = res[:non_separable]
    err = nonseparable_sensitivity_error(r.parametrisation, r.θ, prob.x0, prob.p0, prob.w0)
    err_nc = nonseparable_sensitivity_error(
        r.parametrisation, r.θ, prob.x0, prob.p0, prob.w0;
        centroid_terms = false)
    @printf("  non-separable sensitivities against forward mode: %.1e; without the shift and the boost %.1e\n",
        err, err_nc)
    errs = projected_weight_errors(prob)
    @printf("  projected target weights, rms of w*_h/(c f) − 1 by |v| band: %s; max %.2f\n",
        join([@sprintf("%.3f", b.rms) for b in errs], ", "), maximum(b.max for b in errs))

    @testset "the note's numbers, x_h = $(round(hole_x; digits = 4))" begin
        @test agrees(e.before.empty, "2.75") && agrees(e.before.fill, "0.00")
        @test agrees(e.target_diagnostics.empty, "1.58")
        @test agrees(e.target_diagnostics.fill, "0.86")
        for (v, row) in NOTE[hole_x]
            r = res[v]
            values = (r.diagnostics.empty, r.diagnostics.fill, v === :θ0 ? NaN : r.ratio,
                r.charge_res, r.diagnostics_charge.empty, r.diagnostics_charge.fill, r.dlnw)
            for (x, q) in zip(values, row)
                isempty(q) || @test agrees(x, q)
            end
            @test r.charge_its <= 16
            v === :θ0 || @test r.ΔP < 1e-12
        end
        @test err < 1e-10 && err_nc > 2e-2
    end
    e
end

report(π / 2)
report(Float64(π))
