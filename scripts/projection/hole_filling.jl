# The hole-filling experiment of the Hamiltonian particle projection at both hole positions,
# x_h = π/2 (∂_x W ≠ 0) and x_h = π (∂_x W = 0), each fit by Gauss–Newton with DogLeg: prints the
# two tables of the documentation (the remap with the target weights, and the charge-conserving
# reweighting) with the return code of each fit, the cost of each map, the sensitivities of the
# non-separable map, and the relative error of the projected target weights, and checks every
# number against the value the documentation quotes, a number quoted as a range between its ends.
# It then runs the experiment again with the targets moved by 1e-13 relative and checks every
# quoted number there too, and runs the weight-dependent fit, control A and control B at x_h = π
# from the seeded starts 101:105 and checks that each of their numbers lies in its quoted range.
#
#   julia --startup-file=no --project=scripts scripts/projection/hole_filling.jl
#
# Runtime: about 31 minutes.

using ParticleRemapping
using Printf
using Random
using Test

# `x` agrees with `quoted` to the digits that `quoted` shows; a range `lo–hi` holds `x` between
# its ends.
function agrees(x, quoted::AbstractString)
    if occursin('–', quoted)
        lo, hi = parse.(Float64, split(quoted, '–'))
        return lo <= x <= hi
    end
    m = match(r"^-?(\d+)(?:\.(\d+))?(?:e(-?\d+))?$", quoted)
    decimals = m.captures[2] === nothing ? 0 : length(m.captures[2])
    exponent = m.captures[3] === nothing ? 0 : parse(Int, m.captures[3])
    abs(x - parse(Float64, quoted)) <= 0.5 * 10.0^(exponent - decimals) * (1 + 1e-9)
end

# The return code of a fit by its name, with "stationary" where the stationarity test ended it.
stop(status) = string(status.code) * (status.stationary ? ", stationary" : "")

const ROWS = (:θ0 => "θ = 0, reweighting only", :weight_dependent => "weight-dependent",
    :control_A => "control A", :control_B => "control B", :non_separable => "non-separable",
    :control_C => "control C")

# The numbers of the documentation, row by row: empty disk, fill, objective; charge, disk, fill,
# max |ln|; the return code of the fit.
const NOTE = Dict(
    π / 2 => Dict(:θ0 => ("1.84", "0.51", "", "9.1e-3", "1.84", "0.51", "0.014", ""),
        :weight_dependent =>
            ("1.43", "0.90", "0.258", "8.9e-2", "1.48", "0.78", "0.175–0.176",
                "MAXITERS"),
        :control_A =>
            ("1.67", "0.66", "0.834–0.835", "2.7e-2", "1.66", "0.70", "0.057",
                "MAXITERS"),
        :control_B =>
            ("1.84", "0.62", "0.699", "2.8e-2", "1.83", "0.47", "0.061–0.062",
                "MAXITERS"),
        :non_separable => ("1.56", "0.90", "0.438", "9.1e-2", "1.57", "0.70", "0.166",
            "MAXITERS"),
        :control_C => ("1.88", "0.54–0.55", "0.667–0.668", "3.2e-2", "1.85–1.86",
            "0.54–0.55", "0.071", "LINESEARCH_FAILED")),
    Float64(π) => Dict(:θ0 => ("1.88", "0.47", "", "4.2e-3", "1.91", "0.47", "0.008", ""),
        :weight_dependent =>
            ("1.81", "0.47", "0.927", "4.2e-3", "1.84", "0.47", "0.008",
                "SUCCESS, stationary"),
        :control_A => ("1.88", "0.47", "0.999", "4.2e-3", "1.91", "0.47", "0.008",
            "SUCCESS, stationary"),
        :control_B => ("1.88", "0.55", "0.972", "4.2e-3", "1.86", "0.47", "0.008",
            "SUCCESS, stationary"),
        :non_separable =>
            ("1.51–1.52", "0.93–0.94", "0.520", "6.4e-2", "1.55–1.56", "0.70",
                "0.169", "LINESEARCH_FAILED"),
        :control_C => ("2.04", "0.54–0.55", "0.713–0.714", "1.8e-2", "2.05–2.06",
            "0.54–0.55", "0.053", "LINESEARCH_FAILED")))

# The weight-dependent fit at x_h = π from the seeded starts 101:105: the range over the seeds of
# the empty disk, the fill and the objective, of the charge, the disk, the fill and max |ln| with
# the charge-conserving weights, and of the empty disk of the mirror ellipse at v = −1.2.
const SEEDS = 101:105
const SEEDED = (empty = "1.35–1.45", fill = "0.97–1.06", ratio = "0.285–0.333",
    charge = "5.6e-2–6.2e-2", cempty = "1.32–1.41", cfill = "0.74–0.78",
    dlnw = "0.153–0.168", mirror = "0.76–0.87")

# Control A and control B at x_h = π from the same seeded starts: the range over the seeds of the
# numbers of the tables, and the return codes that the fits end on.
const SEEDED_CONTROLS = Dict(
    :control_A => (
        ("1.49–1.61", "0.81–0.90", "0.654–0.711", "4.5e-2–5.2e-2", "1.50–1.74",
            "0.50–0.55", "0.102–0.116"),
        ("MAXITERS",)),
    :control_B => (
        ("1.51–1.98", "0.62–0.78", "0.749–0.937", "1.6e-2–2.1e-2", "1.45–1.96",
            "0.46–0.78", "0.047–0.071"),
        ("MAXITERS", "SUCCESS, stationary")))

# The values of the documentation's tables of a variant's entry `r` of `results`.
function row_values(v, r)
    (r.diagnostics.empty, r.diagnostics.fill, v === :θ0 ? NaN : r.ratio,
        r.charge_res, r.diagnostics_charge.empty, r.diagnostics_charge.fill, r.dlnw)
end

# Each number of the documentation agrees with the experiment `e`.
function check_note(e, hole_x)
    @test agrees(e.before.empty, "2.75") && agrees(e.before.fill, "0.00")
    @test agrees(e.target_diagnostics.empty, "1.58")
    @test agrees(e.target_diagnostics.fill, "0.86")
    for (v, row) in NOTE[hole_x]
        r = e.results[v]
        for (x, q) in zip(row_values(v, r), row[1:7])
            isempty(q) || @test agrees(x, q)
        end
        @test r.charge_its <= 16
        v === :θ0 || @test r.ΔP < 1e-15
        v === :θ0 || @test stop(r.status) == row[8]
    end
end

function report(hole_x)
    Random.seed!(20261119)
    e = hole_filling_experiment(; hole_x)
    prob, res = e.problem, e.results
    @printf("\nhole at x_h = %.4f: N = %d (%d markers removed); a perfect lattice has the empty disk %.2f\n",
        hole_x, prob.N, prob.nL^2 - prob.N, sqrt(0.5))
    @printf("%-26s %10s %6s %9s %-20s | %8s %6s %6s %8s\n",
        "", "empty disk", "fill", "objective", " stop",
        "charge", "disk", "fill", "max|ln|")
    @printf("%-26s %10.2f %6.2f\n", "before", e.before.empty, e.before.fill)
    @printf("%-26s %10.2f %6.2f\n", "targets", e.target_diagnostics.empty,
        e.target_diagnostics.fill)
    for (v, name) in ROWS
        r = res[v]
        @printf("%-26s %10.2f %6.2f %9s %-20s | %8.1e %6.2f %6.2f %8.3f\n", name,
            r.diagnostics.empty, r.diagnostics.fill,
            v === :θ0 ? "1" : @sprintf("%.3f", r.ratio), v === :θ0 ? "" :
                                                         " " * stop(r.status),
            r.charge_res, r.diagnostics_charge.empty, r.diagnostics_charge.fill, r.dlnw)
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

    @testset "the numbers of the documentation, x_h = $(round(hole_x; digits = 4))" begin
        check_note(e, hole_x)
        @test err < 1e-10 && err_nc > 2e-2
    end
    e
end

# The experiment again, with the targets of `e` moved by 1e-13 relative: every quoted number
# agrees there too, and every fit starts from targets other than those of `e`.
function perturbed(e, hole_x)
    tg = e.targets
    x = tg.x .* (1 .+ 1e-13 .* randn(Xoshiro(1), length(tg.x)))
    Random.seed!(20261119)
    ep = hole_filling_experiment(; hole_x, targets = merge(tg, (; x)))
    @printf("\nx_h = %.4f, targets moved by 1e-13 relative:\n", hole_x)
    for (v, name) in ROWS[2:end]
        r = ep.results[v]
        @printf("%-26s %10.2f %6.2f %9.3f %-20s | %8.1e %6.2f %6.2f %8.3f\n", name,
            row_values(v, r)[1:3]..., " " * stop(r.status), row_values(v, r)[4:7]...)
    end
    @testset "the numbers of the documentation, moved targets, x_h = $(round(hole_x; digits = 4))" begin
        @test ep.targets.x != tg.x
        for (v, _) in ROWS[2:end]
            @test ep.results[v].θ != e.results[v].θ
        end
        check_note(ep, hole_x)
    end
end

# The weight-dependent fit, control A and control B at x_h = π from the seeded starts: each number
# in its quoted range.
function seeded(e)
    prob = e.problem
    mirror = EllipticHole(π, -1.2, 0.45, 0.40; prob.L)
    @printf("\nx_h = π, from θ0 = 1e-3 randn(Xoshiro(seed)):\n")
    @printf("%-6s %-17s %10s %6s %9s %-20s | %8s %6s %6s %8s | %6s\n",
        "seed", "", "empty disk", "fill",
        "objective", " stop", "charge", "disk", "fill", "max|ln|", "mirror")
    @testset "the seeded fits, x_h = π" begin
        for seed in SEEDS
            es = hole_filling_experiment(; hole_x = π,
                variants = (:weight_dependent, :control_A, :control_B), seed,
                targets = e.targets)
            for v in (:weight_dependent, :control_A, :control_B)
                r = es.results[v]
                values = row_values(v, r)
                md = local_empty_disk(mirror, r.X, r.v1; prob.dx, prob.dv)
                @printf("%-6d %-17s %10.2f %6.2f %9.3f %-20s | %8.1e %6.2f %6.2f %8.3f | %6.2f\n",
                    seed, v, values[1:3]..., " " * stop(r.status), values[4:7]..., md)
                @test r.ΔP < 1e-15
                if v === :weight_dependent
                    @test stop(r.status) == "MAXITERS"
                    for (x, q) in zip((values..., md), SEEDED)
                        @test agrees(x, q)
                    end
                else
                    quoted, stops = SEEDED_CONTROLS[v]
                    @test stop(r.status) in stops
                    for (x, q) in zip(values, quoted)
                        @test agrees(x, q)
                    end
                end
            end
        end
    end
end

e2 = report(π / 2)
e1 = report(Float64(π))
perturbed(e2, π / 2)
perturbed(e1, Float64(π))
seeded(e1)
