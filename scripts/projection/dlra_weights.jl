# The target weights after a push: the markers of the hole-filling experiment (x_h = π) follow the
# characteristics with E = −0.15 sin x to T = 3; the target weights come from f(T), from the BUG
# solution at ϑ = 1e-3 … 1e-6 and at the fixed rank one, and from the lumped projection; each set
# enters the fit of the non-separable map and the charge-conserving reweighting. Prints the two
# tables of the documentation and checks every number against the value the note quotes.
#
#   julia --startup-file=no --project=scripts scripts/projection/dlra_weights.jl
#
# Runtime: about 3 minutes.

using ParticleRemapping
using Printf
using Random
using Test

function agrees(x, quoted::AbstractString)
    m = match(r"^-?(\d+)(?:\.(\d+))?(?:e(-?\d+))?$", quoted)
    decimals = m.captures[2] === nothing ? 0 : length(m.captures[2])
    exponent = m.captures[3] === nothing ? 0 : parse(Int, m.captures[3])
    abs(x - parse(Float64, quoted)) <= 0.5 * 10.0^(exponent - decimals) * (1 + 1e-9)
end

const SOURCES = (
    "f(T)", "ϑ = 1e-3", "ϑ = 1e-4", "ϑ = 1e-5", "ϑ = 1e-6", "rank 1", "projection")

const NOTE_WEIGHTS = Dict(
    "ϑ = 1e-3" => ("4", "1.1e-3", "3.0e-2", "7.5e-2", "1.2e-1", "2.0e-1"),
    "ϑ = 1e-4" => ("7", "8.5e-5", "3.1e-3", "1.3e-2", "4.1e-2", "1.2e-1"),
    "ϑ = 1e-5" => ("10", "4.2e-6", "2.7e-4", "1.3e-3", "5.5e-3", "2.0e-2"),
    "ϑ = 1e-6" => ("13", "1.9e-7", "2.3e-5", "8.7e-5", "5.8e-4", "2.9e-3"),
    "rank 1" => ("1", nothing, "1.5e-1", "1.8e-1", "1.4e-1", "1.8e-1"),
    "projection" => ("", "", "4.1e-2", "1.1e-1", "8.4e-2", "2.0e-1"))

const NOTE_REMAP = Dict(
    "f(T)" => ("1.05", "0.82", "0.113", "1.06", "0.74", "2.4e-2", "0.036"),
    "ϑ = 1e-3" => ("1.30", "1.33", "0.364", "1.28", "1.25", "1.8e-2", "0.021"),
    "ϑ = 1e-4" => ("1.46", "1.05", "0.369", "1.42", "1.09", "1.3e-2", "0.021"),
    "ϑ = 1e-5" => ("1.05", "0.82", "0.114", "1.06", "0.74", "2.4e-2", "0.037"),
    "ϑ = 1e-6" => ("1.05", "0.82", "0.113", "1.06", "0.74", "2.4e-2", "0.036"),
    "rank 1" => ("1.10", "0.97", "0.330", "3.33", "0.86", "1.7e-1", "0.271"),
    "projection" => ("2.67", "0.78", "0.510", "2.74", "0.74", "2.1e-2", "0.043"))

Random.seed!(20261119)
e = dlra_weights_experiment()
res = e.results
function rank_of(name)
    name == "rank 1" ? lowrank_rank(e.solutions[name].f) :
    maximum(e.solutions[name].ranks)
end

@printf("push: E = −0.15 sin x, T = %.1f; max |v| = %.2f; pushed markers: empty disk %.2f, fill %.2f, density %.3f\n",
    e.T, maximum(abs, e.vT), e.pushed.largest_empty, e.pushed.fill, e.pushed.discrepancy)
for name in SOURCES[2:6]
    s = e.solutions[name]
    @printf("  BUG %-9s  rank %2d, mass error %.1e, %.2f s\n", name, rank_of(name),
        s.mass_error, s.time)
end
println("\ntarget weights: rank, mass error, rms of w*/(c f(T)) − 1 by |v| band")
for name in SOURCES[2:end]
    r = res[name]
    s = get(e.solutions, name, nothing)
    @printf("  %-11s %4s %9s   %s;  %d not positive\n", name,
        s === nothing ? "" : string(rank_of(name)),
        s === nothing ? "" : @sprintf("%.1e", s.mass_error),
        join([@sprintf("%.1e", b.rms) for b in r.errs], "  "), r.nneg)
end
println("\nremap: with w* (empty disk, fill, density) | charge-conserving (empty disk, fill, charge, max |ln|)")
for name in SOURCES
    r = res[name]
    @printf("  %-11s %5.2f %5.2f %6.3f | %5.2f %5.2f %8.1e %6.3f   (fit %.0f s, |ΔP|/Σ|p| %.1e, charge residual %.1e)\n",
        name, r.diagnostics.largest_empty, r.diagnostics.fill, r.diagnostics.discrepancy,
        r.diagnostics_charge.largest_empty, r.diagnostics_charge.fill, r.charge_res, r.dlnw,
        r.fit_time, r.ΔP, r.charge_res_c)
end

@testset "the note's numbers" begin
    @test agrees(maximum(abs, e.vT), "3.95")
    @test agrees(e.pushed.largest_empty, "1.38") && e.pushed.fill == 0
    @test agrees(e.pushed.discrepancy, "0.099")
    for (name, row) in NOTE_WEIGHTS
        isempty(row[1]) || @test rank_of(name) == parse(Int, row[1])
        # the mass error of the rank one is at the level of rounding: the note quotes its bound
        if row[2] === nothing
            @test e.solutions[name].mass_error < 1e-12
        elseif !isempty(row[2])
            @test agrees(e.solutions[name].mass_error, row[2])
        end
        for (b, q) in zip(res[name].errs, row[3:end])
            @test agrees(b.rms, q)
        end
    end
    for (name, row) in NOTE_REMAP
        r = res[name]
        values = (
            r.diagnostics.largest_empty, r.diagnostics.fill, r.diagnostics.discrepancy,
            r.diagnostics_charge.largest_empty, r.diagnostics_charge.fill, r.charge_res, r.dlnw)
        for (x, q) in zip(values, row)
            @test agrees(x, q)
        end
        @test r.nneg == 0 && r.ΔP < 1e-12 && r.charge_res_c < 1e-12
    end
end
