# The target weights after a push: the markers of the hole-filling experiment (x_h = π) pushed with
# E = −0.15 sin x to T = 3; target weights from f(T), from the BUG solution at four tolerances and at
# the fixed rank one, and from the lumped projection; each followed by the fit of the non-separable
# map and the charge-conserving reweighting. Each number of the tables of the documentation is
# checked with the digits or the range it is quoted with, and each fit by its return code. Case
# that must fail: the fixed rank one misses the sheared density.

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

# The return code of a fit by its name, with "stationary" where the stationarity test ended it.
stop(status) = string(status.code) * (status.stationary ? ", stationary" : "")

# rank, mass error, rms of w*/(c f(T)) − 1 in the bands |v| < 1, [1, 2), [2, 3), [3, 4); the
# mass error of the rank one is at the level of rounding, and `nothing` checks it by its bound
const WEIGHT_TABLE = Dict(
    "ϑ = 1e-3" => (4, "1.1e-3", ("3.0e-2", "7.5e-2", "1.2e-1", "2.0e-1")),
    "ϑ = 1e-4" => (7, "8.5e-5", ("3.1e-3", "1.3e-2", "4.1e-2", "1.2e-1")),
    "ϑ = 1e-5" => (10, "4.2e-6", ("2.7e-4", "1.3e-3", "5.5e-3", "2.0e-2")),
    "ϑ = 1e-6" => (13, "1.9e-7", ("2.3e-5", "8.7e-5", "5.8e-4", "2.9e-3")),
    "rank 1" => (1, nothing, ("1.5e-1", "1.8e-1", "1.4e-1", "1.8e-1")),
    "projection" => (nothing, nothing, ("4.1e-2", "1.1e-1", "8.4e-2", "2.0e-1")))

# with w*: empty disk, fill, density; charge-conserving: empty disk, fill, charge, max |ln|; the
# return code of the fit. No fit ends on the stationarity test, so a cell whose digits move when
# the targets move by 1e-13 relative, or on Julia 1.13 and 1.14, or that lies within a tenth of a
# unit of a rounding boundary, is a range over those runs.
const REMAP_TABLE = Dict(
    "f(T)" =>
        ("1.05", "0.82", "0.113", "1.066", "0.74", "2.4e-2", "0.035–0.036", "MAXITERS"),
    "ϑ = 1e-3" => ("1.43", "1.13", "0.308–0.309", "1.52–1.53", "1.21", "3.3e-2", "0.039",
        "LINESEARCH_FAILED"),
    "ϑ = 1e-4" => ("1.23", "1.01", "0.317", "1.27", "0.936", "1.7e-2", "0.027",
        "LINESEARCH_FAILED"),
    "ϑ = 1e-5" => ("1.05", "0.86", "0.114", "1.064", "0.74", "2.4e-2", "0.032", "MAXITERS"),
    "ϑ = 1e-6" => ("1.05", "0.82", "0.113", "1.066", "0.78", "2.5e-2", "0.035",
        "LINESEARCH_FAILED"),
    "rank 1" => ("1.06", "0.975", "0.332", "2.74", "0.90", "1.75e-1", "0.258",
        "LINESEARCH_FAILED"),
    "projection" => ("2.96–2.97", "0.35", "0.519", "3.05", "0.31", "7.9e-2", "0.0986",
        "LINESEARCH_FAILED"))

Random.seed!(20261119)
e = dlra_weights_experiment()
res = e.results

@testset "Target weights from a dynamical low-rank density" begin
    @testset "push" begin
        @test agrees(maximum(abs, e.vT), "3.95")
        @test agrees(e.pushed.largest_empty, "1.38")
        @test e.pushed.fill == 0
        @test agrees(e.pushed.discrepancy, "0.099")
    end

    @testset "every source" begin
        for (name, r) in res
            @test r.nneg == 0                  # no target weight is negative
            @test r.ΔP < 1e-12                 # the manual: below 1e-12, relative
            @test r.charge_ok && r.charge_res_c < 1e-12 && minimum(r.wc) > 0
        end
        @test res["f(T)"].errs[1].max < 1e-12
    end

    @testset "weight errors" begin
        for (name, (rank, mass, bands)) in WEIGHT_TABLE
            if rank !== nothing
                s = e.solutions[name]
                @test (name == "rank 1" ? ParticleRemapping.lowrank_rank(s.f) :
                       maximum(s.ranks)) == rank
                if mass === nothing
                    @test s.mass_error < 1e-12     # the manual: below 1e-12
                else
                    @test agrees(s.mass_error, mass)
                end
            end
            for (b, q) in zip(res[name].errs, bands)
                @test agrees(b.rms, q)
            end
        end
        # the error falls with the tolerance in every band
        for b in 1:4
            @test issorted(
                [res[n].errs[b].rms
                 for n in ("ϑ = 1e-3", "ϑ = 1e-4", "ϑ = 1e-5",
                    "ϑ = 1e-6")];
                rev = true)
        end
        @test maximum(e.solutions["ϑ = 1e-6"].ranks) < length(e.grid.x) ÷ 2
        # must fail: the fixed rank one misses the sheared density
        @test res["rank 1"].errs[1].rms > 0.05
    end

    @testset "remap" begin
        for (name, row) in REMAP_TABLE
            r = res[name]
            @test agrees(r.diagnostics.largest_empty, row[1])
            @test agrees(r.diagnostics.fill, row[2])
            @test agrees(r.diagnostics.discrepancy, row[3])
            @test agrees(r.diagnostics_charge.largest_empty, row[4])
            @test agrees(r.diagnostics_charge.fill, row[5])
            @test agrees(r.charge_res, row[6])
            @test agrees(r.dlnw, row[7])
            @test stop(r.status) == row[8]
        end
    end
end
