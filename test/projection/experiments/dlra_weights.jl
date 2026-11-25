# The target weights after a push: the markers of the hole-filling experiment (x_h = π) pushed with
# E = −0.15 sin x to T = 3; target weights from f(T), from the BUG solution at four tolerances and at
# the fixed rank one, and from the lumped projection; each followed by the fit of the non-separable
# map and the charge-conserving reweighting. Each number of the tables of the documentation is
# checked with the digits it is quoted with. Case that must fail: the fixed rank one misses the
# sheared density.

using ParticleRemapping
using Random
using Test

function agrees(x, quoted::AbstractString)
    m = match(r"^-?(\d+)(?:\.(\d+))?(?:e(-?\d+))?$", quoted)
    decimals = m.captures[2] === nothing ? 0 : length(m.captures[2])
    exponent = m.captures[3] === nothing ? 0 : parse(Int, m.captures[3])
    abs(x - parse(Float64, quoted)) <= 0.5 * 10.0^(exponent - decimals) * (1 + 1e-9)
end

# rank, mass error, rms of w*/(c f(T)) − 1 in the bands |v| < 1, [1, 2), [2, 3), [3, 4); the
# mass error of the rank one is at the level of rounding, and `nothing` checks it by its bound
const WEIGHT_TABLE = Dict(
    "ϑ = 1e-3" => (4, "1.1e-3", ("3.0e-2", "7.5e-2", "1.2e-1", "2.0e-1")),
    "ϑ = 1e-4" => (7, "8.5e-5", ("3.1e-3", "1.3e-2", "4.1e-2", "1.2e-1")),
    "ϑ = 1e-5" => (10, "4.2e-6", ("2.7e-4", "1.3e-3", "5.5e-3", "2.0e-2")),
    "ϑ = 1e-6" => (13, "1.9e-7", ("2.3e-5", "8.7e-5", "5.8e-4", "2.9e-3")),
    "rank 1" => (1, nothing, ("1.5e-1", "1.8e-1", "1.4e-1", "1.8e-1")),
    "projection" => (nothing, nothing, ("4.1e-2", "1.1e-1", "8.4e-2", "2.0e-1")))

# with w*: empty disk, fill, density; charge-conserving: empty disk, fill, charge, max |ln|
const REMAP_TABLE = Dict(
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
            @test r.ΔP < 1e-12                 # the note quotes ≤ 3.5e-17 relative
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
                    @test s.mass_error < 2e-13     # the manual: below 2e-13
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
        end
    end
end
