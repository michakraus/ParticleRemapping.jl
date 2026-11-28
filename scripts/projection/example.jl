# The two examples of the projection overview page (`docs/src/projection/index.md`): the mildly
# deformed 24 × 24 lattice (`lattice_example()`) and the 24 × 24 clustered markers
# (`clustered_example()`), each remapped by the same 16/20/8 kick–drift map on 12 × 13 hats.
# Prints every quoted number and T1 (`conservation_report` on `PeriodicGrid(2π, 12)`) of the main
# example, remaps the main example again with the positions moved by 1e-13 relative (`Xoshiro(1)`
# and `Xoshiro(2)`) and the clustered example with ten such runs (`Xoshiro(1)` to `Xoshiro(10)`),
# checks that each perturbed `θ` differs from the unperturbed one, checks that every run agrees
# with each quoted number or lies inside its quoted range, and verifies that the lumped projection
# on 12 × 13 hats refuses the clustered markers with a `DomainError`.
#
#   julia --startup-file=no --project=scripts scripts/projection/example.jl
#
# Runtime: about 35 seconds.

using ParticleRemapping
using Printf
using Random
using Test

# One unit of the last digit of the number `quoted`.
function quoted_unit(quoted::AbstractString)
    m = match(r"^-?(\d+)(?:\.(\d+))?(?:e(-?\d+))?$", quoted)
    decimals = m.captures[2] === nothing ? 0 : length(m.captures[2])
    exponent = m.captures[3] === nothing ? 0 : parse(Int, m.captures[3])
    10.0^(exponent - decimals)
end

# `x` agrees with `quoted` to the digits that `quoted` shows; a range `lo–hi` holds `x` between
# its rounded ends.
function agrees(x, quoted::AbstractString)
    if occursin('–', quoted)
        lo, hi = parse.(Float64, split(quoted, '–'))
        return lo <= x <= hi
    end
    abs(x - parse(Float64, quoted)) <= 0.5 * quoted_unit(quoted) * (1 + 1e-9)
end

function median_value(v::AbstractVector)
    s = sort(v)
    n = length(s)
    isodd(n) ? s[cld(n, 2)] : 0.5 * (s[n ÷ 2] + s[n ÷ 2 + 1])
end

const QUOTED_MAIN = (
    code = "MAXITERS",
    iterations = "20",
    K0 = "3.138",
    objective_ratio = "0.220",
    dK = "0.431",
    outside = "7",
    max_v = "4.5"
)

const QUOTED_CLUSTERED = (
    outside = "201–224",
    max_v = "120–180",
    weight_share = "0.090–0.11",
    min_target_w = "2.9e-6",
    median_target_w = "1.2e-4–1.4e-4",
    mean_w = "1.7e-3"
)

function run_main()
    ex = lattice_example()
    mk = ex.markers
    method = KickDriftRemap(; L = ex.L, V = ex.V, projection = HatProjection(12, 13),
        nξ = 16, nv = 20, nw = 8, maxit = 20, density = (xs, vs) -> ex.density.(xs, vs))
    K0 = kinetic_energy(mk)
    r0 = remap(method, mk)
    @printf("== main example: lattice_example() (N = %d, V = %.1f, K0 = %.6f) ==\n",
        ex.N, ex.V, K0)
    println("T1: ", conservation_report(mk, r0.markers; grid = PeriodicGrid(2π, 12)))

    @testset "main example: unperturbed and perturbed runs" begin
        for s in 0:2
            r = s == 0 ? r0 :
                remap(method,
                Markers(mk.x .* (1 .+ 1e-13 .* randn(Xoshiro(s), ex.N)),
                    copy(mk.v), copy(mk.w); m = mk.m))
            i = r.info
            dK = kinetic_energy(r.markers) - K0
            out = abs.(r.markers.v) .> ex.V
            n_out = count(out)
            vmax = maximum(abs, r.markers.v)
            @printf("  run %d: %s, its = %d, F/F0 = %.6f, dK = %.6f, |Δw|/Σw = %.1e, |ΔP|/Σ|p| = %.1e, outside = %d, max |v| = %.4f%s\n",
                s, i.status.code, i.iterations, i.objective_ratio, dK,
                i.weight_change, i.momentum_change, n_out,
                vmax,
                s == 0 ? "" : ", θ moved = $(i.θ != r0.info.θ)")
            s == 0 || @test i.θ != r0.info.θ
            @test string(i.status.code) == QUOTED_MAIN.code
            @test !i.status.stationary
            @test i.iterations == parse(Int, QUOTED_MAIN.iterations)
            @test agrees(K0, QUOTED_MAIN.K0)
            @test agrees(i.objective_ratio, QUOTED_MAIN.objective_ratio)
            @test agrees(dK, QUOTED_MAIN.dK)
            @test n_out == parse(Int, QUOTED_MAIN.outside)
            @test agrees(vmax, QUOTED_MAIN.max_v)
            @test i.weight_change <= 1e-13
            @test i.momentum_change <= 1e-13
        end
    end
end

function run_clustered()
    cl = clustered_example()
    mk = cl.markers
    method = KickDriftRemap(; L = cl.L, V = cl.V, projection = HatProjection(12, 13),
        nξ = 16, nv = 20, nw = 8, maxit = 20, density = (xs, vs) -> cl.density.(xs, vs))
    wmean = sum(mk.w) / cl.N
    r0 = remap(method, mk)
    @printf("\n== second example: clustered_example() (N = %d, V = %.1f, mean w = %.4e) ==\n",
        cl.N, cl.V, wmean)

    @testset "clustered example: unperturbed and perturbed runs" begin
        @test agrees(wmean, QUOTED_CLUSTERED.mean_w)
        for s in 0:10
            r = s == 0 ? r0 :
                remap(method,
                Markers(mk.x .* (1 .+ 1e-13 .* randn(Xoshiro(s), cl.N)),
                    copy(mk.v), copy(mk.w); m = mk.m))
            i = r.info
            out = abs.(r.markers.v) .> cl.V
            w_out = r.markers.w[out]
            n_out = count(out)
            vmax = maximum(abs, r.markers.v)
            share = sum(w_out) / sum(r.markers.w)
            wmin = minimum(w_out)
            wmed = median_value(w_out)
            @printf("  run %d: %s, its = %d, outside = %d, max |v| = %.2f, weight share = %.4f, min w* = %.2e, median w* = %.2e%s\n",
                s, i.status.code, i.iterations, n_out, vmax, share, wmin,
                wmed,
                s == 0 ? "" : ", θ moved = $(i.θ != r0.info.θ)")
            s == 0 || @test i.θ != r0.info.θ
            @test string(i.status.code) == "MAXITERS"
            @test i.iterations == 20
            @test agrees(n_out, QUOTED_CLUSTERED.outside)
            @test agrees(vmax, QUOTED_CLUSTERED.max_v)
            @test agrees(share, QUOTED_CLUSTERED.weight_share)
            @test agrees(wmin, QUOTED_CLUSTERED.min_target_w)
            @test agrees(wmed, QUOTED_CLUSTERED.median_target_w)
        end

        @test_throws r"^DomainError(?s:.)*the target weights must be positive" remap(
            KickDriftRemap(; L = cl.L, V = cl.V, projection = HatProjection(12, 13),
                nξ = 16, nv = 20, nw = 8, maxit = 2),
            mk)
    end
end

run_main()
run_clustered()
