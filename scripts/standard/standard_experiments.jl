# The numbers of the standard experiments: the tables of the pages "Numerical experiments" of every
# chapter and of the comparison.
#
#     julia --startup-file=no --project=scripts scripts/standard/standard_experiments.jl
#
# Runs every construction on the checkerboard and on the long run, with both densities, and the
# runs with a remap after every 100 steps, about 15 minutes. The spline-projection redistribution
# also runs its reduced sizes. Prints the return code of each fit of the projection maps with the
# positive density, and checks every number of these rows that the pages "Numerical experiments" of
# the projection and "Comparison" quote, with the digits they are quoted with, a number quoted as a
# range between its ends.

using ParticleRemapping
using Printf
using Test

short(e) = first(replace(e, "\n" => " "), 100)

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

# The quoted numbers of the projection maps with the positive density, and the return code and
# the iterations of the fit. Checkerboard: the quadrant counts, the weight in Q2 ∪ Q4 and the empty
# cells of 64. Long run: the empty cells, the relative std, the L² difference, the change of the
# kinetic energy and the largest |v|. No fit ends on the stationarity test, so a number whose
# digits move when the targets move by 1e-13 relative, or on Julia 1.13 and 1.14, or that lies
# within a tenth of a unit of a rounding boundary, is a range over those runs.
const CHECKERBOARD = Dict(
    :projection => (("160", "40", "165", "35"), "0.1875", "13", "MAXITERS", 20),
    :nonseparable => (("127", "72", "128", "73"), "0.3625", "19", "MAXITERS", 5))
const LONG_RUN = Dict(
    :projection => (("3–6", "0.808–0.837", "0.494–0.507", "1.36–1.43", "32.4–33.4"),
        "MAXITERS", 20),
    :nonseparable => (("12", "1.01", "0.623", "2.144", "65.16"), "MAXITERS", 5))

# The remaps after every 100 steps with the positive density: the return code and the iterations
# of each fit, the largest |v| after it, and the time and the error of the remap that fails.
const HISTORY = Dict(
    :projection => ((("MAXITERS", 20, "34.5"),), 10.0, "DomainError with 0.0"),
    :nonseparable => ((("MAXITERS", 5, "31.44"), ("MAXITERS", 5, "2.6e303")), 15.0,
        "DomainError with NaN"))

# A remap by `method` at the time `t` that records the return code and the iterations of its fit
# and the largest |v| after it, or the time of the remap that throws.
struct RecordedRemap{M} <: AbstractRemap
    method::M
    t::Float64
    log::Vector{Any}
end
function ParticleRemapping.remap(r::RecordedRemap, mk::Markers)
    res = try
        remap(r.method, mk)
    catch
        push!(r.log, (; r.t, failed = true))
        rethrow()
    end
    push!(r.log,
        (; r.t, failed = false, stop = stop(res.info.status), res.info.iterations,
            max_v = maximum(abs, res.markers.v)))
    res
end

# The long run with a remap by the projection map `key` after every 100 steps, with the positive
# density at the time of each remap, as `standard_remap_history` runs it; returns the record of
# the remaps and the error message.
function recorded_history(key; Nx = 30, Nv = 30, Δt = 0.05)
    log = Any[]
    method(t) = RecordedRemap(
        standard_method(key, :wave;
            density = (x, v) -> electric_wave_density(x, v; variant = :positive,
                steps = round(Int, t / Δt), Nx, Nv)),
        t, log)
    error = try
        electric_wave_simulation(electric_wave_markers(; Nx, Nv); steps = 500,
            record_every = 10, remap_method = method, remap_every = 100)
        nothing
    catch e
        sprint(showerror, e)
    end
    log, error
end

checked = Dict{Symbol, Any}()

println("== the checkerboard (before: quadrants (200, 0, 200, 0), 32 of 64 cells empty)")
for variant in (:represented, :positive), key in STANDARD_METHODS

    r = standard_checkerboard(key; variant)
    if r.after === nothing
        println(rpad("$key, $variant", 28), "fails: ", short(r.error))
    else
        @printf("%-28s quadrants %s, weight in Q2 ∪ Q4 %.3f, empty %d, rstd %.3f, max shift %.2f, ΔW %.1e, ΔP %.1e",
            "$key, $variant", r.quadrants_after, r.hole_weight, r.stats_after.empty,
            r.stats_after.rstd, r.max_shift, r.report.weight, r.report.momentum)
        if variant === :positive && haskey(CHECKERBOARD, key)
            @printf(", fit %s after %d", stop(r.result.info.status),
                r.result.info.iterations)
            checked[Symbol(:checkerboard_, key)] = r
        end
        println()
    end
end

println("\n== the long run (before: 14 empty cells, rstd 0.598)")
for variant in (:represented, :positive), key in STANDARD_METHODS

    r = standard_long_time(key; variant)
    if r.after === nothing
        println(rpad("$key, $variant", 28), "fails: ", short(r.error))
    else
        s = r.stats_after
        @printf("%-28s empty %d, min %d, max %d, rstd %.3f, L2 %.3f, ΔW %.1e, ΔK %.3g, max |v| %.3g",
            "$key, $variant", s.empty, s.min_count, s.max_count, s.rstd, r.l2_difference,
            r.report.weight, r.report.kinetic_energy, r.max_v)
        if variant === :positive && haskey(LONG_RUN, key)
            @printf(", fit %s after %d", stop(r.result.info.status),
                r.result.info.iterations)
            checked[Symbol(:long_, key)] = r
        end
        println()
    end
end

println("\n== the spline-projection redistribution at reduced sizes")
for N in (16, 36, 64)
    r = standard_checkerboard(:spline; N)
    @printf("checkerboard N = %d: quadrants %s → %s, target error %.3g\n", N,
        r.quadrants_before, r.quadrants_after, r.result.info.target_error)
end
for n in (4, 6, 8)
    r = standard_long_time(:spline; Nx = n, Nv = n)
    @printf("long run %d × %d: empty %d → %d, target error %.3g\n", n, n,
        r.stats_before.empty, r.stats_after.empty, r.result.info.target_error)
end

println("\n== remaps after every 100 steps (mean rstd for t ≥ 5; empty cells at t = 0, 12.5, 25)")
for key in (:projection, :nonseparable)
    log, error = recorded_history(key)
    for l in log
        l.failed ? @printf("%-28s t = %4.1f: fails\n", "$key, positive", l.t) :
        @printf("%-28s t = %4.1f: fit %s after %d, max |v| after %.3g\n", "$key, positive",
            l.t, l.stop, l.iterations, l.max_v)
    end
    println(rpad("$key, positive", 28), "fails: ", short(something(error, "no error")))
    checked[Symbol(:history_, key)] = (; log, error)
end
for (key, variant) in ((:gauge, :represented), (:transport, :represented),
    (:transport, :positive), (:sliced, :represented))
    h = standard_remap_history(key; variant)
    if h.error === nothing
        @printf("%-28s %.3f → %.3f, empty %s\n", "$key, $variant", h.rstd_plain,
            h.rstd_remap, h.h_remap.empty[[1, 26, 51]])
    else
        println(rpad("$key, $variant", 28), "fails: ", short(h.error))
    end
end

@testset "the quoted numbers of the projection maps, positive density" begin
    for (key, (quadrants, weight, empty, code, iterations)) in CHECKERBOARD
        r = checked[Symbol(:checkerboard_, key)]
        for (x, q) in zip(r.quadrants_after, quadrants)
            @test agrees(x, q)
        end
        @test agrees(r.hole_weight, weight)
        @test agrees(r.stats_after.empty, empty)
        @test abs(r.report.weight) < 1e-12 && abs(r.report.momentum) < 1e-12
        @test stop(r.result.info.status) == code && r.result.info.iterations == iterations
    end
    for (key, (row, code, iterations)) in LONG_RUN
        r = checked[Symbol(:long_, key)]
        values = (r.stats_after.empty, r.stats_after.rstd, r.l2_difference,
            r.report.kinetic_energy, r.max_v)
        for (x, q) in zip(values, row)
            @test agrees(x, q)
        end
        @test stop(r.result.info.status) == code && r.result.info.iterations == iterations
    end
    for (key, (fits, t_fail, error)) in HISTORY
        h = checked[Symbol(:history_, key)]
        @test length(h.log) == length(fits) + 1
        for (l, (code, iterations, max_v)) in zip(h.log, fits)
            @test !l.failed && l.stop == code && l.iterations == iterations &&
                  agrees(l.max_v, max_v)
        end
        @test last(h.log).failed && last(h.log).t == t_fail
        @test h.error !== nothing && startswith(h.error, error)
    end
end
