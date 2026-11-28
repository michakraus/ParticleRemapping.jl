using ParticleRemapping
using Random
using Test

# The type of the global logger before VlasovMethods loads, which the load must not change.
const LOGGER_BEFORE = typeof(Base.CoreLogging.global_logger())

using VlasovMethods
using ReducedComplexityModeling: Parameter, ParameterSpace

# The message that `f()` throws as an exception of type `T`, as `showerror` prints it.
function shown(T, f)
    try
        f()
    catch e
        e isa T && return sprint(showerror, e)
        rethrow()
    end
    "no exception"
end

const TWO_STREAM_RATE = 0.2592499317  # VlasovMethods' test/helpers/lowrank_bug_rates.jl
const LP2 = 2π / 0.3

# Remap methods of the test: each returns markers of its own and records what it sees.

# `method`, which keeps a copy of the markers it returns on its first call.
struct CopyFirst{M <: AbstractRemap} <: AbstractRemap
    method::M
    seen::Base.RefValue{Any}
end
function ParticleRemapping.remap(c::CopyFirst, mk::Markers)
    r = remap(c.method, mk)
    a = r.markers
    c.seen[] === nothing && (c.seen[] = Markers(copy(a.x), copy(a.v), copy(a.w); m = a.m))
    r
end

# The identity, which keeps the positions it receives; with `kick`, marker 1 leaves with v = 7.
struct Recorder <: AbstractRemap
    kick::Bool
    xs::Vector{Vector{Float64}}
end
Recorder(kick) = Recorder(kick, Vector{Float64}[])
function ParticleRemapping.remap(r::Recorder, mk::Markers)
    push!(r.xs, copy(mk.x))
    v = copy(mk.v)
    r.kick && (v[1] = 7.0)
    RemapResult(Markers(copy(mk.x), v, copy(mk.w); m = mk.m), (;))
end

# A remap that throws an exception outside the failures of a run.
struct Boom <: AbstractRemap end
ParticleRemapping.remap(::Boom, mk::Markers) = error("boom")

# The identity with the velocity of marker 1 set to NaN.
struct NaNVelocity <: AbstractRemap end
function ParticleRemapping.remap(::NaNVelocity, mk::Markers)
    v = copy(mk.v)
    v[1] = NaN
    RemapResult(Markers(copy(mk.x), v, copy(mk.w); m = mk.m), (;))
end

# The identity with every weight times 1e200, whose field energy overflows.
struct HugeWeights <: AbstractRemap end
function ParticleRemapping.remap(::HugeWeights, mk::Markers)
    RemapResult(Markers(copy(mk.x), copy(mk.v), mk.w .* 1e200; m = mk.m), (;))
end

# The identity with the position (`:x`) or the weight (`:w`) of marker 1 set to `value`.
struct OneNonFinite <: AbstractRemap
    component::Symbol
    value::Float64
end
function ParticleRemapping.remap(r::OneNonFinite, mk::Markers)
    x, w = copy(mk.x), copy(mk.w)
    r.component === :x ? (x[1] = r.value) : (w[1] = r.value)
    RemapResult(Markers(x, copy(mk.v), w; m = mk.m), (;))
end

# The markers without the last one, so that N falls by one at each remap.
struct DropLast <: AbstractRemap end
function ParticleRemapping.remap(::DropLast, mk::Markers)
    RemapResult(Markers(mk.x[1:(end - 1)], mk.v[1:(end - 1)], mk.w[1:(end - 1)]; m = mk.m),
        (;))
end

# A remap that overwrites the velocities and the weights it received, and then fails.
struct Scribble <: AbstractRemap end
function ParticleRemapping.remap(::Scribble, mk::Markers)
    fill!(mk.v, NaN)
    fill!(mk.w, NaN)
    throw(ArgumentError("scribbled"))
end

# The same stored markers at every call, as a method that returns a buffer of its own would.
struct Stored{M <: Markers} <: AbstractRemap
    markers::M
end
ParticleRemapping.remap(s::Stored, ::Markers) = RemapResult(s.markers, (;))

# Every weight divided by 1e200, which brings markers of huge weights back to the scale of P1.
struct Shrink <: AbstractRemap end
function ParticleRemapping.remap(::Shrink, mk::Markers)
    RemapResult(Markers(copy(mk.x), copy(mk.v), mk.w ./ 1e200; m = mk.m), (;))
end

@testset "the extension loads with VlasovMethods" begin
    @test Base.get_extension(ParticleRemapping, :ParticleRemappingVlasovMethodsExt) isa
          Module
    @test !isempty(methods(bump_on_tail_case))
    @test !isempty(methods(run_pic))
    @test typeof(Base.CoreLogging.global_logger()) == LOGGER_BEFORE
end

@testset "P2: the scripts' draw, on the case's box" begin
    case = bump_on_tail_case()
    Random.seed!(1234)
    p = initialize!(ParticleDistribution(1, 1, 10_000), BumpOnTail(κ = 0.3)).particles
    x, v, w = vec(p.x), vec(p.v), vec(p.w)
    mk = case.markers
    @test mk.v == v
    @test mk.w == w
    @test mk.x == mod.(x, LP2) .- LP2 / 2
    @test all(==(LP2 / 10_000), mk.w)
    @test isapprox(sum(mk.w), LP2; rtol = 1e-14)
    @test isapprox(maximum(abs, mk.v), 6.1775; rtol = 1e-4)
    @test case isa StandardCase
    @test case.family === :P2
    @test case.experiment === :wave
    @test mk.m == 1
    @test case.L ≈ LP2
    @test case.V == 8.0
    @test case.parameters ==
          (;
        N = 10_000, seed = 1234, V = 8.0, ε = 0.03, κ = 0.3, α = 0.1, σ = 0.5, v₀ = 4.5,
        cells = 16, degree = 3, dt = 0.1, steps = 250, window = (5.0, 17.5), rate = nothing)
    @test case_density(case, :positive) isa Function
    @test case_outside(case, mk).count == 0
end

@testset "P2: the caller's generator, the density and the guards" begin
    Random.seed!(99)
    case = bump_on_tail_case(; N = 100)
    a = rand()
    Random.seed!(99)
    @test a == rand()
    bot = BumpOnTail(κ = 0.3)
    for (x, v) in ((1.3, 4.2), (-9.0, 0.4), (5.0, -1.7), (10.2, 4.6), (-3.3, 6.5))
        ref = VlasovMethods.f(x + LP2 / 2, v, bot)
        @test isapprox(case_density(case, :represented)(x, v), ref; rtol = 1e-14)
        @test isapprox(case_density(case, :positive)(x, v), ref; rtol = 1e-14)
    end
    @test shown(ArgumentError, () -> bump_on_tail_case(; N = 0)) ==
          "ArgumentError: bump_on_tail_case needs N >= 1 and a seed >= 0, got N = 0, seed = 1234"
    @test shown(ArgumentError, () -> bump_on_tail_case(; seed = -1)) ==
          "ArgumentError: bump_on_tail_case needs N >= 1 and a seed >= 0, got N = 10000, seed = -1"
    other = bump_on_tail_case(; N = 2500, seed = 2)
    @test length(other.markers) == 2500
    @test other.V == 8.0
    @test other.parameters.N == 2500 && other.parameters.seed == 2
    wide = bump_on_tail_case(; N = 100, V = 9)
    @test wide.V == 9.0
    @test wide.parameters.V === 9.0
end

const P1_PLAIN = run_pic(two_stream_case())

@testset "the run without remaps: P1 reproduces the two-stream rate" begin
    r = P1_PLAIN
    @test r.status === :success
    @test length(r.W) == 401
    @test abs(r.rate / TWO_STREAM_RATE - 1) < 0.02
    @test r.max_v < 6.0
    @test all(isfinite, r.W) && all(>(0), r.W)
end

@testset "the run without remaps: P2 is VlasovMethods' own march" begin
    case = bump_on_tail_case()
    r = run_pic(case)
    N = 10_000
    p = ParticleDistribution(1, 1, N).particles
    p.x[1, :] .= case.markers.x .+ LP2 / 2
    p.v[1, :] .= case.markers.v
    p.w[1, :] .= case.markers.w
    potential() = VlasovMethods.Potential(VlasovMethods.PeriodicBasisSpline((0.0, LP2), 4, 16))
    nh = length(VlasovMethods.coefficients(PoissonField(potential())))
    SS = Snapshots(Float64, 1, N, nh, 251, ParameterSpace(Parameter(:χ, [1.0])))
    integrate_vp(p, (χ = 1.0,), PoissonField(potential()), SS,
        IntegratorParameters(0.1, 250, 251, nh, N, 1), 1)
    @test r.status === :success
    @test length(r.W) == 251
    @test all(isapprox.(r.W, SS.W[:, 1]; rtol = 1e-12))
    @test isapprox(r.rate, 0.1582; rtol = 1e-3)
    @test 6.8 < r.max_v < 7.0
end

@testset "the row of a run without remaps" begin
    sample = run_pic(two_stream_case(; loader = :sample, seed = 1))
    for r in (P1_PLAIN, sample)
        @test r.label == "none"
        @test r.family === :P1
        @test r.N == 900
        @test r.R === nothing
        @test r.steps == 400
        @test r.variant === :positive
        @test r.status === :success
        @test r.error_class === nothing && r.error_message === nothing
        @test isempty(r.remaps)
        @test r.W_reference === nothing
        @test r.rate_reference === nothing
        @test r.charge_deviation === nothing
        @test r.failed_at === nothing
        @test r.t ≈ (0:400) .* 0.1
        @test r.julia == string(VERSION)
        @test r.package == string(pkgversion(ParticleRemapping))
        @test r.vlasovmethods == string(pkgversion(VlasovMethods))
        @test r.time_push >= 0
        @test r.time_remap == 0
    end
    @test P1_PLAIN.parameters == two_stream_case().parameters
    @test P1_PLAIN.charge_noise < 0.02
    @test 0.05 < sample.charge_noise < 0.12
    @test sample.charge_noise ≈
          case_charge_noise(two_stream_case(; loader = :sample, seed = 1),
        two_stream_case(; loader = :sample, seed = 1).markers, :positive)
    short = run_pic(two_stream_case(); steps = 200)
    @test length(short.W) == 201
    @test all(isapprox.(short.W, P1_PLAIN.W[1:201]; rtol = 8eps()))
    @test short.rate === nothing  # the history ends before the window does
    # the variant enters the charge noise, and only there
    case = two_stream_case()
    represented = case_charge_noise(case, case.markers, :represented)
    @test represented != P1_PLAIN.charge_noise
    rep = run_pic(case; steps = 1, variant = :represented)
    @test rep.variant === :represented
    @test rep.charge_noise == represented
    rec = run_pic(Recorder(false), case; R = 1, steps = 1, variant = :represented)
    @test rec.variant === :represented
    @test rec.charge_noise == represented
    @test rec.W_reference == rep.W
    message = "ArgumentError: variant must be :represented or :positive, got :bogus"
    @test shown(ArgumentError, () -> run_pic(case; variant = :bogus)) == message
    bogus = Recorder(false)
    @test shown(ArgumentError, () -> run_pic(bogus, case; R = 1, variant = :bogus)) ==
          message
    @test isempty(bogus.xs)
end

@testset "the reference remap every 40 steps" begin
    r = run_pic(:grid_lumped, two_stream_case(); R = 40)
    @test r.status === :success
    @test r.label == "grid_lumped"
    @test r.R == 40
    @test length(r.remaps) == 10
    @test [m.step for m in r.remaps] == 40:40:400
    @test [m.time for m in r.remaps] isa Vector{Float64}
    @test all(m -> m.conservation.weight < 1e-12, r.remaps)
    @test all(m -> m.outside.count == 0, r.remaps)
    @test all(m -> m.outside_push.count == 0, r.remaps)
    @test all(iszero, r.charge_deviation[1:40])
    @test r.charge_deviation[41] > 1e-3
    @test length(r.charge_deviation) == 401
    @test r.W_reference == P1_PLAIN.W
    @test r.W[41] == r.remaps[1].W_after
    @test isapprox(r.W[40], P1_PLAIN.W[40]; rtol = 8eps())
    @test isapprox(r.remaps[1].W_before, P1_PLAIN.W[41]; rtol = 8eps())
    @test abs(r.rate_reference / TWO_STREAM_RATE - 1) < 0.02
    @test r.rate_reference == P1_PLAIN.rate
    @test r.time_remap > 0
end

@testset "nothing but x, v and w crosses a remap" begin
    case = two_stream_case()
    method = CopyFirst(case_method(:grid_lumped, case), Ref{Any}(nothing))
    r = run_pic(method, case; R = 40, steps = 79)
    @test r.status === :success
    @test length(r.remaps) == 1
    @test r.label == "CopyFirst"
    fresh = StandardCase(:P1, case.parameters, method.seen[], case.L, case.V,
        case.represented, case.positive, :wave)
    r2 = run_pic(fresh; steps = 39)
    @test r.W[41:80] == r2.W
end

@testset "the remap sees wrapped markers in the box, and its outside is counted after it" begin
    kick = Recorder(true)
    r = run_pic(kick, two_stream_case(); R = 40)
    same = Recorder(false)
    run_pic(same, two_stream_case(); R = 40)
    for rec in (kick, same)
        @test length(rec.xs) == 10
        @test all(xs -> all(x -> -2π <= x < 2π, xs), rec.xs)
    end
    @test r.status === :success
    @test r.remaps[1].outside.count == 1
    @test r.remaps[1].outside.max_v == 7.0
    @test r.remaps[1].outside_push.count == 0
    @test r.max_v >= 7.0
end

@testset "constructions that keep the structure, recorded raw" begin
    for key in (:gauge, :sliced)
        r = run_pic(key, two_stream_case(); R = 40)
        @test r.status === :success
        @test length(r.remaps) == 10
        @test all(m -> m.conservation.weight < 1e-12, r.remaps)
        @test all(m -> m.conservation.momentum < 1e-12, r.remaps)
        @test r.remaps[1].outside.count > 0
        @test r.time_remap > 0
    end
end

@testset "zero weights, and a push that leaves the box" begin
    r = run_pic(:grid_lumped, bump_on_tail_case(); R = 50)
    @test r.status === :success
    @test length(r.remaps) == 5
    @test r.remaps[1].weights.ratio == Inf
    @test r.remaps[1].weights.min == 0
    @test all(m -> m.conservation.weight < 1e-12, r.remaps)
    @test r.max_v > 8.0
end

@testset "a remap that throws in the middle of a run" begin
    r = run_pic(:grid_lumped, two_stream_case(; loader = :sample, N = 901, seed = 1); R = 40)
    @test r.status === :failure
    @test r.failed_at == 40
    @test r.error_class == "ArgumentError"
    @test r.error_message ==
          "ArgumentError: GridResample places an n × n lattice, so N must be a square " *
          "n^2 >= 1, got N = 901"
    @test length(r.W) == length(r.t) == 41
    @test length(r.charge_deviation) == 41
    @test r.W[41] == r.W_reference[41]  # the state before the failed remap
    @test isempty(r.remaps)
    @test r.rate === nothing
    @test length(r.W_reference) == 401
end

@testset "an exception outside the failures, a non-finite result, a non-finite field" begin
    case = two_stream_case()
    @test shown(ErrorException, () -> run_pic(Boom(), case; R = 40)) == "boom"
    @test_throws ErrorException run_pic(Boom(), case; R = 40)
    r = run_pic(NaNVelocity(), case; R = 40)
    @test r.status === :failure
    @test r.failed_at == 40
    @test r.error_class == "NonFiniteResult"
    @test r.error_message ==
          "1 of 900 markers have a non-finite position, velocity or weight after the remap"
    @test length(r.W) == 41 && isempty(r.remaps)
    h = run_pic(HugeWeights(), case; R = 40)
    @test h.status === :failure
    @test h.failed_at == 40
    @test h.error_class == "NonFiniteField"
    @test h.error_message == "the field energy is not finite after step 40"
    @test length(h.W) == 41
    @test h.W[41] == Inf
    @test length(h.t) == length(h.charge_deviation) == 41
end

@testset "the edges of a run" begin
    case = two_stream_case()
    # the fit: half the least-squares slope of log W on the closed window, here by the normal
    # equations of the line, and defined when the history ends exactly at the window's end
    r = run_pic(case; steps = 250)
    @test r.t[end] == 25.0
    sel = findall(s -> 10.0 <= s <= 25.0, r.t)
    @test length(sel) == 151
    coefficients = [ones(length(sel)) r.t[sel]] \ log.(r.W[sel])
    @test isapprox(r.rate, coefficients[2] / 2; rtol = 1e-10)
    @test isapprox(r.rate, P1_PLAIN.rate; rtol = 8eps())
    # the remap's own `|v|` enters max_v, also with no push after it
    kick = run_pic(Recorder(true), case; R = 400)
    @test kick.max_v == 7.0
    @test P1_PLAIN.max_v < 6.0
    # max_v counts the states after the pushes, not the initial markers: two like charges that
    # approach each other both slow down in the first push, below their initial |v| = 5
    pair = StandardCase(:P1, case.parameters, Markers([-1.0, 1.0], [5.0, -5.0], [2π, 2π]),
        case.L, case.V, case.represented, case.positive, :wave)
    slowed = run_pic(pair; steps = 1)
    @test slowed.status === :success
    @test 4.5 < slowed.max_v < 5.0
    # a non-finite position or weight of a remap result
    for (component, value) in ((:x, NaN), (:w, Inf))
        f = run_pic(OneNonFinite(component, value), case; R = 40)
        @test f.status === :failure
        @test f.failed_at == 40
        @test f.error_class == "NonFiniteResult"
        @test f.error_message ==
              "1 of 900 markers have a non-finite position, velocity or weight after the remap"
    end
    # a remap that changes N: the acceleration buffer follows it
    d = run_pic(DropLast(), case; R = 100)
    @test d.status === :success
    @test length(d.remaps) == 4
    @test [m.conservation.weight > 0 for m in d.remaps] == trues(4)
    @test length(d.W) == 401 && all(isfinite, d.W)
    # the march hands the remap copies, so the state before a failed remap is the march's own,
    # and it keeps none of the arrays of a result
    s = run_pic(Scribble(), case; R = 40)
    @test s.status === :failure && s.error_message == "ArgumentError: scribbled"
    @test s.W[41] == s.W_reference[41]
    @test s.charge_deviation[41] == 0
    fixed = Markers(copy(case.markers.x), copy(case.markers.v), copy(case.markers.w))
    v0, w0 = copy(fixed.v), copy(fixed.w)
    st = run_pic(Stored(fixed), case; R = 40, steps = 80)
    @test st.status === :success
    @test fixed.v == v0 && fixed.w == w0
    @test st.W[41] == st.W[81]  # both remaps place the same markers
    # a reference that fails before the run with remaps does: the deviation is NaN past its end
    mk = case.markers
    huge = StandardCase(:P1, case.parameters, Markers(mk.x, mk.v, mk.w .* 1e200), case.L,
        case.V, case.represented, case.positive, :wave)
    plain = run_pic(huge; steps = 3)
    @test plain.status === :failure && plain.error_class == "NonFiniteField"
    @test plain.failed_at == 1
    h = run_pic(Shrink(), huge; R = 1, steps = 3)
    @test h.W_reference == plain.W
    @test h.status === :success
    @test length(h.charge_deviation) == 4
    @test isfinite(h.charge_deviation[2])
    @test all(isnan, h.charge_deviation[3:4])
end

@testset "the guards of the runner" begin
    case = two_stream_case()
    for R in (0, 401)
        rec = Recorder(false)
        @test shown(ArgumentError, () -> run_pic(rec, case; R)) ==
              "ArgumentError: run_pic needs 1 <= R <= steps, got R = $R, steps = 400"
        @test isempty(rec.xs)
    end
    @test shown(ArgumentError, () -> run_pic(case; steps = 0)) ==
          "ArgumentError: run_pic needs steps >= 1, got steps = 0"
    @test shown(ArgumentError, () -> run_pic(:grid_lumped, case; R = 1, steps = 0)) ==
          "ArgumentError: run_pic needs steps >= 1, got steps = 0"
    @test shown(ArgumentError, () -> run_pic(:grid_lumped, case; R = 50, steps = 40)) ==
          "ArgumentError: run_pic needs 1 <= R <= steps, got R = 50, steps = 40"
    wave = wave_case(; N = 16)
    message = "ArgumentError: run_pic needs a case of the family :P1 or :P2, got :W"
    @test shown(ArgumentError, () -> run_pic(wave)) == message
    @test shown(ArgumentError, () -> run_pic(wave; steps = 10)) == message
    @test shown(ArgumentError, () -> run_pic(:grid_lumped, wave; R = 1)) == message
    @test shown(ArgumentError, () -> run_pic(Recorder(false), wave; R = 1, steps = 1)) ==
          message
    @test shown(ArgumentError, () -> run_pic(hole_case(); steps = 10)) ==
          "ArgumentError: run_pic needs a case of the family :P1 or :P2, got :C"
    final = Recorder(false)
    r = run_pic(final, case; R = 400)
    @test length(r.remaps) == 1
    @test r.remaps[1].step == 400
    @test length(final.xs) == 1
end
