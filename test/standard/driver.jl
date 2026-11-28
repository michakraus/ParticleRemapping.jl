using ParticleRemapping
using Test

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

const MEASURES = (:evenness, :outside, :representation, :charge_noise, :weights)
const CODES = (
    "SUCCESS", "STALLED", "MAXITERS", "SINGULAR", "NONFINITE", "LINESEARCH_FAILED")

# Methods of the test: one that throws an error outside the failures of a remap, one that throws
# a DomainError, one that returns a non-finite velocity, and one that counts its calls.
struct Boom <: AbstractRemap end
struct Domain <: AbstractRemap end
struct NotFinite <: AbstractRemap end
struct Counting <: AbstractRemap
    calls::Base.RefValue{Int}
end
ParticleRemapping.remap(::Boom, mk::Markers) = error("boom")
ParticleRemapping.remap(::Domain, mk::Markers) = throw(DomainError(0.0, "x"))
function ParticleRemapping.remap(::NotFinite, mk::Markers)
    v = copy(mk.v)
    v[2] = NaN
    RemapResult(Markers(copy(mk.x), v, copy(mk.w); m = mk.m), (;))
end
# a non-finite position at one marker and a non-finite weight at another
struct NotFiniteXW <: AbstractRemap end
function ParticleRemapping.remap(::NotFiniteXW, mk::Markers)
    x, w = copy(mk.x), copy(mk.w)
    x[3] = NaN
    w[5] = Inf
    RemapResult(Markers(x, copy(mk.v), w; m = mk.m), (;))
end
function ParticleRemapping.remap(c::Counting, mk::Markers)
    c.calls[] += 1
    RemapResult(copy(mk), (;))
end

@testset "the box of the case in every method" begin
    case = wave_case(; N = 144, V = 2.5)
    m = Dict(key => case_method(key, case) for key in STANDARD_METHODS)
    for key in (:projection, :nonseparable)
        @test m[key] isa ShiftedRemap
        @test m[key].method.V == 2.5
        @test m[key].method.L ≈ 2π
        @test m[key].shift ≈ π
        @test m[key].L ≈ 2π
        @test m[key].method.projection.nX == 16 && m[key].method.projection.nV == 17
    end
    @test m[:projection].method isa KickDriftRemap
    @test m[:projection].method.maxit == 20
    @test m[:nonseparable].method isa NonSeparableRemap
    @test (m[:nonseparable].method.nξ, m[:nonseparable].method.nv,
        m[:nonseparable].method.nℓ, m[:nonseparable].method.maxit) == (6, 7, 3, 5)
    spline = m[:spline].method
    @test m[:spline] isa ShiftedRemap && m[:spline].shift ≈ π
    @test spline.Lx ≈ 2π
    @test spline.Vmax == 2.5
    @test spline.relax_step ≈ 0.006 * 16 / 144
    @test spline.σx ≈ 2π / 4 && spline.σv == 1.2
    @test m[:gauge] isa GaugeRezoning && m[:gauge].grid.Ng == 16 && m[:gauge].grid.L ≈ 2π
    @test m[:transport] isa CurrentQuantileMatching
    @test m[:transport].target.V[end] == 2.5
    @test m[:transport].grid.Ng == 16
    @test m[:sliced].Lx ≈ π
    @test m[:sliced].Vmax == 2.5
    @test m[:sliced].hx == 0.45 && m[:sliced].hv == 0.35
    @test case_method(:projection, wave_case(; N = 144)).method.V == 3.0
    @test case_method(:transport, wave_case(; N = 144, T = 0)).target.V[end] == 3.0
    @test case_method(:projection, stress_case(:S1; V = 3.5)).method.V == 3.5
    @test case_method(:projection, lattice_case()).method.V == 4.0
    hole = hole_case()
    s = case_method(:sliced, hole)
    @test (s.Lx, s.Vmax, s.hx, s.hv, s.B_v, s.B_x) == (1, 1, 0.15, 0.15, 10, 20)
    @test case_method(:spline, wave_case(; N = 900)).method.relax_step == 0.006
    sp = case_method(:spline, hole).method
    @test sp.relax_step == 0.006 && sp.Vmax == 1.0 && sp.σv == 0.5 && sp.Lx == 2.0
    p = case_method(:projection, hole)
    @test p.method.V == 1.0 && p.method.L == 2.0 && p.shift == 1.0
    @test p.method.projection.nX == 12 && p.method.projection.nV == 13
    @test case_method(:gauge, hole).grid.Ng == 8
    @test case_method(:transport, hole).target.V[end] == 1.0
    # the density of the case, shifted by L/2 as standard_method shifts it
    f = case_density(case, :represented)
    pr = case_method(:projection, case; variant = :represented).method
    @test pr.density([0.5, 4.0], [0.3, -1.0]) == f.([0.5 - π, 4.0 - π], [0.3, -1.0])
    tr = case_method(:transport, case; variant = :represented)
    @test tr.target.F == DensityTarget(PeriodicGrid(2π, 16),
        (x, v) -> f(ParticleRemapping.periodic_wrap(x, π), v); vmax = 2.5).F
    # the references
    gl = case_method(:grid_lumped, case)
    @test gl isa GridResample && gl.density === nothing
    @test gl.L ≈ 2π && gl.V == 2.5 && gl.projection.nX == 16 && gl.projection.nV == 17
    ge = case_method(:grid_exact, case; variant = :represented)
    @test ge.density === case.represented
    @test case_method(:grid_exact, case).density === case.positive
    @test case_method(:grid_lumped, hole).projection.nX == 12
end

@testset "the slice counts of the sliced map" begin
    for (N, B) in ((16, (2, 8)), (144, (6, 36)), (576, (12, 72)), (900, (15, 90)))
        s = case_method(:sliced, wave_case(; N))
        @test (s.B_v, s.B_x) == B
        @test N % s.B_v == 0 && N % s.B_x == 0
    end
    s = case_method(:sliced, hole_case(; k = 4, N = 401))
    @test (s.B_v, s.B_x) == (1, 1)
end

@testset "an unknown key" begin
    @test shown(ArgumentError, () -> case_method(:hole, wave_case())) ==
          "ArgumentError: unknown construction :hole; the keys are (:projection, " *
          ":nonseparable, :spline, :gauge, :transport, :sliced, :grid_lumped, :grid_exact)"
    @test REFERENCE_METHODS == (:grid_lumped, :grid_exact)
end

@testset "the reference of a row" begin
    case = wave_case(; N = 144, T = 5)
    ref = case_reference(case, :represented)
    a = run_case(:transport, case; variant = :represented)
    b = run_case(:transport, case; variant = :represented, reference = ref)
    after = remap(case_method(:transport, case; variant = :represented), case.markers).markers
    for row in (a, b)
        @test row.before.representation ==
              case_representation(case, case.markers, :represented)
        @test row.before.charge_noise == case_charge_noise(case, case.markers, :represented)
        @test row.after.representation == case_representation(case, after, :represented)
        @test row.after.charge_noise == case_charge_noise(case, after, :represented)
    end
    @test shown(ArgumentError, () -> run_case(:transport, case; reference = ref)) ==
          "ArgumentError: the reference is of the variant :represented, the measure asks " *
          "for :positive"
    # no remap runs before the refusal
    calls = Ref(0)
    @test_throws ArgumentError run_case(Counting(calls), case; reference = ref)
    @test calls[] == 0
end

@testset "rows of constructions that succeed" begin
    case = wave_case(; N = 144, T = 5)
    ref = case_reference(case, :positive)
    rows = Dict(key => run_case(key, case; reference = ref)
    for key in (:gauge, :transport, :sliced, :grid_lumped, :grid_exact))
    for (key, row) in rows
        @test row.label == string(key)
        @test row.status === :success
        @test row.error_class === nothing && row.error_message === nothing
        @test row.N == 144 && row.family === :W && row.variant === :positive
        @test row.parameters == case.parameters
        @test row.time >= 0
        @test row.julia == string(VERSION)
        @test row.package == string(pkgversion(ParticleRemapping))
        @test keys(row.before) == MEASURES && keys(row.after) == MEASURES
        @test row.before.outside.count == 0
        @test row.solver.code === missing
        @test row.solver.iterations === missing && row.solver.objective_ratio === missing
        @test keys(row.conservation) ==
              (:weight, :momentum, :kinetic_energy, :charge, :current)
        @test keys(row.report) == (:weight, :momentum, :kinetic_energy)
        @test keys(row.displacement) ==
              (:max_dx, :max_dv, :max_weight_change, :kinetic_energy_change)
    end
    g = rows[:gauge]
    @test g.displacement.max_dx == 0
    @test g.conservation.weight < 1e-10
    @test g.conservation.charge < 1e-10
    @test g.conservation.momentum < 1e-10
    @test g.solver.tau == 1
    t = rows[:transport]
    @test t.displacement.max_dx == 0
    @test t.conservation.momentum < 1e-12
    @test t.solver.tau === missing
    s = rows[:sliced]
    @test s.conservation.weight < 1e-12
    @test s.conservation.momentum < 1e-12
    @test s.status === :success
    @test s.after.outside.count > 0
    @test s.after.outside.max_v > 3.0
    @test rows[:grid_lumped].conservation.weight < 1e-12
    @test rows[:grid_exact].conservation.weight < 1e-12
    @test fieldnames(typeof(g)) ==
          (:label, :family, :parameters, :N, :variant, :status, :error_class,
        :error_message, :time, :julia, :package, :before, :after, :conservation, :report,
        :displacement, :solver)
end

@testset "the row of a fit" begin
    row = run_case(:nonseparable, wave_case(; N = 144, T = 5))
    @test row.status === :success
    @test row.solver.code in CODES
    @test row.solver.iterations in 1:5
    @test isfinite(row.solver.objective_ratio) && row.solver.objective_ratio > 0
    @test row.solver.tau === missing
end

@testset "a failed remap is recorded" begin
    case = wave_case(; N = 144, T = 5)
    row = run_case(:projection, case; variant = :represented)
    @test row.status === :failure
    @test row.error_class == "DomainError"
    @test occursin("the target weights must be positive", row.error_message)
    @test row.variant === :represented
    @test row.after === nothing && row.solver === nothing
    @test row.conservation === nothing && row.report === nothing
    @test row.displacement === nothing
    @test keys(row.before) == MEASURES
    row = run_case(SlicedQuantileRemap(; Lx = π, Vmax = 3.0), case)
    @test row.status === :failure
    @test row.label == "SlicedQuantileRemap"
    @test row.error_class == "ArgumentError"
    @test row.error_message ==
          "ArgumentError: 144 markers do not split into 15 equal slices"
end

@testset "an exception outside the failures, and a non-finite result" begin
    case = wave_case(; N = 16, T = 0)
    @test_throws ErrorException("boom") run_case(Boom(), case)
    row = run_case(Domain(), case)
    @test row.status === :failure && row.error_class == "DomainError"
    @test row.error_message == sprint(showerror, DomainError(0.0, "x"))
    row = run_case(NotFinite(), case)
    @test row.status === :failure
    @test row.error_class == "NonFiniteResult"
    @test row.error_message ==
          "1 of 16 markers have a non-finite position, velocity or weight after the remap"
    @test row.after === nothing && row.solver === nothing
    row = run_case(NotFiniteXW(), case)
    @test row.error_class == "NonFiniteResult"
    @test row.error_message ==
          "2 of 16 markers have a non-finite position, velocity or weight after the remap"
end

@testset "the warm-up and the exported driver" begin
    case = wave_case(; N = 16, T = 0)
    calls = Ref(0)
    row = run_case(Counting(calls), case; warmup = true)
    @test calls[] == 2
    @test row.status === :success && row.label == "Counting"
    calls[] = 0
    run_case(Counting(calls), case)
    @test calls[] == 1
    @test run_case(:gauge, wave_case(; N = 144, T = 5); warmup = true).status === :success
    for name in (:case_method, :run_case, :GridResample, :REFERENCE_METHODS)
        @test Base.isexported(ParticleRemapping, name)
    end
end
