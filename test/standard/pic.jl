using ParticleRemapping
using Random: Xoshiro
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

# VlasovMethods' two-stream reference f₀ = (1 + α cos(k x)) v² e^{−v²/2}/√(2π) on [0, 4π).
f0(x, v) = (1 + 1e-3 * cos(0.5 * x)) * v^2 * exp(-v^2 / 2) / sqrt(2π)

# The n × n cell centres of [−2π, 2π) × [−5, 5], i outer, j inner, and f₀(x + 2π, v) dx dv.
function lattice(n)
    L, VL = 4π, 5.0
    dx, dv = L / n, 2VL / n
    x = [-L / 2 + (i - 0.5) * dx for i in 1:n for j in 1:n]
    v = [-VL + (j - 0.5) * dv for i in 1:n for j in 1:n]
    x, v, f0.(x .+ L / 2, v) .* (dx * dv)
end

# The equal-weight sample of P1, drawn by the algorithm of its definition.
function sample(N, seed)
    L, VL, α, k = 4π, 5.0, 1e-3, 0.5
    rng = Xoshiro(seed)
    x, v = Float64[], Float64[]
    while length(x) < N
        xa = -L / 2 + L * rand(rng)
        (1 + α) * rand(rng) > 1 + α * cos(k * (xa + L / 2)) && continue
        s = sqrt(randn(rng)^2 + randn(rng)^2 + randn(rng)^2)
        s > VL && continue
        push!(x, xa)
        push!(v, rand(rng) < 0.5 ? -s : s)
    end
    x, v
end

# Two constructions are equal field by field, down to the arrays and the numbers they hold.
function samefields(a, b)
    typeof(a) == typeof(b) || return false
    a isa Union{Number, Symbol, AbstractString, Nothing, AbstractArray} &&
        return isequal(a, b)
    a isa Function && return true  # the type, compared above, holds what a closure captures
    fieldcount(typeof(a)) == 0 && return a === b
    all(samefields(getfield(a, f), getfield(b, f)) for f in fieldnames(typeof(a)))
end

const KEYS = (:projection, :nonseparable, :spline, :gauge, :sliced, :grid_lumped)

@testset "P1: the lattice and its weights" begin
    for N in (900, 3600)
        case = N == 900 ? two_stream_case() : two_stream_case(; N)
        n = isqrt(N)
        x, v, raw = lattice(n)
        mk = case.markers
        @test mk.x == x
        @test mk.v == v
        c = mk.w[1] / raw[1]
        @test all(isapprox.(mk.w, c .* raw; rtol = 1e-14))
        @test isapprox(sum(mk.w), 4π; rtol = 1e-14)
        @test maximum(abs, mk.v) ≈ 5 * (1 - 1 / n)
        @test case isa StandardCase
        @test case.family === :P1
        @test case.experiment === :wave
        @test mk.m == 1
        @test case.L ≈ 4π
        @test case.V == 6.0
        @test case.parameters ==
              (; N, loader = :lattice, seed = 0, V = 6.0, VL = 5.0, α = 1e-3, k = 0.5,
            cells = 16, degree = 3, dt = 0.1, steps = 400, window = (10.0, 25.0),
            rate = 0.2592499317)
        @test case_density(case, :positive) isa Function
        @test case_outside(case, mk).count == 0
        N == 900 && @test isapprox(c, 1.000013946050365; rtol = 1e-12)
    end
    @test maximum(abs, two_stream_case().markers.v) ≈ 4.833333333333333
end

@testset "P1: the two densities" begin
    case = two_stream_case()
    mk = case.markers
    for f in (case_density(case, :represented), case_density(case, :positive))
        @test all(isapprox.(f.(mk.x, mk.v), f0.(mk.x .+ 2π, mk.v); rtol = 1e-14))
    end
    for (x, v) in ((0.0, 5.0), (1.0, -5.0))  # the cut keeps |v| = VL
        @test isapprox(case_density(case, :represented)(x, v), f0(x + 2π, v); rtol = 1e-14)
    end
    for (x, v) in ((0.0, 5.5), (1.0, -5.2))
        @test case_density(case, :represented)(x, v) == 0.0
        p = case_density(case, :positive)(x, v)
        @test p > 0
        @test isapprox(p, f0(x + 2π, v); rtol = 1e-14)
    end
end

@testset "P1: the guards" begin
    @test shown(ArgumentError, () -> two_stream_case(; N = 901)) ==
          "ArgumentError: two_stream_case :lattice needs N = n^2 with n >= 1, got N = 901"
    @test shown(ArgumentError, () -> two_stream_case(; N = 0)) ==
          "ArgumentError: two_stream_case :lattice needs N = n^2 with n >= 1, got N = 0"
    @test shown(ArgumentError, () -> two_stream_case(; seed = 1)) ==
          "ArgumentError: two_stream_case :lattice is deterministic; its only seed is 0, got seed = 1"
    @test shown(ArgumentError, () -> two_stream_case(; loader = :grid)) ==
          "ArgumentError: two_stream_case needs loader :lattice or :sample, got :grid"
    vmax = maximum(abs, lattice(30)[2])
    @test shown(ArgumentError, () -> two_stream_case(; V = 4.8)) ==
          "ArgumentError: StandardCase needs every |v| <= V = 4.8, got max |v| = $vmax"
    @test shown(ArgumentError, () -> two_stream_case(; loader = :sample, N = 0)) ==
          "ArgumentError: two_stream_case :sample needs N >= 1 and a seed >= 0, got N = 0, seed = 0"
    @test shown(ArgumentError, () -> two_stream_case(; loader = :sample, seed = -1)) ==
          "ArgumentError: two_stream_case :sample needs N >= 1 and a seed >= 0, got N = 900, seed = -1"
    @test two_stream_case(; V = 5.0).V == 5.0
    @test two_stream_case(; V = 5).parameters.V === 5.0
    small = two_stream_case(; loader = :sample, N = 7)
    @test length(small.markers) == 7
end

@testset "P1: the sample, its truncation and its seeds" begin
    cases = Dict{Tuple{Int, Int}, Any}()
    for (N, seed) in ((900, 1), (900, 2), (200_000, 1))
        case = two_stream_case(; loader = :sample, N, seed)
        cases[(N, seed)] = case
        mk = case.markers
        x, v = sample(N, seed)
        @test mk.x == x
        @test mk.v == v
        @test all(==(4π / N), mk.w)
        @test all(s -> abs(s) <= 5, mk.v)
        @test all(s -> -2π <= s < 2π, mk.x)
        @test case.family === :P1
        @test case.markers.m == 1
        @test case.parameters ==
              (;
            N, loader = :sample, seed, V = 6.0, VL = 5.0, α = 1e-3, k = 0.5, cells = 16,
            degree = 3, dt = 0.1, steps = 400, window = (10.0, 25.0), rate = 0.2592499317)
    end
    @test cases[(900, 1)].markers.x != cases[(900, 2)].markers.x
    again = two_stream_case(; loader = :sample, N = 900, seed = 1).markers
    once = cases[(900, 1)].markers
    @test again.x == once.x && again.v == once.v && again.w == once.w
    lat = two_stream_case()
    smp = cases[(900, 1)]
    for (x, v) in ((0.3, 1.1), (0.0, 5.5), (-6.0, -2.0))
        @test case_density(smp, :represented)(x, v) == case_density(lat, :represented)(x, v)
        @test case_density(smp, :positive)(x, v) == case_density(lat, :positive)(x, v)
    end
end

@testset "pic_method: the keys, the box and the density" begin
    @test PIC_METHODS == KEYS
    case = two_stream_case()
    built = Dict(key => pic_method(key, case) for key in PIC_METHODS)
    for key in (:projection, :nonseparable)
        m = built[key]
        @test m isa ShiftedRemap
        @test m.shift ≈ 2π
        @test m.L ≈ 4π
        @test m.method.L ≈ 4π
        @test m.method.V == 6.0
        @test m.method.density === nothing
        @test m.method.projection.nX == 16 && m.method.projection.nV == 17
    end
    @test built[:projection].method isa KickDriftRemap
    @test built[:projection].method.maxit == 20
    ns = built[:nonseparable].method
    @test ns isa NonSeparableRemap
    @test (ns.nξ, ns.nv, ns.nℓ, ns.maxit) == (6, 7, 3, 5)
    for key in (:spline, :gauge, :sliced, :grid_lumped)
        @test samefields(built[key], case_method(key, case))
    end
    sl = built[:sliced]
    @test sl.Lx ≈ 2π
    @test sl.Vmax == 6.0
    @test (sl.B_v, sl.B_x) == (15, 90)
    gr = built[:grid_lumped]
    @test gr isa GridResample
    @test gr.V == 6.0
    @test gr.density === nothing
    wide = pic_method(:projection, two_stream_case(; V = 7.0))
    @test wide.method.V == 7.0
    @test wide.method.density === nothing
    # the comparison helper tells two boxes apart
    @test !samefields(pic_method(:sliced, two_stream_case(; V = 7.0)), sl)
end

@testset "pic_method: the keys a run refuses" begin
    case = two_stream_case()
    keys = "(:projection, :nonseparable, :spline, :gauge, :sliced, :grid_lumped)"
    for key in (:transport, :grid_exact)
        @test shown(ArgumentError, () -> pic_method(key, case)) ==
              "ArgumentError: pic_method(:$key) needs the exact density, which a " *
              "particle-in-cell run knows only at t = 0; the keys are $keys"
    end
    @test shown(ArgumentError, () -> pic_method(:hole, case)) ==
          "ArgumentError: unknown construction :hole; the keys are $keys"
end
