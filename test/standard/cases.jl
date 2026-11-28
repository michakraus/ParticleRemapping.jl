using ParticleRemapping
using Random: MersenneTwister, Xoshiro
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

same_markers(a, b) = a.x == b.x && a.v == b.v && a.w == b.w && a.m == b.m

@testset "the box guard of StandardCase" begin
    mk = wave_case(; N = 16, T = 1.25).markers
    one(x, v) = 1.0
    case(x, v, w, L, V) = StandardCase(:W, (;), Markers(x, v, w; m = mk.m), L, V, one, one,
        :wave)
    @test case(mk.x, mk.v, mk.w, 2π, 2.5) isa StandardCase
    v = copy(mk.v)
    v[3] = 2.6
    @test occursin(
        r"^ArgumentError: StandardCase needs every \|v\| <= V = 2\.5, got max \|v\| = 2\.6$",
        shown(ArgumentError, () -> case(mk.x, v, mk.w, 2π, 2.5)))
    v[3] = NaN
    @test occursin(
        r"^ArgumentError: StandardCase needs every \|v\| <= V = 2\.5, got max \|v\| = NaN$",
        shown(ArgumentError, () -> case(mk.x, v, mk.w, 2π, 2.5)))
    x = copy(mk.x)
    x[5] = π
    @test occursin(
        r"^ArgumentError: StandardCase needs every x in \[-L/2, L/2\) with L = 6\.28318\d*$",
        shown(ArgumentError, () -> case(x, mk.v, mk.w, 2π, 2.5)))
    @test shown(ArgumentError, () -> case(Float64[], Float64[], Float64[], 2π, 2.5)) ==
          "ArgumentError: StandardCase needs at least one marker"
    for (V, s) in ((0.0, "0\\.0"), (NaN, "NaN"), (Inf, "Inf"))
        @test occursin(
            Regex("^ArgumentError: StandardCase needs a finite L > 0 and V > 0, " *
                  "got L = 6\\.28318\\d*, V = $s\$"),
            shown(ArgumentError, () -> case(mk.x, mk.v, mk.w, 2π, V)))
    end
    for (L, s) in ((0.0, "0\\.0"), (-2π, "-6\\.28318\\d*"), (NaN, "NaN"), (Inf, "Inf"))
        @test occursin(
            Regex("^ArgumentError: StandardCase needs a finite L > 0 and V > 0, " *
                  "got L = $s, V = 2\\.5\$"),
            shown(ArgumentError, () -> case(mk.x, mk.v, mk.w, L, 2.5)))
    end
    # the box comes first: with V = NaN no marker passes |v| <= V
    @test occursin("finite L > 0",
        shown(ArgumentError, () -> case(mk.x, v, mk.w, 2π, NaN)))
    # L and V are stored as Float64 and print as such
    c = case(mk.x, mk.v, mk.w, 2π, 3)
    @test c.V === 3.0 && c.L === Float64(2π)
    @test occursin("V = 1.0, got",
        shown(ArgumentError, () -> case(mk.x, mk.v, mk.w, 2π, 1)))
    @test occursin("unknown family",
        shown(ArgumentError,
            () -> StandardCase(:X, (;), mk, 2π, 3.0, one, one, :wave)))
    @test occursin("unknown experiment",
        shown(ArgumentError,
            () -> StandardCase(:W, (;), mk, 2π, 3.0, one, one, :hole)))
end

@testset "the two densities of a case" begin
    w = wave_case(; T = 5)
    xc = [-π + (i - 0.5) * 2π / 10 for i in 1:10, j in 1:10]
    vc = [-3 + (j - 0.5) * 6 / 10 for i in 1:10, j in 1:10]
    r = case_density(w, :represented).(xc, vc)
    p = case_density(w, :positive).(xc, vc)
    @test count((r .== 0) .& (p .> 0)) == 52
    @test case_density(w, :represented) === w.represented
    @test case_density(w, :positive) === w.positive
    c = hole_case(; k = 4)
    # the centre of the empty cell i = 1, j = 0
    @test case_density(c, :represented)(-1 + 3 / 4, -1 + 1 / 4) == 0
    @test case_density(c, :positive)(-1 + 3 / 4, -1 + 1 / 4) == 0.25
    for case in (lattice_case(), (stress_case(s) for s in (:S1, :S2, :S3, :S4))...)
        @test case.represented === case.positive
    end
    for case in (w, c, lattice_case())
        @test shown(ArgumentError, () -> case_density(case, :uniform)) ==
              "ArgumentError: variant must be :represented or :positive, got :uniform"
    end
end

@testset "the ladders" begin
    @test STANDARD_LADDERS == (
        W = (N = (16, 144, 576, 900), T = (0.0, 1.25, 2.5, 5.0, 10.0, 15.0, 25.0),
            seed = 0:4),
        C = (k = (1, 2, 4, 8), N = 400, seed = 505:509),
        E = (N = (144, 576), deformed = (false, true)),
        S = (name = (:S1, :S2, :S3, :S4), N = (144, 576), seed = 0:4))
    for name in (:StandardCase, :STANDARD_LADDERS, :wave_case, :hole_case, :lattice_case,
        :stress_case, :case_density)
        @test Base.isexported(ParticleRemapping, name)
    end
end

@testset "W: seed 0 is the old lattice and the old long run" begin
    for (case, ref, steps, n) in (
        (wave_case(; N = 900, T = 25), standard_wave_markers(), 500, 30),
        (wave_case(; N = 144, T = 0), electric_wave_markers(; Nx = 12, Nv = 12), 0, 12),
        (
        wave_case(; N = 16, T = 1.25), standard_wave_markers(; Nx = 4, Nv = 4, steps = 25),
        25, 4))
        @test same_markers(case.markers, ref)
        @test case.markers.m == 1.2
        @test case.L ≈ 2π
        @test case.experiment === :wave && case.family === :W
        @test case.parameters.steps == steps
        @test case.V == 3.0
        mk = case.markers
        for variant in (:represented, :positive)
            f = case_density(case, variant)
            g(x, v) = electric_wave_density(x, v; variant, steps, Nx = n, Nv = n)
            @test isapprox(f.(mk.x, mk.v), g.(mk.x, mk.v); rtol = 1e-14)
            @test isapprox(f(0.0, 2.4), g(0.0, 2.4); rtol = 1e-14)
        end
    end
end

@testset "W: one box for every rung" begin
    for T in (0, 1.25, 25)
        case = wave_case(; N = 144, T)
        @test case.V == 3.0
        @test case.parameters.V == 3.0
    end
    case = wave_case(; N = 144, T = 0)
    @test case_evenness(case, case.markers).empty == 40
    @test case.represented(0.0, 2.4) == 0
    @test case.positive(0.0, 2.4) > 0
end

@testset "W: the box as a keyword, and the longest rungs" begin
    vmax = 0.0
    for N in (576, 900), seed in 0:4

        narrow = wave_case(; N, T = 25, seed, V = 2.5)
        @test narrow.V == 2.5
        wide = wave_case(; N, T = 25, seed)
        @test wide.V == 3.0
        @test same_markers(narrow.markers, wide.markers)
        vmax = max(vmax, maximum(abs, wide.markers.v))
    end
    @test 2.3 < vmax < 2.5
    @test occursin(
        r"^ArgumentError: StandardCase needs every \|v\| <= V = 2\.0, got max \|v\| = 2\.15107836\d*$",
        shown(ArgumentError, () -> wave_case(; N = 144, T = 1.25, V = 2.0)))
end

@testset "W: the stratified seeds" begin
    n = 12
    lattice = electric_wave_markers(; Nx = n, Nv = n)
    dx0, dv0 = 2π / n, 2.8 / n
    f0(x, v) = (1 + 0.35 * cos(x)) * (0.40 + v^2) * exp(-v^2 / 1.4)
    sets = []
    for s in 1:4
        case = wave_case(; N = 144, T = 0, seed = s)
        mk = case.markers
        rng = Xoshiro(s)
        x, v = zeros(n^2), zeros(n^2)
        k = 0
        for ix in 1:n, iv in 1:n

            k += 1
            ux = rand(rng)
            uv = rand(rng)
            x[k] = -π + (ix - 1 + ux) * dx0
            v[k] = -1.4 + (iv - 1 + uv) * dv0
            # inside its own cell of the lattice, and not at its centre
            @test abs(mk.x[k] - lattice.x[k]) < dx0 / 2
            @test abs(mk.v[k] - lattice.v[k]) < dv0 / 2
            @test (mk.x[k], mk.v[k]) != (lattice.x[k], lattice.v[k])
        end
        w = f0.(x, v) .* dx0 .* dv0
        w ./= sum(w)
        @test isapprox(mk.x, x; rtol = 1e-14)
        @test isapprox(mk.v, v; rtol = 1e-14)
        @test isapprox(mk.w, w; rtol = 1e-14)
        @test abs(sum(mk.w) - 1) < 1e-14
        @test same_markers(mk, wave_case(; N = 144, T = 0, seed = s).markers)
        push!(sets, mk.x)
    end
    @test allunique(sets)
end

@testset "W: the guards of N, T and the seed" begin
    for (kw, msg) in (
        ((; N = 0), "wave_case needs N = n^2 with n >= 1, got N = 0"),
        ((; N = 15), "wave_case needs N = n^2 with n >= 1, got N = 15"),
        (
        (; T = -0.05), "wave_case needs T >= 0 on the step grid T = 0.05 k, got T = -0.05"),
        ((; T = 0.07), "wave_case needs T >= 0 on the step grid T = 0.05 k, got T = 0.07"),
        ((; T = NaN), "wave_case needs T >= 0 on the step grid T = 0.05 k, got T = NaN"),
        ((; T = Inf), "wave_case needs T >= 0 on the step grid T = 0.05 k, got T = Inf"),
        (
        (; T = 1e20), "wave_case needs T >= 0 on the step grid T = 0.05 k, got T = 1.0e20"),
        ((; seed = -1), "wave_case needs a seed >= 0, got seed = -1"))
        @test shown(ArgumentError, () -> wave_case(; kw...)) == "ArgumentError: " * msg
    end
    @test occursin(r"^ArgumentError: wave_case needs N = n\^2 with n >= 1, got N = 15$",
        shown(ArgumentError, () -> wave_case(; N = 15)))
    @test wave_case(; T = 1.25).parameters.steps == 25
    @test wave_case(; T = 2.5).parameters.steps == 50
    mk = wave_case(; seed = 7).markers
    @test all(abs.(mk.v) .<= 3.0) && all(-π .<= mk.x .< π)
end

@testset "C: k = 2 is the old checkerboard" begin
    for seed in 505:509
        case = hole_case(; k = 2, seed)
        @test same_markers(case.markers, checkerboard_markers(; seed, N = 400))
        @test case.L == 2.0 && case.V == 1.0 && case.experiment === :checkerboard
        @test case.family === :C
    end
    case = hole_case()
    for variant in (:represented, :positive), x in -1:0.05:1, v in -1:0.05:1
        @test case_density(case, variant)(x, v) == checkerboard_density(x, v; variant)
    end
    r = case.represented
    @test r(1.0, 0.5) == 0.5 && r(0.5, 1.0) == 0.5
    @test (r(1.0, 1.0), r(-1.0, -1.0), r(1.0, -0.5), r(1.05, 0.5)) == (0.5, 0.5, 0.0, 0.0)
    for k in (1, 4)
        p = hole_case(; k).positive
        @test (p(1.0, 1.0), p(1.05, 0.5), p(0.5, -1.05)) == (0.25, 0.0, 0.0)
    end
end

@testset "C: the other orders" begin
    cell(z, k) = min(floor(Int, (z + 1) * k / 2), k - 1)
    for k in (1, 4, 8)
        case = hole_case(; k, seed = 505)
        mk = case.markers
        @test length(mk) == 400
        @test all(mk.w .== 1 / 400)
        counts = zeros(Int, k, k)
        for (x, v) in zip(mk.x, mk.v)
            counts[cell(x, k) + 1, cell(v, k) + 1] += 1
        end
        occupied = [(i, j) for i in 0:(k - 1) for j in 0:(k - 1) if iseven(i + j)]
        empty = [(i, j) for i in 0:(k - 1) for j in 0:(k - 1) if isodd(i + j)]
        @test all(counts[i + 1, j + 1] == 0 for (i, j) in empty)
        c = [counts[i + 1, j + 1] for (i, j) in occupied]
        k == 1 && @test c == [400]
        k == 4 && @test all(==(50), c)
        k == 8 && @test c == vcat(fill(13, 16), fill(12, 16))
        @test all(case.represented.(mk.x, mk.v) .> 0)
        @test all(case.positive.(mk.x, mk.v) .== 0.25)
        for (i, j) in empty
            @test case.represented(-1 + (2i + 1) / k, -1 + (2j + 1) / k) == 0
        end
    end
    # the draws of k = 4 as the definition states them: 50 markers in each occupied cell
    rng = MersenneTwister(505)
    x, v = Float64[], Float64[]
    for i in 0:3, j in 0:3

        isodd(i + j) && continue
        for _ in 1:50
            push!(x, -1 + 2i / 4 + (2 / 4) * rand(rng))
            push!(v, -1 + 2j / 4 + (2 / 4) * rand(rng))
        end
    end
    mk = hole_case(; k = 4).markers
    @test isapprox(mk.x, x; rtol = 1e-14) && isapprox(mk.v, v; rtol = 1e-14)
    for k in (1, 2, 4, 8)
        f = hole_case(; k).represented
        h = 2 / 400
        q = sum(f(-1 + (i - 0.5) * h, -1 + (j - 0.5) * h) for i in 1:400, j in 1:400) * h^2
        @test abs(q - 1) < 1e-12
    end
end

@testset "C: the guards of k, N and the seed" begin
    for (kw, msg) in (
        ((; k = 0), "hole_case needs k >= 1 and N >= 1, got k = 0, N = 400"),
        ((; N = 0), "hole_case needs k >= 1 and N >= 1, got k = 2, N = 0"),
        ((; k = 2, N = 401), "hole_case needs an even N at k = 2, got N = 401"),
        ((; seed = -1), "hole_case needs a seed >= 0, got seed = -1"))
        @test shown(ArgumentError, () -> hole_case(; kw...)) == "ArgumentError: " * msg
    end
    @test length(hole_case(; k = 4, N = 401).markers) == 401
    case = hole_case(; k = 3)
    mk = case.markers
    @test length(mk) == 400
    cell(z) = min(floor(Int, (z + 1) * 3 / 2), 2)
    @test length(unique(zip(cell.(mk.x), cell.(mk.v)))) == 5
    @test all(isapprox.(case.represented.(mk.x, mk.v), 0.45; rtol = 1e-14))
end

@testset "E: the plain and the deformed lattice" begin
    for N in (144, 576), deformed in (false, true)

        case = lattice_case(; N, deformed)
        a, b = deformed ? (0.2, 0.3) : (0.0, 0.0)
        ex = lattice_example(; nL = isqrt(N), a, b)
        mk = case.markers
        @test mk.x == ex.markers.x .- π
        @test mk.v == ex.markers.v && mk.w == ex.markers.w
        @test all(-π .<= mk.x .< π)
        @test case.L ≈ 2π && case.V == 4.0 && case.markers.m == 1.0
        @test case.family === :E && case.experiment === :wave
        @test case.parameters == (; N, deformed)
        if !deformed
            @test isapprox(case.represented.(mk.x, mk.v) .* (2 * 2π * 4 / N), mk.w;
                rtol = 1e-14)
            @test mk.x != lattice_case(; N, deformed = true).markers.x
        end
    end
end

@testset "E: a size that is no square" begin
    for N in (150, 0)
        @test shown(ArgumentError, () -> lattice_case(; N)) ==
              "ArgumentError: lattice_case needs N = n^2 with n >= 1, got N = $N"
    end
end

@testset "S1 and S3, the samples" begin
    for N in (144, 576), seed in (0, 1)

        s1 = stress_case(:S1; N, seed)
        s3 = stress_case(:S3; N, seed)
        for case in (s1, s3)
            mk = case.markers
            @test length(mk) == N
            @test all(abs.(mk.v) .<= case.V) && all(-π .<= mk.x .< π)
        end
        @test s1.V == 4.0 && s3.V == 3.0
        @test all(isapprox.(s1.markers.w, 2π * 0.9999366575 / N; rtol = 1e-9))
        @test all(s3.markers.w .== 1 / N)
        @test same_markers(s1.markers, stress_case(:S1; N, seed).markers)
        @test same_markers(s3.markers, stress_case(:S3; N, seed).markers)
        @test s1.markers.x != stress_case(:S1; N, seed = 1 - seed).markers.x
        @test s3.markers.x != stress_case(:S3; N, seed = 1 - seed).markers.x
    end
    narrow = stress_case(:S1; N = 144, V = 2.5)
    @test length(narrow.markers) == 144 && narrow.V == 2.5
    @test all(abs.(narrow.markers.v) .<= 2.5)
    @test all(isapprox.(narrow.markers.w, 2π * 0.9875806703 / 144; rtol = 1e-9))
    f = stress_case(:S3).represented
    σx, σv = 0.1 * 2π, 0.6
    images = sum(exp(-(3.0 + k * 2π)^2 / (2σx^2)) for k in -1:1) / (2π * σx * σv)
    central = exp(-3.0^2 / (2σx^2)) / (2π * σx * σv)
    @test isapprox(f(3.0, 0.0), images; rtol = 1e-14)
    @test abs(f(3.0, 0.0) - central) > 0.05 * central
    @test isapprox(f(-π, 0.4), f(π - 1e-12, 0.4); rtol = 1e-9)
    # the draws as the definitions state them
    for (V, seed) in ((4.0, 0), (2.5, 1))
        rng = Xoshiro(seed)
        x, v = Float64[], Float64[]
        while length(x) < 144
            xa = -π + 2π * rand(rng)
            1.1 * rand(rng) > 1 + 0.1 * cos(xa) && continue
            va = randn(rng)
            abs(va) > V && continue
            push!(x, xa)
            push!(v, va)
        end
        mk = stress_case(:S1; N = 144, seed, V).markers
        @test mk.x == x && mk.v == v
    end
    rng = Xoshiro(1)
    x, v = Float64[], Float64[]
    while length(x) < 576
        push!(x, ParticleRemapping.periodic_wrap(σx * randn(rng), π))
        push!(v, 0.1 * (2 * 3.0) * randn(rng))
        abs(v[end]) > 3.0 && (pop!(x); pop!(v))
    end
    mk = stress_case(:S3; N = 576, seed = 1).markers
    @test mk.x == x && mk.v == v
end

@testset "S2, the weights over decades" begin
    for (N, wmin, wmax, total, beam) in ((144, 4.8790e-08, 1.7232e-01, 6.274146, 0.618734),
        (576, 2.9860e-09, 4.7360e-02, 6.283175, 0.660346))
        mk = stress_case(:S2; N).markers
        @test isapprox(minimum(mk.w), wmin; rtol = 1e-4)
        @test isapprox(maximum(mk.w), wmax; rtol = 1e-4)
        @test isapprox(sum(mk.w), total; rtol = 1e-6)
        @test isapprox(sum(mk.w[mk.v .> 2.5]), beam; rtol = 1e-5)
    end
end

@testset "S4, the clustered example on the case box" begin
    for N in (576, 144)
        n = isqrt(N)
        case = stress_case(:S4; N)
        ex = clustered_example(; nx0 = n, nv0 = n)
        mk = case.markers
        @test mk.x == ex.markers.x .- π
        @test mk.v == ex.markers.v && mk.w == ex.markers.w
        @test case.V == 3.5 && mk.m == 1.5
        @test isapprox(case.represented.(mk.x, mk.v) .* (2π / n * 4.8 / n), mk.w;
            rtol = 1e-12)
    end
    wide = stress_case(:S4; N = 576, V = 4.0)
    @test wide.V == 4.0
    @test same_markers(wide.markers, stress_case(:S4; N = 576).markers)
    @test startswith(shown(ArgumentError, () -> stress_case(:S4; N = 576, V = 3.0)),
        "ArgumentError: clustered_example needs V ≥ 3.01488")
    @test startswith(shown(ArgumentError, () -> stress_case(:S4; N = 144, V = 2.5)),
        "ArgumentError: clustered_example needs V ≥ 2.84217")
end

@testset "the guards of the stress cases" begin
    for (name, kw, msg) in (
        (:S5, (;), "unknown stress case :S5; the names are (:S1, :S2, :S3, :S4)"),
        (:S2, (; N = 150), "stress_case :S2 needs N = n^2 with n >= 1, got N = 150"),
        (:S4, (; N = 150), "stress_case :S4 needs N = n^2 with n >= 1, got N = 150"),
        (:S2, (; seed = 1),
        "stress_case :S2 is deterministic; its only seed is 0, got seed = 1"),
        (:S4, (; seed = 3),
        "stress_case :S4 is deterministic; its only seed is 0, got seed = 3"),
        (:S1, (; N = 0), "stress_case :S1 needs N >= 1, got N = 0"),
        (:S3, (; seed = -1), "stress_case :S3 needs a seed >= 0, got seed = -1"),
        (:S2, (; V = 7.0),
        "stress_case :S2 keeps the box of its density; V is a keyword of :S1 and :S4 only"),
        (:S3, (; V = 3.0),
        "stress_case :S3 keeps the box of its density; V is a keyword of :S1 and :S4 only"))
        @test shown(ArgumentError, () -> stress_case(name; kw...)) ==
              "ArgumentError: " * msg
    end
    # a box that is not finite and positive is refused before any marker is drawn
    for (name, V, shownV) in ((:S1, 0.0, "0\\.0"), (:S1, -1.0, "-1\\.0"), (
        :S1, -Inf, "-Inf"),
        (:S1, 0, "0\\.0"), (:S4, 0.0, "0\\.0"), (:S4, NaN, "NaN"))
        @test occursin(
            Regex("^ArgumentError: StandardCase needs a finite L > 0 and V > 0, " *
                  "got L = 6\\.28318\\d*, V = $shownV\$"),
            shown(ArgumentError, () -> stress_case(name; V)))
    end
    @test length(stress_case(:S1; N = 7).markers) == 7
end
