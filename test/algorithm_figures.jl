using JSON
using ParticleRemapping
using Test

count_of(s, pattern) = length(collect(eachmatch(pattern, s)))

# the number of markers, of panels and of animated player steps of each example
const EXAMPLES = Dict(
    :projection => (100, 7, 7),
    :spline => (16, 5, 4),
    :gauge => (8, 5, 5),
    :transport => (48, 4, 4),
    :sliced => (40, 6, 6))

@testset "algorithm_svg($key)" for (key, (N, panels, nsteps)) in EXAMPLES
    s = algorithm_svg(key)
    @test s.key == key
    p = s.panels.svg
    @test count_of(p, r"<clipPath") == panels
    # every panel draws every marker, plus the open circles of the targets
    @test count_of(p, r"fill-opacity=\"0.85\"/>") == N * panels
    @test !occursin("NaN", p)

    @test sizeof(s.player) <= 150_000
    d = JSON.parse(s.player)
    @test d["version"] == 1
    @test d["id"] == string(key)
    @test d["markers"] == N
    @test length(d["window"]) == 4 && all(isfinite, d["window"])
    steps = d["steps"]
    @test length(steps) == nsteps
    for st in steps
        @test !isempty(st["caption"])
        @test length(st["x"]) == N && all(isfinite, st["x"])
        @test length(st["v"]) == N && all(isfinite, st["v"])
        @test length(st["r"]) == N && all(isfinite, st["r"])
        @test length(st["c"]) == N && all(c -> occursin(r"^#[0-9a-f]{6}$", c), st["c"])
        tg = st["targets"]
        if tg !== nothing
            @test length(tg["x"]) == N && all(isfinite, tg["x"])
            @test length(tg["v"]) == N && all(isfinite, tg["v"])
        end
    end
    if key === :spline
        @test [st["caption"] for st in steps] == [
            "1. density on a spline basis",
            "2. relaxed targets",
            "4. symplectic step, old weights",
            "5. gauge step: KKT weights"]
        @test steps[1]["targets"] === nothing
        @test steps[2]["targets"] !== nothing
        @test length(steps[2]["targets"]["x"]) == 16
        @test length(steps[2]["targets"]["v"]) == 16
    end
end

@test_throws r"^ArgumentError: unknown construction :other" algorithm_svg(:other)

@testset "write_algorithm_player" begin
    s = algorithm_svg(:spline)
    mktempdir() do dir
        path = write_algorithm_player(dir, s)
        @test path == joinpath(dir, "algorithm-spline.json")
        raw = read(path, String)
        @test raw == s.player
        @test JSON.parse(raw)["id"] == "spline"
        @test_throws SystemError write_algorithm_player(joinpath(dir, "missing"), s)
    end
end

@testset "pages" begin
    root = pkgdir(ParticleRemapping)
    pages = [
        "projection" => (:projection,),
        "redistribution" => (:spline,),
        "remapping" => (:gauge, :transport),
        "sliced" => (:sliced,)]
    mktempdir() do dir
        for (chapter, keys) in pages
            path = joinpath(root, "docs", "src", chapter, "index.md")
            text = read(path, String)
            @test count("<script setup>", text) == 1
            @test count("write_algorithm_player(\".\", ", text) == length(keys)
            for key in keys
                fname = basename(write_algorithm_player(dir, (; key, player = "{}")))
                @test count("import $key from './$fname'", text) == 1
                @test count("<AlgorithmPlayer :player=\"$key\" />", text) == 1
                @test count("algorithm_svg(:$key)", text) == 1
            end
        end
    end
end

# the steps of the projection, which the three testsets below read
const STEPS = ParticleRemapping._projection_steps()

# the kicks and drifts of the projection leave the box |v| ≤ 4: the window holds every stage
@testset "the window of the projection's steps" begin
    s = STEPS
    @test s.window[3] == -s.window[4]
    @test s.window[4] >= 4
    @test all(abs(v) < s.window[4] for st in s.stages for v in st.markers.v)
    vmax = maximum(abs(v) for st in s.stages for v in st.markers.v)
    @test s.window[4] ≈ max(4, 1.1 * vmax)
end

# the hat grid in panel 1, an arrow from each marker to its target in panel 2, no hole outline
@testset "the overlays of the projection's steps" begin
    s = STEPS
    @test length(s.overlays[1]) == 9 + 9
    @test all(o -> o.kind === :line, s.overlays[1])
    @test length(s.overlays[2]) == 100
    @test all(o -> o.kind === :arrow, s.overlays[2])
    @test all(isempty, s.overlays[3:7])
end

# the last stage is the remap of the deformed 10 × 10 lattice by the method of the steps
@testset "the method of the projection's steps" begin
    s = STEPS
    ex = lattice_example(; nL = 10)
    method = KickDriftRemap(;
        L = ex.L, V = ex.V, projection = HatProjection(8, 9), nξ = 16,
        nv = 20, nw = 8, maxit = 20, solver = :gauss_newton,
        density = (xs, vs) -> ex.density.(xs, vs))
    r = remap(method, ex.markers)
    @test s.stages[1].markers.x == ex.markers.x
    @test isapprox(s.stages[7].markers.x, r.markers.x; rtol = 1e-10)
    @test isapprox(s.stages[7].markers.v, r.markers.v; rtol = 1e-10)
    @test isapprox(s.stages[7].markers.w, r.markers.w; rtol = 1e-10)
end
