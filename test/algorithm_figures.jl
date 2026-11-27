using ParticleRemapping
using Test

count_of(s, pattern) = length(collect(eachmatch(pattern, s)))

# the number of markers and of panels of each example
const EXAMPLES = Dict(:projection => (100, 7), :spline => (16, 5), :gauge => (8, 5),
    :transport => (48, 4), :sliced => (40, 6))

@testset "algorithm_svg($key)" for (key, (N, panels)) in EXAMPLES
    s = algorithm_svg(key)
    p, a = s.panels.svg, s.animation.svg
    @test count_of(p, r"<clipPath") == panels
    # every panel draws every marker, plus the open circles of the targets
    @test count_of(p, r"fill-opacity=\"0.85\"/>") == N * panels
    @test count_of(a, r"attributeName=\"cx\"") == N
    @test !occursin("NaN", p) && !occursin("NaN", a)
end

@test_throws r"^ArgumentError: unknown construction :other" algorithm_svg(:other)

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
