using ParticleRemapping
using Test

count_of(s, pattern) = length(collect(eachmatch(pattern, s)))

# the number of markers and of panels of each example
const EXAMPLES = Dict(:projection => (96, 7), :spline => (16, 5), :gauge => (8, 5),
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
