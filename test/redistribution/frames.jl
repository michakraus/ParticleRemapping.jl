# The states and the frames of the standard pictures of the redistribution: before, relaxed
# targets, symplectic step(s) at the old weights, gauge step with the KKT weights.

using ParticleRemapping
using LinearAlgebra
using Test

ex = redistribution_example()
mk = ex.markers
method = SplineRedistribution()
res = remap(method, mk)
states = remap_states(method, mk, res)
x̄, v̄, w̄ = res.info.targets_x, res.info.targets_v, res.info.target_weights

@testset "States: implicit midpoint in one step" begin
    @test [st.label for st in states] == ["before", "relaxed targets",
        "symplectic step, old weights", "gauge step, KKT weights"]
    # the first state is the input, the last the result
    @test states[1].markers.x == mk.x && states[1].markers.v == mk.v &&
          states[1].markers.w == mk.w
    @test states[end].markers === res.markers
    # the targets are computed, the markers do not move
    @test states[2].markers === mk
    # the canonical stage keeps the weights and reaches the targets in (x, p)
    step = states[3].markers
    @test step.w == mk.w
    @test all(0 .<= step.x .< ex.Lx)
    @test step.x ≈ mod.(x̄, ex.Lx) rtol = 1e-12
    @test momentum(step) ≈ w̄ .* mk.m .* v̄ rtol = 1e-12
    # the canonical stage keeps the total momentum
    @test abs(total_momentum(step) - total_momentum(mk)) < 1e-15
    # the gauge step keeps (x, p) and changes the weights and the velocities
    @test states[end].markers.x == step.x
    @test momentum(states[end].markers) ≈ momentum(step) rtol = 8 * eps()
    @test states[end].markers.w == w̄
    # control: without the gauge step the velocities are those of the old weights, far from v̄
    @test norm(step.v - v̄) > 0.1
end

@testset "States: two symplectic Euler substeps" begin
    m2 = SplineRedistribution(; scheme = SymplecticEulerStep(), substeps = 2)
    r2 = remap(m2, mk)
    s2 = remap_states(m2, mk, r2)
    @test [st.label for st in s2][3:4] ==
          ["symplectic step 1 of 2, old weights", "symplectic step 2 of 2, old weights"]
    @test length(s2) == 5
    @test all(st.markers.w == mk.w for st in s2[1:4])
    # the first substep ends half way to the targets
    z0 = vcat(mk.x, mk.w .* mk.m .* mk.v)
    z̄ = vcat(r2.info.targets_x, r2.info.target_weights .* mk.m .* r2.info.targets_v)
    zh = (z0 .+ z̄) ./ 2
    N = length(mk)
    @test s2[3].markers.x ≈ mod.(zh[1:N], ex.Lx) rtol = 1e-12
    @test momentum(s2[3].markers) ≈ zh[(N + 1):end] rtol = 1e-12
    @test s2[end].markers === r2.markers
    @test momentum(s2[end].markers) ≈ momentum(s2[4].markers) rtol = 8 * eps()
end

@testset "Frames: the pseudo-time path of the step" begin
    s = range(0, 1; length = 5)
    hold, fps = 3, 4
    frames = remap_frames(method, mk, res; s = s, hold = hold, frames_per_stage = fps)
    @test length(frames) == 2hold + (length(s) - 1) + hold + fps + hold
    @test frames[1] === states[1]
    @test frames[hold + 1] === states[2]
    @test frames[end] === states[end]
    # the path runs at the old weights and ends at the end of the step
    path = frames[(2hold + 1):(2hold + length(s) - 1)]
    @test all(f.markers.w == mk.w for f in path)
    @test all(all(0 .<= f.markers.x .< ex.Lx) for f in path)
    @test last(path).markers.x == states[3].markers.x
    @test last(path).markers.v == states[3].markers.v
    @test first(path).label == "symplectic step, old weights, s = 0.25"
    # the markers move along the path: not a straight line in (x, p)
    zs = [vcat(f.markers.x, momentum(f.markers)) for f in path]
    @test norm(zs[2] - (zs[1] + zs[3]) / 2) > 1e-6
    # the gauge frames keep (x, p)
    gauge = frames[(end - fps - hold + 1):end]
    @test all(f.markers.x == states[3].markers.x for f in gauge)
    @test all(isapprox(momentum(f.markers), momentum(states[3].markers); rtol = 1e-14)
    for f in gauge)
    @test_throws r"^ArgumentError: the pseudo-time values s must run from 0 to 1" remap_frames(
        method, mk, res; s = range(0, 0.5; length = 3))
end
