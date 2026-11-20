using ParticleRemapping
using Test

@testset "interpolate_stages" begin
    m = 2.0
    a = Markers([0.1, 3.0, 6.1], [1.0, -0.5, 2.0], [0.5, 1.0, 2.0]; m)
    # a canonical stage: positions and momenta move, the weights stay
    b = Markers([0.4, 2.0, 0.2], [0.5, -1.0, 2.5], copy(a.w); m)
    # a gauge step: positions and momenta stay, the weights change
    c = gauge_step(b, [2.0, 0.5, 1.0])
    stages = [RemapStage(a, "a"), RemapStage(b, "b"), RemapStage(c, "c")]
    frames = interpolate_stages(stages; frames_per_stage = 4, hold = 2, period = (0, 2π))
    @test length(frames) == 2 + 2 * (4 + 2)
    @test frames[1].markers === a
    @test frames[end].markers === c
    @test [f.label for f in frames[3:6]] == fill("b", 4)
    # the frames end on the stages
    @test frames[6].markers.x ≈ b.x
    @test frames[6].markers.v ≈ b.v
    # the gauge step keeps x and p in every frame
    for f in frames[9:12]
        @test f.markers.x ≈ b.x
        @test momentum(f.markers) ≈ momentum(b)
    end
    # periodic positions move by the shorter way round: 6.1 → 0.2 passes 2π, not 3
    mid = interpolate_stages(stages[1:2]; frames_per_stage = 2, hold = 0, period = (0, 2π))
    @test mid[1].markers.x[3] ≈ mod(6.1 + (0.2 + 2π - 6.1) / 2, 2π)
    @test all(0 .<= mid[1].markers.x .< 2π)
    # without a period the same marker crosses the whole interval (the control)
    flat = interpolate_stages(stages[1:2]; frames_per_stage = 2, hold = 0)
    @test flat[1].markers.x[3] ≈ (6.1 + 0.2) / 2
    @test_throws r"^ArgumentError: no stages to interpolate" interpolate_stages(RemapStage[])
    @test_throws r"^DimensionMismatch: stage 1 has 3 markers, stage 2 has 1" interpolate_stages(
        [stages[1], RemapStage(Markers([0.0], [0.0], [1.0]), "d")])
end
