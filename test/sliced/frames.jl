# The states and frames of a sliced remap: remap_states and remap_frames.

using ParticleRemapping
using Random: MersenneTwister
using Test

ex = sliced_stages_example()
before, result = ex.before, ex.result
method = SlicedQuantileRemap(; B_v = 5, B_x = 12)
info = result.info
N, m, Lx = length(before), before.m, method.Lx
inbox(x) = all(-Lx .<= x .< Lx)
function same(a, b)
    a.label == b.label && a.markers.x == b.markers.x && a.markers.v == b.markers.v &&
        a.markers.w == b.markers.w
end

@testset "states" begin
    st = remap_states(method, before, result)
    @test [s.label for s in st] ==
          ["before", "Ψ_q: lift in velocity slices", "Ψ_p: lift in position slices",
        "T_μ: weight squeeze", "gauge step: new weights"]
    # the first state is the input, the last the result with its positions in the box
    @test st[1].markers === before
    last = st[end].markers
    @test inbox(last.x)
    @test maximum(abs, periodic_wrap.(last.x .- result.markers.x, Lx)) < 8 * eps(Lx)
    @test last.v == result.markers.v && last.w == result.markers.w
    # the canonical stages keep the old weights and carry the momenta of the stages
    for (k, p) in zip(2:4, (info.p̃, info.p̂, info.p̄))
        @test st[k].markers.w == before.w
        @test inbox(st[k].markers.x)
        @test isapprox(momentum(st[k].markers), p; rtol = 8 * eps())
        @test isapprox(sum(momentum(st[k].markers)), sum(momentum(before));
            atol = 8 * eps() * N * maximum(abs, p))
    end
    # the gauge step keeps (x, p) and sets the weights w_ref μ
    @test st[5].markers.x == st[4].markers.x
    @test isapprox(momentum(st[5].markers), momentum(st[4].markers); rtol = 8 * eps())
    @test st[5].markers.w ≈ info.w_ref .* info.μ rtol = 8 * eps()
    # control: the velocities at the old weights after the squeeze are not those of the output
    @test maximum(abs, st[4].markers.v .- st[5].markers.v) > 0.1
    @test_throws r"^DimensionMismatch: before has 100 markers, the result 120" remap_states(
        method, electric_wave_markers(; Nx = 10, Nv = 10), result)
end

@testset "frames along the canonical path" begin
    st = remap_states(method, before, result)
    fps, hold = 4, 2
    fr = remap_frames(method, before, result; frames_per_stage = fps, hold = hold)
    @test length(fr) == 5hold + 4fps
    # every state is held, and each stage ends exactly at its state
    for k in 1:5
        first_hold = (k - 1) * (hold + fps) + 1
        @test all(same(fr[j], st[k]) for j in first_hold:(first_hold + hold - 1))
        k > 1 && @test fr[first_hold - 1].markers.x == st[k].markers.x &&
              fr[first_hold - 1].markers.v == st[k].markers.v &&
              fr[first_hold - 1].markers.w == st[k].markers.w
    end
    @test [f.label for f in fr[(hold + 1):(hold + fps)]] == fill(st[2].label, fps)
    # halfway through Ψ_q: the lift with β_mean and β_intra halved, at the old weights
    plan = info.plan
    half = fr[hold + fps ÷ 2].markers
    q, p = split_quantile_ot_lift(before.x, momentum(before), plan.bins_v, plan.taus_x;
        anchors = plan.p_anc, β_mean = plan.β_mean_q / 2, β_intra = plan.β_intra_q / 2,
        ε = plan.ε_q)
    @test isapprox(momentum(half), p; rtol = 8 * eps())
    @test half.w == before.w
    # every frame of the canonical stages keeps the total momentum and the old weights
    canonical = [fr[(k - 1) * (hold + fps) + hold + j] for k in 1:3 for j in 1:fps]
    P = sum(momentum(before))
    @test all(f.markers.w == before.w for f in canonical)
    @test all(abs(sum(momentum(f.markers)) - P) < 64 * eps() * N for f in canonical)
    # the gauge frames keep (x, p̄) and move the weights linearly
    gauge = fr[(3 * (hold + fps) + hold + 1):(3 * (hold + fps) + hold + fps)]
    @test all(f.markers.x == st[4].markers.x for f in gauge)
    @test all(isapprox(momentum(f.markers), info.p̄; rtol = 8 * eps()) for f in gauge)
    @test gauge[2].markers.w ≈ (before.w .+ result.markers.w) ./ 2 rtol = 8 * eps()
    # control: the path of the lift is not the straight line between the two states
    line = interpolate_stages(st[1:2]; frames_per_stage = fps, hold = 0)
    @test maximum(abs, momentum(line[fps ÷ 2].markers) .- momentum(half)) >
          1e-3 * maximum(abs, momentum(before))
end

@testset "states of 2D-2V markers" begin
    rng = MersenneTwister(31)
    x = 2π .* rand(rng, 2, 36) .- π
    v = 0.8 .* randn(rng, 2, 36)
    w = 0.5 .+ rand(rng, 36)
    mk = Markers(x, v, w; m = 1.1)
    method2 = SlicedQuantileRemap(; B_v = 3, B_x = 6)
    r2 = remap(method2, mk)
    st = remap_states(method2, mk, r2)
    @test all(size(s.markers.x) == (2, 36) for s in st)
    @test st[end].markers.v == r2.markers.v
    @test isapprox(momentum(st[4].markers), r2.info.p̄; rtol = 8 * eps())
    fr = remap_frames(method2, mk, r2; frames_per_stage = 2, hold = 1)
    @test length(fr) == 5 + 8
    @test same(fr[end], st[end])
    @test fr[3 * 3 + 1].markers.v == st[4].markers.v
end
