# The states and frames of the remaps of the Hamiltonian particle projection, and the data of the
# other animations, on small runs.

using ParticleRemapping
using Random
using Test

Random.seed!(20261126)

prob = hole_filling_problem(; nL = 16, hole_a = 0.9, hole_b = 0.8)
mk = prob.markers
hp = HatProjection(12, 13)
L, m = prob.L, prob.m

# two stages with the same markers and label
function same(a, b)
    a.label == b.label && a.markers.x == b.markers.x && a.markers.v == b.markers.v &&
        a.markers.w == b.markers.w
end

@testset "Frames" begin
    methods = (
        KickDriftRemap(; nξ = 8, nv = 10, nw = 4, projection = hp, maxit = 2),
        NonSeparableRemap(; nξ = 6, nv = 7, nℓ = 3, projection = hp, maxit = 2),
        # 8 splines: on 16, the cubic B-splines at the 16 marker columns miss the alternating mode
        KickDriftRemap(; nξ = 8, nv = 10, nw = 4, projection = hp, maxit = 2, stages = 3,
            charge = ChargeConservingReweighting(L, 8)))
    labels = (
        ["before", "kick 1", "drift 1", "kick 2", "drift 2", "gauge: target weights"],
        ["before", "midpoint stage 1", "midpoint stage 2", "gauge: target weights"],
        ["before", "kick 1", "drift 1", "kick 2", "drift 2", "kick 3", "drift 3",
            "gauge: target weights", "gauge: charge-conserving"])
    @testset "remap states, $(nameof(typeof(method))) $k" for (k, method) in enumerate(methods)
        r = remap(method, mk)
        i = r.info
        states = remap_states(method, mk, r)
        ngauge = haskey(i, :wc) ? 2 : 1
        ψ = states[2:(end - ngauge)]                       # the stages of the canonical map
        @test [st.label for st in states] == labels[k]
        # the first state is the input, the last the result
        @test states[1].markers.x == mk.x && states[1].markers.v == mk.v &&
              states[1].markers.w == mk.w
        @test states[end].markers.x == r.markers.x &&
              states[end].markers.v == r.markers.v &&
              states[end].markers.w == r.markers.w
        # every stage of ψ moves the markers at the old weights and keeps P
        for (a, b) in zip(states[1:(end - ngauge - 1)], ψ)
            @test b.markers.w == mk.w
            @test b.markers.v != a.markers.v || b.markers.x != a.markers.x
            @test sum(momentum(b.markers)) ≈ sum(prob.p0) atol = 8eps() * sum(abs, prob.p0)
            @test all(0 .<= b.markers.x .< L)
        end
        # a kick keeps the positions, a drift the velocities
        if method isa KickDriftRemap
            @test states[2].markers.x == mk.x
            @test states[3].markers.v == states[2].markers.v
        end
        # the stages end at the image of the map with the fitted θ
        @test ψ[end].markers.x ≈ mod.(i.X, L) rtol = 8eps()
        @test momentum(ψ[end].markers) ≈ i.P rtol = 8eps()
        # each gauge step keeps (x, p) and changes the weights
        for (a, b) in zip(states[(end - ngauge):(end - 1)], states[(end - ngauge + 1):end])
            @test b.markers.x == a.markers.x
            @test momentum(b.markers) ≈ momentum(a.markers) rtol = 8eps()
            @test b.markers.w != a.markers.w
        end
        @test states[end - ngauge + 1].markers.w == i.targets.w
        haskey(i, :wc) && @test i.charge_converged && states[end].markers.w == i.wc
        # the case that must fail: markers that are not the input of the remap
        moved = Markers(mk.x .+ 0.1, mk.v, mk.w; m)
        @test_throws r"^ArgumentError: the markers are not the input of the remap" remap_states(
            method, moved, r)
    end

    @testset "remap frames, $(nameof(typeof(method))) $k" for (k, method) in enumerate(methods)
        r = remap(method, mk)
        states = remap_states(method, mk, r)
        nf, hold = 3, 2
        frames = remap_frames(method, mk, r; frames_per_stage = nf, hold)
        @test length(frames) == hold + (length(states) - 1) * (nf + hold)
        @test same(frames[1], states[1]) && same(frames[end], states[end])
        # the path passes through every state: the frame at the end of a sub-step is the state
        for j in 2:length(states)
            e = hold + (j - 1) * (nf + hold)
            @test same(frames[e], states[j])
            @test frames[e - hold].markers.x ≈ states[j].markers.x rtol = 8eps()
            @test frames[e - hold].markers.v ≈ states[j].markers.v rtol = 8eps()
            @test frames[e - hold].label == states[j].label
        end
        # the canonical frames keep the weights, the gauge frames keep (x, p)
        nψ = length(states) - (haskey(r.info, :wc) ? 2 : 1)
        canonical = frames[1:(hold + (nψ - 1) * (nf + hold))]
        @test all(f.markers.w == mk.w for f in canonical)
        gauge = frames[(length(canonical) + 1):end]
        @test all(f.markers.x == states[nψ].markers.x for f in gauge)
        @test all(isapprox(momentum(f.markers), momentum(states[nψ].markers); rtol = 8eps())
        for f in gauge)
        # the frames inside a sub-step lie between its states, at a fraction of its parameters
        @test frames[hold + 1].markers.v != states[1].markers.v
        @test frames[hold + 1].markers.v != states[2].markers.v
    end

    @testset "charge-conserving reweighting" begin
        r0 = remap(KickDriftRemap(; nξ = 8, nv = 10, nw = 4, projection = hp, maxit = 2), mk)
        method = ChargeConservingReweighting(L, 8)
        # a reference charge that positive weights at these positions reproduce
        x0 = r0.markers.x
        ρ = spline_grid_charge(method.basis, x0, r0.markers.w .* exp.(0.1 .* cos.(x0)))
        r = remap(method, r0.markers; reference_charge = ρ)
        states = remap_states(method, r0.markers, r)
        @test [st.label for st in states] == ["before", "gauge: charge-conserving"]
        @test states[1].markers === r0.markers && states[end].markers === r.markers
        @test states[2].markers.x == states[1].markers.x
        @test momentum(states[2].markers) ≈ momentum(states[1].markers) rtol = 8eps()
        @test_throws r"^ArgumentError: the markers are not the input of the remap" remap_states(
            method, Markers(r0.markers.x .+ 0.1, r0.markers.v, r0.markers.w; m), r)
    end

    @testset "remaps of the hole-filling experiment" begin
        e = hole_filling_experiment(; problem = prob, projection = hp,
            variants = (:weight_dependent, :non_separable), maxit = 2)
        for (v, n) in ((:weight_dependent, 4), (:non_separable, 2))
            res = e.results[v]
            method, before, result = experiment_remap(e, v)
            @test method isa (v === :weight_dependent ? KickDriftRemap : NonSeparableRemap)
            @test before === prob.markers
            states = remap_states(method, before, result)
            @test length(states) == 1 + n + 2
            @test states[end - 1].markers.v ≈ res.v1 rtol = 8eps()
            @test states[end - 1].markers.w == e.targets.w
            @test states[end].markers.v == res.vc && states[end].markers.w == res.wc
            @test states[end].markers.x ≈ mod.(res.X, L) rtol = 8eps()
        end
    end

    @testset "density-interpolation flow" begin
        hist = moser_flow_history(prob; projection = hp, steps = 4)
        @test length(hist) == 5
        @test [h.label for h in hist[[1, end]]] ==
              ["density-interpolation flow, t = 0.0", "density-interpolation flow, t = 1.0"]
        @test hist[1].markers.x ≈ prob.x0 rtol = 8eps()
        @test all(h.markers.w == prob.w0 for h in hist)
    end

    @testset "DLRA push" begin
        d = dlra_weights_experiment(; problem = prob, fit = false, tolerances = (1e-4,))
        frames = dlra_push_frames(d; nsub = 50)
        @test length(frames) == 1 + round(Int, d.T / (50 * d.dt))
        @test frames[end].t ≈ d.T
        @test frames[end].x ≈ d.xT rtol = 1e-12
        @test all(f.rank >= 1 for f in frames)
        @test size(frames[end].F) == (length(d.grid.x), length(d.grid.v))
    end
end
