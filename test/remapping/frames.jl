# The states and animation frames of the gauge rezoning and of the canonical transport.

using ParticleRemapping
using Test

@testset "Frames" begin
    @testset "gauge rezoning: states along the homotopy" begin
        ex = eight_marker_example()
        mk, grid, method = ex.markers, ex.grid, ex.method
        r = remap(method, mk)
        relres(after, k) = maximum(abs.(grid_moment(grid, after, k) .-
                                        grid_moment(grid, mk, k)) ./
                                   grid_moment(grid, mk.x, abs.(mk.v), mk.w, k))
        # the homotopy accepts τ = 0.25, 0.5, 0.75, 1; the two steps nearest τ = 1/3 and 2/3
        s = remap_states(method, mk, r)
        @test [st.label for st in s] ==
              ["before", "homotopy, τ = 0.25", "homotopy, τ = 0.75",
            "gauge step: new weights, τ = 1.0"]
        @test s[1].markers === mk
        @test s[end].markers === r.markers
        for st in s
            # a gauge step keeps the positions and the momenta
            @test st.markers.x == mk.x
            @test momentum(st.markers) ≈ momentum(mk) rtol = 8eps()
        end
        # an accepted point of the homotopy meets the charge and energy constraints
        for st in s[2:3], k in (0, 2)

            @test relres(st.markers, k) < 1e-13
        end
        # control: the intermediate states are not the input
        @test maximum(abs.(s[2].markers.w .- mk.w)) > 1e-3
        @test maximum(abs.(s[3].markers.w .- s[2].markers.w)) > 1e-3
        @test length(remap_states(method, mk, r; steps = 5)) == 2 + 3   # three inner points
        @test_throws r"^ArgumentError: the result does not belong to these markers" remap_states(
            method, Markers(mk.x .+ 0.1, mk.v, mk.w; m = mk.m), r)
    end

    @testset "gauge rezoning: frames along the homotopy" begin
        ex = eight_marker_example()
        mk, grid, method = ex.markers, ex.grid, ex.method
        r = remap(method, mk)
        relres(after, k) = maximum(abs.(grid_moment(grid, after, k) .-
                                        grid_moment(grid, mk, k)) ./
                                   grid_moment(grid, mk.x, abs.(mk.v), mk.w, k))
        frames = remap_frames(method, mk, r; frames_per_step = 3, hold = 2)
        @test length(frames) == 2 + 4 * 3 + 2
        @test frames[1].markers === mk
        @test frames[end].markers === r.markers
        @test frames[end].label == "gauge step: new weights, τ = 1.0"
        @test frames[2 + 3].label == "homotopy, τ = 0.25"
        for f in frames
            @test f.markers.x == mk.x
            @test momentum(f.markers) ≈ momentum(mk) rtol = 8eps()
        end
        # the last frame of each homotopy step is the accepted point, the last one the result
        @test relres(frames[2 + 3].markers, 0) < 1e-13
        @test frames[2 + 4 * 3].markers.w ≈ r.markers.w rtol = 8eps()
        # control: between two accepted points the stretches are interpolated, not solved, and
        # the constraints do not hold there
        @test relres(frames[2 + 3 + 2].markers, 0) > 1e-6
        @test_throws r"^ArgumentError: frames_per_step must be positive, got 0" remap_frames(
            method, mk, r; frames_per_step = 0)
    end

    @testset "canonical transport" begin
        setup = quantile_matching_setup(; N = 4000, Ng = 8)
        mk, grid, method = setup.markers, setup.grid, setup.admissible
        r = remap(method, mk)
        s = remap_states(method, mk, r)
        @test [st.label for st in s] ==
              ["before", "ψ at the old weights", "gauge step: matched weights"]
        @test s[1].markers === mk
        @test s[end].markers === r.markers
        # the canonical map keeps the weights and gives (X, P)
        @test s[2].markers.w == mk.w
        @test s[2].markers.x ≈ mod.(r.info.X, grid.L)
        @test momentum(s[2].markers) ≈ r.info.P rtol = 8eps()
        # the matching keeps (X, P) of the canonical map
        @test s[3].markers.x == s[2].markers.x
        @test momentum(s[3].markers) ≈ momentum(s[2].markers) rtol = 8eps()
        # control: at the old weights the moved markers misrepresent f, after the matching not
        err(m) = binned_weight_error(setup.target, grid, m.x, m.v, m.w)
        @test err(s[2].markers) > 3err(s[3].markers)
        frames = remap_frames(method, mk, r; frames_per_stage = 4, hold = 1)
        @test length(frames) == 1 + 2 * (4 + 1)
        @test frames[end].markers === r.markers
        @test all(all(0 .<= f.markers.x .< grid.L) for f in frames)
        @test_throws r"^ArgumentError: the result holds 4000 markers, the input 3" remap_states(
            method, Markers(mk.x[1:3], mk.v[1:3], mk.w[1:3]), r)
    end
end
