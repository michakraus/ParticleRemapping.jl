# The case P2: the markers of VlasovMethods' bump-on-tail draw on the box of a case.

function ParticleRemapping.bump_on_tail_case(; N::Integer = 10_000, seed::Integer = 1234,
        V::Real = 8.0)
    N >= 1 && seed >= 0 || throw(ArgumentError(
        "bump_on_tail_case needs N >= 1 and a seed >= 0, got N = $N, seed = $seed"))
    ε, κ, α, σ, v₀ = 0.03, 0.3, 0.1, 0.5, 4.5
    L = 2π / κ
    # VlasovMethods' draw takes its accept–reject numbers from the default generator, which is
    # seeded with `seed` here; the caller's stream is restored after the draw
    state = copy(Random.default_rng())
    particles = try
        Random.seed!(seed)
        initialize!(ParticleDistribution(1, 1, Int(N)), BumpOnTail(; ε, κ, α, σ, v₀)).particles
    finally
        copy!(Random.default_rng(), state)
    end
    markers = Markers(
        mod.(vec(particles.x), L) .- L / 2, vec(particles.v), vec(particles.w);
        m = 1.0)
    function density(x, v)
        xs = x + L / 2  # the position on VlasovMethods' interval [0, L)
        (1 - ε * cos(κ * xs)) *
        ((1 - α) * exp(-v^2 / 2) + α / σ * exp(-(v - v₀)^2 / (2σ^2))) / sqrt(2π)
    end
    parameters = (;
        N, seed, V = Float64(V), ε, κ, α, σ, v₀, cells = 16, degree = 3, dt = 0.1,
        steps = 250, window = (5.0, 17.5), rate = nothing)
    StandardCase(:P2, parameters, markers, L, V, density, density, :wave)
end
