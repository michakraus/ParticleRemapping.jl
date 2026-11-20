# The set-ups of the two numerical examples of the remapping note.

"""
    two_stream_density(x, v; amplitude = 0.3, beam = 2)

The two-stream density `(1 + amplitude cos x) (e^{−(v−beam)²/2} + e^{−(v+beam)²/2}) / (2√(2π))`
of the weights of [`gauge_rezoning_setup`](@ref).
"""
function two_stream_density(x::Real, v::Real; amplitude::Real = 0.3, beam::Real = 2)
    (1 + amplitude * cos(x)) * 0.5 * (exp(-(v - beam)^2 / 2) + exp(-(v + beam)^2 / 2)) /
    sqrt(2π)
end

"""
    gauge_rezoning_setup(; Ng = 16, markers_per_node = 400, L = 2π, m = 1, ncentres = 6,
                         vmax = 4, spread = 0.08, weight_floor = 1e-6, seed = 20261119,
                         K = (0, 2), dτ0 = 0.05, dτmin = 1e-3, solver = :kkt)

The numerical example of the gauge rezoning: `N = Ng · markers_per_node` markers with uniform
positions on `[0, L)`; in each cell the velocities cluster around `ncentres` random centres in
`[−vmax, vmax]` with the spread `spread`; the weights are
`w_a = f(x_a, v_a) L 2vmax / N + weight_floor` with the [`two_stream_density`](@ref) `f`. The
random numbers come from `Xoshiro(seed)`. Returns `(; markers, grid, method)` with the
[`GaugeRezoning`](@ref) of the example, which solves the full sparse KKT system.
"""
function gauge_rezoning_setup(; Ng::Integer = 16, markers_per_node::Integer = 400,
        L::Real = 2π, m::Real = 1.0, ncentres::Integer = 6, vmax::Real = 4.0,
        spread::Real = 0.08, weight_floor::Real = 1e-6, seed::Integer = 20261119,
        K::Tuple = (0, 2), dτ0::Real = 0.05, dτmin::Real = 1e-3, solver::Symbol = :kkt)
    rng = Xoshiro(seed)
    grid = PeriodicGrid(L, Ng)
    N = Ng * markers_per_node
    x = L .* rand(rng, N)
    centres = [sort(2vmax .* rand(rng, ncentres) .- vmax) for _ in 1:Ng]
    v = [centres[cellof(grid, xa)][rand(rng, 1:ncentres)] + spread * randn(rng) for xa in x]
    w = [two_stream_density(x[a], v[a]) * L * 2vmax / N + weight_floor for a in 1:N]
    method = GaugeRezoning(grid; K, s = 1.0, dτ0, dτmin, solver)
    (; markers = Markers(x, v, w; m), grid, method)
end

"""
    quantile_matching_setup(; N = 200_000, α = 0.4, L = 2π, Ng = 16, m = 1, t = 1.5,
                            dt = 1e-2, seed = 20261119)

The numerical example of the canonical transport: `N` markers with equal weights `w0 = L/N`
sampled from `f(x, v) = (1 + α cos x) e^{−v²/2} / √(2π)` on `[0, L)` (positions by rejection,
velocities normal, from `Xoshiro(seed)`), `Ng` cells and the mass `m`. Returns
`(; markers, grid, target, α, t, identity, admissible, inadmissible, a_admissible,
a_inadmissible)`:
the [`SeparableMaxwellianTarget`](@ref), and three [`CurrentQuantileMatching`](@ref) methods with
`ψ = id`, the admissible cotangent lift of `a = J^{−1/2}` and the inadmissible lift of
`a = 1 + 0.5 cos x` at time `t` (Runge–Kutta step `dt`); `a_…` are the pairs `(a, a')`.
"""
function quantile_matching_setup(; N::Integer = 200_000, α::Real = 0.4, L::Real = 2π,
        Ng::Integer = 16, m::Real = 1.0, t::Real = 1.5, dt::Real = 1e-2,
        seed::Integer = 20261119)
    rng = Xoshiro(seed)
    x = Float64[]
    while length(x) < N
        u = L * rand(rng)
        rand(rng) * (1 + α) < cosine_profile(u, α) && push!(x, u)
    end
    v = randn(rng, N)
    markers = Markers(x, v, fill(L / N, N); m)
    grid = PeriodicGrid(L, Ng)
    target = SeparableMaxwellianTarget(grid, α)
    J(x) = cosine_profile(x, α)
    a_ok = (x -> J(x)^(-1 / 2), x -> -0.5 * J(x)^(-3 / 2) * (-α * sin(x)))
    a_bad = (x -> 1 + 0.5cos(x), x -> -0.5sin(x))
    lift(a) = (x, p) -> cotangent_lift(a..., x, p, t; L, dt)
    (; markers, grid, target, α, t,
        identity = CurrentQuantileMatching(grid, (x, p) -> (x, p), target),
        admissible = CurrentQuantileMatching(grid, lift(a_ok), target),
        inadmissible = CurrentQuantileMatching(grid, lift(a_bad), target),
        a_admissible = a_ok, a_inadmissible = a_bad)
end

"""
    quantile_matching_experiment(setup = quantile_matching_setup())

The configurations of the numerical example of the canonical transport: the `original` markers,
the matching after `ψ = id` (`identity`), after the `admissible` lift and after the
`inadmissible` lift (each a `RemapResult`), and `kept`, the markers after the admissible lift
with their old weights, `W = w`, `V = P/(m w)`.
"""
function quantile_matching_experiment(setup = quantile_matching_setup())
    original = setup.markers
    identity = remap(setup.identity, original)
    admissible = remap(setup.admissible, original)
    inadmissible = remap(setup.inadmissible, original)
    kept = from_canonical(admissible.info.X, admissible.info.P, original.w; m = original.m)
    (; original, identity, admissible, inadmissible, kept)
end

"""
    lift_matching_frames(setup = quantile_matching_setup(); nframes = 75, substeps = 4,
                         every = 50)

The admissible cotangent lift of [`quantile_matching_setup`](@ref) from `0` to `setup.t` in
`nframes` frames of `substeps` Runge–Kutta steps each, with the current-quantile matching applied
at every frame. Each frame holds the time `t`, every `every`-th marker (`X`, `V`, `W`), and the
marker counts `counts` and the summed weights `charges` of the cells; the first frame is the
initial state.
"""
function lift_matching_frames(setup = quantile_matching_setup(); nframes::Integer = 75,
        substeps::Integer = 4, every::Integer = 50)
    mk, grid = setup.markers, setup.grid
    a, da = setup.a_admissible
    dt = setup.t / nframes / substeps
    p = momentum(mk)
    sub = 1:every:length(mk)
    function frame(t, X, V, W)
        cells = cellof.(Ref(grid), X)
        (; t, X = X[sub], V = V[sub], W = W[sub],
            counts = [count(==(c), cells) for c in 1:(grid.Ng)],
            charges = [sum(W[cells .== c]) for c in 1:(grid.Ng)])
    end
    Xa, ξ = copy(mk.x), ones(length(mk))
    frames = [frame(0.0, mod.(Xa, grid.L), mk.v, mk.w)]
    for i in 1:nframes
        for k in eachindex(Xa)
            x, z = Xa[k], ξ[k]
            for _ in 1:substeps
                x, z = _rk4_lift(a, da, x, z, dt)
            end
            Xa[k], ξ[k] = x, z
        end
        Xm = mod.(Xa, grid.L)
        V, W, _ = current_quantile_match(grid, Xm, p ./ ξ, setup.target; m = mk.m)
        push!(frames, frame(i * setup.t / nframes, Xm, V, W))
    end
    frames
end
