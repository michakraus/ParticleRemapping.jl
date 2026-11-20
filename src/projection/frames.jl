# The data of the animations and of the pictures of the algorithm steps: computed here, drawn by the
# Makie extension.

# The frames of one remap of the markers of `prob`: the two stages of `par` with `θ` scaled from
# zero to one, then the reweighting from the old weights to `wt` at the mapped state `(X, P)`,
# then, for `wc !== nothing`, the charge-conserving reweighting from `wt` to `wc`.
function _remap_frames(prob, par, θ, X, P, wt, wc, label; n, nreweight)
    m = prob.m
    h = length(θ) ÷ 2
    frames = Tuple{Vector{Float64}, Vector{Float64}, String}[(prob.x0, prob.v0, "before")]
    for s in 1:2, t in range(0, 1; length = n + 1)[2:end]

        θt = s == 1 ? vcat(t .* θ[1:h], zero(θ[(h + 1):end])) :
             vcat(θ[1:h], t .* θ[(h + 1):end])
        Xt, Pt = remap_stages(par, θt, prob.x0, prob.p0, prob.w0; m)
        push!(frames, (Xt, Pt ./ (m .* prob.w0), "stage $s$label"))
    end
    for t in range(0, 1; length = nreweight + 1)[2:end]
        push!(frames, (X, P ./ (m .* prob.w0 .^ (1 - t) .* wt .^ t),
            "reweighting to the target weights"))
    end
    wc === nothing && return frames
    for t in range(0, 1; length = nreweight + 1)[2:end]
        push!(frames, (X, P ./ (m .* wt .^ (1 - t) .* wc .^ t),
            "charge-conserving reweighting"))
    end
    frames
end

"""
    remap_animation_frames(experiment, variant; n = 30, nreweight = 20) -> Vector of (x, v, title)
    remap_animation_frames(problem, result; n = 30, nreweight = 20) -> Vector of (x, v, title)

The frames of the animation of one remap: the markers before the map, then each of the two stages
with its parameters scaled from zero to one in `n` frames (velocities `p/(m w)` with the old
weights), then the reweighting to the target weights and the charge-conserving reweighting, each
interpolated geometrically in the weights in `nreweight` frames.

The first form animates the remap `variant` of [`hole_filling_experiment`](@ref). The second form
animates the [`RemapResult`](@ref ParticleRemapping.RemapResult) `result` of [`remap`](@ref
ParticleRemapping.remap) with a [`KickDriftRemap`](@ref) or a [`NonSeparableRemap`](@ref) on the
markers of `problem`, a [`hole_filling_problem`](@ref); without a charge-conserving reweighting
its frames end with the reweighting to the target weights.
"""
function remap_animation_frames(experiment, variant::Symbol; n = 30, nreweight = 20)
    prob = experiment.problem
    r = experiment.results[variant]
    _remap_frames(prob, r.parametrisation, r.θ, r.X, r.P, experiment.targets.w, r.wc,
        ", $variant"; n, nreweight)
end

function remap_animation_frames(problem, result::RemapResult; n = 30, nreweight = 20)
    i = result.info
    _remap_frames(problem, i.parametrisation, i.θ, i.X, i.P, i.targets.w,
        get(i, :wc, nothing), ""; n, nreweight)
end

"""
    moser_flow_history(problem = hole_filling_problem(); projection = HatProjection(32, 33),
        steps = 40) -> Vector of (x, v)

The markers of `problem` after each Runge–Kutta step of the density-interpolation flow
([`moser_targets`](@ref)), in `(x, v)`, starting with the markers themselves.
"""
function moser_flow_history(
        problem = hole_filling_problem(); projection = HatProjection(32, 33),
        steps = 40)
    history = Tuple{Vector{Float64}, Vector{Float64}}[]
    moser_targets(
        projection, problem.x0 ./ problem.L, problem.v0 ./ problem.V; steps, history)
    [(problem.L .* s, problem.V .* r) for (s, r) in history]
end

"""
    remap_stage_states(experiment, variant) -> Vector of (x, v, title)

The markers of [`hole_filling_experiment`](@ref) before the map, after each kick and each drift (or
each non-separable stage), with the velocities `p/(m w)` of the old weights, then after the
reweighting with the target weights and after the charge-conserving reweighting.
"""
function remap_stage_states(experiment, variant::Symbol)
    prob = experiment.problem
    r = experiment.results[variant]
    m = prob.m
    history = Tuple{Vector{Float64}, Vector{Float64}}[]
    remap_stages(r.parametrisation, r.θ, prob.x0, prob.p0, prob.w0; m, history)
    names = r.parametrisation isa KickDriftParametrisation ?
            [isodd(k) ? "kick $(cld(k, 2))" : "drift $(cld(k, 2))"
             for k in eachindex(history)] :
            ["stage $k" for k in eachindex(history)]
    [(prob.x0, prob.v0, "before");
     [(X, P ./ (m .* prob.w0), name) for ((X, P), name) in zip(history, names)];
     (r.X, r.v1, "reweighting");
     (r.X, r.vc, "charge-conserving")]
end

"""
    dlra_push_frames(experiment; tol = 1e-4, rmin = 4, nsub = 5, rng = Random.default_rng())
        -> Vector of (; F, x, v, hole_x, hole_v, t, rank)

The frames of the animation of [`dlra_weights_experiment`](@ref): the low-rank density at the
tolerance `tol` on the grid, the markers and the outline of the hole, every `nsub` steps of the
BUG integrator and of the characteristics up to the time of the push.
"""
function dlra_push_frames(experiment; tol = 1e-4, rmin = 4, nsub = 5,
        rng::AbstractRNG = Random.default_rng())
    g, E, dt, prob = experiment.grid, experiment.E, experiment.dt, experiment.problem
    f = lowrank_factorise(g, prob.density.(g.x, g.v'); tol, rmin, rng)
    x, v = prob.x0, prob.v0
    hx, hv = hole_outline(prob.hole)
    frames = [(; F = lowrank_full(f), x, v, hole_x = hx, hole_v = hv, t = 0.0,
        rank = lowrank_rank(f))]
    t = 0.0
    for _ in 1:round(Int, experiment.T / (nsub * dt))
        rn = 0
        for _ in 1:nsub
            f, rn = bug_step(g, f, dt; tol, rmin)
        end
        x, v = vlasov_characteristics(E, x, v, nsub * dt)
        hx, hv = vlasov_characteristics(E, hx, hv, nsub * dt)
        t += nsub * dt
        push!(frames, (; F = lowrank_full(f), x, v, hole_x = hx, hole_v = hv, t, rank = rn))
    end
    frames
end
