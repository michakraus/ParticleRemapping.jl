# The states and frames of a remap of the Hamiltonian particle projection, and the data of the
# other animations of the chapter: computed here, drawn by the Makie extension.

const _ProjectionRemap = Union{KickDriftRemap, NonSeparableRemap}

# The sub-steps of the map `par`: the range of the parameters of each sub-step in `θ` and its label.
# A kick–drift stage has two sub-steps, a kick (the coefficients λ) and a drift (μ); a non-separable
# stage has one.
function _substeps(par::KickDriftParametrisation)
    nk, nd = length(par.Bξ) * _nw(par), length(par.Bv) * _nw(par)
    reduce(vcat,
        [[((s - 1) * (nk + nd)) .+ (1:nk) => "kick $s",
             ((s - 1) * (nk + nd) + nk) .+ (1:nd) => "drift $s"] for s in 1:(par.stages)])
end

function _substeps(par::NonSeparableParametrisation)
    n = _stage_count(par)
    [((s - 1) * n) .+ (1:n) => "midpoint stage $s" for s in 1:(par.stages)]
end

# The markers with the lifted positions `X` (wrapped into [0, L)), the momenta `P` and the weights
# `w`, with the velocities `P / (m w)`.
_stage_markers(X, P, w, m, L) = Markers(mod.(X, L), P ./ (m .* w), w; m)

"""
    remap_states(method::Union{KickDriftRemap, NonSeparableRemap}, before, result)
    remap_states(method::ChargeConservingReweighting, before, result)

The states of a remap of the Hamiltonian particle projection, as a vector of
[`RemapStage`](@ref ParticleRemapping.RemapStage)s: the markers `before` the remap; after each
kick and each drift of [`KickDriftRemap`](@ref) (labels `"kick 1"`, `"drift 1"`, …) or after each
implicit-midpoint stage of [`NonSeparableRemap`](@ref) (`"midpoint stage 1"`, …), at the old
weights with the velocities `p / (m w)`; after the gauge step to the target weights; and, when
`result` holds the charge-conserving weights `wc`, after the gauge step to these weights. The last
state is `result.markers`. Positions are wrapped into `[0, L)`.

For a [`ChargeConservingReweighting`](@ref) the states are the markers before it and after it.

The stages are recomputed from the parameters `θ` in `result.info`; markers that are not the input
of the remap throw an `ArgumentError`.
"""
function ParticleRemapping.remap_states(
        method::_ProjectionRemap, before::Markers, result::RemapResult)
    i = result.info
    par, L, m, w = i.parametrisation, method.L, before.m, before.w
    history = Tuple{Vector{Float64}, Vector{Float64}}[]
    remap_stages(par, i.θ, before.x, momentum(before), w; m, history)
    isapprox(history[end][1], i.X; rtol = 1e-12) || throw(ArgumentError(
        "the markers are not the input of the remap: the stages do not reach the positions of `result`"))
    states = [RemapStage(before, "before");
              [RemapStage(_stage_markers(X, P, w, m, L), label)
               for ((X, P), (_, label)) in zip(history, _substeps(par))]]
    if haskey(i, :wc)
        push!(states,
            RemapStage(_stage_markers(i.X, i.P, i.targets.w, m, L), "gauge: target weights"))
        push!(states, RemapStage(result.markers, "gauge: charge-conserving"))
    else
        push!(states, RemapStage(result.markers, "gauge: target weights"))
    end
    states
end

function ParticleRemapping.remap_states(
        ::ChargeConservingReweighting, before::Markers, result::RemapResult)
    before.x == result.markers.x || throw(ArgumentError(
        "the markers are not the input of the remap: the reweighting keeps the positions"))
    [RemapStage(before, "before"), RemapStage(result.markers, "gauge: charge-conserving")]
end

"""
    remap_frames(method::Union{KickDriftRemap, NonSeparableRemap}, before, result;
        frames_per_stage = 30, hold = 8)

The frames of an animation of the remap along its natural path: the parameters of each kick and
each drift, or of each non-separable stage, scaled from zero to their fitted values in
`frames_per_stage` frames, one sub-step after the other, at the old weights; then the gauge steps
of [`remap_states`](@ref ParticleRemapping.remap_states) at fixed `(x, p)`, interpolated by
[`interpolate_stages`](@ref ParticleRemapping.interpolate_stages). Each state is held for `hold`
frames. A sub-step with its parameters at zero is the identity, so the path passes through every
state of `remap_states`.
"""
function ParticleRemapping.remap_frames(method::_ProjectionRemap, before::Markers,
        result::RemapResult; frames_per_stage::Integer = 30, hold::Integer = 8)
    states = remap_states(method, before, result)
    i = result.info
    par, L, m, w, p = i.parametrisation, method.L, before.m, before.w, momentum(before)
    steps = _substeps(par)
    frames = fill(states[1], hold)
    θt = zero(i.θ)
    for (k, (r, label)) in enumerate(steps)
        for t in range(0, 1; length = frames_per_stage + 1)[2:end]
            θt[r] .= t .* i.θ[r]
            X, P = remap_stages(par, θt, before.x, p, w; m)
            push!(frames, RemapStage(_stage_markers(X, P, w, m, L), label))
        end
        append!(frames, fill(states[k + 1], hold))
    end
    gauge = interpolate_stages(states[(length(steps) + 1):end]; frames_per_stage, hold)
    append!(frames, gauge[(hold + 1):end])
end

"""
    experiment_remap(experiment, variant) -> (method, markers, result)

The remap `variant` of a [`hole_filling_experiment`](@ref) run as the arguments of
[`remap_states`](@ref ParticleRemapping.remap_states) and
[`remap_frames`](@ref ParticleRemapping.remap_frames): the method of the variant, the markers of
the problem, and a [`RemapResult`](@ref ParticleRemapping.RemapResult) with the markers after the
charge-conserving reweighting and the parameters, the lifted image `(X, P)`, the targets and the
charge-conserving weights `wc` of the run.
"""
function experiment_remap(experiment, variant::Symbol)
    prob, r = experiment.problem, experiment.results[variant]
    info = (; r.θ, r.X, r.P, targets = experiment.targets, r.parametrisation, r.wc)
    r.method, prob.markers,
    RemapResult(Markers(mod.(r.X, prob.L), r.vc, r.wc; prob.m), info)
end

"""
    moser_flow_history(problem = hole_filling_problem(); projection = HatProjection(32, 33),
        steps = 40) -> Vector{RemapStage}

The markers of `problem` after each Runge–Kutta step of the density-interpolation flow
([`moser_targets`](@ref)), starting with the markers themselves, at their old weights and labelled
by the flow time `t`. The flow moves the points `(x, v)` and changes volume; it only provides the
targets.
"""
function moser_flow_history(
        problem = hole_filling_problem(); projection = HatProjection(32, 33),
        steps = 40)
    history = Tuple{Vector{Float64}, Vector{Float64}}[]
    moser_targets(
        projection, problem.x0 ./ problem.L, problem.v0 ./ problem.V; steps, history)
    n = length(history) - 1
    [RemapStage(Markers(problem.L .* s, problem.V .* r, problem.w0; problem.m),
         "density-interpolation flow, t = $(round((k - 1) / n; digits = 3))")
     for (k, (s, r)) in enumerate(history)]
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
