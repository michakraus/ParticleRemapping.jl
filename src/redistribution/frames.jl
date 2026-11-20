# The states and the frames of the standard pictures of one redistribution: the markers before,
# with the relaxed targets, after each symplectic step at the old weights and after the gauge step
# with the KKT weights.

# The packed state `(q, p)` of the markers, with the momenta computed as `remap` computes them.
_packed_state(mk::Markers) = vcat(mk.x, mk.w .* mk.m .* mk.v)

# The markers of the packed state `z` at the weights `w`, positions wrapped into `[0, Lx)`.
function _step_markers(z::AbstractVector, w::AbstractVector, m::Real, Lx::Real)
    N = length(w)
    Markers(mod.(z[1:N], Lx), z[(N + 1):(2N)] ./ (m .* w), w; m = m)
end

function _step_label(k, K)
    K == 1 ? "symplectic step, old weights" : "symplectic step $k of $K, old weights"
end

"""
    remap_states(method::SplineRedistribution, before::Markers, result::RemapResult)
        -> Vector{RemapStage}

The states of the redistribution of `before` that gave `result`: the markers before; the same
markers, the stage at which the relaxed targets `(x̄, v̄)` and their KKT weights are computed;
the end point of each symplectic step (one per substep) at the old weights, `v = p / (m w)`;
and the markers of `result`, after the gauge step that assigns the KKT weights `w̄` at fixed
`(x, p)`. The end points are recomputed from the parameters in `result.info`. Positions after
the first state are wrapped into `[0, Lx)`.
"""
function remap_states(method::SplineRedistribution, before::Markers, result::RemapResult)
    info = result.info
    K = length(info.θ)
    states = [RemapStage(before, "before"), RemapStage(before, "relaxed targets")]
    z = _packed_state(before)
    N = length(before)
    for (k, θ) in enumerate(info.θ)
        z = forward_redistribution_step(info.hamiltonian, method.scheme, z[1:N],
            z[(N + 1):(2N)], θ; τ = method.τ).z
        push!(states,
            RemapStage(_step_markers(z, before.w, before.m, method.Lx), _step_label(k, K)))
    end
    push!(states, RemapStage(result.markers, "gauge step, KKT weights"))
end

"""
    remap_frames(method::SplineRedistribution, before::Markers, result::RemapResult;
        s = range(0, 1; length = 41), hold = 8, frames_per_stage = 30) -> Vector{RemapStage}

The frames of an animation of the redistribution of `before` that gave `result`. The markers
hold before and at the stage of the relaxed targets; each symplectic step then follows its
pseudo-time path `Ψ_{s τ}^θ` ([`pseudotime_path`](@ref)) for the values `s` (from 0 to 1) at the
old weights; the gauge step to the KKT weights is interpolated at fixed `(x, p)`
([`interpolate_stages`](@ref ParticleRemapping.interpolate_stages), `frames_per_stage`). Every
stage holds for `hold` frames.
"""
function remap_frames(method::SplineRedistribution, before::Markers, result::RemapResult;
        s = range(0, 1; length = 41), hold::Integer = 8, frames_per_stage::Integer = 30)
    (first(s) == 0 && last(s) == 1) ||
        throw(ArgumentError("the pseudo-time values s must run from 0 to 1, got $(extrema(s))"))
    info = result.info
    states = remap_states(method, before, result)
    K = length(info.θ)
    N = length(before)
    frames = [fill(states[1], hold); fill(states[2], hold)]
    z = _packed_state(before)
    for (k, θ) in enumerate(info.θ)
        path = pseudotime_path(info.hamiltonian, method.scheme, z[1:N], z[(N + 1):(2N)], θ;
            s = s, τ = method.τ)
        for (si, zs) in zip(s[2:end], path[2:end])
            label = _step_label(k, K) * ", s = " * string(round(si; digits = 2))
            push!(frames,
                RemapStage(_step_markers(zs, before.w, before.m, method.Lx), label))
        end
        append!(frames, fill(states[2 + k], hold))
        z = path[end]
    end
    gauge = interpolate_stages(states[(end - 1):end]; frames_per_stage = frames_per_stage,
        hold = hold)
    append!(frames, gauge[(hold + 1):end])
end
