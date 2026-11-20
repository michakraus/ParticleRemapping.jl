# The states and the frames of a sliced remap, for the standard figures of the chapter: the
# markers before, after each stage of the canonical map ψ = T_μ ∘ Ψ_p ∘ Ψ_q at the old weights,
# and after the gauge step F.

const _STAGE_LABELS = ("before", "Ψ_q: lift in velocity slices",
    "Ψ_p: lift in position slices", "T_μ: weight squeeze", "gauge step: new weights")

# A position mapped into [−Lx, Lx); a position inside is returned as it is.
_into_box(x::Real, Lx::Real) = -Lx <= x < Lx ? x : periodic_wrap(x, Lx)

# The markers at the canonical state (q, p), with the positions in the box, at the weights w.
function _state_markers(q, p, w, m, Lx)
    Markers(_into_box.(q, Lx), p ./ (m .* _per_marker(w, p)), w; m = m)
end

"""
    remap_states(method::SlicedQuantileRemap, before::Markers, result::RemapResult)
        -> Vector{RemapStage}

The five states of the sliced remap of `before` that gave `result`: the markers before; after
`Ψ_q`, after `Ψ_p` and after the squeeze `T_μ`, each at the old weights `w` with `v = p / (m w)`;
and after the gauge step `F` (the markers of `result`, with `w̄ = w_ref μ`). The positions are
mapped into `[−Lx, Lx)`. For `d × N` markers the states are `d × N` as well.
"""
function ParticleRemapping.remap_states(
        method::SlicedQuantileRemap, before::Markers, result::RemapResult)
    i, m, w, Lx = result.info, before.m, before.w, method.Lx
    after = result.markers
    length(before) == length(after) || throw(DimensionMismatch(
        "before has $(length(before)) markers, the result $(length(after))"))
    markers = (
        before, _state_markers(i.q̃, i.p̃, w, m, Lx), _state_markers(i.q̂, i.p̂, w, m, Lx),
        _state_markers(i.q̄, i.p̄, w, m, Lx),
        Markers(_into_box.(after.x, Lx), after.v, after.w; m = m))
    [RemapStage(mk, label) for (mk, label) in zip(markers, _STAGE_LABELS)]
end

# The rows of a vector (one) or of a d × N matrix (d), and back.
_axis_rows(x::AbstractVector) = [x]
_axis_rows(x::AbstractMatrix) = [x[ℓ, :] for ℓ in axes(x, 1)]
_from_rows(rows, ::AbstractVector) = only(rows)
_from_rows(rows, ::AbstractMatrix) = permutedims(reduce(hcat, rows))

_plan_list(plan::SlicedLiftPlan) = [plan]
_plan_list(plans::AbstractVector{<:SlicedLiftPlan}) = plans

# The lifts Ψ_q and Ψ_p of the frozen plans with the relaxation parameters β scaled by s: the
# identity for s = 0, the lifts of the remap for s = 1.
function _scaled_lift_q(plans, q, p, s)
    rows = [split_quantile_ot_lift(qℓ, pℓ, pl.bins_v, pl.taus_x; anchors = pl.p_anc,
                β_mean = s * pl.β_mean_q, β_intra = s * pl.β_intra_q, ε = pl.ε_q)
            for (pl, qℓ, pℓ) in zip(plans, _axis_rows(q), _axis_rows(p))]
    _from_rows(first.(rows), q), _from_rows(last.(rows), p)
end

function _scaled_lift_p(plans, q, p, s)
    rows = [dual_quantile_ot_lift(qℓ, pℓ, pl.bins_x, pl.taus_p; anchors = pl.q_anc,
                β_mean = s * pl.β_mean_p, β_intra = s * pl.β_intra_p, ε = pl.ε_p)
            for (pl, qℓ, pℓ) in zip(plans, _axis_rows(q), _axis_rows(p))]
    _from_rows(first.(rows), q), _from_rows(last.(rows), p)
end

"""
    remap_frames(method::SlicedQuantileRemap, before::Markers, result::RemapResult;
                 frames_per_stage = 30, hold = 8) -> Vector{RemapStage}

The frames of an animation of the sliced remap of `before` that gave `result`, along a path of
canonical maps: `frames_per_stage` frames of `Ψ_q` and of `Ψ_p` with the relaxation parameters
`β_mean`, `β_intra` of the frozen plan scaled by `s ∈ (0, 1]` (the identity at `s = 0`); of the
squeeze `T_μ` with the multipliers `(1 − s) + s μ`, which keep `Σ μ = N`; all at the old weights,
`v = p / (m w)`. Then the gauge step at the fixed state `(x̄, p̄)` with the weights
`(1 − s) w + s w̄`. Each state of [`remap_states`](@ref) is held for `hold` frames, and the label
of a frame is that of the state it moves towards. The positions are mapped into `[−Lx, Lx)`.
"""
function ParticleRemapping.remap_frames(
        method::SlicedQuantileRemap, before::Markers, result::RemapResult;
        frames_per_stage::Integer = 30, hold::Integer = 8)
    i, m, w, Lx = result.info, before.m, before.w, method.Lx
    plans = _plan_list(i.plan)
    q_anc = _from_rows([pl.q_anc for pl in plans], i.q̂)
    states = remap_states(method, before, result)
    path = (s -> _scaled_lift_q(plans, before.x, momentum(before), s),
        s -> _scaled_lift_p(plans, i.q̃, i.p̃, s),
        s -> _squeeze_axes(i.q̂, i.p̂, (1 - s) .+ s .* i.μ, q_anc))
    frames = fill(states[1], hold)
    steps = range(0, 1; length = frames_per_stage + 1)[2:end]
    for (k, stage) in enumerate(path)
        for s in steps
            q, p = stage(s)
            push!(frames, RemapStage(_state_markers(q, p, w, m, Lx), states[k + 1].label))
        end
        append!(frames, fill(states[k + 1], hold))
    end
    for s in steps
        push!(frames,
            RemapStage(
                _state_markers(i.q̄, i.p̄, (1 - s) .* w .+ s .* result.markers.w, m, Lx),
                states[5].label))
    end
    append!(frames, fill(states[5], hold))
end
