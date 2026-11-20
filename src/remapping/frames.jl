# The states and the animation frames of the two remaps of the chapter, for `plot_stages` and
# `animate_stages`: computed here, drawn by the Makie extension.

# the markers after the gauge step with the log-stretches μ: x and p kept, W = w e^{-μ}
function _gauge_markers(before::Markers, μ::AbstractVector)
    from_canonical(before.x, momentum(before), before.w .* exp.(-μ); m = before.m)
end

function _check_gauge_result(before::Markers, result::RemapResult)
    length(before) == length(result.markers) && before.x == result.markers.x ||
        throw(ArgumentError("the result does not belong to these markers: the gauge " *
                            "rezoning keeps every position"))
end

_τlabel(τ) = "homotopy, τ = $(round(τ; digits = 2))"
_gauge_label(τ) = "gauge step: new weights, τ = $(round(τ; digits = 3))"

"""
    remap_states(method::GaugeRezoning, before::Markers, result::RemapResult; steps = 2)
        -> Vector{RemapStage}

The states of the gauge rezoning of `before` that gave `result`. The canonical map is the
identity, so the states are the markers before the remap, the gauge steps `W = w e^{−μ}` of
`steps` accepted points of the homotopy (those closest to equispaced fractions of the parameter
`τ` reached), and the result. Every state has the positions and the momenta of `before`.
"""
function remap_states(method::GaugeRezoning, before::Markers,
        result::RemapResult; steps::Integer = 2)
    _check_gauge_result(before, result)
    τend = result.info.τ
    inner = [h for h in result.info.history if h.accepted && 0 < h.τ < τend]
    picks = unique([inner[argmin([abs(h.τ - j * τend / (steps + 1)) for h in inner])]
                    for j in 1:steps if !isempty(inner)])
    vcat(RemapStage(before, "before"),
        [RemapStage(_gauge_markers(before, h.μ), _τlabel(h.τ)) for h in picks],
        RemapStage(result.markers, _gauge_label(τend)))
end

"""
    remap_frames(method::GaugeRezoning, before::Markers, result::RemapResult;
        frames_per_step = 6, hold = 8) -> Vector{RemapStage}

The frames of an animation of the gauge rezoning that follow its homotopy: `hold` copies of the
markers before the remap, then `frames_per_step` frames between each pair of consecutive accepted
points of the homotopy in `result.info.history`, with the log-stretches `μ` interpolated linearly
(each frame is a gauge step of `before`), and `hold` copies of the result.
"""
function remap_frames(method::GaugeRezoning, before::Markers,
        result::RemapResult; frames_per_step::Integer = 6, hold::Integer = 8)
    _check_gauge_result(before, result)
    frames_per_step ≥ 1 ||
        throw(ArgumentError("frames_per_step must be positive, got $frames_per_step"))
    acc = [h for h in result.info.history if h.accepted]
    frames = fill(RemapStage(before, "before"), hold)
    for k in 2:length(acc)
        a, b = acc[k - 1], acc[k]
        for s in range(0, 1; length = frames_per_step + 1)[2:end]
            μ = (1 - s) .* a.μ .+ s .* b.μ
            push!(frames,
                RemapStage(_gauge_markers(before, μ), _τlabel((1 - s) * a.τ + s * b.τ)))
        end
    end
    append!(frames, fill(RemapStage(result.markers, _gauge_label(result.info.τ)), hold))
end

"""
    remap_states(method::CurrentQuantileMatching, before::Markers, result::RemapResult)
        -> Vector{RemapStage}

The states of the canonical transport of `before` that gave `result`: the markers before the
remap, after the canonical map `ψ` at the old weights (`V = P/(m w)`, positions wrapped into the
periodic domain of the grid), and after the current-quantile matching, which is the gauge step
that sets the new weights.
"""
function remap_states(method::CurrentQuantileMatching, before::Markers,
        result::RemapResult)
    length(before) == length(result.markers) ||
        throw(ArgumentError("the result holds $(length(result.markers)) markers, the " *
                            "input $(length(before))"))
    moved = from_canonical(
        mod.(result.info.X, method.grid.L), result.info.P, before.w; m = before.m)
    [RemapStage(before, "before"), RemapStage(moved, "ψ at the old weights"),
        RemapStage(result.markers, "gauge step: matched weights")]
end

"""
    remap_frames(method::CurrentQuantileMatching, before::Markers, result::RemapResult;
        period = (0, L), kwargs...) -> Vector{RemapStage}

The frames of an animation of the canonical transport: its states (`remap_states`) interpolated
by `interpolate_stages`, with the positions periodic on the domain `[0, L)` of the grid; `kwargs`
go to `interpolate_stages`.
"""
function remap_frames(method::CurrentQuantileMatching, before::Markers,
        result::RemapResult; period = (0, method.grid.L), kwargs...)
    interpolate_stages(remap_states(method, before, result); period, kwargs...)
end
