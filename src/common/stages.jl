"""
    RemapStage(markers, label)

One state of a remap: the markers after a stage of the map, with a label for a figure panel or an
animation frame. Every construction returns the states of its remap as a vector of stages, from
the markers before the remap to the markers after it, for [`plot_stages`](@ref) and
[`animate_stages`](@ref).
"""
struct RemapStage{M <: Markers}
    markers::M
    label::String
end

RemapStage(markers::Markers, label::AbstractString) = RemapStage(markers, String(label))

"""
    interpolate_stages(stages; frames_per_stage = 30, hold = 8, period = nothing)

The frames of an animation through `stages`: `frames_per_stage` frames between each pair of
consecutive stages and `hold` copies of every stage. Between two stages the positions, the
momenta `p = m w v` and the weights change linearly, and the velocities are `v = p / (m w)`. A
canonical stage, at fixed weights, thus moves the markers on straight lines, and a gauge step, at
fixed `(x, p)`, moves them along `v ∝ 1/w`. With `period = (a, b)` the positions are periodic on
`[a, b)` and move by the shorter way round. The label of each frame is that of the stage it moves
towards.
"""
function interpolate_stages(
        stages::AbstractVector{<:RemapStage}; frames_per_stage::Integer = 30,
        hold::Integer = 8, period = nothing)
    isempty(stages) && throw(ArgumentError("no stages to interpolate"))
    frames = eltype(stages)[]
    append!(frames, fill(stages[1], hold))
    for k in 2:length(stages)
        a, b = stages[k - 1].markers, stages[k].markers
        length(a) == length(b) || throw(DimensionMismatch(
            "stage $(k - 1) has $(length(a)) markers, stage $k has $(length(b))"))
        pa, pb = momentum(a), momentum(b)
        dx = b.x .- a.x
        if period !== nothing
            L = period[2] - period[1]
            dx = mod.(dx .+ L / 2, L) .- L / 2
        end
        for s in range(0, 1; length = frames_per_stage + 1)[2:end]
            w = (1 - s) .* a.w .+ s .* b.w
            x = a.x .+ s .* dx
            period === nothing ||
                (x = mod.(x .- period[1], period[2] - period[1]) .+ period[1])
            p = (1 - s) .* pa .+ s .* pb
            push!(frames, RemapStage(from_canonical(x, p, w; m = a.m), stages[k].label))
        end
        append!(frames, fill(stages[k], hold))
    end
    frames
end

"""
    remap_states(method::AbstractRemap, before::Markers, result::RemapResult) -> Vector{RemapStage}

The states of the remap of `before` by `method` that gave `result`: the markers before the remap,
after each stage of the canonical map `ψ` (at the old weights, `v = p / (m w)`), and after the
gauge step that sets the new weights, which are the markers of `result`. Each construction adds
a method for its remap types. [`plot_stages`](@ref) draws the states.
"""
function remap_states end

"""
    remap_frames(method::AbstractRemap, before::Markers, result::RemapResult; kwargs...)
        -> Vector{RemapStage}

The frames of an animation of the remap of `before` by `method` that gave `result`, for
[`animate_stages`](@ref). The generic method interpolates the [`remap_states`](@ref) with
[`interpolate_stages`](@ref) and passes `kwargs` to it; a construction whose map has a natural
path, such as a homotopy or the flow of a Hamiltonian, adds a method that follows that path.
"""
function remap_frames(method::AbstractRemap, before::Markers, result::RemapResult; kwargs...)
    interpolate_stages(remap_states(method, before, result); kwargs...)
end

export RemapStage, interpolate_stages, remap_states, remap_frames
