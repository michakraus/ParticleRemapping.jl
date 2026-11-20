# Plotting functions of the sliced redistribution. The methods are in the Makie extension, which
# loads with any Makie backend (`using CairoMakie`).

"""
    plot_sliced_stages(before::Markers, result::RemapResult; Lx = π, size = (1150, 640))

The stages of the sliced map `R = F ∘ T_μ ∘ Ψ_p ∘ Ψ_q` on one example, in six panels: (a) the
input markers `(x, v)` with the boundaries of their velocity slices; (b) `Ψ_q` in `(q̃, p̃)` with
the position targets; (c) the position slices and the momentum targets of `Ψ_p` and its result
`(q̂, p̂)`; (d) the squeezed state `(q̄, p̄)`; (e) the multipliers `μ` against `p̂`; (f) the output
`(x̄, v̄)`. The positions of (d) and (f) are wrapped into `[−Lx, Lx)`. Momenta are shown in
velocity units `p / (m w_ref)`. Panels (a) and (b) colour each marker by its velocity slice and
panel (c) by its position slice, because the slices are what these panels show; (d) and (f)
colour and size the markers by the new weight `w̄ = w_ref μ`, on one scale. The weight-coded
stages at the old weights are [`plot_stages`](@ref ParticleRemapping.plot_stages) of
[`remap_states`](@ref ParticleRemapping.remap_states). `result` is the [`remap`](@ref) of `before`
with a [`SlicedQuantileRemap`](@ref). Needs a Makie backend.
"""
function plot_sliced_stages end

"""
    plot_checkerboard_filling(ex; size = (1150, 430))

The checkerboard of [`checkerboard_experiment`](@ref): (a) the markers before, (b) after the
sliced lifts and (c) after the one-slice (marginal) control, with the quadrant counts in the
titles and the quadrant boundaries dashed. The markers have equal weights, so the colour shows the
quadrant each marker comes from (`Q1` or `Q3`). Needs a Makie backend.
"""
function plot_checkerboard_filling end

"""
    plot_cell_count_history(histories; labels, remap_times = Float64[], size = (900, 600))

(a) The number of empty cells and (b) the relative standard deviation of the cell counts of
[`phase_cell_statistics`](@ref) over time, one line per history of [`cell_count_history`](@ref),
with dashed vertical lines at the remap times. Needs a Makie backend.
"""
function plot_cell_count_history end

"""
    plot_wave_remap(ex; Lx = π, Vmax = 2.25, size = (1150, 820))

The long-time run of [`long_time_remap_experiment`](@ref): (a), (b) the markers at `T` before and
after the remap, coloured and sized by weight on one scale, on the `10 × 10` cell grid of the
statistics on `[−Lx, Lx) × [−Vmax, Vmax]`; (c), (d) the cell counts before and after as heat
maps. Needs a Makie backend.
"""
function plot_wave_remap end

"""
    animate_sliced_stages(path, before::Markers, result::RemapResult; framerate = 20,
                          frames_per_stage = 40, hold = 10, Lx = π)

An animation of the stages of the sliced map in the plane of [`plot_sliced_stages`](@ref): the
markers move from the input through `Ψ_q`, `Ψ_p` and `T_μ` (in `(q, p / (m w_ref))`) to the
output `(x̄, v̄)`, by linear interpolation between consecutive stages, and each stage is held for
`hold` frames. The positions are wrapped into `[−Lx, Lx)` and move the shorter way round. The
colour shows the velocity slice of each marker. The weight-coded animation at the old weights is
[`animate_stages`](@ref ParticleRemapping.animate_stages) of
[`remap_frames`](@ref ParticleRemapping.remap_frames). Written to `path` (MP4). Needs a Makie
backend.
"""
function animate_sliced_stages end

"""
    animate_long_time_run(path, run, history; framerate = 10, Lx = π, vlim = 3.0)

An animation of a run of [`electric_wave_simulation`](@ref) with recorded snapshots: the markers
in `[−Lx, Lx) × [−vlim, vlim]`, coloured and sized by weight on one scale for all frames, and,
below, the number of empty cells and ten times the relative standard deviation of the cell
counts up to the current time (`history` from [`cell_count_history`](@ref)); the remap times are
marked. Written to `path` (MP4). Needs a Makie backend.
"""
function animate_long_time_run end
