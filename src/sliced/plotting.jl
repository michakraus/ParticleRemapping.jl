# Plotting functions of the sliced redistribution. The methods are in the Makie extension, which
# loads with any Makie backend (`using CairoMakie`).

"""
    plot_sliced_stages(before::Markers, result::RemapResult; Lx = π, size = (1500, 820))

The stages of the sliced map `R = F ∘ T_μ ∘ Ψ_p ∘ Ψ_q` on one example, in six panels: the input
markers `(x, v)` with their velocity slices (bands) and the position targets; `Ψ_q` in `(q̃, p̃)`,
coloured by velocity slice, with the targets on the axis; the position slices and momentum targets
of `Ψ_p` and its result `(q̂, p̂)`; the squeezed state `(q̄, p̄)` coloured by `μ`; the multipliers
`μ` against `p̂`; and the output `(x̄, v̄)` coloured by weight. Momenta are shown in velocity units
`p / (m w_ref)`; the output positions are wrapped into `[−Lx, Lx)`. `result` is the
[`remap`](@ref) of `before` with a [`SlicedQuantileRemap`](@ref). Needs a Makie backend.
"""
function plot_sliced_stages end

"""
    plot_checkerboard_filling(ex; size = (1300, 430))

The checkerboard of [`checkerboard_experiment`](@ref): the markers before, after the sliced
lifts and after the one-slice (marginal) control, with the quadrant counts in the titles and the
quadrant boundaries. Needs a Makie backend.
"""
function plot_checkerboard_filling end

"""
    plot_cell_count_history(histories; labels, remap_times = Float64[], size = (900, 600))

The number of empty cells and the relative standard deviation of the cell counts of
[`phase_cell_statistics`](@ref) over time, one line per history of
[`cell_count_history`](@ref), with vertical lines at the remap times. Needs a Makie backend.
"""
function plot_cell_count_history end

"""
    plot_wave_remap(ex; Lx = π, Vmax = 2.25, size = (1300, 880))

The long-time run of [`long_time_remap_experiment`](@ref): the markers at `T` before and after the
remap on the `10 × 10` cell grid of the statistics on `[−Lx, Lx) × [−Vmax, Vmax]` (colour:
weight), and the cell counts before and after as heat maps. Needs a Makie backend.
"""
function plot_wave_remap end

"""
    animate_sliced_stages(path, before::Markers, result::RemapResult; framerate = 30,
                          frames_per_stage = 40)

An animation of the stages of the sliced map: the markers move from the input through `Ψ_q`,
`Ψ_p` and `T_μ` (in `(q, p / (m w_ref))`) to the output `(x̄, v̄)`, by linear interpolation between
consecutive stages. Written to `path` (MP4). Needs a Makie backend.
"""
function animate_sliced_stages end

"""
    animate_long_time_run(path, run, history; framerate = 10, Lx = π, vlim = 3.0)

An animation of a run of [`electric_wave_simulation`](@ref) with recorded snapshots: the markers
in `[−Lx, Lx) × [−vlim, vlim]` (colour: weight) and, below, the number of empty cells and the relative standard deviation of the
cell counts up to the current time (`history` from [`cell_count_history`](@ref)); the remap
times are marked. Written to `path` (MP4). Needs a Makie backend.
"""
function animate_long_time_run end
