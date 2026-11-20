# Plotting functions of the common layer. The methods are in the Makie extension, which loads with
# any Makie backend (`using CairoMakie`).

"""
    plot_markers(markers; color = :weight, title = "")
    plot_markers!(ax, markers; color = :weight, markersize = :weight, wref, colorrange)

Scatter plot of 1D-1V markers in `(x, v)`, in the common style of the package: the colour shows
the weight (`color = :weight`, colour map `viridis`, on a logarithmic scale when the weights span
more than two decades) and the marker area is proportional to the weight
(`markersize = :weight`), relative to the reference weight `wref` (default: the largest weight).
A number for `color` or `markersize` fixes it. Needs a Makie backend.
"""
function plot_markers end
@doc (@doc plot_markers) function plot_markers! end

"""
    plot_remap(before, after; titles = ("before", "after"), targets = nothing, period = nothing)

Two panels of [`plot_markers`](@ref) side by side, with a shared colour scale and a shared size
scale of the weights. `targets = (x, v)` draws target points as open circles in the second panel;
`period = (a, b)` wraps the positions into `[a, b)`. Needs a Makie backend.
"""
function plot_remap end

"""
    plot_stages(stages; columns = 4, targets = nothing, window = nothing, period = nothing,
        decorate! = nothing)

One panel per [`RemapStage`](@ref), titled by its label, with the markers coloured and sized by
weight on one scale for all panels. `targets = (x, v)` draws target points as open circles in the
panels of the stages whose index is in `target_panels` (default: all but the first);
`window = (x₀, x₁, v₀, v₁)` restricts the view; `decorate!(ax, k)` draws more into panel `k`, for
example the outline of a hole. Needs a Makie backend.
"""
function plot_stages end

"""
    animate_stages(path, frames; framerate = 20, targets = nothing, window = nothing,
        period = nothing, decorate! = nothing)

Writes an animation of the [`RemapStage`](@ref) frames to `path` (an `.mp4` file), with the
markers coloured and sized by weight on one scale for all frames and the label of each frame as
the title. [`interpolate_stages`](@ref) makes frames from the stages of a remap. The keywords are
those of [`plot_stages`](@ref). Needs a Makie backend.
"""
function animate_stages end

"""
    remap_theme()

The Makie theme of every figure of the package: font sizes, the colour cycle (Wong's palette),
the default colour map and the marker style. `set_theme!(remap_theme())` applies it; the
documentation and the figure scripts do so. Needs a Makie backend.
"""
function remap_theme end

export plot_markers, plot_markers!, plot_remap, plot_stages, animate_stages, remap_theme
