# Plotting functions of the common layer. The methods are in the Makie extension, which loads with
# any Makie backend (`using CairoMakie`).

"""
    plot_markers(markers; color = :weight, title = "")
    plot_markers!(ax, markers; color = :weight, markersize = :weight)

Scatter plot of 1D-1V markers in `(x, v)`. The colour shows the weight (`color = :weight`) or is
fixed; the marker area is proportional to the weight (`markersize = :weight`) or fixed. Needs a
Makie backend.
"""
function plot_markers end
@doc (@doc plot_markers) function plot_markers! end

"""
    plot_remap(before, after; titles = ("before", "after"))

Two panels of [`plot_markers`](@ref) side by side, with a shared colour scale of the weights.
Needs a Makie backend.
"""
function plot_remap end

export plot_markers, plot_markers!, plot_remap
