# Plots of the common layer.

_markersizes(w, ms::Real) = ms
_markersizes(w, ::Symbol) = 3 .+ 9 .* sqrt.(w ./ maximum(w))

function ParticleRemapping.plot_markers!(
        ax, mk::Markers; color = :weight, markersize = :weight,
        colorrange = extrema(mk.w), colormap = :viridis, kwargs...)
    mk.x isa AbstractVector ||
        throw(ArgumentError("plot_markers! draws 1D-1V markers; x is a $(size(mk.x)) matrix"))
    c = color === :weight ? mk.w : color
    scatter!(ax, mk.x, mk.v; color = c, colorrange = colorrange, colormap = colormap,
        markersize = _markersizes(mk.w, markersize), strokewidth = 0, kwargs...)
end

function ParticleRemapping.plot_markers(mk::Markers; title = "", size = (520, 400), kwargs...)
    fig = Figure(; size = size)
    ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "v", title = title)
    sc = ParticleRemapping.plot_markers!(ax, mk; kwargs...)
    Colorbar(fig[1, 2], sc; label = "weight")
    fig
end

function ParticleRemapping.plot_remap(before::Markers, after::Markers;
        titles = ("before", "after"), size = (980, 400), kwargs...)
    fig = Figure(; size = size)
    cr = extrema(vcat(before.w, after.w))
    axes = [Axis(fig[1, k]; xlabel = "x", ylabel = "v", title = titles[k]) for k in 1:2]
    linkaxes!(axes...)
    ParticleRemapping.plot_markers!(axes[1], before; colorrange = cr, kwargs...)
    sc = ParticleRemapping.plot_markers!(axes[2], after; colorrange = cr, kwargs...)
    Colorbar(fig[1, 3], sc; label = "weight")
    fig
end
