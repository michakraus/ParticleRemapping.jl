# The common style of every figure and animation, and the plots of the common layer.
#
# Style rules, used by the plots of every construction:
#
# - markers in phase space: colour and area by weight, `_WEIGHT_CMAP`, one colour range and one
#   size reference for all panels or frames of a figure (`_weight_scale`), robust to a few extreme
#   weights, and a linear colour scale unless a figure asks for `colorscale = log10`;
# - target points: open circles in `_TARGET_COLOR`; displacements: thin grey segments;
#   outlines of holes and regions: `_OUTLINE_COLOR`, dashed;
# - curves and series: Wong's palette, in order (`_SERIES`);
# - fields: `_FIELD_CMAP` for a non-negative field (a density, counts, a magnitude),
#   `_SIGNED_CMAP` with a symmetric colour range for a signed field;
# - categories (slices): `_category_colors`;
# - panels of a multi-panel figure are titled "(a) …", "(b) …"; axis labels are LaTeX strings.

const _SERIES = Makie.wong_colors()
const _WEIGHT_CMAP = :viridis
const _FIELD_CMAP = :cividis
const _SIGNED_CMAP = :balance
const _TARGET_COLOR = _SERIES[6]
const _OUTLINE_COLOR = :gray25
const _SEGMENT_COLOR = (:gray45, 0.6)

_category_colors(n) = [_SERIES[mod1(k, 7)] for k in 1:n]

# The largest marker size for an ensemble of `N` markers: large for a few markers, small for many.
_max_markersize(N) = clamp(110 / sqrt(N), 2.5, 16.0)

# Marker sizes with the area proportional to the weight, relative to `wref`, at most the largest
# size; never below 35 % of it, so that light markers stay visible.
function _weight_markersize(w, wref, N)
    smax = _max_markersize(N)
    max.(max(0.35smax, 1.5), smax .* sqrt.(clamp.(w ./ wref, 0, 1)))
end

# The colour range, the colour scale and the size reference of the weights of several marker sets.
# For 200 weights or more the range runs from the 0.5 % to the 99.5 % quantile, so that a few
# extreme weights do not wash out the colours of the others; weights outside it take the colour and
# the size of its ends. A logarithmic `colorscale` needs positive weights.
function _weight_scale(ws; colorscale = identity)
    all = sort(reduce(vcat, ws))
    n = length(all)
    lo, hi = n < 200 ? (all[1], all[end]) :
             (all[max(1, ceil(Int, 0.005n))], all[ceil(Int, 0.995n)])
    if lo == hi
        lo, hi = lo - abs(lo) / 10 - eps(), hi + abs(hi) / 10 + eps()
    end
    colorscale === identity || lo > 0 ||
        throw(ArgumentError("a logarithmic colour scale needs positive weights"))
    (colorrange = (lo, hi), colorscale, wref = maximum(abs, (lo, hi)))
end

_wrap(x, ::Nothing) = x
_wrap(x, period) = mod.(x .- period[1], period[2] - period[1]) .+ period[1]

function ParticleRemapping.remap_theme()
    Theme(;
        fontsize = 14,
        palette = (color = _SERIES,),
        colormap = _WEIGHT_CMAP,
        Axis = (titlesize = 14, titlefont = :regular, xgridvisible = false,
            ygridvisible = false),
        Scatter = (strokewidth = 0,),
        Legend = (framevisible = false,),
        Colorbar = (width = 10,))
end

function ParticleRemapping.plot_markers!(
        ax, mk::Markers; color = :weight, markersize = :weight, wref = maximum(abs, mk.w),
        colorrange = nothing, colorscale = nothing, period = nothing, kwargs...)
    mk.x isa AbstractVector ||
        throw(ArgumentError("plot_markers! draws 1D-1V markers; x is a $(size(mk.x)) matrix"))
    sc = _weight_scale([mk.w])
    cr = colorrange === nothing ? sc.colorrange : colorrange
    cs = colorscale === nothing ? sc.colorscale : colorscale
    ms = markersize === :weight ? _weight_markersize(mk.w, wref, length(mk)) : markersize
    if color === :weight
        scatter!(ax, _wrap(mk.x, period), mk.v; color = mk.w, colorrange = cr,
            colorscale = cs, colormap = _WEIGHT_CMAP, markersize = ms, strokewidth = 0,
            kwargs...)
    else
        scatter!(ax, _wrap(mk.x, period), mk.v; color, markersize = ms, strokewidth = 0,
            kwargs...)
    end
end

# Target points as open circles; thinner and smaller for many markers, so that they do not hide
# the markers.
function _targets!(ax, targets, N; period = nothing, kwargs...)
    dense = N > 100
    scatter!(ax, _wrap(targets[1], period), targets[2]; color = :transparent,
        strokecolor = dense ? (_TARGET_COLOR, 0.7) : _TARGET_COLOR,
        strokewidth = dense ? 0.6 : 1.0,
        markersize = (dense ? 1.25 : 1.6) * _max_markersize(N), kwargs...)
end

# The x-limits of a periodic domain, with a margin for the markers and circles at its ends.
function _period_xlims!(ax, period)
    pad = (period[2] - period[1]) / 40
    xlims!(ax, period[1] - pad, period[2] + pad)
end

# A title with the subscripts of `a_b` or `a_{bc}` set as subscripts, such as "Ψ_q" or "w_ref".
function _richtitle(s::AbstractString)
    parts = Any[]
    i = firstindex(s)
    for m in eachmatch(r"_(\{[^}]*\}|[^\s:,;()=]+)", s)
        m.offset > i && push!(parts, s[i:prevind(s, m.offset)])
        push!(parts, subscript(String(strip(m.captures[1], ['{', '}']))))
        i = m.offset + ncodeunits(m.match)
    end
    i <= lastindex(s) && push!(parts, s[i:end])
    rich(parts...)
end

function ParticleRemapping.plot_markers(mk::Markers; title = "", size = (520, 400),
        kwargs...)
    fig = Figure(; size = size)
    ax = Axis(fig[1, 1]; xlabel = L"x", ylabel = L"v", title = title)
    sc = ParticleRemapping.plot_markers!(ax, mk; kwargs...)
    Colorbar(fig[1, 2], sc; label = "weight w")
    fig
end

function ParticleRemapping.plot_remap(before::Markers, after::Markers;
        titles = ("before", "after"), targets = nothing, period = nothing,
        colorscale = identity, size = (980, 400), kwargs...)
    fig = Figure(; size = size)
    s = _weight_scale([before.w, after.w]; colorscale)
    axes = [Axis(fig[1, k]; xlabel = L"x", ylabel = k == 1 ? L"v" : "", title = titles[k])
            for k in 1:2]
    linkaxes!(axes...)
    hideydecorations!(axes[2]; ticks = false, ticklabels = false)
    ParticleRemapping.plot_markers!(axes[1], before; s.colorrange, s.colorscale, s.wref,
        period, kwargs...)
    targets === nothing || _targets!(axes[2], targets, length(after); period)
    sc = ParticleRemapping.plot_markers!(
        axes[2], after; s.colorrange, s.colorscale, s.wref,
        period, kwargs...)
    period === nothing || _period_xlims!(axes[1], period)
    Colorbar(fig[1, 3], sc; label = "weight w")
    fig
end

function ParticleRemapping.plot_stages(stages::AbstractVector{<:RemapStage};
        columns = min(4, length(stages)), targets = nothing,
        target_panels = 2:length(stages), window = nothing, period = nothing,
        decorate! = nothing, ylabel = L"v", panel_size = (240, 230), colorscale = identity,
        kwargs...)
    n = length(stages)
    rows = cld(n, columns)
    s = _weight_scale([st.markers.w for st in stages]; colorscale)
    fig = Figure(; size = (panel_size[1] * columns + 80, panel_size[2] * rows + 20))
    local sc
    axes = Axis[]
    for (k, st) in enumerate(stages)
        r, c = fldmod1(k, columns)
        # x labels on every panel with no panel below it
        bottom = k + columns > n
        ax = Axis(fig[r, c]; title = _richtitle("($(Char('a' + k - 1))) " * st.label),
            xlabel = bottom ? L"x" : "", ylabel = c == 1 ? ylabel : "",
            yticklabelsvisible = c == 1, xticklabelsvisible = bottom)
        push!(axes, ax)
        decorate! === nothing || decorate!(ax, k)
        (targets !== nothing && k in target_panels) &&
            _targets!(ax, targets, length(st.markers); period)
        sc = ParticleRemapping.plot_markers!(ax, st.markers; s.colorrange, s.colorscale,
            s.wref, period, kwargs...)
        window === nothing || limits!(ax, window...)
    end
    linkaxes!(axes...)
    (window === nothing && period !== nothing) && _period_xlims!(axes[1], period)
    Colorbar(fig[1:rows, columns + 1], sc; label = "weight w")
    fig
end

function ParticleRemapping.animate_stages(path, frames::AbstractVector{<:RemapStage};
        framerate = 20, targets = nothing,
        target_frames = something(findfirst(f -> f.label != frames[1].label, frames), 1):length(frames),
        window = nothing, period = nothing, decorate! = nothing, ylabel = L"v",
        colorscale = identity, size = (720, 480))
    s = _weight_scale([f.markers.w for f in frames]; colorscale)
    N = length(frames[1].markers)
    k = Observable(1)
    mk = lift(i -> frames[i].markers, k)
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; xlabel = L"x", ylabel,
        title = lift(i -> _richtitle(frames[i].label), k))
    decorate! === nothing || decorate!(ax, 1)
    targets === nothing ||
        _targets!(ax, targets, N; period, visible = lift(i -> i in target_frames, k))
    sc = scatter!(ax, lift(m -> _wrap(m.x, period), mk), lift(m -> m.v, mk);
        color = lift(m -> m.w, mk), colorrange = s.colorrange, colorscale = s.colorscale,
        colormap = _WEIGHT_CMAP, strokewidth = 0,
        markersize = lift(m -> _weight_markersize(m.w, s.wref, N), mk))
    if window === nothing
        xs = reduce(vcat, [_wrap(f.markers.x, period) for f in frames])
        vs = reduce(vcat, [f.markers.v for f in frames])
        period === nothing ? xlims!(ax, extrema(xs)...) : _period_xlims!(ax, period)
        ylims!(ax, extrema(vs)...)
    else
        limits!(ax, window...)
    end
    Colorbar(fig[1, 2], sc; label = "weight w")
    record(fig, path, eachindex(frames); framerate) do i
        k[] = i
    end
    path
end
