# SVG pictures of the steps of a remap: one panel per state, and one animated panel that moves
# the markers from state to state (SMIL, no script). They need no plotting package.

"""
    SVGFigure(svg::String)

An SVG picture. It shows as `image/svg+xml`, so that a notebook or a Documenter `@example` block
draws it, and `write(file, fig)` writes it.
"""
struct SVGFigure
    svg::String
end

Base.show(io::IO, ::MIME"image/svg+xml", f::SVGFigure) = print(io, f.svg)
Base.write(io::IO, f::SVGFigure) = write(io, f.svg)

"""
    svg_line(x1, v1, x2, v2; color = "#404040", width = 1, dash = false)
    svg_rect(x, v, Δx, Δv; color = "#56B4E9", opacity = 0.15)
    svg_text(x, v, text; color = "#404040", size = 11)
    svg_arrow(x1, v1, x2, v2; color = "#808080", width = 1)

The overlays of [`svg_stage_panels`](@ref) and [`svg_stage_animation`](@ref), in the
coordinates `(x, v)` of the markers: a line, a filled rectangle with the corner `(x, v)`, a text
centred at `(x, v)`, and a line with an arrow head at `(x2, v2)`.
"""
function svg_line(x1, v1, x2, v2; color = "#404040", width = 1, dash = false)
    (; kind = :line,
        x1, v1, x2, v2, color, width, dash)
end
function svg_rect(x, v, Δx, Δv; color = "#56B4E9", opacity = 0.15)
    (; kind = :rect, x, v, Δx, Δv,
        color, opacity)
end
function svg_text(x, v, text; color = "#404040", size = 11)
    (; kind = :text, x, v, text, color, size)
end
function svg_arrow(x1, v1, x2, v2; color = "#808080", width = 1)
    (; kind = :arrow, x1, v1, x2, v2,
        color, width)
end
@doc (@doc svg_line) svg_rect
@doc (@doc svg_line) svg_text
@doc (@doc svg_line) svg_arrow

# viridis at eleven equidistant points
const _VIRIDIS = [(0.267, 0.005, 0.329), (0.283, 0.141, 0.458), (0.254, 0.265, 0.530),
    (0.207, 0.372, 0.553), (0.164, 0.471, 0.558), (0.128, 0.567, 0.551),
    (0.135, 0.659, 0.518), (0.267, 0.749, 0.441), (0.478, 0.821, 0.318),
    (0.741, 0.873, 0.150), (0.993, 0.906, 0.144)]

function _viridis(t)
    s = clamp(t, 0, 1) * (length(_VIRIDIS) - 1)
    i = min(floor(Int, s), length(_VIRIDIS) - 2)
    θ = s - i
    c = (1 - θ) .* _VIRIDIS[i + 1] .+ θ .* _VIRIDIS[i + 2]
    "rgb(" * join(round.(Int, 255 .* c), ",") * ")"
end

_fmt(x) = string(round(x; digits = 2))

# a text with the characters that XML reserves replaced by their entities
_escape(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")

# The frame of one panel: data → pixel maps and the drawing of its box, labels and title.
struct _Panel
    x0::Float64
    y0::Float64
    w::Float64
    h::Float64
    window::NTuple{4, Float64}
end

_px(p::_Panel, x) = p.x0 + (x - p.window[1]) / (p.window[2] - p.window[1]) * p.w
_py(p::_Panel, v) = p.y0 + (p.window[4] - v) / (p.window[4] - p.window[3]) * p.h

function _frame(io, p::_Panel)
    print(io, "<rect x=\"$(_fmt(p.x0))\" y=\"$(_fmt(p.y0))\" width=\"$(_fmt(p.w))\" ",
        "height=\"$(_fmt(p.h))\" fill=\"white\" stroke=\"#404040\"/>")
    print(io, "<text x=\"$(_fmt(p.x0 + p.w / 2))\" y=\"$(_fmt(p.y0 + p.h + 16))\" ",
        "font-size=\"12\" font-style=\"italic\" text-anchor=\"middle\">x</text>")
    print(io, "<text x=\"$(_fmt(p.x0 - 10))\" y=\"$(_fmt(p.y0 + p.h / 2 + 4))\" ",
        "font-size=\"12\" font-style=\"italic\" text-anchor=\"middle\">v</text>")
end

function _overlay(io, p::_Panel, s)
    if s.kind === :line || s.kind === :arrow
        dash = s.kind === :line && s.dash ? " stroke-dasharray=\"4,3\"" : ""
        head = s.kind === :arrow ? " marker-end=\"url(#arrowhead)\"" : ""
        print(io, "<line x1=\"$(_fmt(_px(p, s.x1)))\" y1=\"$(_fmt(_py(p, s.v1)))\" ",
            "x2=\"$(_fmt(_px(p, s.x2)))\" y2=\"$(_fmt(_py(p, s.v2)))\" stroke=\"$(s.color)\" ",
            "stroke-width=\"$(s.width)\"$dash$head/>")
    elseif s.kind === :rect
        x1, x2 = minmax(_px(p, s.x), _px(p, s.x + s.Δx))
        y1, y2 = minmax(_py(p, s.v), _py(p, s.v + s.Δv))
        print(io, "<rect x=\"$(_fmt(x1))\" y=\"$(_fmt(y1))\" width=\"$(_fmt(x2 - x1))\" ",
            "height=\"$(_fmt(y2 - y1))\" fill=\"$(s.color)\" fill-opacity=\"$(s.opacity)\"/>")
    else
        print(io, "<text x=\"$(_fmt(_px(p, s.x)))\" y=\"$(_fmt(_py(p, s.v)))\" ",
            "font-size=\"$(s.size)\" fill=\"$(s.color)\" text-anchor=\"middle\">",
            _escape(s.text), "</text>")
    end
end

const _SVG_DEFS = "<defs><marker id=\"arrowhead\" markerWidth=\"6\" markerHeight=\"6\" " *
                  "refX=\"5\" refY=\"3\" orient=\"auto\"><path d=\"M0,0 L6,3 L0,6 z\" " *
                  "fill=\"#808080\"/></marker></defs>" *
                  # a white ground, so that the black texts also show on a dark page
                  "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>"

# Ids of clip paths, unique among the SVG pictures of one page.
const _CLIP_COUNTER = Ref(0)

# Open a group clipped to the box of the panel `p`; close it with `</g>`.
function _clip_open(io, p::_Panel)
    id = "pr-clip-$(_CLIP_COUNTER[] += 1)"
    print(io, "<clipPath id=\"$id\"><rect x=\"$(_fmt(p.x0))\" y=\"$(_fmt(p.y0))\" ",
        "width=\"$(_fmt(p.w))\" height=\"$(_fmt(p.h))\"/></clipPath>",
        "<g clip-path=\"url(#$id)\">")
end

# The colour and the radius of every marker of every state: one colour range and one size
# reference for all states, the marker area proportional to the weight.
function _marker_styles(stages, colors, rmax)
    ws = [s.markers.w for s in stages]
    lo, hi = extrema(reduce(vcat, ws))
    wref = maximum(abs, (lo, hi))
    iszero(wref) && (wref = one(wref))
    hi == lo && (hi = lo + 1)
    map(enumerate(ws)) do (k, w)
        c = _colors(colors, k)
        fill = c === nothing ? [_viridis((wa - lo) / (hi - lo)) for wa in w] : c
        (fill, max.(0.35rmax, rmax .* sqrt.(clamp.(w ./ wref, 0, 1))))
    end
end

_colors(::Nothing, k) = nothing
_colors(c::AbstractVector, k) = c
_colors(f::Function, k) = f(k)

"""
    svg_category_colors(categories)

One colour of Wong's palette for each entry of `categories` (integers from 1), cycling through
the seven colours, for the keyword `colors` of [`svg_stage_panels`](@ref) and
[`svg_stage_animation`](@ref).
"""
function svg_category_colors(categories::AbstractVector{<:Integer})
    palette = ("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442")
    [palette[mod1(c, length(palette))] for c in categories]
end

_overlays(::Nothing, k) = ()
_overlays(f::Function, k) = f(k)
_overlays(v::AbstractVector, k) = v[k]

_targets(::Nothing, k) = nothing
_targets(t::Tuple, k) = t
_targets(t::AbstractVector, k) = t[k]

function _title(p::_Panel, label)
    "<text x=\"$(_fmt(p.x0 + p.w / 2))\" y=\"$(_fmt(p.y0 - 6))\" " *
    "font-size=\"12\" text-anchor=\"middle\">" * _escape(label) * "</text>"
end

"""
    svg_stage_panels(stages; window, columns = length(stages), panel = (220, 160),
                     overlays = nothing, targets = nothing, colors = nothing, rmax = 6,
                     titles = [s.label for s in stages])

One SVG panel of the markers in `(x, v)` for each state of `stages` (a vector of
[`RemapStage`](@ref)), in rows of `columns` panels of `panel` pixels, over the box
`window = (xmin, xmax, vmin, vmax)`. The markers are coloured (viridis) and sized (area) by
weight, with one colour range and one size reference for all panels; `colors`, one colour per
marker (or `colors(k)`, such a vector or `nothing` for panel `k`), colours them by a category
instead ([`svg_category_colors`](@ref)). `overlays(k)` (or `overlays[k]`) gives the
overlays of panel `k` ([`svg_line`](@ref), [`svg_rect`](@ref), [`svg_text`](@ref),
[`svg_arrow`](@ref)), drawn below the markers; `targets = (x, v)` (or a vector of such pairs or
`nothing`, one per panel) draws target points as open circles. Returns an [`SVGFigure`](@ref).
"""
function svg_stage_panels(stages::AbstractVector{<:RemapStage}; window,
        columns::Integer = length(stages), panel = (220, 160), overlays = nothing,
        targets = nothing, colors = nothing, rmax::Real = 6,
        titles = [s.label for s in stages])
    n = length(stages)
    rows = cld(n, columns)
    mx, mtop, mbottom = 24.0, 22.0, 26.0
    W = columns * (panel[1] + mx) + mx
    H = rows * (panel[2] + mtop + mbottom)
    styles = _marker_styles(stages, colors, rmax)
    io = IOBuffer()
    print(io,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 $(_fmt(W)) $(_fmt(H))\" ",
        "width=\"100%\" style=\"max-width:$(round(Int, W))px\" font-family=\"sans-serif\">",
        _SVG_DEFS)
    for k in 1:n
        r, c = fldmod1(k, columns)
        p = _Panel(
            mx + (c - 1) * (panel[1] + mx), (r - 1) * (panel[2] + mtop + mbottom) + mtop,
            panel[1], panel[2], Float64.(Tuple(window)))
        _frame(io, p)
        _clip_open(io, p)
        foreach(s -> _overlay(io, p, s), _overlays(overlays, k))
        mk = stages[k].markers
        fillc, rad = styles[k]
        for a in eachindex(mk.x)
            print(io,
                "<circle cx=\"$(_fmt(_px(p, mk.x[a])))\" cy=\"$(_fmt(_py(p, mk.v[a])))\" ",
                "r=\"$(_fmt(rad[a]))\" fill=\"$(fillc[a])\" fill-opacity=\"0.85\"/>")
        end
        t = _targets(targets, k)
        if t !== nothing
            for (xt, vt) in zip(t...)
                print(io, "<circle cx=\"$(_fmt(_px(p, xt)))\" cy=\"$(_fmt(_py(p, vt)))\" ",
                    "r=\"$(_fmt(0.6rmax))\" fill=\"none\" stroke=\"#D55E00\"/>")
            end
        end
        print(io, "</g>", _title(p, titles[k]))
    end
    print(io, "</svg>")
    SVGFigure(String(take!(io)))
end

"""
    svg_stage_animation(stages; window, panel = (420, 300), hold = 1.2, move = 1.6,
                        overlays = nothing, targets = nothing, colors = nothing, rmax = 7,
                        titles = [s.label for s in stages])

An animated SVG of the markers through the states of `stages`: each state shows for `hold`
seconds, and between two states the markers move for `move` seconds, with the positions, the
velocities, the colours and the radii interpolated linearly; the animation repeats. The overlays,
the targets and the title of a state show while its state is held and while the markers move
towards it. The other keywords are those of [`svg_stage_panels`](@ref). Returns an
[`SVGFigure`](@ref).
"""
function svg_stage_animation(
        stages::AbstractVector{<:RemapStage}; window, panel = (420, 300),
        hold::Real = 1.2, move::Real = 1.6, overlays = nothing, targets = nothing,
        colors = nothing, rmax::Real = 7, titles = [s.label for s in stages])
    n = length(stages)
    N = length(stages[1].markers)
    all(s -> length(s.markers) == N, stages) ||
        throw(DimensionMismatch("every state must hold the same markers"))
    p = _Panel(30.0, 24.0, Float64(panel[1]), Float64(panel[2]), Float64.(Tuple(window)))
    W, H = p.w + 50, p.h + 54
    # the last move returns to the first state
    dur = n * (hold + move)
    # the times at which state k starts and ends to be held
    t_on = [(k - 1) * (hold + move) for k in 1:n]
    t_off = t_on .+ hold
    kt = Float64[]
    for k in 1:n
        push!(kt, t_on[k] / dur, t_off[k] / dur)
    end
    push!(kt, 1.0)
    keytimes = join(_fmt4.(kt), ";")
    styles = _marker_styles(stages, colors, rmax)
    anim(attr, vals) = "<animate attributeName=\"$attr\" dur=\"$(_fmt(dur))s\" " *
                       "repeatCount=\"indefinite\" keyTimes=\"$keytimes\" values=\"" *
                       join(vcat([[v, v] for v in vals]..., vals[1]), ";") * "\"/>"
    # without animation, only the first state shows
    group(k) = k == 1 ? "<g>" : "<g opacity=\"0\">"
    # a group visible while the markers move towards state k and while state k is held
    function visible(k)
        n == 1 && return ""
        ts, vs = k == 1 ? ([0, t_off[1], t_off[n]] ./ dur, (1, 0, 1)) :
                 ([0, t_off[k - 1], t_off[k]] ./ dur, (0, 1, 0))
        "<animate attributeName=\"opacity\" dur=\"$(_fmt(dur))s\" repeatCount=\"indefinite\" " *
        "calcMode=\"discrete\" keyTimes=\"$(join(_fmt4.(ts), ";"))\" values=\"$(join(vs, ";"))\"/>"
    end
    io = IOBuffer()
    print(io,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 $(_fmt(W)) $(_fmt(H))\" ",
        "width=\"100%\" style=\"max-width:$(round(Int, W))px\" font-family=\"sans-serif\">",
        _SVG_DEFS)
    _frame(io, p)
    _clip_open(io, p)
    for k in 1:n
        print(io, group(k), visible(k))
        foreach(s -> _overlay(io, p, s), _overlays(overlays, k))
        t = _targets(targets, k)
        if t !== nothing
            for (xt, vt) in zip(t...)
                print(io, "<circle cx=\"$(_fmt(_px(p, xt)))\" cy=\"$(_fmt(_py(p, vt)))\" ",
                    "r=\"$(_fmt(0.6rmax))\" fill=\"none\" stroke=\"#D55E00\"/>")
            end
        end
        print(io, "</g>")
    end
    for a in 1:N
        xs = [_fmt(_px(p, s.markers.x[a])) for s in stages]
        vs = [_fmt(_py(p, s.markers.v[a])) for s in stages]
        rs = [_fmt(styles[k][2][a]) for k in 1:n]
        cs = [styles[k][1][a] for k in 1:n]
        print(
            io, "<circle cx=\"$(xs[1])\" cy=\"$(vs[1])\" r=\"$(rs[1])\" fill=\"$(cs[1])\" ",
            "fill-opacity=\"0.85\">", anim("cx", xs), anim("cy", vs), anim("r", rs),
            anim("fill", cs), "</circle>")
    end
    print(io, "</g>")
    for k in 1:n
        print(io, group(k), visible(k), _title(p, titles[k]), "</g>")
    end
    print(io, "</svg>")
    SVGFigure(String(take!(io)))
end

_fmt4(t) = string(round(t; digits = 4))

export SVGFigure, svg_line, svg_rect, svg_text, svg_arrow, svg_category_colors,
       svg_stage_panels, svg_stage_animation
