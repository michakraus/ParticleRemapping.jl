# SVG pictures of the steps of a remap: one panel per state, and JSON data for the Vue/SVG
# algorithm player. They need no plotting package.

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

The overlays of [`svg_stage_panels`](@ref) and [`svg_player_json`](@ref), in the
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
[`svg_player_json`](@ref).
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

# JSON serialisation helpers for `svg_player_json` (no JSON dependency in the package).
_json_num(x::Real) = string(round(Float64(x); digits = 4))

function _json_escape(io::IO, s::AbstractString)
    print(io, '"')
    for c in s
        if c == '"'
            print(io, "\\\"")
        elseif c == '\\'
            print(io, "\\\\")
        elseif c == '\b'
            print(io, "\\b")
        elseif c == '\f'
            print(io, "\\f")
        elseif c == '\n'
            print(io, "\\n")
        elseif c == '\r'
            print(io, "\\r")
        elseif c == '\t'
            print(io, "\\t")
        elseif UInt32(c) < 0x20
            print(io, "\\u", string(UInt32(c); base = 16, pad = 4))
        else
            print(io, c)
        end
    end
    print(io, '"')
end

function _json_hex_color(c::AbstractString, k::Int)
    s = String(c)
    if occursin(r"^#[0-9a-fA-F]{6}\z", s)
        return lowercase(s)
    end
    m = match(r"^rgb\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)\z", s)
    if m !== nothing
        rgb = [tryparse(Int, cap) for cap in m.captures]
        if all(v -> v !== nothing && 0 <= v <= 255, rgb)
            return "#" * join(string.(rgb; base = 16, pad = 2))
        end
    end
    throw(ArgumentError(
        "svg_player_json: the colour $(repr(s)) of state $k is not #rrggbb or rgb(r,g,b)"))
end

function _json_numbers(io::IO, xs)
    print(io, '[')
    for (i, x) in enumerate(xs)
        i > 1 && print(io, ',')
        print(io, _json_num(x))
    end
    print(io, ']')
end

function _json_strings(io::IO, xs)
    print(io, '[')
    for (i, x) in enumerate(xs)
        i > 1 && print(io, ',')
        _json_escape(io, x)
    end
    print(io, ']')
end

function _check_overlay(o, k::Int)
    if o.kind === :line || o.kind === :arrow
        all(isfinite, (o.x1, o.v1, o.x2, o.v2, o.width)) ||
            throw(ArgumentError("svg_player_json: a non-finite number in state $k"))
    elseif o.kind === :rect
        all(isfinite, (o.x, o.v, o.Δx, o.Δv, o.opacity)) ||
            throw(ArgumentError("svg_player_json: a non-finite number in state $k"))
    elseif o.kind === :text
        all(isfinite, (o.x, o.v, o.size)) ||
            throw(ArgumentError("svg_player_json: a non-finite number in state $k"))
    else
        throw(ArgumentError(
            "svg_player_json: an overlay of state $k has the unknown kind $(repr(o.kind))"))
    end
end

function _json_overlay(io::IO, o)
    if o.kind === :line
        print(io, "{\"kind\":\"line\",\"x1\":", _json_num(o.x1), ",\"v1\":",
            _json_num(o.v1), ",\"x2\":", _json_num(o.x2), ",\"v2\":", _json_num(o.v2),
            ",\"color\":")
        _json_escape(io, string(o.color))
        print(io, ",\"width\":", _json_num(o.width), ",\"dash\":",
            Bool(o.dash) ? "true" : "false", '}')
    elseif o.kind === :arrow
        print(io, "{\"kind\":\"arrow\",\"x1\":", _json_num(o.x1), ",\"v1\":",
            _json_num(o.v1), ",\"x2\":", _json_num(o.x2), ",\"v2\":", _json_num(o.v2),
            ",\"color\":")
        _json_escape(io, string(o.color))
        print(io, ",\"width\":", _json_num(o.width), '}')
    elseif o.kind === :rect
        print(io, "{\"kind\":\"rect\",\"x\":", _json_num(o.x), ",\"v\":", _json_num(o.v),
            ",\"dx\":", _json_num(o.Δx), ",\"dv\":", _json_num(o.Δv), ",\"color\":")
        _json_escape(io, string(o.color))
        print(io, ",\"opacity\":", _json_num(o.opacity), '}')
    else
        print(io, "{\"kind\":\"text\",\"x\":", _json_num(o.x), ",\"v\":", _json_num(o.v),
            ",\"text\":")
        _json_escape(io, string(o.text))
        print(io, ",\"color\":")
        _json_escape(io, string(o.color))
        print(io, ",\"size\":", _json_num(o.size), '}')
    end
end

"""
    svg_player_json(stages; window, id, overlays = nothing, targets = nothing,
                    colors = nothing, rmax = 7, titles = [s.label for s in stages]) -> String

Serialise the states of `stages` (a vector of [`RemapStage`](@ref)) as a UTF-8 JSON object for the
Vue/SVG `AlgorithmPlayer` component, over the box `window = (xmin, xmax, vmin, vmax)` with key
`id`. The markers share one colour range and one size reference across all states (radii in pixels
for a `420 × 300` panel up to `rmax`, colours normalised to lowercase `#rrggbb`). `overlays`,
`targets`, `colors` and `titles` accept the same forms as [`svg_stage_panels`](@ref). Every numeric
value is checked finite and rounded to four decimals; invalid inputs throw `ArgumentError` or
`DimensionMismatch`.
"""
function svg_player_json(
        stages::AbstractVector{<:RemapStage}; window, id::AbstractString,
        overlays = nothing, targets = nothing, colors = nothing, rmax::Real = 7,
        titles = [s.label for s in stages])
    Base.require_one_based_indexing(stages)
    n = length(stages)
    n > 0 || throw(ArgumentError("svg_player_json: no states"))
    win = Float64.(Tuple(window))
    (length(win) == 4 && all(isfinite, win) && win[1] < win[2] && win[3] < win[4]) ||
        throw(ArgumentError(
            "svg_player_json: the window must hold finite bounds with xmin < xmax and vmin < vmax"))
    rmax_f = Float64(rmax)
    (isfinite(rmax_f) && rmax_f > 0) ||
        throw(ArgumentError("svg_player_json: rmax must be finite and positive, got $rmax"))
    N = length(stages[1].markers)
    all(s -> length(s.markers) == N, stages) ||
        throw(DimensionMismatch("svg_player_json: every state must hold the same markers"))
    N > 0 || throw(ArgumentError("svg_player_json: the states hold no markers"))
    length(titles) == n ||
        throw(DimensionMismatch("svg_player_json: titles must hold one caption per state"))
    if targets isa AbstractVector && length(targets) != n
        throw(DimensionMismatch("svg_player_json: targets must hold one entry per state"))
    end
    if overlays isa AbstractVector && length(overlays) != n
        throw(DimensionMismatch("svg_player_json: overlays must hold one entry per state"))
    end

    captions = Vector{String}(undef, n)
    step_targets = Vector{Any}(undef, n)
    step_overlays = Vector{Any}(undef, n)
    for k in 1:n
        cap = string(titles[k])
        isempty(strip(cap)) &&
            throw(ArgumentError("svg_player_json: the caption of state $k is empty"))
        captions[k] = cap
        mk = stages[k].markers
        (all(isfinite, mk.x) && all(isfinite, mk.v) && all(isfinite, mk.w)) ||
            throw(ArgumentError("svg_player_json: a non-finite number in state $k"))
        t = _targets(targets, k)
        if t !== nothing
            (length(t) == 2 && length(t[1]) == N && length(t[2]) == N) ||
                throw(DimensionMismatch(
                    "svg_player_json: the targets of state $k must hold one point per marker"))
            (all(isfinite, t[1]) && all(isfinite, t[2])) ||
                throw(ArgumentError("svg_player_json: a non-finite number in state $k"))
        end
        step_targets[k] = t
        ovs = collect(_overlays(overlays, k))
        for o in ovs
            _check_overlay(o, k)
        end
        step_overlays[k] = ovs
    end

    styles = _marker_styles(stages, colors, rmax_f)
    hex_colors = Vector{Vector{String}}(undef, n)
    for k in 1:n
        fillc = styles[k][1]
        length(fillc) == N ||
            throw(DimensionMismatch(
                "svg_player_json: the colours of state $k must hold one colour per marker"))
        hex_colors[k] = [_json_hex_color(c, k) for c in fillc]
    end

    io = IOBuffer()
    print(io, "{\"version\":1,\"id\":")
    _json_escape(io, id)
    print(io, ",\"window\":")
    _json_numbers(io, win)
    print(io, ",\"markers\":", N, ",\"steps\":[")
    for k in 1:n
        k > 1 && print(io, ',')
        mk = stages[k].markers
        print(io, "{\"caption\":")
        _json_escape(io, captions[k])
        print(io, ",\"x\":")
        _json_numbers(io, mk.x)
        print(io, ",\"v\":")
        _json_numbers(io, mk.v)
        print(io, ",\"r\":")
        _json_numbers(io, styles[k][2])
        print(io, ",\"c\":")
        _json_strings(io, hex_colors[k])
        print(io, ",\"targets\":")
        t = step_targets[k]
        if t === nothing
            print(io, "null")
        else
            print(io, "{\"x\":")
            _json_numbers(io, t[1])
            print(io, ",\"v\":")
            _json_numbers(io, t[2])
            print(io, '}')
        end
        print(io, ",\"overlays\":[")
        for (i, o) in enumerate(step_overlays[k])
            i > 1 && print(io, ',')
            _json_overlay(io, o)
        end
        print(io, "]}")
    end
    print(io, "]}")
    return String(take!(io))
end

export SVGFigure, svg_line, svg_rect, svg_text, svg_arrow, svg_category_colors,
       svg_stage_panels, svg_player_json
