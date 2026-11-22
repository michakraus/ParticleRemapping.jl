# The graphical introductions of the overview pages: the steps of each algorithm on a small
# example, as one SVG panel per step and an SVG animation through the steps.

using .Projection: EllipticHole, hole_filling_problem
using .Redistribution: redistribution_example
using .Remapping: eight_marker_example, quantile_matching_setup
using .SlicedRedistribution: sliced_stages_example

# arrows from the markers of `a` to the points `(x, v)`, for the markers that move visibly; with
# a `period`, each arrow takes the shorter way round in x
function _arrows(a::Markers, x, v; tol = 1e-3, period = nothing)
    dx = period === nothing ? x .- a.x : mod.(x .- a.x .+ period / 2, period) .- period / 2
    [svg_arrow(a.x[k], a.v[k], a.x[k] + dx[k], v[k])
     for k in eachindex(x)
     if hypot(dx[k], v[k] - a.v[k]) > tol]
end

# the outline of an elliptic hole as a closed polygon of lines
function _ellipse(h::EllipticHole; n = 48)
    t = range(0, 2π; length = n + 1)
    xs, vs = h.x .+ h.a .* cos.(t), h.v .+ h.b .* sin.(t)
    [svg_line(xs[k], vs[k], xs[k + 1], vs[k + 1]; color = "#404040", width = 1.2)
     for k in 1:n]
end

_vlines(xs, v1, v2; kw...) = [svg_line(x, v1, x, v2; kw...) for x in xs]
_hlines(vs, x1, x2; kw...) = [svg_line(x1, v, x2, v; kw...) for v in vs]

# the slice of each marker of a vector of index sets
function _slice_of(bins, N)
    s = zeros(Int, N)
    for (b, idx) in enumerate(bins)
        s[idx] .= b
    end
    s
end

function _projection_steps()
    prob = hole_filling_problem(; hole_x = π / 2, nL = 10, hole_a = 0.9, hole_b = 0.8)
    hats = HatProjection(8, 9)
    method = KickDriftRemap(; projection = hats, maxit = 20,
        density = (xs, vs) -> prob.density.(xs, vs))
    mk = prob.markers
    r = remap(method, mk)
    S = remap_states(method, mk, r)
    tg = (r.info.targets.x, r.info.targets.v)
    grid = vcat(
        _vlines(range(0, 2π; length = hats.nX + 1), -4, 4; color = "#b0b0b0",
            dash = true),
        _hlines(range(-4, 4; length = hats.nV), 0, 2π; color = "#b0b0b0",
            dash = true))
    stages = [RemapStage(mk, "1. densities on hat functions"),
        RemapStage(mk, "2. targets of the interpolation flow"),
        RemapStage(S[2].markers, "3. fitted map: kick 1"),
        RemapStage(S[3].markers, "3. fitted map: drift 1"),
        RemapStage(S[4].markers, "3. fitted map: kick 2"),
        RemapStage(S[5].markers, "3. fitted map: drift 2"),
        RemapStage(S[6].markers, "4. gauge step: target weights")]
    hole = _ellipse(prob.hole)
    overlays = [vcat(grid, hole), vcat(hole, _arrows(mk, tg...; period = 2π)), hole, hole,
        hole, hole, hole]
    targets = [nothing, tg, tg, tg, tg, tg, tg]
    (; stages, overlays, targets, colors = nothing, window = (0, 2π, -4, 4), columns = 4,
        animated = 1:7)
end

function _spline_steps()
    ex = redistribution_example()
    method = SplineRedistribution()
    mk = ex.markers
    r = remap(method, mk)
    S = remap_states(method, mk, r)
    tg = (r.info.targets_x, r.info.targets_v)
    V = method.Vmax
    grid = vcat(
        _vlines(range(0, ex.Lx; length = method.nx + 1), -V, V; color = "#b0b0b0",
            dash = true),
        _hlines(range(-V, V; length = method.nv - method.degree + 1), 0, ex.Lx;
            color = "#b0b0b0", dash = true))
    at_targets = Markers(tg[1], tg[2], r.info.target_weights; m = mk.m)
    stages = [RemapStage(mk, "1. density on a spline basis"),
        RemapStage(mk, "2. relaxed targets"),
        RemapStage(at_targets, "3. KKT weights at the targets"),
        RemapStage(S[3].markers, "4. symplectic step, old weights"),
        RemapStage(S[end].markers, "5. gauge step: KKT weights")]
    overlays = [
        grid, _arrows(mk, tg...), [], _arrows(mk, S[3].markers.x, S[3].markers.v), []]
    targets = [nothing, tg, tg, tg, tg]
    (; stages, overlays, targets, colors = nothing,
        window = (0, ex.Lx, -V, V), columns = 5,
        animated = [1, 2, 4, 5])
end

function _gauge_steps()
    ex = eight_marker_example()
    mk = ex.markers
    r = remap(ex.method, mk)
    S = remap_states(ex.method, mk, r)
    vstar = mk.v .* exp.(r.info.μstar)
    vmax = 1.1 * maximum(maximum(abs, s.markers.v) for s in S)
    cells = _vlines(0:(ex.grid.Δx):(ex.grid.L), 0, vmax; color = "#404040")
    rays = [svg_line(x, 0, x, vmax; color = "#d0d0d0") for x in mk.x]
    base = vcat(cells, rays)
    stages = [RemapStage(mk, "before"),
        RemapStage(mk, "1. target speeds per cell"),
        [RemapStage(s.markers, "2. " * s.label) for s in S[2:(end - 1)]]...,
        RemapStage(S[end].markers, "3. gauge step W = w e^(−μ)")]
    n = length(stages)
    targets = vcat([nothing], fill((mk.x, vstar), n - 1))
    overlays = fill(base, n)
    (; stages, overlays, targets, colors = nothing, window = (0, ex.grid.L, 0, vmax),
        columns = n, animated = 1:n)
end

function _transport_steps()
    q = quantile_matching_setup(; N = 48, Ng = 2)
    method, mk = q.admissible, q.markers
    r = remap(method, mk)
    S = remap_states(method, mk, r)
    moved, matched = S[2].markers, S[3].markers
    L = q.grid.L
    groups = [(mod(moved.x[a], L) < L / 2 ? 0 : 2) + (moved.v[a] > 0 ? 1 : 2)
              for a in eachindex(moved.x)]
    lines = vcat(_vlines((0, L / 2, L), -3, 3; color = "#404040"),
        _hlines((0.0,), 0, L; color = "#404040"))
    stages = [RemapStage(mk, "before"),
        RemapStage(moved, "1. canonical map ψ, old weights"),
        RemapStage(moved, "2. groups of cell and sign of P"),
        RemapStage(matched, "3. gauge step: quantile velocities")]
    overlays = [lines, vcat(lines, _arrows(mk, moved.x, moved.v; tol = 0.05, period = L)),
        lines, vcat(lines, _arrows(moved, matched.x, matched.v; tol = 0.05))]
    colors = k -> k == 3 ? svg_category_colors(groups) : nothing
    (; stages, overlays, targets = nothing, colors, window = (0, L, -3, 3), columns = 4,
        animated = 1:4)
end

function _sliced_steps()
    ex = sliced_stages_example(; Nx = 8, Nv = 5, B_v = 4, B_x = 8)
    method = SlicedQuantileRemap(; B_v = 4, B_x = 8)
    mk = ex.before
    S = remap_states(method, mk, ex.result)
    plan = ex.result.info.plan
    N = length(mk)
    sv = _slice_of(plan.bins_v, N)
    sx = _slice_of(plan.bins_x, N)
    Lx, V = method.Lx, 2.5
    # the slice boundaries: between the largest and the smallest coordinate of two slices
    vb = [(maximum(mk.v[plan.bins_v[b]]) + minimum(mk.v[plan.bins_v[b + 1]])) / 2
          for b in 1:(length(plan.bins_v) - 1)]
    q̃ = S[2].markers.x
    xb = [(maximum(q̃[plan.bins_x[b]]) + minimum(q̃[plan.bins_x[b + 1]])) / 2
          for b in 1:(length(plan.bins_x) - 1)]
    xtarget = zeros(N)
    for (b, idx) in enumerate(plan.bins_v)
        xtarget[idx] .= plan.taus_x[b]
    end
    vlines_ = _hlines(vb, -Lx, Lx; color = "#808080", dash = true)
    xlines_ = _vlines(xb, -V, V; color = "#808080", dash = true)
    stages = [RemapStage(mk, "1. velocity slices"),
        RemapStage(S[2].markers, "2. Ψ_q: positions to quantiles"),
        RemapStage(S[2].markers, "3. position slices"),
        RemapStage(S[3].markers, "4. Ψ_p: momenta to quantiles"),
        RemapStage(S[4].markers, "5. T_μ: weight squeeze"),
        RemapStage(S[5].markers, "6. gauge step w̄ = w_ref μ")]
    overlays = [vlines_, vlines_, xlines_, xlines_, [], []]
    targets = [(xtarget, mk.v), (xtarget, mk.v), nothing, nothing, nothing, nothing]
    colors = k -> k <= 2 ? svg_category_colors(sv) :
                  k <= 4 ? svg_category_colors(sx) :
                  nothing
    (; stages, overlays, targets, colors, window = (-Lx, Lx, -V, V), columns = 3,
        animated = 1:6)
end

"""
    algorithm_svg(key) -> (; panels, animation)

The graphical introduction to the algorithm of the construction `key` on a small example: one SVG
panel for each step of the algorithm (`panels`) and an SVG animation through the steps
(`animation`), both [`SVGFigure`](@ref)s. The keys and the examples:

- `:projection`: the kick–drift map on 96 markers of a ``10 \\times 10`` lattice with an
  elliptic hole ([`hole_filling_problem`](@ref)); hat grid, target points, the two kicks and
  drifts, the gauge step;
- `:spline`: the spline-projection redistribution of [`redistribution_example`](@ref); spline
  grid, relaxed targets, KKT weights, the symplectic step, the gauge step;
- `:gauge`: the gauge rezoning of [`eight_marker_example`](@ref); cells and rays, target speeds,
  homotopy, gauge step;
- `:transport`: the canonical transport of 48 markers ([`quantile_matching_setup`](@ref) with two
  cells and the admissible lift); the canonical map, the groups of cell and sign, the gauge step;
- `:sliced`: the sliced redistribution of 40 markers after the long run
  ([`sliced_stages_example`](@ref) with 4 velocity and 8 position slices); the slices, the two
  lifts, the squeeze and the gauge step.
"""
function algorithm_svg(key::Symbol)
    s = key === :projection ? _projection_steps() :
        key === :spline ? _spline_steps() :
        key === :gauge ? _gauge_steps() :
        key === :transport ? _transport_steps() :
        key === :sliced ? _sliced_steps() :
        throw(ArgumentError("unknown construction :$key; the keys are :projection, :spline, " *
                            ":gauge, :transport and :sliced"))
    panels = svg_stage_panels(
        s.stages; s.window, s.columns, s.overlays, s.targets, s.colors)
    a = s.animated
    animation = svg_stage_animation(s.stages[a]; s.window, overlays = s.overlays[a],
        targets = s.targets isa AbstractVector ? s.targets[a] : s.targets,
        colors = s.colors isa Function ? (k -> s.colors(a[k])) : s.colors)
    (; panels, animation)
end

export algorithm_svg
