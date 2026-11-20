# Plots of `ParticleRemapping.Projection`, in the common style of `common.jl`.

const _Pr = ParticleRemapping.Projection

function _Pr.plot_hole_outline!(ax, hole::_Pr.EllipticHole)
    hx, hv = _Pr.hole_outline(hole)
    lines!(ax, hx, hv; color = _OUTLINE_COLOR, linestyle = :dash, linewidth = 1)
end

# The indices of the markers that lie within `a` of the centre of `hole` (periodic in x) in some
# stage.
function _near_hole(stages, hole, a)
    L = hole.L
    near(mk) = (abs.(mod.(mk.x .- hole.x .+ L / 2, L) .- L / 2) .<= a) .&
               (abs.(mk.v .- hole.v) .<= a)
    findall(reduce(.|, [near(st.markers) for st in stages]))
end

function _Pr.plot_hole_stages(stages::AbstractVector{<:RemapStage}, hole::_Pr.EllipticHole;
        window = 1.2, columns = 4, kwargs...)
    k = _near_hole(stages, hole, 1.1window)
    sub(mk) = Markers(mk.x[k], mk.v[k], mk.w[k]; mk.m)
    ParticleRemapping.plot_stages([RemapStage(sub(st.markers), st.label) for st in stages];
        columns, period = (0, hole.L),
        window = (hole.x - window, hole.x + window, hole.v - window, hole.v + window),
        decorate! = (ax, _) -> _Pr.plot_hole_outline!(ax, hole), kwargs...)
end

function _Pr.plot_mean_field_step(data = _Pr.mean_field_step_example(); size = (660, 360))
    fig = Figure(; size)
    axes = Axis[]
    for (j, (title, Z, c)) in enumerate((("(a) before", data.before, data.circle),
        ("(b) after the step", data.after, data.image)))
        ax = Axis(fig[1, j]; title, aspect = DataAspect(), xlabel = L"\xi",
            ylabel = j == 1 ? L"\eta" : "", yticklabelsvisible = j == 1)
        push!(axes, ax)
        poly!(ax, Point2f.(first.(c), last.(c)); color = (_OUTLINE_COLOR, 0.15),
            strokecolor = _OUTLINE_COLOR, strokewidth = 1, linestyle = :dash)
        scatter!(ax, first.(Z), last.(Z); markersize = 5, color = _SERIES[1])
    end
    linkaxes!(axes...)
    fig
end

function _Pr.plot_implicit_weight_fold(; vmax = 4, levels = (0.5, 1.5), size = (520, 360))
    fM(r) = exp(-r^2 / 2) / sqrt(2π)
    r = range(0, vmax; length = 400)
    pmax = _Pr.maxwellian_fold_peak()
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; xlabel = L"|v| / v_{\mathrm{th}}", ylabel = L"p / p_{\max}",
        title = "Maxwellian, p = m c |v| f(v)")
    lines!(ax, r, r .* fM.(r) ./ pmax; color = _SERIES[1])
    for (k, q) in enumerate(levels)
        n = _Pr.maxwellian_fold_roots(q * pmax)
        hlines!(ax, q; color = _SERIES[k + 1], linestyle = :dash,
            label = _richtitle("p = $q p_max: $n roots"))
    end
    vlines!(ax, 1.0; color = _OUTLINE_COLOR, linestyle = :dot)
    axislegend(ax; position = :rt)
    fig
end

function _Pr.plot_pm_v_coupling(data = _Pr.pm_v_coupling(); size = (520, 360))
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; xlabel = L"v", ylabel = L"\det DN - 1",
        title = "one marker step at x = $(round(data.x; digits = 4))")
    hlines!(ax, 0; color = _OUTLINE_COLOR, linewidth = 0.5)
    lines!(ax, data.v, data.separable; color = _SERIES[1], label = "separable")
    lines!(ax, data.v, data.nonseparable; color = _SERIES[2], label = "non-separable")
    axislegend(ax; position = :ct)
    fig
end

function _Pr.plot_hole_filling(e; window = 1.2, kwargs...)
    prob, res = e.problem, e.results
    m, L, wt = prob.m, prob.L, e.targets.w
    mk(x, v, w) = Markers(mod.(x, L), v, w; m)
    wd, ns, cc = res[:weight_dependent], res[:non_separable], res[:control_C]
    stages = [RemapStage(prob.markers, "before"),
        RemapStage(mk(e.targets.x, e.targets.v, wt), "targets"),
        RemapStage(mk(prob.x0, res[:θ0].v1, wt), "reweighting only"),
        RemapStage(mk(cc.X, cc.v1, wt), "control C"),
        RemapStage(mk(wd.X, wd.v1, wt), "kicks and drifts"),
        RemapStage(mk(ns.X, ns.v1, wt), "non-separable"),
        RemapStage(mk(wd.X, wd.vc, wd.wc), "(e), charge-conserving"),
        RemapStage(mk(ns.X, ns.vc, ns.wc), "(f), charge-conserving")]
    _Pr.plot_hole_stages(stages, prob.hole; window, columns = 4, kwargs...)
end

function _Pr.plot_charge_correction(experiments...; variant = :non_separable,
        labels = ["x_h = $(round(e.problem.hole.x; digits = 3))" for e in experiments],
        size = (980, 380))
    fig = Figure(; size)
    ax1 = Axis(fig[1, 1]; xlabel = L"x", ylabel = L"(\rho(w) - \rho) / \max \rho",
        title = "(a) grid charge after the map")
    ax2 = Axis(fig[1, 2]; xlabel = L"x", ylabel = L"\ln (w' / w^*)", title = "(b) correction")
    for (k, (e, lab)) in enumerate(zip(experiments, labels))
        b, ρ0, r = e.charge.basis, e.ρ0, e.results[variant]
        xg = (0:(b.n - 1)) .* b.h
        rel(w) = (_Pr.spline_grid_charge(b, r.X, w) .- ρ0) ./ maximum(ρ0)
        lines!(ax1, xg, rel(e.targets.w); color = _SERIES[k], label = lab)
        scatter!(ax1, xg, rel(r.wc); color = _SERIES[k], markersize = 6)
        scatter!(
            ax2, mod.(r.X, e.problem.L), log.(r.wc ./ e.targets.w); color = _SERIES[k],
            markersize = 3, label = lab)
    end
    axislegend(ax1; position = :rt)
    axislegend(ax2; position = :rb)
    fig
end

function _Pr.plot_dlra_weights(e; tolerance = 1e-6, error_tolerance = 1e-4,
        tolerances = (1e-3, 1e-4, 1e-5, 1e-6), size = (980, 760))
    name(ϑ) = "ϑ = 1e" * string(round(Int, log10(ϑ)))
    prob, g, L = e.problem, e.grid, e.problem.L
    r = e.results[name(tolerance)]
    pushed = Markers(mod.(e.xT, L), e.vT, prob.w0; prob.m)
    remapped = Markers(mod.(r.X, L), r.vc, r.wc; prob.m)
    s = _weight_scale([pushed.w, remapped.w])
    hx, hv = e.hole_image
    fig = Figure(; size)
    top = fig[1, 1] = GridLayout()
    local sc
    for (j, (title, mkr)) in enumerate((("(a) pushed markers", pushed),
        ("(b) remap, $(name(tolerance))", remapped)))
        ax = Axis(top[1, j]; title, xlabel = L"x", ylabel = j == 1 ? L"v" : "",
            yticklabelsvisible = j == 1, limits = (0, L, -4, 4))
        sc = ParticleRemapping.plot_markers!(ax, mkr; s.colorrange, s.colorscale, s.wref)
        scatter!(ax, mod.(hx, L), hv; color = _OUTLINE_COLOR, markersize = 2)
    end
    Colorbar(top[1, 3], sc; label = "weight w")
    ax = Axis(top[1, 4]; title = "(c) f_r − f(T), $(name(error_tolerance))", xlabel = L"x",
        limits = (0, L, -4, 4), yticklabelsvisible = false)
    Δ = _Pr.lowrank_full(e.solutions[name(error_tolerance)].f) .- e.Fe
    mΔ = maximum(abs, Δ)
    hm = heatmap!(ax, g.x, g.v, Δ; colormap = _SIGNED_CMAP, colorrange = (-mΔ, mΔ))
    Colorbar(top[1, 5], hm)
    bottom = fig[2, 1] = GridLayout()
    ax = Axis(bottom[1, 1]; title = "(d) numerical rank", xlabel = L"t", ylabel = L"r")
    for (k, ϑ) in enumerate(tolerances)
        rk = e.solutions[name(ϑ)].ranks
        stairs!(ax, (1:length(rk)) .* e.dt, rk; color = _SERIES[k])
    end
    ax = Axis(bottom[1, 2]; title = "(e) target weights", xlabel = L"|v^*|",
        ylabel = L"|w^* / (c f) - 1|", yscale = log10, limits = (0, 4, 1e-7, 2),
        yticks = (10.0 .^ (-6:2:0), [L"10^{%$k}" for k in -6:2:0]))
    vt, wex = e.targets.v, e.results["f(T)"].wt
    for (k, ϑ) in enumerate(tolerances)
        ϑ == 1e-5 && continue
        wr = e.results[name(ϑ)].wt
        scatter!(ax, abs.(vt), max.(abs.(wr ./ wex .- 1), 1e-7); markersize = 2,
            color = _SERIES[k])
    end
    pc = _SERIES[6]
    wp = e.results["projection"].wt
    scatter!(ax, abs.(vt), max.(abs.(wp ./ wex .- 1), 1e-7); markersize = 2, color = pc)
    Legend(bottom[1, 3],
        vcat(
            [[LineElement(color = _SERIES[k]),
                 MarkerElement(color = _SERIES[k], marker = :circle, markersize = 6)]
             for k in eachindex(tolerances)],
            [MarkerElement(color = pc, marker = :circle, markersize = 6)]),
        vcat([name(ϑ) for ϑ in tolerances], ["projection"]))
    rowsize!(fig.layout, 1, Relative(0.5))
    fig
end

function _Pr.plot_lumped_projection(problem = _Pr.hole_filling_problem();
        projection = _Pr.HatProjection(32, 33), size = (1060, 330))
    hp, L, V = projection, problem.L, problem.V
    s, r = problem.x0 ./ L, problem.v0 ./ V
    fh = _Pr.lumped_projection(hp, s, r, problem.w0) ./ (L * V)
    gh = _Pr.lumped_projection(hp, s, r, ones(length(s))) ./ (L * V)
    xs = (0:(hp.nX - 1)) .* (L / hp.nX)
    vs = range(-V, V; length = hp.nV)
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; title = "(a) markers", xlabel = L"x", ylabel = L"v",
        limits = (0, L, -V, V))
    sc = ParticleRemapping.plot_markers!(ax, problem.markers)
    _Pr.plot_hole_outline!(ax, problem.hole)
    Colorbar(fig[1, 2], sc; label = "weight w")
    for (k, (title, G)) in enumerate(((L"(b) $f_h$, lumped projection of $f$", fh),
        (L"(c) $g_h$, lumped marker density", gh)))
        ax = Axis(fig[1, 2k + 1]; title, xlabel = L"x", yticklabelsvisible = false,
            limits = (0, L, -V, V))
        hm = heatmap!(ax, xs, vs, G; colormap = _FIELD_CMAP)
        _Pr.plot_hole_outline!(ax, problem.hole)
        Colorbar(fig[1, 2k + 2], hm)
    end
    colgap!(fig.layout, 2, 30)
    fig
end

function _Pr.plot_target_flow(problem, targets; window = nothing, size = (620, 480))
    L, V = problem.L, problem.V
    lim = window === nothing ? (0, L, -V, V) : window
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; title = "density-interpolation flow: markers to targets",
        xlabel = L"x", ylabel = L"v", limits = lim)
    seg = [Point2f(problem.x0[k], problem.v0[k]) => Point2f(targets.x[k], targets.v[k])
           for k in eachindex(problem.x0)]
    linesegments!(ax, seg; color = _SEGMENT_COLOR, linewidth = 0.7)
    # the marker sizes suit the markers in view
    inview = count(
        k -> lim[1] <= problem.x0[k] <= lim[2] && lim[3] <= problem.v0[k] <= lim[4],
        eachindex(problem.x0))
    _targets!(ax, (targets.x, targets.v), inview)
    w = problem.w0
    sc = ParticleRemapping.plot_markers!(ax, problem.markers;
        markersize = _weight_markersize(w, maximum(w), inview))
    _Pr.plot_hole_outline!(ax, problem.hole)
    Colorbar(fig[1, 2], sc; label = "weight w")
    fig
end

function _Pr.animate_dlra(path, frames, grid; framerate = 20, V = 4.0, size = (640, 480))
    L = length(grid.x) * grid.hx
    F = Observable(frames[1].F)
    pts = Observable(Point2f.(mod.(frames[1].x, L), frames[1].v))
    hole = Observable(Point2f.(mod.(frames[1].hole_x, L), frames[1].hole_v))
    ttl = Observable("t = 0")
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; title = ttl, xlabel = L"x", ylabel = L"v", limits = (0, L, -V, V))
    hm = heatmap!(ax, grid.x, grid.v, F; colormap = _FIELD_CMAP,
        colorrange = (0, maximum(F[])))
    scatter!(ax, pts; markersize = 2, color = (:white, 0.8))
    scatter!(ax, hole; markersize = 2, color = _TARGET_COLOR)
    Colorbar(fig[1, 2], hm; label = "f_r")
    record(fig, path, frames; framerate) do fr
        F[] = fr.F
        pts[] = Point2f.(mod.(fr.x, L), fr.v)
        hole[] = Point2f.(mod.(fr.hole_x, L), fr.hole_v)
        ttl[] = "low-rank density: t = $(round(fr.t; digits = 2)), rank $(fr.rank); markers"
    end
    path
end
