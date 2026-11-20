# Plots of `ParticleRemapping.Projection`.

const _PC = Makie.wong_colors()

function _dots!(ax, x, v; L = 2π, markersize = 2.5, kwargs...)
    scatter!(ax, mod.(x, L), v;
        markersize, color = :black, kwargs...)
end

function _outline!(ax, hole)
    hx, hv = ParticleRemapping.Projection.hole_outline(hole)
    lines!(ax, hx, hv; color = _PC[6], linewidth = 1)
end

function ParticleRemapping.Projection.plot_mean_field_step(
        data = ParticleRemapping.Projection.mean_field_step_example(); size = (660, 360))
    fig = Figure(; size)
    for (j, (title, Z, c)) in enumerate((("(a) before", data.before, data.circle),
        ("(b) after the step", data.after, data.image)))
        ax = Axis(fig[1, j]; title, aspect = DataAspect(), xlabel = "ξ",
            ylabel = j == 1 ? "η" : "", yticklabelsvisible = j == 1)
        poly!(
            ax, Point2f.(first.(c), last.(c)); color = (_PC[6], 0.3), strokecolor = _PC[6],
            strokewidth = 1)
        scatter!(ax, first.(Z), last.(Z); markersize = 4, color = :black)
    end
    fig
end

function ParticleRemapping.Projection.plot_implicit_weight_fold(;
        vmax = 4, levels = (0.5, 1.5),
        size = (520, 360))
    fM(r) = exp(-r^2 / 2) / sqrt(2π)
    r = range(0, vmax; length = 400)
    pmax = ParticleRemapping.Projection.maxwellian_fold_peak()
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; xlabel = "|v| / v_th", ylabel = "p / p_max",
        title = "p = m c |v| f_M(v), Maxwellian")
    lines!(ax, r, r .* fM.(r) ./ pmax; color = :black)
    for (k, q) in enumerate(levels)
        n = ParticleRemapping.Projection.maxwellian_fold_roots(q * pmax)
        hlines!(ax, q; color = _PC[k], linestyle = :dash, label = "p = $q p_max: $n roots")
    end
    vlines!(ax, 1.0; color = :gray60, linestyle = :dot)
    axislegend(ax; position = :rt)
    fig
end

function ParticleRemapping.Projection.plot_pm_v_coupling(
        data = ParticleRemapping.Projection.pm_v_coupling(); size = (520, 360))
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; xlabel = "v", ylabel = "det DN − 1",
        title = "one marker step at x = $(round(data.x; digits = 4))")
    hlines!(ax, 0; color = :gray70, linewidth = 0.5)
    lines!(ax, data.v, data.separable; color = _PC[1], label = "separable")
    lines!(ax, data.v, data.nonseparable; color = _PC[2], label = "non-separable")
    axislegend(ax; position = :ct)
    fig
end

function ParticleRemapping.Projection.plot_hole_filling(e; window = 1.2, size = (980, 520))
    prob, res = e.problem, e.results
    xh = prob.hole.x
    panels = [("(a) before", prob.x0, prob.v0), ("(b) targets", e.targets.x, e.targets.v),
        ("(c) reweighting only", prob.x0, res[:θ0].v1),
        ("(d) control C", res[:control_C].X, res[:control_C].v1),
        ("(e) kicks and drifts", res[:weight_dependent].X, res[:weight_dependent].v1),
        ("(f) non-separable", res[:non_separable].X, res[:non_separable].v1),
        ("(g) (e) + charge", res[:weight_dependent].X, res[:weight_dependent].vc),
        ("(h) (f) + charge", res[:non_separable].X, res[:non_separable].vc)]
    fig = Figure(; size)
    for (j, (title, x, v)) in enumerate(panels)
        r, c = fldmod1(j, 4)
        ax = Axis(fig[r, c]; title, limits = (xh - window, xh + window, 0.0, 2.4),
            xlabel = r == 2 ? "x" : "", ylabel = c == 1 ? "v" : "",
            xticklabelsvisible = r == 2, yticklabelsvisible = c == 1)
        _outline!(ax, prob.hole)
        _dots!(ax, x, v; L = prob.L, markersize = 4)
    end
    fig
end

function ParticleRemapping.Projection.plot_charge_correction(experiments...;
        variant = :non_separable,
        labels = ["x_h = $(round(e.problem.hole.x; digits = 3))" for e in experiments],
        size = (980, 380))
    fig = Figure(; size)
    ax1 = Axis(fig[1, 1]; xlabel = "x", ylabel = "(ρ(w) − ρ) / max ρ",
        title = "(a) grid charge after the map")
    ax2 = Axis(fig[1, 2]; xlabel = "x", ylabel = "ln(w′ / w*)", title = "(b) correction")
    for (k, (e, lab)) in enumerate(zip(experiments, labels))
        b, ρ0, r = e.charge.basis, e.ρ0, e.results[variant]
        xg = (0:(b.n - 1)) .* b.h
        rel(w) = (ParticleRemapping.Projection.spline_grid_charge(b, r.X, w) .- ρ0) ./
                 maximum(ρ0)
        lines!(ax1, xg, rel(e.targets.w); color = _PC[k], label = lab)
        scatter!(ax1, xg, rel(r.wc); color = _PC[k], markersize = 6)
        scatter!(ax2, mod.(r.X, e.problem.L), log.(r.wc ./ e.targets.w); color = _PC[k],
            markersize = 3, label = lab)
    end
    axislegend(ax1; position = :rt)
    axislegend(ax2; position = :rb)
    fig
end

function ParticleRemapping.Projection.plot_dlra_weights(e; tolerance = 1e-6,
        error_tolerance = 1e-4, tolerances = (1e-3, 1e-4, 1e-5, 1e-6), size = (980, 760))
    name(ϑ) = "ϑ = 1e" * string(round(Int, log10(ϑ)))
    prob, g, L = e.problem, e.grid, e.problem.L
    fig = Figure(; size)
    top = fig[1, 1] = GridLayout()
    hx, hv = e.hole_image
    ax = Axis(top[1, 1]; title = "(a) pushed markers", xlabel = "x", ylabel = "v",
        limits = (0, L, -4, 4))
    _dots!(ax, e.xT, e.vT; L, markersize = 2)
    scatter!(ax, mod.(hx, L), hv; color = _PC[6], markersize = 2)
    r = e.results[name(tolerance)]
    ax = Axis(top[1, 2]; title = "(b) remap, $(name(tolerance))", xlabel = "x",
        limits = (0, L, -4, 4), yticklabelsvisible = false)
    _dots!(ax, r.X, r.vc; L, markersize = 2)
    scatter!(ax, mod.(hx, L), hv; color = _PC[6], markersize = 2)
    ax = Axis(top[1, 3]; title = "(c) f_r − f(T), $(name(error_tolerance))", xlabel = "x",
        limits = (0, L, -4, 4), yticklabelsvisible = false)
    Δ = ParticleRemapping.Projection.lowrank_full(e.solutions[name(error_tolerance)].f) .-
        e.Fe
    mΔ = maximum(abs, Δ)
    hm = heatmap!(ax, g.x, g.v, Δ; colormap = :balance, colorrange = (-mΔ, mΔ))
    Colorbar(top[1, 4], hm; width = 8)
    bottom = fig[2, 1] = GridLayout()
    ax = Axis(bottom[1, 1]; title = "(d) numerical rank", xlabel = "t", ylabel = "r")
    for (k, ϑ) in enumerate(tolerances)
        rk = e.solutions[name(ϑ)].ranks
        stairs!(ax, (1:length(rk)) .* e.dt, rk; color = _PC[k])
    end
    ax = Axis(bottom[1, 2]; title = "(e) target weights", xlabel = "|v*|",
        ylabel = "|w* / (c f) − 1|", yscale = log10, limits = (0, 4, 1e-7, 2))
    vt, wex = e.targets.v, e.results["f(T)"].wt
    for (k, ϑ) in enumerate(tolerances)
        ϑ == 1e-5 && continue
        wr = e.results[name(ϑ)].wt
        scatter!(ax, abs.(vt), max.(abs.(wr ./ wex .- 1), 1e-7); markersize = 2,
            color = _PC[k])
    end
    wp = e.results["projection"].wt
    scatter!(
        ax, abs.(vt), max.(abs.(wp ./ wex .- 1), 1e-7); markersize = 2, color = :gray50)
    Legend(bottom[1, 3],
        vcat(
            [[LineElement(color = _PC[k]),
                 MarkerElement(color = _PC[k], marker = :circle, markersize = 6)]
             for k in eachindex(tolerances)],
            [MarkerElement(color = :gray50, marker = :circle, markersize = 6)]),
        vcat([name(ϑ) for ϑ in tolerances], ["projection"]); framevisible = false)
    rowsize!(fig.layout, 1, Relative(0.5))
    fig
end

function ParticleRemapping.Projection.plot_lumped_projection(
        problem = ParticleRemapping.Projection.hole_filling_problem();
        projection = ParticleRemapping.Projection.HatProjection(32, 33), size = (980, 330))
    P = ParticleRemapping.Projection
    hp, L, V = projection, problem.L, problem.V
    s, r = problem.x0 ./ L, problem.v0 ./ V
    fh = P.lumped_projection(hp, s, r, problem.w0) ./ (L * V)
    gh = P.lumped_projection(hp, s, r, ones(length(s))) ./ (L * V)
    xs = (0:(hp.nX - 1)) .* (L / hp.nX)
    vs = range(-V, V; length = hp.nV)
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; title = "markers (area ∝ weight)", xlabel = "x", ylabel = "v",
        limits = (0, L, -V, V))
    scatter!(ax, problem.x0, problem.v0;
        markersize = 1 .+ 4 .* sqrt.(problem.w0 ./
                                     maximum(problem.w0)),
        color = :black)
    _outline!(ax, problem.hole)
    for (k, (title, G)) in enumerate((("f_h, lumped projection of f", fh),
        ("g_h, lumped marker density", gh)))
        ax = Axis(fig[1, 2k]; title, xlabel = "x", yticklabelsvisible = false,
            limits = (0, L, -V, V))
        hm = heatmap!(ax, xs, vs, G; colormap = :viridis)
        _outline!(ax, problem.hole)
        Colorbar(fig[1, 2k + 1], hm; width = 8)
    end
    fig
end

function ParticleRemapping.Projection.plot_target_flow(problem, targets; window = nothing,
        size = (620, 480))
    L, V = problem.L, problem.V
    lim = window === nothing ? (0, L, -V, V) : window
    fig = Figure(; size)
    ax = Axis(
        fig[1, 1]; title = "density-interpolation flow: markers → targets", xlabel = "x",
        ylabel = "v", limits = lim)
    seg = [Point2f(problem.x0[k], problem.v0[k]) => Point2f(targets.x[k], targets.v[k])
           for k in eachindex(problem.x0)]
    linesegments!(ax, seg; color = (:gray40, 0.7), linewidth = 0.7)
    scatter!(ax, problem.x0, problem.v0; markersize = 3, color = :black, label = "markers")
    scatter!(ax, targets.x, targets.v; markersize = 3, color = _PC[1], label = "targets")
    _outline!(ax, problem.hole)
    axislegend(ax; position = :rb)
    fig
end

function ParticleRemapping.Projection.plot_remap_stages(states, hole; window = nothing,
        columns = 4, size = (980, 260 * cld(length(states), columns)))
    lim = window === nothing ? (hole.x - 1.2, hole.x + 1.2, 0.0, 2.4) : window
    fig = Figure(; size)
    for (j, (x, v, title)) in enumerate(states)
        r, c = fldmod1(j, columns)
        ax = Axis(fig[r, c]; title, limits = lim, xlabel = "x", ylabel = c == 1 ? "v" : "",
            yticklabelsvisible = c == 1)
        _outline!(ax, hole)
        _dots!(ax, x, v; L = hole.L, markersize = 4)
    end
    fig
end

function ParticleRemapping.Projection.animate_remap(path, frames, hole; framerate = 20,
        V = 4.0, size = (640, 480))
    pts = Observable(Point2f.(mod.(frames[1][1], hole.L), frames[1][2]))
    ttl = Observable(frames[1][3])
    fig = Figure(; size)
    ax = Axis(
        fig[1, 1]; title = ttl, xlabel = "x", ylabel = "v", limits = (0, hole.L, -V, V))
    _outline!(ax, hole)
    scatter!(ax, pts; markersize = 3, color = :black)
    record(fig, path, frames; framerate) do (x, v, t)
        pts[] = Point2f.(mod.(x, hole.L), v)
        ttl[] = t
    end
    path
end

function ParticleRemapping.Projection.animate_moser_targets(
        path, history, hole; framerate = 15,
        V = 4.0, size = (640, 480))
    pts = Observable(Point2f.(mod.(history[1][1], hole.L), history[1][2]))
    ttl = Observable("t = 0")
    fig = Figure(; size)
    ax = Axis(
        fig[1, 1]; title = ttl, xlabel = "x", ylabel = "v", limits = (0, hole.L, -V, V))
    _outline!(ax, hole)
    scatter!(ax, pts; markersize = 3, color = :black)
    n = length(history)
    record(fig, path, eachindex(history); framerate) do i
        pts[] = Point2f.(mod.(history[i][1], hole.L), history[i][2])
        ttl[] = "density-interpolation flow, t = $(round((i - 1) / (n - 1); digits = 3))"
    end
    path
end

function ParticleRemapping.Projection.animate_dlra(path, frames, grid; framerate = 20,
        V = 4.0, size = (640, 480))
    L = length(grid.x) * grid.hx
    F = Observable(frames[1].F)
    pts = Observable(Point2f.(mod.(frames[1].x, L), frames[1].v))
    hole = Observable(Point2f.(mod.(frames[1].hole_x, L), frames[1].hole_v))
    ttl = Observable("t = 0")
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; title = ttl, xlabel = "x", ylabel = "v", limits = (0, L, -V, V))
    hm = heatmap!(ax, grid.x, grid.v, F; colormap = :Blues, colorrange = (0, maximum(F[])))
    scatter!(ax, pts; markersize = 2, color = :black)
    scatter!(ax, hole; markersize = 2, color = _PC[6])
    Colorbar(fig[1, 2], hm; width = 8)
    record(fig, path, frames; framerate) do fr
        F[] = fr.F
        pts[] = Point2f.(mod.(fr.x, L), fr.v)
        hole[] = Point2f.(mod.(fr.hole_x, L), fr.hole_v)
        ttl[] = "low-rank density: t = $(round(fr.t; digits = 2)), rank $(fr.rank); markers"
    end
    path
end
