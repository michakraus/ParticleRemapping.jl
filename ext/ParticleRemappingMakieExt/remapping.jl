# Plots of `ParticleRemapping.Remapping`, in the common style of `common.jl`.

const _RM = ParticleRemapping.Remapping

const _RM_PITICKS = (0:(π / 2):(2π), ["0", "π/2", "π", "3π/2", "2π"])

# the colours of the weights `w` on the scale `s` of `_weight_scale`, for patches
function _rm_weight_colors(w, s)
    lo, hi = s.colorscale.(s.colorrange)
    g = Makie.cgrad(_WEIGHT_CMAP)
    [get(g, (s.colorscale(wa) - lo) / (hi - lo)) for wa in w]
end

# ------------------------------------------------------------------ the charge-constraint example

function _RM.plot_charge_example(ex = _RM.charge_constraint_example(); size = (880, 660))
    βf = _RM.charge_fold(ex)
    e1 = _RM.charge_stationary_point(ex, 1.0)
    v, w, μs = ex.v, ex.w, ex.μstar
    n = length(v)
    ca = _category_colors(n)
    wref = maximum(vcat(w, e1.W, w .* exp.(-μs)))
    ms(ww) = _weight_markersize(ww, wref, n)
    vmax = 1.1maximum(vcat(v, ex.Vstar, e1.V))

    f = Figure(; size)
    ax = Axis(f[1, 1:2]; xlabel = L"|v|", title = "(a) speeds and weights",
        yticks = (1:3, ["result, β = 1", "targets", "before"]),
        limits = (0, vmax, 0.4, 3.6))
    linesegments!(ax, vec([Point2f(v[a], 3) => Point2f(ex.Vstar[a], 2) for a in 1:n]);
        color = _SEGMENT_COLOR, linestyle = :dot)
    linesegments!(ax, vec([Point2f(v[a], 3) => Point2f(e1.V[a], 1) for a in 1:n]);
        color = _SEGMENT_COLOR)
    scatter!(ax, v, fill(3, n); markersize = ms(w), color = ca)
    scatter!(ax, ex.Vstar, fill(2, n); markersize = ms(w .* exp.(-μs)),
        color = :transparent, strokecolor = _TARGET_COLOR, strokewidth = 1.5)
    scatter!(ax, e1.V, fill(1, n); markersize = ms(e1.W), color = ca)

    ax2 = Axis(f[2, 1]; xlabel = L"A", ylabel = L"c^{(0)}", title = "(b) constraint")
    for (j, β) in enumerate((1.0, βf, 1.6))
        As = range(_RM.charge_multiplier_limit(ex, β), 0.3; length = 300)
        lines!(ax2, As, [_RM.charge_constraint_residual(ex, A, β) for A in As];
            color = _SERIES[4 + j],
            label = β == βf ? "β = $(round(β; digits = 3)) (fold)" :
                    "β = $(round(β; digits = 2))")
    end
    hlines!(ax2, 0; color = _OUTLINE_COLOR, linewidth = 0.5)
    vlines!(ax2, -exp(-1) / ex.s2; color = _OUTLINE_COLOR, linestyle = :dash)
    scatter!(ax2, [e1.A], [0.0]; color = _SERIES[5], markersize = 9)
    axislegend(ax2; position = :rt)

    βs = range(0, βf; length = 200)
    sols = [_RM.charge_stationary_point(ex, β) for β in βs]
    ax3 = Axis(f[2, 2]; xlabel = L"\beta", ylabel = L"H_a / w_a", title = "(c) curvature")
    for a in 1:n
        lines!(ax3, βs, [s.Hw[a] for s in sols]; color = ca[a], label = "a = $a")
    end
    vlines!(ax3, βf; color = _OUTLINE_COLOR, linestyle = :dash)
    ylims!(ax3, 0, 1.05 / ex.s2)
    axislegend(ax3; position = :lb)
    f
end

function _RM.animate_charge_fold(path, ex = _RM.charge_constraint_example(); framerate = 20)
    βf = _RM.charge_fold(ex)
    v, w, μs = ex.v, ex.w, ex.μstar
    n, ca = length(v), _category_colors(length(v))
    βs = range(0, βf; length = 200)
    sols = [_RM.charge_stationary_point(ex, β) for β in βs]
    wref = maximum(maximum(s.W) for s in sols)
    ms(ww) = _weight_markersize(ww, wref, n)
    βa = Observable(0.0)
    sol = lift(β -> _RM.charge_stationary_point(ex, β), βa)
    fa = Figure(; size = (720, 450))
    axa = Axis(fa[1, 1]; xlabel = L"|v|",
        yticks = (1:2, ["stationary point", "targets β μ*"]),
        limits = (0, 1.4maximum(v .* exp.(βf .* μs)), 0.4, 2.6),
        title = lift(β -> "charge constraint only, β = $(round(β; digits = 3))", βa))
    scatter!(axa, lift(β -> v .* exp.(β .* μs), βa), fill(2, n);
        markersize = lift(β -> ms(w .* exp.(-β .* μs)), βa),
        color = :transparent, strokecolor = _TARGET_COLOR, strokewidth = 1.5)
    scatter!(axa, lift(s -> s.V, sol), fill(1, n); markersize = lift(s -> ms(s.W), sol),
        color = ca)
    axb = Axis(fa[2, 1]; xlabel = L"\beta", ylabel = L"H_a/w_a",
        limits = (0, βf * 1.02, 0, 1.05 / ex.s2))
    for a in 1:n
        lines!(axb, βs, [s.Hw[a] for s in sols]; color = (ca[a], 0.35))
    end
    scatter!(
        axb, lift(β -> fill(β, n), βa), lift(s -> s.Hw, sol); color = ca, markersize = 9)
    vlines!(axb, βf; color = _OUTLINE_COLOR, linestyle = :dash)
    record(fa, path, vcat(fill(0.0, 8), range(0, βf; length = 120), fill(βf, 30));
        framerate) do β
        βa[] = β
    end
    path
end

# ----------------------------------------------------------------------- the eight-marker example

# the eight markers before (open circles, sized by weight) and at `V`, `W` (coloured and sized by
# weight), with the targets `Vs` and the nodes
function _rm_eight_markers!(ax, mk, grid, V, W, Vs, s)
    x = mk.x
    vlines!(ax, collect(range(0, grid.L; step = grid.Δx)); color = (_OUTLINE_COLOR, 0.3),
        linestyle = :dash)
    scatter!(ax, x, Vs; marker = :hline, markersize = 18, color = _TARGET_COLOR)
    linesegments!(ax,
        lift(
            V -> vec([Point2f(x[a], mk.v[a]) => Point2f(x[a], V[a]) for a in eachindex(x)]),
            V);
        color = _SEGMENT_COLOR)
    scatter!(ax, x, mk.v; markersize = _weight_markersize(mk.w, s.wref, length(mk)),
        color = :transparent, strokecolor = _OUTLINE_COLOR, strokewidth = 1)
    scatter!(ax, x, V; color = W, colorrange = s.colorrange, colorscale = s.colorscale,
        colormap = _WEIGHT_CMAP,
        markersize = lift(W -> _weight_markersize(W, s.wref, length(mk)), W))
end

function _RM.plot_eight_marker_example(data = _RM.eight_marker_homotopy(); size = (
        880, 720))
    mk, grid = data.example.markers, data.example.grid
    info = data.result.info
    V, W = data.result.markers.v, data.result.markers.w
    s = _weight_scale([mk.w, W])
    vlims = (0.8minimum(vcat(mk.v, V)), 1.15maximum(vcat(mk.v, V)))

    f = Figure(; size)
    ax = Axis(f[1, 1:2]; xlabel = L"x", ylabel = L"v",
        title = "(a) markers before (open) and after (filled)",
        limits = (-0.05grid.L, 1.05grid.L, vlims...))
    sc = _rm_eight_markers!(ax, mk, grid, Observable(V), Observable(W),
        mk.v .* exp.(info.μstar), s)
    Colorbar(f[1, 3], sc; label = "weight w")

    ax2 = Axis(f[2, 1]; xlabel = "Newton iteration", ylabel = "max |c| / scale",
        yscale = log10, title = "(b) Newton per homotopy step",
        limits = (nothing, (1e-18, 1)))
    for (j, h) in enumerate(filter(h -> h.accepted && h.τ > 0, info.history))
        its = [l.it for l in h.iterations][2:end]
        res = max.([l.residual for l in h.iterations][2:end], 1e-17)
        scatterlines!(ax2, its, res; color = _SERIES[mod1(j, 7)], markersize = 7,
            label = "τ = $(round(h.τ; digits = 2))")
    end
    axislegend(ax2; position = :rt)

    ax3 = Axis(f[2, 2:3]; xlabel = L"\tau", ylabel = L"\min_a H_a/w_a",
        title = "(c) along the homotopy", limits = (-0.02, 1.04, -0.15, 1.05))
    hlines!(ax3, 0; color = _OUTLINE_COLOR, linewidth = 0.5)
    lines!(ax3, data.path.τ, data.path.min_curvature; color = _SERIES[3])
    minHd = data.direct.iterations[end].min_curvature
    scatter!(ax3, [1.0], [minHd]; color = _SERIES[6], marker = :xcross, markersize = 12)
    text!(ax3, 0.95, 0.05; text = "no homotopy:\nmin H/w = $(round(minHd; digits = 3))",
        align = (:right, :bottom))
    f
end

function _RM.animate_eight_marker_homotopy(path, data = _RM.eight_marker_homotopy();
        framerate = 20)
    mk, grid = data.example.markers, data.example.grid
    p = momentum(mk)
    μstar = data.result.info.μstar
    τs, μs, minH = data.path.τ, data.path.μ, data.path.min_curvature
    Wall = [mk.w .* exp.(-μ) for μ in μs]
    s = _weight_scale(vcat([mk.w], Wall))
    vlims = (0.8minimum(mk.v), 1.15maximum(p ./ (mk.m .* minimum(Wall))))
    k = Observable(1)
    W = lift(i -> Wall[i], k)
    V = lift(W -> p ./ (mk.m .* W), W)
    fa = Figure(; size = (720, 480))
    axa = Axis(fa[1, 1]; xlabel = L"x", ylabel = L"v",
        limits = (-0.05grid.L, 1.05grid.L, vlims...),
        title = lift(i -> "gauge rezoning, K = {0, 2}, τ = $(round(τs[i]; digits = 3))", k))
    sc = _rm_eight_markers!(
        axa, mk, grid, V, W, lift(i -> mk.v .* exp.(τs[i] .* μstar), k),
        s)
    Colorbar(fa[1, 2], sc; label = "weight w")
    axb = Axis(fa[2, 1]; xlabel = L"\tau", ylabel = L"\min_a H_a/w_a",
        limits = (0, 1, 0, 1.05))
    lines!(axb, τs, minH; color = (_SERIES[3], 0.4))
    scatter!(axb, lift(i -> [τs[i]], k), lift(i -> [minH[i]], k); color = _SERIES[3])
    rowsize!(fa.layout, 2, Relative(0.3))
    record(fa, path, vcat(fill(1, 8), 1:length(τs), fill(length(τs), 25)); framerate) do i
        k[] = i
    end
    path
end

# ------------------------------------------------------- the numerical example of the gauge step

function _RM.plot_gauge_experiment(mk::Markers, result; size = (880, 760))
    grid = result.info.problem.grid
    after = result.markers
    λ = exp.(result.info.μ)
    gs = [idx for (_, idx) in _RM.velocity_groups(grid, mk.x, mk.v) if length(idx) ≥ 3]
    s = _weight_scale([mk.w, after.w])
    vl = 1.15maximum(abs, vcat(mk.v, after.v))
    f = Figure(; size)
    ax1 = Axis(f[1, 1]; xlabel = L"x", ylabel = L"v", title = "(a) before",
        limits = (0, grid.L, -vl, vl))
    ParticleRemapping.plot_markers!(ax1, mk; s.colorrange, s.colorscale, s.wref)
    ax2 = Axis(f[1, 2]; xlabel = L"x",
        title = "(b) after, τ = $(round(result.info.τ; digits = 3))",
        limits = (0, grid.L, -vl, vl))
    sc = ParticleRemapping.plot_markers!(ax2, after; s.colorrange, s.colorscale, s.wref)
    hideydecorations!(ax2; ticks = false)
    Colorbar(f[1, 3], sc; label = "weight w")
    ax3 = Axis(f[2, 1]; xlabel = L"\lambda_a", ylabel = "markers", xscale = log10,
        title = "(c) stretches")
    lo, hi = log10.(extrema(λ))
    hist!(ax3, λ; bins = 10 .^ range(lo - 0.05, hi + 0.05; length = 40), color = _SERIES[1])
    ax4 = Axis(f[2, 2:3]; xlabel = "largest/mean gap before", ylabel = "after",
        xscale = log10, yscale = log10, title = "(d) per cell and sign")
    before = [_RM.speed_gap_ratio(mk.v[g]) for g in gs]
    gaps = [_RM.speed_gap_ratio(after.v[g]) for g in gs]
    scatter!(ax4, before, gaps; markersize = 6, color = _SERIES[2])
    lines!(ax4, [1, maximum(before)], [1, maximum(before)]; color = _OUTLINE_COLOR,
        linestyle = :dash)
    f
end

# -------------------------------------------------------------- the cotangent lift and matching

function _RM.plot_lift_example(setup, experiment; t = setup.t, cell = 5, size = (880, 700))
    α, grid, target = setup.α, setup.grid, setup.target
    J(x) = _RM.cosine_profile(x, α)
    f = Figure(; size)
    ax1 = Axis(f[1, 1]; xlabel = L"x", ylabel = L"t", title = "(a) characteristics",
        limits = (0, 2π, 0, t), xticks = _RM_PITICKS)
    ts = range(0, t; length = 31)
    for x0 in range(0, 2π; length = 21)[1:(end - 1)]
        X = [_RM.admissible_lift(x0, s; α)[1] for s in ts]
        lines!(ax1, X, ts; color = _SERIES[1], linewidth = 0.8)
        lines!(ax1, X .- 2π, ts; color = _SERIES[1], linewidth = 0.8)
    end

    Xs = range(0, 2π; length = 200)
    xback = [_RM.lift_theta_inverse(_RM.lift_theta(X; α) - t; α) for X in Xs]
    ax2 = Axis(f[1, 2]; xlabel = L"X", title = "(b) densities at t = $t",
        limits = (0, 2π, 0, 2.3), xticks = _RM_PITICKS)
    lines!(ax2, Xs, J.(Xs); color = _SERIES[2], label = L"represented $J(X)$")
    lines!(ax2, Xs, sqrt.(J.(xback) .* J.(Xs)); color = _SERIES[1], linestyle = :dash,
        label = L"markers $\sqrt{J(x)J(X)}$")
    lines!(ax2, Xs, sqrt.(J.(Xs) ./ J.(xback)); color = _SERIES[3], linestyle = :dot,
        label = L"W/w = 1/\varphi'")
    adm = experiment.admissible.markers
    cells = cellof.(Ref(grid), adm.x)
    xc = ((1:(grid.Ng)) .- 0.5) .* grid.Δx
    cnt = [count(==(c), cells) for c in 1:(grid.Ng)]
    chg = [sum(adm.w[cells .== c]) for c in 1:(grid.Ng)]
    scatter!(ax2, xc, cnt ./ (length(adm) * grid.Δx / grid.L); color = _SERIES[1],
        markersize = 8)
    scatter!(ax2, xc, chg ./ grid.Δx; color = _SERIES[2], markersize = 8,
        marker = :utriangle)
    axislegend(ax2; position = :ct)

    ax3 = Axis(f[2, 1]; xlabel = "cell", ylabel = L"R_c^s", title = "(c) compatibility residual")
    Ra, Rb = experiment.admissible.info.R, experiment.inadmissible.info.R
    for (j, s) in enumerate(("+", "−"))
        mk = j == 1 ? :circle : :utriangle
        scatterlines!(ax3, 1:(grid.Ng), Ra[:, j]; color = _SERIES[1], marker = mk,
            markersize = 7, label = "admissible, s = $s")
        scatterlines!(ax3, 1:(grid.Ng), Rb[:, j]; color = _SERIES[6], marker = mk,
            markersize = 7, label = "inadmissible, s = $s")
    end
    axislegend(ax3; position = :rt)

    vb = range(-4, 4; length = 25)
    vmid = (vb[1:(end - 1)] .+ vb[2:end]) ./ 2
    marginal(mk) = _RM.binned_weights(grid, mk.x, mk.v, mk.w, vb)[cell, :] ./ step(vb)
    ax4 = Axis(f[2, 2]; xlabel = L"V", title = "(d) velocity marginal, cell $cell",
        limits = (-4, 4, 0, nothing))
    vf = range(-4, 4; step = 0.02)
    lines!(ax4, vf, target.Jbar[cell] * grid.Δx .* exp.(-vf .^ 2 ./ 2) ./ sqrt(2π);
        color = _OUTLINE_COLOR, label = "exact")
    stairs!(ax4, vmid, marginal(adm); step = :center, color = _SERIES[1], label = "matched")
    stairs!(ax4, vmid, marginal(experiment.kept); step = :center, color = _SERIES[3],
        label = "W = w kept")
    stairs!(ax4, vmid, marginal(experiment.inadmissible.markers); step = :center,
        color = _SERIES[6], label = "inadmissible")
    axislegend(ax4; position = :lt)
    f
end

function _RM.animate_lift_matching(path, frames, setup; framerate = 20)
    grid, α = setup.grid, setup.α
    N = length(setup.markers)
    n = length(frames[1].X)
    s = _weight_scale([fr.W for fr in frames])
    Xs = range(0, grid.L; length = 200)
    xc = ((1:(grid.Ng)) .- 0.5) .* grid.Δx
    k = Observable(1)
    fr = lift(i -> frames[i], k)
    fa = Figure(; size = (760, 520))
    axa = Axis(fa[1, 1]; xlabel = L"X", ylabel = L"V", limits = (0, grid.L, -3.5, 3.5),
        title = lift(
            s -> "cotangent lift of a = 1/√J with matching, t = " *
                 string(round(s.t; digits = 2)),
            fr))
    sca = scatter!(axa, lift(s -> s.X, fr), lift(s -> s.V, fr);
        color = lift(s -> s.W, fr), colorrange = s.colorrange, colorscale = s.colorscale,
        colormap = _WEIGHT_CMAP, markersize = lift(f -> _weight_markersize(f.W, s.wref, n), fr))
    Colorbar(fa[1, 2], sca; label = "weight W")
    axb = Axis(fa[2, 1]; xlabel = L"X", ylabel = "per cell", limits = (0, grid.L, 0, 1.7))
    lines!(axb, Xs, _RM.cosine_profile.(Xs, α); color = _SERIES[2])
    scatter!(axb, xc, lift(s -> s.counts ./ (N * grid.Δx / grid.L), fr);
        color = _SERIES[1], markersize = 8, label = "markers")
    scatter!(axb, xc, lift(s -> s.charges ./ grid.Δx, fr); color = _SERIES[2],
        marker = :utriangle, markersize = 8, label = L"\sum W")
    axislegend(axb; position = :rb, orientation = :horizontal)
    rowsize!(fa.layout, 2, Relative(0.3))
    record(fa, path, vcat(fill(1, 8), 1:length(frames), fill(length(frames), 20));
        framerate) do i
        k[] = i
    end
    path
end

# --------------------------------------------------------------------- the five-marker matching

function _rm_draw_matching!(ax, ex, r)
    C = _category_colors(length(ex.v))
    wref = max(ex.w0, maximum(ex.W))
    n = length(ex.v)
    Vs = range(0, 1.75maximum(ex.V); length = 300)
    lines!(ax, Vs, ex.Jbar .* _RM.maxwellian_current.(Vs); color = _OUTLINE_COLOR)
    hlines!(ax, ex.Jbar * _RM.maxwellian_current(Inf); color = _OUTLINE_COLOR,
        linestyle = :dot)
    for j in 1:r
        linesegments!(ax,
            [Point2f(0, ex.Πc[j]) => Point2f(ex.V[j], ex.Πc[j]),
                Point2f(ex.V[j], ex.Πc[j]) => Point2f(ex.V[j], 0)];
            color = C[j], linestyle = :dash)
        scatter!(ax, [ex.V[j]], [ex.Πc[j]]; color = C[j],
            markersize = _weight_markersize([ex.W[j]], wref, n))
        scatter!(ax, [ex.v[j]], [ex.Jbar * _RM.maxwellian_current(ex.v[j])];
            color = :transparent, strokecolor = C[j], strokewidth = 1.2,
            markersize = _weight_markersize([ex.w0], wref, n))
    end
end

function _RM.plot_matching_example(ex = _RM.quantile_matching_example(); size = (880, 400))
    Vmax = 1.75maximum(ex.V)
    Fmax = 1.05ex.Jbar * _RM.maxwellian_current(Inf)
    f = Figure(; size)
    ax1 = Axis(f[1, 1]; xlabel = L"V", ylabel = L"\hat{F}(V)",
        title = "(a) quantile matching", limits = (0, Vmax, 0, Fmax))
    _rm_draw_matching!(ax1, ex, length(ex.v))
    ax2 = Axis(f[1, 2]; xlabel = L"V", ylabel = "weight", title = "(b) weights",
        limits = (0, Vmax, 0, 1.2max(ex.w0, maximum(ex.W))))
    stem!(ax2, ex.v, fill(ex.w0, length(ex.v)); color = :transparent,
        strokecolor = _OUTLINE_COLOR, strokewidth = 1, stemcolor = _SEGMENT_COLOR,
        trunkcolor = :transparent, markersize = 9)
    stem!(ax2, ex.V, ex.W; color = _category_colors(length(ex.v)),
        stemcolor = _OUTLINE_COLOR, trunkcolor = :transparent, markersize = 10)
    f
end

function _RM.animate_matching_example(path, ex = _RM.quantile_matching_example();
        framerate = 20)
    n = length(ex.v)
    fa = Figure(; size = (720, 400))
    axa = Axis(fa[1, 1]; xlabel = L"V", ylabel = L"\hat{F}(V)",
        limits = (0, 1.75maximum(ex.V), 0, 1.05ex.Jbar * _RM.maxwellian_current(Inf)))
    record(fa, path, vcat(fill(0, 10), repeat(1:n; inner = 24), fill(n, 20)); framerate) do r
        empty!(axa)
        _rm_draw_matching!(axa, ex, r)
        axa.title = r == 0 ? "cumulative current of f" :
                    "marker $r: Π = $(round(ex.Πc[r]; digits = 4)) → V = " *
                    "$(round(ex.V[r]; digits = 3)), W = $(round(ex.W[r]; digits = 4))"
    end
    path
end

# ------------------------------------------------------------------- pictures of the algorithms

function _RM.plot_kkt_structure(prob, μ, α; size = (880, 420))
    H = _RM.gauge_hessian(prob, μ, α)
    _, G = _RM.gauge_constraints(prob, μ)
    K = Matrix(_RM.gauge_kkt_matrix(H, G))
    S = Matrix(_RM.gauge_schur_matrix(H, G))
    pattern(A) = [iszero(a) ? NaN : log10(abs(a)) for a in A]
    PK, PS = pattern(K), pattern(S)
    cr = extrema(filter(!isnan, vcat(vec(PK), vec(PS))))
    N, nc = length(H), length(S[:, 1])
    f = Figure(; size)
    ax1 = Axis(f[1, 1]; title = "(a) KKT matrix, $(N + nc) × $(N + nc)", yreversed = true,
        aspect = DataAspect(), xlabel = "column", ylabel = "row")
    hm = heatmap!(ax1, 1:(N + nc), 1:(N + nc), permutedims(PK); colormap = _FIELD_CMAP,
        colorrange = cr)
    hlines!(ax1, N + 0.5; color = _OUTLINE_COLOR, linewidth = 0.8)
    vlines!(ax1, N + 0.5; color = _OUTLINE_COLOR, linewidth = 0.8)
    ax2 = Axis(f[1, 2]; title = "(b) Schur complement, $nc × $nc", yreversed = true,
        aspect = DataAspect(), xlabel = "(k, i)", ylabel = "(k, i)")
    heatmap!(ax2, 1:nc, 1:nc, permutedims(PS); colormap = _FIELD_CMAP, colorrange = cr)
    Colorbar(f[1, 3], hm; label = L"\log_{10} |\text{entry}|")
    f
end

function _RM.plot_gauge_orbits(; m = 1, momenta = (0.25, 0.5, 1, 2), λs = (0.5, 1, 2),
        size = (880, 400))
    v0, w0 = 1.0, 1.0 / m
    ws = [w0 / λ for λ in λs]
    s = _weight_scale([ws])
    f = Figure(; size)
    ax1 = Axis(f[1, 1]; xlabel = L"v", ylabel = L"w", title = "(a) orbits w v = p/m",
        limits = (0, 4, 0, 4))
    vs = range(0.02, 4; length = 400)
    for (j, p) in enumerate(momenta)
        lines!(ax1, vs, p ./ (m .* vs); color = _SERIES[j], label = "p = $p")
    end
    sc = scatter!(ax1, [λ * v0 for λ in λs], ws; color = ws, colorrange = s.colorrange,
        colorscale = s.colorscale, colormap = _WEIGHT_CMAP,
        markersize = _weight_markersize(ws, s.wref, length(ws)))
    for λ in λs
        text!(ax1, λ * v0, w0 / λ; text = "λ = $λ", offset = (8, 6))
    end
    axislegend(ax1; position = :rt)
    ax2 = Axis(f[1, 2]; xlabel = L"x", ylabel = L"v",
        title = "(b) area × λ, weight × 1/λ",
        limits = (0, 1, 0, 2.6))
    δx, δv, vb = 0.5, 0.3, 0.5
    for (λ, c) in zip(λs, _rm_weight_colors(ws, s))
        poly!(ax2, Rect2f(0.25, λ * vb, δx, λ * δv); color = (c, 0.6), strokecolor = c,
            strokewidth = 1)
        text!(ax2, 0.25 + δx + 0.02, λ * (vb + δv / 2); text = "λ = $λ: w/λ",
            align = (:left, :center))
    end
    Colorbar(f[1, 3], sc; label = "weight w")
    f
end

function _RM.plot_velocity_momentum_sheet(; vth = 1, m = 1, size = (880, 400))
    vs = range(-3vth, 3vth; length = 400)
    wM(v) = exp(-v^2 / (2vth^2))
    f = Figure(; size)
    ax1 = Axis(f[1, 1]; xlabel = L"v", ylabel = L"p = m\, w(v)\, v",
        title = "(a) the sheet map at fixed x")
    lines!(ax1, vs, m .* 0.6 .* vs ./ vth; color = _SERIES[1], label = "equal weights")
    lines!(ax1, vs, m .* wM.(vs) .* vs; color = _SERIES[2],
        label = L"w \propto \text{Maxwellian}")
    vlines!(ax1, [-vth, vth]; color = _OUTLINE_COLOR, linestyle = :dash)
    axislegend(ax1; position = :lt)
    ax2 = Axis(f[1, 2]; xlabel = L"v", ylabel = L"\partial_v (w v)",
        title = "(b) the single-sheet condition")
    lines!(ax2, vs, fill(0.6 / vth, length(vs)); color = _SERIES[1])
    lines!(ax2, vs, (1 .- vs .^ 2 ./ vth^2) .* wM.(vs); color = _SERIES[2])
    hlines!(ax2, 0; color = _OUTLINE_COLOR, linewidth = 0.5)
    vlines!(ax2, [-vth, vth]; color = _OUTLINE_COLOR, linestyle = :dash)
    f
end
