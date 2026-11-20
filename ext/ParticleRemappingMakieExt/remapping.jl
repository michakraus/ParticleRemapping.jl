# Plots of `ParticleRemapping.Remapping`.

const _RM = ParticleRemapping.Remapping

_rm_colors(n) = [Makie.wong_colors()[mod1(a, 7)] for a in 1:n]
_rm_msize(w, wref; s = 20) = s .* sqrt.(w ./ wref)                # marker area ∝ weight
const _RM_PITICKS = (0:(π / 2):(2π), ["0", "π/2", "π", "3π/2", "2π"])

# ------------------------------------------------------------------ the charge-constraint example

function _RM.plot_charge_example(ex = _RM.charge_constraint_example(); size = (880, 660))
    βf = _RM.charge_fold(ex)
    e1 = _RM.charge_stationary_point(ex, 1.0)
    v, w, μs = ex.v, ex.w, ex.μstar
    n, wref = length(v), maximum(w)
    C, ca = Makie.wong_colors(), _rm_colors(length(v))
    vmax = 1.1maximum(vcat(v, ex.Vstar, e1.V))

    f = Figure(; size)
    ax = Axis(f[1, 1:2]; xlabel = L"speed $|v|$", title = "(a) speeds and weights",
        yticks = (1:3, ["result, β = 1", "targets", "before"]), limits = (
            0, vmax, 0.4, 3.6))
    linesegments!(ax, vec([Point2f(v[a], 3) => Point2f(ex.Vstar[a], 2) for a in 1:n]);
        color = :gray70, linestyle = :dot)
    linesegments!(ax, vec([Point2f(v[a], 3) => Point2f(e1.V[a], 1) for a in 1:n]);
        color = :gray40)
    scatter!(ax, v, fill(3, n); markersize = _rm_msize(w, wref), color = ca)
    scatter!(ax, ex.Vstar, fill(2, n); markersize = _rm_msize(w .* exp.(-μs), wref),
        color = :white, strokecolor = ca, strokewidth = 1.5)
    scatter!(ax, e1.V, fill(1, n); markersize = _rm_msize(e1.W, wref), color = ca)

    ax2 = Axis(f[2, 1]; xlabel = L"multiplier $A$", ylabel = L"c^{(0)}",
        title = "(b) constraint")
    for (j, β) in enumerate((1.0, βf, 1.6))
        As = range(_RM.charge_multiplier_limit(ex, β), 0.3; length = 300)
        lines!(ax2, As, [_RM.charge_constraint_residual(ex, A, β) for A in As];
            color = C[4 + j],
            label = β == βf ? "β = $(round(β; digits = 3)) (fold)" :
                    "β = $(round(β; digits = 2))")
    end
    hlines!(ax2, 0; color = :black, linewidth = 0.5)
    vlines!(ax2, -exp(-1) / ex.s2; color = :gray50, linestyle = :dash)
    scatter!(ax2, [e1.A], [0.0]; color = C[5], markersize = 9)
    axislegend(ax2; position = :rt, framevisible = false)

    βs = range(0, βf; length = 200)
    sols = [_RM.charge_stationary_point(ex, β) for β in βs]
    ax3 = Axis(f[2, 2]; xlabel = L"target scale $\beta$", ylabel = L"H_a / w_a",
        title = "(c) curvature")
    for a in 1:n
        lines!(ax3, βs, [s.Hw[a] for s in sols]; color = ca[a], label = "a = $a")
    end
    vlines!(ax3, βf; color = :gray50, linestyle = :dash)
    ylims!(ax3, 0, 1.05 / ex.s2)
    axislegend(ax3; position = :lb, framevisible = false)
    f
end

function _RM.animate_charge_fold(path, ex = _RM.charge_constraint_example(); framerate = 30)
    βf = _RM.charge_fold(ex)
    v, w, μs = ex.v, ex.w, ex.μstar
    n, wref, ca = length(v), maximum(w), _rm_colors(length(v))
    βs = range(0, βf; length = 200)
    sols = [_RM.charge_stationary_point(ex, β) for β in βs]
    βa = Observable(0.0)
    sol = lift(β -> _RM.charge_stationary_point(ex, β), βa)
    fa = Figure(; size = (720, 450))
    axa = Axis(fa[1, 1]; xlabel = L"speed $|v|$",
        yticks = (1:2, ["stationary point", "targets β μ*"]),
        limits = (0, 1.4maximum(v .* exp.(βf .* μs)), 0.4, 2.6),
        title = lift(β -> "charge constraint only, β = $(round(β; digits = 3))", βa))
    scatter!(axa, lift(β -> v .* exp.(β .* μs), βa), fill(2, n);
        markersize = lift(β -> _rm_msize(w .* exp.(-β .* μs), wref; s = 26), βa),
        color = :white, strokecolor = ca, strokewidth = 1.5)
    scatter!(axa, lift(s -> s.V, sol), fill(1, n);
        markersize = lift(s -> _rm_msize(s.W, wref; s = 26), sol), color = ca)
    axb = Axis(fa[2, 1]; xlabel = L"\beta", ylabel = L"H_a/w_a",
        limits = (0, βf * 1.02, 0, 1.05 / ex.s2))
    for a in 1:n
        lines!(axb, βs, [s.Hw[a] for s in sols]; color = (ca[a], 0.35))
    end
    scatter!(
        axb, lift(β -> fill(β, n), βa), lift(s -> s.Hw, sol); color = ca, markersize = 9)
    vlines!(axb, βf; color = :gray50, linestyle = :dash)
    record(fa, path, vcat(range(0, βf; length = 120), fill(βf, 30)); framerate) do β
        βa[] = β
    end
    path
end

# ----------------------------------------------------------------------- the eight-marker example

function _RM.plot_eight_marker_example(data = _RM.eight_marker_homotopy(); size = (
        880, 720))
    mk, grid = data.example.markers, data.example.grid
    info = data.result.info
    x, v, w = mk.x, mk.v, mk.w
    V, W = data.result.markers.v, data.result.markers.w
    Vs = v .* exp.(info.μstar)
    wref = maximum(W)
    C = Makie.wong_colors()
    nodes = collect(range(0, grid.L; step = grid.Δx))
    vlims = (0.8minimum(vcat(v, V)), 1.15maximum(vcat(v, V)))

    f = Figure(; size)
    ax = Axis(f[1, 1:2]; xlabel = L"x", ylabel = L"v",
        title = "(a) markers before (open) and after (filled)",
        limits = (-0.05grid.L, 1.05grid.L, vlims...))
    vlines!(ax, nodes; color = :gray80, linestyle = :dash)
    scatter!(ax, x, Vs; marker = :hline, markersize = 18, color = :gray50)
    linesegments!(
        ax, vec([Point2f(x[a], v[a]) => Point2f(x[a], V[a]) for a in eachindex(x)]);
        color = :gray40)
    scatter!(ax, x, v; markersize = _rm_msize(w, wref), color = :white, strokecolor = C[1],
        strokewidth = 1.2)
    scatter!(ax, x, V; markersize = _rm_msize(W, wref), color = C[2])

    ax2 = Axis(f[2, 1]; xlabel = "Newton iteration", ylabel = "max |c| / scale",
        yscale = log10, title = "(b) Newton per homotopy step", limits = (
            nothing, (1e-18, 1)))
    for (j, h) in enumerate(filter(h -> h.accepted && h.τ > 0, info.history))
        its = [l.it for l in h.iterations][2:end]
        res = max.([l.residual for l in h.iterations][2:end], 1e-17)
        scatterlines!(ax2, its, res; color = C[mod1(j, 7)], markersize = 7,
            label = "τ = $(round(h.τ; digits = 2))")
    end
    axislegend(ax2; position = :rt, framevisible = false)

    ax3 = Axis(f[2, 2]; xlabel = L"\tau", ylabel = L"\min_a H_a/w_a",
        title = "(c) along the homotopy", limits = (-0.02, 1.04, -0.15, 1.05))
    hlines!(ax3, 0; color = :black, linewidth = 0.5)
    lines!(ax3, data.path.τ, data.path.min_curvature; color = C[3])
    minHd = data.direct.iterations[end].min_curvature
    scatter!(ax3, [1.0], [minHd]; color = :red, marker = :xcross, markersize = 12)
    text!(ax3, 0.95, 0.05; text = "no homotopy:\nmin H/w = $(round(minHd; digits = 3))",
        align = (:right, :bottom))
    f
end

function _RM.animate_eight_marker_homotopy(path, data = _RM.eight_marker_homotopy();
        framerate = 25)
    mk, grid = data.example.markers, data.example.grid
    x, v, w = mk.x, mk.v, mk.w
    p = mk.m .* w .* v
    μstar = data.result.info.μstar
    τs, μs, minH = data.path.τ, data.path.μ, data.path.min_curvature
    Wall = [w .* exp.(-μ) for μ in μs]
    wref = maximum(maximum, Wall)
    C = Makie.wong_colors()
    vlims = (0.8minimum(v), 1.15maximum(p ./ (mk.m .* minimum(Wall))))
    k = Observable(1)
    fa = Figure(; size = (720, 480))
    axa = Axis(fa[1, 1]; xlabel = L"x", ylabel = L"v",
        limits = (-0.05grid.L, 1.05grid.L, vlims...),
        title = lift(i -> "gauge rezoning, K = {0, 2}, τ = $(round(τs[i]; digits = 3))", k))
    vlines!(axa, collect(range(0, grid.L; step = grid.Δx)); color = :gray80, linestyle = :dash)
    scatter!(axa, x, lift(i -> v .* exp.(τs[i] .* μstar), k); marker = :hline,
        markersize = 18, color = :gray50)
    scatter!(axa, x, v; markersize = _rm_msize(w, wref; s = 24), color = :white,
        strokecolor = C[1], strokewidth = 1)
    scatter!(axa, x, lift(i -> p ./ (mk.m .* Wall[i]), k);
        markersize = lift(i -> _rm_msize(Wall[i], wref; s = 24), k), color = C[2])
    axb = Axis(fa[2, 1]; xlabel = L"\tau", ylabel = L"\min_a H_a/w_a",
        limits = (0, 1, 0, 1.05))
    lines!(axb, τs, minH; color = (C[3], 0.4))
    scatter!(axb, lift(i -> [τs[i]], k), lift(i -> [minH[i]], k); color = C[3])
    rowsize!(fa.layout, 2, Relative(0.3))
    record(fa, path, vcat(1:length(τs), fill(length(τs), 25)); framerate) do i
        k[] = i
    end
    path
end

# ------------------------------------------------------- the numerical example of the gauge step

function _RM.plot_gauge_experiment(mk::Markers, result; size = (880, 760))
    grid = result.info.problem.grid
    x, v, w = mk.x, mk.v, mk.w
    V, W = result.markers.v, result.markers.w
    λ = exp.(result.info.μ)
    C = Makie.wong_colors()
    gs = [idx for (_, idx) in _RM.velocity_groups(grid, x, v) if length(idx) ≥ 3]
    crange = extrema(log10.(w))
    vl = 1.15maximum(abs, vcat(v, V))
    f = Figure(; size)
    ax1 = Axis(f[1, 1]; xlabel = L"x", ylabel = L"v", title = "(a) before",
        limits = (0, grid.L, -vl, vl))
    sc = scatter!(ax1, x, v; markersize = 2, color = log10.(w), colorrange = crange,
        colormap = :viridis)
    ax2 = Axis(f[1, 2]; xlabel = L"x",
        title = "(b) after, τ = $(round(result.info.τ; digits = 3))",
        limits = (0, grid.L, -vl, vl))
    scatter!(ax2, x, V; markersize = 2, color = log10.(W), colorrange = crange,
        colormap = :viridis)
    hideydecorations!(ax2; ticks = false)
    Colorbar(f[1, 3], sc; label = L"\log_{10} w")
    ax3 = Axis(f[2, 1]; xlabel = L"\lambda_a", ylabel = "markers", xscale = log10,
        title = "(c) stretches")
    lo, hi = log10.(extrema(λ))
    hist!(ax3, λ; bins = 10 .^ range(lo - 0.05, hi + 0.05; length = 40), color = C[1])
    ax4 = Axis(f[2, 2:3]; xlabel = "largest/mean gap before", ylabel = "after",
        xscale = log10, yscale = log10, title = "(d) per cell and sign")
    before = [_RM.speed_gap_ratio(v[g]) for g in gs]
    after = [_RM.speed_gap_ratio(V[g]) for g in gs]
    scatter!(ax4, before, after; markersize = 6, color = C[2])
    lines!(
        ax4, [1, maximum(before)], [1, maximum(before)]; color = :gray60, linestyle = :dash)
    f
end

function _RM.animate_gauge_homotopy(path, mk::Markers, result; framerate = 15)
    grid = result.info.problem.grid
    x, v, w = mk.x, mk.v, mk.w
    p = mk.m .* w .* v
    steps = [h for h in result.info.history if h.accepted]
    crange = extrema(log10.(w))
    vl = 1.15maximum(abs, v .* exp.(maximum(abs, steps[end].μ)))
    vl = min(vl, 2maximum(abs, v))
    sel = findall(<(2grid.Δx), x)
    k = Observable(1)
    Vk = lift(i -> p ./ (mk.m .* w .* exp.(-steps[i].μ)), k)
    fa = Figure(; size = (760, 480))
    axa = Axis(fa[1, 1]; xlabel = L"x", ylabel = L"v", limits = (0, 2grid.Δx, -vl, vl),
        title = lift(i -> "cells 1 and 2, τ = $(round(steps[i].τ; digits = 3))", k))
    scatter!(axa, x[sel], lift(V -> V[sel], Vk);
        markersize = lift(
            i -> 3 .+ 8 .* sqrt.(w[sel] .* exp.(-steps[i].μ[sel]) ./ maximum(w)), k),
        color = log10.(w[sel]), colorrange = crange, colormap = :viridis)
    axb = Axis(fa[1, 2]; xlabel = L"x", ylabel = L"v", limits = (0, grid.L, -vl, vl),
        title = "all cells")
    scatter!(axb, x, Vk; markersize = 1.5, color = log10.(w), colorrange = crange,
        colormap = :viridis)
    record(fa, path, vcat(1:length(steps), fill(length(steps), 20)); framerate) do i
        k[] = i
    end
    path
end

# -------------------------------------------------------------- the cotangent lift and matching

function _RM.plot_lift_example(setup, experiment; t = setup.t, cell = 5, size = (880, 700))
    α, grid, target = setup.α, setup.grid, setup.target
    J(x) = _RM.cosine_profile(x, α)
    C = Makie.wong_colors()
    f = Figure(; size)
    ax1 = Axis(f[1, 1]; xlabel = L"x", ylabel = L"t", title = "(a) characteristics",
        limits = (0, 2π, 0, t), xticks = _RM_PITICKS)
    ts = range(0, t; length = 31)
    for x0 in range(0, 2π; length = 21)[1:(end - 1)]
        X = [_RM.admissible_lift(x0, s; α)[1] for s in ts]
        lines!(ax1, X, ts; color = C[1], linewidth = 0.8)
        lines!(ax1, X .- 2π, ts; color = C[1], linewidth = 0.8)
    end

    Xs = range(0, 2π; length = 200)
    xback = [_RM.lift_theta_inverse(_RM.lift_theta(X; α) - t; α) for X in Xs]
    ax2 = Axis(f[1, 2]; xlabel = L"X", title = "(b) densities at t = $t",
        limits = (0, 2π, 0, 2.3), xticks = _RM_PITICKS)
    lines!(ax2, Xs, J.(Xs); color = C[2], label = L"represented $J(X)$")
    lines!(ax2, Xs, sqrt.(J.(xback) .* J.(Xs)); color = C[1], linestyle = :dash,
        label = L"markers $\sqrt{J(x)J(X)}$")
    lines!(ax2, Xs, sqrt.(J.(Xs) ./ J.(xback)); color = C[3], linestyle = :dot,
        label = L"W/w = 1/\varphi'")
    adm = experiment.admissible.markers
    cells = cellof.(Ref(grid), adm.x)
    xc = ((1:(grid.Ng)) .- 0.5) .* grid.Δx
    cnt = [count(==(c), cells) for c in 1:(grid.Ng)]
    chg = [sum(adm.w[cells .== c]) for c in 1:(grid.Ng)]
    scatter!(ax2, xc, cnt ./ (length(adm) * grid.Δx / grid.L); color = C[1], markersize = 8)
    scatter!(ax2, xc, chg ./ grid.Δx; color = C[2], markersize = 8, marker = :utriangle)
    axislegend(ax2; position = :ct, framevisible = false)

    ax3 = Axis(f[2, 1]; xlabel = "cell", ylabel = L"R_c^s", title = "(c) compatibility residual")
    Ra, Rb = experiment.admissible.info.R, experiment.inadmissible.info.R
    for (j, s) in enumerate(("+", "−"))
        mk = j == 1 ? :circle : :utriangle
        scatterlines!(
            ax3, 1:(grid.Ng), Ra[:, j]; color = C[1], marker = mk, markersize = 7,
            label = "admissible, s = $s")
        scatterlines!(
            ax3, 1:(grid.Ng), Rb[:, j]; color = C[6], marker = mk, markersize = 7,
            label = "inadmissible, s = $s")
    end
    axislegend(ax3; position = :rt, framevisible = false)

    vb = range(-4, 4; length = 25)
    vmid = (vb[1:(end - 1)] .+ vb[2:end]) ./ 2
    marginal(mk) = _RM.binned_weights(grid, mk.x, mk.v, mk.w, vb)[cell, :] ./ step(vb)
    ax4 = Axis(f[2, 2]; xlabel = L"V", title = "(d) velocity marginal, cell $cell",
        limits = (-4, 4, 0, nothing))
    vf = range(-4, 4; step = 0.02)
    lines!(ax4, vf, target.Jbar[cell] * grid.Δx .* exp.(-vf .^ 2 ./ 2) ./ sqrt(2π);
        color = :black, label = "exact")
    stairs!(ax4, vmid, marginal(adm); step = :center, color = C[1], label = "matched")
    stairs!(ax4, vmid, marginal(experiment.kept); step = :center, color = C[3],
        label = "W = w kept")
    stairs!(ax4, vmid, marginal(experiment.inadmissible.markers); step = :center,
        color = C[6], label = "inadmissible")
    axislegend(ax4; position = :lt, framevisible = false)
    f
end

function _RM.animate_lift_matching(path, frames, setup; framerate = 20)
    grid, α = setup.grid, setup.α
    N, w0 = length(setup.markers), setup.markers.w[1]
    C = Makie.wong_colors()
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
    sca = scatter!(axa, lift(s -> s.X, fr), lift(s -> s.V, fr); markersize = 2.5,
        color = lift(s -> s.W ./ w0, fr), colorrange = (0.6, 1.6), colormap = :coolwarm)
    Colorbar(fa[1, 2], sca; label = L"W/w")
    axb = Axis(fa[2, 1]; xlabel = L"X", ylabel = "per cell", limits = (0, grid.L, 0, 1.7))
    lines!(axb, Xs, _RM.cosine_profile.(Xs, α); color = C[2])
    scatter!(axb, xc, lift(s -> s.counts ./ (N * grid.Δx / grid.L), fr); color = C[1],
        markersize = 8, label = "markers")
    scatter!(
        axb, xc, lift(s -> s.charges ./ grid.Δx, fr); color = C[2], marker = :utriangle,
        markersize = 8, label = L"\sum W")
    axislegend(axb; position = :rb, framevisible = false, orientation = :horizontal)
    rowsize!(fa.layout, 2, Relative(0.3))
    record(fa, path, vcat(1:length(frames), fill(length(frames), 20)); framerate) do i
        k[] = i
    end
    path
end

# --------------------------------------------------------------------- the five-marker matching

function _rm_draw_matching!(ax, ex, r)
    C = _rm_colors(length(ex.v))
    Vs = range(0, 1.75maximum(ex.V); length = 300)
    lines!(ax, Vs, ex.Jbar .* _RM.maxwellian_current.(Vs); color = :black)
    hlines!(ax, ex.Jbar * _RM.maxwellian_current(Inf); color = :gray60, linestyle = :dot)
    for j in 1:r
        linesegments!(ax,
            [Point2f(0, ex.Πc[j]) => Point2f(ex.V[j], ex.Πc[j]),
                Point2f(ex.V[j], ex.Πc[j]) => Point2f(ex.V[j], 0)];
            color = C[j], linestyle = :dash)
        scatter!(ax, [ex.V[j]], [ex.Πc[j]]; color = C[j], markersize = 10)
        scatter!(
            ax, [ex.v[j]], [ex.Jbar * _RM.maxwellian_current(ex.v[j])]; color = :white,
            strokecolor = C[j], strokewidth = 1.2, markersize = 10)
    end
end

function _RM.plot_matching_example(ex = _RM.quantile_matching_example(); size = (880, 400))
    Vmax = 1.75maximum(ex.V)
    Fmax = 1.05ex.Jbar * _RM.maxwellian_current(Inf)
    f = Figure(; size)
    ax1 = Axis(
        f[1, 1]; xlabel = L"V", ylabel = L"\hat{F}(V)", title = "(a) quantile matching",
        limits = (0, Vmax, 0, Fmax))
    _rm_draw_matching!(ax1, ex, length(ex.v))
    ax2 = Axis(f[1, 2]; xlabel = L"V", ylabel = "weight", title = "(b) weights",
        limits = (0, Vmax, 0, 1.2max(ex.w0, maximum(ex.W))))
    stem!(ax2, ex.v, fill(ex.w0, length(ex.v)); color = :gray60, stemcolor = :gray60,
        trunkcolor = :transparent, markersize = 9)
    stem!(ax2, ex.V, ex.W; color = _rm_colors(length(ex.v)), stemcolor = :black,
        trunkcolor = :transparent, markersize = 10)
    f
end

function _RM.animate_matching_example(path, ex = _RM.quantile_matching_example();
        framerate = 15)
    n = length(ex.v)
    fa = Figure(; size = (720, 400))
    axa = Axis(fa[1, 1]; xlabel = L"V", ylabel = L"\hat{F}(V)",
        limits = (0, 1.75maximum(ex.V), 0, 1.05ex.Jbar * _RM.maxwellian_current(Inf)))
    record(fa, path, vcat(fill(0, 10), repeat(1:n; inner = 18), fill(n, 20)); framerate) do r
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
    hm = heatmap!(ax1, 1:(N + nc), 1:(N + nc), permutedims(PK); colormap = :viridis,
        colorrange = cr)
    hlines!(ax1, N + 0.5; color = :red, linewidth = 0.8)
    vlines!(ax1, N + 0.5; color = :red, linewidth = 0.8)
    ax2 = Axis(f[1, 2]; title = "(b) Schur complement, $nc × $nc", yreversed = true,
        aspect = DataAspect(), xlabel = "(k, i)", ylabel = "(k, i)")
    heatmap!(ax2, 1:nc, 1:nc, permutedims(PS); colormap = :viridis, colorrange = cr)
    Colorbar(f[1, 3], hm; label = L"\log_{10} |\text{entry}|")
    f
end

function _RM.plot_gauge_orbits(; m = 1, momenta = (0.25, 0.5, 1, 2), λs = (0.5, 1, 2),
        size = (880, 400))
    C = Makie.wong_colors()
    f = Figure(; size)
    ax1 = Axis(f[1, 1]; xlabel = L"v", ylabel = L"w", title = L"(a) orbits $w v = p/m$",
        limits = (0, 4, 0, 4))
    vs = range(0.02, 4; length = 400)
    for (j, p) in enumerate(momenta)
        lines!(ax1, vs, p ./ (m .* vs); color = (C[j], 0.8), label = "p = $p")
    end
    v0, w0 = 1.0, 1.0 / m
    scatter!(
        ax1, [λ * v0 for λ in λs], [w0 / λ for λ in λs]; color = :black, markersize = 10)
    for λ in λs
        text!(ax1, λ * v0, w0 / λ; text = "λ = $λ", offset = (6, 6))
    end
    axislegend(ax1; position = :rt, framevisible = false)
    ax2 = Axis(f[1, 2]; xlabel = L"x", ylabel = L"v",
        title = L"(b) area $\times \lambda$, weight $\times 1/\lambda$", limits = (
            0, 1, 0, 2.6))
    δx, δv, vb = 0.5, 0.3, 0.5
    for (j, λ) in enumerate(λs)
        poly!(
            ax2, Rect2f(0.25, λ * vb, δx, λ * δv); color = (C[j], 0.35), strokecolor = C[j],
            strokewidth = 1)
        text!(ax2, 0.25 + δx + 0.02, λ * (vb + δv / 2); text = "λ = $λ: w/λ",
            align = (:left, :center))
    end
    f
end

function _RM.plot_velocity_momentum_sheet(; vth = 1, m = 1, size = (880, 400))
    C = Makie.wong_colors()
    vs = range(-3vth, 3vth; length = 400)
    wM(v) = exp(-v^2 / (2vth^2))
    f = Figure(; size)
    ax1 = Axis(f[1, 1]; xlabel = L"v", ylabel = L"p = m\, w(v)\, v",
        title = "(a) the sheet map at fixed x")
    lines!(ax1, vs, m .* 0.6 .* vs ./ vth; color = C[1], label = "equal weights")
    lines!(
        ax1, vs, m .* wM.(vs) .* vs; color = C[2], label = L"w \propto \text{Maxwellian}")
    vlines!(ax1, [-vth, vth]; color = :gray50, linestyle = :dash)
    axislegend(ax1; position = :lt, framevisible = false)
    ax2 = Axis(f[1, 2]; xlabel = L"v", ylabel = L"\partial_v (w v)",
        title = "(b) the single-sheet condition")
    lines!(ax2, vs, fill(0.6 / vth, length(vs)); color = C[1])
    lines!(ax2, vs, (1 .- vs .^ 2 ./ vth^2) .* wM.(vs); color = C[2])
    hlines!(ax2, 0; color = :black, linewidth = 0.5)
    vlines!(ax2, [-vth, vth]; color = :gray50, linestyle = :dash)
    f
end
