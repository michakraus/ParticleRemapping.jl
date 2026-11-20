# Plots of `ParticleRemapping.Redistribution`. The standard pictures of the remap (S1, S2, A1) are
# the common `plot_remap`, `plot_stages` and `animate_stages`.

const _RD = ParticleRemapping.Redistribution

function _RD.plot_velocity_momentum_fold(v::AbstractVector, P::AbstractVector,
        D::AbstractVector; vth::Real, pmax::Real, size = (900, 360))
    fig = Figure(; size = size)
    ax1 = Axis(fig[1, 1]; xlabel = L"\bar{v}",
        ylabel = L"\bar{p} = m \, \Delta\Omega \, \bar{v} f_h(\bar{v})",
        title = "(a) canonical momentum of a density-tied weight")
    ax2 = Axis(fig[1, 2]; xlabel = L"\bar{v}", ylabel = L"\det D_{\bar{v}} P",
        title = "(b) Jacobian determinant")
    vmax = maximum(abs, v)
    for ax in (ax1, ax2)
        vspan!(ax, [-vmax, vth], [-vth, vmax]; color = (_OUTLINE_COLOR, 0.08))
        vlines!(ax, [-vth, vth]; color = _OUTLINE_COLOR, linestyle = :dash, linewidth = 1)
    end
    lines!(ax1, v, P; color = _SERIES[1], linewidth = 2)
    hlines!(ax1, [-pmax, pmax]; color = _SERIES[2], linestyle = :dot, linewidth = 1.5,
        label = L"\pm p_{\max}")
    axislegend(ax1; position = :rb)
    lines!(ax2, v, D; color = _SERIES[1], linewidth = 2)
    hlines!(ax2, [0]; color = _OUTLINE_COLOR, linewidth = 1)
    fig
end

function _RD.plot_spline_projections(l2, kkt; markers = nothing, targets = nothing,
        size = (1200, 400))
    fig = Figure(; size = size)
    xs, vs, F2 = l2
    _, _, Fk = kkt
    Lx = xs[end]
    lim = max(maximum(abs, F2), maximum(abs, Fk))
    ax1 = Axis(fig[1, 1]; xlabel = L"x", ylabel = L"v",
        title = L"\text{(a) } L^2 \text{ projection } f_h^{L^2}")
    heatmap!(ax1, xs, vs, F2; colormap = _SIGNED_CMAP, colorrange = (-lim, lim))
    contour!(ax1, xs, vs, F2; levels = [0.0], color = _OUTLINE_COLOR, linewidth = 1.5)
    if markers !== nothing
        ParticleRemapping.plot_markers!(ax1, markers; color = _OUTLINE_COLOR,
            period = (zero(Lx), Lx))
    end
    ax2 = Axis(fig[1, 2]; xlabel = L"x", ylabel = "",
        title = L"\text{(b) KKT projection } f_h^{\mathrm{KKT}}")
    hm = heatmap!(ax2, xs, vs, Fk; colormap = _SIGNED_CMAP, colorrange = (-lim, lim))
    contour!(ax2, xs, vs, Fk; levels = [0.0], color = _OUTLINE_COLOR, linewidth = 1.5)
    targets === nothing ||
        _targets!(ax2, targets, length(targets[1]); period = (zero(Lx), Lx))
    Colorbar(fig[1, 3], hm; label = L"f_h")
    linkaxes!(ax1, ax2)
    hideydecorations!(ax2; ticks = false)
    ax3 = Axis(fig[1, 4]; xlabel = L"\min_x f_h", ylabel = L"v",
        title = "(c) smallest value at each v")
    lines!(ax3, vec(minimum(F2; dims = 1)), vs; color = _SERIES[1], label = L"L^2")
    lines!(ax3, vec(minimum(Fk; dims = 1)), vs; color = _SERIES[2], label = "KKT")
    vlines!(ax3, [0]; color = _OUTLINE_COLOR, linestyle = :dash, linewidth = 1)
    axislegend(ax3; position = :rb)
    fig
end

function _RD.plot_target_relaxation(history; Lx::Real, cv = nothing, size = (900, 400))
    fig = Figure(; size = size)
    ax = Axis(fig[1, 1]; xlabel = L"x", ylabel = L"v",
        title = (cv === nothing ? "" : "(a) ") * "relaxation of the target points")
    N = length(history[1][1])
    for a in 1:N
        xa = [mod(h[1][a], Lx) for h in history]
        va = [h[2][a] for h in history]
        # break the path where it crosses the periodic boundary
        for i in reverse(2:length(xa))
            abs(xa[i] - xa[i - 1]) > Lx / 2 && (insert!(xa, i, NaN); insert!(va, i, NaN))
        end
        lines!(ax, xa, va; color = _SEGMENT_COLOR, linewidth = 1.5)
    end
    scatter!(ax, mod.(history[1][1], Lx), history[1][2]; color = _OUTLINE_COLOR,
        markersize = 0.6 * _max_markersize(N), label = "start")
    _targets!(ax, (history[end][1], history[end][2]), N; period = (zero(Lx), Lx))
    xlims!(ax, 0, Lx)
    if cv !== nothing
        ax2 = Axis(fig[1, 2]; xlabel = "relaxation step", ylabel = "kernel-density CV",
            title = "(b) clustering")
        scatterlines!(ax2, 0:(length(cv) - 1), cv; color = _SERIES[1])
        colsize!(fig.layout, 2, Relative(0.35))
    end
    fig
end

function _RD.plot_newton_convergence(
        histories...; labels = ["$k" for k in eachindex(histories)],
        size = (760, 460))
    fig = Figure(; size = size)
    ax = Axis(fig[1, 1]; xlabel = "Newton iterate (all homotopy stages)",
        ylabel = "scaled residual", yscale = log10, title = "forward Newton solves")
    for (k, h) in enumerate(histories)
        r = [max(e[2], 1e-17) for e in h]
        scatterlines!(ax, eachindex(r), r; color = _SERIES[k], label = labels[k])
        # the first iterate of each homotopy stage
        starts = [i for i in eachindex(h) if i == 1 || h[i][1] != h[i - 1][1]]
        vlines!(ax, starts; color = (_SERIES[k], 0.3), linestyle = :dash)
    end
    hlines!(ax, [1e-14]; color = _OUTLINE_COLOR, linestyle = :dot)
    Legend(fig[2, 1], ax; orientation = :horizontal)
    fig
end

function _RD.plot_closed_loop_defect(lhs::AbstractMatrix, rhs::AbstractMatrix;
        variational = nothing, size = (1200, 330))
    fig = Figure(; size = size)
    lim = maximum(abs, lhs)
    panels = [(L"\text{(a) } D\Phi^T J \, D\Phi - J", lhs),
        (L"\text{(b) } -\tau \, (Dg^T D\Theta - D\Theta^T Dg)", rhs)]
    variational === nothing || push!(panels, ("(c) variational closure", variational))
    letter(k) = "($(Char('a' + k - 1))) "
    for (k, (title, M)) in enumerate(panels)
        ax = Axis(fig[1, k]; title = title, yreversed = true, aspect = DataAspect())
        hidedecorations!(ax)
        heatmap!(ax, transpose(M); colormap = _SIGNED_CMAP, colorrange = (-lim, lim))
    end
    Colorbar(fig[1, length(panels) + 1]; colormap = _SIGNED_CMAP, limits = (-lim, lim))
    d = lhs - rhs
    k = length(panels) + 1
    ax = Axis(fig[1, k + 1]; title = letter(k) * "difference of (a) and (b)",
        yreversed = true, aspect = DataAspect())
    hidedecorations!(ax)
    dl = maximum(abs, d)
    heatmap!(ax, transpose(d); colormap = _SIGNED_CMAP, colorrange = (-dl, dl))
    Colorbar(fig[1, k + 2]; colormap = _SIGNED_CMAP, limits = (-dl, dl))
    fig
end

function _RD.animate_target_relaxation(path::AbstractString, history; Lx::Real,
        framerate::Integer = 20, hold::Integer = 8, size = (620, 460))
    N = length(history[1][1])
    order = [fill(1, hold); eachindex(history); fill(lastindex(history), hold)]
    k = Observable(1)
    x = @lift mod.(history[$k][1], Lx)
    v = @lift history[$k][2]
    fig = Figure(; size = size)
    title = lift(i -> "de-clustering relaxation, step $(i - 1)", k)
    ax = Axis(fig[1, 1]; xlabel = L"x", ylabel = L"v", title = title)
    scatter!(ax, mod.(history[1][1], Lx), history[1][2]; color = _OUTLINE_COLOR,
        markersize = 0.6 * _max_markersize(N))
    scatter!(ax, x, v; color = :transparent, strokecolor = _TARGET_COLOR, strokewidth = 1,
        markersize = 1.6 * _max_markersize(N))
    vs = reduce(vcat, [h[2] for h in history])
    limits!(ax, 0, Lx, extrema(vs)...)
    record(fig, path, order; framerate = framerate) do i
        k[] = i
    end
    path
end
