# Plots of `ParticleRemapping.Redistribution`.

const _RD = ParticleRemapping.Redistribution

function _RD.plot_velocity_momentum_fold(v::AbstractVector, P::AbstractVector,
        D::AbstractVector; vth::Real, pmax::Real, size = (900, 360))
    fig = Figure(; size = size)
    ax1 = Axis(fig[1, 1]; xlabel = "v̄", ylabel = "p̄ = m ΔΩ v̄ f_h(v̄)",
        title = "canonical momentum of a density-tied weight")
    ax2 = Axis(fig[1, 2]; xlabel = "v̄", ylabel = "det D_v̄ P", title = "Jacobian determinant")
    for ax in (ax1, ax2)
        vspan!(ax, [-maximum(abs, v), vth], [-vth, maximum(abs, v)];
            color = (:grey, 0.12))
        vlines!(ax, [-vth, vth]; color = :grey, linestyle = :dash)
    end
    lines!(ax1, v, P; color = :black, linewidth = 2)
    hlines!(ax1, [-pmax, pmax]; color = :firebrick, linestyle = :dot, label = "±p_max")
    axislegend(ax1; position = :rb)
    lines!(ax2, v, D; color = :black, linewidth = 2)
    hlines!(ax2, [0]; color = :grey)
    fig
end

function _RD.plot_spline_projections(l2, kkt; markers = nothing, targets = nothing,
        size = (1200, 400))
    fig = Figure(; size = size)
    xs, vs, F2 = l2
    _, _, Fk = kkt
    lim = maximum(abs, F2)
    ax1 = Axis(fig[1, 1]; xlabel = "x", ylabel = "v", title = "L² projection f_h^{L²}")
    hm1 = heatmap!(ax1, xs, vs, F2; colormap = :RdBu, colorrange = (-lim, lim))
    contour!(ax1, xs, vs, F2; levels = [0.0], color = :black, linewidth = 1.5)
    markers === nothing ||
        scatter!(ax1, mod.(markers.x, xs[end]), markers.v; color = :black, markersize = 6)
    Colorbar(fig[1, 2], hm1)
    ax2 = Axis(fig[1, 3]; xlabel = "x", ylabel = "v", title = "KKT projection f_h^KKT")
    hm2 = heatmap!(ax2, xs, vs, Fk; colormap = :RdBu, colorrange = (-lim, lim))
    contour!(ax2, xs, vs, Fk; levels = [0.0], color = :black, linewidth = 1.5)
    targets === nothing || scatter!(ax2, mod.(targets[1], xs[end]), targets[2];
        color = :white, strokecolor = :black, strokewidth = 1, markersize = 7)
    Colorbar(fig[1, 4], hm2)
    linkaxes!(ax1, ax2)
    ax3 = Axis(fig[1, 5]; xlabel = "min over x", ylabel = "v",
        title = "smallest value at each v")
    lines!(ax3, vec(minimum(F2; dims = 1)), vs; color = :firebrick, label = "L²")
    lines!(ax3, vec(minimum(Fk; dims = 1)), vs; color = :navy, label = "KKT")
    vlines!(ax3, [0]; color = :grey, linestyle = :dash)
    axislegend(ax3; position = :rb)
    fig
end

function _RD.plot_target_relaxation(history; Lx::Real, cv = nothing, size = (900, 400))
    fig = Figure(; size = size)
    ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "v", title = "relaxation of the target points")
    N = length(history[1][1])
    for a in 1:N
        xa = [mod(h[1][a], Lx) for h in history]
        va = [h[2][a] for h in history]
        lines!(ax, xa, va; color = (:steelblue, 0.8))
    end
    scatter!(ax, mod.(history[1][1], Lx), history[1][2]; color = :grey, markersize = 8,
        label = "start")
    scatter!(
        ax, mod.(history[end][1], Lx), history[end][2]; color = :black, markersize = 8,
        label = "targets")
    axislegend(ax; position = :rt)
    if cv !== nothing
        ax2 = Axis(fig[1, 2]; xlabel = "relaxation step", ylabel = "kernel-density CV",
            title = "clustering")
        scatterlines!(ax2, 0:(length(cv) - 1), cv; color = :black)
        colsize!(fig.layout, 2, Relative(0.35))
    end
    fig
end

function _RD.plot_redistribution_markers(
        before::Markers, after::Markers; targets = nothing,
        Lx::Real = 2π, size = (980, 400))
    fig = ParticleRemapping.plot_remap(before, after; size = size,
        titles = ("before: clustered markers", "after: redistributed markers"))
    if targets !== nothing
        ax = contents(fig[1, 2])[1]
        scatter!(ax, mod.(targets[1], Lx), targets[2]; color = :transparent,
            strokecolor = :red, strokewidth = 1, markersize = 16)
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
        scatterlines!(ax, eachindex(r), r; label = labels[k])
        # the first iterate of each homotopy stage
        starts = [i for i in eachindex(h) if i == 1 || h[i][1] != h[i - 1][1]]
        vlines!(ax, starts; color = (Makie.wong_colors()[k], 0.3), linestyle = :dash)
    end
    hlines!(ax, [1e-14]; color = :grey, linestyle = :dot)
    Legend(fig[2, 1], ax; orientation = :horizontal, framevisible = false)
    fig
end

function _RD.plot_closed_loop_defect(lhs::AbstractMatrix, rhs::AbstractMatrix;
        variational = nothing, size = (1200, 330))
    fig = Figure(; size = size)
    lim = maximum(abs, lhs)
    panels = [("DΦᵀ J DΦ − J", lhs), ("−τ (Dgᵀ DΘ − DΘᵀ Dg)", rhs)]
    variational === nothing || push!(panels, ("variational closure", variational))
    for (k, (t, M)) in enumerate(panels)
        ax = Axis(fig[1, k]; title = t, yreversed = true, aspect = DataAspect())
        hidedecorations!(ax)
        heatmap!(ax, transpose(M); colormap = :RdBu, colorrange = (-lim, lim))
    end
    Colorbar(fig[1, length(panels) + 1]; colormap = :RdBu, limits = (-lim, lim))
    d = lhs - rhs
    ax = Axis(fig[1, length(panels) + 2]; title = "difference of the first two",
        yreversed = true, aspect = DataAspect())
    hidedecorations!(ax)
    dl = maximum(abs, d)
    heatmap!(ax, transpose(d); colormap = :RdBu, colorrange = (-dl, dl))
    Colorbar(fig[1, length(panels) + 3]; colormap = :RdBu, limits = (-dl, dl))
    fig
end

function _RD.animate_redistribution(path::AbstractString, states, targets; p_scale::Real,
        framerate::Integer = 12, size = (620, 460))
    n = length(targets) ÷ 2
    k = Observable(1)
    q = @lift states[$k][1:n]
    p = @lift states[$k][(n + 1):(2n)] ./ p_scale
    fig = Figure(; size = size)
    title = lift(k) do i
        "symplectic step, pseudo-time s = $(round((i - 1) / (length(states) - 1); digits = 2))"
    end
    ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "p / p_s", title = title)
    scatter!(ax, targets[1:n], targets[(n + 1):(2n)] ./ p_scale; color = :transparent,
        strokecolor = :red, strokewidth = 1.5, markersize = 14, label = "targets")
    scatter!(ax, states[1][1:n], states[1][(n + 1):(2n)] ./ p_scale; color = (:grey, 0.5),
        markersize = 7, label = "start")
    scatter!(ax, q, p; color = :black, markersize = 9, label = "markers")
    axislegend(ax; position = :rt)
    record(fig, path, eachindex(states); framerate = framerate) do i
        k[] = i
    end
    path
end

function _RD.animate_target_relaxation(path::AbstractString, history; Lx::Real,
        framerate::Integer = 6, size = (620, 460))
    k = Observable(1)
    x = @lift mod.(history[$k][1], Lx)
    v = @lift history[$k][2]
    fig = Figure(; size = size)
    title = lift(i -> "de-clustering relaxation, step $(i - 1)", k)
    ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "v", title = title)
    scatter!(
        ax, mod.(history[1][1], Lx), history[1][2]; color = (:grey, 0.5), markersize = 8)
    scatter!(ax, x, v; color = :black, markersize = 10)
    xlims!(ax, 0, Lx)
    record(fig, path, eachindex(history); framerate = framerate) do i
        k[] = i
    end
    path
end
