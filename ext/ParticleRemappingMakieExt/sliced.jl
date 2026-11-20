# Plots of `ParticleRemapping.SlicedRedistribution`.

const _SR = ParticleRemapping.SlicedRedistribution

# The index of the slice of each marker.
function _slice_index(bins, N)
    idx = zeros(Int, N)
    for (b, members) in pairs(bins), a in members

        idx[a] = b
    end
    idx
end

# The boundaries between consecutive slices of the values `key`.
function _slice_boundaries(bins, key)
    [(maximum(key[bins[b]]) + minimum(key[bins[b + 1]])) / 2
     for b in 1:(length(bins) - 1)]
end

# The colour of the slice of each marker.
_slice_colors(bins, N) = _category_colors(length(bins))[_slice_index(bins, N)]

const _PU_LABEL = L"p / (m \, w_{\mathrm{ref}})"

function _SR.plot_sliced_stages(before::Markers, result::RemapResult; Lx = π,
        size = (1150, 640))
    info = result.info
    plan = info.plan
    N = length(before)
    m = before.m
    pu = m * info.w_ref
    ms = _max_markersize(N)
    col_v = _slice_colors(plan.bins_v, N)
    col_x = _slice_colors(plan.bins_x, N)
    fig = Figure(; size = size)
    axis(pos, title; xlabel = L"x", ylabel = _PU_LABEL) = Axis(
        fig[pos...]; xlabel, ylabel, title = _richtitle(title))

    ax = axis((1, 1), "(a) before: velocity slices"; ylabel = L"v")
    hlines!(ax, _slice_boundaries(plan.bins_v, before.v); color = _OUTLINE_COLOR,
        linestyle = :dash, linewidth = 1)
    scatter!(ax, before.x, before.v; color = col_v, markersize = ms)

    ax = axis((1, 2), "(b) Ψ_q: positions to the targets")
    vlines!(ax, plan.taus_x[1]; color = (:gray, 0.3), linewidth = 1)
    scatter!(ax, info.q̃, info.p̃ ./ pu; color = col_v, markersize = ms)
    scatter!(ax, plan.taus_x[1], fill(minimum(info.p̃ ./ pu), length(plan.taus_x[1]));
        marker = :vline, color = _TARGET_COLOR, markersize = 10)

    ax = axis((1, 3), "(c) Ψ_p: momenta to the targets")
    vlines!(ax, _slice_boundaries(plan.bins_x, info.q̃); color = _OUTLINE_COLOR,
        linestyle = :dash, linewidth = 1)
    hlines!(ax, plan.taus_p[1] ./ pu; color = (:gray, 0.3), linewidth = 1)
    scatter!(ax, info.q̂, info.p̂ ./ pu; color = col_x, markersize = ms)

    after = result.markers
    s = _weight_scale([after.w])
    ax = axis((2, 1), "(d) T_μ: weight squeeze, w̄ = w_ref μ")
    # the squeeze moves some offsets far out of the box (K6): positions wrapped into it
    ParticleRemapping.plot_markers!(ax, Markers(info.q̄, info.p̄ ./ pu, after.w; m = m);
        s.colorrange, s.colorscale, s.wref, period = (-Lx, Lx))

    ax = axis(
        (2, 2), "(e) slot multipliers"; xlabel = L"\hat{p} / (m \, w_{\mathrm{ref}})",
        ylabel = L"\mu")
    hlines!(ax, [1.0]; color = _OUTLINE_COLOR, linestyle = :dash, linewidth = 1)
    scatter!(ax, info.p̂ ./ pu, info.μ; color = _SERIES[1], markersize = ms)

    ax = axis((2, 3), "(f) gauge step: v̄ = p̄ / (m w̄)"; xlabel = L"\bar{x}",
        ylabel = L"\bar{v}")
    sc = ParticleRemapping.plot_markers!(ax, after; s.colorrange, s.colorscale, s.wref,
        period = (-Lx, Lx))
    Colorbar(fig[2, 4], sc; label = "weight w̄")
    fig
end

function _SR.plot_checkerboard_filling(ex; size = (1150, 430))
    fig = Figure(; size = size)
    from_q1 = ex.before.x .>= 0
    ms = _max_markersize(length(ex.before))
    panels = ((ex.before, ex.quadrants_before, "(a) before"),
        (ex.sliced, ex.quadrants_sliced,
            "(b) sliced lifts, B_v = $(length(ex.plan.bins_v)), " *
            "B_x = $(length(ex.plan.bins_x))"),
        (ex.marginal, ex.quadrants_marginal, "(c) one slice (marginal control)"))
    local q1, q3
    for (k, (mk, quad, label)) in pairs(panels)
        ax = Axis(fig[1, k]; xlabel = L"x", ylabel = k == 1 ? L"v" : "", aspect = 1,
            title = "$label\nQ1–Q4: $(join(quad, ", "))")
        hlines!(ax, [0.0]; color = _OUTLINE_COLOR, linestyle = :dash, linewidth = 1)
        vlines!(ax, [0.0]; color = _OUTLINE_COLOR, linestyle = :dash, linewidth = 1)
        q1 = scatter!(ax, mk.x[from_q1], mk.v[from_q1]; color = _SERIES[1], markersize = ms)
        q3 = scatter!(ax, mk.x[.!from_q1], mk.v[.!from_q1]; color = _SERIES[2],
            markersize = ms)
    end
    Legend(fig[2, 1:3], [q1, q3], ["markers from Q1", "markers from Q3"];
        orientation = :horizontal)
    fig
end

function _SR.plot_cell_count_history(
        histories; labels = ["run $k" for k in eachindex(histories)],
        remap_times = Float64[], size = (900, 600))
    fig = Figure(; size = size)
    ax1 = Axis(fig[1, 1]; ylabel = "empty cells", title = "(a) empty cells",
        xticklabelsvisible = false)
    ax2 = Axis(fig[2, 1]; xlabel = L"t", ylabel = "relative std",
        title = "(b) relative standard deviation of the cell counts")
    linkxaxes!(ax1, ax2)
    for (k, (h, label)) in enumerate(zip(histories, labels))
        lines!(ax1, h.times, h.empty; label = label, color = _SERIES[k])
        lines!(ax2, h.times, h.rstd; label = label, color = _SERIES[k])
    end
    if !isempty(remap_times)
        vlines!(ax1, remap_times; color = (:gray, 0.5), linestyle = :dash)
        vlines!(ax2, remap_times; color = (:gray, 0.5), linestyle = :dash)
    end
    axislegend(ax1; position = :rt)
    fig
end

function _SR.plot_wave_remap(ex; Lx = π, Vmax = 2.25, size = (1150, 820))
    before, after = ex.before, ex.after
    fig = Figure(; size = size)
    s = _weight_scale([before.w, after.w])
    titles = ("(a) T = 25, before the remap: $(ex.stats_before.empty) empty cells",
        "(b) after the remap: $(ex.stats_after.empty) empty cells")
    local sc
    for (k, (mk, title)) in enumerate(zip((before, after), titles))
        ax = Axis(fig[1, k]; xlabel = L"x", ylabel = k == 1 ? L"v" : "", title = title)
        vlines!(ax, range(-Lx, Lx; length = 11); color = (:gray, 0.4), linewidth = 1)
        hlines!(ax, range(-Vmax, Vmax; length = 11); color = (:gray, 0.4), linewidth = 1)
        sc = ParticleRemapping.plot_markers!(ax, mk; s.colorrange, s.colorscale, s.wref,
            period = (-Lx, Lx))
        limits!(ax, -Lx, Lx, -1.15Vmax, 1.15Vmax)
    end
    Colorbar(fig[1, 3], sc; label = "weight w")
    crange = (0, max(maximum(ex.stats_before.counts), maximum(ex.stats_after.counts)))
    local hm
    for (k, (st, label)) in enumerate(zip((ex.stats_before, ex.stats_after),
        ("(c) cell counts before", "(d) cell counts after")))
        ax = Axis(fig[2, k]; xlabel = "x cell", ylabel = k == 1 ? "v cell" : "",
            title = "$label: relative std $(round(st.rstd; digits = 3))")
        hm = heatmap!(ax, st.counts; colorrange = crange, colormap = _FIELD_CMAP)
    end
    Colorbar(fig[2, 3], hm; label = "markers per cell")
    fig
end

function _SR.animate_sliced_stages(path, before::Markers, result::RemapResult;
        framerate = 20, frames_per_stage = 40, hold = 10, Lx = π)
    info = result.info
    pu = before.m * info.w_ref
    p = ParticleRemapping.momentum(before)
    after = result.markers
    period = (-Lx, Lx)
    stages = [(_wrap(x, period), y)
              for (x, y) in ((before.x, p ./ pu), (info.q̃, info.p̃ ./ pu),
        (info.q̂, info.p̂ ./ pu), (info.q̄, info.p̄ ./ pu), (after.x, after.v))]
    # the displacement of each stage the shorter way round the periodic box
    shifts = [mod.(stages[k + 1][1] .- stages[k][1] .+ Lx, 2Lx) .- Lx
              for k in 1:(length(stages) - 1)]
    names = ["before: (x, p / (m w_ref))", "Ψ_q: lift in velocity slices",
        "Ψ_p: lift in position slices", "T_μ: weight squeeze",
        "gauge step: (x̄, v̄), v̄ = p̄ / (m w̄)"]
    xs = Observable(copy(stages[1][1]))
    ys = Observable(copy(stages[1][2]))
    title = Observable(_richtitle(names[1]))
    fig = Figure(; size = (720, 520))
    ax = Axis(fig[1, 1]; xlabel = L"x", ylabel = L"p / (m \, w_{\mathrm{ref}}),\ \bar{v}",
        title = title)
    limits!(ax, -Lx, Lx, extrema(reduce(vcat, last.(stages)))...)
    scatter!(ax, xs, ys; color = _slice_colors(info.plan.bins_v, length(before)),
        markersize = _max_markersize(length(before)))
    frames = [(1, 0.0) for _ in 1:hold]
    for k in 1:(length(stages) - 1)
        append!(frames, [(k, s) for s in range(0, 1; length = frames_per_stage + 1)[2:end]])
        append!(frames, [(k, 1.0) for _ in 1:hold])
    end
    record(fig, path, frames; framerate = framerate) do (k, s)
        xs[] = _wrap(stages[k][1] .+ s .* shifts[k], period)
        ys[] = (1 - s) .* stages[k][2] .+ s .* stages[k + 1][2]
        title[] = _richtitle(s == 0 ? names[k] : names[k + 1])
    end
    path
end

function _SR.animate_long_time_run(path, run, history; framerate = 10, Lx = π, vlim = 3.0)
    snaps = run.snapshots
    N = length(snaps[1])
    s = _weight_scale([mk.w for mk in snaps])
    fig = Figure(; size = (820, 720))
    title = Observable("t = 0")
    ax = Axis(fig[1, 1]; xlabel = L"x", ylabel = L"v", title = title)
    k_obs = Observable(1)
    mk = lift(k -> snaps[k], k_obs)
    sc = scatter!(ax, lift(m -> m.x, mk), lift(m -> m.v, mk); color = lift(m -> m.w, mk),
        colorrange = s.colorrange, colorscale = s.colorscale, colormap = _WEIGHT_CMAP,
        markersize = lift(m -> _weight_markersize(m.w, s.wref, N), mk), strokewidth = 0)
    Colorbar(fig[1, 2], sc; label = "weight w")
    limits!(ax, -Lx, Lx, -vlim, vlim)
    ax2 = Axis(fig[2, 1]; xlabel = L"t", ylabel = "cells")
    tmax = history.times[end]
    limits!(ax2, 0, tmax, 0, 1.1 * max(maximum(history.empty), 10 * maximum(history.rstd)))
    isempty(run.remap_times) ||
        vlines!(ax2, run.remap_times; color = (:gray, 0.5), linestyle = :dash)
    lines!(ax2, lift(k -> history.times[1:k], k_obs), lift(k -> history.empty[1:k], k_obs);
        label = "empty cells", color = _SERIES[1])
    lines!(ax2, lift(k -> history.times[1:k], k_obs),
        lift(k -> 10 .* history.rstd[1:k], k_obs);
        label = "10 × relative std", color = _SERIES[2])
    axislegend(ax2; position = :rt)
    rowsize!(fig.layout, 2, Relative(0.3))
    record(fig, path, eachindex(snaps); framerate = framerate) do k
        k_obs[] = k
        title[] = "t = $(round(run.times[k]; digits = 2))"
    end
    path
end
