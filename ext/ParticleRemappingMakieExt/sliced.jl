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

function _stage_axis(fig, pos, title; ylabel = "p / (m w_ref)")
    Axis(fig[pos...]; xlabel = "x", ylabel = ylabel, title = title, titlesize = 13)
end

function _SR.plot_sliced_stages(before::Markers, result::RemapResult; Lx = π,
        size = (1500, 820))
    info = result.info
    plan = info.plan
    N = length(before)
    m = before.m
    pu = m * info.w_ref
    slice_v = _slice_index(plan.bins_v, N)
    slice_x = _slice_index(plan.bins_x, N)
    B_v, B_x = length(plan.bins_v), length(plan.bins_x)
    fig = Figure(; size = size)

    ax = _stage_axis(fig, (1, 1), "input and velocity slices"; ylabel = "v")
    hlines!(ax, _slice_boundaries(plan.bins_v, before.v); color = (:gray, 0.6), linestyle = :dash)
    scatter!(
        ax, before.x, before.v; color = slice_v, colormap = :turbo, colorrange = (1, B_v),
        markersize = 7)

    ax = _stage_axis(fig, (1, 2), "Ψ_q: positions to the targets of each velocity slice")
    vlines!(ax, plan.taus_x[1]; color = (:gray, 0.3))
    scatter!(ax, info.q̃, info.p̃ ./ pu; color = slice_v, colormap = :turbo,
        colorrange = (1, B_v), markersize = 7)
    scatter!(ax, plan.taus_x[1], fill(minimum(info.p̃ ./ pu), length(plan.taus_x[1]));
        marker = :vline, color = :black, markersize = 10)

    ax = _stage_axis(fig, (1, 3), "Ψ_p: momenta to the targets of each position slice")
    vlines!(ax, _slice_boundaries(plan.bins_x, info.q̃); color = (:gray, 0.5), linestyle = :dash)
    hlines!(ax, plan.taus_p[1] ./ pu; color = (:gray, 0.3))
    scatter!(ax, info.q̂, info.p̂ ./ pu; color = slice_x, colormap = :phase,
        colorrange = (1, B_x), markersize = 7)

    ax = _stage_axis(fig, (2, 1), "T_μ: canonical weight squeeze")
    sc = scatter!(ax, info.q̄, info.p̄ ./ pu; color = info.μ, colormap = :viridis,
        markersize = 7)
    Colorbar(fig[2, 1][1, 2], sc; label = "μ")

    ax = Axis(
        fig[2, 2]; xlabel = "p̂ / (m w_ref)", ylabel = "μ", title = "slot multipliers",
        titlesize = 13)
    hlines!(ax, [1.0]; color = (:gray, 0.6), linestyle = :dash)
    scatter!(ax, info.p̂ ./ pu, info.μ; color = :black, markersize = 6)

    ax = Axis(fig[2, 3]; xlabel = "x̄", ylabel = "v̄",
        title = "F: weights w̄ = w_ref μ, v̄ = p̄/(m w̄)",
        titlesize = 13)
    after = result.markers
    shown = Markers(ParticleRemapping.periodic_wrap.(after.x, Lx), after.v, after.w; m = m)
    sc = ParticleRemapping.plot_markers!(ax, shown)
    Colorbar(fig[2, 3][1, 2], sc; label = "weight")
    fig
end

function _SR.plot_checkerboard_filling(ex; size = (1300, 430))
    fig = Figure(; size = size)
    panels = ((ex.before, ex.quadrants_before, "before"),
        (ex.sliced,
            ex.quadrants_sliced,
            "sliced lifts (B_v = $(length(ex.plan.bins_v)), " *
            "B_x = $(length(ex.plan.bins_x)))"),
        (ex.marginal, ex.quadrants_marginal, "one slice (marginal control)"))
    for (k, (mk, quad, label)) in pairs(panels)
        ax = Axis(fig[1, k]; xlabel = "x", ylabel = "v", aspect = 1,
            title = "$label\nQ1–Q4: $(join(quad, ", "))", titlesize = 13)
        hlines!(ax, [0.0]; color = :gray)
        vlines!(ax, [0.0]; color = :gray)
        scatter!(ax, mk.x, mk.v; color = :black, markersize = 4)
    end
    fig
end

function _SR.plot_cell_count_history(
        histories; labels = ["run $k" for k in eachindex(histories)],
        remap_times = Float64[], size = (900, 600))
    fig = Figure(; size = size)
    ax1 = Axis(fig[1, 1]; ylabel = "empty cells")
    ax2 = Axis(fig[2, 1]; xlabel = "t", ylabel = "relative std of cell counts")
    linkxaxes!(ax1, ax2)
    for (h, label) in zip(histories, labels)
        lines!(ax1, h.times, h.empty; label = label)
        lines!(ax2, h.times, h.rstd; label = label)
    end
    if !isempty(remap_times)
        vlines!(ax1, remap_times; color = (:gray, 0.5), linestyle = :dash)
        vlines!(ax2, remap_times; color = (:gray, 0.5), linestyle = :dash)
    end
    axislegend(ax1; position = :rt)
    fig
end

function _SR.plot_wave_remap(ex; Lx = π, Vmax = 2.25, size = (1300, 880))
    before = ex.before
    after = Markers(
        ParticleRemapping.periodic_wrap.(ex.after.x, Lx), ex.after.v, ex.after.w;
        m = ex.after.m)
    fig = Figure(; size = size)
    cr = extrema(vcat(before.w, after.w))
    titles = ("T = 25, before the remap: $(ex.stats_before.empty) empty cells",
        "after the remap: $(ex.stats_after.empty) empty cells")
    local sc
    for (k, (mk, title)) in enumerate(zip((before, after), titles))
        ax = Axis(fig[1, k]; xlabel = "x", ylabel = "v", title = title)
        vlines!(ax, range(-Lx, Lx; length = 11); color = (:gray, 0.4))
        hlines!(ax, range(-Vmax, Vmax; length = 11); color = (:gray, 0.4))
        sc = ParticleRemapping.plot_markers!(ax, mk; colorrange = cr)
        limits!(ax, -Lx, Lx, -1.15Vmax, 1.15Vmax)
    end
    Colorbar(fig[1, 3], sc; label = "weight")
    crange = (0, max(maximum(ex.stats_before.counts), maximum(ex.stats_after.counts)))
    local hm
    for (k, (s, label)) in enumerate(zip((ex.stats_before, ex.stats_after), (
        "before", "after")))
        ax = Axis(fig[2, k]; xlabel = "x cell", ylabel = "v cell",
            title = "cell counts $label: relative std $(round(s.rstd; digits = 3))")
        hm = heatmap!(ax, s.counts; colorrange = crange, colormap = :magma)
    end
    Colorbar(fig[2, 3], hm; label = "markers per cell")
    fig
end

function _SR.animate_sliced_stages(path, before::Markers, result::RemapResult;
        framerate = 30, frames_per_stage = 40)
    info = result.info
    pu = before.m * info.w_ref
    p = ParticleRemapping.momentum(before)
    after = result.markers
    stages = [(before.x, p ./ pu), (info.q̃, info.p̃ ./ pu), (info.q̂, info.p̂ ./ pu),
        (info.q̄, info.p̄ ./ pu), (after.x, after.v)]
    names = ["input (x, p/(m w_ref))", "Ψ_q: equalise positions in velocity slices",
        "Ψ_p: equalise momenta in position slices", "T_μ: weight squeeze",
        "F: final weights, (x̄, v̄)"]
    colour = _slice_index(info.plan.bins_v, length(before))
    xs = Observable(copy(stages[1][1]))
    ys = Observable(copy(stages[1][2]))
    title = Observable(names[1])
    fig = Figure(; size = (760, 560))
    ax = Axis(fig[1, 1]; xlabel = "position", ylabel = "momentum / (m w_ref), velocity",
        title = title)
    allx = reduce(vcat, first.(stages))
    ally = reduce(vcat, last.(stages))
    limits!(ax, extrema(allx)..., extrema(ally)...)
    scatter!(ax, xs, ys; color = colour, colormap = :turbo, markersize = 7)
    frames = [(k, s) for k in 1:(length(stages) - 1)
              for s in range(0, 1; length = frames_per_stage)]
    record(fig, path, frames; framerate = framerate) do (k, s)
        xs[] = (1 - s) .* stages[k][1] .+ s .* stages[k + 1][1]
        ys[] = (1 - s) .* stages[k][2] .+ s .* stages[k + 1][2]
        title[] = s < 1 ? "$(names[k]) → $(names[k + 1])" : names[k + 1]
    end
    path
end

function _SR.animate_long_time_run(path, run, history; framerate = 10, Lx = π, vlim = 3.0)
    snaps = run.snapshots
    fig = Figure(; size = (900, 760))
    title = Observable("t = 0")
    ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "v", title = title)
    cr = extrema(reduce(vcat, [mk.w for mk in snaps]))
    xs = Observable(snaps[1].x)
    vs = Observable(snaps[1].v)
    ws = Observable(snaps[1].w)
    sc = scatter!(
        ax, xs, vs; color = ws, colorrange = cr, colormap = :viridis, markersize = 5)
    Colorbar(fig[1, 2], sc; label = "weight")
    limits!(ax, -Lx, Lx, -vlim, vlim)
    ax2 = Axis(fig[2, 1]; xlabel = "t", ylabel = "empty cells / relative std × 10")
    tmax = history.times[end]
    limits!(ax2, 0, tmax, 0, 1.1 * max(maximum(history.empty), 10 * maximum(history.rstd)))
    isempty(run.remap_times) || vlines!(ax2, run.remap_times; color = (:gray, 0.5))
    k_obs = Observable(1)
    lines!(ax2, lift(k -> history.times[1:k], k_obs), lift(k -> history.empty[1:k], k_obs);
        label = "empty cells")
    lines!(ax2, lift(k -> history.times[1:k], k_obs),
        lift(k -> 10 .* history.rstd[1:k], k_obs);
        label = "10 × relative std")
    axislegend(ax2; position = :rt)
    rowsize!(fig.layout, 2, Relative(0.3))
    record(fig, path, eachindex(snaps); framerate = framerate) do k
        xs[] = snaps[k].x
        vs[] = snaps[k].v
        ws[] = snaps[k].w
        title[] = "t = $(round(run.times[k]; digits = 2))"
        k_obs[] = k
    end
    path
end
