# [Long-time run: variants and repeated remaps](@id sliced-long-time-page)

```@meta
CurrentModule = ParticleRemapping.SlicedRedistribution
```

This page gives the sliced-specific parts of the long-time run: the variants of the parameters
and the run with repeated remaps. The [standard experiments](@ref sliced-experiments) compare the
same run across all constructions.

## [The remap after the long run](@id sliced-long-time)

[`long_time_remap_experiment`](@ref): ``30 \times 30`` markers at the cell centres of
``[-\pi, \pi) \times [-1.4, 1.4]`` with weights ``\propto f_0(x, v)``,
``f_0 = (1 + 0.35 \cos x)(0.40 + v^2) e^{-v^2/1.4}`` (``\sum w = 1``, ``m = 1.2``), pushed by 500
drift–kick–drift steps of ``\Delta t = 0.05`` through the field
``E(t, x) = -0.85 \sin x - 0.75 \sin(2x - 1.3 t)`` to ``T = 25``
([`electric_wave_markers`](@ref), [`electric_wave_simulation`](@ref)), then one
[`remap`](@ref) with the default [`SlicedQuantileRemap`](@ref). The statistics count the markers
in the ``10 \times 10`` cells of ``[-\pi, \pi) \times [-2.25, 2.25]``
([`phase_cell_statistics`](@ref)); the ``L^2`` difference compares kernel-smoothed densities on a
``30 \times 30`` grid ([`smoothed_density_difference`](@ref), ``h_x = 0.45``, ``h_v = 0.35``).

| quantity | before | after | threshold |
|:--|--:|--:|:--|
| empty cells | 14 | 0 | before ``\ge 10``, after ``= 0`` |
| fewest, most markers in a cell | 0, 20 | 1, 15 | |
| relative std of the cell counts | 0.598 | 0.367 | after ``< 0.65 \times`` before |
| relative ``L^2`` difference of the smoothed densities | | 0.104 | ``< 0.12`` |
| total momentum error | | ``1.1 \times 10^{-16}`` | ``< 10^{-15}`` |
| total weight error | | 0 | ``< 10^{-15}`` |
| smallest eigenvalue of ``\nabla^2 S_1``, ``\nabla^2 S_2`` | | 0.2226, 0.2187 | ``\ge 0.22``, ``\ge 0.20`` |

The multipliers of this remap depend on the state; the map is canonical for the frozen
multipliers ([theorem](@ref sliced-composed-map)).

```@example sliced-experiments
using ParticleRemapping, CairoMakie
set_theme!(remap_theme())
ex = long_time_remap_experiment()
plot_wave_remap(ex)
```

### Variants

`scripts/sliced/long_time_experiment.jl` runs the same remap with other parameters. Columns:
empty cells, fewest and most markers per cell and relative std after the remap; the ``L^2``
difference; the smallest eigenvalues ``\lambda_q``, ``\lambda_p`` of the two Hessians; the largest
distance of ``\hat{x}`` and of ``\bar{x}`` from the anchors (the fields `q̂`, `q̄` and `q_anc`);
the largest ``|\bar{v}|``.

| variant | empty | min | max | rstd | ``L^2`` | ``\lambda_q`` | ``\lambda_p`` | ``\lvert\hat{x} - x^{\mathrm{anc}}\rvert`` | ``\lvert\bar{x} - x^{\mathrm{anc}}\rvert`` | ``\lvert\bar{v}\rvert`` |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| defaults | 0 | 1 | 15 | 0.367 | 0.104 | 0.2226 | 0.2187 | 0.361 | 3.46 | 2.87 |
| no squeeze (control) | 3 | 0 | 18 | 0.740 | 0.337 | 0.2226 | 0.2187 | 0.361 | 3.46 | 31.9 |
| ``\beta_{\mathrm{mean}} = \beta_{\mathrm{intra}} = 0.88`` | 0 | 3 | 18 | 0.337 | 0.096 | 0.1200 | 0.1200 | 0.921 | 13.2 | 2.88 |
| ``\varepsilon`` = 0.1 spacings | 1 | 0 | 17 | 0.398 | 0.097 | 0.2200 | 0.2000 | 0.559 | 2.48 | 2.74 |
| ``\varepsilon`` = 0.5 spacings | 1 | 0 | 16 | 0.369 | 0.096 | 0.2202 | 0.2006 | 0.438 | 1.98 | 2.72 |
| ``\varepsilon`` = 1 spacing | 0 | 1 | 16 | 0.353 | 0.097 | 0.2206 | 0.2034 | 0.374 | 2.48 | 2.77 |

In the control without the squeeze the weights ``w_{\mathrm{ref}} \mu`` are set at the fixed
momenta ``\hat{p}``, and ``\bar{v} = \hat{p}/(m\bar{w})`` stretches where the weights are small
(the [failures](@ref sliced-failures) page). The other variants keep every cell occupied or leave
one empty; they move the positions further from the anchors (``\beta = 0.88``) or the lifted
positions (small ``\varepsilon``). The squeeze moves some markers by up to 3.46 from their anchors
in the default run: the slots with small ``\mu_a``.

### Remaps during the run

The same run with a remap after every 100 steps (``t = 5, 10, 15, 20``), statistics every 10
steps. At ``t = 0`` the lattice covers ``|v| \le 1.4`` only, hence the empty cells of the
statistics box ``|v| \le 2.25``.

| ``t`` | empty cells, no remap | rstd, no remap | empty cells, remaps | rstd, remaps |
|--:|--:|--:|--:|--:|
| 0 | 20 | 0.699 | 20 | 0.699 |
| 2.5 | 16 | 0.664 | 16 | 0.664 |
| 5 | 17 | 0.663 | 0 | 0.389 |
| 10 | 17 | 0.592 | 0 | 0.304 |
| 15 | 13 | 0.569 | 0 | 0.349 |
| 20 | 13 | 0.600 | 0 | 0.315 |
| 25 | 14 | 0.598 | 0 | 0.427 |

The mean relative standard deviation for ``t \ge 5`` is 0.587 without and 0.365 with the remaps
(`scripts/sliced/long_time_experiment.jl`).

```@example sliced-experiments
initial = electric_wave_markers()
plain = electric_wave_simulation(initial; record_every = 10)
remapped = electric_wave_simulation(initial; record_every = 10,
    remap_method = SlicedQuantileRemap(), remap_every = 100)
h_plain = cell_count_history(plain.snapshots, plain.times)
h_remap = cell_count_history(remapped.snapshots, remapped.times)
plot_cell_count_history([h_plain, h_remap]; labels = ["no remap", "remap every 100 steps"],
    remap_times = remapped.remap_times)
```

The animations show the two runs (colour and size: weight; below: empty cells and ten times the
relative standard deviation).

```@example sliced-experiments
animate_long_time_run("long-time-run.mp4", remapped, h_remap)
animate_long_time_run("long-time-run-no-remap.mp4", plain, h_plain)
nothing # hide
```

With a remap after every 100 steps:

```@raw html
<video src="./long-time-run.mp4" controls loop muted playsinline style="width:100%"></video>
```

Without remaps:

```@raw html
<video src="./long-time-run-no-remap.mp4" controls loop muted playsinline style="width:100%"></video>
```
