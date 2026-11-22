# [Numerical experiments](@id sliced-experiments)

```@meta
CurrentModule = ParticleRemapping.SlicedRedistribution
```

The [standard experiments](@ref standard-experiments) with the sliced redistribution (`:sliced`
of [`standard_method`](@ref ParticleRemapping.standard_method)). The
[comparison](@ref comparison-standard) sets the results beside those of the other constructions.
The construction computes everything from the markers, so the two density variants give the same
result. The structural checks are on the page [structural suites](@ref sliced-suites), the lifts
alone on the checkerboard on the page [hole filling](@ref sliced-holes), and the variants of the
parameters of the long run on the page
[long-time run: variants and repeated remaps](@ref sliced-long-time-page).

## The checkerboard

The full remap, with the lifts, the squeeze and the gauge step, on ``[-1, 1)^2`` with
``B_v = 10``, ``B_x = 20`` and the kernel widths ``h_x = h_v = 0.15``:

| | ``Q_1`` | ``Q_2`` | ``Q_3`` | ``Q_4`` | weight in ``Q_2 \cup Q_4`` | empty cells of 64 | relative std |
|:--|--:|--:|--:|--:|--:|--:|--:|
| before | 200 | 0 | 200 | 0 | 0 | 32 | 1.022 |
| after | 134 | 67 | 128 | 71 | 0.205 | 2 | 0.606 |

The lifts couple the markers through their slices and move markers into the empty quadrants
([hole filling](@ref sliced-holes)). The squeeze then sets the weights from a kernel estimate of
the density of the input, which is small in ``Q_2`` and ``Q_4``: the markers there carry the weight
``0.205`` together, against ``1/2`` for the positive density and zero for the represented one.
The squeeze moves some positions by up to the half width of the box (K6). The total weight and
the total momentum are kept to round-off.

```@example sliced-standard
using ParticleRemapping, CairoMakie
set_theme!(remap_theme())
board = standard_checkerboard(:sliced)
(quadrants = board.quadrants_after, hole_weight = board.hole_weight,
    stats = board.stats_after[(:empty, :rstd)], max_shift = board.max_shift,
    report = board.report)
```

```@example sliced-standard
plot_remap(board.before, board.after; period = (-1, 1),
    titles = ("(a) before", "(b) after the sliced remap"))
```

## The long run

The remap at ``T = 25`` with the default [`SlicedQuantileRemap`](@ref) is suite 5b of the source
script:

| | empty cells | fewest, most markers in a cell | relative std | ``L^2`` difference | change of the kinetic energy |
|:--|--:|:--|--:|--:|--:|
| before | 14 | 0, 20 | 0.598 | | |
| after | 0 | 1, 15 | 0.367 | 0.104 | 0.088 |

```@example sliced-standard
long = standard_long_time(:sliced)
(stats = long.stats_after[(:empty, :min_count, :max_count, :rstd)],
    l2_difference = long.l2_difference, report = long.report)
```

```@example sliced-standard
plot_remap(long.before, long.after; period = (-π, π),
    titles = ("(a) after 500 steps", "(b) after the sliced remap"))
```

## Remaps during the run

With a remap after every 100 steps, no cell is empty right after a remap, and at most two cells
empty again between two remaps; the mean relative standard deviation of the cell counts for
``t \ge 5`` falls from ``0.587`` to ``0.365``. The
[long-time page](@ref sliced-long-time-page) shows the run as an animation.

```@example sliced-standard
h = standard_remap_history(:sliced)
plot_cell_count_history([h.h_plain, h.h_remap]; labels = ["no remap", "sliced remap"],
    remap_times = h.remapped.remap_times)
```
