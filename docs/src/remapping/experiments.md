# [Numerical experiments](@id remap-experiments)

```@meta
CurrentModule = ParticleRemapping.Remapping
```

The [standard experiments](@ref standard-experiments) with the gauge rezoning (`:gauge`) and the
canonical transport with the identity map (`:transport`) of
[`standard_method`](@ref ParticleRemapping.standard_method). The
[comparison](@ref comparison-standard) sets the results beside those of the other constructions.
The two numerical examples of the chapter are on the page
[clustered markers and a separable density](@ref remap-examples).

Both constructions keep every position: the gauge rezoning by construction, the transport because
its canonical map is the identity. Both keep every momentum and every sign of ``v``. They can
spread the velocities of the markers in a cell, but they cannot move a marker into another cell
of ``x`` or across ``v = 0``.

## The checkerboard

| construction | quadrants after | weight in ``Q_2 \cup Q_4`` | change of the total weight |
|:--|:--|--:|--:|
| gauge rezoning | ``(200, 0, 200, 0)`` | 0 | ``-2 \cdot 10^{-16}`` |
| transport, represented density | ``(200, 0, 200, 0)`` | 0 | ``1.6 \cdot 10^{-2}`` |
| transport, positive density | ``(200, 0, 200, 0)`` | 0 | ``1.6 \cdot 10^{-2}`` |

The empty quadrants stay empty for both densities. In a cell with ``x \ge 0`` every marker has
``v \ge 0``, so the transport matches only the branch ``v \ge 0`` of the target density there,
whatever weight the target puts on the other branch. The gauge rezoning keeps the total weight,
the total momentum and the kinetic energy to round-off, and the grid charge by its constraint.
The transport keeps every momentum; its weights follow the quantiles of the target, and the total
weight changes by the quantile error of 50 markers per cell and branch.

```@example remap-standard
using ParticleRemapping, CairoMakie
set_theme!(remap_theme())
gauge = standard_checkerboard(:gauge)
transport = standard_checkerboard(:transport; variant = :positive)
(gauge = (gauge.quadrants_after, gauge.report),
    transport = (transport.quadrants_after, transport.report))
```

```@example remap-standard
plot_remap(gauge.after, transport.after;
    titles = ("(a) after the gauge rezoning", "(b) after the transport, positive density"))
```

## The long run

| construction | empty cells (before: 14) | relative std (before: 0.598) | ``L^2`` difference | change of the total weight |
|:--|--:|--:|--:|--:|
| gauge rezoning | 13 | 0.559 | 0.069 | ``6 \cdot 10^{-15}`` |
| transport, represented density | 15 | 0.572 | 0.066 | ``7.5 \cdot 10^{-3}`` |
| transport, positive density | 0 | 0.392 | 0.328 | ``-0.20`` |

The gauge rezoning spreads the speeds in each cell and sign towards an even spacing, as far as the
homotopy reaches under the constraints of the grid charge and the kinetic energy; the number of
empty cells falls from 14 to 13. The transport with the represented density gives each marker the velocity at
its quantile of the exact density of the run; the result stays close to the markers before.
With the positive density the quantiles reach into the tails of ``f_0``, which the markers never
sampled, so the velocities spread over the whole box and every cell holds a marker. The closure
of the quantiles then puts the weight of the unsampled tails nowhere: the total weight falls by
``0.20``.

```@example remap-standard
long_gauge = standard_long_time(:gauge)
long_transport = standard_long_time(:transport; variant = :positive)
(gauge = long_gauge.stats_after[(:empty, :rstd)], gauge_l2 = long_gauge.l2_difference,
    transport = long_transport.stats_after[(:empty, :rstd)],
    transport_weight = long_transport.report.weight)
```

```@example remap-standard
plot_remap(long_gauge.after, long_transport.after; period = (-π, π),
    titles = ("(a) after the gauge rezoning", "(b) after the transport, positive density"))
```

## Remaps during the run

A remap after every 100 steps; the mean relative standard deviation of the cell counts for
``t \ge 5`` is ``0.587`` without remaps.

| construction | mean relative std, ``t \ge 5`` | empty cells at ``t = 0``, ``12.5``, ``25`` |
|:--|--:|:--|
| gauge rezoning | 0.550 | 20, 10, 8 |
| transport, represented density | 0.597 | 20, 12, 15 |
| transport, positive density | 0.384 | 20, 0, 0 |

```@example remap-standard
h_gauge = standard_remap_history(:gauge)
h_transport = standard_remap_history(:transport; variant = :positive)
plot_cell_count_history([h_gauge.h_plain, h_gauge.h_remap, h_transport.h_remap];
    labels = ["no remap", "gauge rezoning", "transport, positive density"],
    remap_times = h_gauge.remapped.remap_times)
```
