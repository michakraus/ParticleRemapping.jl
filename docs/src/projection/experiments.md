# [Numerical experiments](@id projection-experiments)

```@meta
CurrentModule = ParticleRemapping.Projection
```

The [standard experiments](@ref standard-experiments) with the kick–drift map (`:projection`)
and the non-separable map (`:nonseparable`) of [`standard_method`](@ref ParticleRemapping.standard_method).
The [comparison](@ref comparison-standard) sets the results beside those of the other
constructions. The structural checks of the chapter are on the page
[structural statements and target weights after a push](@ref projection-checks), the hole of the
chapter's own example on the page [hole filling](@ref projection-hole-filling).

## The checkerboard

With the represented density, both maps stop before the fit: the target points of the
density-interpolation flow cover the whole box, and the target weights ``c f`` vanish at the
targets in ``Q_2`` and ``Q_4``, where the represented density is zero. The construction needs a
density that is positive on its box.

```@example projection-standard
using ParticleRemapping, CairoMakie
set_theme!(remap_theme())
represented = standard_checkerboard(:projection; variant = :represented)
represented.error
```

With the positive density ``1/4``, both maps fill the empty quadrants:

| map | ``Q_1`` | ``Q_2`` | ``Q_3`` | ``Q_4`` | weight in ``Q_2 \cup Q_4`` | empty cells of 64 |
|:--|--:|--:|--:|--:|--:|--:|
| before | 200 | 0 | 200 | 0 | 0 | 32 |
| kick–drift | 131 | 73 | 130 | 66 | 0.348 | 8 |
| non-separable | 123 | 76 | 120 | 81 | 0.393 | 16 |

The exact weight in ``Q_2 \cup Q_4`` is ``1/2``. The total weight and the total momentum are kept
to round-off. The markers have equal weights, so the stage Hamiltonians do not depend on the
weight of a marker, and each map is the same area-preserving map for every marker
([one map for all markers](@ref one-map-for-all-markers)). It cannot shrink the area of the empty
region, but it can fold the occupied region into filaments that reach into every quadrant, which
fills cells but leaves empty cells between the filaments. The target weights are those of the
positive density, so every marker keeps the weight ``1/400``.

```@example projection-standard
positive = standard_checkerboard(:projection; variant = :positive)
(quadrants = positive.quadrants_after, hole_weight = positive.hole_weight,
    empty = positive.stats_after.empty, report = positive.report)
```

```@example projection-standard
plot_remap(positive.before, positive.after;
    titles = ("(a) before", "(b) after the kick–drift map, positive density"))
```

```@example projection-standard
nonseparable = standard_checkerboard(:nonseparable; variant = :positive)
(quadrants = nonseparable.quadrants_after, hole_weight = nonseparable.hole_weight,
    empty = nonseparable.stats_after.empty)
```

## The long run

With the represented density, which vanishes beyond the velocities that the initial lattice
covers, both maps stop for the same reason as on the checkerboard. With the positive density:

| map | empty cells (before: 14) | relative std (before: 0.598) | ``L^2`` difference | change of the kinetic energy | largest ``\lvert v \rvert`` |
|:--|--:|--:|--:|--:|--:|
| kick–drift | 0 | 0.475 | 0.310 | 0.289 | 6.48 |
| non-separable | 7 | 0.937 | 0.631 | 2.10 | 56.5 |

The kick–drift map leaves no cell empty. Its target weights come from the density at the targets,
and where the targets reach the tails of the density the weights are small; a small weight at a
fixed momentum gives a large velocity, up to ``6.48``. The relative ``L^2`` difference of the
smoothed densities, ``0.310``, is that of the positive density, whose tails the markers did not
sample before. The non-separable map with its five fit iterations and ``6 \times 7 \times 3``
splines reaches its targets less closely and moves some markers to velocities up to ``56.5``.

```@example projection-standard
long = standard_long_time(:projection; variant = :positive)
(stats = long.stats_after[(:empty, :min_count, :max_count, :rstd)],
    l2_difference = long.l2_difference, max_v = long.max_v, report = long.report)
```

```@example projection-standard
plot_remap(long.before, long.after; period = (-π, π),
    titles = ("(a) after 500 steps", "(b) after the kick–drift map, positive density"))
```

## Remaps during the run

With a remap after every 100 steps, the kick–drift map fails at the third remap (``t = 15``).
The first remap leaves some markers outside the velocity box ``|v| \le 3`` of the construction.
For these markers the density-interpolation flow sends the targets out of the box: the largest
target speed is ``2.9`` at ``t = 5``, ``12.3`` at ``t = 10`` and ``7.5 \cdot 10^{29}`` at
``t = 15``, where the density underflows to zero and the target weights vanish. The construction
assumes that every marker lies in its box; a repeated remap needs a box that grows with the
velocities, or a limit on the velocities that the gauge step gives
([limitations](@ref projection-limitations)). The non-separable map fails with the same error.
[`standard_remap_history`](@ref ParticleRemapping.standard_remap_history) reports the failure:

```julia
h = standard_remap_history(:projection; variant = :positive)
h.error   # "DomainError with 0.0: the target weights must be positive; …"
```

The run takes about two minutes, so the manual does not repeat it.
