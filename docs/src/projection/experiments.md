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

With the positive density ``1/4``, both maps move markers into the empty quadrants:

| map | ``Q_1`` | ``Q_2`` | ``Q_3`` | ``Q_4`` | weight in ``Q_2 \cup Q_4`` | empty cells of 64 | stop |
|:--|--:|--:|--:|--:|--:|--:|:--|
| before | 200 | 0 | 200 | 0 | 0 | 32 | |
| kick–drift | 160 | 40 | 165 | 35 | 0.1875 | 13 | `MAXITERS` |
| non-separable | 127 | 72 | 128 | 73 | 0.3625 | 19 | `MAXITERS` |

The fit is GeometricSolvers' Gauss–Newton method with the `DogLeg` globalisation, with the
iteration budget of [`standard_method`](@ref ParticleRemapping.standard_method): 20 iterations for
the kick–drift map and 5 for the non-separable map. The column stop gives its return code;
`MAXITERS` means that the fit ends on its budget, so the numbers are those after the budget, not
at a minimum. A number of this page is quoted with its digits where, with the targets moved by
``10^{-13}`` relative and on Julia 1.13 and 1.14, it does not move at these digits and stays at
least a tenth of a unit of its last digit from a rounding boundary; otherwise it is quoted as the
range over these runs.

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

| map | empty cells (before: 14) | relative std (before: 0.598) | ``L^2`` difference | change of the kinetic energy | largest ``\lvert v \rvert`` | stop |
|:--|--:|--:|--:|--:|--:|:--|
| kick–drift | 3–6 | 0.808–0.837 | 0.494–0.507 | 1.36–1.43 | 32.4–33.4 | `MAXITERS` |
| non-separable | 12 | 1.01 | 0.623 | 2.144 | 65.16 | `MAXITERS` |

Both fits end on their budgets, and both maps leave cells empty. The velocity of a marker after
the map is ``p^S / (m w^*)``, its momentum after the map over its mass and its target weight, so a
momentum that misses its target at a small target weight gives a large velocity: up to
``32.4``–``33.4`` after the kick–drift map, and up to ``65.16`` after the non-separable map with
its ``6 \times 7 \times 3`` splines and five fit iterations.

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

With a remap after every 100 steps, the kick–drift map fails at the second remap (``t = 10``) on
Julia 1.13 and at the third remap (``t = 15``) on Julia 1.14. Its fit at ``t = 5`` ends on its
20 iterations (`MAXITERS`), and the first remap leaves some markers outside the velocity box
``|v| \le 3`` of the construction: the largest speed after it is 26.8–36.2. The run then fails
through a target weight that is zero or that underflows. On Julia 1.13 a target weight at
``t = 10`` is zero, and that remap throws `DomainError with 0.0`. On Julia 1.14 the remap at
``t = 10`` passes, and its fit ends on `MAXITERS` after 20 iterations or on `SUCCESS` after one;
where it ends on `SUCCESS`, target weights below ``10^{-200}`` pass the check of positive target
weights (`KNOWN_ISSUES.md`, K18). The remap at ``t = 15`` then throws `DomainError with NaN`.
The construction assumes that every marker lies in its box; a repeated remap needs a box that
grows with the velocities, or a limit on the velocities that the gauge step gives
([limitations](@ref projection-limitations)). The fit of the non-separable map ends on its five
iterations (`MAXITERS`) at ``t = 5`` and at ``t = 10``. The largest speed is ``31.44`` after the
first remap and ``2.6 \cdot 10^{303}`` after the second, where a target weight near the underflow
threshold passes the check of positive target weights (`KNOWN_ISSUES.md`, K18); the remap at
``t = 15`` throws `DomainError with NaN`.
[`standard_remap_history`](@ref ParticleRemapping.standard_remap_history) reports the failure:

```julia
h = standard_remap_history(:projection; variant = :positive)
h.error   # "DomainError with 0.0: …" on Julia 1.13, "DomainError with NaN: …" on 1.14
```

The run takes about two minutes, so the manual does not repeat it.
