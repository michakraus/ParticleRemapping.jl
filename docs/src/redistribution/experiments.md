# [Numerical experiments](@id rd-experiments)

```@meta
CurrentModule = ParticleRemapping.Redistribution
```

The [standard experiments](@ref standard-experiments) with the spline-projection redistribution
(`:spline` of [`standard_method`](@ref ParticleRemapping.standard_method)). The
[comparison](@ref comparison-standard) sets the results beside those of the other constructions.
The construction computes everything from the markers, so the two density variants give the same
result. Its own example is on the page [the verification set-up](@ref rd-setup), and a lattice with
an elliptic hole on the page [hole filling](@ref rd-holes).

## The full sizes

At 400 and at 900 markers the construction stops in the KKT solve, as for the 900 markers of
[K2](@ref rd-experiments-failures): the relaxation sums the repulsion over all pairs, its default
step fits 16 markers, and the relaxed targets leave the velocity box of the spline basis.

```@example rd-standard
using ParticleRemapping, CairoMakie
set_theme!(remap_theme())
(checkerboard = standard_checkerboard(:spline).error,
    long_run = standard_long_time(:spline).error)
```

With the relaxation step scaled to ``0.006 \cdot 16 / N``, the 400 markers of the checkerboard
reach the end of the remap after 316 s, with the scaled distance ``3.1 \cdot 10^9`` of the end
point from the targets: the symplectic step misses its targets, as in [K3](@ref rd-holes).
The 900 markers of the long run did not finish within 20 minutes.

## Reduced sizes

The reduced runs scale the relaxation step to ``0.006 \cdot 16 / N``. The checkerboard holds
``N`` markers from the same seed, half in ``Q_1`` and half in ``Q_3``; the long run starts from a
coarser lattice.

| marker set | ``N`` | quadrants or empty cells after | scaled distance from the targets |
|:--|--:|:--|--:|
| checkerboard | 16 | ``(5, 1, 5, 5)`` | 116 |
| checkerboard | 36 | ``(9, 11, 10, 6)`` | ``4.4 \cdot 10^3`` |
| checkerboard | 64 | ``(23, 12, 22, 7)`` | ``4.6 \cdot 10^6`` |
| long run, ``4 \times 4`` lattice | 16 | 85 of 100 cells empty (before: 86) | 1.15 |
| long run, ``6 \times 6`` lattice | 36 | 90 of 100 cells empty (before: 70) | 152 |
| long run, ``8 \times 8`` lattice | 64 | 63 of 100 cells empty (before: 58) | ``3.3 \cdot 10^4`` |

The pairwise Hamiltonian couples the markers, so the step can move markers into the empty
quadrants. It reaches its targets only for the 16 markers of the long run, the size and the kind
of distribution of the [verification set-up](@ref rd-setup); for every other set the distance from
the targets grows with ``N``, and the end point is not the step that the parameters describe. The
total weight is kept to round-off in every run, the total momentum where the Newton solve of the
step converges.

```@example rd-standard
small = standard_checkerboard(:spline; N = 16)
(quadrants = small.quadrants_after, target_error = small.result.info.target_error,
    report = small.report)
```

```@example rd-standard
plot_remap(small.before, small.after;
    targets = (small.result.info.targets_x .- 1, small.result.info.targets_v),
    titles = ("(a) 16 checkerboard markers", "(b) after the redistribution"))
```

```@example rd-standard
coarse = standard_long_time(:spline; Nx = 4, Nv = 4)
(target_error = coarse.result.info.target_error, report = coarse.report)
```

The run with repeated remaps needs the full long run, so the construction has no result there.
