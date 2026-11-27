# [Hole filling](@id sliced-holes)

```@meta
CurrentModule = ParticleRemapping.SlicedRedistribution
```

This page shows the sliced lifts alone, without the squeeze, on the markers of the standard
checkerboard with the one-slice control, and a 2D-2V hole between the axes that the per-axis lifts
do not see. The full remap of the checkerboard is part of the
[standard experiments](@ref sliced-experiments).

## [The lifts on the checkerboard](@id sliced-checkerboard)

[`checkerboard_experiment`](@ref): 400 markers of weight ``1/400`` (``m = 1``), uniform in the
quadrants ``Q_1 = [0, 1)^2`` and ``Q_3 = [-1, 0)^2``; ``Q_2`` and ``Q_4`` are empty, and both
marginals are uniform on ``[-1, 1]``. Only the lifts act, the weights do not change.

| map | ``Q_1`` | ``Q_2`` | ``Q_3`` | ``Q_4`` | largest ``\lvert x \rvert`` | largest ``\lvert v \rvert`` |
|:--|--:|--:|--:|--:|--:|--:|
| before | 200 | 0 | 200 | 0 | 1.00 | 1.00 |
| sliced lifts, ``B_v = 10``, ``B_x = 20`` | 98 | 101 | 112 | 89 | 1.28 | 1.04 |
| one slice for each lift (control) | 200 | 0 | 200 | 0 | 6.77 | 1.79 |

The thresholds are ``Q_2 + Q_4 < 20`` for the control and ``65 \le Q_k \le 135`` for the sliced
lifts. Within each velocity slice the positions lie on one side of ``x = 0``; ``\Psi_q`` spreads
them over ``[-1, 1]`` and so fills ``Q_2`` and ``Q_4``. The one-slice lifts see the uniform
marginals and leave the markers in their quadrants. They also stretch the conjugate variables:
without anchors the cotangent part ``(\nabla^2 S)^{-1}`` amplifies them by up to
``1 / (1 - \beta_{\mathrm{intra}}) = 6.7`` where the Laplacian of the potential is weak. The
markers have equal weights, so the figure colours each marker by the quadrant it comes from.

```@example sliced-experiments
using ParticleRemapping, CairoMakie
set_theme!(remap_theme())
plot_checkerboard_filling(checkerboard_experiment())
```

## [A hole between the axes](@id sliced-hole-between-axes)

Each lift reads its slices and targets from the two-dimensional marginal ``(x_\ell, v_\ell)``
only, so the lifts cannot see a hole in the joint distribution of different axes, such as in
``(x_1, x_2)``, ``(v_1, v_2)`` or ``(x_1, v_2)``. The analogue of the
[checkerboard](@ref sliced-checkerboard) one dimension up shows it: 360 equal-weight markers in
2D-2V fill the quadrants ``Q_1`` and ``Q_3`` of the ``(x_1, x_2)`` plane, with velocities
uniform on ``[-1, 1]^2`` and independent of the positions, so every marginal
``(x_\ell, v_\ell)`` is even. The checkerboard of the same markers in ``(x_1, v_1)`` is the
control:

```@example sliced-hole-between-axes
using ParticleRemapping
using Random: MersenneTwister
rng = MersenneTwister(7)
N, h = 360, 180
method = SlicedQuantileRemap(; B_v = 6, B_x = 36)
x = hcat(π .* rand(rng, 2, h), -π .* rand(rng, 2, h))
v = 2 .* rand(rng, 2, N) .- 1
r = remap(method, Markers(x, v, fill(1 / N, N)))
x1 = vcat(π .* rand(rng, h), -π .* rand(rng, h))
v1 = vcat(rand(rng, h), -rand(rng, h))
r1 = remap(SlicedQuantileRemap(; B_v = 6, B_x = 36, Vmax = 1.0),
    Markers(x1, v1, fill(1 / N, N)))
x̄ = periodic_wrap.(r.markers.x, π)
(; x1_v1_before = phase_quadrant_counts(x[1, :], v[1, :]),
    x1_x2_before = phase_quadrant_counts(x[1, :], x[2, :]),
    x1_x2_after_lifts = phase_quadrant_counts(r.info.q̂[1, :], r.info.q̂[2, :]),
    x1_x2_after_remap = phase_quadrant_counts(x̄[1, :], x̄[2, :]),
    control_x1_v1_after_lifts = phase_quadrant_counts(r1.info.q̂, r1.info.p̂))
```

The per-axis lifts leave most markers in the quadrants ``Q_1`` and ``Q_3`` of ``(x_1, x_2)``,
while the same lifts fill the checkerboard in ``(x_1, v_1)``. After the squeeze the counts are
more even, but not by the lifts: the multipliers come from the density ratio in all four
coordinates, which is small in the empty quadrants, and the squeeze moves the markers of those
slots far from their anchors and across the periodic box (K6). `test/sliced/peraxis.jl` checks
the counts after the lifts.
