# [Hole filling](@id sliced-holes)

```@meta
CurrentModule = ParticleRemapping.SlicedRedistribution
```

This page shows the sliced lifts alone, without the squeeze, on the markers of the standard
checkerboard, and the control with one slice for each lift. The full remap of the same markers is
part of the [standard experiments](@ref sliced-experiments).

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
