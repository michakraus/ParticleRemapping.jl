# [Hole filling](@id projection-hole-filling)

```@meta
CurrentModule = ParticleRemapping.Projection
```

This page gives the hole-filling experiment of the chapter. The script
`scripts/projection/hole_filling.jl` (about 31 minutes) prints the tables of the experiment and
checks every number against the value quoted here; the test
`test/projection/experiments/hole_filling.jl` (about 18 minutes) checks each of them.

[`hole_filling_experiment`](@ref) runs the whole construction for ``d = 1``. The markers sit on a
``48 \times 48`` lattice in ``[0, 2\pi) \times [-4, 4]``, with weights ``w = c f`` for
``f = (1 + 0.3 \cos x) \, (2\pi)^{-1/2} e^{-v^2/2}``. An elliptic hole with semi-axes ``0.45`` and
``0.40`` around ``(x_h, 1.2)`` removes 24 markers, which leaves ``N = 2280``. The projection uses
``32 \times 33`` hat functions, the density-interpolation flow 40 steps of the classical
Runge–Kutta method. The map has two stages with 32 periodic splines in ``\xi``, 40 splines in
``v`` and 16 splines in ``\ln w``. The fit is GeometricSolvers' Gauss–Newton method with the
`DogLeg` globalisation, with ``\beta = 10^{-3}``, at most 40 iterations and the stationarity test
of GeometricSolvers; it starts at ``\theta = 0``, and each table gives its return code (stop):
`SUCCESS, stationary` where the stationarity test ends it, `MAXITERS` after the 40 iterations,
and `LINESEARCH_FAILED` where `DogLeg` finds no acceptable step. A row that ends on
`MAXITERS` gives the values after the 40 iterations, not at a minimum. The target weights are
``c f`` at the targets, so that the run measures the map and not the projection. Five
parametrisations are compared ([`hole_remap_method`](@ref)):

- the weight-dependent one of the [parametrisation](@ref projection-parametrisation);
- control A, which does not depend on the weight in ``(x,p)``;
- control B, ``\kappa = w \sum_i \lambda_i B_i(\xi)`` and ``\eta = w \sum_j \mu_j \hat{B}_j(v)``,
  which is one area-preserving map of ``(x,v)`` by the [stage
  Jacobians](@ref projection-stage-jacobians);
- the non-separable map, two [non-separable stages](@ref projection-nonseparable) on 16 periodic
  splines in ``\xi``, 17 splines in ``v`` and 6 splines in ``\ln w``, 3264 parameters;
- control C, the same stages with ``\phi_s`` evaluated at one value of ``\ln w`` for all markers,
  which is one area-preserving map of ``(x,v)`` by item 4 of the proposition on non-separable
  stages.

For the non-separable stages, the fit rejects every ``\theta`` for which the fixed-point
iteration of the velocity form needs more than 50 iterations, which keeps its contraction rate
below about ``0.55``. The Jacobian of the fit comes from the linearised step equations. The table
gives the largest empty disk near the hole in units of the lattice spacing (a perfect lattice has
``0.71``; [`local_empty_disk`](@ref)), the number of markers in the hole relative to an
equidistribution ([`hole_fill`](@ref)), and the objective after the fit relative to
``\theta = 0``. All numbers of this section come from `scripts/projection/hole_filling.jl`; the
test `test/projection/experiments/hole_filling.jl` checks each of them. A number is quoted with
its digits where it does not move when the targets move by ``10^{-13}`` relative. The numbers of
control A at ``x_h = \pi`` move, so its row gives the range over the run and two runs with moved
targets; the row "seeded" gives the range over the starts of the seeds 101 to 105 (item 2).

|  | ``x_h = \pi/2``, ``\partial_x W \neq 0``: empty disk | fill | objective | stop | ``x_h = \pi``, ``\partial_x W = 0``: empty disk | fill | objective | stop |
|:--|--:|--:|--:|:--|--:|--:|--:|:--|
| before | 2.75 | 0.00 | | | 2.75 | 0.00 | | |
| targets | 1.58 | 0.86 | | | 1.58 | 0.86 | | |
| ``\theta = 0``, reweighting only | 1.84 | 0.51 | 1 | | 1.88 | 0.47 | 1 | |
| weight-dependent | 1.43 | 0.90 | 0.258 | `MAXITERS` | 1.81 | 0.47 | 0.927 | `SUCCESS`, stationary (saddle) |
| weight-dependent, seeded | | | | | 1.35–1.45 | 0.97–1.06 | 0.285–0.333 | `MAXITERS` |
| control A | 1.67 | 0.66 | 0.834 | `MAXITERS` | 1.64–1.66 | 0.81–0.86 | 0.697–0.704 | `MAXITERS` |
| control B | 1.84 | 0.62 | 0.699 | `MAXITERS` | 1.88 | 0.55 | 0.972 | `SUCCESS`, stationary |
| non-separable | 1.56 | 0.90 | 0.438 | `MAXITERS` | 1.52 | 0.94 | 0.520 | `LINESEARCH_FAILED` |
| control C | 1.88 | 0.55 | 0.668 | `LINESEARCH_FAILED` | 2.04 | 0.55 | 0.713 | `LINESEARCH_FAILED` |

After each map, the [charge-conserving reweighting](@ref projection-charge) on 32 periodic cubic
B-splines in ``x`` replaces ``w^*``; its reference is the grid charge of the markers before the
map, which carries the deficit of the hole. The second table gives the grid-charge residual of
``w^*`` relative to ``\max_i \rho_i``, the empty disk (disk) and the fill with the corrected
weights, and ``\max_a |\ln (w'_a / w^*_a)|``. Newton's method converges in at most 16 iterations,
and the residual after the correction is at most ``8.5 \cdot 10^{-14}``.

|  | ``x_h = \pi/2``: charge | disk | fill | ``\max\lvert\ln\rvert`` | ``x_h = \pi``: charge | disk | fill | ``\max\lvert\ln\rvert`` |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| ``\theta = 0`` | ``9.1 \cdot 10^{-3}`` | 1.84 | 0.51 | 0.014 | ``4.2 \cdot 10^{-3}`` | 1.91 | 0.47 | 0.008 |
| weight-dependent | ``8.9 \cdot 10^{-2}`` | 1.48 | 0.78 | 0.175 | ``4.2 \cdot 10^{-3}`` | 1.84 | 0.47 | 0.008 |
| weight-dependent, seeded | | | | | ``5.6 \cdot 10^{-2}``–``6.2 \cdot 10^{-2}`` | 1.32–1.41 | 0.74–0.78 | 0.153–0.168 |
| control A | ``2.7 \cdot 10^{-2}`` | 1.66 | 0.70 | 0.057 | ``4.0 \cdot 10^{-2}``–``4.4 \cdot 10^{-2}`` | 1.67–1.79 | 0.55 | 0.091–0.096 |
| control B | ``2.8 \cdot 10^{-2}`` | 1.83 | 0.47 | 0.061 | ``4.2 \cdot 10^{-3}`` | 1.86 | 0.47 | 0.008 |
| non-separable | ``9.1 \cdot 10^{-2}`` | 1.57 | 0.70 | 0.166 | ``6.4 \cdot 10^{-2}`` | 1.55 | 0.70 | 0.169 |
| control C | ``3.2 \cdot 10^{-2}`` | 1.86 | 0.55 | 0.071 | ``1.8 \cdot 10^{-2}`` | 2.06 | 0.55 | 0.053 |

In every run, ``|\Delta P| / \sum_a |p_a| \leq 1.1 \cdot 10^{-17}`` and the total weight is exact.
The results are:

1. At ``x_h = \pi/2`` the weight-dependent map fills the hole: after the 40 iterations
   (`MAXITERS`, not at a minimum) its empty disk is ``1.43``, below the ``1.58`` of the targets,
   and the fill ``0.90``. Control B leaves the empty disk at the value of the reweighting alone.
2. At ``x_h = \pi``, a point with ``\partial_x W = 0``, the ``\pm v`` coupling of the [stage
   Jacobians](@ref projection-stage-jacobians) makes the fit of the weight-dependent map from
   ``\theta = 0`` lead to a mirror-symmetric saddle of its objective: the fit ends on the
   stationarity test within 5 iterations (`SUCCESS`), and the map leaves the hole open (``1.81``
   and ``0.47``). The kicks and drifts fill the hole once the fit breaks the symmetry: from the
   starts ``\theta_0 = 10^{-3} \, \mathrm{randn}`` of the seeds 101 to 105 (row "seeded", the range
   over the five seeds), the fit ends on `MAXITERS`, and after the 40 iterations, not at a
   minimum, the empty disk is ``1.35`` to ``1.45`` and the fill ``0.97`` to ``1.06``. No hole
   opens at ``-v``: the empty disk of the mirror ellipse around ``(\pi, -1.2)`` is ``0.76`` to
   ``0.87``.
3. The reweighting alone moves markers towards the targets, because ``v' = v \, w / w^*``.
4. Control A shrinks the empty disk at ``x_h = \pi/2``. It keeps the area of the hole in ``(x,p)``
   by the corollary [holes persist](@ref holes-persist), but it deforms the hole, and
   in ``(x,v)`` it acts on each weight differently through ``v = p / (m w)``.
5. With the lumped projection on the same grid, the relative error of the target weights at
   ``x_h = \pi`` has the root mean square ``0.10``, ``0.17``, ``0.13`` and ``0.21`` in the velocity
   bands ``[0,1)``, ``[1,2)``, ``[2,3)`` and ``[3,4)``, and maxima up to ``0.96``
   ([`projected_weight_errors`](@ref)); at ``x_h = \pi/2`` the root mean squares are ``0.10``,
   ``0.18``, ``0.13`` and ``0.22``. By the [reweighting](@ref projection-final-reweighting), the
   velocities inherit errors of the same size.
6. At ``x_h = \pi`` the non-separable map fills the hole from ``\theta = 0``, with the empty disk
   ``1.52`` and the fill ``0.94``, where the kicks and drifts give ``1.81`` and ``0.47`` at the
   saddle and ``1.35`` to ``1.45`` and ``0.97`` to ``1.06`` from the seeded starts. The
   non-separable stages do not have the ``\pm v`` coupling, and their fit from ``\theta = 0``
   ends on `LINESEARCH_FAILED` with the hole filled. Control C does
   not fill the hole (``2.04`` and ``0.55``), so the filling comes from the dependence of
   ``\phi_s`` on ``w``, as the proposition on [non-separable stages](@ref projection-nonseparable)
   states. At ``x_h = \pi/2`` the non-separable map fills the hole as well (``1.56`` and
   ``0.90``, after the 40 iterations).
7. The non-separable map has more parameters than the kicks and drifts (3264 against 2304) and a
   coarser basis in ``\xi`` (16 against 32 splines). Its fixed-point iteration needs 50 iterations
   at the fitted ``\theta`` in both runs, the limit that the fit allows. In the run of the script,
   one application of the map takes 0.21 to 0.34 s, against 0.004 to 0.007 s for the kicks and
   drifts; the fit takes 411 s (40 iterations) and 67 s (25 iterations). The linearised step equations agree with forward-mode
   derivatives along random directions at the fitted ``\theta`` to a relative
   ``2.3 \cdot 10^{-13}`` ([`nonseparable_sensitivity_error`](@ref)); without the terms of ``c``
   and ``b`` they miss by at least ``2.9 \cdot 10^{-2}``.
8. The target weights ``c f`` describe the density without the hole, and the markers before the
   map carry its charge deficit. A map that fills the hole with these weights adds charge in the
   columns of the hole: up to ``9 \%`` of ``\max_i \rho_i`` for the weight-dependent and the
   non-separable map. The charge-conserving reweighting removes it and changes the weights by up
   to ``16 \%``; the velocities ``p / (m w')`` change by the inverse factor, and the fill drops
   to between ``0.70`` and ``0.78``. For these two maps the empty disk changes by less than
   ``0.06`` with the correction. At ``\theta = 0`` the correction changes the weights by
   at most ``1.4 \%``.

The figures below run the experiment with the weight-dependent map, the non-separable map and
control C (the other variants do not enter the figures).

```@example experiments
using ParticleRemapping, CairoMakie, Random
Random.seed!(20261119)
variants = (:weight_dependent, :non_separable, :control_C)
pi2 = hole_filling_experiment(; hole_x = π / 2, variants)
pi1 = hole_filling_experiment(; hole_x = π, variants)
plot_hole_filling(pi2)
```

Hole filling at ``x_h = \pi/2``: markers near the hole (dashed outline) in ``(x, v)``, coloured
and sized by weight, with ``v = p / (m w)`` for the weights of each panel. Panels (g) and (h) are
(e) and (f) after the charge-conserving reweighting.

```@example experiments
plot_hole_filling(pi1)
```

Hole filling at ``x_h = \pi``, as above.

```@example experiments
plot_charge_correction(pi2, pi1; labels = ["x_h = π/2", "x_h = π"])
```

Charge-conserving reweighting after the non-separable map. (a) Grid charge of the target weights
(lines) and of the corrected weights (dots) at the positions after the map, against the grid
charge before the map. (b) The correction ``\ln (w' / w^*)`` per marker.

The stages of the map, each with its parameters scaled from zero to their fitted values, then the
reweighting and the charge-conserving reweighting: the non-separable map at ``x_h = \pi`` and the
kicks and drifts at ``x_h = \pi/2`` ([`experiment_remap`](@ref),
[`remap_frames`](@ref ParticleRemapping.remap_frames)).

```@example experiments
for (e, v, name) in ((pi1, :non_separable, "remap-non-separable-pi.mp4"),
    (pi2, :weight_dependent, "remap-kick-drift-pi2.mp4"))
    animate_stages(name, remap_frames(experiment_remap(e, v)...); period = (0, 2π),
        decorate! = (ax, _) -> plot_hole_outline!(ax, e.problem.hole))
end
nothing # hide
```

```@raw html
<video src="./remap-non-separable-pi.mp4" controls loop muted playsinline style="width:100%"></video>
<video src="./remap-kick-drift-pi2.mp4" controls loop muted playsinline style="width:100%"></video>
```
