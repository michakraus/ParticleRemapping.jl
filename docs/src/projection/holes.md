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
test `test/projection/experiments/hole_filling.jl` checks each of them. A number of a fit that
ends on the stationarity test is quoted with its digits. A number of any other fit is quoted with
its digits where, with the targets moved by ``10^{-13}`` relative and on Julia 1.13 and 1.14, it
does not move at these digits and stays at least a tenth of a unit of its last digit from a
rounding boundary; otherwise it is quoted as the range over these runs. The rows "seeded" start
at ``\theta_0 = 10^{-3} \, \mathrm{randn}`` and give the range over the seeds 101 to 105
(item 2).

|  | ``x_h = \pi/2``, ``\partial_x W \neq 0``: empty disk | fill | objective | stop | ``x_h = \pi``, ``\partial_x W = 0``: empty disk | fill | objective | stop |
|:--|--:|--:|--:|:--|--:|--:|--:|:--|
| before | 2.75 | 0.00 | | | 2.75 | 0.00 | | |
| targets | 1.58 | 0.86 | | | 1.58 | 0.86 | | |
| ``\theta = 0``, reweighting only | 1.84 | 0.51 | 1 | | 1.88 | 0.47 | 1 | |
| weight-dependent | 1.43 | 0.90 | 0.258 | `MAXITERS` | 1.81 | 0.47 | 0.927 | `SUCCESS`, stationary (saddle) |
| weight-dependent, seeded | | | | | 1.35–1.45 | 0.97–1.06 | 0.285–0.333 | `MAXITERS` |
| control A | 1.67 | 0.66 | 0.834–0.835 | `MAXITERS` | 1.88 | 0.47 | 0.999 | `SUCCESS`, stationary |
| control A, seeded | | | | | 1.49–1.61 | 0.81–0.90 | 0.654–0.711 | `MAXITERS` |
| control B | 1.84 | 0.62 | 0.699 | `MAXITERS` | 1.88 | 0.55 | 0.972 | `SUCCESS`, stationary |
| control B, seeded | | | | | 1.51–1.98 | 0.62–0.78 | 0.749–0.937 | `MAXITERS` or `SUCCESS`, stationary |
| non-separable | 1.56 | 0.90 | 0.438 | `MAXITERS` | 1.51–1.52 | 0.93–0.94 | 0.520 | `LINESEARCH_FAILED` |
| control C | 1.88 | 0.54–0.55 | 0.667–0.668 | `LINESEARCH_FAILED` | 2.04 | 0.54–0.55 | 0.713–0.714 | `LINESEARCH_FAILED` |

After each map, the [charge-conserving reweighting](@ref projection-charge) on 32 periodic cubic
B-splines in ``x`` replaces ``w^*``; its reference is the grid charge of the markers before the
map, which carries the deficit of the hole. The second table gives the grid-charge residual of
``w^*`` relative to ``\max_i \rho_i``, the empty disk (disk) and the fill with the corrected
weights, and ``\max_a |\ln (w'_a / w^*_a)|``. Newton's method converges in at most 16 iterations,
and the residual after the correction is at most ``8.5 \cdot 10^{-14}``.

|  | ``x_h = \pi/2``: charge | disk | fill | ``\max\lvert\ln\rvert`` | ``x_h = \pi``: charge | disk | fill | ``\max\lvert\ln\rvert`` |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| ``\theta = 0`` | ``9.1 \cdot 10^{-3}`` | 1.84 | 0.51 | 0.014 | ``4.2 \cdot 10^{-3}`` | 1.91 | 0.47 | 0.008 |
| weight-dependent | ``8.9 \cdot 10^{-2}`` | 1.48 | 0.78 | 0.175–0.176 | ``4.2 \cdot 10^{-3}`` | 1.84 | 0.47 | 0.008 |
| weight-dependent, seeded | | | | | ``5.6 \cdot 10^{-2}``–``6.2 \cdot 10^{-2}`` | 1.32–1.41 | 0.74–0.78 | 0.153–0.168 |
| control A | ``2.7 \cdot 10^{-2}`` | 1.66 | 0.70 | 0.057 | ``4.2 \cdot 10^{-3}`` | 1.91 | 0.47 | 0.008 |
| control A, seeded | | | | | ``4.5 \cdot 10^{-2}``–``5.2 \cdot 10^{-2}`` | 1.50–1.74 | 0.50–0.55 | 0.102–0.116 |
| control B | ``2.8 \cdot 10^{-2}`` | 1.83 | 0.47 | 0.061–0.062 | ``4.2 \cdot 10^{-3}`` | 1.86 | 0.47 | 0.008 |
| control B, seeded | | | | | ``1.6 \cdot 10^{-2}``–``2.1 \cdot 10^{-2}`` | 1.45–1.96 | 0.46–0.78 | 0.047–0.071 |
| non-separable | ``9.1 \cdot 10^{-2}`` | 1.57 | 0.70 | 0.166 | ``6.4 \cdot 10^{-2}`` | 1.55–1.56 | 0.70 | 0.169 |
| control C | ``3.2 \cdot 10^{-2}`` | 1.85–1.86 | 0.54–0.55 | 0.071 | ``1.8 \cdot 10^{-2}`` | 2.05–2.06 | 0.54–0.55 | 0.053 |

In every run, ``|\Delta P| / \sum_a |p_a| < 10^{-15}`` and the total weight is exact.
The results are:

1. At ``x_h = \pi/2`` the weight-dependent map fills the hole: after the 40 iterations
   (`MAXITERS`, not at a minimum) its empty disk is ``1.43``, below the ``1.58`` of the targets,
   and the fill ``0.90``. Control B leaves the empty disk at the value of the reweighting alone.
2. At ``x_h = \pi``, a point with ``\partial_x W = 0``, the ``\pm v`` coupling of the [stage
   Jacobians](@ref projection-stage-jacobians) makes the fit of the weight-dependent map from
   ``\theta = 0`` lead to a mirror-symmetric saddle of its objective, and the map leaves the hole
   open (``1.81`` and ``0.47``). On a released Julia the fit ends there on the stationarity test
   (`SUCCESS`); on a Julia with a prerelease tag it ends there on the stationarity test or on
   `LINESEARCH_FAILED`, with the same numbers (`KNOWN_ISSUES.md`, K17). The
   kicks and drifts fill the hole once the fit breaks the symmetry: from the starts
   ``\theta_0 = 10^{-3} \, \mathrm{randn}`` of the seeds 101 to 105 (row "seeded", the range over
   the five seeds), the fit ends on `MAXITERS`, and after the 40 iterations, not at a minimum, the
   empty disk is ``1.35`` to ``1.45`` and the fill ``0.97`` to ``1.06``. No hole opens at
   ``-v``: the empty disk of the mirror ellipse around ``(\pi, -1.2)`` is ``0.76`` to ``0.87``.
   From ``\theta = 0`` the fits of control A and control B end on the stationarity test as well,
   with the empty disk of the reweighting alone (``1.88``); from the seeded starts their empty
   disks are ``1.49`` to ``1.61`` and ``1.51`` to ``1.98``, and their fills ``0.81`` to ``0.90``
   and ``0.62`` to ``0.78``. Control A ends on `MAXITERS` for every seed, control B on the
   stationarity test for three of the five seeds and on `MAXITERS` for the other two.
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
   ``1.51`` to ``1.52`` and the fill ``0.93`` to ``0.94``, where the kicks and drifts give
   ``1.81`` and ``0.47`` at the saddle and ``1.35`` to ``1.45`` and ``0.97`` to ``1.06`` from the
   seeded starts: at this point with ``\partial_x W = 0`` the ``\pm v`` coupling makes their fit
   from ``\theta = 0`` lead to a mirror-symmetric saddle, and they fill the hole once the fit
   breaks the symmetry. The non-separable stages do not have the ``\pm v`` coupling, and their
   fit from ``\theta = 0`` ends on `LINESEARCH_FAILED` with the hole filled. Control C does not
   fill the hole (``2.04`` and ``0.54`` to ``0.55``), so the filling comes from the dependence of
   ``\phi_s`` on ``w``, as the proposition on [non-separable stages](@ref projection-nonseparable)
   states. At ``x_h = \pi/2`` the non-separable map fills the hole as well (``1.56`` and
   ``0.90``, after the 40 iterations).
7. The non-separable map has more parameters than the kicks and drifts (3264 against 2304) and a
   coarser basis in ``\xi`` (16 against 32 splines). Its fixed-point iteration needs 50 iterations
   at the fitted ``\theta`` in both runs, the limit that the fit allows. In the run of the script,
   one application of the map takes 0.20 to 0.34 s, against 0.004 to 0.007 s for the kicks and
   drifts; the fit takes 411 s (40 iterations) and 67 s (25 iterations). The linearised step
   equations agree with forward-mode derivatives along random directions at the fitted
   ``\theta`` to a relative ``2.3 \cdot 10^{-13}`` ([`nonseparable_sensitivity_error`](@ref));
   without the terms of ``c`` and ``b`` they miss by at least ``2.9 \cdot 10^{-2}``.
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

```@example experiments
plot_target_flow(pi1.problem, pi1.targets; window = (π - 1.2, π + 1.2, 0, 2.4))
```

Density-interpolation flow near the hole at ``x_h = \pi``: each marker, coloured and sized by
weight, moves to its target (open circle) along the segment and carries markers into the hole
(row *targets*, fill ``0.86``; [`plot_target_flow`](@ref)).

```@example experiments
plot_hole_stages(remap_states(experiment_remap(pi2, :weight_dependent)...), pi2.problem.hole)
```

Stages of the kick–drift map at ``x_h = \pi/2``: markers near the hole before the map, after each
kick and each drift (velocities ``p / (m w)`` with the old weights), after the reweighting with
the target weights and after the charge-conserving reweighting ([`experiment_remap`](@ref),
[`remap_states`](@ref ParticleRemapping.remap_states), [`plot_hole_stages`](@ref)).

```@example experiments
plot_hole_stages(remap_states(experiment_remap(pi1, :non_separable)...), pi1.problem.hole;
    columns = 5)
```

Stages of the non-separable map at ``x_h = \pi``: markers near the hole before the map, after each
of the two stages (at the old weights), after the reweighting and after the charge-conserving
reweighting.

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

## [Target weights after a push](@id projection-push-test)

[`dlra_weights_experiment`](@ref) pushes the markers of the hole-filling experiment, with the
hole at ``x_h = \pi``, along the characteristics of
``\partial_t f + v \, \partial_x f + E(x) \, \partial_v f = 0`` with ``E = -0.15 \sin x`` up to
``T = 3``. The weights stay constant, the largest speed is ``3.95``, and the image of the hole is
sheared and contains no marker. The density ``f(T) = f_0 \circ \Phi_{-T}`` follows from backward
characteristics. The [BUG integrator](@ref projection-bug) solves the same equation on
``64 \times 256`` points of ``[0, 2\pi) \times [-6, 6]`` with the step ``0.01``, from rank four.
The target weights come from ``f(T)``, from ``f_r`` at four tolerances ``\vartheta``, from ``f_r``
at the fixed rank one, and from the lumped projection of the markers. Each set of weights enters
the fit of the non-separable map, and the charge-conserving reweighting follows against the grid
charge of the pushed markers. The first table gives the rank, the mass error of ``f_r``, and the
root mean square of ``w^* / (c f(T)) - 1`` in four velocity bands
([`weight_errors_by_band`](@ref)); no target weight is negative in any run. The script
`scripts/projection/dlra_weights.jl` (about 14 minutes) prints the tables of the push experiment
and checks every number against the value quoted here; the test
`test/projection/experiments/dlra_weights.jl` checks each of them.

| source | rank | mass error | ``\lvert v\rvert < 1`` | ``[1, 2)`` | ``[2, 3)`` | ``[3, 4)`` |
|:--|--:|--:|--:|--:|--:|--:|
| ``\vartheta = 10^{-3}`` | 4 | ``1.1 \cdot 10^{-3}`` | ``3.0 \cdot 10^{-2}`` | ``7.5 \cdot 10^{-2}`` | ``1.2 \cdot 10^{-1}`` | ``2.0 \cdot 10^{-1}`` |
| ``\vartheta = 10^{-4}`` | 7 | ``8.5 \cdot 10^{-5}`` | ``3.1 \cdot 10^{-3}`` | ``1.3 \cdot 10^{-2}`` | ``4.1 \cdot 10^{-2}`` | ``1.2 \cdot 10^{-1}`` |
| ``\vartheta = 10^{-5}`` | 10 | ``4.2 \cdot 10^{-6}`` | ``2.7 \cdot 10^{-4}`` | ``1.3 \cdot 10^{-3}`` | ``5.5 \cdot 10^{-3}`` | ``2.0 \cdot 10^{-2}`` |
| ``\vartheta = 10^{-6}`` | 13 | ``1.9 \cdot 10^{-7}`` | ``2.3 \cdot 10^{-5}`` | ``8.7 \cdot 10^{-5}`` | ``5.8 \cdot 10^{-4}`` | ``2.9 \cdot 10^{-3}`` |
| rank one | 1 | below ``10^{-12}`` | ``1.5 \cdot 10^{-1}`` | ``1.8 \cdot 10^{-1}`` | ``1.4 \cdot 10^{-1}`` | ``1.8 \cdot 10^{-1}`` |
| projection | | | ``4.1 \cdot 10^{-2}`` | ``1.1 \cdot 10^{-1}`` | ``8.4 \cdot 10^{-2}`` | ``2.0 \cdot 10^{-1}`` |

The second table gives, after the map and the reweighting with ``w^*``, the largest empty disk
over the whole domain ([`largest_empty_disk`](@ref)), the fill of the image of the hole, and the
root mean square deviation of the marker density from its mean ([`density_discrepancy`](@ref));
then the same with the charge-conserving weights, the grid-charge residual of ``w^*``, and
``\max_a |\ln (w'_a / w^*_a)|``. The pushed markers before the map have the empty disk ``1.38``,
the fill ``0`` and the density deviation ``0.099``. The fit is GeometricSolvers' Gauss–Newton
method with the `DogLeg` globalisation, at most 40 iterations, from ``\theta = 0``; the column
stop gives its return code: `MAXITERS` after the 40 iterations, and `LINESEARCH_FAILED` where
`DogLeg` finds no acceptable step. No fit ends on the stationarity test, so each row gives the
values where the fit stops, not at a minimum. A number is quoted with its digits where, with the
targets moved by ``10^{-13}`` relative and on Julia 1.13 and 1.14, it does not move at these
digits and stays at least a tenth of a unit of its last digit from a rounding boundary; otherwise
it is quoted as the range over these runs.

|  | with ``w^*``: empty disk | fill | density | charge-conserving: empty disk | fill | charge | ``\max\lvert\ln\rvert`` | stop |
|:--|--:|--:|--:|--:|--:|--:|--:|:--|
| ``f(T)`` | 1.05 | 0.82 | 0.113 | 1.066 | 0.74 | ``2.4 \cdot 10^{-2}`` | 0.035–0.036 | `MAXITERS` |
| ``\vartheta = 10^{-3}`` | 1.43 | 1.13 | 0.308–0.309 | 1.52–1.53 | 1.21 | ``3.3 \cdot 10^{-2}`` | 0.039 | `LINESEARCH_FAILED` |
| ``\vartheta = 10^{-4}`` | 1.23 | 1.01 | 0.317 | 1.27 | 0.936 | ``1.7 \cdot 10^{-2}`` | 0.027 | `LINESEARCH_FAILED` |
| ``\vartheta = 10^{-5}`` | 1.05 | 0.86 | 0.114 | 1.064 | 0.74 | ``2.4 \cdot 10^{-2}`` | 0.032 | `MAXITERS` |
| ``\vartheta = 10^{-6}`` | 1.05 | 0.82 | 0.113 | 1.066 | 0.78 | ``2.5 \cdot 10^{-2}`` | 0.035 | `LINESEARCH_FAILED` |
| rank one | 1.06 | 0.975 | 0.332 | 2.74 | 0.90 | ``1.75 \cdot 10^{-1}`` | 0.258 | `LINESEARCH_FAILED` |
| projection | 2.96–2.97 | 0.35 | 0.519 | 3.05 | 0.31 | ``7.9 \cdot 10^{-2}`` | 0.0986 | `LINESEARCH_FAILED` |

In every run ``|\Delta P| / \sum_a |p_a|`` is below ``10^{-12}``, and so is the grid-charge
residual after the correction. The results are:

1. The weight error of ``f_r`` falls with ``\vartheta`` in every band, at ranks between 4 and 13:
   by about one order of magnitude per order of ``\vartheta`` for ``|v| < 1``, by less in the
   tails between ``\vartheta = 10^{-3}`` and ``10^{-4}``. At every tolerance it grows with
   ``|v|``, by a factor of about 100 from ``|v| < 1`` to ``|v| \in [3, 4)`` at
   ``\vartheta = 10^{-6}``: the truncation bounds an absolute error, and ``f`` is small in the
   tails.
2. At ``\vartheta \leq 10^{-5}`` the empty disks and the density deviation of the remap agree
   with those with the weights of ``f(T)`` to within two units of the last digit of the table,
   and the fills differ by up to ``0.04``. At ``\vartheta = 10^{-4}`` and ``10^{-3}`` the
   velocities inherit the weight errors, and the density deviation nearly triples, although the
   error in the core is ``3 \cdot 10^{-3}`` at ``\vartheta = 10^{-4}``; the error in the tails
   dominates.
3. The lumped projection gives weight errors of the size of those at ``\vartheta = 10^{-3}`` and
   leaves the largest empty disk at ``2.96``–``2.97``.
4. The charge-conserving reweighting changes the weights by at most ``4 \%`` for every source
   except rank one and the lumped projection, which it changes by up to ``10 \%``. With the
   weights of ``f(T)`` it lowers the fill of the image of the hole from ``0.82`` to ``0.74``. For
   rank one it changes the weights by up to ``29 \%``, and the largest empty disk grows from
   ``1.06`` to ``2.74``.
5. In the run of the script, one BUG solve to ``T = 3`` takes about ``0.1`` s at every tolerance
   once compiled.

The field is weak and fixed, and the dimension is one, so the ranks are small. A self-consistent
field and ``d = 3`` are not tested.

The figure below runs the experiment with the fit of the map for the weights of ``f_r`` at
``\vartheta = 10^{-6}`` only (`fit = ("ϑ = 1e-6",)`), a reduced version of the script run; the
weight errors are those of the full run.

```@example experiments
using ParticleRemapping, CairoMakie, Random
Random.seed!(20261119)
push = dlra_weights_experiment(; fit = ("ϑ = 1e-6",))
plot_dlra_weights(push)
```

Target weights after a push. (a) Pushed markers and the image of the hole. (b) The non-separable
map with the weights of ``f_r`` at ``\vartheta = 10^{-6}``, after the charge-conserving
reweighting; (a) and (b) are coloured and sized by weight on one scale. (c) ``f_r - f(T)`` at ``\vartheta = 10^{-4}``. (d) Numerical rank of ``f_r`` over
time. (e) Relative error of the target weights against ``|v^*|``.

```@example experiments
animate_dlra("dlra.mp4", dlra_push_frames(push), push.grid)
nothing # hide
```

```@raw html
<video src="./dlra.mp4" controls loop muted playsinline style="width:100%"></video>
```

The low-rank density at ``\vartheta = 10^{-4}`` and the markers along the push
([`dlra_push_frames`](@ref)).
