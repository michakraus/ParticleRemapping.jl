# [Numerical experiments](@id projection-experiments)

```@meta
CurrentModule = ParticleRemapping.Projection
```

Every statement of this chapter has a test under `test/projection/`, and each test runs the case
that must fail. Derivatives are computed by forward-mode automatic differentiation, so the
tolerances are set by rounding error. The script `scripts/projection/checks.jl` runs the tests of
the structural statements and prints the measured numbers of the first table; the scripts
`scripts/projection/hole_filling.jl` and `scripts/projection/dlra_weights.jl` (about 3 minutes
each) print the tables of the two experiments and check every number against the
value quoted here.

## The structural statements

The first three rows check statements of
[What a symplectic redistribution can do](@ref obstructions) with the functions of this package;
the other rows check the statements of this chapter.

| statement | result | case that must fail |
|:--|:--|:--|
| [compressibility](@ref weight-dependent-hamiltonians), ``d = 1, 2, 3`` | ``\lvert\operatorname{div} u - \{W, \partial_w \Psi\}\rvert \leq 1.1 \cdot 10^{-15}`` at random points, for random ``\Psi`` and ``W`` | opposite sign: ``\geq 7.8 \cdot 10^{-2}``; frozen weight: ``\geq 3.9 \cdot 10^{-2}`` |
| [fold](@ref velocity-momentum-fold), ``d = 1, 2, 3`` | relative error of the determinant ``\leq 1.7 \cdot 10^{-15}`` for random positive ``f``; Maxwellian: sign change at ``\lvert v\rvert = v_{\mathrm{th}}``, two roots at ``p = p_{\max}/2``, none at ``p = 3 p_{\max}/2`` | ``f^d`` in place of ``f^{d-1}``: ``\geq 2.4 \cdot 10^{-3}`` (``d \geq 2``) |
| [mean-field Hamiltonian](@ref one-map-for-all-markers), [holes persist](@ref holes-persist) | symplectic Euler and implicit midpoint on 274 markers: ``\lvert\psi(z_a) - z'_a\rvert \leq 6.8 \cdot 10^{-16}`` while markers move up to ``0.39``; ``\lvert\det D\psi - 1\rvert \leq 2.2 \cdot 10^{-16}``; hole area ``0.384845`` unchanged, no marker in the image of the hole | two coincident markers with weights ``0.5`` and ``2`` under a weight-dependent Hamiltonian: images differ by ``0.14`` |
| [kick–drift map](@ref projection-hamiltonian), ``d = 1, 2`` | ``\max \lvert DR_\theta^\top \mathbb{J} DR_\theta - \mathbb{J}\rvert \leq 2.9 \cdot 10^{-15}``; ``\lvert\Delta P\rvert \leq 8.9 \cdot 10^{-16}``; boost spread ``2.8 \cdot 10^{-17}``; kick and drift agree with symplectic Euler of ``H^K + H^D`` to ``3.1 \cdot 10^{-15}`` | naive kick: defect ``0.86`` and ``1.6`` with ``\lvert\Delta P\rvert \leq 8.9 \cdot 10^{-16}``; ``\theta = \theta(z)`` ([parameters that depend on the state](@ref state-dependent-parameters)): defect ``0.17`` and ``0.16``, against ``8.6 \cdot 10^{-15}`` with ``\theta`` frozen |
| [stage Jacobians](@ref projection-stage-jacobians) | error ``\leq 2.2 \cdot 10^{-16}``; at ``\pm v``: ``\det DK - 1 = \mp 0.047`` | none; without ``\ell``-dependence ``\det DK = 1`` exactly |
| [non-separable stages](@ref projection-nonseparable), ``d = 1, 2`` | the velocity form agrees with the implicit midpoint rule for ``H^N_s`` to ``4.4 \cdot 10^{-16}``; ``\max \lvert DR^\top \mathbb{J} DR - \mathbb{J}\rvert \leq 1.1 \cdot 10^{-15}`` for two stages; ``\lvert\Delta P\rvert \leq 4.4 \cdot 10^{-16}``; error of the determinant ``\leq 4.4 \cdot 10^{-15}``, of the divergence ``\leq 8.9 \cdot 10^{-16}``; at ``x = \pi``, ``\pm v``: ``\det DN - 1 = \mp 0.120`` for a separable, ``-0.144`` at both for a non-separable ``\phi_s`` | trapezoidal rule: defect ``7.5 \cdot 10^{-4}`` and ``3.6 \cdot 10^{-4}`` with ``\Delta P = 0``; no centroid shift: ``\lvert\Delta P\rvert = 0.29`` and ``1.7``; first-order determinant: ``\geq 0.315``; sign of ``\tfrac{1}{2} \mathbb{J} K`` reversed: ``\geq 0.475`` |
| [charge-conserving reweighting](@ref projection-charge), ``d = 1``, 400 markers, 16 shape functions | after a kick and a drift, Newton converges in 5 iterations; grid charge to ``7.1 \cdot 10^{-16}``; total weight and ``P`` exact; ``\min w' > 0``; gradient of the relative entropy orthogonal to the null space of ``A`` to ``9.1 \cdot 10^{-16}``, and 50 feasible perturbations increase it | ``w^*`` alone: charge error ``0.29``; additive correction: exact charge, ``\min w' = -0.015``; markers compressed to ``[0, \pi)``: no solution, reported |
| [BUG integrator](@ref projection-bug), ``64 \times 256`` grid, ``T = 3`` | free streaming from rank four: numerical rank 3, error against the full grid ``7.6 \cdot 10^{-7}``, ``1.9 \cdot 10^{-7}`` at half the step; full grid against characteristics for ``E = -0.5 \sin x``: ``2.7 \cdot 10^{-5}`` (``3.6 \cdot 10^{-4}`` on half the velocity points); ``\vartheta = 10^{-4}, 10^{-6}``: rank ``13, 31``, error ``3.6 \cdot 10^{-3}, 2.8 \cdot 10^{-5}`` below ``n_{\mathrm{steps}} \vartheta``; mass error ``6.8 \cdot 10^{-5}, 3.2 \cdot 10^{-7}`` | free streaming from rank one: stays at rank one, error ``0.28``; fixed rank one with the field: error ``0.16`` |

The [verification](@ref projection-verification) lists the test file and the script of every
row. All checks pass. The figures of three of the statements — the image of a hole under one
weight-independent step, the fold of the implicit weight for a Maxwellian, and the ``\pm v``
coupling at ``x = \pi`` — are in [one map for all markers](@ref one-map-for-all-markers),
[the velocity–momentum fold](@ref velocity-momentum-fold) and
[the non-separable stages](@ref projection-nonseparable).

## [Hole filling](@id projection-hole-filling)

[`hole_filling_experiment`](@ref) runs the whole construction for ``d = 1``. The markers sit on a
``48 \times 48`` lattice in ``[0, 2\pi) \times [-4, 4]``, with weights ``w = c f`` for
``f = (1 + 0.3 \cos x) \, (2\pi)^{-1/2} e^{-v^2/2}``. An elliptic hole with semi-axes ``0.45`` and
``0.40`` around ``(x_h, 1.2)`` removes 24 markers, which leaves ``N = 2280``. The projection uses
``32 \times 33`` hat functions, the density-interpolation flow 40 steps of the classical
Runge–Kutta method. The map has two stages with 32 periodic splines in ``\xi``, 40 splines in
``v`` and 16 splines in ``\ln w``; the fit is Levenberg–Marquardt with ``\beta = 10^{-3}`` and at
most 40 iterations. The target weights are ``c f`` at the targets, so that the run measures the
map and not the projection. Five parametrisations are compared ([`hole_remap_method`](@ref)):

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
test `test/projection/experiments/hole_filling.jl` checks each of them.

|  | ``x_h = \pi/2``, ``\partial_x W \neq 0``: empty disk | fill | objective | ``x_h = \pi``, ``\partial_x W = 0``: empty disk | fill | objective |
|:--|--:|--:|--:|--:|--:|--:|
| before | 2.75 | 0.00 | | 2.75 | 0.00 | |
| targets | 1.58 | 0.86 | | 1.58 | 0.86 | |
| ``\theta = 0``, reweighting only | 1.84 | 0.51 | 1 | 1.88 | 0.47 | 1 |
| weight-dependent | 1.58 | 0.94 | 0.264 | 1.81 | 0.47 | 0.927 |
| control A | 1.65 | 0.82 | 0.795 | 1.88 | 0.47 | 0.999 |
| control B | 1.84 | 0.62 | 0.699 | 1.88 | 0.55 | 0.972 |
| non-separable | 1.54 | 0.94 | 0.439 | 1.50 | 0.94 | 0.521 |
| control C | 1.80 | 0.55 | 0.663 | 2.00 | 0.55 | 0.716 |

After each map, the [charge-conserving reweighting](@ref projection-charge) on 32 periodic cubic
B-splines in ``x`` replaces ``w^*``; its reference is the grid charge of the markers before the
map, which carries the deficit of the hole. The second table gives the grid-charge residual of
``w^*`` relative to ``\max_i \rho_i``, the empty disk (disk) and the fill with the corrected
weights, and ``\max_a |\ln (w'_a / w^*_a)|``. Newton's method converges in at most 16 iterations,
and the residual after the correction is at most ``8.5 \cdot 10^{-14}``.

|  | ``x_h = \pi/2``: charge | disk | fill | ``\max\lvert\ln\rvert`` | ``x_h = \pi``: charge | disk | fill | ``\max\lvert\ln\rvert`` |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| ``\theta = 0`` | ``9.1 \cdot 10^{-3}`` | 1.84 | 0.51 | 0.014 | ``4.2 \cdot 10^{-3}`` | 1.91 | 0.47 | 0.008 |
| weight-dependent | ``9.0 \cdot 10^{-2}`` | 1.56 | 0.70 | 0.173 | ``4.2 \cdot 10^{-3}`` | 1.84 | 0.47 | 0.008 |
| control A | ``5.4 \cdot 10^{-2}`` | 1.90 | 0.78 | 0.142 | ``4.2 \cdot 10^{-3}`` | 1.91 | 0.47 | 0.008 |
| control B | ``2.8 \cdot 10^{-2}`` | 1.83 | 0.47 | 0.061 | ``4.2 \cdot 10^{-3}`` | 1.86 | 0.47 | 0.008 |
| non-separable | ``9.2 \cdot 10^{-2}`` | 1.56 | 0.70 | 0.165 | ``6.3 \cdot 10^{-2}`` | 1.53 | 0.70 | 0.165 |
| control C | ``3.1 \cdot 10^{-2}`` | 1.83 | 0.55 | 0.069 | ``1.6 \cdot 10^{-2}`` | 2.03 | 0.55 | 0.050 |

In every run, ``|\Delta P| / \sum_a |p_a| \leq 1.1 \cdot 10^{-17}`` and the total weight is exact.
The results are:

1. At ``x_h = \pi/2`` the weight-dependent map reaches the empty disk of the targets and fills
   the hole. Control B leaves the empty disk at the value of the reweighting alone.
2. At ``x_h = \pi`` the weight-dependent map does not fill the hole, as the [stage
   Jacobians](@ref projection-stage-jacobians) predict for a point with ``\partial_x W = 0``.
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
6. At ``x_h = \pi`` the non-separable map fills the hole, with the empty disk ``1.50`` and the
   fill ``0.94``, where the kicks and drifts give ``1.81`` and ``0.47``. Control C does not fill
   it (``2.00`` and ``0.55``), so the filling comes from the dependence of ``\phi_s`` on ``w``, as
   the proposition on [non-separable stages](@ref projection-nonseparable) states. At
   ``x_h = \pi/2`` the non-separable map fills the hole as well (``1.54`` and ``0.94``).
7. The non-separable map has more parameters than the kicks and drifts (3264 against 2304) and a
   coarser basis in ``\xi`` (16 against 32 splines). Its fixed-point iteration needs 50 iterations
   at the fitted ``\theta`` in both runs, the limit that the fit allows. In the run of the script,
   one application of the map takes 0.16 to 0.19 s, against 0.004 to 0.007 s for the kicks and
   drifts; the fit takes 17 and 24 s. The linearised step equations agree with forward-mode
   derivatives along random directions at the fitted ``\theta`` to a relative
   ``2.8 \cdot 10^{-13}`` ([`nonseparable_sensitivity_error`](@ref)); without the terms of ``c``
   and ``b`` they miss by at least ``2.4 \cdot 10^{-2}``.
8. The target weights ``c f`` describe the density without the hole, and the markers before the
   map carry its charge deficit. A map that fills the hole with these weights adds charge in the
   columns of the hole: up to ``9 \%`` of ``\max_i \rho_i`` for the two maps that fill it. The
   charge-conserving reweighting removes it and changes the weights by up to ``16 \%``; the
   velocities ``p / (m w')`` change by the inverse factor, and the fill drops from ``0.94`` to
   ``0.70``. For these two maps the empty disk stays within ``0.03`` of its value without the
   correction; for control A at ``x_h = \pi/2`` it grows from ``1.65`` to ``1.90``. At
   ``\theta = 0`` the correction changes the weights by at most ``1.4 \%``.

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

Hole filling at ``x_h = \pi/2``: markers near the hole (outline) in ``(x, v)``, with
``v = p / (m w)`` for the weights of each panel. Panels (g) and (h) are (e) and (f) after the
charge-conserving reweighting.

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

The two stages of the map, each with its parameters scaled from zero to one, then the reweighting
and the charge-conserving reweighting: the non-separable map at ``x_h = \pi`` and the kicks and
drifts at ``x_h = \pi/2`` ([`remap_animation_frames`](@ref)).

```@example experiments
animate_remap("remap-non-separable-pi.mp4", remap_animation_frames(pi1, :non_separable),
    pi1.problem.hole)
animate_remap("remap-kick-drift-pi2.mp4", remap_animation_frames(pi2, :weight_dependent),
    pi2.problem.hole)
nothing # hide
```

```@raw html
<video src="remap-non-separable-pi.mp4" controls loop muted playsinline style="width:100%"></video>
<video src="remap-kick-drift-pi2.mp4" controls loop muted playsinline style="width:100%"></video>
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
([`weight_errors_by_band`](@ref)); no target weight is negative in any run. All numbers of this
section come from `scripts/projection/dlra_weights.jl`; the test
`test/projection/experiments/dlra_weights.jl` checks each of them.

| source | rank | mass error | ``\lvert v\rvert < 1`` | ``[1, 2)`` | ``[2, 3)`` | ``[3, 4)`` |
|:--|--:|--:|--:|--:|--:|--:|
| ``\vartheta = 10^{-3}`` | 4 | ``1.1 \cdot 10^{-3}`` | ``3.0 \cdot 10^{-2}`` | ``7.5 \cdot 10^{-2}`` | ``1.2 \cdot 10^{-1}`` | ``2.0 \cdot 10^{-1}`` |
| ``\vartheta = 10^{-4}`` | 7 | ``8.5 \cdot 10^{-5}`` | ``3.1 \cdot 10^{-3}`` | ``1.3 \cdot 10^{-2}`` | ``4.1 \cdot 10^{-2}`` | ``1.2 \cdot 10^{-1}`` |
| ``\vartheta = 10^{-5}`` | 10 | ``4.2 \cdot 10^{-6}`` | ``2.7 \cdot 10^{-4}`` | ``1.3 \cdot 10^{-3}`` | ``5.5 \cdot 10^{-3}`` | ``2.0 \cdot 10^{-2}`` |
| ``\vartheta = 10^{-6}`` | 13 | ``1.9 \cdot 10^{-7}`` | ``2.3 \cdot 10^{-5}`` | ``8.7 \cdot 10^{-5}`` | ``5.8 \cdot 10^{-4}`` | ``2.9 \cdot 10^{-3}`` |
| rank one | 1 | ``1.1 \cdot 10^{-13}`` | ``1.5 \cdot 10^{-1}`` | ``1.8 \cdot 10^{-1}`` | ``1.4 \cdot 10^{-1}`` | ``1.8 \cdot 10^{-1}`` |
| projection | | | ``4.1 \cdot 10^{-2}`` | ``1.1 \cdot 10^{-1}`` | ``8.4 \cdot 10^{-2}`` | ``2.0 \cdot 10^{-1}`` |

The second table gives, after the map and the reweighting with ``w^*``, the largest empty disk
over the whole domain ([`largest_empty_disk`](@ref)), the fill of the image of the hole, and the
root mean square deviation of the marker density from its mean ([`density_discrepancy`](@ref));
then the same with the charge-conserving weights, the grid-charge residual of ``w^*``, and
``\max_a |\ln (w'_a / w^*_a)|``. The pushed markers before the map have the empty disk ``1.38``,
the fill ``0`` and the density deviation ``0.099``.

|  | with ``w^*``: empty disk | fill | density | charge-conserving: empty disk | fill | charge | ``\max\lvert\ln\rvert`` |
|:--|--:|--:|--:|--:|--:|--:|--:|
| ``f(T)`` | 1.05 | 0.82 | 0.113 | 1.06 | 0.74 | ``2.4 \cdot 10^{-2}`` | 0.036 |
| ``\vartheta = 10^{-3}`` | 1.30 | 1.33 | 0.364 | 1.28 | 1.25 | ``1.8 \cdot 10^{-2}`` | 0.021 |
| ``\vartheta = 10^{-4}`` | 1.46 | 1.05 | 0.369 | 1.42 | 1.09 | ``1.3 \cdot 10^{-2}`` | 0.021 |
| ``\vartheta = 10^{-5}`` | 1.05 | 0.82 | 0.114 | 1.06 | 0.74 | ``2.4 \cdot 10^{-2}`` | 0.037 |
| ``\vartheta = 10^{-6}`` | 1.05 | 0.82 | 0.113 | 1.06 | 0.74 | ``2.4 \cdot 10^{-2}`` | 0.036 |
| rank one | 1.10 | 0.97 | 0.330 | 3.33 | 0.86 | ``1.7 \cdot 10^{-1}`` | 0.271 |
| projection | 2.67 | 0.78 | 0.510 | 2.74 | 0.74 | ``2.1 \cdot 10^{-2}`` | 0.043 |

In every run ``|\Delta P| / \sum_a |p_a| \leq 3.5 \cdot 10^{-17}``, and the grid charge after the
correction is exact to ``2.9 \cdot 10^{-14}``. The results are:

1. The weight error of ``f_r`` falls with ``\vartheta`` in every band, at ranks between 4 and 13:
   by about one order of magnitude per order of ``\vartheta`` for ``|v| < 1``, by less in the
   tails between ``\vartheta = 10^{-3}`` and ``10^{-4}``. At every tolerance it grows with
   ``|v|``, by a factor of about 100 from ``|v| < 1`` to ``|v| \in [3, 4)`` at
   ``\vartheta = 10^{-6}``: the truncation bounds an absolute error, and ``f`` is small in the
   tails.
2. At ``\vartheta \leq 10^{-5}`` the remap agrees with the remap with the weights of ``f(T)`` to
   the digits of the table. At ``\vartheta = 10^{-4}`` and ``10^{-3}`` the velocities inherit the
   weight errors, and the density deviation triples, although the error in the core is
   ``3 \cdot 10^{-3}`` at ``\vartheta = 10^{-4}``; the error in the tails dominates.
3. The lumped projection gives weight errors of the size of those at ``\vartheta = 10^{-3}`` and
   leaves the largest empty disk at ``2.67``.
4. The charge-conserving reweighting changes the weights by at most ``4 \%`` for every source
   except rank one. As in the hole-filling experiment, it lowers the fill of the image of the
   hole, here from ``0.82`` to ``0.74`` with the weights of ``f(T)``. For rank one it changes the
   weights by up to ``31 \%``, and the largest empty disk grows to ``3.33``.
5. In the run of the script, one BUG solve to ``T = 3`` takes about ``0.1`` s at every tolerance
   once compiled; the fit of the map takes 16 to 36 s.

The field is weak and fixed, and the dimension is one, so the ranks are small. A self-consistent
field and ``d = 3`` are not tested.

The figure below runs the experiment with the fit of the map for the weights of ``f_r`` at
``\vartheta = 10^{-6}`` only (`fit = ("ϑ = 1e-6",)`), a reduced version of the script run; the
weight errors are those of the full run.

```@example experiments
Random.seed!(20261119)
push = dlra_weights_experiment(; fit = ("ϑ = 1e-6",))
plot_dlra_weights(push)
```

Target weights after a push. (a) Pushed markers and the image of the hole. (b) The non-separable
map with the weights of ``f_r`` at ``\vartheta = 10^{-6}``, after the charge-conserving
reweighting. (c) ``f_r - f(T)`` at ``\vartheta = 10^{-4}``. (d) Numerical rank of ``f_r`` over
time. (e) Relative error of the target weights against ``|v^*|``.

```@example experiments
animate_dlra("dlra.mp4", dlra_push_frames(push), push.grid)
nothing # hide
```

```@raw html
<video src="dlra.mp4" controls loop muted playsinline style="width:100%"></video>
```

The low-rank density at ``\vartheta = 10^{-4}`` and the markers along the push
([`dlra_push_frames`](@ref)).

The script `scripts/projection/make_figures.jl` writes all figures as PDF files and the four
animations as MP4 files to `results/projection/`.
