# [Structural statements and target weights after a push](@id projection-checks)

```@meta
CurrentModule = ParticleRemapping.Projection
```

This page gives the checks of the structural statements and the target weights after a push.
Every statement of this chapter has a test under `test/projection/`, and each test runs the case
that must fail. Derivatives are computed by forward-mode automatic differentiation, so the
tolerances are set by rounding error. The script `scripts/projection/checks.jl` runs the tests of
the structural statements and prints the measured numbers of the first table; the script
`scripts/projection/dlra_weights.jl` (about 3 minutes) prints the tables of the push experiment and
checks every number against the value quoted here.

## The structural statements

The first three rows check statements of
[Symplectic redistribution](@ref obstructions) with the functions of this package;
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
<video src="dlra.mp4" controls loop muted playsinline style="width:100%"></video>
```

The low-rank density at ``\vartheta = 10^{-4}`` and the markers along the push
([`dlra_push_frames`](@ref)).
