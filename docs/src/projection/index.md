# [Hamiltonian particle projection](@id projection-index)

By the [factorisation lemma](@ref factorisation-lemma), every symplectic redistribution is a
canonical map ``\psi`` of ``(x, p)`` followed by a [gauge step](@ref gauge-step). In this chapter,
``\psi`` is the map ``R_\theta`` of a parametrised Hamiltonian that depends on the weights of the
markers: a composition of exactly integrable kicks and drifts, or of non-separable stages that the
implicit midpoint rule integrates. A least-squares fit chooses the parameters ``\theta`` so that
the markers reach the target points of a density-interpolation flow. The gauge step then gives the
markers their new weights at fixed ``(x, p)``: the values of a projected density at the target
points, the minimiser of a relative entropy under the constraint of the grid charge, or the values
of a dynamical low-rank density. For every fixed ``\theta``, the map preserves ``\omega`` exactly,
and it conserves the total momentum ``P`` and the total weight exactly. The charge-conserving
variant also keeps the grid charge. Consistency (R3) and equidistribution (R4) hold up to the fit
residual and the resolution of the projection ([Properties](@ref projection-properties)).

The chapter [What a symplectic redistribution can do](@ref obstructions) states four facts that
shape the construction:

- **The gauge step.** A change of weight at fixed ``(x, p)`` is the identity on ``(x, p)``
  ([gauge step](@ref gauge-step)). The construction puts all its symplectic content into the
  motion of ``(x, p)``, and it assigns the new weights in a final gauge step, the
  [reweighting](@ref projection-final-reweighting).
- **One map for all markers.** A Hamiltonian that does not depend on the weights moves all
  markers with one volume-preserving map ([one map for all markers](@ref one-map-for-all-markers)),
  so [holes persist](@ref holes-persist). The stage Hamiltonians of this chapter therefore depend
  on the weights, through splines in ``\ln w``.
- **Weight-dependent Hamiltonians.** Such a Hamiltonian compresses the marker density at the rate
  of the [compressibility bracket](@ref weight-dependent-hamiltonians), that is only where the old
  weights vary. The [stage Jacobians](@ref projection-stage-jacobians) refine this in ``(x, v)``:
  where ``\partial_x W = 0``, a kick that compresses the markers at ``+v`` expands them at ``-v``.
  The [non-separable stages](@ref projection-nonseparable) do not have this coupling.
- **The velocity–momentum fold.** The implicit equation ``w' = c f_h(x', p' / (m w'))`` for a new
  weight at fixed ``(x', p')`` has two solutions below a bound ``p_{\max}`` and none above it
  ([the velocity–momentum fold](@ref velocity-momentum-fold)). The construction never solves this
  equation: it takes the new weight from the target point and computes the velocity from it.

The fit makes ``\theta`` a function of the state, so the map is symplectic for each fixed
``\theta``, in the sense of [parameters that depend on the state](@ref state-dependent-parameters).
The [spline-projection redistribution](@ref redistribution-overview) also moves the markers to
target points; there the Hamiltonian couples the markers instead of depending on their weights.

The positions lie on the torus ``\mathbb{T}^d = \mathbb{R}^d / L \mathbb{Z}^d``, and the markers
lie in the box ``\Omega = \mathbb{T}^d \times [-V, V]^d`` of ``(x, v)``. An equidistribution has
the marker density ``\bar{g} = N / |\Omega|`` and the weights ``w_a = c \, f(x_a, v_a)`` with
``c = |\Omega| / N`` ([requirements](@ref requirements)). A point of the ensemble phase space
``\mathcal{Z}`` is ``z = (z_1, \dots, z_N)`` with ``z_a = (x_a, p_a)``, and a point of ``(x, v)``
is ``y = (x, v)``. The chapter marks the steps of the map as follows: ``y^*_a = (x^*_a, v^*_a)``
and ``w^*_a`` are the target points and the target weights, ``(x^S_a, p^S_a)`` is the image of
marker ``a`` under ``R_\theta``, and ``x'_a``, ``p'_a``, ``w'_a`` and ``v'_a = p'_a / (m w'_a)``
describe the markers after the final reweighting. The statements hold in ``d`` dimensions; the
code implements the whole algorithm for ``d = 1``.

## The algorithm in short

1. Project the represented density and the marker density onto a B-spline basis, ``f_h`` and
   ``g_h``.
2. Integrate the density-interpolation flow from ``g_h`` to ``\bar{g}``. It gives the target
   points ``y^*_a``, the target weights ``w^*_a`` and the target momenta
   ``p^*_a = m w^*_a v^*_a``.
3. Fit ``\theta`` so that ``R_\theta`` moves each marker close to its target.
4. Apply ``R_\theta``.
5. Reweight at fixed ``(x, p)``: ``(x'_a, p'_a, w'_a) = (x^S_a \bmod L, p^S_a, w^*_a)``, or the
   charge-conserving weights in place of ``w^*_a``.

The [algorithm](@ref projection-algorithm) of the construction page gives the steps in full.

## Example

The markers sit on a ``24 \times 24`` lattice in ``[0, 2\pi) \times [-4, 4]``, with the weights
``w = c f`` for ``f = (1 + 0.3 \cos x)(2\pi)^{-1/2} e^{-v^2/2}``, and an elliptic hole around
``(\pi/2, 1.2)`` removes some of them
([`hole_filling_problem`](@ref ParticleRemapping.Projection.hole_filling_problem)). The
kick–drift map remaps them, with the target weights of the lumped projection on ``12 \times 13``
hat functions and a fit that stops after 10 iterations. The result reports the objective after the
fit relative to ``\theta = 0`` and the changes of the total momentum and the total weight.

```@example projection-index
using ParticleRemapping, CairoMakie, Random
Random.seed!(20261119)
small = hole_filling_problem(; hole_x = π / 2, nL = 24, hole_a = 0.9, hole_b = 0.8)
result = remap(KickDriftRemap(; projection = HatProjection(12, 13), maxit = 10), small.markers)
(; objective_ratio = result.info.objective_ratio,
    momentum_change = result.info.momentum_change, weight_change = result.info.weight_change)
```

```@example projection-index
plot_remap(small.markers, result.markers; titles = ("before", "after the remap"))
```

## The pages of this chapter

- [Construction](@ref projection-construction): the projection, the target points, the
  kick–drift Hamiltonian with the centroid shift, its parametrisation, the stage Jacobians, the
  fit, the reweighting and the algorithm.
- [Non-separable stages](@ref projection-nonseparable): stages generated by one Hamiltonian in
  ``x`` and ``p`` together, integrated with the implicit midpoint rule.
- [Charge-conserving reweighting](@ref projection-charge): the final weights that keep the grid
  charge exactly.
- [Target weights from a dynamical low-rank density](@ref projection-lowrank): target weights from
  a low-rank solution of the Vlasov equation.
- [Properties](@ref projection-properties): what the map keeps exactly and what approximately.
- [Numerical experiments](@ref projection-experiments): the checks of every statement, the
  hole-filling experiment and the target weights after a push.
- [Library: Projection](@ref projection-library): the functions of the module.

## The code

The module `ParticleRemapping.Projection` implements every construction of this chapter for
``d = 1`` (and the stages in any dimension). Three method types implement
[`remap`](@ref ParticleRemapping.remap):

- [`KickDriftRemap`](@ref ParticleRemapping.Projection.KickDriftRemap): the whole algorithm with
  the kick–drift map;
- [`NonSeparableRemap`](@ref ParticleRemapping.Projection.NonSeparableRemap): the whole algorithm
  with the non-separable implicit-midpoint stages;
- [`ChargeConservingReweighting`](@ref ParticleRemapping.Projection.ChargeConservingReweighting):
  the charge-conserving choice of the final weights.

The experiments are functions with the parameters of this chapter as defaults:
[`hole_filling_experiment`](@ref ParticleRemapping.Projection.hole_filling_experiment) and
[`dlra_weights_experiment`](@ref ParticleRemapping.Projection.dlra_weights_experiment). The scripts
`scripts/projection/hole_filling.jl`, `scripts/projection/dlra_weights.jl` and
`scripts/projection/checks.jl` print every number of the
[numerical experiments](@ref projection-experiments), and `scripts/projection/make_figures.jl`
writes the figures and animations to `results/projection/`. The
[verification](@ref projection-verification) lists every test file and script of the chapter.
