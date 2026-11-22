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
residual and the resolution of the projection ([properties](@ref projection-properties)).

The chapter [Symplectic redistribution](@ref obstructions) states four facts that
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

## The algorithm step by step

The figures follow the steps on 96 markers of a ``10 \times 10`` lattice with an elliptic hole
(dark outline), remapped by the kick–drift map with two stages on ``8 \times 9`` hat functions
([`algorithm_svg`](@ref ParticleRemapping.algorithm_svg)). The markers are coloured and sized by
weight; open circles are the target points.

1. **Project.** Project the represented density and the marker density onto a B-spline basis,
   ``f_h`` and ``g_h`` (panel 1, the dashed lines are the hat grid).
2. **Targets.** Integrate the density-interpolation flow from ``g_h`` to ``\bar{g}``. It gives
   the target points ``y^*_a``, the target weights ``w^*_a`` and the target momenta
   ``p^*_a = m w^*_a v^*_a`` (panel 2, the arrows point from each marker to its target). The
   targets cover the hole.
3. **Fit and map.** Fit ``\theta`` so that ``R_\theta`` moves each marker close to its target,
   and apply ``R_\theta``, the canonical map ``\psi``, at the old weights: each kick moves the
   markers in ``v`` only, each drift in ``x`` only (panels 3 to 6). The stages depend on the
   weights, so markers of different weight move differently.
4. **Reweight.** Set the weights at fixed ``(x, p)``, the gauge step:
   ``(x'_a, p'_a, w'_a) = (x^S_a \bmod L, p^S_a, w^*_a)``, or the charge-conserving weights in
   place of ``w^*_a`` (panel 7).

```@example projection-steps
using ParticleRemapping
steps = algorithm_svg(:projection)
steps.panels
```

The animation runs through the same steps; between two steps the positions, the velocities and
the weights change linearly:

```@example projection-steps
steps.animation
```

The [algorithm](@ref projection-algorithm) page gives the steps in full.

## [Example](@id projection-example)

The markers sit on a ``24 \times 24`` lattice in ``[0, 2\pi) \times [-4, 4]``, with the weights
``w = c f`` for ``f = (1 + 0.3 \cos x)(2\pi)^{-1/2} e^{-v^2/2}``, and an elliptic hole with the
semi-axes ``0.9`` and ``0.8`` around ``(\pi/2, 1.2)`` removes some of them
([`hole_filling_problem`](@ref ParticleRemapping.Projection.hole_filling_problem)). The
kick–drift map with two stages remaps them
([`KickDriftRemap`](@ref ParticleRemapping.Projection.KickDriftRemap)): the targets come from the
density-interpolation flow on ``12 \times 13`` hat functions, the target weights are ``c f`` at the
targets, so that the example shows the map and not the error of a projection, and the fit stops
after 20 iterations. The result reports the number of markers, the objective after the fit
relative to ``\theta = 0`` and the relative changes of the total momentum and the total weight.
The same example gives the figures `projection-remap.pdf`, `projection-stages.pdf` and the
animation `projection-stages.mp4` of `scripts/projection/make_figures.jl`.

```@example projection-index
using ParticleRemapping, CairoMakie, Random
Random.seed!(20261119)
small = hole_filling_problem(; hole_x = π / 2, nL = 24, hole_a = 0.9, hole_b = 0.8)
method = KickDriftRemap(; projection = HatProjection(12, 13), maxit = 20,
    density = (xs, vs) -> small.density.(xs, vs))
result = remap(method, small.markers)
(; N = small.N, objective_ratio = result.info.objective_ratio,
    momentum_change = result.info.momentum_change, weight_change = result.info.weight_change)
```

The markers before and after the remap, coloured and sized by weight, with the target points as
open circles:

```@example projection-index
targets = (result.info.targets.x, result.info.targets.v)
plot_remap(small.markers, result.markers; titles = ("(a) before", "(b) after the remap"),
    targets, period = (0, 2π))
```

The states of the remap ([`remap_states`](@ref ParticleRemapping.remap_states)): the markers
after each kick and each drift at the old weights, with the velocities ``p / (m w)``, and after
the gauge step to the target weights, with the targets in the last panel. A kick moves the markers
only in ``v``, a drift only in ``x``; the gauge step keeps ``(x, p)`` and changes each velocity by
the factor ``w / w^*``. The view is ``|v| \leq 5``.

```@example projection-index
states = remap_states(method, small.markers, result)
hole! = (ax, _) -> plot_hole_outline!(ax, small.hole)
plot_stages(states; columns = 3, period = (0, 2π), window = (0, 2π, -5, 5), decorate! = hole!,
    targets, target_panels = [length(states)])
```

The animation follows the natural path of the map
([`remap_frames`](@ref ParticleRemapping.remap_frames)): the parameters of each kick and each
drift scaled from zero to their fitted values, one after the other, then the gauge step at fixed
``(x, p)``.

```@example projection-index
animate_stages("projection-stages.mp4", remap_frames(method, small.markers, result);
    period = (0, 2π), decorate! = hole!)
nothing # hide
```

```@raw html
<video src="projection-stages.mp4" controls loop muted playsinline style="width:100%"></video>
```

The changes of the total weight, the total momentum and the kinetic energy, and the largest
relative changes of the grid moments on 12 linear shape functions in ``x``
([`conservation_report`](@ref ParticleRemapping.conservation_report)):

```@example projection-index
conservation_report(small.markers, result.markers; grid = PeriodicGrid(2π, 12))
```

The total weight and the total momentum are kept to round-off. The kinetic energy, the grid charge
and the grid current change: the map moves the markers and the target weights differ from the old
weights. The [charge-conserving reweighting](@ref projection-charge) keeps the grid charge on its
own spline basis.

## The pages of this chapter

- [Construction](@ref projection-construction): the projection, the target points, the
  kick–drift Hamiltonian with the centroid shift, its parametrisation, the stage Jacobians, the
  fit and the reweighting.
- [Non-separable stages](@ref projection-nonseparable): stages generated by one Hamiltonian in
  ``x`` and ``p`` together, integrated with the implicit midpoint rule.
- [Charge-conserving reweighting](@ref projection-charge): the final weights that keep the grid
  charge exactly.
- [Target weights from a dynamical low-rank density](@ref projection-lowrank): target weights from
  a low-rank solution of the Vlasov equation.
- [Algorithm](@ref projection-algorithm): the steps of the remap in order, a run with the target
  weights from the lumped projection, and the states of the kick–drift map.
- [Numerical experiments](@ref projection-experiments): the standard experiments, which every
  construction runs.
- [Structural statements and target weights after a push](@ref projection-checks): the checks of
  every statement and the target weights from a low-rank density after a push.
- [Hole filling](@ref projection-hole-filling): the hole-filling experiment with five
  parametrisations and the charge-conserving reweighting.
- [Properties, limitations and higher dimensions](@ref projection-properties): what the map keeps
  exactly and what approximately, what it cannot do, the known issues, the open questions, and the
  sizes and costs in 3d3v.
- [Library: Projection](@ref projection-library): the functions of the module.

## The code

The module `ParticleRemapping.Projection` implements every construction of this chapter for
``d = 1`` (and the stages in any dimension). Three method types implement
[`remap`](@ref ParticleRemapping.remap), [`remap_states`](@ref ParticleRemapping.remap_states) and,
for the first two, [`remap_frames`](@ref ParticleRemapping.remap_frames):

- [`KickDriftRemap`](@ref ParticleRemapping.Projection.KickDriftRemap): the whole algorithm with
  the kick–drift map;
- [`NonSeparableRemap`](@ref ParticleRemapping.Projection.NonSeparableRemap): the whole algorithm
  with the non-separable implicit-midpoint stages;
- [`ChargeConservingReweighting`](@ref ParticleRemapping.Projection.ChargeConservingReweighting):
  the charge-conserving choice of the final weights.

The experiments are functions with the parameters of this chapter as defaults:
[`hole_filling_experiment`](@ref ParticleRemapping.Projection.hole_filling_experiment) and
[`dlra_weights_experiment`](@ref ParticleRemapping.Projection.dlra_weights_experiment);
[`experiment_remap`](@ref ParticleRemapping.Projection.experiment_remap) gives the states and the
frames of each map of the hole-filling experiment. The scripts
`scripts/projection/hole_filling.jl`, `scripts/projection/dlra_weights.jl` and
`scripts/projection/checks.jl` print every number of the
[structural statements](@ref projection-checks) and of the
[hole filling](@ref projection-hole-filling), `scripts/projection/higher_dimensions.jl`
prints the sizes of the [higher dimensions](@ref projection-higher-dimensions), and
`scripts/projection/make_figures.jl` writes the figures and animations to `results/projection/`.
The [verification](@ref projection-verification) lists every test file and script of the chapter.
