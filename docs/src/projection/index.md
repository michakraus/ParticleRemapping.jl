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

The figures follow the steps on 100 markers of the ``10 \times 10`` lattice whose cells a mild
smooth map moved ([`lattice_example`](@ref ParticleRemapping.lattice_example)`(; nL = 10)`),
remapped by the kick–drift map with two stages on ``8 \times 9`` hat functions
([`algorithm_svg`](@ref ParticleRemapping.algorithm_svg)). The markers are coloured and sized by
weight; open circles are the target points.

1. **Project.** Project the represented density and the marker density onto a B-spline basis,
   ``f_h`` and ``g_h`` (panel 1, the dashed lines are the hat grid).
2. **Targets.** Integrate the density-interpolation flow from ``g_h`` to ``\bar{g}``. It gives
   the target points ``y^*_a``, the target weights ``w^*_a`` and the target momenta
   ``p^*_a = m w^*_a v^*_a`` (panel 2, the arrows point from each marker to its target). The
   targets move the markers toward an even spacing.
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

The markers start at the cell centres ``(\xi, \eta)`` of a ``24 \times 24`` lattice in
``[0, 2\pi) \times [-4, 4]`` with the density
``f = (1 + 0.3 \cos x)(2\pi)^{-1/2} e^{-v^2/2}``, moved by the mild smooth map
``x = (\xi + 0.2 \sin(\xi - 0.5)) \bmod 2\pi``,
``v = \eta + 0.3 \sin(\pi (\eta + 4) / 8)`` with the weights ``w = c f J``, where ``J`` is the
Jacobian of the map ([`lattice_example`](@ref ParticleRemapping.lattice_example)). The kick–drift
map with two stages remaps them
([`KickDriftRemap`](@ref ParticleRemapping.Projection.KickDriftRemap)): the targets come from the
density-interpolation flow on ``12 \times 13`` hat functions, the target weights are ``c f`` from
the exact density at the targets, so that the example shows the map and not the error of a
projection, and the fit ends after its budget of ``20`` iterations, with the return code
`MAXITERS`. The result reports the number of markers, the objective after the fit relative to
``\theta = 0``, the return code of the fit, the number of iterations, and the relative changes of
the total momentum and the total weight. The same example gives the figures `projection-remap.pdf`,
`projection-stages.pdf` and the animation `projection-stages.mp4` of
`scripts/projection/make_figures.jl`.

```@example projection-index
using ParticleRemapping, CairoMakie, Random
Random.seed!(20261119)
ex = lattice_example()
method = KickDriftRemap(; L = ex.L, V = ex.V, projection = HatProjection(12, 13),
    nξ = 16, nv = 20, nw = 8, maxit = 20, density = (xs, vs) -> ex.density.(xs, vs))
result = remap(method, ex.markers)
(; N = ex.N, objective_ratio = result.info.objective_ratio, code = result.info.status.code,
    iterations = result.info.iterations, momentum_change = result.info.momentum_change,
    weight_change = result.info.weight_change)
```

The markers before and after the remap, coloured and sized by weight, with the target points as
open circles:

```@example projection-index
targets = (result.info.targets.x, result.info.targets.v)
plot_remap(ex.markers, result.markers; titles = ("(a) before", "(b) after the remap"),
    targets, period = (0, 2π))
```

The states of the remap ([`remap_states`](@ref ParticleRemapping.remap_states)): the markers
after each kick and each drift at the old weights, with the velocities ``p / (m w)``, and after
the gauge step to the target weights, with the targets in the last panel. A kick moves the markers
only in ``v``, a drift only in ``x``; the gauge step keeps ``(x, p)`` and changes each velocity by
the factor ``w / w^*``. The view is ``|v| \leq 6``.

```@example projection-index
states = remap_states(method, ex.markers, result)
plot_stages(states; columns = 3, period = (0, 2π), window = (0, 2π, -6, 6),
    targets, target_panels = [length(states)])
```

The animation follows the natural path of the map
([`remap_frames`](@ref ParticleRemapping.remap_frames)): the parameters of each kick and each
drift scaled from zero to their fitted values, one after the other, then the gauge step at fixed
``(x, p)``.

```@example projection-index
animate_stages("projection-stages.mp4", remap_frames(method, ex.markers, result);
    period = (0, 2π))
nothing # hide
```

```@raw html
<video src="./projection-stages.mp4" controls loop muted playsinline style="width:100%"></video>
```

The changes of the total weight, the total momentum and the kinetic energy, and the largest
relative changes of the grid moments on 12 linear shape functions in ``x``
([`conservation_report`](@ref ParticleRemapping.conservation_report)):

```@example projection-index
conservation_report(ex.markers, result.markers; grid = PeriodicGrid(2π, 12))
```

The total weight and the total momentum are kept to round-off. The fit ends with `MAXITERS` after
``20`` iterations at ``F / F_0 =`` ``0.220``. The kinetic energy, the grid charge and the grid
current change: the kinetic energy changes by ``0.431`` from ``3.138``, because the map moves the
markers and the target weights differ from the old weights. After the remap, ``7`` of the 576
markers lie slightly outside the velocity box ``|v| \leq 4``, with largest ``|v|`` equal to
``4.5`` (`scripts/projection/example.jl`). The
[charge-conserving reweighting](@ref projection-charge) keeps the grid charge on its own spline
basis.

## [A second example: clustered markers](@id projection-clustered-example)

The second example takes the clustered markers of the spline chapter
([`redistribution_example`](@ref ParticleRemapping.Redistribution.redistribution_example)) on a
``24 \times 24`` lattice (576 markers, mass ``m = 1.5``) in ``[0, 2\pi) \times [-2.4, 2.4]``,
evolved by 60 leapfrog steps of length ``0.05`` in the static force
``1.4 \sin(x - 0.5) + 0.8 \sin 2x``, embedded in the velocity box ``[-3.5, 3.5]`` with the exact
density ``f_0 \circ \Phi^{-1}`` given by the inverse of the leapfrog steps
([`clustered_example`](@ref ParticleRemapping.clustered_example)). Remapping them with the same
kick–drift method on their own box and density sends the velocity tails far outside the box, so
that the axes of the figure and the animation span the escaped markers and the bulk appears as a
narrow band near ``v = 0``:

```@example projection-index
cl = clustered_example()
method2 = KickDriftRemap(; L = cl.L, V = cl.V, projection = HatProjection(12, 13),
    nξ = 16, nv = 20, nw = 8, maxit = 20, density = (xs, vs) -> cl.density.(xs, vs))
r2 = remap(method2, cl.markers)
targets2 = (r2.info.targets.x, r2.info.targets.v)
plot_remap(cl.markers, r2.markers;
    titles = ("(a) before: clustered markers", "(b) after the remap"),
    targets = targets2, period = (0, cl.L))
```

```@example projection-index
animate_stages("projection-clustered.mp4", remap_frames(method2, cl.markers, r2);
    period = (0, cl.L))
nothing # hide
```

```@raw html
<video src="./projection-clustered.mp4" controls loop muted playsinline style="width:100%"></video>
```

Over the run, two runs with the positions moved by ``10^{-13}`` relative, and the CI runners,
``210``–``224`` of the 576 markers leave ``|v| \leq 3.5`` after the remap, reaching ``|v|`` up to
``130``–``180`` and carrying ``0.094``–``0.11`` of the total weight
(`scripts/projection/example.jl`).

The cause is the gauge step at small target weights. The density-interpolation flow spreads the
clustered markers over the whole box ``[0, 2\pi) \times [-3.5, 3.5]``, where the targets in the
velocity tails receive target weights down to ``2.9 \cdot 10^{-6}`` with a median of
``1.2 \cdot 10^{-4}``–``1.4 \cdot 10^{-4}`` among the escaped markers, against a mean weight of
``1.7 \cdot 10^{-3}`` before the remap. The gauge step ``v' = p' / (m w^*)`` divides the momentum
error of the fit by those small weights and amplifies it into large velocities; this behaviour is
recorded as K21 under [known issues](@ref projection-known-issues). Without the exact density, the
lumped projection on ``12 \times 13`` hat functions vanishes at a target of these markers, so
`remap` refuses them with a `DomainError`.

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
- [Structural statements](@ref projection-checks): the checks of every statement.
- [Hole filling](@ref projection-hole-filling): the hole-filling experiment with five
  parametrisations, the charge-conserving reweighting, and the target weights from a low-rank
  density after a push.
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

The two examples of the overview are
[`lattice_example`](@ref ParticleRemapping.lattice_example) and
[`clustered_example`](@ref ParticleRemapping.clustered_example), and
`scripts/projection/example.jl` prints and checks their quoted numbers. The experiments are
functions with the parameters of this chapter as defaults:
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
