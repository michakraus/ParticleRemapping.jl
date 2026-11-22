# [Spline-projection redistribution](@id redistribution-overview)

```@meta
CurrentModule = ParticleRemapping.Redistribution
```

In this chapter the canonical map ``\psi`` is one symplectic step of an ensemble Hamiltonian
``H_\theta(x, p)`` that couples the markers in pairs. The Hamiltonian is a linear combination of
spline functions of the relative positions ``x_a - x_b`` and of the momenta. It is invariant under
translations and permutations of the markers, so the step keeps the total momentum ``P``. Its
parameters ``\theta`` solve a linear system, so that the step moves the markers onto relaxed,
approximately equidistributed targets ``(\bar{x}_a, \bar{v}_a)``. The new weights come from a
positive spline projection of the distribution, corrected by Karush–Kuhn–Tucker (KKT) conditions
so that the total weight and the total momentum are exact. By the
[factorisation lemma](@ref factorisation-lemma), the step followed by this
[gauge step](@ref gauge-step) keeps ``\omega`` exactly for fixed parameters, and it keeps ``P``
and the total weight exactly.

A canonical map applied to each marker on its own moves a hole but does not shrink it
([uncoupled maps](@ref uncoupled-maps)). The pairwise Hamiltonian is free of this obstruction
because it couples the markers. The new weights depend on the projected density at the targets, so
the [velocity–momentum fold](@ref velocity-momentum-fold) rules out a choice of ``\bar{p}_a``
followed by a solve for ``\bar{v}_a``: the chapter specifies the targets in ``(\bar{x}, \bar{v})``
and evaluates ``\bar{w}_a`` and ``\bar{p}_a = m \bar{w}_a \bar{v}_a`` forward. The parameters are
computed from the markers, so the map is symplectic for their frozen values, in the open-loop
reading of [two readings of symplecticity](@ref two-readings); the
[closed-loop identity](@ref closed-loop-identity) states what a closed loop needs. The
[Hamiltonian particle projection](@ref projection-index) also fits the parameters of a Hamiltonian
to target points; here the Hamiltonian does not depend on the weights, and the parameters solve
the linear system of one step. The [sliced redistribution](@ref sliced-overview) also couples the
markers, by an explicit map without parameters.

The ``N`` markers have the positions ``x_a`` on the torus ``\mathbb{T}^d``, the velocities
``v_a`` in the box ``\Omega_v = [-V_{\max}, V_{\max}]^d``, the weights ``w_a`` and the momenta
``p_a = m w_a v_a``; ``z = (x, p) \in \mathbb{R}^{2 N d}`` is the packed state. Bars denote the
targets and the redistributed markers, ``(\bar{x}_a, \bar{v}_a, \bar{w}_a)`` and ``\bar{p}_a``.
Hats denote the end point of a step that does not reach the targets exactly. ``\theta`` are the
``N_\theta`` parameters of the Hamiltonian, and ``p_{\mathrm{scale}}`` its momentum scale. The
statements hold in ``d`` dimensions; the code implements ``d = 1``.

## The algorithm step by step

The figures follow the steps on the 16 markers of the standard example below
([`algorithm_svg`](@ref ParticleRemapping.algorithm_svg)). The markers are coloured and sized by
weight; open circles are the relaxed targets.

1. **Project.** Project the markers onto a tensor-product B-spline basis, periodic in ``x`` and
   clamped in ``v``, by the positive lumped projection (panel 1, the dashed lines are the knots).
2. **Relax.** A few damped steps of repulsive gradient descent move the points ``(x_a, v_a)`` to
   targets ``(\bar{x}_a, \bar{v}_a)`` (panel 2, arrows from the markers to their targets).
3. **Weigh.** The weighted KKT projection gives the target weights
   ``\bar{w}_a = \Delta\bar{\Omega}_a f_h^{\mathrm{KKT}}(\bar{x}_a, \bar{v}_a)``, which keep the
   discrete mass and momentum (panel 3 draws the targets with these weights).
4. **Step.** One linear solve gives the parameters of a pairwise spline Hamiltonian whose
   symplectic step moves ``(x, p)`` to ``(\bar{x}, \bar{p})``; at the old weights the end point
   lies on the targets in ``(x, p)`` but not in ``(x, v)`` (panel 4).
5. **Gauge step.** The end point gets the weights ``\bar{w}_a`` at fixed ``(x, p)`` (panel 5).

```@example rd-steps
using ParticleRemapping
steps = algorithm_svg(:spline)
steps.panels
```

The animation runs through steps 1, 2, 4 and 5 (step 3 moves no marker):

```@example rd-steps
steps.animation
```

The [algorithm](@ref rd-algorithm) page gives the steps with the functions of the package.

## Example

The standard example is the verification set-up [`redistribution_example`](@ref): ``4 \times 4``
markers at the cell centres of a lattice on ``[0, 2\pi) \times [-2.4, 2.4]``, weighted by
``w = \Delta x \Delta v f_0`` with ``f_0 = (1 + 0.45 \cos x) e^{-v^2/2} / (2\pi \sqrt{2\pi})``,
mass ``m = 1.5``, and clustered by 60 leapfrog steps of length ``0.05`` in the static force
``1.4 \sin(x - 0.5) + 0.8 \sin 2x``. [`SplineRedistribution`](@ref) with its defaults
redistributes them by one implicit midpoint step. The figure shows the markers before and after
the remap, coloured and sized by weight, with the relaxed targets as open circles:

```@example rd-index
using ParticleRemapping, CairoMakie
set_theme!(remap_theme())
ex = redistribution_example()
method = SplineRedistribution()
res = remap(method, ex.markers)
targets = (res.info.targets_x, res.info.targets_v)
plot_remap(ex.markers, res.markers; targets = targets, period = (0, ex.Lx),
    titles = ("(a) before: clustered markers", "(b) after: redistributed markers"))
```

The stages of the remap: the markers before; the same markers with the relaxed targets, at which
the KKT weights are computed; the end point of the symplectic step at the old weights,
``v = p / (m w)``, which lies on the targets in ``(x, p)`` but not in ``(x, v)``; and the markers
after the gauge step, which assigns the KKT weights at fixed ``(x, p)``
([`remap_states`](@ref ParticleRemapping.remap_states)):

```@example rd-index
plot_stages(remap_states(method, ex.markers, res); targets = targets, period = (0, ex.Lx))
```

The animation follows the symplectic step in pseudo-time, ``\Psi^{\theta}_{s\tau}(z)`` for
``s \in [0, 1]`` at the old weights ([`pseudotime_path`](@ref)), and then the gauge step
([`remap_frames`](@ref ParticleRemapping.remap_frames)):

```@example rd-index
animate_stages("redistribution-stages.mp4", remap_frames(method, ex.markers, res);
    targets = targets, period = (0, ex.Lx))
nothing # hide
```

```@raw html
<video src="redistribution-stages.mp4" controls loop muted playsinline style="width:100%"></video>
```

The changes of the conserved quantities
([`conservation_report`](@ref ParticleRemapping.conservation_report)), and the scaled distance of
the end point of the step from the targets:

```@example rd-index
(conservation_report(ex.markers, res.markers)..., target_error = res.info.target_error)
```

The total weight and the total momentum are kept to round-off. The kinetic energy is not a
constraint of the construction and changes by ``0.1323`` of ``1.0812``
(`scripts/redistribution/failure_cases.jl`). The [verification set-up](@ref rd-setup) gives the
example in more detail and shows where the relaxation fails; [hole filling](@ref rd-holes) shows
where the symplectic step fails.

## The pages of this chapter

- [Positive, moment-compatible spline projection](@ref rd-projection): the ``L^2`` projection
  keeps mass and momentum but overshoots to negative values in sparse regions; the de-clustering
  relaxation gives the targets; the weighted KKT projection gives positive weights with the exact
  discrete mass and momentum that the parameter solver needs.
- [Symmetric ensemble Hamiltonians](@ref rd-hamiltonians): uncoupled Hamiltonians cannot exchange
  momentum while they keep the total momentum, and the subtraction of the mean force breaks
  symplecticity; the pairwise, Fourier-moment and centroid-frame Hamiltonians are invariant under
  permutations and translations.
- [Symplectic step and parameter solver](@ref rd-solver): for ``H_\theta = \sum_k \theta_k
  \Psi_k``, one step of the implicit midpoint or the symplectic Euler method gives a linear system
  for ``\theta`` on the ``(2 N d - d)``-dimensional subspace of zero net momentum; the forward
  Newton solve reaches the identity-connected branch when
  ``\tfrac{\tau}{2} \| \mathbb{J} \nabla^2 H_\theta \| < 1``; caustics, substeps and the coarse
  basis.
- [Closed loops](@ref rd-closed-loop): the [closed-loop identity](@ref closed-loop-identity)
  evaluated on the parameter maps of this chapter, the fixed-grid rule, the closed-loop cotangent
  lifts and the two closed loops of the verification.
- [Algorithm](@ref rd-algorithm): the five steps of [`SplineRedistribution`](@ref), from the
  positive projection to the gauge step.
- [Numerical experiments](@ref rd-experiments): the standard experiments, which every
  construction runs.
- [The verification set-up](@ref rd-setup): the standard example in detail, the schemes and
  substeps, and the relaxation outside the verification set-up (K2).
- [Hole filling](@ref rd-holes): on a marker set with a hole, the symplectic step misses its
  targets and `remap` does not report it (K3).
- [Properties, limitations and higher dimensions](@ref rd-limitations-page): what the map keeps,
  what it cannot do, the known issues, the open questions, and the cost in 3d3v.
- [Library: Redistribution](@ref rd-library): the functions of the module.

## The code

The module `ParticleRemapping.Redistribution` implements the construction for ``d = 1``, with
positions on ``[0, 2\pi)`` and velocities in ``[-3.5, 3.5]`` by default. The method type
[`SplineRedistribution`](@ref) implements [`remap`](@ref ParticleRemapping.remap),
[`remap_states`](@ref ParticleRemapping.remap_states) and
[`remap_frames`](@ref ParticleRemapping.remap_frames); its keywords set the spline basis, the
relaxation, the Hamiltonian basis, the scheme and the substeps. The building blocks are public:
the projections ([`l2_spline_projection`](@ref), [`kkt_spline_projection`](@ref)), the
relaxation ([`declustering_relaxation`](@ref)), the Hamiltonians
([`PairwiseSplineHamiltonian`](@ref), [`FourierMomentHamiltonian`](@ref),
[`CentroidFrameHamiltonian`](@ref)), the parameter solve and the step
([`solve_redistribution_parameters`](@ref), [`forward_redistribution_step`](@ref)) and the
closed loops ([`closed_loop_defect`](@ref)).

The script `scripts/redistribution/verify_symplectic_redistribution.jl` checks every statement of
the note and prints its numbers; `scripts/redistribution/failure_cases.jl` prints the numbers of
the [verification set-up](@ref rd-setup), of [hole filling](@ref rd-holes) and of the
[limitations](@ref rd-limitations-page); `scripts/redistribution/make_figures.jl` writes the
figures and the animations to `results/redistribution/`. The
[verification](@ref rd-verification) lists every test file and script of the chapter.
