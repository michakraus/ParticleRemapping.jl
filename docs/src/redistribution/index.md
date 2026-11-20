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

## The algorithm in short

The method [`SplineRedistribution`](@ref) runs the whole redistribution of 1D-1V markers,
through [`remap`](@ref ParticleRemapping.remap):

1. Project the markers onto a tensor-product B-spline basis, periodic in ``x`` and clamped in
   ``v``, by the positive lumped projection ``c^{\mathrm{pos}}`` ([`lumped_spline_projection`](@ref)).
2. Relax the points ``(x_a, v_a)`` by a few damped steps of repulsive gradient descent to targets
   ``(\bar{x}_a, \bar{v}_a)`` ([`declustering_relaxation`](@ref)).
3. Weigh the targets by the weighted KKT projection, ``\bar{w}_a = \Delta\bar{\Omega}_a
   f_h^{\mathrm{KKT}}(\bar{x}_a, \bar{v}_a)``, with local volumes ``\Delta\bar{\Omega}_a`` from the
   number density of the targets, so that ``\bar{p}_a = m \bar{w}_a \bar{v}_a`` conserves the
   discrete mass and momentum ([`kkt_redistribution_weights`](@ref)).
4. Solve the linear system ``G(z_{\mathrm{stage}}) \theta = \bar{z} - z``, with
   ``z = (x, p)``, for the parameters of a pairwise spline Hamiltonian
   ([`solve_redistribution_parameters`](@ref)), and advance ``(x, p)`` by one step of the implicit
   midpoint or the symplectic Euler method, or by a few substeps
   ([`forward_redistribution_step`](@ref), [`redistribution_substeps`](@ref)). This is the
   canonical map ``\psi``.
5. Assign the weights ``\bar{w}_a`` to the end point of the step, and the velocities
   ``p_a / (m \bar{w}_a)`` with the momenta ``p_a`` of the end point: the gauge step. By the
   [factorisation lemma](@ref factorisation-lemma), every positive split of the momenta keeps
   ``\omega``.

The [positive projection](@ref rd-projection), the [Hamiltonians](@ref rd-hamiltonians) and the
[parameter solver](@ref rd-solver) give each step in full.

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
(`scripts/redistribution/failure_cases.jl`). The [numerical experiments](@ref rd-experiments)
give the example in more detail and show where the method fails.

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
- [Numerical experiments](@ref rd-experiments): the standard example in detail, the schemes and
  substeps, and the failures outside the verification set-up.
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
the [numerical experiments](@ref rd-experiments) and of the
[limitations](@ref rd-limitations-page); `scripts/redistribution/make_figures.jl` writes the
figures and the animations to `results/redistribution/`. The
[verification](@ref rd-verification) lists every test file and script of the chapter.
