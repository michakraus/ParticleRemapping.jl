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

Bars denote the targets and the redistributed markers, ``(\bar{x}_a, \bar{v}_a, \bar{w}_a)`` and
``\bar{p}_a``. Hats denote the end point of a step that does not reach the targets exactly.

## The algorithm

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
   ([`forward_redistribution_step`](@ref), [`redistribution_substeps`](@ref)).
5. Assign the weights ``\bar{w}_a`` to the end point of the step, and the velocities
   ``p_a / (m \bar{w}_a)`` with the momenta ``p_a`` of the end point. By the
   [factorisation lemma](@ref factorisation-lemma), every positive split of the momenta keeps
   ``\omega``.

## The verification set-up

The verification set-up [`redistribution_example`](@ref) produces a clustered ensemble:
``4 \times 4`` markers at the cell centres of a lattice on ``[0, 2\pi) \times [-2.4, 2.4]``,
weighted by ``w = \Delta x \Delta v f_0`` with
``f_0 = (1 + 0.45 \cos x) e^{-v^2/2} / (2\pi \sqrt{2\pi})``, mass ``m = 1.5``, advanced by 60
leapfrog steps of length ``0.05`` in the static force ``1.4 \sin(x - 0.5) + 0.8 \sin 2x``.

```@example rd-index
using ParticleRemapping, CairoMakie
ex = redistribution_example()
plot_remap(ex.initial, ex.markers; titles = ("lattice", "after 60 leapfrog steps"))
```

The example below runs the redistribution on these 16 markers: the markers move onto the relaxed
targets (red circles) and take the KKT weights, while the total weight and momentum are kept.

```@example rd-index
res = remap(SplineRedistribution(), ex.markers)
(target_error = res.info.target_error,
 weight_change = total_weight(res.markers) - total_weight(ex.markers),
 momentum_change = total_momentum(res.markers) - total_momentum(ex.markers))
```

```@example rd-index
plot_redistribution_markers(ex.markers, res.markers;
    targets = (res.info.targets_x, res.info.targets_v), Lx = ex.Lx)
```

## The pages

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

The checks, the controls and the comparison with the numbers of the sources are in
[the verification](@ref rd-verification).

The package implements the construction for ``d = 1``, with positions on ``[0, 2\pi)`` and
velocities in ``[-3.5, 3.5]``; the statements of this chapter hold in ``d`` dimensions.
