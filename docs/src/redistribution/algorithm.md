# [Algorithm](@id rd-algorithm)

```@meta
CurrentModule = ParticleRemapping.Redistribution
```

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
