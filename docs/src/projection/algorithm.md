# [Algorithm](@id projection-algorithm)

```@meta
CurrentModule = ParticleRemapping.Projection
```

The algorithm runs the steps of the [construction](@ref projection-construction) in this order:

1. Project ``f_h`` and ``g_h``.
2. Solve the Poisson problem and integrate the density-interpolation flow from ``t = 0`` to ``1``
   for every marker to get ``y^*_a``. Compute ``w^*_a`` and ``p^*_a``.
3. Lift the positions to ``\mathbb{R}^d`` and compute ``w_{\mathrm{tot}}`` and ``\alpha_a``.
4. Minimise ``J`` over ``\theta``, starting from ``\theta = 0``.
5. Apply ``R_\theta`` with the optimal ``\theta``.
6. Apply the reweighting ``(x'_a, p'_a, w'_a) = (x^S_a \bmod L, p^S_a, w^*_a)``.

[`remap`](@ref) with a [`KickDriftRemap`](@ref) runs these steps for 1D-1V
[`Markers`](@ref ParticleRemapping.Markers); [`NonSeparableRemap`](@ref) replaces the kicks and
drifts by non-separable stages, and the option `charge` adds the [charge-conserving
reweighting](@ref projection-charge). The [example](@ref projection-example) of the overview runs
the algorithm on a lattice of ``24 \times 24`` markers with the target weights ``c f``. The example
below runs it on the same markers with the target weights from the lumped projection (the
default) on ``12 \times 13`` hats, coarser than the markers, and a fit stopped after 10
iterations. On hats finer than the marker spacing the lumped projection vanishes between the
markers, and `remap` refuses the zero target weights.

```@example construction
using ParticleRemapping, CairoMakie, Random
Random.seed!(20261119)
small = hole_filling_problem(; hole_x = π / 2, nL = 24, hole_a = 0.9, hole_b = 0.8)
result = remap(KickDriftRemap(; projection = HatProjection(12, 13), maxit = 10), small.markers)
(; objective_ratio = result.info.objective_ratio, iterations = result.info.iterations,
    momentum_change = result.info.momentum_change, weight_change = result.info.weight_change)
```

```@example construction
plot_remap(small.markers, result.markers;
    titles = ("(a) before", "(b) after the remap, lumped projection"), period = (0, 2π))
```

The errors of the lumped projection pass to the velocities through ``v' = p' / (m w^*)``
([reweighting](@ref projection-final-reweighting)), and the markers of the tails scatter.

## The steps of the map

The kick–drift map of the [hole-filling experiment](@ref projection-hole-filling) at
``x_h = \pi/2``, here with the fit stopped after 10 iterations: the markers near the hole after
each kick and each drift (velocities ``p/(m w)`` with the old weights), after the reweighting with
the target weights and after the charge-conserving reweighting
([`experiment_remap`](@ref), [`remap_states`](@ref ParticleRemapping.remap_states),
[`plot_hole_stages`](@ref)).

```@example construction
demo = hole_filling_experiment(; hole_x = π / 2, variants = (:weight_dependent,), maxit = 10)
plot_hole_stages(remap_states(experiment_remap(demo, :weight_dependent)...), demo.problem.hole)
```

A kick moves the markers only in ``v``, a drift only in ``x``; the reweighting changes the
velocities by ``w/w^*``.
