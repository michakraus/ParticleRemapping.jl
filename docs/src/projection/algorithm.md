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
the algorithm on the deformed ``24 \times 24`` lattice ([`lattice_example`](@ref
ParticleRemapping.lattice_example)) with the target weights ``c f``. The example below runs it on
the same markers with the target weights from the lumped projection (the default) on
``12 \times 13`` hats, coarser than the markers, and a fit that ends after its budget of 10
iterations, with the return code `MAXITERS`. On hats finer than the marker spacing the lumped
projection vanishes between the markers, and `remap` refuses the zero target weights.

```@example construction
using ParticleRemapping, CairoMakie, Random
Random.seed!(20261119)
ex = lattice_example()
result = remap(KickDriftRemap(; L = ex.L, V = ex.V, projection = HatProjection(12, 13),
    nξ = 16, nv = 20, nw = 8, maxit = 10), ex.markers)
(; objective_ratio = result.info.objective_ratio, iterations = result.info.iterations,
    code = result.info.status.code, momentum_change = result.info.momentum_change,
    weight_change = result.info.weight_change)
```

```@example construction
plot_remap(ex.markers, result.markers;
    titles = ("(a) before", "(b) after the remap, lumped projection"), period = (0, ex.L))
```

The errors of the lumped projection pass to the velocities through ``v' = p' / (m w^*)``
([reweighting](@ref projection-final-reweighting)): the lumped target weights in the tails are
larger than ``c f`` there, so the gauge step pulls the tail markers inward.

## The steps of the map

The kick–drift map with the charge-conserving reweighting on 12 periodic cubic B-splines on the
deformed lattice, here with a fit that ends after its budget of 10 iterations (`MAXITERS`): the
markers before the map, after each kick and each drift (velocities ``p / (m w)`` with the old
weights), after the reweighting with the target weights and after the charge-conserving
reweighting ([`remap_states`](@ref ParticleRemapping.remap_states); view ``|v| \leq 6``).

```@example construction
method = KickDriftRemap(; L = ex.L, V = ex.V, projection = HatProjection(12, 13),
    nξ = 16, nv = 20, nw = 8, maxit = 10, density = (xs, vs) -> ex.density.(xs, vs),
    charge = ChargeConservingReweighting(ex.L, 12))
r_charge = remap(method, ex.markers)
plot_stages(remap_states(method, ex.markers, r_charge); columns = 4, period = (0, ex.L),
    window = (0, ex.L, -6, 6))
```

A kick moves the markers only in ``v``, a drift only in ``x``; the reweighting changes the
velocities by ``w/w^*``.
