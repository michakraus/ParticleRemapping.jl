# [Gauge rezoning and canonical transport](@id remap-overview)

```@meta
CurrentModule = ParticleRemapping.Remapping
```

By the [factorisation lemma](@ref factorisation-lemma), every symplectic redistribution is a
canonical map ``\psi`` of ``(x, p)``, symplectic for ``\omega = \sum_a \mathrm{d} x_a \wedge
\mathrm{d} p_a``, followed by a [gauge step](@ref gauge-step). This chapter works out the two
extreme cases of this factorisation in 1D-1V, on a periodic grid with ``N_g`` nodes and linear
shape functions ``S_i``. Both are exactly symplectic for ``\omega_N``.

- The **gauge rezoning** takes ``\psi = \mathrm{id}``. It moves no marker in ``(x, p)``: it
  stretches each velocity along its ray, ``V_a = \lambda_a v_a``, and divides the weight by the
  same ``\lambda_a > 0``. A constrained optimisation chooses the stretches. The positions, the
  momenta, the total momentum ``P`` and the grid current are exact by construction, and the grid
  charge and the kinetic energy are exact by constraints. The gauge step preserves ``\omega_N``
  also for stretches that depend on the whole ensemble, so the remap itself needs no
  closed-loop argument ([parameters that depend on the state](@ref state-dependent-parameters));
  only the composition with the flow is subject to the [two readings](@ref two-readings).
- The **canonical transport** takes a given canonical map ``\psi``, the same for all markers, and
  chooses the new weights by *current-quantile matching*, so that the markers represent the same
  density ``f``. In the continuum limit these weights are unique, and they exist if and only if
  ``\psi`` preserves the sign-resolved current profile. One map for all markers
  [moves holes but keeps their area](@ref holes-persist); the
  [limit of equidistribution](@ref remap-limit) states this obstruction for consistent weights.

Throughout the chapter, lower-case letters ``x_a, v_a, w_a, p_a`` denote a marker before the
remap and capitals ``X_a, V_a, W_a, P_a`` after it; a star marks a target, ``\mu^*_a`` or
``V^*_a``.

## The algorithms in short

**Gauge rezoning** ([`GaugeRezoning`](@ref)), with log-stretches ``\mu_a = \ln \lambda_a``:

1. In each cell and for each sign of ``v``, set targets ``\mu^*_a`` that spread the speeds evenly.
2. Minimise the weighted distance ``\sum_a w_a (\mu_a - \mu^*_a)^2`` subject to unchanged grid
   charge and grid kinetic energy, by a Newton method on the KKT system, and follow the
   minimiser by a homotopy from ``\mu^* = 0`` to the full targets. The homotopy stops early where
   the solution branch has a fold.
3. Set ``W_a = w_a e^{-\mu_a}`` and ``V_a = p_a / (m W_a)``.

**Canonical transport** ([`CurrentQuantileMatching`](@ref)):

1. Map ``(X_a, P_a) = \psi(x_a, p_a)``.
2. In each cell and for each sign of ``P``, sort the markers by ``|P_a|`` and form the cumulative
   momenta.
3. Give each marker the velocity ``V_a`` at the same quantile of the cumulative current of ``f``,
   and the weight ``W_a = P_a / (m V_a)``.

Steps 2 and 3 of the transport are a gauge step. The velocities they give can serve as targets for
the gauge rezoning, which then makes the grid charge and the kinetic energy exact
([combining the two constructions](@ref remap-comparison)).

```@example overview
using ParticleRemapping, CairoMakie
setup = gauge_rezoning_setup()
result = remap(setup.method, setup.markers)
plot_remap(setup.markers, result.markers; titles = ("clustered markers",
    "after the gauge rezoning, τ = $(round(result.info.τ; digits = 3))"), markersize = 2)
```

The figure shows the [numerical example of the gauge rezoning](@ref remap-gauge-example): 6400
markers whose velocities cluster in each cell, before and after the rezoning, coloured by weight.
The positions and the momenta are the same in both panels.

## Pages of this chapter

- [Velocity rezoning by the gauge step](@ref remap-gauge): the map and its exact properties, the
  constrained optimisation, its optimality conditions, the Newton step and the homotopy, two
  worked examples and a numerical example.
- [Canonical transport with consistent weights](@ref remap-transport): the matching proposition,
  the compatibility condition, the limit of equidistribution, admissible canonical maps, the
  discrete algorithm, two worked examples and a numerical example.
- [Combining the two constructions](@ref remap-comparison): what each construction keeps exactly,
  and the matched velocities as targets of the gauge rezoning.
- [Library: Remapping](library.md): the docstrings of the module [`Remapping`](@ref).

The tests and scripts that check every claim and every number of this chapter are listed in the
[verification](@ref remap-verification).
