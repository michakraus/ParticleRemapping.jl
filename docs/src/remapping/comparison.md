# [Combining the two constructions](@id remap-comparison)

The gauge rezoning and the canonical transport differ in what they move and in what they keep
exactly; the [properties](@ref remapping-properties) of the chapter list both side by side, and
the [comparison of all constructions](../comparison.md) sets them beside the other chapters.

The two constructions combine. The current-quantile matching gives velocities that are consistent
with an estimate of ``f``; used as [targets](@ref remap-gauge-algorithm) for the gauge rezoning,
they give a remap whose grid charge and kinetic energy are exact. Its grid current is exact only if
``\psi`` preserves the grid current, the discrete form of the compatibility condition.

`remap(::GaugeRezoning, markers; targets)` takes such targets as log-stretches
``\mu^*_a = \ln(V_a / v_a)``. Here the matched velocities of ``\psi = \mathrm{id}`` serve as
targets for the charge- and energy-conserving rezoning; the rezoning keeps the momenta, so the
matched velocities with the sign of ``v`` are admissible targets:

```@example comparison
using ParticleRemapping
grid = PeriodicGrid(2π, 8)
target = SeparableMaxwellianTarget(grid, 0.4)
mk = quantile_matching_setup(; N = 8000, Ng = 8).markers
matched = remap(CurrentQuantileMatching(grid, (x, p) -> (x, p), target), mk).markers
μstar = log.(matched.v ./ mk.v)
r = remap(GaugeRezoning(grid; K = (0, 2), dτ0 = 0.25), mk; targets = μstar)
(τ = r.info.τ, report = conservation_report(mk, r.markers; grid))
```
