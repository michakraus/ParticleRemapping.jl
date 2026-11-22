# [Hole filling](@id rd-holes)

```@meta
CurrentModule = ParticleRemapping.Redistribution
```

This page runs [`SplineRedistribution`](@ref) on a marker set with a hole, outside the
[verification set-up](@ref rd-setup). Here the symplectic step misses its targets
([known issues](@ref redistribution-known-issues) K3). The script
`scripts/redistribution/failure_cases.jl` prints every number of this page.

## A hole: the step misses its targets (K3)

The 244 markers of [`hole_filling_problem`](@ref ParticleRemapping.Projection.hole_filling_problem)
with `nL = 16`, `hole_a = 0.9` and `hole_b = 0.8` (a ``16 \times 16`` lattice on
``[0, 2\pi) \times [-4, 4]`` with an elliptic hole), redistributed by
`SplineRedistribution(; Vmax = 4.0, relax_step = 0.006 * 16 / 256)`, the relaxation step scaled
down with the number of pairs. The remap takes 44 s and returns without an error:

| quantity | value |
|:--|:--|
| target error (scaled) | ``1.98 \cdot 10^{8}`` |
| contraction number | ``4.85 \cdot 10^{8}`` |
| Newton residual at the end point | ``3.96 \cdot 10^{8}`` |
| change of the total weight | ``1.8 \cdot 10^{-15}`` |
| change of the total momentum | ``0.165`` |
| change of the kinetic energy | ``3394`` (of ``2.933``) |

The contraction number ``4.85 \cdot 10^{8}`` puts ``\theta^*`` far outside the regime of the
[parameter solver](@ref rd-linear-parameter-solver), and the forward Newton solve does not
converge: the end point is not a step of the scheme, so neither the targets
nor the total momentum are reached, and only the KKT weights keep the total weight. `remap` does
not report the failure; only `info.target_error` and `info.residual` show it (K3).
