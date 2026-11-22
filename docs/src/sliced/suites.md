# [Structural suites](@id sliced-suites)

```@meta
CurrentModule = ParticleRemapping.SlicedRedistribution
```

Four suites verify the [symplecticity condition](@ref symplecticity-conditions) and the
statements of [the stages and their structure](@ref sliced-structure); suites 5a and 5b, which
test the redistribution on two phase-space distributions, are on the pages
[hole filling](@ref sliced-holes) and [long-time run](@ref sliced-long-time-page). The
symplecticity defect of a map ``F`` of the
packed state ``z = (x, p)`` is the largest entry of
``\mathrm{D}F^\top \mathbb{J} \, \mathrm{D}F - \mathbb{J}`` with the canonical matrix
``\mathbb{J}`` ([`canonical_J`](@ref ParticleRemapping.canonical_J)); the Jacobians come from
ForwardDiff. The thresholds in the tables are those of the source verification script
([verification](@ref sliced-verification)). The numbers come from
`scripts/sliced/verify_sliced_redistribution.jl` (suites 1–5) and
`scripts/sliced/long_time_experiment.jl` (the variants and the runs with remaps); the tests
`test/sliced/experiments.jl` and `test/sliced/experiments_long_time.jl` check them.

## [Suite 1: factorisation on the full state space](@id sliced-factorisation-suite)

[`factorisation_pullback_experiment`](@ref): eight random markers (``m = 1.35``). The pull-back of
``\sum_a \mathrm{d}\bar{x}_a \wedge \mathrm{d}\bar{p}_a`` under each map of ``(x, v, w)`` is
compared with ``\omega_N``.

| map | defect | threshold |
|:--|--:|:--|
| sliced lift of ``(x, p)``, weights ``\bar{w} \propto w (1 + 0.3 \cos\bar{x}\, e^{-v^2/4})``, ``\bar{v} = \bar{p}/(m\bar{w})`` | ``1.1 \times 10^{-15}`` | ``< 10^{-13}`` |
| the same with ``\bar{v} = \bar{p}/(m w)`` (no recoil, control) | ``0.507`` | ``> 10^{-1}`` |
| ``\bar{p} = p - w \sin x`` on the leaf ``\mathrm{d}w = 0`` | ``0`` | ``< 10^{-14}`` |
| ``\bar{p} = p - w \sin x`` across ``\mathrm{d}w`` (control) | ``0.974`` | ``> 10^{-1}`` |

The first row is the full-space form of the [condition](@ref symplecticity-conditions)
with a weight rule that depends on ``v`` and ``w``; the last two rows show the leaf form: a
canonical map that reads ``w`` is symplectic on each leaf and not across leaves.

## [Suite 2: pseudo-time flow](@id sliced-flow-suite)

[`pseudo_time_flow_experiment`](@ref): six markers (``m = 1.1``), one implicit-midpoint step
``\Delta\tau = 0.12`` of the generator
``\mathcal{G}(x, p) = \frac{1}{4} \sum_{a,b} e^{-(x_a - x_b)^2/2} (p_a - p_b)^2`` (25 fixed-point
sweeps), weights ``\bar{w}_a = w_a \exp(0.25 \, \Delta\tau \cos\bar{x}_a)``.

| velocities | defect | threshold |
|:--|--:|:--|
| ``\bar{v} = \bar{p}/(m\bar{w})`` (with recoil) | ``3.2 \times 10^{-15}`` | ``< 10^{-13}`` |
| ``\bar{v} = \bar{p}/(m w)`` (no recoil, control) | ``0.0306`` | ``> 10^{-2}`` |
| with recoil, 3 fixed-point sweeps only | ``0.0113`` | (not converged) |

The implicit midpoint rule is symplectic for any Hamiltonian [HairerLubichWanner:2006](@cite) once
its fixed-point equation is solved; the last row shows that the check sees an unconverged step.

## [Suite 3: the autonomous map in one and two dimensions](@id sliced-autonomous-suite)

[`autonomous_map_experiment`](@ref). In ``d = 1``: 16 random markers on ``[-\pi, \pi) \times \mathbb{R}``
(``m = 1.25``), a frozen plan with ``B_v = B_x = 4``, multipliers that depend on the slot only,
``\mu_a \propto 0.25 + \exp(-v_j^2 / 2.5)``, and final weights from a kernel ratio limited to
``w_{\mathrm{ref}} \mu_a e^{\pm 0.25}``. In ``d = 2``: 12 random markers in 2D-2V, the per-axis
lift with two fixed index slices, one squeeze per axis with multipliers from ``|v|^2``, weights
``\propto \mu (1 + 0.15 \cos(\bar{x}_1 + \bar{x}_2))``.

| case | value | threshold |
|:--|--:|:--|
| ``d = 1``: symplecticity defect | ``1.0 \times 10^{-14}`` | ``< 10^{-13}`` |
| ``d = 1``: total momentum error | ``1.1 \times 10^{-16}`` | ``< 10^{-15}`` |
| ``d = 1``: lifts without the cotangent parts (control) | ``0.760`` | ``> 10^{-1}`` |
| ``d = 1``: multipliers ``\mu \propto 0.3 + e^{-\hat{x}^2/4 - (\hat{p}/L_p)^2/2}`` (control) | ``10.8`` | ``> 10^{-2}`` |
| ``d = 2``: symplecticity defect | ``6.7 \times 10^{-16}`` | ``< 10^{-13}`` |
| ``d = 2``: momentum error (largest component) | ``5.6 \times 10^{-17}`` | ``< 10^{-15}`` |

The first control drops ``(\nabla^2 S)^{-1}``: the positions follow the quantile map and the
momenta stay. The second lets the multipliers of the squeeze depend on the squeezed state; the
squeeze is then not the linear cotangent lift of the [proposition](@ref sliced-squeeze-structure).

## [Suite 4: composition with a later flow step](@id sliced-composition-suite)

[`composition_defect_experiment`](@ref): two markers (``m = 1.4``), the canonical shear
``\bar{p}_a = p_a - \cos x_a`` with the weights
``\bar{w} = (1 + 0.4 \bar{x}_1^2 + 0.2 \bar{p}_2, \ 1.2 + 0.3 \bar{x}_2 \bar{p}_1)``, then the
free flow for ``t = 0.7`` at frozen weights.

| quantity | value | threshold |
|:--|--:|:--|
| defect of the composition | ``0.0357`` | ``> 10^{-2}`` |
| distance to ``-\sum_a \mathrm{d}A_a \wedge \mathrm{d}\bar{w}_a`` | ``1.1 \times 10^{-16}`` | ``< 10^{-14}`` |
| defect with constant weights | ``4.3 \times 10^{-17}`` | ``< 10^{-15}`` |
| extended flow, ``\dot{\theta} = -\partial H/\partial\bar{w}`` | ``1.7 \times 10^{-16}`` | ``< 10^{-14}`` |
| extended flow, ``\dot{\theta} = +\partial H/\partial\bar{w}`` (control) | ``0.781`` | ``> 10^{-1}`` |

This checks the [composition defect](@ref sliced-composition-defect) with
``A_a = -t \, \bar{p}_a^2 / (2 m \bar{w}_a^2)``.
