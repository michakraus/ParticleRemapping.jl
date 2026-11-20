# [Experiments](@id sliced-experiments)

```@meta
CurrentModule = ParticleRemapping.SlicedRedistribution
```

Five suites verify the [symplecticity condition](@ref symplecticity-conditions) and the
statements of [the stages and their structure](@ref sliced-structure), and test the
redistribution on two phase-space distributions. The symplecticity defect of a map ``F`` of the
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

## [Suite 5a: hole filling on a checkerboard](@id sliced-checkerboard)

[`checkerboard_experiment`](@ref): 400 markers of weight ``1/400`` (``m = 1``), uniform in the
quadrants ``Q_1 = [0, 1)^2`` and ``Q_3 = [-1, 0)^2``; ``Q_2`` and ``Q_4`` are empty, and both
marginals are uniform on ``[-1, 1]``. Only the lifts act, the weights do not change.

| map | ``Q_1`` | ``Q_2`` | ``Q_3`` | ``Q_4`` | largest ``\lvert x \rvert`` | largest ``\lvert v \rvert`` |
|:--|--:|--:|--:|--:|--:|--:|
| before | 200 | 0 | 200 | 0 | 1.00 | 1.00 |
| sliced lifts, ``B_v = 10``, ``B_x = 20`` | 98 | 101 | 112 | 89 | 1.28 | 1.04 |
| one slice for each lift (control) | 200 | 0 | 200 | 0 | 6.77 | 1.79 |

The thresholds are ``Q_2 + Q_4 < 20`` for the control and ``65 \le Q_k \le 135`` for the sliced
lifts. Within each velocity slice the positions lie on one side of ``x = 0``; ``\Psi_q`` spreads
them over ``[-1, 1]`` and so fills ``Q_2`` and ``Q_4``. The one-slice lifts see the uniform
marginals and leave the markers in their quadrants. They also stretch the conjugate variables:
without anchors the cotangent part ``(\nabla^2 S)^{-1}`` amplifies them by up to
``1 / (1 - \beta_{\mathrm{intra}}) = 6.7`` where the Laplacian of the potential is weak.

```@example sliced-experiments
using ParticleRemapping, CairoMakie
plot_checkerboard_filling(checkerboard_experiment())
```

## [Suite 5b: long-time run in an electric field](@id sliced-long-time)

[`long_time_remap_experiment`](@ref): ``30 \times 30`` markers at the cell centres of
``[-\pi, \pi) \times [-1.4, 1.4]`` with weights ``\propto f_0(x, v)``,
``f_0 = (1 + 0.35 \cos x)(0.40 + v^2) e^{-v^2/1.4}`` (``\sum w = 1``, ``m = 1.2``), pushed by 500
drift–kick–drift steps of ``\Delta t = 0.05`` through the field
``E(t, x) = -0.85 \sin x - 0.75 \sin(2x - 1.3 t)`` to ``T = 25``
([`electric_wave_markers`](@ref), [`electric_wave_simulation`](@ref)), then one
[`remap`](@ref) with the default [`SlicedQuantileRemap`](@ref). The statistics count the markers
in the ``10 \times 10`` cells of ``[-\pi, \pi) \times [-2.25, 2.25]``
([`phase_cell_statistics`](@ref)); the ``L^2`` difference compares kernel-smoothed densities on a
``30 \times 30`` grid ([`smoothed_density_difference`](@ref), ``h_x = 0.45``, ``h_v = 0.35``).

| quantity | before | after | threshold |
|:--|--:|--:|:--|
| empty cells | 14 | 0 | before ``\ge 10``, after ``= 0`` |
| fewest, most markers in a cell | 0, 20 | 1, 15 | |
| relative std of the cell counts | 0.598 | 0.367 | after ``< 0.65 \times`` before |
| relative ``L^2`` difference of the smoothed densities | | 0.104 | ``< 0.12`` |
| total momentum error | | ``1.1 \times 10^{-16}`` | ``< 10^{-15}`` |
| total weight error | | 0 | ``< 10^{-15}`` |
| smallest eigenvalue of ``\nabla^2 S_1``, ``\nabla^2 S_2`` | | 0.2226, 0.2187 | ``\ge 0.22``, ``\ge 0.20`` |

The multipliers of this remap depend on the state; the map is canonical for the frozen
multipliers ([theorem](@ref sliced-composed-map)).

```@example sliced-experiments
ex = long_time_remap_experiment()
plot_wave_remap(ex)
```

### Variants

`scripts/sliced/long_time_experiment.jl` runs the same remap with other parameters. Columns:
empty cells, fewest and most markers per cell and relative std after the remap; the ``L^2``
difference; the smallest eigenvalues ``\lambda_q``, ``\lambda_p`` of the two Hessians; the largest
distance of ``\hat{x}`` and of ``\bar{x}`` from the anchors (the fields `q̂`, `q̄` and `q_anc`);
the largest ``|\bar{v}|``.

| variant | empty | min | max | rstd | ``L^2`` | ``\lambda_q`` | ``\lambda_p`` | ``\lvert\hat{x} - x^{\mathrm{anc}}\rvert`` | ``\lvert\bar{x} - x^{\mathrm{anc}}\rvert`` | ``\lvert\bar{v}\rvert`` |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| defaults | 0 | 1 | 15 | 0.367 | 0.104 | 0.2226 | 0.2187 | 0.361 | 3.46 | 2.87 |
| no squeeze (control) | 3 | 0 | 18 | 0.740 | 0.337 | 0.2226 | 0.2187 | 0.361 | 3.46 | 31.9 |
| ``\beta_{\mathrm{mean}} = \beta_{\mathrm{intra}} = 0.88`` | 0 | 3 | 18 | 0.337 | 0.096 | 0.1200 | 0.1200 | 0.921 | 13.2 | 2.88 |
| ``\varepsilon`` = 0.1 spacings | 1 | 0 | 17 | 0.398 | 0.097 | 0.2200 | 0.2000 | 0.559 | 2.48 | 2.74 |
| ``\varepsilon`` = 0.5 spacings | 1 | 0 | 16 | 0.369 | 0.096 | 0.2202 | 0.2006 | 0.438 | 1.98 | 2.72 |
| ``\varepsilon`` = 1 spacing | 0 | 1 | 16 | 0.353 | 0.097 | 0.2206 | 0.2034 | 0.374 | 2.48 | 2.77 |

In the control without the squeeze the weights ``w_{\mathrm{ref}} \mu`` are set at the fixed
momenta ``\hat{p}``, and ``\bar{v} = \hat{p}/(m\bar{w})`` stretches where the weights are small
(the [failures](@ref sliced-failures) page). The other variants keep every cell occupied or leave
one empty; they move the positions further from the anchors (``\beta = 0.88``) or the lifted
positions (small ``\varepsilon``). The squeeze moves some markers by up to 3.46 from their anchors
in the default run: the slots with small ``\mu_a``.

### Remaps during the run

The same run with a remap after every 100 steps (``t = 5, 10, 15, 20``), statistics every 10
steps. At ``t = 0`` the lattice covers ``|v| \le 1.4`` only, hence the empty cells of the
statistics box ``|v| \le 2.25``.

| ``t`` | empty cells, no remap | rstd, no remap | empty cells, remaps | rstd, remaps |
|--:|--:|--:|--:|--:|
| 0 | 20 | 0.699 | 20 | 0.699 |
| 2.5 | 16 | 0.664 | 16 | 0.664 |
| 5 | 17 | 0.663 | 0 | 0.389 |
| 10 | 17 | 0.592 | 0 | 0.304 |
| 15 | 13 | 0.569 | 0 | 0.349 |
| 20 | 13 | 0.600 | 0 | 0.315 |
| 25 | 14 | 0.598 | 0 | 0.427 |

The mean relative standard deviation for ``t \ge 5`` is 0.587 without and 0.365 with the remaps
(`scripts/sliced/long_time_experiment.jl`).

```@example sliced-experiments
initial = electric_wave_markers()
plain = electric_wave_simulation(initial; record_every = 10)
remapped = electric_wave_simulation(initial; record_every = 10,
    remap_method = SlicedQuantileRemap(), remap_every = 100)
h_plain = cell_count_history(plain.snapshots, plain.times)
h_remap = cell_count_history(remapped.snapshots, remapped.times)
plot_cell_count_history([h_plain, h_remap]; labels = ["no remap", "remap every 100 steps"],
    remap_times = remapped.remap_times)
```

The animations show the two runs (colour: weight; below: empty cells and ten times the relative
standard deviation).

```@example sliced-experiments
animate_long_time_run("long-time-run.mp4", remapped, h_remap)
animate_long_time_run("long-time-run-no-remap.mp4", plain, h_plain)
nothing # hide
```

With a remap after every 100 steps:

```@raw html
<video src="long-time-run.mp4" controls loop muted playsinline style="width:100%"></video>
```

Without remaps:

```@raw html
<video src="long-time-run-no-remap.mp4" controls loop muted playsinline style="width:100%"></video>
```
