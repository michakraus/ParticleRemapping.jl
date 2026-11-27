# [Verification: sliced redistribution](@id sliced-verification)

Every number of the [sliced redistribution](@ref sliced-overview) comes from a test in
`test/sliced/`, a script in `scripts/sliced/` or an `@example` block of its pages. The functions
of [`SlicedRedistribution`](@ref ParticleRemapping.SlicedRedistribution) port the source
verification script of the construction, `verify_symplectic_redistribution.jl` of the sliced
construction in the particle-in-cell project; the package replaces it there, and the project's
note `Sliced-Symplectic-Particle-Redistribution.tex` is written from this chapter. The derivation
that produced the script is the source of the observations marked *derivation* on the page
[Constructions that fail](@ref sliced-failures).
The tables of the suites are on the [structural suites](@ref sliced-suites),
[hole filling](@ref sliced-holes) and [long-time run](@ref sliced-long-time-page) pages.

## Tests

All eleven files run in the `core` group of `test/runtests.jl`.

| file | what it checks | controls and error checks |
|:--|:--|:--|
| `test/sliced/binning.jl` | [`uniform_quantile_targets`](@ref ParticleRemapping.SlicedRedistribution.uniform_quantile_targets): the midpoints of ``n`` equal cells, symmetric, with the spacing ``2L/n``; [`sorted_slices`](@ref ParticleRemapping.SlicedRedistribution.sorted_slices): equal slices that partition the markers, every key of a slice below every key of the next, each slice ordered by the second key; [`slice_anchors`](@ref ParticleRemapping.SlicedRedistribution.slice_anchors): the centre of each marker's slice | ``N = 24`` markers in 5 slices throw an `ArgumentError` |
| `test/sliced/lift.jl` | for eight combinations of ``\beta_{\mathrm{mean}}``, ``\beta_{\mathrm{intra}}`` and ``\varepsilon``: the lift is the ForwardDiff gradient of the potential, the Hessian its ForwardDiff Hessian; ``\nabla^2 S \, \mathbf{1} = \mathbf{1}``, ``\nabla^2 S \succeq (1 - \max\beta) I``, slice-constant vectors with zero mean have the eigenvalue ``1 - \beta_{\mathrm{mean}}``; the cotangent part solves ``\nabla^2 S \, \tilde{v} = v``; both sums are kept; translation equivariance. The lift and the dual lift are symplectic with and without anchors; the anchored lift keeps ``\sum_a v_a``. The quantile limit ``\tilde{u}_a \to \bar{u} + \tau_r`` for ``\beta = 1``, ``\varepsilon = 10^{-9}`` | the lift without its cotangent part has the defect ``> 10^{-1}`` and leaves ``v`` unchanged; a mutant that keeps each slice mean instead of the global mean misses the quantile limit by ``> 10^{-3}`` |
| `test/sliced/squeeze.jl` | [`canonical_weight_squeeze`](@ref ParticleRemapping.SlicedRedistribution.canonical_weight_squeeze) against the dense ``M`` and ``M^{-\top}``; ``\det M = \prod_a \mu_a``; the total momentum; the common shift ``\frac{1}{N}\sum_c (1 - \mu_c) p_c``; the ForwardDiff Jacobian is symplectic and equals the block matrix of the linear cotangent lift | ``\sum_a \mu_a \ne N`` throws an `ArgumentError`; for ``\sum_a \mu_a = N(1 + \delta)`` the matrix ``M`` changes the total momentum by ``\delta \sum_c (1 - \mu_c) p_c``, more than ``10^{-3}`` |
| `test/sliced/reweighting.jl` | [`kernel_density_ratio`](@ref ParticleRemapping.SlicedRedistribution.kernel_density_ratio): one source and one evaluation point give ``w K`` with the periodic kernel; evaluation at equal-weight sources gives the weight; periodicity in ``x``; a ``1 \times N`` matrix gives the vector result; in 2D the kernel is the product over the axes. [`squeeze_multipliers`](@ref ParticleRemapping.SlicedRedistribution.squeeze_multipliers): ``\sum_a \mu_a = N`` and ``\mu > 0`` after one and three passes; one pass is the normalised density ratio plus the floor | the corrector passes change the multipliers |
| `test/sliced/diagnostics.jl` | [`periodic_wrap`](@ref ParticleRemapping.SlicedRedistribution.periodic_wrap), [`phase_quadrant_counts`](@ref ParticleRemapping.SlicedRedistribution.phase_quadrant_counts), [`phase_cell_statistics`](@ref ParticleRemapping.SlicedRedistribution.phase_cell_statistics) (one marker per cell, a moved, a wrapped and a dropped marker, the relative standard deviation), [`cell_count_history`](@ref ParticleRemapping.SlicedRedistribution.cell_count_history), [`smoothed_density_difference`](@ref ParticleRemapping.SlicedRedistribution.smoothed_density_difference) (zero for equal sets, 1 for doubled weights, periodic in ``x``) | |
| `test/sliced/map.jl` | the defaults of [`SlicedQuantileRemap`](@ref ParticleRemapping.SlicedRedistribution.SlicedQuantileRemap); the plan of the 120-marker stages example (5 velocity slices of 24, 12 position slices of 10, their orders, anchors at the slice centres, targets at the uniform quantiles, mollifier widths of two target spacings); the total weight and momentum, ``\bar{w} = w_{\mathrm{ref}} \mu > 0``, the field `q̄` holds the output positions ``\bar{x}``; the frozen map (both lifts and the squeeze at the plan and ``\mu`` of the remap) reproduces the output and is symplectic | the frozen map without the cotangent parts has the defect ``> 10^{-1}`` |
| `test/sliced/peraxis.jl` | [`per_axis_sliced_lift`](@ref ParticleRemapping.SlicedRedistribution.per_axis_sliced_lift) equals the lift of each axis; the remap of 36 markers in 2D-2V (``B_v = 3``, ``B_x = 6``) keeps the total weight and the total momentum, has positive weights and slices each axis by its own velocity component; the frozen per-axis map is symplectic; on a checkerboard in ``(x_1, x_2)`` with even marginals ``(x_\ell, v_\ell)`` (360 markers) the per-axis lifts leave fewer than 36 markers in the empty quadrants ([a hole between the axes](@ref sliced-hole-between-axes)) | a slice list of the wrong length throws a `DimensionMismatch`; the same lifts fill the checkerboard in ``(x_1, v_1)``, more than 72 markers in every quadrant |
| `test/sliced/frames.jl` | [`remap_states`](@ref ParticleRemapping.remap_states) of the standard example: the five labels, the first state is the input, the last the result with its positions in ``[-\pi, \pi)``; the canonical stages keep the old weights, carry the momenta `p̃`, `p̂`, `p̄` and the total momentum; the gauge step keeps ``(x, p)`` and sets ``w_{\mathrm{ref}} \mu``. [`remap_frames`](@ref ParticleRemapping.remap_frames): the number of frames, every state held and reached exactly, the frame halfway through ``\Psi_q`` is the lift with halved ``\beta``, every canonical frame keeps the old weights and the total momentum, the gauge frames keep ``(\bar{x}, \bar{p})`` and move the weights linearly; the states and frames of 2D-2V markers | the velocities after the squeeze at the old weights differ from the output by more than 0.1; the straight line between two states misses the frame of the lift; a result of another number of markers throws a `DimensionMismatch` |
| `test/sliced/forwarddiff.jl` | the port against the source script; see below | |
| `test/sliced/experiments.jl` | suites 1–5a of the [structural suites](@ref sliced-suites) and [hole filling](@ref sliced-holes) pages with the thresholds of the source script (in suite 5a also the total momentum and the unchanged weights), and the probes of [Constructions that fail](@ref sliced-failures): the fold of the repulsion with a wall (first ``\lambda_{\min} < 0`` at pass 3, frame 1; first non-finite Jacobian at pass 4, frame 3; no fold for ``k_w = 0``), the rotated-frame lift (quadrants ``(198, 2, 200, 0)``), the checkerboard reweighted at fixed momenta (smallest weight ratio below 0.2, largest speed above 5, against at most 1.05 after the sliced lifts) | every control of suites 1–4 (see below); the one-slice lifts of suite 5a keep the quadrants ``(200, 0, 200, 0)``; an unconverged midpoint step (3 fixed-point sweeps) has the defect ``> 10^{-6}`` |
| `test/sliced/experiments_long_time.jl` | suite 5b: the thresholds of the source script, the cell counts before and after, the relative standard deviations and the ``L^2`` difference to round-off, the Hessian bounds ``1 - \max\beta`` | without the squeeze, cells empty again, the relative standard deviation exceeds that before the remap and the ``L^2`` difference triples; with ``\beta_{\mathrm{mean}} = \beta_{\mathrm{intra}} = 0.88`` both Hessians have the smallest eigenvalue 0.12; with ``\varepsilon`` = 0.1 target spacings the largest ``\lvert\hat{x} - x^{\mathrm{anc}}\rvert`` exceeds 1.4 times that of the defaults |

## Scripts

| script | what it does |
|:--|:--|
| `scripts/sliced/verify_sliced_redistribution.jl` | runs suites 1–5 and the probes of the constructions that do not work, prints every number and checks it against the threshold of the source script; it also times suite 5b, the run time quoted on the [algorithm](@ref sliced-algorithm) page |
| `scripts/sliced/long_time_experiment.jl` | the variants of the long-time remap (no squeeze, ``\beta = 0.88``, three mollifier widths) and the run with a remap after every 100 steps; prints the tables of [Suite 5b](@ref sliced-long-time) and the mean relative standard deviation for ``t \ge 5``, and checks the claims of the pages: 14 empty cells before and none after the remap, the control without the squeeze, the eigenvalue 0.12, the larger excursion for the narrow mollifier, the weights of the run without remaps unchanged and the total weight of the run with remaps kept to ``10^{-13}`` |
| `scripts/sliced/make_figures.jl` | writes the figures (PDF) and animations (MP4) of the chapter to `results/sliced/` or to `--results-dir`: the standard pictures of the standard example, `sliced-remap.pdf` (S1), `sliced-stages.pdf` (S2) and `sliced-stages.mp4` (A1), the stages coloured by slice, `sliced-slices.pdf` and `sliced-slices.mp4`, and `checkerboard.pdf`, `wave-remap.pdf`, `cell-count-history.pdf`, `long-time-run.mp4`, `long-time-run-no-remap.mp4`; it checks nothing |

## Controls that must fail

The tables of the [structural suites](@ref sliced-suites), [hole filling](@ref sliced-holes) and
[long-time run](@ref sliced-long-time-page) pages give the values; each control must
exceed its threshold.

- Suite 1: the velocity split without the recoil, ``\bar{v} = \bar{p}/(m w)`` (``> 10^{-1}``),
  and a canonical map that reads ``w`` across ``\mathrm{d}w`` (``> 10^{-1}``).
- Suite 2: the flow step without the recoil (``> 10^{-2}``), and the step with three fixed-point
  sweeps only.
- Suite 3: the lifts without their cotangent parts (``> 10^{-1}``), and multipliers that depend
  on the squeezed state (``> 10^{-2}``).
- Suite 4: the composition itself (defect ``> 10^{-2}``), and the extended flow with
  ``\dot{\theta} = +\partial H / \partial \bar{w}`` (``> 10^{-1}``).
- Suite 5a: one slice for each lift (``Q_2 + Q_4 < 20``).
- Suite 5b: the remap without the squeeze (cells empty again). The variants with
  ``\beta = 0.88`` and with a mollifier of 0.1 target spacings must show the soft Hessian and the
  larger excursion from the anchors.

## Comparison with the source script

`test/sliced/forwarddiff.jl` holds, in a module `ScriptReference`, the source script's code
verbatim: its hand-written forward-mode `Dual` number, its linear solve, its Jacobian
`jacobian_dual`, its canonical matrix, its sliced lift and its weight squeeze. At the plan and the
multipliers of suite 3 (16 markers) and a random state it compares

- the values of the port (both lifts and the squeeze) with those of the script, and the lift with
  and without anchors, to the relative tolerance ``8 \varepsilon_{\mathrm{mach}}``;
- the ForwardDiff Jacobians of the port with the script's `Dual` Jacobians, of the composed map
  and of the first lift, to ``8 \varepsilon_{\mathrm{mach}} N \operatorname{cond}(J) \lVert J \rVert``.
  The two differ by the pivoting of the linear solves only.

The thresholds of `test/sliced/experiments.jl`, `test/sliced/experiments_long_time.jl` and
`scripts/sliced/verify_sliced_redistribution.jl` are those of the source script. The reference
numbers in the two test files come from a cold run (Julia 1.13.4, aarch64) of the working copy of
the source script, with its own `Dual` numbers; the copy in the project does not parse (see the
[findings](@ref verification-findings)). The defects are round-off and differ from those of the script
in the last digits; the controls agree with the script to the round-off of the Jacobians
(``\operatorname{cond}(\nabla^2 S) \lesssim 10``), so the tests require a relative agreement of
``10^{-13}``.

| quantity | source script (cold run) | the package |
|:--|--:|:--|
| suite 1, no recoil | 0.5073615519131363 | rtol ``10^{-13}`` |
| suite 1, across ``\mathrm{d}w`` | 0.9740596015484777 | rtol ``10^{-13}`` |
| suite 2, no recoil | 0.030573267254231595 | rtol ``10^{-13}`` |
| suite 3, no cotangent parts | 0.7600914624090467 | rtol ``10^{-13}`` |
| suite 3, multipliers from the state | 10.786509290207738 | rtol ``10^{-13}`` |
| suite 4, defect of the composition | 0.03572970047023824 | rtol ``10^{-13}`` |
| suite 5a, quadrants after the sliced lifts | ``(98, 101, 112, 89)`` | equal |
| suite 5a, quadrants after the one-slice lifts | ``(200, 0, 200, 0)`` | equal |
| suite 5b, empty, fewest, most per cell before; after | ``(14, 0, 20)``; ``(0, 1, 15)`` | equal |
| suite 5b, relative std before | 0.5976086723909207 | rtol ``8 \varepsilon_{\mathrm{mach}}`` |
| suite 5b, relative std after | 0.367197786605922 | rtol ``8 \varepsilon_{\mathrm{mach}} \cdot 900`` |
| suite 5b, ``L^2`` difference | 0.10414395551641199 | rtol ``8 \varepsilon_{\mathrm{mach}} \cdot 900`` |

Two probes reproduce observations of the derivation: [`rotated_frame_lift_probe`](@ref ParticleRemapping.SlicedRedistribution.rotated_frame_lift_probe)
gives the quadrant counts ``(198, 2, 200, 0)`` that the derivation reports, and
[`repulsion_wall_probe`](@ref ParticleRemapping.SlicedRedistribution.repulsion_wall_probe) runs with
the derivation's parameters into a fold and then into non-finite numbers, where the derivation's
run stops with a linear solve on non-finite entries. The other observations marked *derivation*
on [Constructions that fail](@ref sliced-failures) have no reproduction in the package.
