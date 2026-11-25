# [Verification: Hamiltonian particle projection](@id projection-verification)

Every statement of the [Hamiltonian particle projection](@ref projection-index) has a test under
`test/projection/`, and each test runs the case that must fail. The pages
[Structural statements and target weights after a push](@ref projection-checks) and
[Hole filling](@ref projection-hole-filling) hold the tables of the measured numbers:
the structural statements, the hole-filling experiment and the target weights after a push. This
section says which file checks what. Run a test file with

```
julia --startup-file=no --project=test test/projection/<file>.jl
```

and a script with

```
julia --startup-file=no --project=scripts scripts/projection/<file>.jl
```

The files of `test/projection/` are in the group `core` of `test/runtests.jl`; the two files of
`test/projection/experiments/` are in the group `slow`.

## Test files of the group `core`

- **`test/projection/structure.jl`**: the statements of
  [Symplectic redistribution](@ref obstructions) that this package computes. The
  divergence of the effective velocity field against the
  [compressibility bracket](@ref weight-dependent-hamiltonians) for random ``\Psi(x, p, w)`` and
  ``W(x, p)`` in ``d = 1, 2, 3``; the determinant of the implicit weight map against the
  forward-mode Jacobian in ``d = 1, 2, 3``, and the Maxwellian
  [fold](@ref velocity-momentum-fold) with its sign change at ``|v| = v_{\mathrm{th}}``, two roots
  below ``p_{\max}`` and none above; ``\nabla_{z_a} H = \nabla \Psi_\mu(z_a)`` and
  [one map for all markers](@ref one-map-for-all-markers) for one step of symplectic Euler and of
  the implicit midpoint rule on 274 markers, with the area of the hole unchanged and no marker in
  its image ([holes persist](@ref holes-persist)). Cases that must fail: the bracket with the
  opposite sign and the divergence at frozen weight; ``f^d`` in place of ``f^{d-1}`` for
  ``d \geq 2``; two coincident markers of different weights under a weight-dependent Hamiltonian,
  which have no single map.
- **`test/projection/stages.jl`**: the [kick–drift map](@ref projection-hamiltonian) in ``d = 1``
  and ``d = 2`` (symplectic, momentum-conserving, the momentum correction a Galilean boost, a kick
  and a drift one step of symplectic Euler of ``H^K + H^D``); the
  [stage Jacobians](@ref projection-stage-jacobians) in ``(x, v)`` for ``d = 1``; the
  [gauge step](@ref gauge-step) at fixed ``(x, p)`` in ``d = 1`` and ``d = 2`` (the velocity rule
  and ``P``); the [non-separable stages](@ref projection-nonseparable) in ``d = 1`` and ``d = 2``
  (the velocity form against the implicit midpoint rule of ``H^N_s``, symplecticity, ``P``, the
  Jacobian determinant and the divergence in ``(x, v)``) and the ``\pm v`` coupling at
  ``\partial_x W = 0``; the guards of the stage solvers (an unknown rule, a step that does not
  converge). Cases that must fail: the naive kick; parameters that depend on the state
  ([parameters that depend on the state](@ref state-dependent-parameters)); the trapezoidal rule;
  the non-separable stage without the centroid shift; the first-order determinant and the
  determinant with ``(I - \tfrac12 \mathbb{J} K)^{-1}``.
- **`test/projection/charge.jl`**: the [charge-conserving reweighting](@ref projection-charge) on
  400 random markers in ``d = 1`` after a kick and a drift, with 16 shape functions: the grid
  charge, positivity, the total weight and ``P``; the minimal relative entropy, checked on the
  null space of the charge matrix independently of the dual problem; the identity when the weights
  already have the charge; the method type
  [`ChargeConservingReweighting`](@ref ParticleRemapping.Projection.ChargeConservingReweighting).
  Cases that must fail: the target weights alone miss the charge; the additive least-norm
  correction gives negative weights; after positions that leave a shape function with charge
  without markers, no solution exists, and `remap` reports it.
- **`test/projection/lowrank.jl`**: the [BUG integrator](@ref projection-bug) on
  ``64 \times 256`` points up to ``T = 3``: free streaming from rank four keeps the numerical rank
  3 and converges with the step; the full grid against the characteristics for
  ``E = -0.5 \sin x``; the BUG solution against the full grid at full rank, of second order in the
  step; the truncation error below ``n_{\mathrm{steps}} \vartheta`` at ``\vartheta = 10^{-4}`` and
  ``10^{-6}``; the evaluation off the grid; the numerical rank. Cases that must fail: free
  streaming from rank one stays at rank one; the fixed rank one with the field.
- **`test/projection/targets.jl`**: the lumped [projection](@ref projection-projection) on hat
  functions (nonnegativity, the partition of unity and the normalisation, interpolation); the
  five-point Laplacian, whose null space is the constants; the
  [target points](@ref projection-targets) of the hole-filling problem (the total weight is exact,
  the density deviation drops, the targets fill the hole, the target weights are ``c f`` with
  ``c`` close to ``|\Omega| / N``) and the flow history; the target weights from the lumped
  projection, whose relative error exceeds ``0.1``.
- **`test/projection/parametrisation.jl`**: the parametrised maps on a small lattice: the
  parameter counts and the guard of an unknown variant; ``\theta = 0`` is the identity; for random
  ``\theta``, every parametrisation (weight-dependent, controls A, B and C, non-separable) moves
  the markers, conserves ``P`` and is symplectic; the controls B and C move two coincident markers
  of different weights to one point, one map of ``(x, v)``, which the weight-dependent and the
  non-separable maps do not; the history of the kick–drift map; the
  sensitivities of the non-separable map against forward-mode derivatives, and a stage beyond the
  iteration limit returns `NaN`; the [fit](@ref projection-fit) lowers the objective. Case that
  must fail: the sensitivities without the terms of the centroid shift and the boost.
- **`test/projection/diagnostics.jl`**: the diagnostics of the experiments against closed forms:
  the largest and the local empty disk of a perfect lattice (``\sqrt{2}/2``) and of the hole of
  the experiment (24 markers removed, ``2.75`` lattice units); the hole, its outline and the
  fill; the density deviation of a lattice that the hats resolve; the weight errors by velocity
  band.
- **`test/projection/remap.jl`**: the whole algorithm through the method types
  [`KickDriftRemap`](@ref ParticleRemapping.Projection.KickDriftRemap) and
  [`NonSeparableRemap`](@ref ParticleRemapping.Projection.NonSeparableRemap) on a small lattice
  with a hole: the objective drops, ``P`` and the total weight are kept, the result carries the
  positions modulo ``L``, the target weights and the momenta of the map; the option `charge`; the
  default target weights from the lumped projection; the guards of nonpositive target weights, of
  an unknown variant and of markers that are not 1D-1V. Case that must fail: the target weights
  alone miss the grid charge.
- **`test/projection/frames.jl`**: the states and the frames of the remaps on small runs
  ([`remap_states`](@ref ParticleRemapping.remap_states),
  [`remap_frames`](@ref ParticleRemapping.remap_frames)) for the kick–drift map, the
  non-separable map, a kick–drift map of three stages with the charge-conserving reweighting, and
  [`ChargeConservingReweighting`](@ref ParticleRemapping.Projection.ChargeConservingReweighting)
  alone: the labels; the first state is the input and the last one `result.markers`; every stage
  of the canonical map keeps the old weights and ``P``, a kick keeps the positions and a drift the
  velocities, and the stages end at the image of the fitted map; each gauge step keeps
  ``(x, p)``; the frames pass through every state, keep the weights along the canonical stages and
  ``(x, p)`` along the gauge steps. The remaps of the hole-filling experiment
  ([`experiment_remap`](@ref ParticleRemapping.Projection.experiment_remap)) end at its velocities
  and charge-conserving weights; the history of the density-interpolation flow and the frames of
  the push ([`moser_flow_history`](@ref ParticleRemapping.Projection.moser_flow_history),
  [`dlra_push_frames`](@ref ParticleRemapping.Projection.dlra_push_frames)). Case that must fail:
  markers that are not the input of the remap.

## Test files of the group `slow`

- **`test/projection/experiments/hole_filling.jl`** (about 18 minutes): the
  [hole-filling experiment](@ref projection-hole-filling) at ``x_h = \pi/2`` and ``x_h = \pi``.
  It checks the set-up (2280 markers, 24 removed, 2304 and 3264 parameters), ``P`` and the total
  weight, the charge-conserving reweighting (at most 16 Newton iterations), every number of the
  two tables with the digits it is quoted with, the return code of every fit, the sensitivities
  of the non-separable map at the fitted ``\theta`` with the active iteration limit of 50, the
  result of `remap` against the experiment, and the relative error of the projected target
  weights in the four velocity bands. At ``x_h = \pi`` it checks that the fit of the kicks and
  drifts from ``\theta = 0`` ends on the mirror-symmetric saddle (`SUCCESS`, stationary, and the
  stationarity ratio ``\|J^{\mathsf{T}} F\| / (\|J\|_F \|F\|)`` of the fit below ``10^{-8}``),
  that the fit from the seeded start 101 fills the hole (`MAXITERS` after 40 iterations, the
  empty disk below ``1.55`` and below ``0.9`` of the reweighting alone, the fill above ``0.85``,
  the empty disk of the mirror ellipse at ``v = -1.2`` below ``1.0``, the objective between
  ``0.25`` and ``0.40``), and that control A and control B from the same seeded start give every
  number of their rows "seeded" within the quoted range, widened on each side by its width.
  Cases that must fail: controls B and C, one area-preserving map of ``(x, v)``, fitted from
  ``\theta = 0`` do not shrink the hole beyond the reweighting alone; the kicks and drifts at the
  saddle do not fill the hole at ``x_h = \pi``; the target weights alone miss the grid charge;
  the sensitivities without the terms of the shift and the boost miss the derivatives.
- **`test/projection/experiments/dlra_weights.jl`** (about 3 minutes): the
  [target weights after a push](@ref projection-push-test). It checks the push (largest speed,
  empty disk, fill and density deviation of the pushed markers), that no target weight of any
  source is negative, ``P`` and the grid charge after the correction, every number of the two
  tables with the digits it is quoted with, and that the weight error falls with ``\vartheta`` in
  every band. Case that must fail: the fixed rank one misses the sheared density.

## Scripts

- **`scripts/projection/checks.jl`** (about one minute): runs `structure.jl`, `stages.jl`,
  `charge.jl` and `lowrank.jl` of `test/projection/`, each in a module of its own, and prints the
  measured numbers of the table of structural statements.
- **`scripts/projection/hole_filling.jl`** (about 31 minutes): prints the two tables of
  the hole-filling experiment with the return code of each fit, the cost of each map, the
  sensitivities of the non-separable map and the relative error of the projected target weights,
  and checks every number against the value of the
  [hole-filling page](@ref projection-hole-filling). It runs the experiment again with the
  targets moved by ``10^{-13}`` relative and checks every number there too, and runs the seeded
  fits of the weight-dependent map, control A and control B at ``x_h = \pi`` for the seeds 101
  to 105 and checks that each of their numbers lies in its quoted range.
- **`scripts/projection/dlra_weights.jl`** (about 3 minutes): prints the two tables of the target
  weights after a push and checks every number against the value of the source note.
- **`scripts/projection/higher_dimensions.jl`** (seconds): prints the numbers of parameters,
  nodes and basis functions of the [higher dimensions](@ref projection-higher-dimensions) for
  ``d = 1, 2, 3`` and checks the column ``d = 1`` against the code.
- **`scripts/projection/make_figures.jl`** (about ten minutes): writes the figures of the chapter
  as PDF files and the five animations as MP4 files to `results/projection/`, with the standard
  pictures of the [example](@ref projection-example) (`projection-remap.pdf`,
  `projection-stages.pdf`, `projection-stages.mp4`), and prints the changes of the conserved
  quantities of the example.

## The numbers of the source note

The numbers of the note *Hamiltonian particle projection: a symplectic redistribution of weighted
markers* are the numbers of the tables on
[Structural statements and target weights after a push](@ref projection-checks). The tables on
[Hole filling](@ref projection-hole-filling) differ from the note in value: the note's fit is
Levenberg–Marquardt, and the tables give the Gauss–Newton fit with the return code of each row,
the mirror-symmetric saddle that it reaches from ``\theta = 0`` at ``x_h = \pi``, and the ranges of
the seeded fits there.
The two experiment scripts hold them in their tables `NOTE_*` and `NOTE`, and the two experiment
tests in `TABLES`, `WEIGHT_TABLE` and `REMAP_TABLE`; a measured number agrees with a quoted one
when it lies within half a unit of the last quoted digit. A quoted range holds a number of the
script between its ends, and a number of the test within the range widened on each side by its
width. The experiment tests check every number of the tables in this way. Three points differ in
form, not in value:

- The note quotes the relative error of the projected target weights at ``x_h = \pi`` only. At
  ``x_h = \pi/2`` the bands ``[1, 2)`` and ``[3, 4)`` differ in the second digit (``0.18`` and
  ``0.22`` against ``0.17`` and ``0.21``); the [hole-filling page](@ref projection-hole-filling) gives both.
- The tests assert bounds where the note quotes a measured value: ``|\Delta P|`` below
  ``10^{-15}`` in the hole-filling test and below ``10^{-12}`` in the test of the target weights
  (the note quotes ``1.1 \cdot 10^{-17}`` and ``3.5 \cdot 10^{-17}`` relative), the
  error of the sensitivities below ``10^{-10}`` (the note quotes ``2.8 \cdot 10^{-13}``), and the
  cases that must fail above a threshold below the quoted value (for example ``> 0.28`` for the
  charge error ``0.29``, ``> 0.15`` for the rank-one error ``0.16``). `checks.jl` and the
  experiment scripts print the measured values.
- The cost of the maps and of the fit (seconds) is printed by `hole_filling.jl` and
  `dlra_weights.jl` and depends on the machine; no test checks it.
