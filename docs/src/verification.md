# [Verification](@id verification)

```@meta
CurrentModule = ParticleRemapping
```

Every statement and every number of this manual is checked by a test file of the package, by a
script that prints it next to the value of its source, or by an `@example` block of its page.
Each check of a statement runs with a control that must fail, and the control asserts the margin
by which it fails. This chapter says which file checks what. Its pages follow the chapters of the
manual:

- [Gauge rezoning and canonical transport](@ref remap-verification),
- [Hamiltonian particle projection](@ref projection-verification),
- [Spline-projection redistribution](@ref rd-verification),
- [Sliced redistribution](@ref sliced-verification).

The last section of this page lists the findings against the four sources: the statements and
numbers of the notes and scripts that the package does not reproduce.

## Running the checks

From the package root:

```
julia --startup-file=no --project=test test/runtests.jl              # all test files
julia --startup-file=no --project=test test/runtests.jl core         # without the slow files
julia --startup-file=no --project=test test/runtests.jl slow         # the two slow files
julia --startup-file=no --project=test test/<construction>/<file>.jl # one test file
julia --startup-file=no --project=scripts scripts/<construction>/<file>.jl
```

`test/runtests.jl` runs every file in a `@safetestset` of its own. The group `slow` holds the two
experiments of the Hamiltonian particle projection, `test/projection/experiments/hole_filling.jl`
and `test/projection/experiments/dlra_weights.jl`, of about 18 and 3 minutes; the group `core`
holds every other file and runs in about 5 minutes.

A cold run of both groups with Julia 1.13.4 (aarch64) passes every check:

| part | test files | checks |
|:--|--:|--:|
| quality: Aqua, explicit imports | 2 | 13 |
| common layer | 6 | 136 |
| gauge rezoning and canonical transport | 7 | 512 |
| Hamiltonian particle projection, group `core` | 9 | 494 |
| Hamiltonian particle projection, group `slow` | 2 | 309 |
| spline-projection redistribution | 8 | 190 |
| sliced redistribution | 11 | 315 |
| **total** | **45** | **1969** |

## The common layer and the foundations

- **`test/common/markers.jl`**: the constructors of [`Markers`](@ref) and their dimension
  errors; the momenta, totals and kinetic energy for ``d = 1`` and ``d = 2``;
  [`from_canonical`](@ref) as the inverse of [`canonical`](@ref); the [`gauge_step`](@ref): the
  positions and the momenta are kept, the total weight is not, and a non-positive stretch throws.
- **`test/common/symplectic.jl`**: [`pack`](@ref), [`unpack`](@ref) and [`canonical_J`](@ref);
  the [`symplecticity_defect`](@ref) of a coupled kick (round-off) and of a scaling (the control,
  defect 1); [`fd_jacobian`](@ref) of orders 2 and 4 against an exact Jacobian; the
  [`presymplectic_matrix`](@ref), its skew symmetry, its rank and its kernel; the
  [`pullback_residual`](@ref) of a state-dependent gauge step (below ``10^{-8}``) and of free
  streaming (the control, above ``0.1``).
- **`test/common/grids.jl`**: [`PeriodicGrid`](@ref), [`cellof`](@ref) with periodic wrapping,
  the hat functions [`shape`](@ref) and [`shape_matrix`](@ref) (partition of unity) and
  [`grid_moment`](@ref).
- **`test/common/bsplines.jl`**: the cubic B-spline [`bspline3`](@ref) (values, symmetry,
  partition of unity) and its derivative [`dbspline3`](@ref) against ForwardDiff; the
  [`UniformBasis`](@ref) with [`eachnonzero`](@ref) and [`basis_values`](@ref) (partition of
  unity, four non-zero functions at each point).
- **`test/common/interface.jl`**: [`conservation_report`](@ref) of a gauge step, which keeps the
  momenta and the grid current and changes the grid charge and the energy.
- **`test/common/stages.jl`**: [`interpolate_stages`](@ref): the frames start and end on the
  stages, a gauge step keeps ``(x, p)`` in every frame, periodic positions move by the shorter way
  round (and cross the whole interval without a period, the control), and the guards for no
  stages and for stages of different sizes.
- **`test/common/svg.jl`**: [`svg_stage_panels`](@ref) and [`svg_stage_animation`](@ref): the
  number of markers, targets, clip paths and arrows, the position of a marker in pixels, the
  category colours, the key times of the animation, the escaped titles, clip ids that differ
  between two pictures, and the guard for states of different sizes.
- **`test/algorithm_figures.jl`**: [`algorithm_svg`](@ref) for every construction: the panels,
  the markers in every panel, an animated position per marker, no `NaN`, and the guard for an
  unknown key.
- **`test/standard.jl`** (group `slow`): the [standard experiments](@ref standard-experiments):
  the densities (normalisation, constancy along the discrete characteristics, the cut of the
  represented one), the quadrants kept by the gauge rezoning and the transport and filled by the
  sliced remap, the failures of the projection with the represented density and of the
  spline-projection redistribution at full size, the cell statistics of the long run, the loss
  of weight of the transport with the positive density, and the repeated sliced remaps. The tables
  of the chapter pages *Numerical experiments* and of the
  [comparison](@ref comparison-standard) give the numbers of these runs; the test checks bounds,
  and `scripts/standard/standard_experiments.jl` prints every number of the tables (about
  15 minutes).
- **`test/quality/aqua.jl`** and **`test/quality/explicit_imports.jl`**: `Aqua.test_all`, and
  that no module relies on an implicit import.

The statements of the two foundation chapters are checked where a construction uses them: the
factorisation lemma, the gauge step and the streaming pull-back in `test/remapping/canonicity.jl`;
the three forms of the symplecticity condition in suites 1 and 2 of `test/sliced/experiments.jl`;
the presymplectic form in `test/redistribution/canonical.jl`; one map for all markers, the
compressibility bracket and the fold in `test/projection/structure.jl` and
`test/redistribution/canonical.jl`; the closed-loop identity in
`test/redistribution/closedloop.jl`; and frozen against state-dependent parameters in
`test/projection/stages.jl`. The `@example` blocks of the foundation chapters compute their
numbers when the manual is built.

## [Findings against the sources](@id verification-findings)

The package reproduces the sources where they are right and records where they are not. The
findings below are those against the notes and scripts as the package ported them. The notes now
take their checks, numbers and figures from the package, quote the measured values, and no longer
have scripts of their own; the fit section of the projection note (the third item below) is the
one finding that is not corrected there.

**Hamiltonian remapping of weighted particles in 1D-1V** (the note
`Hamiltonian-Particle-Remapping.tex`). Every number reproduces. Two are rounded in the note: the
residual ``2.0 \times 10^{-8}`` of the eight-marker example is ``1.95 \times 10^{-8}``, and the
agreement ``7 \times 10^{-16}`` of the Schur-complement step with the full KKT step is
``5.6 \times 10^{-16}`` ([eight markers on two nodes](@ref remap-eight-example)).

**Hamiltonian particle projection** (the note `Hamiltonian-Particle-Projection.tex`).

- The lower bounds by which the first-order determinant and the determinant with
  ``(I - \tfrac12 \mathbb{J} K)^{-1}`` miss the non-separable stage are quoted as ``\geq 0.32``
  and ``\geq 0.48``; the measured minima are ``0.315`` and ``0.475``, so the bounds are rounded
  up ([non-separable stages](@ref projection-nonseparable)).
- The relative errors of the projected target weights are quoted without saying that they hold
  for the hole at ``x_h = \pi``; at ``x_h = \pi/2`` two bands differ in the second digit
  ([verification of the projection](@ref projection-verification)).
- The section on the fit describes a minimum-image distance, an adjoint recursion and a
  kinetic-energy constraint. The code fits the lifted positions by the Levenberg–Marquardt method
  with forward-mode Jacobians and has no kinetic-energy constraint ([the fit](@ref projection-fit)).

**Symplectic redistribution of weighted particles in phase space** (the note and its script
`verify_symplectic_redistribution.jl`). The script never passes in a cold process: 27 passes, 2
failures and 1 error. The bound of R1 lies below the round-off of its difference quotient, one
symplectic Euler step of R4 reaches another root, and R5 needs more than 4 periodic cells in
SimpleSplines 0.4. The 16-marker numbers of the research log come from a relaxation of eight
steps of size ``0.02`` that the script does not contain; the minimum weight ``0.0227``, the raw
momentum defect ``0.215``, the 36-marker overshoot ``-0.140`` and the two symplectic Euler
Jacobian numbers do not reproduce with either relaxation
([verification of the redistribution](@ref rd-verification)).

**Sliced redistribution** (the verification script of the derivation; the construction had no
written note, and the [sliced chapter](@ref sliced-overview) was its first write-up, from which
the note `Sliced-Symplectic-Particle-Redistribution.tex` is written).

- The script in the project does not parse: a `]` is missing at line 638. A working copy of the
  script passes its 23 checks in a cold process, and the package ports that copy
  ([verification of the sliced redistribution](@ref sliced-verification)).
- The headline numbers of the derivation mix two runs. The package gives 14 empty cells before
  and none after the long-time remap, a relative standard deviation of the cell counts of
  ``0.598 \to 0.367``, and an ``L^2`` difference of ``0.104`` ([suite 5b](@ref sliced-long-time)).
- The derivation states ``\det M = 1`` for the weight squeeze; the determinant is
  ``\prod_a \mu_a`` ([the canonical weight squeeze](@ref sliced-squeeze-structure)).
- The derivation has the wrong sign of ``\dot{\theta}`` in the extended flow of the composition
  defect; ``\dot{\theta} = -\partial H / \partial \bar{w}`` is right, and suite 4 runs the
  opposite sign as a control ([composition with a later flow step](@ref sliced-composition-defect)).
