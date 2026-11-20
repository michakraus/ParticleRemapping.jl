# Changelog

## [Unreleased]

### Added

- The initial package. It consolidates the notes and scripts on the symplectic remapping,
  projection and redistribution of weighted particle-in-cell markers behind one interface:
  `remap(method, markers)` returns a `RemapResult` with the new markers and the diagnostics of the
  method, and `conservation_report` gives the changes of the conserved quantities.
- A common layer: the `Markers` type with canonical momenta and conserved quantities, the
  canonical structure (`canonical_J`, `symplecticity_defect`, `presymplectic_matrix`,
  `pullback_residual`, `gauge_step`), periodic grids, cubic B-splines and spline bases.
- Four modules, one for each construction: `Remapping` (gauge rezoning and canonical transport
  with consistent weights), `Projection` (Hamiltonian particle projection with kick–drift and
  non-separable stages, charge-conserving reweighting and dynamical low-rank target weights),
  `Redistribution` (spline-projection redistribution) and `SlicedRedistribution` (sliced quantile
  cotangent lifts and the canonical weight squeeze). Each module also holds the worked examples and
  experiments of its note.
- The Makie extension `ParticleRemappingMakieExt`, with the figures and animations of every note.
- One manual, built by `docs/make.jl`: foundations, a chapter for each construction in one
  notation, a comparison with the open questions, a verification chapter with the findings against
  the source notes, and a library page for each module.
- A test suite of 45 test files in the groups `core` and `slow`, and scripts in
  `scripts/<construction>/` that reproduce the numbers and figures of each note.
- `KNOWN_ISSUES.md` with the entries K1–K7.
- The states of a remap, common to every construction: `RemapStage`, `remap_states(method,
  before, result)` (the markers before, after each stage of the canonical map at the old weights,
  and after the gauge step), `remap_frames` for an animation (the homotopy of the gauge rezoning,
  the stage parameters of the projection scaled from zero, the pseudo-time of the spline step, the
  canonical path of the sliced map; otherwise `interpolate_stages`, linear in `(x, p, w)`), and
  methods for every remap type, with `test/{common,remapping,projection,redistribution,sliced}/`
  tests of them.
- One figure style: `remap_theme`, `plot_stages` and `animate_stages`, and the style constants of
  the Makie extension. Markers are coloured and sized by weight on one scale per figure, target
  points are open circles, and the colour maps are `viridis` (weights), `cividis` (non-negative
  fields) and `balance` (signed fields).
- The standard pictures S1 (before and after), S2 (the states), A1 (the animation) and T1 (the
  conserved quantities) of every construction on the example of its overview page, in the manual
  and as `results/<construction>/<prefix>-{remap,stages}.{pdf,mp4}` for the notes.
- Every chapter of a construction has the same pages: an overview with the sections *The
  algorithm in short*, *Example*, *The pages of this chapter* and *The code*; the theory pages;
  *Numerical experiments* (new for the remapping and the spline-projection redistribution); and
  *Properties, limitations and higher dimensions* (new; it replaces the projection's
  *Properties*), with the properties table, the limitations, the known issues, the open questions
  and the counts and costs of each step in 3d3v.
- The manual page *Figures and animations*, with the style rules and the inventory of every
  figure and animation of the manual and of the notes; the comparison has a section on what binds
  first in 3d3v.
- Scripts that print the counts of the higher-dimension sections and check them for `d = 1`
  against the code: `scripts/remapping/cost_3d3v.jl`, `scripts/projection/higher_dimensions.jl`,
  and `scripts/redistribution/failure_cases.jl`, which also measures K2 and K3.
- Tests: the current-quantile matching along rays in two velocity dimensions
  (`test/remapping/transport.jl`), the per-axis sliced lifts on a hole between the axes
  (`test/sliced/peraxis.jl`), and the controls B and C as one map of `(x, v)`
  (`test/projection/parametrisation.jl`, which closes K8).
- `KNOWN_ISSUES.md`: K9 (the dense Schur solve of the gauge rezoning) and K10 (the dense
  pseudo-inverse of the density-interpolation flow).

### Changed

- Every figure and animation follows the common style. The figure files of the
  spline-projection redistribution are named `redistribution-*`; the slice-coloured stages of the
  sliced map are `sliced-slices.{pdf,mp4}`, and `sliced-stages.{pdf,mp4}` are its S2 and A1.
- The numerical examples of the gauge rezoning and of the canonical transport move to the page
  *Numerical experiments* of their chapter; the properties table of *Combining the two
  constructions* moves to the new limitations page.
- The evidence of K2, K3 and K4 is re-measured by `scripts/redistribution/failure_cases.jl` on
  the same machine and code: K3 returns after 44 s instead of 46 s (timing noise), with the
  momentum change 0.165.

### Removed

- `remap_animation_frames`, `remap_stage_states`, `plot_remap_stages`, `animate_remap` and
  `animate_moser_targets` of the projection, `plot_redistribution_markers` and
  `animate_redistribution` of the spline-projection redistribution, and `animate_gauge_homotopy`
  of the gauge rezoning: `remap_states`, `remap_frames`, `plot_stages`, `animate_stages` and
  `plot_remap` with `targets` replace them. The projection adds `experiment_remap` for the remaps
  of `hole_filling_experiment`.
