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
- Standard experiments that every construction runs with the same markers and measures
  (`src/standard.jl`): the checkerboard of suite 5a and the long run in a two-wave field of suite
  5b, each with the density the markers represent and a positive density
  (`checkerboard_density`, `electric_wave_density`, the latter exact by the inverse of the
  drift–kick–drift steps); `standard_method`, `standard_checkerboard`, `standard_long_time`,
  `standard_remap_history`, `standard_wave_markers`, `STANDARD_METHODS` and `ShiftedRemap`, which
  applies a construction on `[0, L)` to markers on `[−L/2, L/2)`. A construction that fails on a
  marker set reports the message instead of throwing.
- `DensityTarget`: a target of the current-quantile matching for any density function, tabulated
  per cell and branch.
- `checkerboard_markers` (the markers of suite 5a), and `electric_wave_simulation` takes a
  `remap_method` that is a function of the remap time.
- SVG pictures that need no plotting package: `svg_stage_panels` (one panel per state),
  `svg_stage_animation` (SMIL, no script), the overlays `svg_line`, `svg_rect`, `svg_text`,
  `svg_arrow`, `svg_category_colors` and `SVGFigure`; `algorithm_svg(key)` draws the steps of each
  algorithm on a small example for the overview pages.
- Tests: `test/standard.jl` (group `slow`), `test/common/svg.jl`, `test/algorithm_figures.jl`,
  `DensityTarget` against `SeparableMaxwellianTarget`, and the parametrisations with equal weights.
- `KNOWN_ISSUES.md`: K11 (the target points of markers outside the velocity box of the projection
  diverge, which stops the repeated remaps of the long run).
- Continuous integration, the same generated workflows as the other repositories of the tree. The
  tests run the groups `core` and `slow` on Julia `min` and `1` on Ubuntu, and on fewer events on
  macOS and Windows; a pinned doctest job runs next to them. The manual is built on every change;
  nothing is published, because there is no `deploydocs`.
- An advisory `Downgrade` job tests the `[compat]` lower bounds; its failure does not fail the run.
- Dependabot opens the weekly `[compat]` bumps of the root `Project.toml` and the GitHub Actions
  updates.
- `codecov.yml` sets the project and patch checks to target `auto` with a 1 % threshold.
- `KNOWN_ISSUES.md`: K12 (the charge-conserving reweighting does not converge for markers in the
  columns of its spline cells: the Hessian is singular up to round-off).
- Tests: the remap states (`test/projection/frames.jl`) use 8 instead of 16 charge splines, which
  on 16 marker columns made the case of K12, and check that the charge step converged.
- A second method for the fit of the projection: `fit_remap_parameters(...; solver =
  :gauss_newton)` runs GeometricSolvers' `GaussNewton(; globalization = DogLeg())` on the
  regularised residual `[r; √β θ]`, with the budget `maxit` and GeometricSolvers' default
  tolerances, its stationarity test among them. The Levenberg–Marquardt fit stays the default with
  the same iteration, so no quoted number moves; both methods refuse a `β` that is negative,
  infinite or `NaN`, and a non-integer `maxit`, which the Levenberg–Marquardt fit accepted
  before. The fit takes a start `θ0` and returns
  `status`: the solver's `SolverStatus` for `:gauss_newton`, `nothing` for `:levenberg_marquardt`
  (K13). `KickDriftRemap` and `NonSeparableRemap`
  take `solver` and `seed`; an integer `seed` starts the fit at `1e-3 .* randn(Xoshiro(seed), n)`,
  and the `info` of `remap` holds `status`. Tests: both solvers, the start and the zero budget,
  the invalid arguments and the seed (`test/projection/parametrisation.jl`,
  `test/projection/remap.jl`).

### Fixed

- The kick–drift and the non-separable parametrisation of the projection evaluate to `NaN` for
  markers of equal weight: the spline basis in `ln w` had the spacing zero. It now spans
  `[ln w, ln w + 1]`, and the stages are the same map for every marker, as they must be.
- `charge_reweight`, the Newton iteration of the charge-conserving reweighting, no longer stalls
  short of its tolerance. Near the solution the decrease of Φ in a Newton step (about 1e-25 against
  Φ ≈ 6) is far below the rounding error of Φ (one ulp, 8.9e-16), so whether the Armijo test
  accepted the full step depended on how the sum rounded. Where Φ rounded up, the line search
  halved the step to 1e-12 and the iteration stopped after 50 iterations at a residual near 1e-10
  (`converged = false`, and `remap` of a `ChargeConservingReweighting` threw "no positive weights
  reproduce the reference grid charge"). This happened on Linux with Julia 1.13 in the example of
  the manual page on the charge-conserving reweighting and with Julia 1.14 in
  `test/projection/charge.jl`. The line search now evaluates the change of Φ along the step,
  `Σ w (exp(s Aᵀδ) − 1) − s δᵀρ` with `expm1`, instead of the difference of two values of Φ, so
  its rounding error stays below the decrease of the step. Each of 20 starts within 1e-11 of the
  grid charge converges in one step (`test/projection/charge.jl`); before the fix one of them
  stopped after 50 iterations and one took 29.

### Changed

- Every figure and animation follows the common style. The figure files of the
  spline-projection redistribution are named `redistribution-*`; the slice-coloured stages of the
  sliced map are `sliced-slices.{pdf,mp4}`, and `sliced-stages.{pdf,mp4}` are its S2 and A1.
- The numerical examples of the gauge rezoning and of the canonical transport move to the page
  *Numerical experiments* of their chapter; the properties table of *Combining the two
  constructions* moves to the new limitations page.
- The manual has one structure per chapter: an overview with the algorithm step by step in SVG
  figures and an animation, the pages of the construction, an *Algorithm* page, *Numerical
  experiments* with the standard experiments, the additional experiments (with a page *Hole
  filling* where the construction has a hole experiment), and the limitations. The chapters run
  projection, spline-projection redistribution, gauge rezoning and canonical transport, sliced
  redistribution; *Figures and animations* follows *Library*. The page *Standard experiments*
  closes the foundations, and the comparison gives the results of all constructions side by side.
  Headings are no questions: *Symplectic redistribution*, *Constructions that fail*, and others.
  The numerical examples of the gauge rezoning and the canonical transport are now on the page
  *Clustered markers and a separable density*, and *Combining the two constructions* is a section
  of the *Algorithm* page.
- The evidence of K2, K3 and K4 is re-measured by `scripts/redistribution/failure_cases.jl` on
  the same machine and code: K3 returns after 44 s instead of 46 s (timing noise), with the
  momentum change 0.165.
- The package depends on GeometricSolvers, which is not registered: `[sources]` takes it from
  `https://github.com/JuliaGNI/GeometricSolvers.jl` at the commit `47d817b`. The Julia floor
  rises from 1.12 to 1.13, the floor of GeometricSolvers, and the ForwardDiff floor from 1 to
  1.4, the floor that GeometricSolvers needs.

### Removed

- `remap_animation_frames`, `remap_stage_states`, `plot_remap_stages`, `animate_remap` and
  `animate_moser_targets` of the projection, `plot_redistribution_markers` and
  `animate_redistribution` of the spline-projection redistribution, and `animate_gauge_homotopy`
  of the gauge rezoning: `remap_states`, `remap_frames`, `plot_stages`, `animate_stages` and
  `plot_remap` with `targets` replace them. The projection adds `experiment_remap` for the remaps
  of `hole_filling_experiment`.
