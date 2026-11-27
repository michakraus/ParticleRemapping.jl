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
  the source notes and a table of the checks of every test file, and a library page for each
  module.
- A test suite of 48 test files in the groups `core` and `slow` (the group `slow` holds the two
  projection experiments and `test/standard.jl`, as the README, the home page and the verification
  page say), and scripts in `scripts/<construction>/` that reproduce the numbers and figures of
  each note.
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
  *Properties*), with the properties table, the limitations, the known issues (every entry of
  `KNOWN_ISSUES.md` that concerns the construction), the open questions and the counts and costs
  of each step in 3d3v.
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
  macOS and Windows; a pinned doctest job runs next to them. The manual is built on every change
  and published to the `gh-pages` branch by `DocumenterVitepress.deploydocs` (see *Changed*).
- An advisory `Downgrade` job tests the `[compat]` lower bounds; its failure does not fail the run.
- Dependabot opens the weekly `[compat]` bumps of the root `Project.toml` and the GitHub Actions
  updates.
- `codecov.yml` sets the project and patch checks to target `auto` with a 1 % threshold.
- `KNOWN_ISSUES.md`: K12 (the charge-conserving reweighting does not converge for markers in the
  columns of its spline cells: the Hessian is singular up to round-off).
- Tests: the remap states (`test/projection/frames.jl`) use 8 instead of 16 charge splines, which
  on 16 marker columns made the case of K12, and check that the charge step converged.
- A new method for the fit of the projection: `fit_remap_parameters(...; solver =
  :gauss_newton)` runs GeometricSolvers' `GaussNewton(; globalization = DogLeg())` on the
  regularised residual `[r; √β θ]`, with the budget `maxit` and GeometricSolvers' default
  tolerances, its stationarity test among them. The fit refuses a `β` that is negative, infinite
  or `NaN`, a non-integer or negative `maxit`, which the Levenberg–Marquardt fit accepted before
  (negative `maxit` returned negative `iterations`), and `maxit > typemax(Int32)`. The fit takes a
  start `θ0`, stops at `iterations = 0` when `J(θ0)` is not finite, and returns `iterations` as an
  `Int` and `status`, the solver's `SolverStatus`. `KickDriftRemap` and `NonSeparableRemap` gain a
  fourth type parameter `S <: Union{Nothing, Int}` and take `solver` and `seed`; a non-negative
  integer `seed` starts the fit at `1e-3 .* randn(Xoshiro(seed), n)`, and the `info` of `remap`
  holds `status`. Tests: the solver, the start and the zero budget,
  non-finite starts, `Int` iterations, trace output, the invalid arguments and the seed
  (`test/projection/parametrisation.jl`, `test/projection/remap.jl`).
- `KNOWN_ISSUES.md`: K20 (the kick–drift map does not infer its return type, so the fields `F0`
  and `F` of `fit_remap_parameters` with a `KickDriftParametrisation` infer `Any`; no test
  measures what this costs).

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
- The manual is built with DocumenterVitepress instead of Documenter's HTML writer, in the style
  of the VitePress tutorial on the Kolmogorov n-width: its colours, content width, code themes,
  local search and outline. `docs/make.jl` builds the site into one directory for each deploy
  base: `docs/build/1` for a local build, a push to `main` or a pull request, and `docs/build/1`,
  `docs/build/2`, … for the aliases of a release tag; the animations are embedded as
  `./<name>.mp4`, so that VitePress bundles them with their page.
- The Documentation workflow publishes the manual to the `gh-pages` branch:
  `docs/make.jl` ends in `DocumenterVitepress.deploydocs`, and the deploy decision comes from the
  CI environment. A push to `main` publishes `dev`, a release tag publishes `stable` and the
  version aliases, and a pull request publishes a preview under `previews/PR<n>`. A local build
  deploys nothing.
- The evidence of K2, K3 and K4 is re-measured by `scripts/redistribution/failure_cases.jl` on
  the same machine and code: K3 returns after 44 s instead of 46 s (timing noise), with the
  momentum change 0.165.
- The package depends on GeometricSolvers, which is not registered: `[sources]` takes it from
  `https://github.com/JuliaGNI/GeometricSolvers.jl` at the commit `d41d379`. The Julia floor
  rises from 1.12 to 1.13, the floor of GeometricSolvers, and the ForwardDiff floor from 1 to
  1.4, the floor that GeometricSolvers needs.
- The hole-filling experiment of the projection fits by GeometricSolvers' Gauss–Newton method with
  `DogLeg`: `hole_filling_experiment` takes `solver` (default `:gauss_newton`), `seed` for the
  start `1e-3 .* randn(Xoshiro(seed), n)` of every fit, and precomputed `targets`; each entry of
  `results` holds the `status` of its fit, and `verbose = true` prints its return code. On Linux
  the Levenberg–Marquardt fit at `x_h = π` left the saddle of the weight-dependent fit that it
  stalls at on macOS, and 9 checks of `test/projection/experiments/hole_filling.jl` failed;
  GeometricSolvers at `d41d379` returns the start point of a step that leaves a stationary point,
  with `SUCCESS` and `stationary`, so the fit from `θ = 0` ends on the saddle that way on Julia
  1.13.4, 1.14.0-beta1 and the nightly builds 1.15.0-DEV.996 and 1.15.0-DEV.1022, and with the
  targets moved by `1e-13` relative; on the nightly build 1.15.0-DEV.1003 (AMD EPYC 9V74, znver4)
  it ends at the same numbers with `LINESEARCH_FAILED`. The test and the script check the numbers
  of the saddle row on every Julia. Its return code is `SUCCESS, stationary` on a Julia without a
  prerelease tag (`isempty(VERSION.prerelease)`), and `SUCCESS, stationary` or
  `LINESEARCH_FAILED` on a Julia with one (`accepted_codes`, 1.14.0-beta1 among them); the
  stationarity ratio, below `1e-8`, is checked wherever the fit ends stationary. The tables of the
  manual page *Hole filling*, of the test and of `scripts/projection/hole_filling.jl` are restated
  from the new fit, each row with the return code of its fit. A number of a fit that does not end on the stationarity test is quoted as the
  range over the runs with the targets moved by `1e-13` relative (`Xoshiro(1:10)`) and on Julia
  1.13 and 1.14 where its digits move in a run or a run lies less than a tenth of a unit of its
  last digit from a rounding boundary; these are 14 cells of the
  two tables (at `x_h = π/2` the objective of control A, max |ln| of the weight-dependent map
  and of control B, and the fill, the objective, and the empty disk and the fill after the
  charge-conserving reweighting of control C; at `x_h = π` the empty disk, the fill and the
  empty disk after the charge-conserving reweighting of the non-separable map, and the fill,
  the objective, and the empty disk and the fill after the charge-conserving reweighting of
  control C). The test checks such a range widened on each side by its width, and prints the
  values of every fitted row; `scripts/projection/hole_filling_ranges.jl` rebuilds the 14 ranges
  from the experiment and the ten runs with the targets moved, each end rounded outward, and
  checks that each lies within its quote (14 of 14 on Julia 1.13.4 and 1.14.0-beta1). At `x_h = π/2` every fitted row changes but control B: the
  weight-dependent map gives the empty disk 1.43 and the fill 0.90 (was 1.58 and 0.94), the
  non-separable map 1.56 and 0.90 (was 1.54 and 0.94), and the weight-dependent and the
  non-separable fit end on `MAXITERS`. At `x_h = π` the saddle row, control A (which ends on the
  stationarity test as well) and control B keep their numbers (control A 1.88 and 0.47), the
  non-separable map gives the empty disk 1.51–1.52 and the fill 0.93–0.94 (was 1.50 and 0.94),
  and control C 2.04 and 0.54–0.55 (was 2.00 and 0.55). The claim that the kicks and drifts do
  not fill the hole at `x_h = π` is a property of the saddle: the manual (*Hole filling*,
  *Properties, limitations and higher dimensions*, *Construction*, the verification page) now
  says that the `±v` coupling makes the fit from `θ = 0` lead to a mirror-symmetric saddle at a
  point with `∂_x W = 0`, and that the kicks and drifts fill the hole once a seeded start breaks
  the symmetry: over the seeds 101–105 the empty disk is 1.35–1.45 and the fill 0.97–1.06 after
  the 40 iterations (row "seeded"). New rows "seeded" give control A and control B at `x_h = π`
  from the same starts, quoted as ranges over the five seeds (empty disk 1.49–1.61 and 1.51–1.98,
  fill 0.81–0.90 and 0.62–0.78); control A ends on `MAXITERS`, control B
  on the stationarity test for three seeds and on `MAXITERS` for two. The test checks the saddle
  with the return codes above and, where its fit ends stationary, the stationarity ratio below
  `1e-8`, the seeded weight-dependent fit of the seed
  101 with bounds, and control A and control B from the seed 101 within their ranges widened on
  each side by their width; the script runs the five seeds for the three variants, and runs every
  experiment again with the targets moved by `1e-13` relative. The test and the script bound the
  relative momentum change of every fit by `1e-15` (was `1e-12`), and the manual quotes it as
  below `10^{-15}`. `KNOWN_ISSUES.md`: K16 (the weight-dependent fit does not converge) and K17
  (on a prerelease Julia the saddle fit does not always end stationary; reported as
  JuliaGNI/GeometricSolvers.jl#124).
- The test of the seeded controls at `x_h = π` also compares control A and control B with the
  `θ = 0` result of the same run (the widened fill range of control B admits the 0.47 of
  `θ = 0`): the fill above that of `θ = 0` plus 0.07, the charge residual above 2 times and
  `max |ln(w′/w*)|` above 2.5 times. The empty disk is not compared (control B's 1.51–1.98 holds
  1.88). The manual states it.
- The manual page *Structural statements and target weights after a push* quotes the mass error
  of the BUG solution at the fixed rank one as below `10^{-12}` (was `1.1·10^{-13}`): the
  number is at the level of rounding, 1.1e-13 on macOS with Julia 1.13, 1.6e-13 with Julia 1.14,
  7.1e-14 to 8.0e-14 on the Linux runners of CI, and 2.6e-13 on Windows with Julia 1.13, so
  `test/projection/experiments/dlra_weights.jl` and `scripts/projection/dlra_weights.jl` check that
  bound and no longer its digits.
- Compatibility with `SimpleSplines` moves from `"0.4"` to `"0.5"` (compat-only bump; no
  behavioural change and no change to the Julia `1.12` floor).
- The README now links to the dev manual (`docs-dev` badge), lists the four modules in the
  manual's chapter order (`Projection`, `Redistribution`, `Remapping`, `SlicedRedistribution`),
  and demonstrates `NonSeparableRemap` (`Projection`) in the quick start.
- The numbers of the projection's fits with the Gauss–Newton method: the remap table of the
  target weights after a push (*Structural statements and target weights after a push*), and the
  checkerboard, the long run and the non-separable remap history with the positive density
  (*Numerical experiments* of the projection, *Comparison*). Each table gives the return code of
  its fits, and a number whose digits move with the targets moved by `1e-13` relative or between
  Julia 1.13 and 1.14 is quoted as the range over those runs. The weight in `Q_2 ∪ Q_4` of the
  checkerboard falls from 0.348 to 0.1875 (kick–drift) and from 0.393 to 0.3625
  (non-separable), and both maps leave cells empty in the long run (3–6 and 12, were 0 and 7).
  The kick–drift remap history fails at `t = 10` on Julia 1.13 and at `t = 15` on Julia 1.14
  (was `t = 15`), and the script accepts either outcome.
  `scripts/standard/standard_experiments.jl` checks these numbers and the return codes,
  `test/projection/experiments/dlra_weights.jl` and `scripts/projection/dlra_weights.jl` the
  return code of each fit of the remap table.

### Removed

- `remap_animation_frames`, `remap_stage_states`, `plot_remap_stages`, `animate_remap` and
  `animate_moser_targets` of the projection, `plot_redistribution_markers` and
  `animate_redistribution` of the spline-projection redistribution, and `animate_gauge_homotopy`
  of the gauge rezoning: `remap_states`, `remap_frames`, `plot_stages`, `animate_stages` and
  `plot_remap` with `targets` replace them. The projection adds `experiment_remap` for the remaps
  of `hole_filling_experiment`.
- Breaking: the Levenberg–Marquardt fit of the projection. `fit_remap_parameters`,
  `KickDriftRemap` and `NonSeparableRemap` default to `solver = :gauss_newton`, GeometricSolvers'
  `GaussNewton(; globalization = DogLeg())`, and refuse `solver = :levenberg_marquardt` and every
  other symbol with `ArgumentError("unknown solver :<name>; use :gauss_newton")`. Every fit
  returns a `SolverStatus` as `status`, and `info.status` of `remap` is never `nothing`, which
  closes K13 (the return type of the fit was a union over the two solvers). A call that passed no
  `solver` gets another fit and other numbers. The manual states the method, its QR step and its
  cost per iteration. The remap examples of the home page and of the projection's overview and
  algorithm pages show the return code of their fit, which `test/projection/remap.jl` checks; the
  text beside the figures of the steps of the kick–drift and the non-separable map names it.
