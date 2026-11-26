# Known issues

## Upstream

### K1 · SimpleSplines 0.4 returns NaN derivatives at `x = 0` of a periodic basis

- **location:** `src/redistribution/hamiltonians.jl:42`
- **problem:** A periodic B-spline basis of SimpleSplines 0.4, evaluated at a dual number exactly
  at the knot `x = 0`, returns NaN derivative partials for every basis function; the value is
  correct, and the derivative at another knot is finite. The pairwise Hamiltonian of
  `Redistribution` evaluates the self pair `x_a − x_a` at the plain number 0, whose derivative
  vanishes identically, to avoid it; the workaround can go when SimpleSplines is fixed.
- **evidence:**
  ```julia
  using SimpleSplines, ForwardDiff
  using SimpleSplines: BSplineBasis, UniformMesh, Periodic, evaluate
  bx = BSplineBasis(UniformMesh(8, (0.0, 2π)), 4, Periodic())
  ForwardDiff.derivative(t -> evaluate(bx, 6, t, 0), 0.0)      # NaN
  ForwardDiff.derivative(t -> evaluate(bx, 6, t, 0), 1e-12)    # -0.6366197723683921
  ForwardDiff.derivative(t -> evaluate(bx, 6, t, 0), 2π / 8)   # -0.2122065907891938
  ```
- **kind:** upstream
- **found:** 2026-11-20

### K15 · The Downgrade job cannot resolve the unregistered GeometricSolvers

- **location:** `Project.toml:25`
- **problem:** the Downgrade job (`julia-actions/julia-downgrade-compat@v2`, mode `forcedeps`,
  `julia_version: min`, projects `.` and `test`) resolves the lower bounds of `[compat]` from the
  General registry, which does not hold GeometricSolvers, so it stops before any test runs; the
  job is advisory, and the entry can go when GeometricSolvers is registered in General.
- **evidence:** run 45206146166 of PR 5, job "Downgrade - ubuntu-latest", Julia 1.13.4:
  ```
  Info: Found source package: GeometricSolvers (url=https://github.com/JuliaGNI/GeometricSolvers.jl)
  Info: Adding runtime dependency from direct path source ParticleRemapping: GeometricSolvers
  Info: Running resolver on merged project with --min=@deps
  ERROR: LoadError: unknown package UUID: 3a2b8526-dde5-4630-ad43-3d2dcd69f118
  ```
  The UUID is GeometricSolvers' (`Project.toml:10`).
- **kind:** upstream
- **found:** PR 5

### K17 · On Julia nightly the saddle fit of the hole-filling experiment does not end stationary

- **location:** `test/projection/experiments/hole_filling.jl`, `code_asserted`
- **problem:** the weight-dependent fit of the hole-filling experiment from `θ = 0` at `x_h = π`
  ends on the mirror-symmetric saddle with `SUCCESS` and `stationary` on Julia 1.13.4 and
  1.14.0-beta1; on the nightly build 1.15.0-DEV.1003 on an AMD EPYC 9V74 (znver4) it ends after 2
  iterations with `LINESEARCH_FAILED`, at the quoted numbers but with the stationarity ratio
  `4.03e-5` against the bound `1e-8`. The test and `scripts/projection/hole_filling.jl` therefore
  assert the return code and the stationarity of this row only on a Julia without a prerelease
  tag; on a prerelease Julia, 1.14.0-beta1 among them, they assert its numbers. Reported as
  [JuliaGNI/GeometricSolvers.jl#124](https://github.com/JuliaGNI/GeometricSolvers.jl/issues/124);
  the entry can go when GeometricSolvers keeps the saddle stationary there.
- **evidence:** CI run 45653025090, job "Julia nightly - ubuntu-latest", on commit `cad5bb8`:
  ```
  hole table, x_h = 3.1416, weight_dependent: 1.807122024142029 0.4678362573099415 0.9265457135502585 0.004241059149177599 1.844606460379962 0.4678362573099415 0.007980562889722034, LINESEARCH_FAILED
  Expression: stationarity_ratio(e, :weight_dependent) < 1.0e-8
   Evaluated: 4.0260004156085914e-5 < 1.0e-8
  ```
  On the same run, control A and control B at `x_h = π` end `SUCCESS, stationary`. In CI run
  45596795331, on commit `a00c653` with the same GeometricSolvers commit `d41d379`, the nightly
  build 1.15.0-DEV.996 on an AMD EPYC 7763 (znver3) ended the row `SUCCESS, stationary`.
- **kind:** upstream
- **found:** 2026-11-25

## Remapping

### K9 · The Schur-complement step of the gauge rezoning solves a dense system

- **location:** `src/remapping/gauge.jl:229`
- **problem:** `gauge_newton_step` with `solver = :schur` converts the sparse, periodic
  block-tridiagonal Schur complement `G diag(H)⁻¹ Gᵀ` to a dense `Matrix` and solves it by LU, at
  the cost `O((|K| N_g)³)`; a block-tridiagonal solve costs `O(|K|³ N_g)`. In 3d3v with
  `N_g = 32` per dimension and `K = {charge, energy}` the dense solve takes `1.9e14` operations and
  `3.4e10` bytes.
- **evidence:** `Matrix(gauge_schur_matrix(H, G)) \ (c .- G * (gJ ./ H))`;
  `scripts/remapping/cost_3d3v.jl` prints the counts.
- **kind:** found late
- **found:** 2026-11-20

## Redistribution

### K2 · The relaxation step of `SplineRedistribution` does not scale with the number of markers

- **location:** `src/redistribution/method.jl:62`
- **problem:** The repulsive gradient descent of `declustering_relaxation` sums over all pairs,
  and the default `relax_step = 0.006` fits the 16 markers of the verification; for more markers
  the targets leave the velocity box and `remap` fails in the KKT solve.
- **evidence:** 900 markers with uniform positions on `[0, 2π)`, `v ~ N(0, 1)` and equal weights
  (`MersenneTwister(1)`): 309 of the 900 relaxed targets leave the velocity box, up to
  `|v̄| = 6.16`, outside `Vmax = 3.5`, and `remap(SplineRedistribution(), markers)` throws
  `ArgumentError: matrix contains Infs or NaNs` from the solve in `kkt_spline_projection`;
  `scripts/redistribution/failure_cases.jl` prints it.
- **kind:** defect
- **found:** 2026-11-20

### K3 · `remap` with `SplineRedistribution` does not report a step that misses its targets

- **location:** `src/redistribution/method.jl:108`
- **problem:** `remap` returns without an error when the symplectic step does not reach the
  targets; only `info.target_error` shows it.
- **evidence:** the 244 markers of `hole_filling_problem(; nL = 16, hole_a = 0.9, hole_b = 0.8)`
  with `SplineRedistribution(; Vmax = 4.0, relax_step = 0.006 * 16 / 256)`: `remap` returns after
  44 s with `info.target_error = 1.98e8`, the contraction number 4.85e8 and the Newton residual
  3.96e8, a change of the kinetic energy of 3394 (of 2.933) and of the total momentum of 0.165;
  the total weight is kept to 1.8e-15 (`scripts/redistribution/failure_cases.jl`).
- **kind:** defect
- **found:** 2026-11-20

### K4 · The verification case of the redistribution violates the contraction condition

- **location:** `docs/src/redistribution/solver.md`
- **problem:** The contraction number of the implicit midpoint step at the stage point is 2.266 in
  the verification (2.01 with the second relaxation), above the sufficient bound 1, and one step of
  symplectic Euler reaches another root (`‖∇_x ∇_p̄ H‖ = 1.391`); the verification uses two
  substeps for symplectic Euler, and the midpoint Newton solve still reaches the targets to
  round-off.
- **evidence:** `test/redistribution/solver.jl`, and the table "symplectic Euler target error, two
  substeps" in the verification chapter of the manual.
- **kind:** found late
- **found:** 2026-11-20

## Projection

### K10 · The density-interpolation flow solves its Poisson problem with a dense pseudo-inverse

- **location:** `src/projection/targets.jl:105`
- **problem:** `moser_targets` forms the dense five-point Laplacian on the `nX × nV` nodes of the
  hats and applies its pseudo-inverse, `O(M³)` operations on `M` nodes; this fits `d = 1`
  (`M = 1056`) but not `d = 3` (`M = 1.2·10⁹`); a Fourier transform in `x` and a cosine transform
  in `v` diagonalise the operator at `O(M log M)`.
- **evidence:** `scripts/projection/higher_dimensions.jl`; the section *The cost in 3d3v* of the
  manual's projection chapter.
- **kind:** found late
- **found:** 2026-11-20

### K11 · The target points of markers outside the velocity box of the projection diverge

- **location:** `src/projection/targets.jl`
- **problem:** for a marker with `|v| > V`, the density-interpolation flow of `projection_targets`
  is not defined, and its target moves far out of the box; in a run with repeated remaps the
  speeds that the gauge step gives grow from remap to remap until the density underflows at a
  target and `remap` stops with a zero target weight.
- **evidence:** `standard_remap_history(:projection; variant = :positive)` (the long run of the
  standard experiments, box `V = 3`): after the first remap at `t = 5` the largest speed is 6.48;
  the largest target speed is 2.94 at `t = 5`, 12.3 at `t = 10` and `7.5·10²⁹` at `t = 15`, where
  `remap` throws `DomainError with 0.0: the target weights must be positive`. The same holds for
  `:nonseparable`.
- **kind:** defect
- **found:** 2026-11-22

### K12 · The charge-conserving reweighting fails for markers in the columns of its spline cells

- **location:** `src/projection/charge.jl`
- **problem:** with one marker column at the centre of each of the `n` spline cells, the cubic
  B-splines at the markers do not see the alternating mode, so the Hessian `A diag(w) Aᵀ` of
  `charge_reweight` is singular up to round-off; whether its Cholesky test passes depends on the
  platform, and the iteration stops with `converged = false`, after no step or after one; inside
  `KickDriftRemap` and `NonSeparableRemap` the weights `wc` of such a run replace the target
  weights with only `charge_converged = false` to show it, while `remap` of a
  `ChargeConservingReweighting` throws.
- **evidence:** `hole_filling_problem(; nL = 16, hole_a = 0.9, hole_b = 0.8)` with
  `KickDriftRemap(; nξ = 8, nv = 10, nw = 4, projection = HatProjection(12, 13), maxit = 2,
  stages = 3, charge = ChargeConservingReweighting(2π, 16))`: smallest eigenvalue of the Hessian
  `4.4·10⁻¹⁷`, largest `0.49`; one Newton step on macOS with Julia 1.13, none with 1.14 and on
  Linux. With 8 splines the eigenvalues are `0.034` and `0.91`, and the iteration converges.
- **kind:** defect
- **found:** PR 3

### K13 · The return type of `fit_remap_parameters` is a union over the two solvers

- **location:** `src/projection/fit.jl`
- **problem:** the field `status` is a `SolverStatus` for `:gauss_newton` and `nothing` for
  `:levenberg_marquardt`, so the inferred return type of `fit_remap_parameters`, and the `info` of
  `remap`, is a `Union` of two `NamedTuple` types for a `solver` that is not a constant; no test
  measures what this costs.
- **evidence:** `Base.return_types` of a call with a `Symbol` argument `solver` and a
  `KickDriftParametrisation`, Julia 1.13:
  `Union{NamedTuple{(:θ, :F0, :F, :iterations, :status), <:Tuple{Vector{Float64}, Any, Any, Int64,
  GeometricSolvers.SolverStatus{Float64}}}, NamedTuple{(:θ, :F0, :F, :iterations, :status),
  <:Tuple{Any, Any, Any, Int64, Nothing}}}`. The `Any` fields come from `remap_stages` of the
  kick–drift map, which infers `Tuple{Any, Any}`; with a `NonSeparableParametrisation` both
  members are concrete.
- **kind:** not verified
- **found:** 2026-11-23

### K16 · The weight-dependent fit of the hole-filling experiment does not converge

- **location:** `src/projection/experiments/hole_filling.jl`
- **problem:** with `β = 1e-3`, GeometricSolvers' `GaussNewton(; globalization = DogLeg())` does
  not reach a minimum of the weight-dependent fit at either hole position: after hundreds of
  iterations the gradient and `|θ|` stay large, so every weight-dependent row of the hole-filling
  tables, apart from the saddle at `x_h = π` that the fit from `θ = 0` reaches, is the value after
  the budget of 40 iterations, not at a minimum; a larger `β` or a continuation in `β` may give a
  minimum. The experiment keeps `β = 1e-3` and the budget of 40 iterations.
- **evidence:** the fit of the experiment's weight-dependent variant run past the budget on macOS
  with `BLAS.set_num_threads(2)`, by Levenberg–Marquardt and by GeometricSolvers' `GaussNewton`
  with `Backtracking` or `DogLeg`, from `θ = 0`, from targets moved by `1e-13` relative and from
  `θ0 = 1e-3 .* randn(Xoshiro(100 + k), n)`: after 100–300 iterations `‖g‖/‖g0‖` is 0.045–5.9
  and `|θ|` is 19–48 at `x_h = π/2` and at `x_h = π`. For example, `DogLeg` from `θ0` with
  `k = 1` at `x_h = π`: `MAXITERS` after 300 iterations, `F/F0 = 0.297382`,
  `‖g‖/‖g0‖ = 1.779`, `|θ| = 43.72`.
- **kind:** found late
- **found:** 2026-11-23

## SlicedRedistribution

### K5 · The long-time sliced remap is canonical only for frozen multipliers

- **location:** `src/sliced/map.jl:164`
- **problem:** `remap(::SlicedQuantileRemap, ·)` computes the squeeze multipliers `μ` from the
  state by a predictor–corrector loop, so the map is canonical for the frozen plan and `μ` only.
- **evidence:** `squeeze_multipliers` takes the lifted state; the manual's section *Parameters
  that depend on the state* states the consequence.
- **kind:** docs
- **found:** 2026-11-20

### K6 · The weight squeeze can move markers far from their anchors

- **location:** `src/sliced/squeeze.jl`
- **problem:** In the long-time experiment of 900 markers the squeeze moves markers by up to 3.46
  from their anchors.
- **evidence:** `test/sliced/experiments_long_time.jl`; the variants table of the sliced
  experiments in the manual.
- **kind:** found late
- **found:** 2026-11-20

### K7 · Each sliced lift solves a dense `N × N` system

- **location:** `src/sliced/lift.jl`
- **problem:** Each lift assembles the dense Hessian and solves with an LU factorisation,
  `O(N³)`; the Hessian is block diagonal plus a rank-one term, so a block solve with the
  Sherman–Morrison formula would cost `O(Σ n_b³)`.
- **evidence:** the section *Cost* of the sliced algorithm in the manual.
- **kind:** found late
- **found:** 2026-11-20
