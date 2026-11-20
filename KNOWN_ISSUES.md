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

## Redistribution

### K2 · The relaxation step of `SplineRedistribution` does not scale with the number of markers

- **location:** `src/redistribution/method.jl:62`
- **problem:** The repulsive gradient descent of `declustering_relaxation` sums over all pairs,
  and the default `relax_step = 0.006` fits the 16 markers of the verification; for more markers
  the targets leave the velocity box and `remap` fails in the KKT solve.
- **evidence:** 900 markers with uniform positions on `[0, 2π)`, `v ~ N(0, 1)` and equal weights
  (`MersenneTwister(1)`): the relaxed targets reach `|v̄| ≈ 6.1`, outside `Vmax = 3.5`, and
  `remap(SplineRedistribution(), markers)` throws `ArgumentError: matrix contains Infs or NaNs`
  from the solve in `kkt_spline_projection`.
- **kind:** defect
- **found:** 2026-11-20

### K3 · `remap` with `SplineRedistribution` does not report a step that misses its targets

- **location:** `src/redistribution/method.jl:108`
- **problem:** `remap` returns without an error when the symplectic step does not reach the
  targets; only `info.target_error` shows it.
- **evidence:** the 244 markers of `hole_filling_problem(; nL = 16, hole_a = 0.9, hole_b = 0.8)`
  with `SplineRedistribution(; Vmax = 4.0, relax_step = 0.006 * 16 / 256)`: `remap` returns after
  46 s with `info.target_error = 1.98e8`, a change of the kinetic energy of 3394 and of the total
  momentum of 0.17; the total weight is kept to round-off.
- **kind:** defect
- **found:** 2026-11-20

### K4 · The verification case of the redistribution violates the contraction condition

- **location:** `docs/src/redistribution/solver.md`
- **problem:** The contraction number of the implicit midpoint step at the stage point is 2.27 in
  the verification (2.01 with the second relaxation), above the sufficient bound 1, and one step of
  symplectic Euler reaches another root (`‖∇_x ∇_p̄ H‖ = 1.39`); the verification uses two
  substeps for symplectic Euler, and the midpoint Newton solve still reaches the targets to
  round-off.
- **evidence:** `test/redistribution/solver.jl`, and the table "symplectic Euler target error, two
  substeps" in the verification chapter of the manual.
- **kind:** found late
- **found:** 2026-11-20

## Projection

### K8 · The header of the parametrisation test promises an area check that the file does not run

- **location:** `test/projection/parametrisation.jl:2`
- **problem:** The header says that the controls B and C preserve the area in `(x, v)`, but no
  test in the file computes an area or a Jacobian determinant of these controls; only the slow
  test `test/projection/experiments/hole_filling.jl` checks the consequence, that they do not
  shrink the hole beyond the reweighting alone.
- **evidence:** `grep -n -i "area\|det" test/projection/parametrisation.jl` matches only the
  header line.
- **kind:** missing test
- **found:** 2026-11-20

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
