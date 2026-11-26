# [Comparison and open questions](@id comparison)

```@meta
CurrentModule = ParticleRemapping
```

The four constructions answer the same question with different choices of the canonical map
``\psi`` and of the new weights ([the constructions](@ref constructions)). This chapter sets them
side by side: what each moves, what it keeps exactly, what it needs, and where it stops. It does
not rank them; each table row states a condition or a property, and which one matters depends on
the simulation that uses the redistribution. The overview page of each chapter shows its remap
with the same pictures (the [standard pictures](@ref figures) S1, S2, A1 and T1), each on the
example of its chapter. A section on higher dimensions and the open questions of all chapters
close the chapter.

## [Results of the standard experiments](@id comparison-standard)

Every construction remaps the same markers ([standard experiments](@ref standard-experiments));
the page *Numerical experiments* of each chapter gives the figures and the details. "fails" names
the reason the construction stops; "rep." and "pos." are the represented and the positive density,
which matter only for the projection and the transport.

| | checkerboard: quadrants after (before ``200, 0, 200, 0``) | weight in ``Q_2 \cup Q_4`` (exact: 0 rep., 1/2 pos.) | long run: empty cells (before 14) | long run: relative std (before 0.598) | remaps every 100 steps: mean relative std, ``t \ge 5`` (no remap 0.587) |
|:--|:--|--:|--:|--:|--:|
| [kick–drift projection](@ref projection-experiments), rep. | fails: zero target weights | | fails | fails | |
| kick–drift projection, pos. | ``131, 73, 130, 66`` | 0.348 | 0 | 0.475 | fails at ``t = 15``: markers outside the box |
| non-separable projection, pos. | ``123, 76, 120, 81`` | 0.393 | 7 | 0.937 | fails |
| [spline redistribution](@ref rd-experiments) | fails: relaxation leaves the box (K2) | | fails (K2) | fails (K2) | |
| [gauge rezoning](@ref remap-experiments) | ``200, 0, 200, 0`` | 0 | 13 | 0.559 | 0.550 |
| canonical transport, ``\psi = \mathrm{id}``, rep. | ``200, 0, 200, 0`` | 0 | 15 | 0.572 | 0.597 |
| canonical transport, ``\psi = \mathrm{id}``, pos. | ``200, 0, 200, 0`` | 0 | 0 | 0.392 | 0.384 |
| [sliced redistribution](@ref sliced-experiments) | ``134, 67, 128, 71`` | 0.205 | 0 | 0.367 | 0.365 |

The rows follow the [obstructions](@ref obstructions). The gauge rezoning and the transport with
the identity keep every position and every sign of ``v``, so the empty quadrants stay empty. The
projection with equal weights moves every marker by one area-preserving map: it fills the
quadrants by folding the occupied region, not by shrinking the empty one, and it needs a density
that is positive on its whole box. The sliced redistribution couples the markers and fills the
quadrants with any density, but its weights follow a kernel estimate of the input, not a given
density. The spline-projection redistribution couples the markers too, but its relaxation and its
step fit 16 markers; its reduced runs are on its [experiments page](@ref rd-experiments). With
the positive density, the transport changes the total weight by ``-0.20`` on the long run, the
weight of the tails that the markers never sampled.

## Moved and kept quantities

| | [Hamiltonian particle projection](@ref projection-index) | [spline-projection redistribution](@ref redistribution-overview) | [gauge rezoning](@ref remap-gauge) | [canonical transport](@ref remap-transport) | [sliced redistribution](@ref sliced-overview) |
|:--|:--|:--|:--|:--|:--|
| dimension in the package | 1D-1V (stages in any ``d``) | 1D-1V | 1D-1V | 1D-1V | any ``d``, per axis |
| ``\psi`` | kicks and drifts, or non-separable stages, of a weight-dependent Hamiltonian | one step (or substeps) of a pairwise spline Hamiltonian | identity | a given canonical map | sliced quantile cotangent lifts and a linear weight squeeze |
| positions | moved by ``\psi`` | moved by ``\psi`` | unchanged | moved by ``\psi`` | moved by ``\psi`` |
| new weights | target weights from a projected or low-rank density, or charge-conserving | weighted KKT spline projection at the targets | stretches ``\lambda_a`` from a constrained optimisation | current-quantile matching | ``w_{\mathrm{ref}} \mu_a`` from the squeeze multipliers |
| ``\omega`` | exactly, for fixed parameters | exactly, for fixed parameters | exactly, for any state dependence | exactly, for a fixed ``\psi`` | exactly, for frozen slices and multipliers |
| total momentum ``P`` | exactly | exactly | exactly | exactly when ``\psi`` keeps it | exactly |
| total weight | exactly | exactly | exactly (moment ``k = 0``) | up to the quantile error | exactly |
| grid charge | not kept; exactly with the charge-conserving reweighting | not kept | exactly (constraint) | up to the quantile error | not kept |
| kinetic energy | not kept (the constraint of the fit is not implemented) | not kept | exactly (constraint) | up to the quantile error | not kept |
| represented density | consistent with ``f_h`` up to the fit residual | consistent with the KKT projection at the targets | in the constrained moments | consistent with the target density, up to the quantile error | through the multipliers, approximately |

## Requirements of each construction

| | needs | computed from the markers |
|:--|:--|:--|
| Hamiltonian particle projection | a box, a hat projection, a B-spline parametrisation of the stages, a regularisation | the target points (density-interpolation flow), the parameters (least-squares fit), the target weights |
| spline-projection redistribution | a box, a phase-space spline basis, relaxation parameters, a pairwise spline basis | the relaxed targets, the KKT weights, the parameters (one linear solve per step) |
| gauge rezoning | a grid, the set ``\mathcal{K}`` of constrained moments, target stretches | the stretches, by a Newton method with a homotopy |
| canonical transport | a grid, an admissible canonical map, a target density | the cumulative currents and their quantiles per cell |
| sliced redistribution | a box, the numbers of slices ``B_v`` and ``B_x`` (``N`` a multiple of both), the lift and kernel parameters | the slices and anchors, the multipliers (predictor–corrector loop) |

## Limits of each construction

Each chapter derives its own limits on its page *Properties, limitations and higher dimensions*,
with the known issues of its code and its open questions; this list collects the main limits.

- **Hamiltonian particle projection.** The fit can reduce its residual only where the old weights
  vary ([weight-dependent Hamiltonians](@ref weight-dependent-hamiltonians)); with separable
  stages, a kick that compresses the markers at one velocity expands them at the opposite one,
  which the [non-separable stages](@ref projection-nonseparable) avoid. The target weights must be
  positive, which fails for hats finer than the marker spacing. Details:
  [limitations](@ref projection-limitations).
- **Spline-projection redistribution.** The step reaches its targets on the identity-connected
  branch when the contraction condition holds; the verification case violates the sufficient
  bound and uses substeps ([the parameter solver](@ref rd-solver)). The defaults fit the
  16-marker verification set-up; the pages [the verification set-up](@ref rd-setup),
  [hole filling](@ref rd-holes) and the [standard experiments](@ref rd-experiments) show the
  failures outside it (K2, K3). Details: [limitations](@ref redistribution-limitations).
- **Gauge rezoning.** It moves the velocities along rays only and never moves a position, so it
  cannot fill a hole in ``x``. Its solution branch ends at a fold of the homotopy
  ([one constraint in closed form](@ref remap-charge-example)). The package solves the Schur
  complement as a dense matrix (K9). Details: [limitations](@ref remapping-limitations).
- **Canonical transport.** A consistent weight field exists if and only if ``\psi`` preserves the
  sign-resolved current profile, and a smooth canonical map cannot, in general, make the marker
  distribution uniform ([the limit of equidistribution](@ref remap-limit)). A folded sheet is
  treated one layer at a time ([the fold](@ref velocity-momentum-fold)). Details:
  [limitations](@ref remapping-limitations).
- **Sliced redistribution.** The map is canonical for frozen slices and multipliers, and the
  multipliers come from the state. The squeeze can move markers far from their anchors (K6), and
  each lift solves a dense ``N \times N`` system (K7, [cost](@ref sliced-cost-3d3v)). Details:
  [limitations](@ref sliced-limitations).

Every construction computes something from the markers and is symplectic for its frozen value,
except the gauge rezoning, whose gauge step is symplectic for any stretch
([parameters that depend on the state](@ref state-dependent-parameters)).

The two constructions of the gauge-rezoning chapter also combine: the current-quantile matching gives
targets for the gauge rezoning ([combining the two constructions](@ref remap-comparison)).

## [Higher dimensions](@id comparison-higher-dimensions)

The statements of every chapter hold in ``d`` dimensions; the code and the cost differ. Each
chapter derives the counts on its page *Properties, limitations and higher dimensions*, from a
script of the package that checks the counts for ``d = 1`` against the code. This table collects
what binds first in 3d3v.

| | in the package | the step that binds in 3d3v | a route the chapter derives |
|:--|:--|:--|:--|
| [Hamiltonian particle projection](@ref projection-higher-dimensions) | stages in any ``d``, the rest in 1D-1V | the grid in ``2d`` dimensions (``1.2 \cdot 10^9`` hat nodes), its Poisson problem (K10), and the parameters (``3.1 \cdot 10^6`` kick–drift, ``2.4 \cdot 10^8`` non-separable) with the dense Jacobian of their fit | transform-based Poisson solver, low-rank target weights, adjoint gradients with L-BFGS or matrix-free Gauss–Newton, a sparse or low-rank parametrisation |
| [spline-projection redistribution](@ref redistribution-higher-dimensions) | 1D-1V | the pairwise basis, ``N_\theta \approx 3 \cdot 10^9`` parameters, the ``O(N^2)`` pair sums and the dense ``2Nd \times 2Nd`` Newton matrix | the coarse collective (Fourier-moment) basis with the Tikhonov solve and the velocity recovery along the fold; a relaxation with a kernel cutoff |
| [gauge rezoning](@ref remapping-higher-dimensions) | 1D-1V | the Schur complement on ``n_K N_g^d`` multipliers, solved dense (K9): ``1.9 \cdot 10^{14}`` operations for ``N_g = 32``, ``n_K = 2`` | a sparse solve of the Schur complement (``27``-node stencil); the stretch per marker and the diagonal Hessian stay as in 1D; the stretch cannot change the direction of ``v`` |
| [canonical transport](@ref remapping-higher-dimensions) | 1D-1V | the discrete matching: the continuum matching holds along rays, a discrete one does not exist | direction bins, the radial cumulative current, canonical maps that keep the direction-resolved current |
| [sliced redistribution](@ref sliced-higher-dimensions) | any ``d``, per axis | the dense ``N \times N`` lifts (K7), the ``O(N^2)`` density ratio whose kernel holds ``1.6 \cdot 10^{-4}`` markers at ``N = 900``, and the correlations between the axes that the per-axis lifts do not see | block-diagonal plus rank-one lift solves; kernel widths that grow with ``d``; lifts conditioned on chosen coordinates or in rotated frames |

## [Open questions](@id open-questions)

### Repeated remaps and state dependence

- The long-time behaviour of a simulation with repeated remaps, given that a state-dependent
  remap breaks the symplecticity of the composite map ([two readings](@ref two-readings)). The
  composition of a push and the projection map is symplectic as a sequence of maps for fixed
  weights and parameters, but not as one map of the initial state.
- A choice of the sliced multipliers ``\mu`` that follows the density and keeps the sliced map
  symplectic as a map of the state, for example from a generating function
  ([multipliers that depend on the state](@ref sliced-failure-state-mu)), or more generally a
  [variational closure](@ref variational-closure) that matches targets.
- The injectivity of the sliced map with slices and anchors read from the state, across the sets
  where two markers swap their order ([slices that depend on the state](@ref sliced-lift-structure)).
- Whether the push conserves ``P``: the redistributions conserve it exactly, the push depends on
  the field discretisation.

### Targets and weights

- The choice of admissible generators ``h`` of the canonical transport that reduce clustering
  most, and how much such a remap gains, within the [limit of equidistribution](@ref remap-limit).
- Folded sheets, where the matching must split ``f`` between the layers
  ([the fold](@ref velocity-momentum-fold)).
- The target marker density: an equidistribution is one choice. A density ``g \propto f_h^a``
  with ``0 < a < 1`` moves the fold outwards, to ``|v| = v_{\mathrm{th}} / \sqrt{1 - a}`` for a
  Maxwellian, and reduces the spread of the weights.
- The projection of the weights: a dynamical low-rank density at ``\vartheta \leq 10^{-5}``
  gives a small relative error in the tails in the [test after a push](@ref projection-push-test),
  at rank 13 or less, in a weak, fixed field in ``d = 1``. Open are the rank that a
  self-consistent field and ``d = 3`` need, a truncation that bounds the relative error in the
  tails, and the sign of ``f_r``, which the test never saw negative.
- Charge against filling: the [charge-conserving reweighting](@ref projection-charge) keeps the
  grid charge of the markers before the map, with the charge deficit of a hole, and lowers the
  fill (from between ``0.93`` and ``0.94`` to ``0.70`` in the
  [hole-filling experiment](@ref projection-hole-filling)).
  Which charge is the right reference depends on whether the hole is a defect of the sampling or
  a feature of ``f``. The correction keeps the grid charge only.

### Solvers and stages

- The [non-separable stages](@ref projection-nonseparable) end at the iteration limit of their
  fixed-point solver, so the solvability of the midpoint equation limits the step; a Newton
  solver, more stages with smaller steps, or a composition with kicks and drifts may relax this.
  One application of the non-separable stages takes 0.20 to 0.34 s against 0.004 to 0.007 s for
  the kicks and drifts in the run of the script.
- The [coarse-basis regime](@ref rd-coarse-basis) of the spline-projection redistribution at
  large ``N``: with ``N_\theta \ll 2 N d - d`` the regularised parameter solve reaches the targets
  only approximately, and the velocities are recovered along the branch of the fold. The package
  tests the recovery on 16 markers; it has no experiment of this regime at large ``N``.

### Higher dimensions

- The current-quantile matching in several velocity dimensions: the matching holds along rays in
  the continuum ([higher dimensions](@ref remapping-higher-dimensions)); a discrete matching with
  direction bins, and the admissible canonical maps in ``d`` dimensions, are open.
- A sparse solve of the Schur complement of the gauge rezoning on ``N_g^3`` nodes, and its
  iteration count (K9).
- The tensor product of the non-separable stages has ``n_\xi^d n_v^d n_\ell`` parameters per
  stage, which is too many for ``d = 3``; a sparse or low-rank parametrisation is needed there.
- The Poisson problem of the density-interpolation flow lives on a grid in ``2d`` dimensions,
  ``1.2 \cdot 10^9`` nodes in 3d3v with the hats of the experiment, and the package solves it with
  a dense pseudo-inverse (K10) ([higher dimensions](@ref projection-higher-dimensions)).
- [`SplineRedistribution`](@ref ParticleRemapping.Redistribution.SplineRedistribution)
  implements ``d = 1`` only, while the statements of its chapter hold in ``d`` dimensions; its
  pairwise basis has ``N_\theta \approx 3 \cdot 10^9`` parameters in 3d3v with the defaults, so
  the regime of exact targets is out of reach there, and the coarse collective basis is untested
  at large ``N`` ([higher dimensions](@ref redistribution-higher-dimensions)).
- Correlations between the axes: the per-axis sliced map sees the two-dimensional marginals
  ``(x_\ell, v_\ell)`` only ([per-axis version](@ref sliced-per-axis)); lifts conditioned on a
  few chosen coordinates or lifts in rotated frames lie between the per-axis map and a full
  conditioning with ``B^{2d-1}`` slices ([higher dimensions](@ref sliced-higher-dimensions)).

### Tests

- A dynamical test: the holes of the projection experiments are cut out by hand, and the push of
  the [test after a push](@ref projection-push-test) only shears them. A test on holes that a
  simulation produces, for example in a bump-on-tail instability, is open for every
  construction.
- The mean-field limit of a one-slice quantile lift, and the rate at which slices of finite size
  escape the [obstruction of uncoupled maps](@ref uncoupled-maps)
  ([coupling is necessary, not sufficient](@ref sliced-coupling-not-sufficient)).
