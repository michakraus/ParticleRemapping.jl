# [Comparison and open questions](@id comparison)

```@meta
CurrentModule = ParticleRemapping
```

The four constructions answer the same question with different choices of the canonical map
``\psi`` and of the new weights ([the constructions](@ref constructions)). This chapter sets them
side by side: what each moves, what it keeps exactly, what it needs, and where it stops. It does
not rank them; each table row states a condition or a property, and which one matters depends on
the simulation that uses the redistribution. The open questions of all chapters close the
chapter.

## What moves, and what is kept

| | [gauge rezoning](@ref remap-gauge) | [canonical transport](@ref remap-transport) | [Hamiltonian particle projection](@ref projection-index) | [spline-projection redistribution](@ref redistribution-overview) | [sliced redistribution](@ref sliced-overview) |
|:--|:--|:--|:--|:--|:--|
| dimension in the package | 1D-1V | 1D-1V | 1D-1V (stages in any ``d``) | 1D-1V | any ``d``, per axis |
| ``\psi`` | identity | a given canonical map | kicks and drifts, or non-separable stages, of a weight-dependent Hamiltonian | one step (or substeps) of a pairwise spline Hamiltonian | sliced quantile cotangent lifts and a linear weight squeeze |
| positions | unchanged | moved by ``\psi`` | moved by ``\psi`` | moved by ``\psi`` | moved by ``\psi`` |
| new weights | stretches ``\lambda_a`` from a constrained optimisation | current-quantile matching | target weights from a projected or low-rank density, or charge-conserving | weighted KKT spline projection at the targets | ``w_{\mathrm{ref}} \mu_a`` from the squeeze multipliers |
| ``\omega`` | exactly, for any state dependence | exactly, for a fixed ``\psi`` | exactly, for fixed parameters | exactly, for fixed parameters | exactly, for frozen slices and multipliers |
| total momentum ``P`` | exactly | exactly when ``\psi`` keeps it | exactly | exactly | exactly |
| total weight | exactly (moment ``k = 0``) | up to the quantile error | exactly | exactly | exactly |
| grid charge | exactly (constraint) | up to the quantile error | not kept; exactly with the charge-conserving reweighting | not kept | not kept |
| kinetic energy | exactly (constraint) | up to the quantile error | not kept (the constraint of the fit is not implemented) | not kept | not kept |
| represented density | in the constrained moments | consistent with the target density, up to the quantile error | consistent with ``f_h`` up to the fit residual | consistent with the KKT projection at the targets | through the multipliers, approximately |

## What each construction needs

| | needs | computed from the markers |
|:--|:--|:--|
| gauge rezoning | a grid, the set ``\mathcal{K}`` of constrained moments, target stretches | the stretches, by a Newton method with a homotopy |
| canonical transport | a grid, an admissible canonical map, a target density | the cumulative currents and their quantiles per cell |
| Hamiltonian particle projection | a box, a hat projection, a B-spline parametrisation of the stages, a regularisation | the target points (density-interpolation flow), the parameters (least-squares fit), the target weights |
| spline-projection redistribution | a box, a phase-space spline basis, relaxation parameters, a pairwise spline basis | the relaxed targets, the KKT weights, the parameters (one linear solve per step) |
| sliced redistribution | a box, the numbers of slices ``B_v`` and ``B_x`` (``N`` a multiple of both), the lift and kernel parameters | the slices and anchors, the multipliers (predictor–corrector loop) |

## Where each construction stops

Each chapter derives its own limits; this list collects them.

- **Gauge rezoning.** It moves the velocities along rays only and never moves a position, so it
  cannot fill a hole in ``x``. Its solution branch ends at a fold of the homotopy
  ([one constraint in closed form](@ref remap-charge-example)).
- **Canonical transport.** A consistent weight field exists if and only if ``\psi`` preserves the
  sign-resolved current profile, and a smooth canonical map cannot, in general, make the marker
  distribution uniform ([the limit of equidistribution](@ref remap-limit)). A folded sheet is
  treated one layer at a time ([the fold](@ref velocity-momentum-fold)).
- **Hamiltonian particle projection.** The fit can reduce its residual only where the old weights
  vary ([weight-dependent Hamiltonians](@ref weight-dependent-hamiltonians)); with separable
  stages, a kick that compresses the markers at one velocity expands them at the opposite one,
  which the [non-separable stages](@ref projection-nonseparable) avoid. The target weights must be
  positive, which fails for hats finer than the marker spacing.
- **Spline-projection redistribution.** The step reaches its targets on the identity-connected
  branch when the contraction condition holds; the verification case violates the sufficient
  bound and uses substeps ([the parameter solver](@ref rd-solver)). The defaults fit the
  16-marker verification set-up; `KNOWN_ISSUES.md` records the failures outside it.
- **Sliced redistribution.** The map is canonical for frozen slices and multipliers, and the
  multipliers come from the state. The squeeze can move markers far from their anchors, and each
  lift solves a dense ``N \times N`` system ([cost](@ref sliced-algorithm)).

Every construction computes something from the markers and is symplectic for its frozen value,
except the gauge rezoning, whose gauge step is symplectic for any stretch
([parameters that depend on the state](@ref state-dependent-parameters)).

The two constructions of the first chapter also combine: the current-quantile matching gives
targets for the gauge rezoning ([combining the two constructions](@ref remap-comparison)).

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
  fill (from ``0.94`` to ``0.70`` in the [hole-filling experiment](@ref projection-hole-filling)).
  Which charge is the right reference depends on whether the hole is a defect of the sampling or
  a feature of ``f``. The correction keeps the grid charge only.

### Solvers and stages

- The [non-separable stages](@ref projection-nonseparable) end at the iteration limit of their
  fixed-point solver, so the solvability of the midpoint equation limits the step; a Newton
  solver, more stages with smaller steps, or a composition with kicks and drifts may relax this.
  A non-separable stage costs about 50 times as much as a kick and a drift in the experiment.
- The [coarse-basis regime](@ref rd-coarse-basis) of the spline-projection redistribution at
  large ``N``: with ``N_\theta \ll 2 N d - d`` the regularised parameter solve reaches the targets
  only approximately, and the velocities are recovered along the branch of the fold. The package
  tests the recovery on 16 markers; it has no experiment of this regime at large ``N``.

### Higher dimensions

- The extension of the current-quantile matching to several velocity dimensions, where the gauge
  step acts only on the speed.
- The tensor product of the non-separable stages has ``n_\xi^d n_v^d n_\ell`` parameters per
  stage, which is too many for ``d = 3``; a sparse or low-rank parametrisation is needed there.
- [`SplineRedistribution`](@ref ParticleRemapping.Redistribution.SplineRedistribution)
  implements ``d = 1`` only, while the statements of its chapter hold in ``d`` dimensions.
- Correlations between the axes: the per-axis sliced map sees the two-dimensional marginals
  ``(x_\ell, v_\ell)`` only ([per-axis version](@ref sliced-per-axis)).

### Tests

- A dynamical test: the holes of the projection experiments are cut out by hand, and the push of
  the [test after a push](@ref projection-push-test) only shears them. A test on holes that a
  simulation produces, for example in a bump-on-tail instability, is open for every
  construction.
- The mean-field limit of a one-slice quantile lift, and the rate at which slices of finite size
  escape the [obstruction of uncoupled maps](@ref uncoupled-maps)
  ([coupling is necessary, not sufficient](@ref sliced-coupling-not-sufficient)).
