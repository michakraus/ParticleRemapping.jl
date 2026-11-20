# [Verification: gauge rezoning and canonical transport](@id remap-verification)

Every claim and every number of the chapter
[Gauge rezoning and canonical transport](@ref remap-overview) is checked by a test file in
`test/remapping/`, and each check runs with a control that must fail; the control asserts the
margin by which it fails. The scripts in `scripts/remapping/` print the numbers of the tables and
check them again. Run a test file with

```
julia --startup-file=no --project=test test/remapping/<file>.jl
```

and a script with

```
julia --startup-file=no --project=scripts scripts/remapping/<file>.jl
```

## What is checked where

- **`test/remapping/canonicity.jl`**, script `symplectic_condition.jl`: the
  [factorisation lemma](@ref factorisation-lemma) (sufficiency), the [gauge step](@ref gauge-step)
  with a state-dependent stretch, the [streaming pull-back](@ref two-readings), and the
  canonicity of ``X = \varphi(x)``, ``V = c(x) v + b(x)``, ``W = w / (\varphi' c)`` for a frozen
  ``w``
  ([`symplectic_condition_checks`](@ref ParticleRemapping.Remapping.symplectic_condition_checks)).
  Pull-back residuals are below ``3 \times 10^{-9}`` (central differences); the controls
  ``X = x + v`` and a non-affine ``V`` fail with residuals ``3.4`` and ``1.2``, which the tests
  compare with their closed forms ``m \max|v|`` and ``m w |q(v) - 1|``. 7 of 7 as expected.
- **`test/remapping/gauge.jl`**: the exact properties of the rezoning for random stretches, the
  Jacobian, the stationarity and the diagonal Hessian against automatic differentiation, the
  Schur-complement step against the full KKT step and its block-tridiagonal structure, and the
  rezoning for every combination of solver and constraint set.
- **`test/remapping/gauge_examples.jl`**, script `worked_examples.jl`: every number of the
  [charge-constraint example](@ref remap-charge-example) and of the
  [eight-marker example](@ref remap-eight-example): the solution by the Lambert function and the
  fold, the Newton iterations and the homotopy, with the Schur-complement step checked against the
  full KKT step. Controls: the targets violate the constraints, Newton without homotopy is
  rejected.
- **`test/remapping/transport.jl`**: the identity in the proof of the
  [matching proposition](@ref remap-transport), to a relative error of ``1.8 \times 10^{-14}``
  (the frozen-weight control fails by ``5 \times 10^{-2}``); the cotangent lift of
  ``a = J^{-1/2}`` keeps the current profile, to ``6.4 \times 10^{-15}`` (the control
  ``a = 1 + 0.5 \cos x`` fails by ``5.7``); the matching on a grid; and the identity of the
  matching along rays in two velocity dimensions
  ([higher dimensions](@ref remapping-higher-dimensions)), to below ``10^{-12}``, with the
  determinant of ``V \mapsto m W V`` from automatic differentiation and against its closed form
  ``m^d W^{d-1} \partial_r (W r)`` (the frozen-weight control ``m^d W^d`` fails by more than
  ``10^{-3}``).
- **`test/remapping/transport_examples.jl`**, script `worked_examples.jl`: every number of the
  [closed-form lift](@ref remap-lift-example) and of the
  [five-marker matching](@ref remap-five-example). Controls: the lift of ``a = 1 + 0.5\cos x`` is
  incompatible, the marker density at ``t = 0`` is uneven, and the old velocities miss the
  current quantiles.
- **`test/remapping/experiments.jl`**, scripts `gauge_rezoning.jl` and
  `current_quantile_matching.jl`: the [numerical example of the gauge rezoning](@ref remap-gauge-example)
  and the [numerical example of the canonical transport](@ref remap-transport-example), with their
  tables. Controls: the unconstrained targets violate the charge constraint,
  ``\mathcal{K} = \{0, 2\}`` changes the third moment, the admissible lift without matching
  misrepresents ``f``, and the inadmissible lift has a large compatibility residual.
- **`test/remapping/frames.jl`**: the states and the animation frames of both constructions
  ([`remap_states`](@ref ParticleRemapping.remap_states),
  [`remap_frames`](@ref ParticleRemapping.remap_frames)). The first state is the input and the
  last is the result; every state of the gauge rezoning keeps ``(x, p)``, and its inner states,
  accepted points of the homotopy, meet the charge and energy constraints to ``10^{-13}``; the
  canonical map of the transport keeps the weights, and its matching keeps ``(X, P)``. Controls:
  the inner states differ from the input, the interpolated frames between two accepted points
  miss the constraints by more than ``10^{-6}``, and at the old weights the moved markers
  misrepresent ``f`` by more than three times the error after the matching.
- **`scripts/remapping/cost_3d3v.jl`**: the counts of one Newton step of the gauge rezoning,
  checked in 1D against the matrices of the package (non-zeros of ``\mathsf{G}``, size and
  stencil of the Schur complement, with the periodic coupling as a control), and evaluated in
  ``d`` dimensions for the [cost in 3d3v](@ref remapping-higher-dimensions). 7 of 7 as expected.
- **`scripts/remapping/make_figures.jl`**: writes the figures and animations listed in the last
  subsection; it checks no number.

## Counts against the source note

The worked examples of the source note give 21 checks with their controls, the three
verifications of the symplecticity condition, the gauge rezoning and the matching 7, 10 and 9; the
tests port all 47, with the numbers of the tables at the digits they are quoted with.

The necessity in the [factorisation lemma](@ref factorisation-lemma), the second-order proposition
of the [gauge rezoning](@ref remap-gauge), the [limit of equidistribution](@ref remap-limit) and
the proposition on admissible generators in general rest on the proofs given in the chapter; the
tests check special cases only.

## Figures and animations of the remapping

`scripts/remapping/make_figures.jl` writes the figures (PDF) and the animations (MP4) into
`results/remapping/`:

```
julia --startup-file=no --project=scripts scripts/remapping/make_figures.jl
```

| file                                   | content                                                                 | function                                     |
|:---------------------------------------|:------------------------------------------------------------------------|:---------------------------------------------|
| `remapping-gauge-remap.pdf`            | S1 of the [gauge rezoning](@ref remap-overview)                          | [`plot_remap`](@ref ParticleRemapping.plot_remap) |
| `remapping-gauge-stages.pdf`           | S2 of the gauge rezoning: before, two homotopy steps, after              | [`plot_stages`](@ref ParticleRemapping.plot_stages) of [`remap_states`](@ref ParticleRemapping.remap_states) |
| `remapping-gauge-stages.mp4`           | A1 of the gauge rezoning, along the homotopy                             | [`animate_stages`](@ref ParticleRemapping.animate_stages) of [`remap_frames`](@ref ParticleRemapping.remap_frames) |
| `remapping-transport-remap.pdf`        | S1 of the [canonical transport](@ref remap-overview)                     | [`plot_remap`](@ref ParticleRemapping.plot_remap) |
| `remapping-transport-stages.pdf`       | S2 of the canonical transport: before, ``\psi`` at the old weights, after the matching | [`plot_stages`](@ref ParticleRemapping.plot_stages) of [`remap_states`](@ref ParticleRemapping.remap_states) |
| `remapping-transport-stages.mp4`       | A1 of the canonical transport                                            | [`animate_stages`](@ref ParticleRemapping.animate_stages) of [`remap_frames`](@ref ParticleRemapping.remap_frames) |
| `remapping-ex1.pdf`                    | [charge-constraint example](@ref remap-charge-example)                  | [`plot_charge_example`](@ref ParticleRemapping.Remapping.plot_charge_example) |
| `remapping-ex1-fold.mp4`               | the same, target scale ``\beta`` from ``0`` to the fold                 | [`animate_charge_fold`](@ref ParticleRemapping.Remapping.animate_charge_fold) |
| `remapping-ex2.pdf`                    | [eight-marker example](@ref remap-eight-example)                        | [`plot_eight_marker_example`](@ref ParticleRemapping.Remapping.plot_eight_marker_example) |
| `remapping-ex2-homotopy.mp4`           | the same, homotopy from ``\tau = 0`` to ``1``                           | [`animate_eight_marker_homotopy`](@ref ParticleRemapping.Remapping.animate_eight_marker_homotopy) |
| `remapping-kkt.pdf`                    | the KKT matrix and its Schur complement                                 | [`plot_kkt_structure`](@ref ParticleRemapping.Remapping.plot_kkt_structure) |
| `remapping-gauge.pdf`                  | [numerical example of the gauge rezoning](@ref remap-gauge-example)    | [`plot_gauge_experiment`](@ref ParticleRemapping.Remapping.plot_gauge_experiment) |
| `remapping-gauge-homotopy.mp4`         | the same, cells 1 and 2, homotopy up to the fold                        | [`animate_stages`](@ref ParticleRemapping.animate_stages) of [`remap_frames`](@ref ParticleRemapping.remap_frames) |
| `remapping-ex3.pdf`                    | [closed-form lift](@ref remap-lift-example) and the numerical example   | [`plot_lift_example`](@ref ParticleRemapping.Remapping.plot_lift_example) |
| `remapping-ex3-lift.mp4`               | the lift on ``2 \times 10^5`` markers, ``t`` from ``0`` to ``1.5``       | [`animate_lift_matching`](@ref ParticleRemapping.Remapping.animate_lift_matching) |
| `remapping-ex4.pdf`                    | [five-marker matching](@ref remap-five-example)                         | [`plot_matching_example`](@ref ParticleRemapping.Remapping.plot_matching_example) |
| `remapping-ex4-matching.mp4`           | the same, one marker at a time                                          | [`animate_matching_example`](@ref ParticleRemapping.Remapping.animate_matching_example) |
| `remapping-gauge-orbits.pdf`           | [the gauge direction](@ref kernel-gauge)                                | [`plot_gauge_orbits`](@ref ParticleRemapping.Remapping.plot_gauge_orbits) |
| `remapping-sheet.pdf`                  | [the velocity–momentum fold](@ref velocity-momentum-fold)               | [`plot_velocity_momentum_sheet`](@ref ParticleRemapping.Remapping.plot_velocity_momentum_sheet) |

The figures and animations are build products. The script regenerates them, and the pages of the
manual build them again.
