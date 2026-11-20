# [What does not work](@id sliced-failures)

```@meta
CurrentModule = ParticleRemapping.SlicedRedistribution
```

Each section names a construction that fails, the reason it fails, and the control that shows
it. The controls marked *port* run in the package (`test/sliced/experiments.jl`,
`test/sliced/experiments_long_time.jl`); the observations marked *derivation* come from the
probes of the derivation of the sliced map and have no reproduction in the package.

## Maps that do not couple the markers

A symplectic map applied to each marker on its own keeps every hole with its volume
([uncoupled maps](@ref uncoupled-maps)). This holds for any single-marker map, however it is
built.

## [Lifts of the marginals](@id sliced-failure-marginal)

A quantile lift of all markers in one slice couples them, but in the limit of many markers it
acts on each marker's own coordinate. *Port*, two controls:

- [`checkerboard_experiment`](@ref): one slice in each lift leaves the quadrant counts at
  ``(200, 0, 200, 0)``, because both marginals of the checkerboard are already uniform. The
  conjugate variables stretch by up to ``1/(1 - \beta_{\mathrm{intra}})`` (largest ``|x|`` 6.77).
- [`rotated_frame_lift_probe`](@ref): a symplectic lift of the derivation. A convex quantile
  lift of all 400 markers (``\beta = 0.25``, mollifier ``0.08 \sqrt{m/N}``, no target weights) in
  six rotated frames ``(Q, P) = (\cos\theta \, x s + \sin\theta \, p/s, \, -\sin\theta \, x s + \cos\theta \, p/s)``
  per pass, four passes. It moves two markers into the empty quadrants: ``(198, 2, 200, 0)``,
  the counts that the derivation reports.

The slices remove this: within a slice of small velocity range the positions are spread over the
whole interval whatever their correlation with the velocity.

## [A repulsion with a confining wall](@id sliced-failure-wall)

A repulsive spread is the cotangent lift of ``\Phi(Q) = Q - \alpha \nabla U(Q)``
with a Gaussian pair repulsion and the wall ``k_w Q^4 / 4``. The Jacobian
``\mathrm{D}\Phi = I - \alpha \nabla^2 U`` has the diagonal entry
``1 - \alpha (\dots + 3 k_w Q_a^2)``, which turns negative for large ``|Q_a|``: ``\Phi`` folds, is
not a diffeomorphism, and its cotangent lift ``(\mathrm{D}\Phi)^{-\top} P`` blows up.
*Port*, [`repulsion_wall_probe`](@ref) with the derivation's parameters
(``\alpha = 0.12``, ``\sigma = 0.35``, ``k_w = 0.4`` in units of ``\sqrt{m/N}``, six rotated frames
per pass, the checkerboard of 400 markers):

| quantity | value |
|:--|:--|
| first sub-step with ``\lambda_{\min}(\mathrm{D}\Phi) \le 0`` | pass 3, frame 1, ``\lambda_{\min} = -4.78`` |
| first sub-step with a non-finite ``\mathrm{D}\Phi`` | pass 4, frame 3 |
| the same without the wall (``k_w = 0``), 12 passes | no fold |

In the derivation the run stops with a linear solve on a matrix with non-finite entries. The
sliced potential is strictly convex instead: its Hessian is bounded below by
``1 - \max(\beta_{\mathrm{mean}}, \beta_{\mathrm{intra}}) > 0`` everywhere, and the gradient map is a
global diffeomorphism ([proposition](@ref sliced-lift-structure)).

## [Reweighting at fixed momenta](@id sliced-failure-reweighting)

After the lifts the markers sit at ``(\hat{x}_a, \hat{p}_a)``. Setting the weights from the
density at that point and the velocity from ``\bar{v}_a = \hat{p}_a / (m \bar{w}_a)`` is
symplectic (a [gauge step](@ref gauge-step)), but the velocity grows where the weight is small, and the markers leave
the region they were moved into.

- *Derivation*: on the checkerboard, whose density is a step function, a kernel estimate of the
  input density vanished in the empty quadrants; the weights there went to zero and the velocities
  diverged (that probe also had an error in the units of the density ratio).
- *Port*, [`checkerboard_reweighting_control`](@ref): the sliced checkerboard reweighted by the
  [`kernel_density_ratio`](@ref) of the input (width 0.25, periodic on ``[-1, 1)``). The weights
  range over 0.171–2.35 times the input weight, and the largest speed is 5.69, against 1.04 after
  the lifts.
- *Derivation*: after a pilot pass that left no cell empty, a kernel reweighting at fixed momenta
  emptied cells again. *Port*, the long-time run with `squeeze = false`
  ([`long_time_remap_experiment`](@ref)): 3 empty cells, relative standard deviation 0.740 (0.598
  before the remap, 0.367 with the squeeze), ``L^2`` difference 0.337, largest ``|\bar{v}|`` 31.9.

This is why the [weight squeeze](@ref sliced-squeeze-structure) exists: it scales the momenta by
the same factor ``\mu_a`` as the weights, so the velocity of a slot stays where the lifts put it.

## [Stiff relaxation and narrow mollifiers](@id sliced-failure-stiff)

The cotangent part ``(\nabla^2 S)^{-1}`` amplifies the conjugate variables in the soft directions
of the Hessian. The slice-mean directions have the eigenvalue ``1 - \beta_{\mathrm{mean}}``
exactly, for every ``\varepsilon`` ([proposition](@ref sliced-lift-structure)); the directions
within a slice have ``1 - \beta_{\mathrm{intra}}`` plus the Laplacian of the mollified repulsion,
whose edge weights ``\varepsilon^2 / ((u_r - u_s)^2 + \varepsilon^2)^{3/2}`` vanish for pairs much
further apart than ``\varepsilon``.

- *Derivation*: with ``\beta = 0.88`` for both parts and ``\varepsilon`` below the marker spacing,
  the Hessian had the eigenvalue 0.12, the slice means were amplified by ``1/0.12 = 8.3``, and the
  positions moved far from their anchors. The sliced map therefore separates
  ``\beta_{\mathrm{mean}}`` from ``\beta_{\mathrm{intra}}`` and sets ``\varepsilon`` to 1.5–2
  target spacings ([parameters](@ref sliced-parameters)).
- *Port*, the long-time run (`scripts/sliced/long_time_experiment.jl`): with
  ``\beta_{\mathrm{mean}} = \beta_{\mathrm{intra}} = 0.88`` both Hessians have the smallest
  eigenvalue ``0.1200``, the largest ``|\hat{x} - x^{\mathrm{anc}}|`` is 0.921 against 0.361 with
  the defaults, and ``|\bar{x} - x^{\mathrm{anc}}|`` reaches 13.2. With ``\varepsilon`` = 0.1 target
  spacings ``|\hat{x} - x^{\mathrm{anc}}|`` reaches 0.559 and one cell stays empty. For
  ``\varepsilon`` from 0.1 to 2 spacings the smallest eigenvalues grow only from 0.2200 to 0.2226
  (``\Psi_q``, bound 0.22) and from 0.2000 to 0.2187 (``\Psi_p``, bound 0.20): a wider mollifier
  stiffens most directions within a slice, and the softest ones little.

## [Multipliers that depend on the state](@id sliced-failure-state-mu)

A squeeze whose multipliers are computed from the squeezed state is not the linear cotangent lift
of a constant ``M``. *Port*, suite 3: multipliers
``\mu \propto 0.3 + \exp(-\hat{x}^2/4 - (\hat{p}/L_p)^2/2)`` give the symplecticity defect 10.8.
The density-ratio multipliers of [`SlicedQuantileRemap`](@ref) depend on the state too; that remap
is canonical for frozen multipliers only
([parameters that depend on the state](@ref state-dependent-parameters)). A choice that follows
the density and keeps the map symplectic is an open question.

## Dropping a part of a canonical map

*Port*, the controls of suites 1–3: the lifts without their cotangent parts (defect 0.760), the
split ``\bar{v} = \bar{p}/(m w)`` without the recoil (0.507 for a map, 0.0306 for a flow step), a
map that reads ``w`` at fixed ``(x, p)`` (0.974 across leaves). Each part of the construction is
needed for the canonical form.
