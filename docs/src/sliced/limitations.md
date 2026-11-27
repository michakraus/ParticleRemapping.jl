# [Properties, limitations and higher dimensions](@id sliced-limitations-page)

```@meta
CurrentModule = ParticleRemapping.SlicedRedistribution
```

This page collects what the sliced redistribution keeps, what it cannot do, its known issues and
open questions, and what it needs in ``d`` dimensions, in particular in 3d3v.

## [Properties](@id sliced-properties)

| quantity | status |
|:--|:--|
| symplectic form ``\omega`` | kept for frozen slices, orders, anchors and multipliers ([the composed map](@ref sliced-composed-map); suite 3 of the [experiments](@ref sliced-autonomous-suite), `test/sliced/map.jl`); the remap reads them from the state (K5) |
| total weight | kept to round-off: ``\sum_a w_{\mathrm{ref}} \mu_a = N w_{\mathrm{ref}}`` ([standard example](@ref sliced-overview), [suite 5b](@ref sliced-long-time)) |
| total momentum ``P`` | kept to round-off: both lifts and the squeeze with ``\sum_a \mu_a = N`` keep ``\sum_a p_a`` ([the weight squeeze](@ref sliced-squeeze-structure), [standard example](@ref sliced-overview)) |
| kinetic energy | not kept, not a constraint of the construction ([standard example](@ref sliced-overview)) |
| grid charge | not kept: positions and weights change ([standard example](@ref sliced-overview)) |
| grid current | not kept ([standard example](@ref sliced-overview)) |
| represented density (consistency) | approximate: the multipliers follow a kernel density ratio of the input; relative ``L^2`` difference 0.104 of the smoothed densities in [suite 5b](@ref sliced-long-time) |
| equidistribution and hole filling | fills holes that the slices see, also with uniform marginals: quadrants ``(200, 0, 200, 0) \to (98, 101, 112, 89)`` on the [checkerboard](@ref sliced-checkerboard); 14 empty cells to none, relative standard deviation ``0.598 \to 0.367`` in the [long-time run](@ref sliced-long-time) |
| parameters computed from the state | the slices, their orders and the anchors (piecewise constant in the state), and the multipliers ``\mu`` (a predictor–corrector loop on a kernel density ratio, [multipliers](@ref sliced-multipliers)) |
| cost of one remap | per axis two dense ``N \times N`` solves, ``O(N^3)`` each, and an assembly of ``O(N^2)``; the density ratio ``O(I d N^2)`` for ``I`` passes; the squeeze ``O(N d)`` ([cost](@ref sliced-algorithm), [the cost in 3d3v](@ref sliced-cost-3d3v)) |
| dimensions implemented in the package | 1D-1V, and ``d \times N`` markers along each axis, tested in 2D-2V ([per-axis version](@ref sliced-per-axis)) |

## [Limitations](@id sliced-limitations)

- **Canonical only for frozen parameters.** The slices, their orders and anchors and the
  multipliers come from the state. The map is canonical for each frozen choice
  ([the composed map](@ref sliced-composed-map)); as a map of the state it is not, by the
  [closed-loop identity](@ref closed-loop-identity) (K5). Multipliers that depend on the squeezed
  state break the canonical form of the squeeze by ``10.8`` in suite 3
  ([multipliers that depend on the state](@ref sliced-failure-state-mu)).
- **The squeeze trades velocity for position.** ``T_\mu`` keeps the velocity of a slot at the
  weight ``w_{\mathrm{ref}} \mu_a`` and scales the offset ``\hat{x}_a - x^{\mathrm{anc}}_a`` by
  ``1/\mu_a``. Slots with small multipliers, in the tails of the density, move their markers far
  from their anchors: up to 3.46 in [suite 5b](@ref sliced-long-time) (K6). Skipping the squeeze
  moves the problem into the velocities, ``|\bar{v}|`` up to 31.9
  ([reweighting at fixed momenta](@ref sliced-failure-reweighting)).
- **Relaxation, not target matching.** The lifts relax towards the targets with
  ``\beta_{\mathrm{mean}}, \beta_{\mathrm{intra}} < 1``. The Hessian is bounded below by
  ``1 - \max(\beta_{\mathrm{mean}}, \beta_{\mathrm{intra}})``, and the cotangent part amplifies the
  conjugate variables by up to its inverse; ``\beta \to 1`` reaches the targets and makes the
  lift stiff ([stiff relaxation and narrow mollifiers](@ref sliced-failure-stiff)).
- **A box and a slicing.** The targets are uniform in the box
  ``[-L_x, L_x) \times [-V_{\max}, V_{\max}]``, the slices have equal numbers of markers, and ``N``
  must be a multiple of ``B_v`` and ``B_x``. The weights carry the density; the marker
  distribution is pulled towards uniform in the box, whatever the density.
- **Only what the slices see.** A lift of one slice sees the marginal of its coordinate and cannot
  fill a hole with uniform marginals ([lifts of the marginals](@ref sliced-failure-marginal)). In
  ``d > 1`` the per-axis map slices each pair ``(x_\ell, v_\ell)`` and cannot fill a hole in the
  correlation of two axes ([higher dimensions](@ref sliced-higher-dimensions)).
- **No grid moments and no energy.** The weights ``w_{\mathrm{ref}} \mu_a`` and the velocities
  ``\bar{p}_a / (m \bar{w}_a)`` keep the total weight and the total momentum; nothing in the
  construction constrains the kinetic energy or the grid charge and current, which change in the
  [standard example](@ref sliced-overview).
- **Dense linear algebra.** Each lift assembles and solves a dense ``N \times N`` system (K7).

## [Known issues](@id sliced-known-issues)

The entries of `KNOWN_ISSUES.md` of the package that concern this chapter; the notes of the
particle-in-cell project list none for it.

- **K5** — the long-time sliced remap is canonical only for frozen multipliers: `remap` computes
  ``\mu`` from the state by a predictor–corrector loop.
- **K6** — the weight squeeze can move markers far from their anchors, by up to 3.46 in the
  long-time experiment of 900 markers.
- **K7** — each sliced lift solves a dense ``N \times N`` system with an LU factorisation,
  ``O(N^3)``, although the Hessian is block diagonal plus a rank-one term.

## [Open questions](@id sliced-open-questions)

- **Multipliers.** A choice of the multipliers ``\mu`` that follows the density and keeps the map
  symplectic as a map of the state, for example from a generating function
  ([multipliers that depend on the state](@ref sliced-failure-state-mu)).
- **Slices that depend on the state.** Whether the map with slices, orders and anchors read from
  the state is injective across the sets where two markers swap their order
  ([the sliced quantile lift](@ref sliced-lift-structure)).
- **Many markers in one slice.** The mean-field limit of a one-slice quantile lift, and the rate at
  which slices of finite size escape the [obstruction of uncoupled maps](@ref uncoupled-maps)
  ([coupling is necessary, not sufficient](@ref sliced-coupling-not-sufficient)).
- **Correlations between axes.** Which lifts fill a hole in the correlation of two axes with
  fewer slices than the full conditioning of [the cost in 3d3v](@ref sliced-cost-3d3v): lifts
  conditioned on a few chosen coordinates, or lifts in rotated frames (see
  [obstacles](@ref sliced-route-3d3v)).
- **Repeated remaps.** The [remaps during the run](@ref sliced-long-time) cover ``t \le 25`` with
  four remaps; the long-time behaviour with many remaps, and the effect of the far moves of K6 on
  the dynamics, are not tested.

## [Higher dimensions](@id sliced-higher-dimensions)

### Statements in d dimensions

The canonical pairs ``(x_{a\ell}, p_{a\ell})``, ``\ell = 1, \dots, d``, are independent: the form
is ``\omega = \sum_\ell \sum_a \mathrm{d}x_{a\ell} \wedge \mathrm{d}p_{a\ell}``. The proposition on
[the sliced quantile lift](@ref sliced-lift-structure) is a statement about a vector
``u \in \mathbb{R}^N`` and its conjugate ``v``; with ``(u, v) = (x_{\cdot\ell}, p_{\cdot\ell})`` it
holds for every axis, and the lift is the identity on the other pairs. The slices of such a lift
may be read from any coordinates: for frozen slices, orders and anchors, the lift of one
coordinate is a symplectic diffeomorphism of the whole ensemble that keeps every component of
``P``. The same holds for [the weight squeeze](@ref sliced-squeeze-structure) of each axis, and
[the composed map](@ref sliced-composed-map) holds in ``d`` dimensions with one common set of
multipliers. The common ``\mu`` is needed because each marker has one weight: the gauge step sets
``\bar{v}_{a\ell} = \bar{p}_{a\ell} / (m w_{\mathrm{ref}} \mu_a)`` in every axis, so a squeeze of
axis ``\ell`` with other multipliers ``\mu'`` would scale the velocities of that axis by
``\mu'_a / \mu_a``.

### Implementation in the package

For ``d \times N`` markers, [`remap`](@ref ParticleRemapping.remap) with a
[`SlicedQuantileRemap`](@ref) applies ``\Psi_q`` and ``\Psi_p`` to each pair: velocity slices by
``v_\ell`` for the positions ``x_\ell``, position slices by ``\tilde{x}_\ell`` for the momenta
``p_\ell`` (one [`SlicedLiftPlan`](@ref) per axis). It computes one set of multipliers from the
product kernel in all ``2d`` coordinates and squeezes every axis with them
([per-axis version](@ref sliced-per-axis)). `test/sliced/peraxis.jl` and suite 3 check it in
2D-2V.

Each lift reads its slices and targets from the two-dimensional marginal ``(x_\ell, v_\ell)``
only, so the lifts cannot see a hole in the joint distribution of different axes
([a hole between the axes](@ref sliced-hole-between-axes)).

### [The cost in 3d3v](@id sliced-cost-3d3v)

Let ``N`` be the number of markers, ``d`` the dimension, ``B_v``, ``B_x`` the numbers of slices
and ``n_v = N / B_v``, ``n_x = N / B_x`` their sizes, and ``I`` the number of passes of the
multipliers.

**The per-axis map, as implemented.** One remap has ``2d`` lifts with ``N`` unknowns each:

| step | count | cost |
|:--|:--|:--|
| sort into slices | ``2d`` sorts | ``O(d N \log N)`` |
| gradient ``\nabla S`` | pair sums within the slices | ``O(d N (n_v + n_x))`` |
| Hessian, dense | ``N \times N`` per lift | ``O(d N^2)`` memory and time |
| solve, dense LU | one per lift | ``O(d N^3)`` |
| density ratio | ``2 N^2`` kernel values in ``2d`` coordinates per pass | ``O(I d N^2)`` |
| squeeze | per axis | ``O(d N)`` |

In 3d3v the remap solves six dense systems of size ``N``. The Hessian is block diagonal, with one
block of size ``n_b`` per slice, plus the rank-one term
``\beta_{\mathrm{mean}} \mathbf{1}\mathbf{1}^\top / N``; a block solve with the Sherman–Morrison
formula costs ``O(\sum_b n_b^3) = O(N n^2)`` per lift with the slice size ``n``, and stores
``O(N n)`` entries (K7). With fixed slice sizes the lifts are then linear in ``N``, and the pair
sums of the density ratio dominate.

**Slices that see the joint distribution.** A triangular (Knothe–Rosenblatt) sequence of lifts
moves the ``2d`` coordinates one after the other, each conditioned on the coordinates before it.
With ``B`` slices along each conditioning coordinate, the ``k``-th lift has ``B^{k-1}`` slices,
and the last one ``B^{2d-1}``. In 1D-1V this is ``B`` slices for one lift, as in the sliced map;
in 3d3v the last lift has ``B^5``. With ``n`` markers per slice the remap needs
``N = n B^{2d-1}`` markers. The kernel of the multipliers has the same problem: for a uniform cover
of the box, the expected number of markers within one kernel width in all ``2d`` coordinates is
``N \big( \frac{h_x}{2 L_x} \cdot \frac{h_v}{2 V_{\max}} \big)^d``. The numbers for the defaults
(``h_x = 0.45``, ``h_v = 0.35``, ``L_x = \pi``, ``V_{\max} = 2.25``):

```@example sliced-limitations
using ParticleRemapping
m = SlicedQuantileRemap()
per_axis = m.hx / (2m.Lx) * m.hv / (2m.Vmax)
slicing = [(; d, B, slices = B^(2d - 1), markers_for_10_per_slice = 10 * B^(2d - 1))
           for d in (1, 3), B in (5, 15)]
kernel = [(; d, N, markers_per_kernel = N * per_axis^d) for (d, N) in ((1, 900), (3, 900),
              (3, 10^6))]
(; slicing, kernel, markers_for_5_per_kernel_3d3v = 5 / per_axis^3)
```

### [Obstacles and a route to 3d3v](@id sliced-route-3d3v)

- **The dense solves (K7).** The block-diagonal plus rank-one solve makes each lift
  ``O(N n^2)``; it changes no statement of the chapter, because the Hessian and the map are the
  same.
- **The density ratio.** Its ``O(N^2)`` pair sums dominate once the lifts are linear in ``N``.
  With the default widths, a kernel holds 5.0 markers of a uniform cover at ``N = 900`` in
  1D-1V, but ``1.6 \times 10^{-4}`` in 3d3v, and the same five markers need ``N = 2.9 \times 10^7``
  (the numbers above): the widths, or the estimator, must change with the dimension.
- **Correlations between axes.** The per-axis map sees ``d`` two-dimensional marginals. A full
  conditioning needs ``B^{2d-1}`` slices in its last lift, which is out of reach for any useful
  ``B`` in 3d3v. Between the two lie lifts conditioned on a few chosen coordinates, and lifts in
  rotated frames: an orthogonal rotation ``R`` of the positions with the same rotation of the
  momenta, ``(x, p) \mapsto (R x, R p)``, is canonical, and a per-axis lift in the rotated frame
  sees other two-dimensional marginals. Which of them fill the holes of a kinetic simulation is
  an [open question](@ref sliced-open-questions).
- **The state dependence and the far moves (K5, K6)** remain in every dimension.
