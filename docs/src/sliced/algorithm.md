# [Algorithm](@id sliced-algorithm)

```@meta
CurrentModule = ParticleRemapping.SlicedRedistribution
```

The sliced redistribution acts on 1D-1V markers ``(x_a, v_a, w_a)``, ``a = 1, \dots, N``, in the
box ``[-L_x, L_x) \times [-V_{\max}, V_{\max}]``, periodic in ``x``. In ``d`` dimensions it acts
on each axis (see [below](@ref sliced-per-axis)). Let ``w_{\mathrm{ref}} = \frac{1}{N}\sum_a w_a``
and ``p_a = m w_a v_a``. The numbers of slices ``B_v`` and ``B_x`` divide ``N``;
``n_v = N / B_v`` and ``n_x = N / B_x`` are the slice sizes. The subscripts ``q`` and ``p`` of
the parameters name the lift ``\Psi_q`` or ``\Psi_p`` they belong to, as in the keywords
`β_mean_q` and `ε_p_spacings` of [`SlicedQuantileRemap`](@ref).

## [The map](@id sliced-map-algorithm)

**Algorithm (sliced redistribution).**

1. **Velocity slices.** Sort the markers by ``v`` and cut them into ``B_v`` slices of ``n_v``
   markers; order each slice by ``x``. The position targets are the ``n_v`` uniform quantiles
   ``\tau^x_r = -L_x + (2r - 1) L_x / n_v`` of ``[-L_x, L_x]``. The momentum anchor of a marker
   in slice ``b`` is ``p^{\mathrm{anc}}_a = m w_{\mathrm{ref}} (-V_{\max} + (2b - 1) V_{\max} / B_v)``.
2. **``\Psi_q``.** ``(\tilde{x}, \tilde{p}) = (\nabla S_1(x), \, p^{\mathrm{anc}} + (\nabla^2 S_1(x))^{-1} (p - p^{\mathrm{anc}}))``
   with the [sliced quantile potential](@ref sliced-lift-structure) ``S_1`` of these slices,
   ``\beta_{\mathrm{mean}, q}``, ``\beta_{\mathrm{intra}, q}`` and ``\varepsilon_q``.
3. **Position slices.** Sort the markers by ``\tilde{x}`` and cut them into ``B_x`` slices of
   ``n_x`` markers; order each slice by ``\tilde{p}``. The position anchor of a marker in slice
   ``b`` is ``x^{\mathrm{anc}}_a = -L_x + (2b - 1) L_x / B_x``. The momentum targets are
   ``\tau^p_j = m w_{\mathrm{ref}} \, v_j`` with the slot velocities
   ``v_j = -V_{\max} + (2j - 1) V_{\max} / n_x``.
4. **``\Psi_p``.** ``(\hat{x}, \hat{p}) = (x^{\mathrm{anc}} + (\nabla^2 S_2(\tilde{p}))^{-1} (\tilde{x} - x^{\mathrm{anc}}), \, \nabla S_2(\tilde{p}))``
   with the potential ``S_2`` of these slices, ``\beta_{\mathrm{mean}, p}``,
   ``\beta_{\mathrm{intra}, p}`` and ``\varepsilon_p``.
5. **Multipliers.** Compute ``\mu > 0`` with ``\sum_a \mu_a = N`` from the density of the input
   markers (the [predictor–corrector loop](@ref sliced-multipliers)).
6. **``T_\mu``.** ``(\bar{x}, \bar{p}) = T_\mu(\hat{x}, \hat{p})``, the
   [weight squeeze](@ref sliced-squeeze-structure) around the anchors ``x^{\mathrm{anc}}``.
7. **``\mathcal{F}``.** Keep the positions ``\bar{x}``, set ``\bar{w}_a = w_{\mathrm{ref}} \mu_a``
   and ``\bar{v}_a = \bar{p}_a / (m \bar{w}_a)``.

[`remap`](@ref) with a [`SlicedQuantileRemap`](@ref) runs these steps; the slices and anchors
come from [`sliced_lift_plan`](@ref), the two lifts from [`apply_sliced_lift`](@ref), the
multipliers from [`squeeze_multipliers`](@ref) and the squeeze from
[`canonical_weight_squeeze`](@ref). The mollifier widths are ``\varepsilon_q = 2 \cdot 2L_x / n_v``
and ``\varepsilon_p = 2 \cdot 2L_p / n_x``, two target spacings, with ``L_p = \max_j |\tau^p_j|``.

```@example sliced-algorithm
using ParticleRemapping
mk = electric_wave_simulation(electric_wave_markers()).markers
result = remap(SlicedQuantileRemap(), mk)
(; weight = sum(result.markers.w) - sum(mk.w),
   momentum = total_momentum(result.markers) - total_momentum(mk),
   empty_cells = (phase_cell_statistics(mk.x, mk.v).empty,
                  phase_cell_statistics(result.markers.x, result.markers.v).empty))
```

The slices, orders and anchors depend on the state; steps 2, 4 and 6 are canonical for each
frozen choice of them and of ``\mu`` ([theorem](@ref sliced-composed-map)).

## [The stages coloured by slice](@id sliced-stages-by-slice)

The figure follows the standard example of the [overview](@ref sliced-overview), 120 markers
pushed by the field of the [long-time experiment](@ref sliced-long-time) to ``T = 25``, through
the stages with ``B_v = 5`` velocity and ``B_x = 12`` position slices
([`sliced_stages_example`](@ref)). Momenta are in velocity units ``p / (m w_{\mathrm{ref}})``, the
same unit for every marker, and the markers of panels (a)–(c) are coloured by their slice, which
is what these panels show; the overview draws the same stages at the old weights, coloured and
sized by weight.

```@example sliced-algorithm
using CairoMakie
set_theme!(remap_theme())
ex = sliced_stages_example()
plot_sliced_stages(ex.before, ex.result)
```

- **(a) Input**: the markers and the boundaries of the five velocity slices (dashed), one colour
  per slice.
- **(b) ``\Psi_q``**: each velocity slice is spread over the 24 position targets (ticks at the
  bottom). The momenta change by the cotangent part ``(\nabla^2 S_1)^{-1}`` around the slice
  anchors.
- **(c) ``\Psi_p``**: the twelve position slices (between the dashed lines, one colour each,
  the colours repeat after seven) have their momenta moved towards the ten momentum targets (grey
  lines).
- **(d) ``T_\mu``**: the squeeze scales the momenta by ``\mu_a`` and the position offsets by
  ``1/\mu_a``; markers in slots with small ``\mu_a`` move far from their anchors, out of the box,
  and the panel wraps their positions into it. Colour and size show the new weight
  ``w_{\mathrm{ref}} \mu_a``.
- **(e) Multipliers**: ``\mu`` against the slot momentum. It follows the density ratio of the
  input: about one in the bulk, small in the tails.
- **(f) ``\mathcal{F}``**: the output with the weights ``w_{\mathrm{ref}} \mu_a`` (colour and size)
  and the velocities ``\bar{p}_a / (m \bar{w}_a)``.

The animation moves the markers through the same stages one after the other, by linear
interpolation between consecutive stages in this plane, coloured by velocity slice.

```@example sliced-algorithm
animate_sliced_stages("sliced-slices.mp4", ex.before, ex.result)
nothing # hide
```

```@raw html
<video src="sliced-slices.mp4" controls loop muted playsinline style="width:100%"></video>
```

## [Slices and quantile targets](@id sliced-slices)

[`sorted_slices`](@ref)`(slice_key, order_key, B)` cuts the markers into ``B`` slices by one key
and orders each slice by another. [`uniform_quantile_targets`](@ref)`(L, n)` gives the midpoints
``-L + (2r - 1) L / n`` of ``n`` equal cells of ``[-L, L]``: the targets of a slice and the slice
centres. [`slice_anchors`](@ref) assigns each marker the centre of its slice.

```@example sliced-algorithm
bins = sorted_slices([0.3, -1.2, 0.8, -0.1, 1.5, -0.7], [5.0, 1.0, 3.0, 2.0, 4.0, 6.0], 2)
(; bins, targets = uniform_quantile_targets(1.0, 3),
   anchors = slice_anchors(bins, uniform_quantile_targets(2.0, 2)))
```

## [The lifts](@id sliced-lifts)

**Algorithm (sliced quantile lift ``(u, v) \mapsto (\tilde{u}, \tilde{v})``).**

1. Form the global mean ``\bar{u}`` and the slice means ``\bar{u}_b``.
2. For the marker ``a = a_r`` of slice ``b`` compute the mollified repulsion
   ``\rho_a = \sum_s |\tau_r - \tau_s| (u_a - u_{a_s}) / \sqrt{(u_a - u_{a_s})^2 + \varepsilon^2}``
   and ``\tilde{u}_a = \beta_{\mathrm{mean}} \bar{u} + (\beta_{\mathrm{intra}} - \beta_{\mathrm{mean}}) \bar{u}_b
   + (1 - \beta_{\mathrm{intra}}) u_a + \beta_{\mathrm{intra}} \rho_a / (n_b \kappa_b)``.
3. Assemble ``H = \nabla^2 S(u)``.
4. Solve ``H (\tilde{v} - v^{\mathrm{anc}}) = v - v^{\mathrm{anc}}``.

[`split_quantile_ot_lift`](@ref) implements it, [`dual_quantile_ot_lift`](@ref) applies it with
``(u, v) = (p, x)``. The potential and its Hessian are [`sliced_quantile_potential`](@ref) and
[`sliced_quantile_hessian`](@ref).

```@example sliced-algorithm
u = [-0.9, -0.2, 0.1, 0.25, 0.7, 1.4]
v = [0.3, -0.1, 0.2, 0.0, -0.4, 0.1]
bins = [[1, 2, 3], [4, 5, 6]]
taus = [uniform_quantile_targets(1.0, 3), uniform_quantile_targets(1.0, 3)]
ũ, ṽ = split_quantile_ot_lift(u, v, bins, taus; β_mean = 0.5, β_intra = 0.9, ε = 0.05)
(; ũ, ṽ, sum_u = (sum(u), sum(ũ)), sum_v = (sum(v), sum(ṽ)))
```

Each slice of three is spread towards ``\bar{u}_b + (-\tfrac{2}{3}, 0, \tfrac{2}{3})``, relaxed
by ``\beta_{\mathrm{intra}} = 0.9``, and its mean moves halfway (``\beta_{\mathrm{mean}} = 0.5``)
towards the global mean; both sums are kept.

## [The squeeze multipliers](@id sliced-multipliers)

The multipliers come from the ratio of two kernel density estimates,

```math
r_a = \frac{\sum_c w_c \, K(\xi_a - x_c, \nu_a - v_c)}{\sum_c K(\xi_a - \xi_c, \nu_a - \nu_c)} ,
\qquad
K(\Delta x, \Delta v) = \exp\Big( -\frac{1}{2} \Big( \frac{2 s^2 (1 - \cos(\Delta x / s))}{h_x^2}
  + \frac{\Delta v^2}{h_v^2} \Big) \Big) , \quad s = \frac{L_x}{\pi} ,
```

the density of the input markers at the evaluation points ``(\xi, \nu)`` over the density of the
evaluation points themselves: the weight a marker at ``(\xi_a, \nu_a)`` carries when the
evaluation points sample the input density ([`kernel_density_ratio`](@ref)).

**Algorithm (predictor–corrector multipliers).**

1. Start from ``(\xi, \nu) = (\hat{x}, \hat{p} / (m w_{\mathrm{ref}}))``.
2. Repeat three times: set ``\mu_a \propto r_a + 0.08 \, w_{\mathrm{ref}}`` with
   ``\sum_a \mu_a = N``; squeeze ``(\hat{x}, \hat{p})`` with these ``\mu`` and set
   ``(\xi, \nu) = (\bar{x}, \bar{p} / (m w_{\mathrm{ref}} \mu))``.

The floor ``0.08 \, w_{\mathrm{ref}}`` keeps every multiplier positive. The corrector passes
evaluate the density ratio where the markers end up after the squeeze. [`squeeze_multipliers`](@ref)
implements it.

## [Per-axis version for ``d > 1``](@id sliced-per-axis)

For ``d \times N`` positions and velocities the pairs ``(x_{a\ell}, p_{a\ell})``,
``\ell = 1, \dots, d``, are canonically independent. [`remap`](@ref) applies steps 1–4 to each
axis (velocity slices by ``v_\ell``, position slices by ``\tilde{x}_\ell``), computes one set of
multipliers with the product kernel in all ``2d`` coordinates, and squeezes every axis with these
``\mu``. [`per_axis_sliced_lift`](@ref) is the lift of each axis with given slices. Each lift
sees the two-dimensional marginal ``(x_\ell, v_\ell)`` only: a hole in the correlation of two axes
is not addressed ([higher dimensions](@ref sliced-higher-dimensions)).

```@example sliced-algorithm
using Random: MersenneTwister
rng = MersenneTwister(31)
x2, v2, w2 = 2π .* rand(rng, 2, 36) .- π, 0.8 .* randn(rng, 2, 36), 0.5 .+ rand(rng, 36)
mk2 = Markers(x2, v2, w2; m = 1.1)
r2 = remap(SlicedQuantileRemap(; B_v = 3, B_x = 6), mk2)
(; weight = sum(r2.markers.w) - sum(w2),
   momentum = total_momentum(r2.markers) - total_momentum(mk2))
```

## [Parameters](@id sliced-parameters)

| parameter | [`SlicedQuantileRemap`](@ref) | autonomous map (suite 3) | checkerboard (suite 5a) |
|:--|:--|:--|:--|
| box ``L_x``, ``V_{\max}`` | ``\pi``, 2.25 | ``\pi``, 2.0 | 1, 1 |
| slices ``B_v``, ``B_x`` | 15, 90 (``N = 900``) | 4, 4 (``N = 16``) | 10, 20 (``N = 400``) |
| ``\beta_{\mathrm{mean}}``, ``\beta_{\mathrm{intra}}`` of ``\Psi_q`` | 0.75, 0.78 | 0.60, 0.72 | 0.80, 0.85 |
| ``\beta_{\mathrm{mean}}``, ``\beta_{\mathrm{intra}}`` of ``\Psi_p`` | 0.78, 0.80 | 0.60, 0.72 | 0.80, 0.85 |
| ``\varepsilon_q`` | 2 target spacings | 0.4 | 0.15 |
| ``\varepsilon_p`` | 2 target spacings | 1.5 target spacings | ``0.15 \, m w`` |
| multipliers | density ratio, ``h_x = 0.45``, ``h_v = 0.35`` | ``\propto 0.25 + e^{-v_j^2/2.5}`` | none |
| final weights | ``w_{\mathrm{ref}} \mu`` | limited kernel ratio | unchanged |

## Cost

Each lift assembles the dense ``N \times N`` Hessian, ``O(\sum_b n_b^2)`` non-trivial entries plus
the rank-one term ``\beta_{\mathrm{mean}} \mathbf{1}\mathbf{1}^\top / N``, and solves with a dense
LU factorisation, ``O(N^3)``. The Hessian is block diagonal plus this rank-one term, so the
Sherman–Morrison formula would reduce the solve to ``O(\sum_b n_b^3)``; the package does not do
this (K7 of the [known issues](@ref sliced-known-issues)). The density ratio costs ``O(N^2)`` per
pass. The [higher dimensions](@ref sliced-higher-dimensions) section counts the cost in 3d3v. The long-time run with ``N = 900`` (500 time steps
and one remap) takes 0.6 s, compilation included, in
`scripts/sliced/verify_sliced_redistribution.jl`.
