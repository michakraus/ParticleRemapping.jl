# [Gauge rezoning and canonical transport](@id remap-overview)

```@meta
CurrentModule = ParticleRemapping.Remapping
```

By the [factorisation lemma](@ref factorisation-lemma), every symplectic redistribution is a
canonical map ``\psi`` of ``(x, p)``, symplectic for ``\omega = \sum_a \mathrm{d} x_a \wedge
\mathrm{d} p_a``, followed by a [gauge step](@ref gauge-step). This chapter works out the two
extreme cases of this factorisation in 1D-1V, on a periodic grid with ``N_g`` nodes and linear
shape functions ``S_i``. Both are exactly symplectic for ``\omega_N``.

- The **gauge rezoning** takes ``\psi = \mathrm{id}``. It moves no marker in ``(x, p)``: it
  stretches each velocity along its ray, ``V_a = \lambda_a v_a``, and divides the weight by the
  same ``\lambda_a > 0``. A constrained optimisation chooses the stretches. The positions, the
  momenta, the total momentum ``P`` and the grid current are exact by construction, and the grid
  charge and the kinetic energy are exact by constraints. The gauge step preserves ``\omega_N``
  also for stretches that depend on the whole ensemble, so the remap itself needs no
  closed-loop argument ([parameters that depend on the state](@ref state-dependent-parameters));
  only the composition with the flow is subject to the [two readings](@ref two-readings).
- The **canonical transport** takes a given canonical map ``\psi``, the same for all markers, and
  chooses the new weights by *current-quantile matching*, a gauge step at the mapped ``(X, P)``,
  so that the markers represent the same density ``f``. In the continuum limit these weights are
  unique, and they exist if and only if ``\psi`` preserves the sign-resolved current profile. One
  map for all markers [moves holes but keeps their area](@ref holes-persist); the
  [limit of equidistribution](@ref remap-limit) states this obstruction for consistent weights.

Throughout the chapter, lower-case letters ``x_a, v_a, w_a, p_a = m w_a v_a`` denote a marker
before the remap and capitals ``X_a, V_a, W_a, P_a`` after it; a star marks a target, ``\mu^*_a``
or ``V^*_a``. The log-stretch is ``\mu_a = \ln \lambda_a``, the constrained grid moments are
``\mathcal{M}_i^{(k)} = \sum_a w_a S_i(x_a) v_a^k`` for ``k \in \mathcal{K}``, and the
cumulative current of ``f`` at ``X`` is ``F_X(V) = \int_0^V f(X, V') V' \, \mathrm{d} V'``.

## The algorithm step by step

**Gauge rezoning** ([`GaugeRezoning`](@ref)). The figures follow the steps on the eight markers
on two nodes of the [worked example](@ref remap-eight-example)
([`algorithm_svg`](@ref ParticleRemapping.algorithm_svg)); the vertical grey lines are the rays
along which the markers move, the black lines the cell boundaries, the open circles the target
speeds.

1. **Targets.** In each cell and for each sign of ``v``, set targets ``\mu^*_a`` that spread the
   speeds evenly (panel 2).
2. **Homotopy.** Minimise the weighted distance ``\sum_a w_a (\mu_a - \mu^*_a)^2`` subject to
   unchanged grid charge and grid kinetic energy, by a Newton method on the KKT system, and follow
   the minimiser by a homotopy from ``\mu^* = 0`` to the full targets (panels 3 and 4). The
   homotopy stops early where the solution branch has a fold.
3. **Gauge step.** Set ``W_a = w_a e^{-\mu_a}`` and ``V_a = p_a / (m W_a)`` (panel 5). A marker
   that speeds up loses weight.

```@example remap-steps
using ParticleRemapping
gauge_steps = algorithm_svg(:gauge)
write_algorithm_player(".", gauge_steps) # hide
nothing # hide
```

```@raw html
<script setup>
import gauge from './algorithm-gauge.json'
import transport from './algorithm-transport.json'
</script>
```

The player runs through the steps automatically; pause it or step back and forward through them:

```@raw html
<AlgorithmPlayer :player="gauge" />
```

```@example remap-steps
gauge_steps.panels
```

**Canonical transport** ([`CurrentQuantileMatching`](@ref)). The figures follow the steps on 48
markers in two cells with the admissible cotangent lift of the
[closed-form example](@ref remap-lift-example).

1. **Canonical map.** ``(X_a, P_a) = \psi(x_a, p_a)``, at the old weights (panel 2, arrows from
   the markers before).
2. **Groups.** In each cell and for each sign of ``P``, sort the markers by ``|P_a|`` and form
   the cumulative momenta (panel 3, the colour is the group).
3. **Gauge step.** Give each marker the velocity ``V_a`` at the same quantile of the cumulative
   current of ``f``, and the weight ``W_a = P_a / (m V_a)`` (panel 4, each marker moves along
   ``v \propto 1/w`` at fixed ``x``).

```@example remap-steps
transport_steps = algorithm_svg(:transport)
write_algorithm_player(".", transport_steps) # hide
nothing # hide
```

The player runs through the steps automatically; pause it or step back and forward through them:

```@raw html
<AlgorithmPlayer :player="transport" />
```

```@example remap-steps
transport_steps.panels
```

Steps 2 and 3 of the transport are a gauge step. The velocities they give can serve as targets
for the gauge rezoning, which then makes the grid charge and the kinetic energy exact
([combining the two constructions](@ref remap-comparison)). The [algorithm](@ref remap-algorithm)
page gives both in full, each with a worked example.

## Example

### Gauge rezoning

The standard example of the gauge rezoning is [`gauge_rezoning_setup`](@ref): ``N = 6400``
markers on ``N_g = 16`` nodes of ``[0, 2\pi)``, whose velocities cluster around six random centres
in each cell, with weights from a two-stream density, and the constraints
``\mathcal{K} = \{0, 2\}``. The figure shows the markers before and after the remap, coloured and
sized by weight. The positions and the momenta are the same in both panels; the targets are not
drawn, since an open circle for each of the 6400 markers would hide them.

```@example overview
using ParticleRemapping, CairoMakie
set_theme!(remap_theme())
gauge = gauge_rezoning_setup()
mk, grid = gauge.markers, gauge.grid
r = remap(gauge.method, mk)
plot_remap(mk, r.markers;
    titles = ("(a) before", "(b) after, τ = $(round(r.info.τ; digits = 3))"))
```

The canonical map is the identity, so the states of the remap are points of the homotopy: the
markers before, the gauge steps of two accepted homotopy steps, and the result at the parameter
``\tau`` where the homotopy stops ([`remap_states`](@ref ParticleRemapping.remap_states)). Every
state has the same ``(x, p)``; only the velocities and the weights change.

```@example overview
plot_stages(remap_states(gauge.method, mk, r); columns = 2, panel_size = (400, 280))
```

The animation follows the homotopy through its accepted steps, with the log-stretches
interpolated between them ([`remap_frames`](@ref ParticleRemapping.remap_frames)):

```@example overview
animate_stages("remapping-gauge-stages.mp4", remap_frames(gauge.method, mk, r))
nothing # hide
```

```@raw html
<video src="./remapping-gauge-stages.mp4" controls loop muted playsinline style="width:100%"></video>
```

The changes of the conserved quantities
([`conservation_report`](@ref ParticleRemapping.conservation_report)); the grid moments are
relative to the moments of ``|v|``:

```@example overview
conservation_report(mk, r.markers; grid)
```

The total weight, the total momentum, the kinetic energy and the grid charge, current and energy
are kept to round-off. The [numerical experiments](@ref remap-gauge-example) give this example in
more detail.

### Canonical transport

The standard example of the canonical transport is [`quantile_matching_setup`](@ref) with
``N = 10^4`` markers of equal weight, sampled from ``f = (1 + 0.4 \cos x) e^{-v^2/2}/\sqrt{2\pi}``
on ``16`` cells, and the admissible cotangent lift of ``a = J^{-1/2}`` at ``t = 1.5`` as
``\psi`` ([the closed-form example](@ref remap-lift-example)). The method has no target points.

```@example overview
transport = quantile_matching_setup(; N = 10_000)
tm, method = transport.markers, transport.admissible
rt = remap(method, tm)
plot_remap(tm, rt.markers; titles = ("(a) before", "(b) after"), period = (0, 2π))
```

The states: the markers before; after the cotangent lift at the old weights, ``V = P/(m w)``,
which moves the markers in ``x`` and changes their velocities by ``1/\varphi_t'``; and after the
current-quantile matching, the gauge step that gives each marker its matched velocity and the
weight ``W = P/(m V)``:

```@example overview
plot_stages(remap_states(method, tm, rt); period = (0, 2π), panel_size = (300, 260))
```

```@example overview
animate_stages("remapping-transport-stages.mp4", remap_frames(method, tm, rt);
    period = (0, 2π))
nothing # hide
```

```@raw html
<video src="./remapping-transport-stages.mp4" controls loop muted playsinline style="width:100%"></video>
```

The changes of the conserved quantities, the largest compatibility residual ``|R_c^s|`` of the
matching, and the range of the weight ratio ``W_a / w_a`` with the speed of the markers at its
ends:

```@example overview
ratio = rt.markers.w ./ tm.w
(conservation_report(tm, rt.markers; grid = transport.grid)...,
    max_residual = rt.info.max_residual, weight_ratio = extrema(ratio),
    speed_at_extremes = abs.(rt.markers.v[[argmin(ratio), argmax(ratio)]]))
```

The matching keeps every momentum ``P_a`` of the map; the cotangent lift does not keep the total
momentum, since ``\sum_a P_a = \sum_a p_a / \varphi_t'(x_a)``. The total weight, the kinetic energy
and the grid moments are reproduced to the accuracy of the quantile rule, with about ``300``
markers per cell and branch here. The largest changes of weight are those of the slowest marker
of a branch, where ``P`` and ``V`` are both small. The
[numerical experiments](@ref remap-transport-example) give the example with ``2 \times 10^5``
markers.

## The pages of this chapter

- [Velocity rezoning by the gauge step](@ref remap-gauge): the map and its exact properties, the
  constrained optimisation, its optimality conditions, the Newton step, and a worked example in
  closed form with the fold of the solution branch.
- [Canonical transport with consistent weights](@ref remap-transport): the matching proposition,
  the compatibility condition, the limit of equidistribution, admissible canonical maps, and a
  worked example in closed form.
- [Algorithm](@ref remap-algorithm): the gauge rezoning and the current-quantile matching step by
  step, with worked examples, and their combination.
- [Numerical experiments](@ref remap-experiments): the standard experiments, which every
  construction runs.
- [Clustered markers and a separable density](@ref remap-examples): the two standard examples in
  detail, with the homotopy up to its fold and the admissible and inadmissible lifts.
- [Properties, limitations and higher dimensions](@ref remap-limitations-page): what each
  construction keeps, what it cannot do, the open questions, and the cost in 3d3v.
- [Library: Remapping](library.md): the docstrings of the module [`Remapping`](@ref).

## The code

The module `ParticleRemapping.Remapping` implements both constructions in 1D-1V on a
[`PeriodicGrid`](@ref ParticleRemapping.PeriodicGrid). The method types [`GaugeRezoning`](@ref) and
[`CurrentQuantileMatching`](@ref) implement [`remap`](@ref ParticleRemapping.remap),
[`remap_states`](@ref ParticleRemapping.remap_states) and
[`remap_frames`](@ref ParticleRemapping.remap_frames). The steps of the gauge rezoning are public:
the targets ([`gauge_targets`](@ref)), the constraints and the optimality conditions
([`gauge_constraints`](@ref), [`gauge_stationarity`](@ref), [`gauge_hessian`](@ref)), the Newton
step ([`gauge_newton_step`](@ref)) and the homotopy ([`gauge_homotopy`](@ref)). So are those of the
transport: the cotangent lift ([`cotangent_lift`](@ref)) and the matching of one branch and of a
grid ([`branch_quantile_matching`](@ref), [`current_quantile_match`](@ref)).

The tests in `test/remapping/` check every claim and every number of the chapter, and the scripts
in `scripts/remapping/` print the numbers of its tables; `scripts/remapping/make_figures.jl`
writes the figures and the animations to `results/remapping/`. The
[verification](@ref remap-verification) lists what is checked where.
