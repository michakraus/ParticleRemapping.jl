# [Sliced redistribution](@id sliced-overview)

```@meta
CurrentModule = ParticleRemapping.SlicedRedistribution
```

The sliced redistribution moves clustered markers back towards an even cover of a box
``[-L_x, L_x) \times [-V_{\max}, V_{\max}]`` without a solve for parameters. Its canonical map
``\psi`` of ``(x, p)`` is an explicit composition of two sliced quantile cotangent lifts and a
linear canonical weight squeeze; the [gauge step](@ref gauge-step) that follows sets the weights
``\bar{w}_a = w_{\mathrm{ref}} \mu_a``, with ``w_{\mathrm{ref}} = \frac{1}{N} \sum_a w_a`` and the
multipliers ``\mu_a`` of the squeeze. It answers the [obstruction of uncoupled maps](@ref uncoupled-maps)
by slices: each lift couples the markers of a slice, and within a slice of small velocity range it
spreads the positions over the whole interval, whatever their correlation with the velocity. A
lift of the marginals alone cannot fill a hole whose marginals are uniform; the slices can
([hole filling on a checkerboard](@ref sliced-checkerboard)). For frozen slices, anchors and
multipliers the map satisfies ``\mathcal{R}^* \omega_N = \omega_N``; it keeps the total momentum
``P`` and the total weight to round-off. The multipliers of the remap are computed from the
markers, so the map is symplectic for their frozen values
([parameters that depend on the state](@ref state-dependent-parameters)).

The setting is that of the foundations: the canonical coordinates ``(x, p)`` with
``p_a = m w_a v_a``, the [presymplectic form](@ref presymplectic-form) ``\omega_N``, the
[factorisation lemma](@ref factorisation-lemma) and the
[symplecticity condition](@ref symplecticity-conditions). This chapter marks the image of each
stage by its own decoration: ``(\tilde{x}, \tilde{p})`` after the first lift,
``(\hat{x}, \hat{p})`` after the second, ``(\bar{x}, \bar{p})`` after the squeeze and
``(\bar{x}, \bar{v}, \bar{w})`` for the output markers. ``B_v`` and ``B_x`` are the numbers of
velocity and position slices, ``n_v = N / B_v`` and ``n_x = N / B_x`` their sizes. In the code
the fields `q̃`, `p̃`, `q̂`, `p̂`, `q̄` and `p̄` of the result hold these values, and the field
`q_anc` holds the position anchors ``x^{\mathrm{anc}}``.

## The algorithm in short

The map is a composition of four stages,

```math
\mathcal{R} = \mathcal{F} \circ T_\mu \circ \Psi_p \circ \Psi_q , \qquad
\psi = T_\mu \circ \Psi_p \circ \Psi_q ,
```

1. ``\Psi_q`` sorts the markers into ``B_v`` slices of equal velocity range and moves the
   positions in each slice towards equispaced targets: the cotangent lift of the gradient of a
   strictly convex potential;
2. ``\Psi_p`` does the same with the roles exchanged: it sorts the markers into ``B_x`` slices of
   equal position range and moves the momenta in each slice towards equispaced targets;
3. ``T_\mu`` scales each momentum by a slot multiplier ``\mu_a`` and each position offset by
   ``1/\mu_a``, a linear canonical map that keeps the total momentum; the multipliers follow the
   density of the input markers;
4. ``\mathcal{F}`` sets the weights ``\bar{w}_a = w_{\mathrm{ref}} \mu_a`` and the velocities
   ``\bar{v}_a = \bar{p}_a / (m \bar{w}_a)``.

Stages 1–3 are canonical maps of ``(x, p)``; stage 4 is a [gauge step](@ref gauge-step), which
any choice of weights keeps symplectic. The order of the two lifts is that of a Knothe–Rosenblatt
rearrangement [Villani:2009](@cite): first the positions conditional on the velocity, then the
momenta conditional on the position. The [algorithm](@ref sliced-algorithm) page gives each stage
in full.

## Example

The standard example is [`sliced_stages_example`](@ref): ``12 \times 10`` markers on a lattice of
``[-\pi, \pi) \times [-1.4, 1.4]`` with weights ``\propto f_0(x, v)`` of the
[long-time run](@ref sliced-long-time), pushed by 500 steps of the two-wave field to ``T = 25``
and remapped by [`SlicedQuantileRemap`](@ref) with ``B_v = 5`` velocity and ``B_x = 12`` position
slices (the other parameters at their defaults). The figure shows the markers before and after the
remap, coloured and sized by weight; the positions are wrapped into ``[-\pi, \pi)``:

```@example sliced-overview
using ParticleRemapping, CairoMakie
set_theme!(remap_theme())
ex = sliced_stages_example()
method = SlicedQuantileRemap(; B_v = 5, B_x = 12)
before, result = ex.before, ex.result
plot_remap(before, result.markers; period = (-π, π),
    titles = ("(a) before", "(b) after the sliced remap"))
```

The stages of the remap ([`remap_states`](@ref ParticleRemapping.remap_states)): the markers
before, after ``\Psi_q``, after ``\Psi_p`` and after the squeeze ``T_\mu``, at the old weights
with ``v = p / (m w)``, and after the gauge step with the weights ``w_{\mathrm{ref}} \mu``. The
lifts move the momenta towards anchors and targets of the form ``m w_{\mathrm{ref}} v``, with the
mean weight, so at the old weights the velocities spread; the squeeze scales the momenta by ``\mu_a``, and the gauge step
divides them by the new weights ``w_{\mathrm{ref}} \mu_a``:

```@example sliced-overview
plot_stages(remap_states(method, before, result); columns = 3, period = (-π, π))
```

The animation follows a path of canonical maps through the stages
([`remap_frames`](@ref ParticleRemapping.remap_frames)): the lifts with their relaxation
parameters ``\beta_{\mathrm{mean}}``, ``\beta_{\mathrm{intra}}`` scaled from 0 (the identity) to
their values, the squeeze with the multipliers ``(1 - s) + s \mu``, and then the gauge step at
fixed ``(\bar{x}, \bar{p})`` with the weights moving from ``w`` to ``\bar{w}``:

```@example sliced-overview
animate_stages("sliced-stages.mp4", remap_frames(method, before, result); period = (-π, π))
nothing # hide
```

```@raw html
<video src="sliced-stages.mp4" controls loop muted playsinline style="width:100%"></video>
```

The changes of the conserved quantities
([`conservation_report`](@ref ParticleRemapping.conservation_report)), with the grid moments on
16 cells of ``[-\pi, \pi)``, and the kinetic energy before the remap:

```@example sliced-overview
(conservation_report(before, result.markers; grid = PeriodicGrid(2π, 16))...,
    kinetic_energy_before = kinetic_energy(before))
```

The total weight and the total momentum are kept to round-off. The kinetic energy, the grid
charge and the grid current are not constraints of the construction and change. The
[algorithm](@ref sliced-algorithm) page draws the same stages coloured by slice, and the
[numerical experiments](@ref sliced-experiments) test the remap on a checkerboard and in a
long-time run of 900 markers.

## The pages of this chapter

- [The stages and their structure](@ref sliced-structure): why coupling alone does not fill a
  hole, and the statements about the sliced lift, the weight squeeze, the composed map and the
  composition with a later flow step, with proofs.
- [Algorithm](@ref sliced-algorithm): the four stages with their formulas, the stages of the
  standard example coloured by slice, the parameters and the cost.
- [Numerical experiments](@ref sliced-experiments): the five verification suites, the hole
  filling on a checkerboard and the long-time run of 900 markers in an electric field, with
  tables, figures and animations.
- [What does not work](@ref sliced-failures): the constructions that fail, with the reason and a
  control that shows it.
- [Properties, limitations and higher dimensions](@ref sliced-limitations-page): what the map
  keeps, what it cannot do, the known issues, the open questions, and the cost in 3d3v.
- [Library](@ref sliced-library): the functions of [`SlicedRedistribution`](@ref).

## The code

The module `ParticleRemapping.SlicedRedistribution` implements the construction for 1D-1V markers
and, for ``d \times N`` markers, along each axis ([per-axis version](@ref sliced-per-axis)). The
method type [`SlicedQuantileRemap`](@ref) implements [`remap`](@ref ParticleRemapping.remap),
[`remap_states`](@ref ParticleRemapping.remap_states) and
[`remap_frames`](@ref ParticleRemapping.remap_frames); its keywords set the box, the numbers of
slices, the relaxation parameters and mollifier widths of the two lifts, and the kernel widths of
the multipliers. The building blocks are public: the slices and targets
([`sorted_slices`](@ref), [`uniform_quantile_targets`](@ref)), the lifts
([`split_quantile_ot_lift`](@ref), [`dual_quantile_ot_lift`](@ref), [`per_axis_sliced_lift`](@ref)),
the frozen plan ([`sliced_lift_plan`](@ref), [`apply_sliced_lift`](@ref)), the squeeze and its
multipliers ([`canonical_weight_squeeze`](@ref), [`squeeze_multipliers`](@ref)) and the
diagnostics ([`phase_cell_statistics`](@ref)).

The script `scripts/sliced/verify_sliced_redistribution.jl` runs the verification suites and
prints their numbers; `scripts/sliced/long_time_experiment.jl` runs the variants and the remaps
during the long-time run; `scripts/sliced/make_figures.jl` writes the figures and the animations
to `results/sliced/`. The [verification](@ref sliced-verification) lists every test file and
script of the chapter.
