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
spreads the positions over the whole interval, whatever their correlation with the velocity. For
frozen slices, anchors and multipliers the map satisfies ``\mathcal{R}^* \omega_N = \omega_N``; it
keeps the total momentum ``P`` and the total weight to round-off. The multipliers of the remap are
computed from the markers, so the map is symplectic for their frozen values
([parameters that depend on the state](@ref state-dependent-parameters)).

The setting is that of the foundations: the canonical coordinates ``(x, p)`` with
``p_a = m w_a v_a``, the [presymplectic form](@ref presymplectic-form) ``\omega_N``, the
[factorisation lemma](@ref factorisation-lemma) and the
[symplecticity condition](@ref symplecticity-conditions). This chapter marks the image of each
stage by its own decoration: ``(\tilde{x}, \tilde{p})`` after the first lift,
``(\hat{x}, \hat{p})`` after the second, ``(\bar{x}, \bar{p})`` after the squeeze and
``(\bar{x}, \bar{v}, \bar{w})`` for the output markers. In the code the fields `q̃`, `p̃`, `q̂`,
`p̂`, `q̄` and `p̄` of the result hold these values, and the field `q_anc` holds the position
anchors ``x^{\mathrm{anc}}``.

The map is a composition of four stages,

```math
\mathcal{R} = \mathcal{F} \circ T_\mu \circ \Psi_p \circ \Psi_q , \qquad
\psi = T_\mu \circ \Psi_p \circ \Psi_q ,
```

1. ``\Psi_q`` sorts the markers into slices of equal velocity range and moves the positions in
   each slice towards equispaced targets: the cotangent lift of the gradient of a strictly convex
   potential;
2. ``\Psi_p`` does the same with the roles exchanged: it sorts the markers into slices of equal
   position range and moves the momenta in each slice towards equispaced targets;
3. ``T_\mu`` scales each momentum by a slot multiplier ``\mu_a`` and each position offset by
   ``1/\mu_a``, a linear canonical map that keeps the total momentum;
4. ``\mathcal{F}`` sets the weights ``\bar{w}_a = w_{\mathrm{ref}} \mu_a`` and the velocities
   ``\bar{v}_a = \bar{p}_a / (m \bar{w}_a)``.

Stages 1–3 are canonical maps of ``(x, p)``; stage 4 is a [gauge step](@ref gauge-step), which
any choice of weights keeps symplectic. The order of the two lifts is that of a Knothe–Rosenblatt
rearrangement [Villani:2009](@cite): first the positions conditional on the velocity, then the
momenta conditional on the position. A lift of the marginals alone cannot fill a hole whose
marginals are uniform; the slices can.

```@example sliced-overview
using ParticleRemapping, CairoMakie
ex = checkerboard_experiment()
plot_checkerboard_filling(ex)
```

The figure shows 400 equal-weight markers that fill two of the four quadrants of ``[-1, 1]^2``
(left). Both marginals are uniform. The sliced lifts with 10 velocity and 20 position slices put
about a quarter of the markers into each quadrant (middle); one slice for each lift leaves the
empty quadrants empty (right). The counts are those of [`checkerboard_experiment`](@ref) and
the test `test/sliced/experiments.jl`.

## The pages

- [The stages and their structure](@ref sliced-structure): why coupling alone does not fill a
  hole, and the statements about the sliced lift, the weight squeeze, the composed map and the
  composition with a later flow step, with proofs.
- [Algorithm](@ref sliced-algorithm): the four stages with their formulas, a picture of each stage
  on a small example, the parameters and the cost.
- [Experiments](@ref sliced-experiments): the five verification suites, the hole filling on a
  checkerboard and the long-time run of 900 markers in an electric field, with tables, figures
  and animations.
- [What does not work](@ref sliced-failures): the constructions that fail, with the reason and a
  control that shows it.
- [Library](@ref sliced-library): the functions of [`SlicedRedistribution`](@ref).

The [verification](@ref sliced-verification) chapter lists the tests and scripts behind every
number of this chapter.
