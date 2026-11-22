# [Clustered markers and a separable density](@id remap-examples)

```@meta
CurrentModule = ParticleRemapping.Remapping
```

This page gives the two numerical examples of the chapter, each with its table, the figure of its
experiment and an animation. The first is the standard example of the gauge rezoning on the
[overview page](@ref remap-overview); the second takes the standard example of the canonical
transport to ``2 \times 10^5`` markers and adds the identity and an inadmissible map. The scripts
`scripts/remapping/gauge_rezoning.jl` and `scripts/remapping/current_quantile_matching.jl` print
the numbers of the tables, and `test/remapping/experiments.jl` checks them.

## [The gauge rezoning of clustered markers](@id remap-gauge-example)

[`gauge_rezoning_setup`](@ref) applies the [algorithm](@ref remap-gauge-algorithm) to
``N = 6400`` markers on ``N_g = 16`` nodes with ``L = 2\pi`` and ``m = 1``. It assembles and solves
the full sparse KKT system, with ``\Delta\tau_0 = 0.05`` and ``\Delta\tau_{\min} = 10^{-3}``. The
positions are uniform; in each cell the velocities cluster around six random centres in
``[-4, 4]`` with a spread of ``0.08``, an extreme form of the clustering that a long simulation
produces. The weights are ``w_a \propto f(x_a, v_a) + \text{const}`` with a two-stream density
``f`` ([`two_stream_density`](@ref)). The measure of evenness is the mean, over the groups of
cell and sign, of the ratio of the largest gap to the mean gap between consecutive sorted speeds
([`mean_gap_ratio`](@ref)); it is ``1`` for an equispaced group.

|                                       | ``\mathcal{K} = \{0, 2\}`` | ``\mathcal{K} = \{0, 2, 3\}`` | unconstrained |
|:--------------------------------------|:-------------------------|:----------------------------|:--------------|
| homotopy parameter ``\tau`` reached   | 0.594                    | 0.402                       | 1             |
| largest/mean gap (before: 77.5)       | 27.5                     | 39.2                        | 1.00          |
| charge residual                       | ``2.6 \times 10^{-15}``  | ``2.6 \times 10^{-15}``     | 0.45          |
| kinetic energy residual               | ``1.9 \times 10^{-15}``  | ``3.5 \times 10^{-15}``     | –             |
| current residual                      | ``1.2 \times 10^{-16}``  | –                           | –             |
| change of ``\mathcal{M}_i^{(3)}``     | 0.22                     | ``7.5 \times 10^{-16}``     | –             |
| range of ``\lambda_a``                | ``[0.29, 8.1]``          | –                           | –             |

*Velocity rezoning by the gauge step. Residuals are maxima over the nodes, relative to
``\mathcal{M}_i^{(k)}[w, |v|]``. At the end of the ``\mathcal{K} = \{0, 2\}`` continuation the
smallest value of ``H_a / w_a`` is ``8.5 \times 10^{-4}``, against ``s^{-2} = 1`` at the start.*

The constrained remap reduces the gap ratio by a factor of almost three and meets every
constraint to round-off. Every momentum ``p_a = m W_a V_a`` holds to a relative difference of
``2 \times 10^{-16}``, and every sign of ``v`` is kept. The unconstrained targets give a perfect
spread, but they change the grid charge by ``45\,\%``, which shows that the constraints are not
optional. Each additional constraint reduces the reachable spread. Unlike the
[eight-marker example](@ref remap-eight-example), the homotopy stops at a fold: with hundreds of
markers per node and clusters far from the targets, some marker's curvature reaches zero first.

```@example experiments
using ParticleRemapping, CairoMakie
set_theme!(remap_theme())
setup = gauge_rezoning_setup()
mk, grid = setup.markers, setup.grid
result = remap(setup.method, mk)
relres(after, k) = maximum(abs.(grid_moment(grid, after, k) .- grid_moment(grid, mk, k)) ./
                           grid_moment(grid, mk.x, abs.(mk.v), mk.w, k))
(τ = result.info.τ, min_curvature = result.info.min_curvature,
    gap = (mean_gap_ratio(grid, mk.x, mk.v), mean_gap_ratio(grid, mk.x, result.markers.v)),
    charge = relres(result.markers, 0), energy = relres(result.markers, 2),
    current = relres(result.markers, 1), third_moment = relres(result.markers, 3),
    λ = extrema(exp.(result.info.μ)))
```

The accepted points of the homotopy and the attempts that it rejected; the last rejections halve
``\Delta\tau`` until it falls below ``\Delta\tau_{\min}``:

```@example experiments
h = result.info.history
(accepted = [round(e.τ; digits = 3) for e in h if e.accepted],
    rejected = [round(e.τ; digits = 4) for e in h if !e.accepted])
```

```@example experiments
plot_gauge_experiment(mk, result)
```

*The gauge rezoning of clustered markers, ``\mathcal{K} = \{0, 2\}``. (a), (b) The markers in
``(x, v)`` before and after the remap, coloured and sized by weight. (c) The distribution of the
stretches ``\lambda_a``. (d) The ratio of the largest to the mean gap per group of cell and sign,
before against after; the dashed line is the identity.*

The animation shows the first two cells along the continuation up to the fold, with
``\Delta\tau_0 = 0.01``, which reaches the same ``\tau = 0.594`` in 61 accepted steps. Each frame
is a gauge step of the markers ([`remap_frames`](@ref ParticleRemapping.remap_frames)), so the
markers move along their rays at fixed ``x``; a marker that speeds up loses weight, and one that
slows down gains weight:

```@example experiments
fine = remap(gauge_rezoning_setup(; dτ0 = 0.01).method, mk)
sel = findall(<(2grid.Δx), mk.x)
zoom = [RemapStage(Markers(f.markers.x[sel], f.markers.v[sel], f.markers.w[sel]; m = mk.m),
            "cells 1 and 2, " * f.label)
        for f in remap_frames(setup.method, mk, fine; frames_per_step = 2)]
vl = 1.05maximum(abs, fine.markers.v[sel])
animate_stages("remapping-gauge-homotopy.mp4", zoom; window = (0, 2grid.Δx, -vl, vl))
(accepted = count(h -> h.accepted, fine.info.history) - 1, τ = fine.info.τ)
```

```@raw html
<video src="remapping-gauge-homotopy.mp4" controls loop muted playsinline style="width:100%"></video>
```

## [The canonical transport of a separable density](@id remap-transport-example)

[`quantile_matching_setup`](@ref) samples ``N = 2 \times 10^5`` markers with equal weights from
``f(x, v) = (1 + 0.4 \cos x) \, e^{-v^2/2} / \sqrt{2\pi}`` on ``[0, 2\pi)``, with ``16`` cells and
``m = 1``; this is the density of the [closed-form example](@ref remap-lift-example). The
admissible map is the cotangent lift with ``a = J_+^{-1/2}`` at ``t = 1.5``, which moves markers
by up to ``1.85`` in ``x``. The inadmissible map is the cotangent lift with
``a = 1 + 0.5 \cos x``. The error measure ([`binned_weight_error`](@ref)) is the relative ``L^1``
difference between the binned weights and the exact bin integrals of ``f``, on ``16`` cells times
``24`` velocity bins in ``[-4, 4]``.

| configuration                                  | binned ``L^1`` error | largest ``\lvert R_c^s \rvert`` |
|:-----------------------------------------------|:---------------------|:--------------------|
| original markers (sampling noise)              | 0.0273               | –                   |
| ``\psi = \mathrm{id}``, then matching          | 0.0113               | 0.038               |
| admissible ``\psi``, then matching             | 0.0117               | 0.035               |
| admissible ``\psi``, old weights kept (``W = w``) | 0.2605            | –                   |
| inadmissible ``\psi``, then matching           | 1.277                | 5.4                 |

*Canonical transport with current-quantile matching.*

After an admissible map, the matching restores the representation of ``f`` below the sampling
noise of the original markers. The error falls below that noise because the matching places the
velocities at the quantiles of the exact ``f``. Without the matching, the same map gives an error
ten times the noise. After an inadmissible map, the residual ``R_c^s`` shows the violation, and
the representation fails. The sampled marker counts and charges per cell after the admissible map
follow the closed forms ``(J(x) J(X))^{1/2}`` and ``J(X)`` of the closed-form example.

```@example experiments
transport = quantile_matching_setup()
experiment = quantile_matching_experiment(transport)
err(m) = binned_weight_error(transport.target, transport.grid, m.x, m.v, m.w)
(original = err(experiment.original),
    identity = (err(experiment.identity.markers), experiment.identity.info.max_residual),
    admissible = (err(experiment.admissible.markers), experiment.admissible.info.max_residual),
    kept = err(experiment.kept),
    inadmissible = (err(experiment.inadmissible.markers),
        experiment.inadmissible.info.max_residual))
```

```@example experiments
plot_lift_example(transport, experiment)
```

*The closed-form example and the numerical example. (a) Characteristics of the admissible lift
``a = J^{-1/2}``. (b) At ``t = 1.5``: the represented density ``J(X)``, the marker density
``(J(x)J(X))^{1/2}`` and the weight ratio ``1/\varphi'`` (lines), with the sampled marker counts
(circles) and charges (triangles) per cell, normalised. (c) The discrete compatibility residual
``R_c^s`` per cell for the admissible and the inadmissible lift. (d) The velocity marginal in
cell 5: exact, after matching, with the old weights kept, and after the inadmissible lift.*

The animation runs the admissible lift from ``t = 0`` to ``t = 1.5`` with the algorithm applied at
every frame ([`lift_matching_frames`](@ref)), and shows every 50th marker, coloured and sized by
weight: the markers move in ``x``, their weights change, and the represented charge per cell stays
at ``J``.

```@example experiments
frames = lift_matching_frames(transport)
animate_lift_matching("remapping-ex3-lift.mp4", frames, transport); nothing # hide
```

```@raw html
<video src="remapping-ex3-lift.mp4" controls loop muted playsinline style="width:100%"></video>
```
