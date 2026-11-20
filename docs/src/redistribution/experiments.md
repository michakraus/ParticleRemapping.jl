# [Numerical experiments](@id rd-experiments)

```@meta
CurrentModule = ParticleRemapping.Redistribution
```

The defaults of [`SplineRedistribution`](@ref) fit the verification set-up of 16 markers: there
every step of the construction does what the chapter states, and the symplectic step reaches its
targets to round-off. This page runs that example step by step, compares the schemes and the
substeps on it, and then shows two cases outside it, where the relaxation or the symplectic step
fails ([known issues](@ref redistribution-known-issues) K2 and K3). The script
`scripts/redistribution/failure_cases.jl` prints every number of this page (Julia 1.13.4, about a
minute; the second failure case takes 44 s).

## [The standard example in detail](@id rd-experiments-standard)

[`redistribution_example`](@ref) starts from ``4 \times 4`` markers at the cell centres of a
lattice and clusters them by 60 leapfrog steps in a static force:

```@example rd-experiments
using ParticleRemapping, CairoMakie
set_theme!(remap_theme())
ex = redistribution_example()
mk = ex.markers
plot_remap(ex.initial, mk; period = (0, ex.Lx),
    titles = ("(a) lattice", "(b) after 60 leapfrog steps"))
```

The leapfrog steps keep the weights. One remap with the defaults:

```@example rd-experiments
method = SplineRedistribution()
res = remap(method, mk)
i = res.info
(cv_before = i.cv_before, cv_after = i.cv_after, min_weight = minimum(i.target_weights),
 contraction = only(i.contraction), newton_residual = only(i.residual),
 target_error = i.target_error)
```

The relaxation lowers the coefficient of variation of the kernel density from ``0.2199`` to
``0.1939`` ([de-clustering relaxation](@ref rd-projection)). The KKT weights are positive, from
``0.0198`` to ``0.1232`` (before: ``0.0161`` to ``0.1318``). The contraction number of the
implicit midpoint step at its stage point is ``2.266``, above the sufficient bound ``1`` of the
[parameter solver](@ref rd-linear-parameter-solver) (K4); the Newton homotopy still reaches the
targets, to the scaled distance ``9.87 \cdot 10^{-15}``. The changes of the conserved quantities
are those of the [overview](@ref redistribution-overview): the total weight and the total momentum
are kept to round-off, and the kinetic energy changes by ``0.1323`` of ``1.0812``, because no
constraint of the construction holds it.

## [Schemes and substeps](@id rd-experiments-schemes)

The same markers with the two schemes, in one step and in substeps
([`redistribution_substeps`](@ref)):

```@example rd-experiments
runs = ["midpoint, 1 step" => SplineRedistribution(),
    "midpoint, 3 substeps" => SplineRedistribution(; substeps = 3),
    "symplectic Euler, 1 step" => SplineRedistribution(; scheme = SymplecticEulerStep()),
    "symplectic Euler, 2 substeps" =>
        SplineRedistribution(; scheme = SymplecticEulerStep(), substeps = 2)]
[(run = name, target_error = r.info.target_error, contraction = r.info.contraction,
  momentum_change = conservation_report(mk, r.markers).momentum)
 for (name, r) in ((name, remap(m, mk)) for (name, m) in runs)]
```

| run | target error | contraction number of each step |
|:--|:--|:--|
| midpoint, 1 step | ``9.87 \cdot 10^{-15}`` | ``2.266`` |
| midpoint, 3 substeps | ``1.61 \cdot 10^{-15}`` | ``0.965``, ``0.755``, ``0.619`` |
| symplectic Euler, 1 step | ``0.184`` | ``1.391`` |
| symplectic Euler, 2 substeps | ``2.18 \cdot 10^{-15}`` | ``0.782``, ``0.690`` |

Every run keeps the total momentum to round-off (at most ``1.1 \cdot 10^{-16}``). The single symplectic Euler step converges to
another root of its equations, a caustic of the [remark](@ref rd-coarse-basis) (a); the substeps
keep each contraction number below one and reach the targets. Its stages show where the single
step ends:

```@example rd-experiments
se = SplineRedistribution(; scheme = SymplecticEulerStep())
rse = remap(se, mk)
plot_stages(remap_states(se, mk, rse); targets = (rse.info.targets_x, rse.info.targets_v),
    period = (0, ex.Lx))
```

Panel (c) shows the end point of the step at the old weights; in (d) some markers sit away from
their targets (open circles) after the gauge step.

## [Outside the verification set-up](@id rd-experiments-failures)

### More markers: the relaxation leaves the velocity box (K2)

The relaxation sums the repulsion over all pairs, and its step size ``0.006`` fits 16 markers.
With 900 markers of uniform positions, standard normal velocities and equal weights, the same ten
steps move the points by up to ``3.74`` and push ``309`` of the ``900`` targets out of the
velocity box ``|v| \le 3.5`` of the spline basis, up to ``|\bar{v}| = 6.16``:

```@example rd-experiments
using Random
rng = MersenneTwister(1)
x9, v9 = 2π .* rand(rng, 900), randn(rng, 900)
r9 = declustering_relaxation(x9, v9; Lx = method.Lx, σx = method.σx, σv = method.σv)
fig = plot_target_relaxation(r9.history; Lx = method.Lx)
hlines!(content(fig[1, 1]), [-method.Vmax, method.Vmax]; color = :gray25, linestyle = :dash)
fig
```

The dashed lines bound the velocity box. `remap(SplineRedistribution(), markers)` on these
markers stops in the KKT solve with `ArgumentError: matrix contains Infs or NaNs`.

### A hole: the step misses its targets (K3)

The 244 markers of [`hole_filling_problem`](@ref ParticleRemapping.Projection.hole_filling_problem)
with `nL = 16`, `hole_a = 0.9` and `hole_b = 0.8` (a ``16 \times 16`` lattice on
``[0, 2\pi) \times [-4, 4]`` with an elliptic hole), redistributed by
`SplineRedistribution(; Vmax = 4.0, relax_step = 0.006 * 16 / 256)`, the relaxation step scaled
down with the number of pairs. The remap takes 44 s and returns without an error:

| quantity | value |
|:--|:--|
| target error (scaled) | ``1.98 \cdot 10^{8}`` |
| contraction number | ``4.85 \cdot 10^{8}`` |
| Newton residual at the end point | ``3.96 \cdot 10^{8}`` |
| change of the total weight | ``1.8 \cdot 10^{-15}`` |
| change of the total momentum | ``0.165`` |
| change of the kinetic energy | ``3394`` (of ``2.933``) |

The contraction number ``4.85 \cdot 10^{8}`` puts ``\theta^*`` far outside the regime of the
[parameter solver](@ref rd-linear-parameter-solver), and the forward Newton solve does not
converge: the end point is not a step of the scheme, so neither the targets
nor the total momentum are reached, and only the KKT weights keep the total weight. `remap` does
not report the failure; only `info.target_error` and `info.residual` show it (K3).

## What the experiments show

On the verification set-up, each step of the construction holds: the relaxation reduces the
clustering, the KKT weights are positive and keep the discrete moments, and the symplectic step
reaches the targets, in one step for the midpoint rule and in two substeps for symplectic Euler.
The contraction condition is sufficient, not necessary: the midpoint step reaches the targets with
the contraction number ``2.266``, and the single symplectic Euler step with ``1.391`` does not.
Outside the set-up, the relaxation parameters do not scale with the number of markers, and a
single step to distant targets meets no root on the identity-connected branch. The
[limitations](@ref redistribution-limitations) state the reasons.
