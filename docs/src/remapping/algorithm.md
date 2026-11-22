# [Algorithm](@id remap-algorithm)

```@meta
CurrentModule = ParticleRemapping.Remapping
```

This page states the gauge rezoning and the current-quantile matching step by step, each with a
worked example that follows every step in numbers. It closes with the combination of the two, in
which the matched velocities serve as targets of the gauge rezoning.

## [Gauge rezoning](@id remap-gauge-algorithm)

The point ``\mu = 0``, ``\alpha = 0`` satisfies the constraints and the stationarity conditions
with the targets ``\mu^* = 0``, and there ``H_a = w_a / s^2 > 0``. A homotopy follows the minimiser
from this point towards the full targets.

!!! info "Algorithm (gauge rezoning)"
    Input: markers ``(x_a, v_a, w_a)``, the constrained moments ``\mathcal{K}``, the scale ``s``,
    an initial step ``\Delta\tau_0`` and a smallest step ``\Delta\tau_{\min}``.

    1. Compute ``p_a = m w_a v_a``, the shape-function values ``S_i(x_a)``, the scales
       ``\mathcal{M}_i^{(k)}[w, |v|]``, and the targets ``\mu^*``. Set ``\tau = 0``, ``\mu = 0``,
       ``\alpha = 0``, ``\Delta\tau = \Delta\tau_0``.
    2. Set ``\tau' = \min(1, \tau + \Delta \tau)`` and replace ``\mu^*`` by ``\tau' \mu^*`` in
       ``J``. Starting from ``(\mu, \alpha)``, repeat:
        1. evaluate ``c(\mu)``, ``\mathsf{G}``, ``\nabla J``, and ``H`` with the current
           ``\alpha``;
        2. if some ``H_a \le 0``, reject: go to step 3 with failure;
        3. solve the Schur-complement system for ``\alpha^+`` and form ``\Delta \mu``;
        4. stop with success if
           ``\max_{i,k} |c_i^{(k)}| / \mathcal{M}_i^{(k)}[w, |v|] < 10^{-13}`` and
           ``\max_a |\Delta\mu_a| < 10^{-10}``;
        5. update ``\nu``, find the step length ``t`` on the merit function, and set
           ``\mu \leftarrow \mu + t \Delta\mu``, ``\alpha \leftarrow \alpha^+``.
       Fail if the iteration does not stop within a fixed number of steps.
    3. On success accept ``\tau \leftarrow \tau'`` with its ``(\mu, \alpha)``, and set
       ``\Delta\tau \leftarrow \min(2\Delta\tau, \Delta\tau_0)``. On failure keep the last
       accepted ``(\tau, \mu, \alpha)`` and halve ``\Delta\tau``.
    4. Return to step 2 until ``\tau = 1`` or ``\Delta\tau < \Delta\tau_{\min}``.
    5. Apply the rezoning with the last accepted ``\mu``: ``W_a = w_a e^{-\mu_a}``, and
       ``V_a = p_a / (m W_a)``, so that ``P_a = p_a`` holds to round-off.

The method is [`GaugeRezoning`](@ref) with `remap`; step 2 is [`gauge_newton`](@ref), steps 2–4
are [`gauge_homotopy`](@ref). The tolerances `tol` and `steptol`, the largest number of Newton
iterations `maxit` (30) and the linear solver are fields of the method.

```@example gauge
using ParticleRemapping, CairoMakie, Random
Random.seed!(1)
mk = Markers(2π .* rand(400), randn(400), 0.5 .+ rand(400))
method = GaugeRezoning(PeriodicGrid(2π, 8); K = (0, 2), dτ0 = 0.25)
r = remap(method, mk)
(τ = r.info.τ, min_curvature = r.info.min_curvature,
    gap_before = mean_gap_ratio(method.grid, mk.x, mk.v),
    gap_after = mean_gap_ratio(method.grid, r.markers.x, r.markers.v),
    report = conservation_report(mk, r.markers; grid = method.grid))
```

The continuation stops before ``\tau = 1`` when the solution branch has a fold: some ``H_a`` tends
to zero, and by the argument after the second-order proposition no nearby minimiser continues the
branch. The value of ``\tau`` reached measures how much of the target spread the constraints
admit. The [first worked example](@ref remap-charge-example) shows such a fold in closed form.
This fold belongs to the solution branch of the optimisation; it is not the
[velocity–momentum fold](@ref velocity-momentum-fold) of the sheet map.

### [Worked example: eight markers on two nodes](@id remap-eight-example)

!!! info "Example (the algorithm step by step)"
    Take the periodic interval ``[0, 2)`` with ``N_g = 2`` nodes, ``x_1 = 0`` and ``x_2 = 1``,
    ``\Delta x = 1``, ``m = 1``, ``s = 1`` and ``\mathcal{K} = \{0, 2\}``. Eight markers with
    equal weights ``w_a = 0.25`` and positive velocities sit four in each cell. By the hat
    functions, ``S_1(x) = 1 - x`` on ``[0, 1)`` and ``S_1(x) = x - 1`` on ``[1, 2)``, and
    ``S_2 = 1 - S_1``.

**Step 1: data, targets and scales.**

| marker ``a``    | 1    | 2      | 3      | 4    | 5    | 6       | 7       | 8    |
|:----------------|:-----|:-------|:-------|:-----|:-----|:--------|:--------|:-----|
| ``x_a``         | 0.2  | 0.4    | 0.6    | 0.8  | 1.2  | 1.4     | 1.6     | 1.8  |
| ``v_a``         | 1.00 | 1.05   | 1.10   | 2.50 | 0.50 | 1.80    | 1.85    | 1.90 |
| ``S_1(x_a)``    | 0.8  | 0.6    | 0.4    | 0.2  | 0.2  | 0.4     | 0.6     | 0.8  |
| ``V^*_a``       | 1.0  | 1.5    | 2.0    | 2.5  | 0.5  | 0.9667  | 1.4333  | 1.9  |
| ``\mu^*_a``     | 0    | 0.3567 | 0.5978 | 0    | 0    | −0.6217 | −0.2552 | 0    |

Each cell is one group (all velocities are positive), and the targets spread its speeds evenly
between the slowest and the fastest. The scales are
``\mathcal{M}_1^{(0)} = \mathcal{M}_2^{(0)} = 1`` and the corresponding moments of ``v^2``. At the
targets, the constraints are
``c(\mu^*) = (c^{(0)}_1, c^{(0)}_2, c^{(2)}_1, c^{(2)}_2) = (0.0398, 0.0609, -0.0958, -0.1063)``:
the targets add charge and remove kinetic energy at both nodes.

**Why a homotopy.** A Newton step from ``\mu = 0``, ``\alpha = 0`` with the full targets,
``\tau = 1``, gives ``\max_a |\Delta \mu_a| = 0.56`` and new multipliers for which the smallest
curvature is ``H_a / w_a = -0.065``. Step 2.2 of the algorithm rejects it.

**Steps 2–4: the homotopy.** With ``\Delta\tau_0 = 1/4`` the homotopy reaches ``\tau = 1`` in four
accepted steps, ``\tau = 0.25, 0.5, 0.75, 1``, without a rejection. Each step starts from the
previous solution and converges quadratically. For the last one:

| iteration | ``\max \lvert c \rvert / \text{scale}`` | ``\max_a \lvert\Delta\mu_a\rvert`` | ``\min_a H_a/w_a`` |
|:----------|:----------------------|:-----------------------|:--------|
| 0         | ``8.7 \times 10^{-15}`` | ``1.65 \times 10^{-1}``  | 0.740   |
| 1         | ``3.5 \times 10^{-3}``  | ``1.81 \times 10^{-2}``  | 0.791   |
| 2         | ``5.8 \times 10^{-5}``  | ``3.42 \times 10^{-4}``  | 0.824   |
| 3         | ``1.95 \times 10^{-8}`` | ``1.19 \times 10^{-7}``  | 0.826   |
| 4         | ``3.9 \times 10^{-15}`` | ``1.18 \times 10^{-13}`` | 0.826   |

Every step takes the full Newton length ``t = 1``. The Schur-complement solve, here a
``4 \times 4`` system, gives the same step as the full ``12 \times 12`` KKT system to
``5.6 \times 10^{-16}``, the largest difference over all steps.

**Step 5: the result.** The multipliers are

```math
(\alpha^{(0)}_1, \alpha^{(0)}_2, \alpha^{(2)}_1, \alpha^{(2)}_2)
= (-0.1184, 0.1688, -0.0446, 0.0260) ,
```

and the remap is

| marker ``a``    | 1        | 2       | 3       | 4       | 5       | 6        | 7        | 8       |
|:----------------|:---------|:--------|:--------|:--------|:--------|:---------|:---------|:--------|
| ``\mu_a``       | −0.0336  | 0.3806  | 0.6315  | 0.0312  | 0.0977  | −0.5262  | −0.2145  | 0.0591  |
| ``\lambda_a``   | 0.9669   | 1.4632  | 1.8805  | 1.0317  | 1.1026  | 0.5908   | 0.8070   | 1.0609  |
| ``V_a``         | 0.9669   | 1.5363  | 2.0685  | 2.5792  | 0.5513  | 1.0635   | 1.4929   | 2.0157  |
| ``W_a``         | 0.2585   | 0.1709  | 0.1329  | 0.2423  | 0.2267  | 0.4231   | 0.3098   | 0.2356  |
| ``H_a/w_a``     | 0.908    | 0.971   | 1.024   | 1.185   | 1.104   | 1.087    | 0.951    | 0.826   |

Stationarity holds to ``4 \times 10^{-14}``, the four constraints to ``4 \times 10^{-15}``, the
grid current and every momentum ``p_a = m W_a V_a`` to round-off. Every ``H_a > 0``, so by the
second-order proposition the result is a strict local minimiser. The ratio of the largest to the
mean gap between sorted speeds falls from ``2.80`` to ``1.06`` in cell 1 and from ``2.79`` to
``1.07`` in cell 2. The result is close to the targets but not equal to them: the end markers,
whose targets are ``\mu^* = 0``, move by up to ``0.1`` to balance charge and energy.

The example is [`eight_marker_example`](@ref); [`eight_marker_homotopy`](@ref) runs it, with the
attempt without homotopy and a fine continuation for the figures. The numbers are checked in
`test/remapping/gauge_examples.jl` and printed by `scripts/remapping/worked_examples.jl`:

```@example gauge
data = eight_marker_homotopy()
info = data.result.info
(rejected = data.direct.iterations[end].min_curvature,
    τ = [h.τ for h in info.history], α = info.α, μ = info.μ,
    V = data.result.markers.v, W = data.result.markers.w)
```

```@example gauge
plot_eight_marker_example(data)
```

*The eight-marker example. (a) The markers in ``(x, v)`` before (open) and after (filled); the area
is proportional to the weight, the grey bars are the targets, and the dashed lines are the nodes.
(b) The scaled constraint residual during the Newton iterations of each homotopy step. (c) The
smallest curvature ``H_a/w_a`` along the homotopy; the cross marks the rejected step without
homotopy.*

The animation follows the markers from ``\tau = 0`` to ``\tau = 1``:

```@example gauge
animate_eight_marker_homotopy("remapping-ex2-homotopy.mp4", data); nothing # hide
```

```@raw html
<video src="./remapping-ex2-homotopy.mp4" controls loop muted playsinline style="width:100%"></video>
```

## [Current-quantile matching](@id remap-matching)

The discrete version replaces the integrals of the cumulative currents by cumulative sums over the
markers of a cell. Let ``\hat{f}`` be the density to be preserved: a known function, or an
estimate from the markers before the remap. Define its cell-averaged cumulative current,

```math
\hat{F}_c(V) = \frac{1}{\Delta x} \int_{C_c} \int_0^V \hat{f}(X, V') \, V' \, \mathrm{d} V' \,
\mathrm{d} X .
```

For ``\hat{f} = J(x) M(v)`` it is ``\hat{F}_c(V) = \bar{J}_c \, G(V)`` with the cell mean
``\bar{J}_c`` of ``J`` ([`SeparableMaxwellianTarget`](@ref), [`cell_cumulative_current`](@ref)).
The cell average of ``\Phi_X`` is a sum over markers, because
``\tilde{n} \, \mathrm{d} X \, \mathrm{d} P`` counts markers:

```math
\frac{1}{\Delta x} \int_{C_c} \Phi_X(P) \, \mathrm{d} X
= \frac{1}{m \Delta x} \iint_{C_c \times [0, P]} P' \, \tilde{n} \, \mathrm{d} X \, \mathrm{d} P'
\approx \frac{1}{m \Delta x} \sum_{X_a \in C_c, \; 0 < P_a < P} P_a .
```

Evaluated at the momentum of a marker itself, the sum counts half of that marker (the midpoint
rule). Matching the two cell averages gives the algorithm.

!!! info "Algorithm (current-quantile matching)"
    Input: markers ``(x_a, v_a, w_a)``, a canonical map ``\psi``, and ``\hat{F}_c`` for every
    cell.

    1. Compute ``p_a = m w_a v_a`` and ``(X_a, P_a) = \psi(x_a, p_a)``.
    2. For each cell ``C_c`` and each sign ``s``, sort the markers with ``X_a \in C_c`` and
       ``\operatorname{sgn} P_a = s`` by ``|P_a|``: ``|P_{(1)}| \le \dots \le |P_{(n)}|``.
    3. Form the midpoint cumulative momenta
       ``\Pi_{(r)} = \big( \sum_{j < r} |P_{(j)}| + \tfrac{1}{2} |P_{(r)}| \big) / (m \, \Delta x)``,
       and the discrete compatibility residual
       ``R_c^s = \sum_j |P_{(j)}| / \big( m \, \Delta x \, \hat{F}_c(s \infty) \big) - 1``.
    4. Close the discrete problem by the rescaling ``\Pi_{(r)} \leftarrow \Pi_{(r)} / (1 + R_c^s)``,
       so that the cumulative momentum ends at ``\hat{F}_c(s\infty)``.
    5. Set ``|V_{(r)}| = \hat{F}_c^{-1}(\Pi_{(r)})`` on branch ``s``, ``V_{(r)} = s |V_{(r)}|``
       and ``W_{(r)} = |P_{(r)}| / (m \, |V_{(r)}|)``.

The method is [`CurrentQuantileMatching`](@ref) with `remap`; steps 2–5 are
[`current_quantile_match`](@ref), steps 3–5 on one branch [`branch_quantile_matching`](@ref).

Steps 2–5 change only weights and velocities at fixed ``(X, P)``: they form a gauge step. By the
[factorisation lemma](@ref factorisation-lemma) the whole map is symplectic for ``\omega_N``
whenever ``\psi`` is canonical, and ``P_a = m W_a V_a`` holds to round-off. The rescaling in step
4 absorbs the sampling fluctuation of the cell totals. If ``R_c^s`` is large, ``\psi`` violates
the compatibility condition, and the representation of ``\hat{f}`` fails. The sliced
redistribution meets a gauge step of the same kind, with weights set from a density estimate at
fixed momenta, and shows what it does to the velocities there
([reweighting at fixed momenta](@ref sliced-failure-reweighting)).

On ``4000`` markers sampled from ``f = (1 + 0.4 \cos x) M(v)`` in ``8`` cells, the matching after
``\psi = \mathrm{id}`` has compatibility residuals of the size of the sampling noise and keeps every
momentum to round-off:

```@example transport
using ParticleRemapping
small = quantile_matching_setup(; N = 4000, Ng = 8)
r = remap(small.identity, small.markers)
(max_residual = r.info.max_residual,
    momentum = maximum(abs.(momentum(r.markers) .- momentum(small.markers)) ./
                       abs.(momentum(small.markers))))
```

### [Worked example: five markers in one cell](@id remap-five-example)

!!! info "Example (the matching step by step)"
    Take one cell with ``\Delta x = 1``, ``m = 1``, ``\hat{f} = M(v)`` (so ``\bar{J}_c = 1`` and
    ``\hat{F}_c = G``), and five markers on the positive branch with equal weights ``w_0 = 0.1``
    and velocities ``v = (0.2, 0.5, 0.75, 1.1, 1.45)``, as a canonical map has left them. The
    positive branch represents ``\int_0^\infty M = 1/2``, which five markers of weight ``0.1``
    carry exactly.

1. The momenta are ``|P| = m w_0 v = (0.02, 0.05, 0.075, 0.11, 0.145)``, already sorted.
2. The midpoint cumulative momenta are ``\Pi = (0.01, 0.045, 0.1075, 0.2, 0.3275)``; for example
   ``\Pi_{(3)} = 0.02 + 0.05 + 0.075/2 = 0.1075``.
3. The total is ``\sum |P| / (m \Delta x) = 0.4``, against
   ``\hat{F}_c(\infty) = 1/\sqrt{2\pi} = 0.39894``, so ``R = 0.00265``.
4. The closed quantiles are ``\Pi / (1 + R) = (0.00997, 0.04488, 0.10722, 0.19947, 0.32663)``.
5. By the inverse of ``G``, for example
   ``V_{(5)} = \sqrt{-2 \ln (1 - \sqrt{2\pi} \cdot 0.32663)} = 1.8482``, and in total

| ``r``         | 1       | 2       | 3       | 4       | 5       |
|:--------------|:--------|:--------|:--------|:--------|:--------|
| ``v_{(r)}``   | 0.2     | 0.5     | 0.75    | 1.1     | 1.45    |
| ``V_{(r)}``   | 0.2250  | 0.4886  | 0.7912  | 1.1774  | 1.8482  |
| ``W_{(r)}``   | 0.0889  | 0.1023  | 0.0948  | 0.0934  | 0.0785  |

The markers sit at the current quantiles of ``M``: ``\hat{F}_c(V_{(r)}) = \Pi_{(r)}/(1 + R)``
to ``2 \times 10^{-17}``, while the old velocities miss them by up to ``0.067``. The momenta are
unchanged, so the current is exact. The charge and the kinetic energy are reproduced only to the
accuracy of the quantile rule: ``\sum W = 0.458`` and ``\sum W V^2 = 0.486``, against ``1/2`` for
both, while the old markers give ``0.5`` and ``0.417``. With five markers the quantile rule is
crude; its error decreases with the number of markers per cell and branch, and the
[numerical example](@ref remap-transport-example) has about ``6000``. Where the charge and the energy must be exact, the matched
velocities serve as targets for the [gauge rezoning](@ref remap-gauge).

The example is [`quantile_matching_example`](@ref); the numbers are checked in
`test/remapping/transport_examples.jl` and printed by `scripts/remapping/worked_examples.jl`:

```@example transport
e = quantile_matching_example()
(Π = e.Π, R = e.R, Πc = e.Πc, V = e.V, W = e.W, charge = sum(e.W),
    energy = sum(e.W .* e.V .^ 2),
    miss = maximum(abs.(e.Jbar .* maxwellian_current.(e.v) .- e.Πc)))
```

```@example transport
using CairoMakie
plot_matching_example(e)
```

*The five-marker example. (a) The cumulative current ``\hat{F}(V)``; each marker's closed quantile
(dashed horizontal) is mapped to its new velocity (filled), the open circles are the old
velocities. (b) The weights before (grey) and after the matching.*

The animation adds the markers one at a time:

```@example transport
animate_matching_example("remapping-ex4-matching.mp4", e); nothing # hide
```

```@raw html
<video src="./remapping-ex4-matching.mp4" controls loop muted playsinline style="width:100%"></video>
```

## [Combining the two constructions](@id remap-comparison)

The gauge rezoning and the canonical transport differ in what they move and in what they keep
exactly; the [properties](@ref remapping-properties) of the chapter list both side by side, and
the [comparison of all constructions](../comparison.md) sets them beside the other chapters.

The two constructions combine. The current-quantile matching gives velocities that are consistent
with an estimate of ``f``; used as [targets](@ref remap-gauge-algorithm) for the gauge rezoning,
they give a remap whose grid charge and kinetic energy are exact. Its grid current is exact only if
``\psi`` preserves the grid current, the discrete form of the compatibility condition.

`remap(::GaugeRezoning, markers; targets)` takes such targets as log-stretches
``\mu^*_a = \ln(V_a / v_a)``. Here the matched velocities of ``\psi = \mathrm{id}`` serve as
targets for the charge- and energy-conserving rezoning; the rezoning keeps the momenta, so the
matched velocities with the sign of ``v`` are admissible targets:

```@example comparison
using ParticleRemapping
grid = PeriodicGrid(2π, 8)
target = SeparableMaxwellianTarget(grid, 0.4)
mk = quantile_matching_setup(; N = 8000, Ng = 8).markers
matched = remap(CurrentQuantileMatching(grid, (x, p) -> (x, p), target), mk).markers
μstar = log.(matched.v ./ mk.v)
r = remap(GaugeRezoning(grid; K = (0, 2), dτ0 = 0.25), mk; targets = μstar)
(τ = r.info.τ, report = conservation_report(mk, r.markers; grid))
```
