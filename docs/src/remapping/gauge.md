# [Velocity rezoning by the gauge step](@id remap-gauge)

```@meta
CurrentModule = ParticleRemapping.Remapping
```

The gauge rezoning uses the [gauge step](@ref gauge-step) alone, the factorisation with
``\psi = \mathrm{id}``. It moves no marker in ``(x, p)``; it redistributes the velocities along
their rays and lets the weights follow.

## The map and its exact properties

Write ``\lambda_a = e^{\mu_a}``, and keep the momentum as the primary variable:

```math
X_a = x_a , \qquad P_a = p_a , \qquad W_a = w_a \, e^{-\mu_a} , \qquad
V_a = \frac{p_a}{m W_a} = e^{\mu_a} v_a .
```

Each marker moves along its velocity ray at a fixed position. A marker with ``v_a = 0`` does not
move, and no marker changes the sign of its velocity. A stretch ``\mu_a > 0`` speeds the marker up
and lightens it; ``\mu_a < 0`` slows it down and makes it heavier.

!!! info "Proposition (exact properties of the rezoning)"
    For every choice of ``\mu \in \mathbb{R}^N``, also when ``\mu`` depends on the whole ensemble,
    the rezoning map has the following properties.

    1. It preserves ``\omega_N``, and it leaves every ``x_a`` and ``p_a`` unchanged.
    2. It leaves the grid current ``j_i = (q/m) \sum_a p_a S_i(x_a)`` and the total momentum
       ``\sum_a p_a`` unchanged.
    3. For every test function ``\varphi(x, v)``,
       ```math
       \sum_a W_a \, \varphi(X_a, V_a) - \sum_a w_a \, \varphi(x_a, v_a)
       = \sum_a w_a \Big( e^{-\mu_a} \varphi(x_a, e^{\mu_a} v_a) - \varphi(x_a, v_a) \Big) .
       ```
    4. In particular, the grid moments change by
       ```math
       \mathcal{M}_i^{(k)}[W, V] - \mathcal{M}_i^{(k)}[w, v]
       = \sum_a w_a \, S_i(x_a) \, v_a^k \big( e^{(k-1) \mu_a} - 1 \big) =: c_i^{(k)}(\mu) ,
       ```
       so the first moment, ``k = 1``, does not change for any ``\mu``.

!!! details "Proof"
    (1) is the [factorisation lemma](@ref factorisation-lemma) with ``\psi = \mathrm{id}``.
    (2) follows from (1), since ``\mathcal{M}_i^{(1)} = \frac{1}{m} \sum_a p_a S_i(x_a)`` depends
    only on ``x`` and ``p``. (3): insert ``X_a = x_a``, ``V_a = e^{\mu_a} v_a`` and
    ``W_a = w_a e^{-\mu_a}`` into the first sum. (4): take ``\varphi = S_i(x) \, v^k`` in (3);
    then ``e^{-\mu_a} S_i(x_a) (e^{\mu_a} v_a)^k - S_i(x_a) v_a^k
    = S_i(x_a) v_a^k (e^{(k-1)\mu_a} - 1)``. For ``k = 1`` the bracket is ``e^0 - 1 = 0``.

The change of the moments is [`gauge_constraints`](@ref). For random stretches the current is
kept to round-off and the other moments change by ``c^{(k)}``:

```@example gauge
using ParticleRemapping, LinearAlgebra, Random
Random.seed!(1)
grid = PeriodicGrid(3.0, 4)
mk = Markers(3 .* rand(40), randn(40), 0.2 .+ rand(40); m = 1.3)
μ = 0.4 .* randn(40)
after = from_canonical(mk.x, momentum(mk), mk.w .* exp.(-μ); m = mk.m)
prob = GaugeRezoningProblem(grid, mk; K = (0, 2))
c, G = gauge_constraints(prob, μ)
(current = norm(grid_moment(grid, after, 1) - grid_moment(grid, mk, 1)),
    charge = norm(c[1:4] - (grid_moment(grid, after, 0) - grid_moment(grid, mk, 0))),
    energy = norm(c[5:8] - (grid_moment(grid, after, 2) - grid_moment(grid, mk, 2))))
```

## Constraints

The zeroth and second moments are not automatic. Consider the constraints

```math
\begin{aligned}
c_i^{(0)}(\mu) &= \sum_a w_a S_i(x_a) \big( e^{-\mu_a} - 1 \big) = 0 , \\
c_i^{(2)}(\mu) &= \sum_a w_a S_i(x_a) \, v_a^2 \big( e^{\mu_a} - 1 \big) = 0 ,
\end{aligned}
\qquad i = 1, \dots, N_g .
```

Their consequences, step by step:

1. ``c_i^{(0)} = 0`` for all ``i`` means ``\rho_i`` does not change. The potential ``\phi_h`` solves
   the discrete Poisson equation with right-hand side ``\rho``, so ``\phi_h``, ``E_h`` and
   ``U(\rho)`` do not change either.
2. ``c_i^{(2)} = 0`` for all ``i`` means that each local kinetic energy
   ``\frac{m}{2} \mathcal{M}_i^{(2)}`` does not change. Summing over ``i`` with
   ``\sum_i S_i = 1``, the total kinetic energy
   ``\frac{m}{2} \sum_a w_a v_a^2 = \frac{m}{2} \sum_i \mathcal{M}_i^{(2)}`` does not change.
3. By the form of the Hamiltonian, ``H = \frac{m}{2} \sum_a w_a v_a^2 + U(\rho)``, its value is
   then the same before and after the remap. The kinetic energy of each marker changes,
   ``P_a^2 / (2 m W_a) = e^{\mu_a} p_a^2 / (2 m w_a)``; only the sums over the support of each
   shape function are fixed.

Further moments, for example ``k = 3`` (heat flux), can be added to the constraint set. Let
``\mathcal{K} \subset \{0, 2, 3, \dots\}`` be the set of constrained moments; ``k = 1`` never
needs a constraint. The represented density ``f`` is preserved only through these moments: the
higher velocity moments change, by the identity (3) of the proposition.

## Targets

The aim is an even spread of the velocities. Group the markers by cell
``C_c = [x_c, x_{c+1})`` and by the sign ``s`` of the velocity ([`velocity_groups`](@ref)). In a
group with ``n \ge 3`` markers, sort the speeds, ``|v_{(1)}| \le \dots \le |v_{(n)}|``, and assign
the equispaced targets

```math
|V^*_{(r)}| = |v_{(1)}| + \frac{r - 1}{n - 1} \big( |v_{(n)}| - |v_{(1)}| \big) ,
\qquad
\mu^*_{(r)} = \ln \frac{|V^*_{(r)}|}{|v_{(r)}|} , \qquad r = 1, \dots, n .
```

The assignment keeps the order of the speeds and the range of each group: the slowest and the
fastest marker keep their speed, ``\mu^*_{(1)} = \mu^*_{(n)} = 0``. Markers in groups with fewer
than three members get ``\mu^* = 0``. Other targets, for example quantiles of a smoothed estimate
of ``f``, fit into the same framework. The targets are [`gauge_targets`](@ref).

The targets themselves generally violate the constraints. The remap therefore solves

```math
\min_{\mu \in \mathbb{R}^N} \;
J(\mu) = \sum_{a=1}^{N} \frac{w_a}{2 s^2} \big( \mu_a - \mu_a^* \big)^2
\qquad \text{subject to} \qquad
c_i^{(k)}(\mu) = 0 , \quad i = 1, \dots, N_g , \; k \in \mathcal{K} .
```

The logarithmic variable ``\mu_a`` keeps ``\lambda_a`` positive. The factor ``w_a`` weighs a heavy
marker's deviation more, and it makes the stationarity condition below free of the weight scale.
The parameter ``s`` only scales the multipliers, because the objective has no other term. The
problem is a [`GaugeRezoningProblem`](@ref), its objective [`gauge_objective`](@ref).

## Optimality conditions

Each constraint is a sum of single-marker terms, so its derivative with respect to ``\mu_a`` has
one term:

```math
\mathsf{G}_{(k,i), a} := \frac{\partial c_i^{(k)}}{\partial \mu_a}
= (k - 1) \, w_a \, S_i(x_a) \, v_a^k \, e^{(k-1)\mu_a} .
```

A column of ``\mathsf{G}`` has at most ``2 |\mathcal{K}|`` non-zero entries, from the two nodes of
the marker's cell. With multipliers ``\alpha_i^{(k)}`` the Lagrangian is
``\mathcal{L}(\mu, \alpha) = J(\mu) + \sum_{i, k} \alpha_i^{(k)} c_i^{(k)}(\mu)``. Define the
multipliers interpolated to the markers,

```math
A_a^{(k)} = \sum_i \alpha_i^{(k)} S_i(x_a) ,
```

which by the form of the hat functions involves the two nodes of the marker's cell. Then
``\sum_{i} \alpha_i^{(k)} \mathsf{G}_{(k,i), a} = (k - 1) w_a A_a^{(k)} v_a^k e^{(k-1)\mu_a}``,
and the stationarity conditions decouple marker by marker:

```math
\frac{\partial \mathcal{L}}{\partial \mu_a}
= w_a \Big[ \frac{\mu_a - \mu_a^*}{s^2}
+ \sum_{k \in \mathcal{K}} (k - 1) \, A_a^{(k)} \, v_a^k \, e^{(k-1) \mu_a} \Big] = 0 .
```

Differentiating once more, the Hessian of the Lagrangian with respect to ``\mu`` is diagonal,

```math
H_a = \frac{\partial^2 \mathcal{L}}{\partial \mu_a^2}
= w_a \Big[ \frac{1}{s^2} + \sum_{k \in \mathcal{K}} (k - 1)^2 \, A_a^{(k)} \, v_a^k \,
e^{(k-1) \mu_a} \Big] ,
\qquad \frac{\partial^2 \mathcal{L}}{\partial \mu_a \partial \mu_b} = 0 \quad (a \neq b) .
```

For ``\mathcal{K} = \{0, 2\}`` these read

```math
\frac{\mu_a - \mu_a^*}{s^2} = A_a^{(0)} e^{-\mu_a} - A_a^{(2)} v_a^2 e^{\mu_a} , \qquad
\frac{H_a}{w_a} = \frac{1}{s^2} + A_a^{(0)} e^{-\mu_a} + A_a^{(2)} v_a^2 e^{\mu_a} .
```

The multipliers act as local forces on the log-stretches. A negative ``A^{(0)}`` pushes ``\mu_a``
down, which makes the marker heavier and raises the local charge; a negative ``A^{(2)}`` pushes
``\mu_a`` up, which raises the local kinetic energy. The stationarity is
[`gauge_stationarity`](@ref), the diagonal Hessian [`gauge_hessian`](@ref); the tests compare both
with automatic differentiation of the Lagrangian (`test/remapping/gauge.jl`).

!!! info "Proposition (second-order condition)"
    Let ``(\mu, \alpha)`` satisfy the stationarity conditions and the constraints, and let
    ``\mathsf{G}`` have full row rank. If ``H_a > 0`` for every marker, then ``\mu`` is a strict
    local minimiser of the constrained problem.

!!! details "Proof"
    A positive diagonal Hessian is positive definite on every subspace, in particular on
    ``\ker \mathsf{G}``, and the second-order sufficient condition holds.

The converse holds in practice. Suppose ``H_b < 0`` for one marker ``b``. The direction
``d = e_b + \delta``, where ``\delta`` is the minimum-norm solution of
``\mathsf{G}_{\cdot, \neq b} \, \delta = - \mathsf{G}_{\cdot, b}`` supported on the other markers,
lies in ``\ker \mathsf{G}``, and ``d^T H d = H_b + \sum_{a \neq b} H_a \delta_a^2``. With ``n``
markers on the nodes of marker ``b``, the entries of ``\delta`` are of order ``n^{-1}``, so
``\sum_a H_a \delta_a^2`` is of order ``n^{-1} \max_a H_a``, and the stationary point is a saddle
once ``H_b < 0`` is larger in size than that. One marker can trade its stretch against the many
others of its nodes.

## The Newton step

The stationary points solve the nonlinear system ``\nabla_\mu \mathcal{L}(\mu, \alpha) = 0``,
``c(\mu) = 0``. Newton's method linearises it at the current ``(\mu, \alpha)``, with the unknowns
``\Delta \mu`` and ``\alpha^+ = \alpha + \Delta \alpha``:

```math
\begin{aligned}
\nabla_\mu \mathcal{L}(\mu + \Delta\mu, \alpha + \Delta\alpha)
&\approx \nabla J(\mu) + \mathsf{G}^T \alpha + \operatorname{diag}(H) \, \Delta\mu
+ \mathsf{G}^T \Delta\alpha \\
&= \nabla J(\mu) + \operatorname{diag}(H) \, \Delta\mu + \mathsf{G}^T \alpha^+ , \\
c(\mu + \Delta\mu) &\approx c(\mu) + \mathsf{G} \, \Delta\mu .
\end{aligned}
```

Setting both to zero gives the Karush–Kuhn–Tucker (KKT) system

```math
\begin{pmatrix} \operatorname{diag}(H) & \mathsf{G}^T \\ \mathsf{G} & 0 \end{pmatrix}
\begin{pmatrix} \Delta \mu \\ \alpha^+ \end{pmatrix}
= - \begin{pmatrix} \nabla J(\mu) \\ c(\mu) \end{pmatrix} ,
\qquad \nabla J(\mu)_a = \frac{w_a}{s^2} ( \mu_a - \mu^*_a ) ,
```

a sparse linear system of size ``N + |\mathcal{K}| N_g`` ([`gauge_kkt_matrix`](@ref)). Because the
Hessian is diagonal, the first block row can be solved for the stretches,

```math
\Delta \mu = - \operatorname{diag}(H)^{-1} \big( \nabla J + \mathsf{G}^T \alpha^+ \big) , \qquad
\Delta \mu_a = - \frac{w_a}{H_a} \Big[ \frac{\mu_a - \mu^*_a}{s^2}
+ \sum_{k \in \mathcal{K}} (k - 1) A_a^{(k),+} v_a^k e^{(k-1)\mu_a} \Big] ,
```

and inserting this into the second block row gives the Schur-complement system for the
multipliers alone,

```math
\mathsf{G} \operatorname{diag}(H)^{-1} \mathsf{G}^T \, \alpha^+
= c(\mu) - \mathsf{G} \operatorname{diag}(H)^{-1} \nabla J(\mu) .
```

Its matrix ([`gauge_schur_matrix`](@ref)) has size ``|\mathcal{K}| N_g`` and the entries

```math
\big( \mathsf{G} \operatorname{diag}(H)^{-1} \mathsf{G}^T \big)_{(k,i),(l,j)}
= \sum_a \frac{(k - 1)(l - 1) \, w_a^2 \, S_i(x_a) S_j(x_a) \, v_a^{k + l} \,
e^{(k + l - 2)\mu_a}}{H_a} ,
```

which vanish unless the nodes ``i`` and ``j`` share a cell, ``|i - j| \le 1`` modulo ``N_g``. The
matrix is therefore periodic block-tridiagonal with blocks of size ``|\mathcal{K}|``. It is
symmetric, and positive definite when every ``H_a > 0`` and ``\mathsf{G}`` has full row rank.
Assembling it costs ``O(N)`` operations and solving it ``O(N_g)``; the step for ``\Delta\mu`` is
then explicit per marker. Both solves are [`gauge_newton_step`](@ref), with
`solver = :schur` or `solver = :kkt`.

The structure shows on 40 markers on 8 nodes, sorted by position: the KKT matrix has the
diagonal block ``\operatorname{diag}(H)`` and the sparse blocks ``\mathsf{G}``, ``\mathsf{G}^T``;
the Schur complement couples each node with its two neighbours, periodically, in blocks of size
``|\mathcal{K}| = 2`` (rows ``(k, i)`` ordered by ``k`` first).

```@example gauge
using CairoMakie
small = gauge_rezoning_setup(; Ng = 8, markers_per_node = 5)
o = sortperm(small.markers.x)
sorted = Markers(small.markers.x[o], small.markers.v[o], small.markers.w[o])
sprob = GaugeRezoningProblem(small.method, sorted)
plot_kkt_structure(sprob, zeros(length(sorted)), zeros(16))
```

A full Newton step can overshoot. Measure progress by the ``\ell^1`` merit function
([`gauge_merit`](@ref))

```math
\Psi_\nu(\mu) = J(\mu) + \nu \sum_{i, k} \frac{| c_i^{(k)}(\mu) |}{\mathcal{M}_i^{(k)}[w, |v|]} ,
\qquad \nu \ge 2 \max_{i,k} \big| \alpha^{+}_{(k,i)} \, \mathcal{M}_i^{(k)}[w, |v|] \big| ,
```

in which each constraint is scaled by the corresponding moment of ``|v|``. For this ``\nu`` the
KKT direction is a descent direction of ``\Psi_\nu``, and the step length
``t \in \{1, \frac12, \frac14, \dots\}`` is the first one with
``\Psi_\nu(\mu + t \Delta\mu) \le \Psi_\nu(\mu)``.

## [Algorithm](@id remap-gauge-algorithm)

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

## [Worked example: one constraint in closed form](@id remap-charge-example)

!!! info "Example (charge constraint at one node)"
    Take ``n = 4`` markers on the support of one node with ``S = 1`` for all of them (equivalently,
    a single global charge constraint), ``\mathcal{K} = \{0\}`` and ``s = 1``. The data and
    targets are

    | marker ``a``                       | 1     | 2        | 3        | 4     |
    |:-----------------------------------|:------|:---------|:---------|:------|
    | ``v_a``                            | 0.5   | 0.6      | 1.4      | 3.0   |
    | ``w_a``                            | 0.1   | 0.4      | 0.4      | 0.1   |
    | ``V^*_a`` (equispaced targets)     | 0.5   | 1.3333   | 2.1667   | 3.0   |
    | ``\mu^*_a = \ln (V^*_a / v_a)``    | 0     | 0.7985   | 0.4367   | 0     |

    At the targets the charge is ``\sum_a w_a e^{-\mu^*_a} = 0.6385`` instead of ``1``: the two
    middle markers are sped up and therefore lightened.

**Stationarity in closed form.** With one multiplier ``A = A^{(0)}``, the stationarity for
``\mathcal{K} = \{0\}`` is ``(\mu_a - \mu^*_a)/s^2 = A e^{-\mu_a}``. Put
``u_a = \mu_a - \mu^*_a``; then ``u_a e^{u_a} = s^2 A \, e^{-\mu^*_a} =: z_a``, so

```math
\mu_a(A) = \mu^*_a + W_0\big( s^2 A \, e^{-\mu^*_a} \big) ,
```

where ``W_0`` is the principal branch of the Lambert function ([`lambert_w0`](@ref)), the inverse
of ``u \mapsto u e^u`` on ``[-1, \infty)``. It is real for ``z_a \ge -1/e``.

**Second-order condition.** With ``A e^{-\mu_a} = u_a / s^2`` from stationarity,

```math
\frac{H_a}{w_a} = \frac{1}{s^2} + A e^{-\mu_a} = \frac{1 + W_0(z_a)}{s^2} ,
```

which is positive exactly on the principal branch, ``z_a > -1/e``, and zero at the branch point
``z_a = -1/e``, where ``W_0 = -1``. The other real branch, ``W_{-1} < -1``, gives saddles.

**The multiplier.** Since ``W_0`` increases, ``\mu_a(A)`` increases with ``A``, so the constraint
``c^{(0)}(A) = \sum_a w_a ( e^{-\mu_a(A)} - 1 )`` decreases strictly. As ``A \to \infty`` it tends
to ``- \sum_a w_a < 0``. The principal branch of every marker requires
``A \ge A_{\lim} = - e^{\min_a \mu^*_a} / (s^2 e)``. Hence a minimiser exists if and only if
``c^{(0)}(A_{\lim}) \ge 0``, and then it is unique: the root of ``c^{(0)}``, found by bisection.
Using ``e^{-\mu_a} = W_0(z_a) / (s^2 A)``, the equation for ``A`` can also be written as
``\sum_a w_a W_0(s^2 A e^{-\mu^*_a}) = s^2 A \sum_a w_a``.

**Result.** The root is ``A = -0.35537``, close to ``A_{\lim} = -1/e = -0.36788``. The remap is

| marker ``a``                | 1        | 2       | 3       | 4        |
|:----------------------------|:---------|:--------|:--------|:---------|
| ``\mu_a``                   | −0.7595  | 0.6043  | 0.1222  | −0.7595  |
| ``\lambda_a = e^{\mu_a}``   | 0.4679   | 1.8300  | 1.1300  | 0.4679   |
| ``V_a``                     | 0.2339   | 1.0980  | 1.5820  | 1.4037   |
| ``W_a``                     | 0.2137   | 0.2186  | 0.3540  | 0.2137   |
| ``H_a / w_a``               | 0.2405   | 0.8058  | 0.6855  | 0.2405   |

The charge is ``\sum_a W_a = 1`` to round-off. The two light end markers, whose targets are
``\mu^* = 0``, pay for the constraint: they slow down and become heavier, because the objective
weighs their deviation by their small weight. The kinetic energy is not constrained here and falls
from ``\sum_a w_a v_a^2 = 1.8530`` to ``\sum_a W_a V_a^2 = 1.5823``. The fast end marker even
falls behind marker 3, so a single constraint with few markers need not keep the order of the
speeds.

**The fold.** Scale the targets, ``\mu^* \to \beta \mu^*``. As ``\beta`` grows, the root ``A``
approaches ``A_{\lim}``, and ``H_a/w_a \to 0`` for the end markers. The branch ends at
``\beta_f = 1.3205``, where ``A = A_{\lim}`` and the end markers sit at the branch point,
``\mu_1 = \mu_4 = -1``. For ``\beta > \beta_f`` no stationary point with ``H_a > 0`` exists. This
is the mechanism that stops the homotopy of the algorithm. A fold needs the end markers
(``\mu^* = 0``) to be light enough: with ``w = (0.3, 0.2, 0.4, 0.1)`` there is no fold for any
target scale.

The same reduction works for the energy constraint alone: the stationarity for
``\mathcal{K} = \{2\}`` gives ``\mu_a = \mu^*_a - W_0\big( s^2 A^{(2)} v_a^2 e^{\mu^*_a} \big)``.

The example is [`charge_constraint_example`](@ref), its stationary point
[`charge_stationary_point`](@ref), the fold [`charge_fold`](@ref). Every number above is checked in
`test/remapping/gauge_examples.jl` and printed by `scripts/remapping/worked_examples.jl`:

```@example gauge
ex1 = charge_constraint_example()
s1 = charge_stationary_point(ex1)
βf = charge_fold(ex1)
(A = s1.A, A_lim = charge_multiplier_limit(ex1), μ = s1.μ, V = s1.V, W = s1.W, Hw = s1.Hw,
    charge = sum(s1.W), β_f = βf, μ_fold = charge_stationary_point(ex1, βf).μ)
```

```@example gauge
plot_charge_example(ex1)
```

*The charge-constraint example. (a) Speeds before, at the targets and after the constrained
remap; the marker area is proportional to the weight. (b) The constraint ``c^{(0)}(A)`` for three
target scales ``\beta``; the dashed line is ``A_{\lim} = -1/e``. For ``\beta > \beta_f`` the curve
has no root. (c) The curvatures ``H_a/w_a`` as functions of ``\beta``; they vanish for the end
markers at the fold ``\beta_f = 1.3205``.*

The animation follows the stationary point from ``\beta = 0`` to ``\beta_f``:

```@example gauge
animate_charge_fold("remapping-ex1-fold.mp4", ex1); nothing # hide
```

```@raw html
<video src="remapping-ex1-fold.mp4" controls loop muted playsinline style="width:100%"></video>
```

## [Worked example: eight markers on two nodes](@id remap-eight-example)

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
<video src="remapping-ex2-homotopy.mp4" controls loop muted playsinline style="width:100%"></video>
```

## [A numerical example](@id remap-gauge-example)

[`gauge_rezoning_setup`](@ref) applies the algorithm to ``N = 6400`` markers on ``N_g = 16``
nodes with ``L = 2\pi`` and ``m = 1``. It assembles and solves the full sparse KKT system, with
``\Delta\tau_0 = 0.05`` and ``\Delta\tau_{\min} = 10^{-3}``. The positions are uniform; in each
cell the velocities cluster around six random centres in ``[-4, 4]`` with a spread of ``0.08``,
an extreme form of the clustering that a long simulation produces. The weights are
``w_a \propto f(x_a, v_a) + \text{const}`` with a two-stream density ``f``
([`two_stream_density`](@ref)). The measure of evenness is the mean, over the groups of cell and
sign, of the ratio of the largest gap to the mean gap between consecutive sorted speeds
([`mean_gap_ratio`](@ref)); it is ``1`` for an equispaced group. The table gives the results; the
script `scripts/remapping/gauge_rezoning.jl` prints them, and `test/remapping/experiments.jl`
checks them.

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
optional. Each additional constraint reduces the reachable spread. Unlike the eight-marker
example, the homotopy stops at a fold: with hundreds of markers per node and clusters far from the
targets, some marker's curvature reaches zero first.

```@example gauge
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

```@example gauge
plot_gauge_experiment(mk, result)
```

*The numerical example, ``\mathcal{K} = \{0, 2\}``. (a), (b) The markers in ``(x, v)`` before
and after the remap, coloured by weight. (c) The distribution of the stretches ``\lambda_a``.
(d) The ratio of the largest to the mean gap per group of cell and sign, before against after;
the dashed line is the identity.*

The animation shows the continuation up to the fold, with ``\Delta\tau_0 = 0.01``; it reaches
the same ``\tau = 0.594`` in 61 accepted steps:

```@example gauge
fine = remap(gauge_rezoning_setup(; dτ0 = 0.01).method, mk)
animate_gauge_homotopy("remapping-gauge-homotopy.mp4", mk, fine)
(accepted = count(h -> h.accepted, fine.info.history) - 1, τ = fine.info.τ)
```

```@raw html
<video src="remapping-gauge-homotopy.mp4" controls loop muted playsinline style="width:100%"></video>
```

## Limits and extensions

- **Rays.** A marker moves only along its velocity ray. It cannot cross ``v = 0``, and it does
  not move in ``x``. A hole that lies across ``v = 0``, or a hole in ``x``, needs a canonical map
  in addition. A per-marker momentum translation ``p_a \mapsto p_a + \kappa_a`` is canonical for
  fixed ``\kappa_a``; it keeps the grid current if ``\sum_a \kappa_a S_i(x_a) = 0`` for every
  node, and the kinetic-energy constraint then includes the translation.
- **Higher dimension.** In ``d`` velocity dimensions the gauge step scales the velocity vector,
  ``V_a = \lambda_a v_a``, and keeps its direction. The moment identity and the constraints carry
  over, with the tensor moments ``v^{\otimes k}``; filling holes in angle needs a canonical map.
- **Global energy.** The local constraints ``c_i^{(2)} = 0`` can be replaced by the single global
  constraint ``\sum_i c_i^{(2)} = 0``, which is enough for exact energy conservation and leaves
  more freedom.
- **Representation of ``f``.** The weight follows from the move: ``W_a = p_a / (m V_a)``. The
  construction does not ensure that ``W_a`` is close to ``f(x_a, V_a)`` times the phase-space area
  of the marker. Targets taken from an estimate of ``f``, for example those of the
  [current-quantile matching](@ref remap-matching), combine the two aims.
