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
Assembling it costs ``O(N)`` operations, and a periodic block-tridiagonal solve costs
``O(N_g)``; the step for ``\Delta\mu`` is then explicit per marker. Both solves are
[`gauge_newton_step`](@ref), with `solver = :schur` or `solver = :kkt`; the package factorises
the Schur complement as a dense matrix, at the cost ``O((|\mathcal{K}| N_g)^3)``
([known issue K9](@ref remapping-known-issues)).

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

## The numerical example and the limits

The [additional experiments](@ref remap-gauge-example) apply the
[algorithm](@ref remap-gauge-algorithm) to 6400 clustered
markers on 16 nodes, with two constraint sets and the unconstrained targets as a control, and
follow the homotopy up to its fold. What the gauge rezoning cannot do, and why, is collected under
[limitations](@ref remapping-limitations), and its form in ``d`` dimensions under
[higher dimensions](@ref remapping-higher-dimensions).
