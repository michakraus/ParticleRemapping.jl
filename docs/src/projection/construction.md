# [Construction](@id projection-construction)

```@meta
CurrentModule = ParticleRemapping.Projection
```

The redistribution map is the composition

```math
\mathcal{R} = \varrho \circ R_\theta, \qquad
R_\theta = D_S \circ K_S \circ \dots \circ D_1 \circ K_1 ,
```

of ``S`` kick–drift stages ``K_s, D_s`` and the reweighting ``\varrho``, a
[gauge step](@ref gauge-step). The parameters ``\theta`` follow from a fit to target points. The
steps are described in the order of the [algorithm](@ref projection-algorithm). The pictures on this page use the markers of the
[hole-filling experiment](@ref projection-hole-filling): a ``48 \times 48`` lattice in
``[0, 2\pi) \times [-4, 4]`` with the weights ``w = c f``,
``f = (1 + 0.3 \cos x)(2\pi)^{-1/2} e^{-v^2/2}``, and an elliptic hole around ``(x_h, 1.2)``
([`hole_filling_problem`](@ref)).

```@example construction
using ParticleRemapping, CairoMakie, Random
Random.seed!(20261119)
problem = hole_filling_problem(; hole_x = π)
nothing # hide
```

## [Projection](@id projection-projection)

Let ``\{\Lambda_\alpha\}`` be a tensor-product B-spline basis on ``\Omega``, periodic in ``x`` and
clamped in ``v``, with ``\Lambda_\alpha \geq 0`` and ``\sum_\alpha \Lambda_\alpha = 1``. With
``y_a = (x_a, v_a)``, the lumped projections of the represented density and of the marker
density are

```math
f_h = \sum_\alpha \frac{\sum_a w_a \Lambda_\alpha(y_a)}{\int_\Omega \Lambda_\alpha} \, \Lambda_\alpha , \qquad
g_h = \sum_\alpha \frac{\sum_a \Lambda_\alpha(y_a)}{\int_\Omega \Lambda_\alpha} \, \Lambda_\alpha .
```

Both are nonnegative, and the partition of unity gives ``\int_\Omega f_h = \sum_a w_a`` and
``\int_\Omega g_h = N``. An ``L^2`` projection with the mass matrix is more accurate but can take
negative values, which the next step does not tolerate in ``g_h``.

The code uses ``32 \times 33`` bilinear hat functions in the scaled variables ``(x/L, v/V)``
([`HatProjection`](@ref), [`lumped_projection`](@ref)).

```@example construction
plot_lumped_projection(problem)
```

(a) The markers, coloured and sized by weight, (b) the lumped projection ``f_h`` of the represented
density and (c) the lumped marker density ``g_h``; the hole (dashed) shows in both.

## [Target points](@id projection-targets)

The target points come from a density-interpolation flow in ``(x,v)``, a construction due to
[Moser:1965](@citet). Let ``\bar{g} = N / |\Omega|`` and

```math
g_t = (1 - t) \, g_h + t \, \bar{g} , \qquad t \in [0, 1] .
```

A velocity field ``u_t`` whose flow ``\chi_t`` pushes ``g_h`` forward to ``g_t`` satisfies
``\partial_t g_t + \operatorname{div} (g_t u_t) = 0``, that is
``\operatorname{div} (g_t u_t) = g_h - \bar{g}``. The ansatz ``g_t u_t = \nabla \phi`` with a
time-independent potential gives the Poisson problem

```math
\Delta \phi = g_h - \bar{g} \quad \text{in } \Omega, \qquad
\phi \text{ periodic in } x, \qquad
\partial_n \phi = 0 \text{ on } |v^i| = V .
```

Its compatibility condition ``\int_\Omega (g_h - \bar{g}) = 0`` holds by the normalisation of the
projection. The Laplacian acts in scaled variables ``(x / L, v / V)``, so that the flow treats
both directions alike. The velocity field is

```math
u_t = \frac{\nabla \phi}{(1-t) \, g_h + t \, \bar{g}} ,
```

which is defined for ``t > 0`` because ``g_t \geq t \bar{g}``, and at ``t = 0`` at every point
where ``g_h > 0``, in particular at every marker. The markers flow from dense to sparse regions:
where ``g_h > \bar{g}``, the flux ``g_t u_t = \nabla \phi`` has positive divergence, so markers
leave the region. The target points are

```math
y^*_a = (x^*_a, v^*_a) = \chi_1(y_a),
```

computed with a few steps of an explicit Runge–Kutta method. In the continuum, ``\chi_1`` pushes
``g_h`` exactly to ``\bar{g}``. For the markers, the targets are equidistributed at the resolution
of the spline basis. The map ``\chi_1`` is not symplectic; it changes volume by design. It only
provides the targets.

The target weights and target momenta are

```math
w^*_a = c \, f_h(y^*_a), \qquad
c = \frac{\sum_j w_j}{\sum_j f_h(y^*_j)}, \qquad
p^*_a = m \, w^*_a \, v^*_a .
```

The normalisation ``c`` makes the total weight exact, and ``c \approx |\Omega| / N``.

The code solves the Poisson problem with the five-point Laplacian on the nodes of the hats
(periodic in ``x``, Neumann in ``v`` through a ghost node; [`moser_laplacian`](@ref)) and
integrates the flow with 40 steps of the classical Runge–Kutta method ([`moser_targets`](@ref)).
[`projection_targets`](@ref) returns the target points and the target weights, from ``f_h`` or
from a given density.

```@example construction
targets = projection_targets(HatProjection(32, 33), problem.x0, problem.v0, problem.w0;
    problem.L, problem.V, density = (xs, vs) -> problem.density.(xs, vs))
plot_target_flow(problem, targets; window = (problem.hole.x - 1.2, problem.hole.x + 1.2, 0, 2.4))
```

Each marker, coloured and sized by weight, moves to its target (open circle) along the segment;
near the hole the markers flow into it.

```@example construction
history = moser_flow_history(problem)
animate_stages("moser-targets.mp4", [fill(history[1], 10); history; fill(history[end], 10)];
    period = (0, 2π), decorate! = (ax, _) -> plot_hole_outline!(ax, problem.hole))
nothing # hide
```

```@raw html
<video src="./moser-targets.mp4" controls loop muted playsinline style="width:100%"></video>
```

The density-interpolation flow from ``t = 0`` to ``t = 1``
([`moser_flow_history`](@ref)), with the old weights; the flow moves the points ``(x, v)`` and
changes volume, and it only provides the targets.

## [The remap Hamiltonian](@id projection-hamiltonian)

During the remap, the positions are lifted to ``\mathbb{R}^d``: each marker keeps the
representative that it has at the start, and no position is wrapped back into the periodic box
until the remap ends. Let

```math
w_{\mathrm{tot}} = \sum_{j} w_j, \qquad \alpha_j = \frac{w_j}{w_{\mathrm{tot}}}, \qquad
X_{\mathrm{c}}(x) = \sum_{j} \alpha_j \, x_j
```

be the total weight, the weight fractions and the weighted centroid of the lifted positions. For
each stage ``s``, let ``\kappa_s(\xi, w)`` be periodic in ``\xi \in \mathbb{T}^d``, and let
``\eta_s(v, w)`` be defined on ``[-V,V]^d``. The kick and drift Hamiltonians are

```math
H^K_s(x) = \sum_{a} \kappa_s\big(x_a - X_{\mathrm{c}}(x), w_a\big) , \qquad
H^D_s(p) = \sum_{a} \eta_s\big(p_a / (m w_a), w_a\big) .
```

The kick Hamiltonian depends only on the positions and the drift Hamiltonian only on the momenta.
Their exact time-one flows are therefore explicit:

```math
\begin{aligned}
K_s: \quad p_a &\mapsto p_a - \nabla_\xi \kappa_s(x_a - X_{\mathrm{c}}, w_a)
  + \alpha_a \sum_{j} \nabla_\xi \kappa_s(x_j - X_{\mathrm{c}}, w_j) , \\
D_s: \quad x_a &\mapsto x_a + \frac{1}{m w_a} \nabla_v \eta_s\big(p_a / (m w_a), w_a\big) .
\end{aligned}
```

The kick follows from ``\partial X_{\mathrm{c}} / \partial x_a = \alpha_a I``, so

```math
\nabla_{x_a} H^K_s = \nabla_\xi \kappa_s(x_a - X_{\mathrm{c}}, w_a) - \alpha_a \sum_j \nabla_\xi \kappa_s(x_j - X_{\mathrm{c}}, w_j) .
```

The code implements the two flows as [`centroid_kick`](@ref) and [`velocity_drift`](@ref), for
`d × N` matrices of positions and momenta.

!!! info "Proposition (the kick–drift map)"
    For every fixed ``\theta``:

    1. ``R_\theta`` is symplectic on ``(\mathcal{Z}, \omega)``;
    2. ``R_\theta`` conserves the total momentum ``P`` exactly;
    3. in terms of velocities, the last term of the kick adds the same vector

       ```math
       \frac{1}{m w_{\mathrm{tot}}} \sum_j \nabla_\xi \kappa_s(x_j - X_{\mathrm{c}}, w_j)
       ```

       to every velocity, a Galilean boost;
    4. a pair ``D_s \circ K_s`` is one step of symplectic Euler for ``H^K_s + H^D_s``; the
       artificial time step is absorbed into ``\kappa_s`` and ``\eta_s``.

!!! details "Proof"
    (1) Each stage is the exact flow of a Hamiltonian, hence symplectic, and so is their
    composition. (2) Under ``x_a \mapsto x_a + \delta`` for all ``a``, the centroid moves to
    ``X_{\mathrm{c}} + \delta``, so ``H^K_s`` is invariant; ``H^D_s`` does not depend on ``x``. By
    Noether's theorem, both flows conserve ``P``. Directly,
    ``\sum_a \nabla_{x_a} H^K_s = (1 - \sum_a \alpha_a) \sum_j \nabla_\xi \kappa_s = 0``, and the
    drift does not change ``p``. (3) The term ``\alpha_a \sum_j \nabla_\xi \kappa_s`` divided by
    ``m w_a`` is ``(m w_{\mathrm{tot}})^{-1} \sum_j \nabla_\xi \kappa_s``, independent of ``a``.
    (4) For a separable Hamiltonian ``H^K(x) + H^D(p)``, symplectic Euler with step ``\tau`` is
    the kick with ``\tau H^K`` followed by the drift with ``\tau H^D``.

For two random stages on six markers in ``d = 1`` and ``d = 2``, the largest entry of
``DR_\theta^\top \mathbb{J} DR_\theta - \mathbb{J}`` is at most ``2.9 \cdot 10^{-15}``,
``|\Delta P| \le 8.9 \cdot 10^{-16}``, the boost differs between the markers by
``2.8 \cdot 10^{-17}``, and a kick and a drift agree with symplectic Euler of ``H^K + H^D`` to
``3.1 \cdot 10^{-15}`` (`test/projection/stages.jl`).

!!! note "Remark (the centroid shift is necessary)"
    The kick ``p_a \mapsto p_a - F_a(x)`` with the mean force removed at unshifted positions,
    ``F_a = \nabla \kappa(x_a) - \alpha_a \sum_j \nabla \kappa(x_j)``, also conserves ``P``. It is
    not symplectic. The map ``(x, p) \mapsto (x, p - F(x))`` has the Jacobian
    ``\begin{pmatrix} I & 0 \\ -DF & I \end{pmatrix}``, which is symplectic if and only if ``DF``
    is symmetric. For ``d = 1``,
    ``\partial F_a / \partial x_j = \kappa''(x_a) \, \delta_{aj} - \alpha_a \, \kappa''(x_j)``, and
    the symmetry requires ``\alpha_a \kappa''(x_j) = \alpha_j \kappa''(x_a)`` for all ``a \neq j``,
    which fails in general. With the centroid shift, ``F = \nabla H^K`` and ``DF`` is a Hessian,
    hence symmetric.

The naive kick ([`naive_kick`](@ref)) is the case that must fail: on the same stages its
symplecticity defect is ``0.86`` (``d = 1``) and ``1.6`` (``d = 2``), while
``|\Delta P| \le 8.9 \cdot 10^{-16}`` (`test/projection/stages.jl`).

The kick is of the weight-dependent form of
[weight-dependent Hamiltonians](@ref weight-dependent-hamiltonians), with the single-particle
Hamiltonian

```math
\Psi^K_s(\xi, w) = \kappa_s(\xi - X_{\mathrm{c}}, w) - \frac{w}{w_{\mathrm{tot}}} \, \xi \cdot \sum_j \nabla_\xi \kappa_s(x_j - X_{\mathrm{c}}, w_j) .
```

By the compressibility bracket, its compressibility in ``(x,p)`` is

```math
\{W, \partial_w \Psi^K_s\}
= \{W, \partial_w \kappa_s\} + \frac{1}{w_{\mathrm{tot}}} \nabla_p W \cdot \sum_j \nabla_\xi \kappa_s(x_j - X_{\mathrm{c}}, w_j) .
```

The second term comes from the momentum correction. It changes the marker density in ``(x,p)``,
but by item 3 of the proposition it is a uniform shift in ``(x,v)`` and does not fill holes
there. The redistribution comes from the dependence of ``\kappa_s`` and ``\eta_s`` on ``w``. By
the corollary [holes persist](@ref holes-persist), without this dependence the stages could not
close a hole in ``(x,p)``.

## [Parametrisation](@id projection-parametrisation)

The generating functions are tensor-product B-splines,

```math
\kappa_s(\xi, w) = w \sum_{i, l} \lambda^s_{il} \, B_i(\xi) \, C_l(\ln w) , \qquad
\eta_s(v, w) = w \sum_{j, l} \mu^s_{jl} \, \hat{B}_j(v) \, C_l(\ln w) ,
```

with periodic splines ``B_i`` on ``\mathbb{T}^d``, splines ``\hat{B}_j`` on ``[-V,V]^d``, and
splines ``C_l`` on an interval that contains the logarithms of all weights. The logarithm is used
because the weights span several orders of magnitude between the bulk and the tail of the
distribution. The prefactor ``w`` makes the velocity change of a stage independent of the size of
the weight: with ``\ell = \ln w``, the kick and the drift become

```math
v_a \mapsto v_a - \frac{1}{m} \sum_{i,l} \lambda^s_{il} \nabla B_i(x_a - X_{\mathrm{c}}) \, C_l(\ell_a) + b_s ,
\qquad
x_a \mapsto x_a + \frac{1}{m} \sum_{j,l} \mu^s_{jl} \nabla \hat{B}_j(v_a) \, C_l(\ell_a) ,
```

where ``b_s`` is the Galilean boost of the proposition. Without the prefactor, the velocity change
of a marker scales with ``1/w_a``, and the markers in the tail restrict the fit. The parameter
vector is ``\theta = (\lambda^s_{il}, \mu^s_{jl})_{s=1}^{S}``. The functions that are constant in
``\xi``, respectively in ``v``, generate no motion, so these combinations are removed or left to
the regularisation. The parameter ``\theta = 0`` gives ``R_\theta = \mathrm{id}``.

[`KickDriftParametrisation`](@ref) holds the cubic B-spline bases (32 periodic splines in
``\xi``, 40 in ``v``, 16 in ``\ln w`` by default, two stages: 2304 parameters) and the two
controls of the experiments; [`kick_drift_map`](@ref) applies ``R_\theta``.

## [Compressibility of the stages in ``(x,v)``](@id projection-stage-jacobians)

The [compressibility bracket](@ref weight-dependent-hamiltonians) describes the compressibility
in ``(x,p)``. The marker density that is to become uniform lives in ``(x,v)``. For ``d = 1``,
write the stages as ``v \mapsto v + \Delta v(x, \ell)`` and ``x \mapsto x + \Delta x(v, \ell)``,
and let the weights be the values of a smooth field, ``\ell_a = \ln W(x_a, v_a)``.

!!! info "Proposition (stage Jacobians)"
    With the centroid and the boost held fixed, the kick and the drift act on ``(x,v)`` with the
    Jacobian determinants

    ```math
    \det D K = 1 + \partial_\ell \Delta v \; \partial_v \ln W , \qquad
    \det D D = 1 + \partial_\ell \Delta x \; \partial_x \ln W .
    ```

!!! details "Proof"
    The kick changes only ``v``, by ``\Delta v(x, \ln W(x,v))``, so its Jacobian is triangular
    with the diagonal entries ``1`` and ``1 + \partial_\ell \Delta v \, \partial_v \ln W``. The
    drift is analogous.

Three consequences follow.

1. Without a dependence on ``\ell``, every stage is area-preserving in ``(x,v)``, and the map
   cannot change the marker density in ``(x,v)``; only the reweighting can.
2. A kick changes the velocities of all markers with the same ``(x, \ell)`` by the same amount,
   and compresses them at the rate ``\partial_v \ln W``. If ``W`` is even in ``v``, as for
   ``W = c f`` with a Maxwellian, a kick that compresses at ``+v`` expands at ``-v`` by the same
   amount.
3. A drift compresses only where ``\partial_x W \neq 0``.

At a point where ``\partial_x W = 0`` and ``W`` is even in ``v``, the separable stages can
therefore fill a hole only by opening its mirror image at ``-v``. In the
[experiment](@ref projection-hole-filling), the fit from ``\theta = 0`` at such a point ends on a
mirror-symmetric saddle of its objective, where the hole stays open; a fit from a seeded start
that breaks the symmetry fills the hole. A
stage Hamiltonian that depends on ``x`` and ``p`` together does not have this coupling, at the
price of an implicit integrator ([Non-separable stages](@ref projection-nonseparable)).

The formulas ([`kick_drift_jacobian_determinants`](@ref)) agree with the forward-mode Jacobians
to ``2.2 \cdot 10^{-16}``; at ``\pm v`` the kick of the test gives ``\det DK - 1 = \mp 0.047``;
without the dependence on ``\ell``, ``\det DK = 1`` exactly (`test/projection/stages.jl`).

## [The fit](@id projection-fit)

Let ``(x^S_a, p^S_a)`` be the image of marker ``a`` under ``R_\theta``. The parameters minimise

```math
J(\theta) = \frac{1}{2} \sum_a \bigg( \frac{|x^S_a - x^*_a|_L^2}{\ell_x^2}
  + \frac{|p^S_a / (m w^*_a) - v^*_a|^2}{\ell_v^2} \bigg) + \frac{\beta}{2} |\theta|^2 ,
```

where ``|\cdot|_L`` is the minimum-image distance on ``\mathbb{T}^d`` and ``\ell_x``, ``\ell_v``
are the spline widths in ``x`` and ``v``. The momentum residual is measured in velocity units
with the target weight, so the objective does not involve the implicit weight equation. The
regularisation keeps ``R_\theta`` close to the identity and fixes the directions of ``\theta``
that the targets do not determine.

The map ``R_\theta`` is a composition of ``2S`` explicit stages,
``z^{(r)} = \Phi_r(z^{(r-1)}; \theta_r)`` for ``r = 1, \dots, 2S``. The gradient follows from the
adjoint recursion

```math
a^{(2S)} = \frac{\partial J}{\partial z^{(2S)}}, \qquad
a^{(r-1)} = \big( D_z \Phi_r \big)^{\top} a^{(r)}, \qquad
\frac{\partial J}{\partial \theta_r} = \big( \partial_{\theta_r} \Phi_r \big)^{\top} a^{(r)} + \beta \, \theta_r .
```

For a kick, ``D_z \Phi_r = \begin{pmatrix} I & 0 \\ -\nabla^2 H^K & I \end{pmatrix}``. The Hessian
``\nabla^2 H^K`` is block diagonal up to a term of rank ``d`` that the centroid contributes, so
each adjoint stage costs ``O(N)`` operations. A drift is analogous with ``\nabla^2 H^D``, which is
block diagonal. Gauss–Newton, or L-BFGS with a continuation in ``\beta``, solves the problem; the
map is linear in ``\theta`` within each stage and nonlinear only through the composition. One to
three stages are a reasonable start.

If the kinetic energy is to be kept, the fit takes the equality constraint

```math
\sum_a \frac{|p^S_a|^2}{2 m w^*_a} = \sum_a \frac{|p_a|^2}{2 m w_a}
```

with a Lagrange multiplier. The field energy changes with the charge density in any case.

The code ([`fit_remap_parameters`](@ref)) minimises ``J`` from ``\theta = 0`` with
``\beta = 10^{-3}`` and at most 40 iterations, with the lifted positions (so
``|x^S_a - x^*_a|_L`` is the plain difference of the lifted positions and targets) and
``\ell_x = L/32``, ``\ell_v = 2V/32``, the widths of the hats. Its Jacobian is the forward-mode
Jacobian of the residual for the kicks and drifts, and comes from the linearised step equations
for the [non-separable stages](@ref projection-nonseparable). The kinetic-energy constraint is not
implemented. Two methods are available, chosen by `solver`. The default, `:levenberg_marquardt`,
is the Levenberg–Marquardt method. It stops after the budget, after 20 rejected trial steps, or
when one step lowers ``J`` by less than ``10^{-6}`` relative. `:gauss_newton` is the Gauss–Newton
method with a dogleg trust region of GeometricSolvers on the residual
``[r; \sqrt\beta\,\theta]``. It stops on GeometricSolvers' convergence tests, among them a
stationarity test, after the budget, or when its step rule fails, and returns the solver's status
with its return code. The methods [`KickDriftRemap`](@ref) and [`NonSeparableRemap`](@ref) pass
`solver` to the fit; with an integer `seed` they start it at ``\theta_0 = 10^{-3}\,\xi`` with
``\xi`` a standard normal draw of `Xoshiro(seed)` instead of at ``\theta = 0``.

## [Reweighting](@id projection-final-reweighting)

The final step is the reweighting ``\varrho`` with the target weights:

```math
x'_a = x^S_a \bmod L, \qquad p'_a = p^S_a, \qquad w'_a = w^*_a, \qquad
v'_a = \frac{p'_a}{m w'_a} .
```

This step is a [gauge step](@ref gauge-step) with the stretch ``w_a / w^*_a`` of marker ``a``: it
leaves every ``(x_a, p_a)`` unchanged, so it keeps ``\omega`` and ``P``.

The reweighting at fixed momentum passes an error of the weight on to the velocity. If
``w^*_a = c f(y^*_a)(1 + \varepsilon_a)`` with a relative error ``\varepsilon_a`` of the
projection, then ``v'_a = p'_a / (m w^*_a)`` carries the relative error
``-\varepsilon_a / (1 + \varepsilon_a)``. The lumped projection with piecewise linear functions
has a large relative error where ``f`` varies fast compared with the grid, in particular in the
tails of a Maxwellian ([hole filling](@ref projection-hole-filling)). The target weights
therefore need a projection with a small relative error, for example splines of higher degree or
a projection of ``\ln f``.
