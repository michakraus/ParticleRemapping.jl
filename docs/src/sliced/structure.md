# [The stages and their structure](@id sliced-structure)

```@meta
CurrentModule = ParticleRemapping.SlicedRedistribution
```

The foundations state when a redistribution of weighted markers is symplectic, in the
[three forms of the condition](@ref symplecticity-conditions) (full space, leaf of constant
weights, pseudo-time flow), and why a map that moves each marker on its own keeps every hole
([uncoupled maps](@ref uncoupled-maps)). Suites 1 and 2 of the
[experiments](@ref sliced-experiments) check the three forms with ForwardDiff Jacobians and with
controls that drop the recoil or read ``w`` at fixed ``(x, p)``. This page adds what the stages of
the sliced map do: coupling the markers is not enough to fill a hole, the sliced lift is a global
symplectic diffeomorphism, the weight squeeze is a linear canonical map that keeps the total
momentum, the composed map satisfies ``\mathcal{R}^* \omega_N = \omega_N``, and a later flow step
with weights read from the state has a composition defect.

## [Coupling is necessary, not sufficient](@id sliced-coupling-not-sufficient)

A canonical map of all ``N`` markers keeps the volume of ``\mathbb{R}^{2Nd}``, not that of the
single-marker space, so a map that couples the markers is free of the
[obstruction of uncoupled maps](@ref uncoupled-maps). Coupling is necessary but not sufficient:
the one-slice lifts of the [checkerboard](@ref sliced-checkerboard) couple all markers and still
leave the empty quadrants empty. For a large number of markers in one slice, the quantile lift
approaches a map of each marker's own coordinate (its rank becomes its quantile), and the
obstruction applies in that limit; this is a heuristic reading, not a statement proved here. The
slices break it: within a slice of small velocity range the positions are spread over the whole
interval, whatever their correlation with the velocity.

## [The sliced quantile lift](@id sliced-lift-structure)

Fix a partition of the markers into slices ``S_1, \dots, S_B`` with ``n_b`` markers each, an
order ``a_1, \dots, a_{n_b}`` of each slice, targets ``\tau^{(b)}_1 < \dots < \tau^{(b)}_{n_b}``,
parameters ``\beta_{\mathrm{mean}}, \beta_{\mathrm{intra}} \in [0, 1)``, a mollifier width
``\varepsilon > 0`` and positive constants ``\kappa_b``. With the means
``\bar{u} = \frac{1}{N} \sum_c u_c`` and ``\bar{u}_b = \frac{1}{n_b} \sum_{a \in S_b} u_a``, the
**sliced quantile potential** of ``u \in \mathbb{R}^N`` is

```math
\begin{aligned}
S(u) ={}& \frac{1}{2N} \Big( \sum_{c=1}^{N} u_c \Big)^2
       + \frac{1 - \beta_{\mathrm{mean}}}{2} \sum_{b=1}^{B} n_b \, (\bar{u}_b - \bar{u})^2
       + \frac{1 - \beta_{\mathrm{intra}}}{2} \sum_{b=1}^{B} \sum_{a \in S_b} (u_a - \bar{u}_b)^2 \\
      &+ \sum_{b=1}^{B} \frac{\beta_{\mathrm{intra}}}{2 n_b \kappa_b} \sum_{r, s = 1}^{n_b}
         \big| \tau^{(b)}_r - \tau^{(b)}_s \big| \sqrt{(u_{a_r} - u_{a_s})^2 + \varepsilon^2} .
\end{aligned}
```

Its gradient is, for the marker ``a = a_r`` of slice ``b``,

```math
\tilde{u}_a = \frac{\partial S}{\partial u_a}
  = \beta_{\mathrm{mean}} \, \bar{u} + (\beta_{\mathrm{intra}} - \beta_{\mathrm{mean}}) \, \bar{u}_b
  + (1 - \beta_{\mathrm{intra}}) \, u_a
  + \frac{\beta_{\mathrm{intra}}}{n_b \kappa_b} \sum_{s=1}^{n_b} \big| \tau^{(b)}_r - \tau^{(b)}_s \big|
    \frac{u_a - u_{a_s}}{\sqrt{(u_a - u_{a_s})^2 + \varepsilon^2}} ,
```

and its Hessian is

```math
\nabla^2 S = \frac{\beta_{\mathrm{mean}}}{N} \mathbf{1} \mathbf{1}^\top
  + \sum_{b} \frac{\beta_{\mathrm{intra}} - \beta_{\mathrm{mean}}}{n_b} \mathbf{1}_b \mathbf{1}_b^\top
  + (1 - \beta_{\mathrm{intra}}) \, I + \sum_b \Lambda_b ,
```

where ``\mathbf{1}_b`` is the indicator vector of slice ``b`` and ``\Lambda_b`` is the weighted
graph Laplacian of slice ``b`` with the edge weights
``\frac{\beta_{\mathrm{intra}}}{n_b \kappa_b} |\tau^{(b)}_r - \tau^{(b)}_s| \,
\varepsilon^2 / ((u_{a_r} - u_{a_s})^2 + \varepsilon^2)^{3/2} \ge 0``.
[`sliced_quantile_potential`](@ref), [`split_quantile_ot_lift`](@ref) and
[`sliced_quantile_hessian`](@ref) compute ``S``, ``\nabla S`` and ``\nabla^2 S``; the test
`test/sliced/lift.jl` checks the last two against ForwardDiff derivatives of the first.

!!! info "Proposition (the sliced lift is a global symplectic diffeomorphism)"
    For a fixed partition, order, targets and parameters as above:
    1. ``\nabla^2 S \succeq (1 - \max(\beta_{\mathrm{mean}}, \beta_{\mathrm{intra}})) \, I`` and
       ``\nabla^2 S \, \mathbf{1} = \mathbf{1}``; every slice-constant vector orthogonal to
       ``\mathbf{1}`` is an eigenvector with the eigenvalue ``1 - \beta_{\mathrm{mean}}``.
    2. ``u \mapsto \nabla S(u)`` is a diffeomorphism of ``\mathbb{R}^N`` with
       ``\sum_a \tilde{u}_a = \sum_a u_a`` and ``\nabla S(u + c \mathbf{1}) = \nabla S(u) + c \mathbf{1}``.
    3. For any constant anchors ``v^{\mathrm{anc}} \in \mathbb{R}^N`` the anchored cotangent lift
       ```math
       \tilde{u} = \nabla S(u) , \qquad
       \tilde{v} = v^{\mathrm{anc}} + \big( \nabla^2 S(u) \big)^{-1} (v - v^{\mathrm{anc}})
       ```
       is a symplectic diffeomorphism of ``(u, v) \in \mathbb{R}^{2N}`` with
       ``\sum_a \tilde{v}_a = \sum_a v_a``, and ``\tilde{v}`` does not change under
       ``u \mapsto u + c \mathbf{1}``.
    4. The same holds for the dual lift ``\hat{p} = \nabla S(p)``,
       ``\hat{x} = x^{\mathrm{anc}} + (\nabla^2 S(p))^{-1} (x - x^{\mathrm{anc}})``.

    With ``(u, v) = (x, p)`` the lift ``\Psi_q`` keeps ``\sum_a x_a`` and the total momentum; with
    ``(u, v) = (p, x)`` the lift ``\Psi_p`` keeps the total momentum and ``\sum_a x_a``.

!!! details "Proof"
    *1.* The three subspaces ``E_0 = \mathrm{span}\{\mathbf{1}\}``, ``E_1`` (slice-constant
    vectors orthogonal to ``\mathbf{1}``) and ``E_2`` (vectors with zero mean on every slice) are
    orthogonal and span ``\mathbb{R}^N``. Each ``\Lambda_b`` is symmetric with
    ``\Lambda_b \mathbf{1}_b = 0``, so it vanishes on ``E_0 \oplus E_1`` and maps ``E_2`` into
    itself; ``\mathbf{1}\mathbf{1}^\top`` vanishes on ``E_1 \oplus E_2``;
    ``\sum_b \mathbf{1}_b \mathbf{1}_b^\top / n_b`` is the identity on ``E_0 \oplus E_1`` and zero
    on ``E_2``. Hence ``\nabla^2 S`` is block diagonal with
    ``\nabla^2 S|_{E_0} = \beta_{\mathrm{mean}} + (\beta_{\mathrm{intra}} - \beta_{\mathrm{mean}})
    + (1 - \beta_{\mathrm{intra}}) = 1``,
    ``\nabla^2 S|_{E_1} = (\beta_{\mathrm{intra}} - \beta_{\mathrm{mean}}) + (1 - \beta_{\mathrm{intra}})
    = 1 - \beta_{\mathrm{mean}}`` and
    ``\nabla^2 S|_{E_2} = (1 - \beta_{\mathrm{intra}}) I + \sum_b \Lambda_b \succeq
    (1 - \beta_{\mathrm{intra}}) I``, since the Laplacians are positive semidefinite (non-negative
    edge weights). All three bounds are at least ``1 - \max(\beta_{\mathrm{mean}},
    \beta_{\mathrm{intra}}) > 0``.

    *2.* With ``c_0 = 1 - \max(\beta_{\mathrm{mean}}, \beta_{\mathrm{intra}})``, ``S`` is strongly
    convex: ``(\nabla S(u) - \nabla S(u'))^\top (u - u') \ge c_0 |u - u'|^2``, so ``\nabla S`` is
    injective. For ``y \in \mathbb{R}^N`` the function ``S(u) - y^\top u`` grows at least like
    ``\frac{c_0}{2}|u|^2`` minus a linear function, so it has a minimiser, where
    ``\nabla S(u) = y``: ``\nabla S`` is onto. Its Jacobian ``\nabla^2 S`` is invertible, so the
    inverse is smooth. The sum: ``\mathbf{1}^\top \nabla S(u) = \frac{\mathrm{d}}{\mathrm{d}t}
    S(u + t \mathbf{1})|_{t=0}``; only the first term of ``S`` changes under the shift, by
    ``t \sum_c u_c + N t^2 / 2``, so ``\mathbf{1}^\top \nabla S(u) = \sum_c u_c``. The same
    expansion gives ``\nabla S(u + c \mathbf{1}) = \nabla S(u) + c \mathbf{1}`` and
    ``\nabla^2 S(u + c \mathbf{1}) = \nabla^2 S(u)``.

    *3.* The cotangent lift ``(u, v) \mapsto (\phi(u), \mathrm{D}\phi(u)^{-\top} v)`` of a
    diffeomorphism ``\phi`` is symplectic [HairerLubichWanner:2006](@cite); here
    ``\mathrm{D}\phi = \nabla^2 S`` is symmetric. The anchored lift is this lift between the
    translations ``v \mapsto v - v^{\mathrm{anc}}`` and ``v \mapsto v + v^{\mathrm{anc}}``, which
    are symplectic because the anchors are constants. By 1, ``(\nabla^2 S)^{-1} \mathbf{1} = \mathbf{1}``
    and, by symmetry, ``\mathbf{1}^\top (\nabla^2 S)^{-1} = \mathbf{1}^\top``, so
    ``\mathbf{1}^\top \tilde{v} = \mathbf{1}^\top v^{\mathrm{anc}} + \mathbf{1}^\top (v - v^{\mathrm{anc}})
    = \mathbf{1}^\top v``. The invariance of ``\tilde{v}`` follows from that of ``\nabla^2 S``.

    *4.* The exchange ``(x, p) \mapsto (p, -x)`` is symplectic. The dual lift is this exchange,
    the lift of 3 with ``(u, v) = (p, -x)`` and anchors ``-x^{\mathrm{anc}}``, and the exchange
    back.

!!! note "Slices that depend on the state"
    The sliced map reads its slices, their order and the anchors from the input state: velocity
    slices from the sorted velocities, position slices from the sorted output of ``\Psi_q``. These
    are piecewise constant functions of the state. On every open set where no two markers swap
    their order the map is the smooth map of the proposition, so it is symplectic there. It jumps
    where the order changes. Whether the map with state-dependent slices is injective across these
    sets is an open question. A canonical map whose parameters are read from the data is a weak
    statement for a finite marker set ([finite marker sets](@ref finite-marker-sets)); the content
    of the proposition is that each lift is one smooth map of all markers for its frozen data.

**Quantile limit.** For ``\beta_{\mathrm{mean}} = \beta_{\mathrm{intra}} = 1`` and
``\varepsilon \to 0`` with a slice ordered like its centred targets, the mollified sign tends to
``\mathrm{sign}(r - s)``, ``\kappa_b \to 1`` and the sum becomes
``\sum_s (\tau_r - \tau_s) = n_b \tau_r``: ``\tilde{u}_{a_r} \to \bar{u} + \tau^{(b)}_r``, the
quantile map onto the targets, shifted to the mean. This one-dimensional quantile map is the
optimal transport map of the slice for a convex cost [Villani:2009](@cite). For
``\beta < 1`` the lift relaxes towards it: the slice means by the factor ``\beta_{\mathrm{mean}}``
towards the global mean, the deviations within a slice by ``\beta_{\mathrm{intra}}`` towards the
targets. The constant ``\kappa_b`` is the positive root of ``\kappa^2 + \rho\kappa - 1 = 0`` with
``\rho = \varepsilon / (\tfrac{1}{2}(\tau^{(b)}_{n_b} - \tau^{(b)}_1))``; it compensates in part
the attenuation of the mollified sign and tends to 1 as ``\varepsilon \to 0``. Any positive
constant keeps the proposition true. The test `test/sliced/lift.jl` checks the limit to
``10 \varepsilon / L`` for ``\varepsilon = 10^{-9}``.

## [The canonical weight squeeze](@id sliced-squeeze-structure)

!!! info "Proposition (weight squeeze)"
    Let ``\mu \in (0, \infty)^N`` and ``x^{\mathrm{anc}} \in \mathbb{R}^N`` be constants and
    ``M = \mathrm{diag}(\mu) \big( I + \tfrac{1}{N} \mathbf{1} (\mathbf{1} - \mu)^\top \big)``.
    1. ``\det M = \prod_a \mu_a`` if ``\sum_a \mu_a = N``.
    2. If ``\sum_a \mu_a = N``, then
       ``M^{-\top} = \mathrm{diag}(\mu)^{-1} \big( I - \tfrac{1}{N} (\mathbf{1} - \mu) \mathbf{1}^\top \big)``.
    3. The map ``T_\mu : (x, p) \mapsto (\bar{x}, \bar{p})``,
       ``\bar{p} = M p``, ``\bar{x} - x^{\mathrm{anc}} = M^{-\top} (x - x^{\mathrm{anc}})``, is
       symplectic for every invertible ``M``. In components, for ``\sum_a \mu_a = N``,
       ```math
       \bar{p}_a = \mu_a \Big( p_a + \frac{1}{N} \sum_{c=1}^{N} (1 - \mu_c) \, p_c \Big) , \qquad
       \bar{x}_a = x^{\mathrm{anc}}_a + \frac{1}{\mu_a} \Big( (x_a - x^{\mathrm{anc}}_a)
         - (1 - \mu_a) \frac{1}{N} \sum_{c=1}^{N} (x_c - x^{\mathrm{anc}}_c) \Big) ,
       ```
       ``O(N)`` operations.
    4. ``\sum_a \bar{p}_a = \sum_a p_a`` for all ``p`` if and only if ``\sum_a \mu_a = N``.

!!! details "Proof"
    *1.* ``\det(I + u v^\top) = 1 + v^\top u``; with ``u = \mathbf{1}/N`` and
    ``v = \mathbf{1} - \mu``, ``v^\top u = (N - \sum_a \mu_a)/N = 0``.

    *2.* ``M^\top = (I + \tfrac{1}{N} (\mathbf{1} - \mu) \mathbf{1}^\top) \, \mathrm{diag}(\mu)``.
    By the Sherman–Morrison formula
    ``(I + u v^\top)^{-1} = I - u v^\top / (1 + v^\top u)`` with ``u = (\mathbf{1} - \mu)/N``,
    ``v = \mathbf{1}`` and ``v^\top u = 0``, so
    ``(M^\top)^{-1} = \mathrm{diag}(\mu)^{-1} (I - \tfrac{1}{N}(\mathbf{1} - \mu) \mathbf{1}^\top)``.

    *3.* The linear map ``(x, p) \mapsto (A x, B p)`` keeps
    ``\sum_a \mathrm{d}x_a \wedge \mathrm{d}p_a = \mathrm{d}x^\top \wedge \mathrm{d}p`` if and only
    if ``A^\top B = I``; here ``A = M^{-\top}``, ``B = M``. The constant anchors are a
    translation. The components are 2 applied to ``x - x^{\mathrm{anc}}`` and the definition of
    ``M``.

    *4.* ``\mathbf{1}^\top M = \mu^\top + \tfrac{1}{N} (\sum_a \mu_a) (\mathbf{1} - \mu)^\top``,
    which equals ``\mathbf{1}^\top`` if and only if ``\sum_a \mu_a = N``. For
    ``\sum_a \mu_a = N(1 + \delta)`` the total momentum changes by
    ``\delta \sum_c (1 - \mu_c) p_c``.

The squeeze scales each momentum by ``\mu_a``. With the final weight ``\bar{w}_a = w_{\mathrm{ref}} \mu_a``
the velocity ``\bar{v}_a = \bar{p}_a / (m \bar{w}_a) = (p_a + \text{const}) / (m w_{\mathrm{ref}})``
is the velocity of the slot before the squeeze, shifted by a common constant. The price is in the
positions: ``\bar{x}_a - x^{\mathrm{anc}}_a`` is ``(x_a - x^{\mathrm{anc}}_a)/\mu_a`` up to a common
term, so a slot with a small multiplier moves its marker far from its anchor. The test
`test/sliced/squeeze.jl` checks the closed forms against the dense ``M``, its determinant and the
momentum defect for ``\sum_a \mu_a \ne N``.

## [The composed map](@id sliced-composed-map)

!!! info "Theorem (sliced redistribution)"
    Freeze the slices, orders, targets and anchors of both lifts and the multipliers
    ``\mu > 0`` with ``\sum_a \mu_a = N``. Let ``\bar{W}(x, v, w) > 0`` be any weight rule. Then
    ``\mathcal{R} = \mathcal{F} \circ T_\mu \circ \Psi_p \circ \Psi_q`` with
    ``\mathcal{F} : \bar{w} = \bar{W}, \ \bar{v} = \bar{p} / (m \bar{w})`` satisfies
    ``\mathcal{R}^* \omega_N = \omega_N`` on ``(x, v, w)`` and keeps the total momentum. With
    ``\bar{W} = w_{\mathrm{ref}} \mu``, ``w_{\mathrm{ref}} = \frac{1}{N} \sum_a w_a``, it keeps the
    total weight. In ``d`` dimensions the same holds for the lifts and squeezes of each axis with
    one common ``\mu``.

!!! details "Proof"
    ``\Psi_q``, ``\Psi_p`` and ``T_\mu`` are symplectic diffeomorphisms of ``\mathbb{R}^{2N}`` that
    keep ``\sum_a p_a`` (the two propositions above), so their composition ``\psi`` is one. The
    output ``(\bar{x}, \bar{p}) = \psi(x, p)`` depends on the state only through ``(x, p)``, and
    ``\mathcal{F}`` is the split of the full-space form of the
    [condition](@ref symplecticity-conditions); so ``\mathcal{R}^* \omega_N = \omega_N``.
    The momentum is ``\sum_a m \bar{w}_a \bar{v}_a = \sum_a \bar{p}_a = \sum_a p_a``, and
    ``\sum_a w_{\mathrm{ref}} \mu_a = w_{\mathrm{ref}} N = \sum_a w_a``. In ``d`` dimensions the
    pairs ``(x_{a\ell}, p_{a\ell})`` of different axes are canonically independent, the lifts and
    squeezes of axis ``\ell`` act on the pairs of that axis only, and the weight rule is common.

The remap of the long-time experiment computes ``\mu`` from the state itself (the
[`squeeze_multipliers`](@ref) of a kernel density ratio). That map is the map of the theorem
evaluated at the multipliers of its input: it is canonical for the frozen ``\mu``, not as a map
of the state ([parameters that depend on the state](@ref state-dependent-parameters)). The
control of suite 3 with multipliers that depend on ``(\hat{x}, \hat{p})`` breaks the
symplecticity by ``10.8`` (`test/sliced/experiments.jl`).

## [Composition with a later flow step](@id sliced-composition-defect)

A redistribution is followed by time steps of the marker dynamics, whose Hamiltonian depends on
the weights, ``H(x, p; \bar{w}) = \sum_a |p_a|^2 / (2 m \bar{w}_a) + \sum_a \bar{w}_a \, q \, \phi(x_a)``,
and which keep the weights (the second of the [two readings](@ref two-readings)). Symplectic
integrators keep the canonical form of ``(x, p)`` at fixed parameters
[MarsdenWest:2001, HairerLubichWanner:2006](@cite); here the parameters come from the remap. The
proposition below has the form of the [closed-loop identity](@ref closed-loop-identity), with the
weights in the role of the parameters and the exact flow in place of one step.

!!! info "Proposition (composition defect)"
    Let ``\Phi^{\lambda}_t`` be the flow of ``H(\cdot\,; \lambda)`` at fixed parameters
    ``\lambda \in (0, \infty)^N``, let ``\psi`` be canonical, and let
    ``C(z) = \Phi^{\bar{W}(z)}_t(\psi(z))`` with ``\bar{W}(z) = W(\psi(z))`` the weights of the
    remap. Then
    ```math
    C^* \omega - \omega = -\sum_{a=1}^{N} \mathrm{d}A_a \wedge \mathrm{d}\bar{W}_a , \qquad
    A_a = \int_0^t \frac{\partial H}{\partial \bar{w}_a}\big(\Phi^{\bar{W}}_s(\psi(z)); \bar{W}\big) \, \mathrm{d}s .
    ```
    The defect vanishes if the weights do not depend on the state
    (``\mathrm{d}\bar{W} = 0``). On the extended space ``(x, p, \bar{w}, \theta)`` with the form
    ``\sum_a (\mathrm{d}x_a \wedge \mathrm{d}p_a + \mathrm{d}\bar{w}_a \wedge \mathrm{d}\theta_a)``
    the flow with ``\dot{\bar{w}}_a = 0`` and ``\dot{\theta}_a = -\partial H / \partial \bar{w}_a``
    is symplectic.

!!! details "Proof"
    On ``(x, p, \lambda, \theta)`` with ``\tilde{\omega} = \sum_a (\mathrm{d}x_a \wedge \mathrm{d}p_a
    + \mathrm{d}\lambda_a \wedge \mathrm{d}\theta_a)`` the Hamiltonian ``\tilde{H} = H(x, p; \lambda)``
    has the equations of ``H`` for ``(x, p)``, ``\dot{\lambda}_a = \partial \tilde{H} / \partial \theta_a = 0``
    and ``\dot{\theta}_a = -\partial \tilde{H} / \partial \lambda_a``. Its flow is
    ``\tilde{\Phi}_t(x, p, \lambda, \theta) = (G(x, p, \lambda), \lambda, \theta - A(x, p, \lambda))``
    with ``G(x, p, \lambda) = \Phi^\lambda_t(x, p)`` and ``A_a`` the integral of
    ``\partial H / \partial \lambda_a`` along the orbit. A Hamiltonian flow keeps
    ``\tilde{\omega}``:
    ``G^* \omega + \sum_a \mathrm{d}\lambda_a \wedge \mathrm{d}(\theta_a - A_a) =
    \omega + \sum_a \mathrm{d}\lambda_a \wedge \mathrm{d}\theta_a``, so
    ``G^* \omega - \omega = \sum_a \mathrm{d}\lambda_a \wedge \mathrm{d}A_a
    = -\sum_a \mathrm{d}A_a \wedge \mathrm{d}\lambda_a`` as forms on ``(x, p, \lambda)``. The
    composition is ``C = G \circ (\psi, \bar{W})``; pulling back by
    ``z \mapsto (\psi(z), \bar{W}(z))`` and using ``\psi^* \omega = \omega`` gives the formula.
    The last claim is the Hamiltonian character of the extended flow.

For the free flow ``\partial H / \partial \bar{w}_a = -|p_a|^2 / (2 m \bar{w}_a^2)`` is constant
along the orbit, and ``A_a = -t |\bar{p}_a|^2 / (2 m \bar{w}_a^2)``. Suite 4 of the
[experiments](@ref sliced-composition-suite) checks the formula, its vanishing for constant
weights, and the sign of ``\dot{\theta}``: with ``\dot{\theta}_a = +\partial H / \partial \bar{w}_a``
the extended flow is not symplectic.
