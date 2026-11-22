# [Canonical transport with consistent weights](@id remap-transport)

```@meta
CurrentModule = ParticleRemapping.Remapping
```

The canonical transport takes a given canonical map ``\psi`` of ``(x, p)``, the same for all
markers, and then chooses the weights so that the markers represent the same density ``f``. The
analysis uses the continuum limit, in which the markers form a density in phase space.

## Transport

The continuum description is that of
[Liouville's theorem for the marker density](@ref liouville-holes): the markers have a smooth
density ``g(x, v) > 0`` and a smooth weight field ``w(x, v) > 0``, they represent ``f = g \, w``,
and the sheet map ``\sigma(x, v) = (x, m \, w(x, v) \, v)`` carries them to the canonical marker
density

```math
n = \frac{g}{m \, \partial_v (w v)} \circ \sigma^{-1} ,
```

under the single-sheet condition ``\partial_v (w v) > 0`` of the
[velocity–momentum fold](@ref velocity-momentum-fold). This chapter treats one layer of a folded
sheet. The canonical map ``\psi`` transports the canonical density to
``\tilde{n} = n \circ \psi^{-1}``.

After the transport every marker receives a new weight. Write the new weight as a field
``W(X, V)`` of the new coordinates, with ``V = P / (m W)`` and the single-sheet condition
``\partial_V (W V) > 0``. The same counting argument as for ``n``, this time from ``(X, V)`` to
``(X, P)`` with ``\mathrm{d} P = m \, \partial_V (W V) \, \mathrm{d} V``, gives the new marker
density in ``(X, V)`` and the density that the markers represent:

```math
g' = \tilde{n}(X, m W V) \, m \, \partial_V (W V) , \qquad f' = g' \, W .
```

This is [`represented_density`](@ref).

## Consistency of the weights

!!! info "Proposition (current-quantile matching)"
    Define, at fixed ``X`` and on each sign branch,
    ```math
    F_X(V) = \int_0^V f(X, V') \, V' \, \mathrm{d} V' , \qquad
    \Phi_X(P) = \frac{1}{m} \int_0^P \tilde{n}(X, P') \, P' \, \mathrm{d} P' .
    ```
    Then ``f' = f`` if and only if, for every ``X`` and ``V``,
    ```math
    F_X(V) = \Phi_X\big( m \, W(X, V) \, V \big) .
    ```

!!! details "Proof"
    Fix ``X`` and substitute ``P = m W V``. Then
    ``\mathrm{d} P = m \, \partial_V (W V) \, \mathrm{d} V`` and ``W V = P / m``. By the
    expression for ``f'``,
    ```math
    f' \, V \, \mathrm{d} V = \tilde{n} \, m \, \partial_V(W V) \, W \, V \, \mathrm{d} V
    = \tilde{n} \cdot \frac{P}{m} \cdot m \, \partial_V(W V) \, \mathrm{d} V
    = \tilde{n} \, \frac{P}{m} \, \mathrm{d} P .
    ```
    At ``V = 0`` also ``P = 0``. Integrating from ``0`` gives
    ``\int_0^V f' V' \mathrm{d} V' = \Phi_X(m W V)`` for every ``V``. Hence the matching holds if
    and only if ``\int_0^V (f' - f) V' \mathrm{d} V' = 0`` for all ``V``. Differentiating in ``V``
    gives ``(f' - f) V = 0``, that is ``f' = f`` for ``V \neq 0``, and by continuity everywhere.
    The converse direction is the same computation read backwards.

The proof uses the exact Jacobian ``\partial_V(W V)``. The cruder relation
``g' = \tilde{n} \, m \, W``, which treats the weight of each marker as frozen, is wrong as soon as
``W`` depends on ``V``. The identity of the proof holds to a relative error of
``1.8 \times 10^{-14}`` (Simpson's rule) for a canonical density and a weight field that depends on
``V``, and the frozen-weight density fails it by ``5 \times 10^{-2}``
(`test/remapping/transport.jl`):

```@example transport
using ParticleRemapping
m = 1.3
ñ(X, P) = (1 + 0.3sin(X)) * exp(-P^2 / 2) * (1 + 0.2P^2)
W(X, V) = 1 + 0.3cos(X) + 0.2tanh(V)
dWV(X, V) = W(X, V) + 0.2V * sech(V)^2                  # ∂_V (W V)
pts = [(X, V) for X in (0.4, 2.1) for V in (-2.5, -0.7, 0.9, 3.0)]
exact = represented_density(ñ, W, dWV; m)
frozen = represented_density(ñ, W, dWV; m, frozen = true)
(exact = maximum(cumulative_current_error(exact, ñ, W, X, V; m) for (X, V) in pts),
    frozen = minimum(cumulative_current_error(frozen, ñ, W, X, V; m) for (X, V) in pts))
```

On each branch, ``F_X`` is strictly increasing in ``|V|`` where ``f > 0``, because
``\partial_V F_X = f V`` has the sign of ``V``, and ``\Phi_X`` is strictly increasing in ``|P|``
where ``\tilde{n} > 0``. The matching therefore defines the new velocity by a
*current-quantile matching*,

```math
V(X, P) = F_X^{-1}\big( \Phi_X(P) \big) , \qquad W = \frac{P}{m \, V} ,
```

on each branch. The marker whose canonical momentum lies at a given quantile of the cumulative
momentum gets the velocity at the same quantile of the cumulative current of ``f``.

!!! info "Corollary (compatibility)"
    A weight field with ``f' = f`` exists if and only if, for every ``X``,
    ```math
    J_\pm^\psi(X) := \frac{1}{m} \int_{\pm P > 0} P \, \tilde{n}(X, P) \, \mathrm{d} P
    = \int_{\pm V > 0} V f(X, V) \, \mathrm{d} V =: J_\pm^f(X) .
    ```
    For ``\psi = \mathrm{id}`` both sides agree. The condition says that ``\psi`` preserves the
    sign-resolved current profile of the markers.

!!! details "Proof"
    On the branch ``V > 0``, ``F_X`` maps ``[0, \infty)`` onto ``[0, J_+^f(X))`` and ``\Phi_X``
    maps ``[0, \infty)`` onto ``[0, J_+^\psi(X))``. The matching needs every value of ``\Phi_X``
    to be a value of ``F_X`` and conversely, since every velocity must be reached; this is
    ``J_+^\psi = J_+^f``. The same holds on the branch ``V < 0``. Conversely, under the
    compatibility condition the map ``P \mapsto V`` of the matching is a strictly increasing
    bijection on each branch with ``V P > 0``, so ``W = P / (m V) > 0``. The derivative is
    ``\partial_P V = \tilde{n} P / (m f V)``, so
    ``\partial_V (W V) = \partial_V P / m = f V / (\tilde{n} P) > 0``, and the single-sheet
    condition holds.

The physical meaning is direct. The gauge step cannot change any momentum, and in a consistent
representation the momentum density at ``X`` is ``m`` times the current of ``f``; the canonical
map must therefore deliver the correct current profile on each branch.

!!! info "Corollary (uniqueness)"
    For a given ``\psi``, the weight field is unique. For ``\psi = \mathrm{id}`` it is ``W = w``,
    and ``V = v``. All the freedom of the construction is in ``\psi``.

!!! details "Proof"
    Uniqueness follows from the strict monotonicity of ``F_X`` in the matching. For
    ``\psi = \mathrm{id}``, ``\tilde{n} = n``. Substituting ``p' = m w v'``, with
    ``\mathrm{d} p' = m \, \partial_{v'}(w v') \, \mathrm{d} v'``, in ``\Phi_X`` and using the
    expression for ``n``,
    ```math
    \Phi_X(m w v) = \frac{1}{m} \int_0^v \frac{g}{m \, \partial_{v'}(w v')} \, m w v' \,
    m \, \partial_{v'}(w v') \, \mathrm{d} v' = \int_0^v g \, w \, v' \, \mathrm{d} v' = F_X(v) .
    ```
    So ``V = v`` solves the matching, and ``W = p / (m v) = w``.

!!! info "Corollary (the limit of equidistribution)"
    The new marker density is ``g' = f / W = m f V / P``. A uniform density
    ``g' \equiv \bar{g}`` requires ``W = f / \bar{g}`` and ``P = m f V / \bar{g}``. This needs
    ``\partial_V (f V) > 0`` on the support, and the transported canonical density must be
    ``\tilde{n}(X, m f V / \bar{g}) = \bar{g}^2 / \big( m \, \partial_V (f V) \big)``. Since
    ``\tilde{n} = n \circ \psi^{-1}``, a canonical ``\psi`` reaches it only if ``n`` and this
    target are equimeasurable: the areas ``|\{ n > c \}|`` and ``|\{ \tilde{n} > c \}|`` agree
    for every ``c``.

!!! details "Proof"
    ``g' = f / W`` is the expression for ``g'`` with ``f' = f``. Setting ``g' = \bar{g}`` gives
    ``W = f/\bar{g}`` and ``P = m W V = m f V / \bar{g}``, and the single-sheet condition
    ``\partial_V(WV) > 0`` becomes ``\partial_V(fV) > 0``. Solving the expression for ``g'`` for
    ``\tilde{n}`` gives
    ``\tilde{n} = g' / (m \partial_V(WV)) = \bar{g}^2 / (m \partial_V(fV))``. A map that preserves
    area preserves the area of every level set of the density it transports, which is the
    equimeasurability.

### [The limit](@id remap-limit)

For a Maxwellian, ``\partial_V (f V) < 0`` for ``|V| > v_{\mathrm{th}}``, so a uniform marker
distribution with consistent weights does not even form a single sheet. In general the
equimeasurability fails as well. A smooth canonical map with consistent weights can therefore
rearrange the canonical marker density, but it cannot make the marker distribution uniform. In
particular it cannot fill a region in which ``n`` is small: it can only move that region, for
example into the tails of ``f``, where few markers are needed. This is the obstruction
[holes persist](@ref holes-persist) of [one map for all markers](@ref one-map-for-all-markers),
stated for consistent weights.

## Admissible canonical maps

The compatibility condition restricts ``\psi``. For a map generated by a Hamiltonian ``h(x, p)``,
the restriction is linear in ``h``.

!!! info "Proposition (admissible generators)"
    Let ``\psi_t`` be the flow of ``h(x, p)`` with ``\partial_x h(x, 0) = 0``, so that the line
    ``p = 0`` is invariant and no marker changes its branch. Then
    ```math
    \frac{\mathrm{d}}{\mathrm{d} t} J_\pm^{\psi_t} \Big|_{t = 0}
    = - \frac{1}{m} \bigg[ \partial_x \int_{\pm p > 0} p \, n \, \partial_p h \, \mathrm{d} p
    + \int_{\pm p > 0} n \, \partial_x h \, \mathrm{d} p \bigg] .
    ```
    The flow keeps the compatibility to first order if and only if the bracket vanishes for both
    branches; along a finite flow the same condition applies to the transported density
    ``\tilde{n}_t``.

!!! details "Proof"
    For a test function ``\varphi(x)``, the definition of ``J_+`` and the substitution
    ``(X, P) = \psi_t(x, p)``, which preserves area, give
    ```math
    \int \varphi \, J_+^{\psi_t} \mathrm{d} X = \frac{1}{m} \iint_{P > 0} \varphi(X) \, P \,
    \tilde{n}_t \, \mathrm{d} X \, \mathrm{d} P = \frac{1}{m} \iint_{p > 0} \varphi(x_t) \, p_t \,
    n(x, p) \, \mathrm{d} x \, \mathrm{d} p ,
    ```
    where ``(x_t, p_t) = \psi_t(x, p)``. By Hamilton's equations ``\dot{x}_t = \partial_p h`` and
    ``\dot{p}_t = - \partial_x h``, so at ``t = 0`` the derivative of the integrand is
    ``\varphi'(x) \, \partial_p h \, p - \varphi(x) \, \partial_x h``. Integration by parts in
    ``x`` moves the derivative from ``\varphi`` to ``\int p \, n \, \partial_p h \, \mathrm{d} p``
    and gives the formula. The invariance of ``p = 0`` keeps the domain ``p > 0`` fixed, so there
    is no boundary term.

!!! info "Example (cotangent lift)"
    Let ``h = p \, a(x)``, whose flow is ``X = \varphi_t(x)``, ``P = p / \varphi_t'(x)`` with
    ``\dot{\varphi}_t = a(\varphi_t)``. Here ``\partial_p h = a`` and ``\partial_x h = p a'``, so
    the bracket becomes ``\partial_x ( a \, m J_\pm ) + a' \, m J_\pm``, using
    ``\frac{1}{m} \int p \, n \, \mathrm{d} p = J_\pm``. Up to the factor ``m`` this is
    ``(a J_\pm)' + a' J_\pm = (a^2 J_\pm)' / a``. The lift is admissible if and only if
    ``a^2 J_+`` and ``a^2 J_-`` are constant, that is ``a = c \, J_+^{-1/2}``, which requires
    ``J_+ / |J_-|`` to be constant, as for an ``f`` that is even in ``v``. The flow then keeps the
    profile exactly, ``J_\pm(\varphi_t(x)) \, \varphi_t'(x)^2 = J_\pm(x)`` (the
    [closed-form example](@ref remap-lift-example) derives this). This map moves markers in
    ``x``, faster where the current is small. It is not designed to equidistribute.

The cotangent lift is [`cotangent_lift`](@ref), its flow [`cotangent_lift_flow`](@ref) (the
classical Runge–Kutta method on ``(x, \varphi')``). Cotangent lifts are among the
[elementary canonical maps](@ref elementary-canonical-maps); the sliced redistribution composes
[sliced quantile cotangent lifts](@ref sliced-lift-structure).

The set of admissible generators is a linear space. A practical remap can choose ``h`` from a
finite basis, impose the linear constraints of the proposition for both branches, and pick the
coefficients that decrease a measure of marker clustering.

## [Worked example: a separable density in closed form](@id remap-lift-example)

!!! info "Example (cotangent lift with matching)"
    Let ``f(x, v) = J(x) \, M(v)`` with ``M(v) = e^{-v^2/2} / \sqrt{2\pi}``,
    ``J(x) = 1 + 0.4 \cos x`` on ``[0, 2\pi)`` ([`cosine_profile`](@ref)), and markers with equal
    weights ``w_0``, so ``g = f / w_0`` and ``p = m w_0 v``. Write
    ```math
    G(V) = \int_0^{|V|} v M(v) \, \mathrm{d} v = \frac{1 - e^{-V^2/2}}{\sqrt{2\pi}} , \qquad
    G^{-1}(y) = \sqrt{-2 \ln \big( 1 - \sqrt{2\pi} \, y \big)} ,
    ```
    with ``G(\infty) = 1 / \sqrt{2\pi}`` ([`maxwellian_current`](@ref),
    [`maxwellian_current_inverse`](@ref)).

**Step 1: the densities.** The cumulative current of ``f`` on the positive branch is
``F_X(V) = J(X) \, G(V)``, so ``J_\pm^f = J(X) / \sqrt{2\pi}``, and ``J_+ / |J_-| = 1``. By the
expression for ``n`` with ``\partial_v(w_0 v) = w_0``, the canonical marker density is
``n(x, p) = J(x) \, M(p / (m w_0)) / (m w_0^2)``.

**Step 2: the canonical map.** Take the cotangent lift of ``\dot{x} = a(x)``,
``X = \varphi_t(x)``, ``P = p / \varphi_t'(x)``. Its inverse is ``x = \varphi_t^{-1}(X)``,
``p = P \, \varphi_t'(x)``, so by ``\tilde{n} = n \circ \psi^{-1}``

```math
\tilde{n}(X, P) = \frac{J(x) \, M\big( P \varphi_t'(x) / (m w_0) \big)}{m w_0^2} ,
\qquad x = \varphi_t^{-1}(X) .
```

**Step 3: the cumulative momentum.** Substituting ``u = P' \varphi_t' / (m w_0)``, so that
``P' = m w_0 u / \varphi_t'`` and ``\mathrm{d} P' = m w_0 \, \mathrm{d} u / \varphi_t'``,

```math
\Phi_X(P) = \frac{1}{m} \int_0^P \frac{J(x) \, M(u)}{m w_0^2} \, P' \, \mathrm{d} P'
= \frac{J(x)}{\varphi_t'(x)^2} \, G\Big( \frac{P \varphi_t'(x)}{m w_0} \Big) .
```

**Step 4: compatibility.** As ``P \to \infty``,
``J_+^\psi(X) = J(x) / (\sqrt{2\pi} \, \varphi_t'^2)``, and the compatibility condition becomes
``J(\varphi_t(x)) \, \varphi_t'(x)^2 = J(x)``. For a general ``a`` the relative compatibility
residual ([`lift_compatibility_residual`](@ref)) is

```math
R(x) = \frac{J(x)}{\varphi_t'(x)^2 \, J(\varphi_t(x))} - 1 .
```

**Step 5: the admissible flow in closed form.** For ``a = J^{-1/2}``, separate the variables in
``\dot{\varphi} = J(\varphi)^{-1/2}``: with ``\Theta(x) = \int_0^x J(y)^{1/2} \, \mathrm{d} y``
([`lift_theta`](@ref)),

```math
\Theta(\varphi_t(x)) = \Theta(x) + t .
```

Differentiating in ``x`` gives ``J(\varphi_t)^{1/2} \varphi_t' = J(x)^{1/2}``, that is
``\varphi_t'(x) = \big( J(x) / J(\varphi_t(x)) \big)^{1/2}`` and ``R \equiv 0``
([`admissible_lift`](@ref)).

**Step 6: the matching.** With ``R = 0``, ``J(x) / \varphi_t'^2 = J(X)``, and the matching reads
``J(X) \, G(V) = J(X) \, G(P \varphi_t' / (m w_0))``. Hence

```math
V = \frac{P \, \varphi_t'(x)}{m w_0} = \frac{p}{m w_0} = v , \qquad
W = \frac{P}{m V} = \frac{w_0}{\varphi_t'(x)} .
```

The remap moves every marker in ``x`` at a fixed velocity and scales its weight by the inverse
stretch of the map ([`admissible_lift_matching`](@ref)). It is exactly the family
``X = \varphi(x)``, ``V = c(x) v + b(x)``, ``W = w / (\varphi' c)`` with ``c = 1``, ``b = 0``,
which [`symplectic_condition_checks`](@ref) confirms to be canonical.

**Step 7: what it does to the markers.** Before the remap the spatial marker density is
proportional to ``J``. After it, it is proportional to
``J(x) / \varphi_t'(x) = \big( J(x) J(X) \big)^{1/2}``, the geometric mean of ``J`` at the
departure and the arrival point, while ``\sum W`` per cell still represents ``J(X)``. The ratio of
the largest to the smallest marker density is ``7/3 = 2.33`` at ``t = 0``. At half a period of
``\Theta``, ``t = \Theta(2\pi)/2 = 3.1089``, every marker travels from near a maximum of ``J`` to
near a minimum or conversely, and the ratio falls to ``1.18``.

**Numbers at ``t = 1.5``.**

| ``x``                                        | 0       | ``\pi/2`` | ``\pi``  | 4        |
|:---------------------------------------------|:--------|:----------|:---------|:---------|
| ``X = \varphi_{1.5}(x)``                     | 1.3197  | 3.3328    | 4.8586   | 5.4954   |
| ``\varphi_{1.5}'(x)``                        | 1.1285  | 1.2832    | 0.7530   | 0.7590   |
| ``W / w_0 = 1/\varphi_{1.5}'``               | 0.8862  | 0.7793    | 1.3281   | 1.3176   |
| ``R(x)`` for ``a = 1 + 0.5\cos x``           | 3.70    | 4.03      | −0.53    | −0.85    |

For ``a = J^{-1/2}``, ``J(X) \varphi'^2 / J(x) = 1`` holds to ``3 \times 10^{-16}`` and the
matched velocity equals ``v`` to ``2 \times 10^{-15}``. The last row shows that the lift of
``a = 1 + 0.5 \cos x``, which is not admissible, violates the compatibility by several hundred per
cent. The numbers are checked in `test/remapping/transport_examples.jl`, which also compares the
closed form with the Runge–Kutta flow of ``\dot{x} = J^{-1/2}``, and printed by
`scripts/remapping/worked_examples.jl`:

```@example transport
xs = [0.0, π / 2, π, 4.0]
lifts = [admissible_lift(x, 1.5) for x in xs]
J(x) = cosine_profile(x)
R = [lift_compatibility_residual(J, x,
         cotangent_lift_flow(x -> 1 + 0.5cos(x), x -> -0.5sin(x), x, 1.5)...) for x in xs]
(X = first.(lifts), dφ = last.(lifts), W_over_w0 = 1 ./ last.(lifts), R_inadmissible = R,
    V_minus_v = maximum(abs(admissible_lift_matching(x, v, 1.5, 0.01).V - v)
                        for x in xs, v in (-1.3, 0.4, 2.2)),
    half_period = lift_theta(2π) / 2)
```

## The numerical example and the limits

The [additional experiments](@ref remap-transport-example) apply the
[algorithm](@ref remap-matching) to ``2 \times 10^5``
markers sampled from the density of the [closed-form example](@ref remap-lift-example), after the
identity, the admissible lift and an inadmissible lift. What the canonical transport cannot do,
and why, is collected under [limitations](@ref remapping-limitations), and the matching along
rays in ``d`` velocity dimensions under [higher dimensions](@ref remapping-higher-dimensions).
