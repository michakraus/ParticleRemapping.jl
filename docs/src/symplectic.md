# Canonical structure of weighted markers

This chapter sets out the canonical structure of weighted particle-in-cell markers and
characterises every redistribution of the markers that preserves it. Its central statement is the
[factorisation lemma](@ref factorisation-lemma): a redistribution is symplectic if and only if it
is a canonical map of positions and momenta followed by a free choice of the new weights. The
[next chapter](@ref obstructions) asks what such maps can achieve, and the four constructions of
the manual are four ways to choose the canonical map and the weights.

```@meta
CurrentModule = ParticleRemapping
```

## [Weighted markers](@id weighted-markers)

Consider ``N`` markers, ``a = 1, \dots, N``, of one species with mass ``m > 0`` and charge
``q``. Marker ``a`` has a position ``x_a`` in the spatial domain (``\mathbb{R}^d`` or the torus
``\mathbb{T}^d``), a velocity ``v_a \in \mathbb{R}^d`` and a weight ``w_a > 0``. The markers
represent the density ``f_N = \sum_a w_a \, \delta(x - x_a) \, \delta(v - v_a)``, with the total
weight ``\sum_a w_a`` and the total momentum ``m \sum_a w_a v_a``.

In the package a set of markers is a [`Markers`](@ref) object: vectors `x`, `v`, `w` for
``d = 1``, or ``d \times N`` matrices `x`, `v` and a vector `w`. [`momentum`](@ref),
[`canonical`](@ref), [`from_canonical`](@ref), [`total_weight`](@ref),
[`total_momentum`](@ref) and [`kinetic_energy`](@ref) give the quantities of this chapter.

```@example symplectic
using ParticleRemapping
mk = Markers([0.3, 1.2, 2.0], [0.5, -1.0, 2.0], [0.2, 0.5, 0.3]; m = 2.0)
(p = momentum(mk), P = total_momentum(mk), kinetic = kinetic_energy(mk))
```

### Grid moments in one dimension

The constructions in one spatial and one velocity dimension (1D-1V) measure the markers on a
grid. The positions lie in the periodic interval ``[0, L)``, and the grid has ``N_g`` nodes
``x_i = (i - 1) \Delta x``, ``\Delta x = L / N_g``, with linear B-spline shape functions
``S_i``. On the cell ``[x_i, x_{i+1})`` the only non-zero shape functions are

```math
S_i(x) = 1 - \frac{x - x_i}{\Delta x} , \qquad S_{i+1}(x) = \frac{x - x_i}{\Delta x} ,
```

with the index ``i + 1`` taken modulo ``N_g``. They form a partition of unity,
``\sum_i S_i \equiv 1``. The grid moments of the marker set are

```math
\mathcal{M}_i^{(k)} = \sum_{a=1}^{N} w_a \, S_i(x_a) \, v_a^k ,
\qquad i = 1, \dots, N_g , \quad k = 0, 1, 2, \dots
```

The grid charge is ``\rho_i = q \mathcal{M}_i^{(0)}``, the grid current is
``j_i = q \mathcal{M}_i^{(1)}``, and ``\frac{m}{2} \mathcal{M}_i^{(2)}`` is the kinetic energy that
the shape function ``S_i`` collects. Because of the partition of unity,
``\sum_i \mathcal{M}_i^{(k)} = \sum_a w_a v_a^k``. The grid is a [`PeriodicGrid`](@ref), and the
moments are [`grid_moment`](@ref):

```@example symplectic
using Random
Random.seed!(1)
grid = PeriodicGrid(2π, 8)
mk = Markers(2π .* rand(200), randn(200), 0.5 .+ rand(200); m = 1.0)
# the partition of unity: the moments sum to the marker sums
sum(grid_moment(grid, mk, 2)) ≈ sum(mk.w .* mk.v .^ 2)
```

## [The Lagrangian and the canonical coordinates](@id canonical-coordinates)

Inserting ``f_N`` into the phase-space action of [Low:1958](@citet) gives the particle
Lagrangian of one species with an electrostatic potential ``\phi``
[BottinoSonnendrucker:2015, Kraus:2017gempic](@cite),

```math
L(x, \dot{x}, v; w) = \sum_{a=1}^{N} w_a \Big( m \, v_a \cdot \dot{x}_a
- \tfrac{1}{2} m |v_a|^2 - q \, \phi(x_a) \Big) .
```

Its Poincaré–Cartan one-form ``\vartheta = \sum_a (\partial L / \partial \dot{x}_a) \cdot
\mathrm{d}x_a`` identifies the canonical coordinates

```math
x_a , \qquad p_a = m \, w_a \, v_a , \qquad a = 1, \dots, N ,
```

and the canonical symplectic form

```math
\omega = -\mathrm{d}\vartheta = \sum_{a=1}^{N} \mathrm{d}x_a \wedge \mathrm{d}p_a
= \sum_{a=1}^{N} \sum_{i=1}^{d} \mathrm{d}x_a^i \wedge \mathrm{d}p_{a,i} .
```

The weight enters the canonical momentum, as it enters the Poisson structure of geometric
particle-in-cell methods [Kraus:2017gempic](@cite). The Hamiltonian is

```math
H(x, p; w) = \sum_{a=1}^{N} \frac{|p_a|^2}{2 m w_a} + U(x; w) ,
```

where ``U`` is the potential and field energy; it depends on the weights through the charge
density. In a 1D-1V particle-in-cell method with the grid potential
``\phi_h = \sum_i \phi_i S_i``, ``U`` depends on the markers only through the grid charge
``\rho``, with ``\partial U / \partial \rho_i = \phi_i``. Hamilton's equations are

```math
\dot{x}_a = \frac{\partial H}{\partial p_a} = \frac{p_a}{m w_a} = v_a , \qquad
\dot{p}_a = - \frac{\partial H}{\partial x_a} = q \, w_a \, E(x_a) ,
```

with ``E = -\nabla \phi``. Dividing the second equation by ``m w_a`` gives
``\dot{v}_a = (q/m) E(x_a)``: in ``(x, v)`` the motion does not depend on the weights. For fixed
weights the flow preserves ``\omega``.

!!! info "Definition (ensemble phase space)"
    The *single-marker phase space* is ``\mathcal{M}`` (``\mathbb{T}^d \times \mathbb{R}^d`` or
    ``\mathbb{R}^{2d}``) with points ``\zeta = (\xi, \eta)``. The *ensemble phase space* is
    ``\mathcal{Z} = \mathcal{M}^N`` with points ``z = (z_1, \dots, z_N)``,
    ``z_a = (x_a, p_a)``, and the form ``\omega``. The weights are parameters of the ensemble, not
    coordinates of ``\mathcal{Z}``. In matrix form, Hamilton's equations read
    ``\dot{z}_a = \mathbb{J} \nabla_{z_a} H`` with
    ``\mathbb{J} = \begin{pmatrix} 0 & I \\ -I & 0 \end{pmatrix}``
    ([`canonical_J`](@ref)).

The group ``\mathbb{R}^d`` acts on ``\mathcal{Z}`` by the translation ``x_a \mapsto x_a + c`` of
every position. Its momentum map is the total momentum

```math
P = \sum_{a=1}^{N} p_a = m \sum_{a=1}^{N} w_a v_a ,
```

which is the physical momentum of the represented density. By Noether's theorem, the flow of a
Hamiltonian that is invariant under the translation conserves ``P``. Requirement (R2) of the
[background](@ref requirements) asks a redistribution to keep it.

## [The presymplectic form on weighted markers](@id presymplectic-form)

A redistribution may change positions, velocities and weights. Its state space is the
``(2d + 1)N``-dimensional space of the states ``(x_a, v_a, w_a)``, ``w_a > 0``. Let

```math
\Pi(x, v, w) = (x, p) , \qquad p_a = m \, w_a \, v_a ,
```

be the *kinematic projection* onto ``\mathcal{Z}``. With Leibniz's rule,
``\mathrm{d}p_a = m ( w_a \, \mathrm{d}v_a + v_a \, \mathrm{d}w_a )``, so the canonical form,
written on the states, is the pull-back

```math
\omega_N = \Pi^* \omega
= \sum_{a=1}^{N} m \big( w_a \, \mathrm{d}x_a \wedge \mathrm{d}v_a
+ (v_a \cdot \mathrm{d}x_a) \wedge \mathrm{d}w_a \big) .
```

It is closed, because it is a pull-back of a closed form, and it is degenerate. A redistribution
``\mathcal{R} : (x, v, w) \mapsto (\bar{x}, \bar{v}, \bar{w})`` with ``\bar{w} > 0`` is
*symplectic* when ``\mathcal{R}^* \omega_N = \omega_N``.

The matrix of ``\omega_N`` for 1D-1V states ``(x_1, v_1, w_1, x_2, \dots)`` is
[`presymplectic_matrix`](@ref). It agrees with ``D\Pi^{\mathsf{T}} \mathbb{J} \, D\Pi``
([`kinematic_pullback`](@ref ParticleRemapping.Redistribution.kinematic_pullback)), exactly with
automatic differentiation and to round-off with fourth-order central differences of step
``10^{-5}`` (the map is quadratic, so the truncation error vanishes):

```@example symplectic
using LinearAlgebra
y0 = [0.4, -0.8, 0.25, 1.9, 1.3, 0.40]          # (x₁, v₁, w₁, x₂, v₂, w₂)
Ω = presymplectic_matrix(y0; m = 1.7)
Ω_fd = kinematic_pullback(y0; m = 1.7, jacobian = (F, y) -> fd_jacobian(F, y; order = 4))
(automatic = norm(kinematic_pullback(y0; m = 1.7) - Ω), differences = norm(Ω_fd - Ω),
 rank = rank(Ω))
```

### [The kernel: the gauge direction](@id kernel-gauge)

To find the kernel of ``\omega_N``, insert a tangent vector
``\xi = \sum_a ( \xi^x_a \cdot \nabla_{x_a} + \xi^v_a \cdot \nabla_{v_a} + \xi^w_a \partial_{w_a} )``:

```math
\iota_\xi \omega_N = \sum_a m \Big( \xi^x_a \cdot ( w_a \, \mathrm{d}v_a + v_a \, \mathrm{d}w_a )
- ( w_a \xi^v_a + v_a \xi^w_a ) \cdot \mathrm{d}x_a \Big) .
```

This vanishes if and only if ``\xi^x_a = 0`` and ``w_a \xi^v_a + v_a \xi^w_a = 0`` for every
``a``. Hence ``\ker \omega_N = \ker \mathrm{d}\Pi`` has dimension ``N`` and is spanned by the
Euler scaling fields

```math
\xi_a = w_a \frac{\partial}{\partial w_a} - v_a \cdot \nabla_{v_a} ,
\qquad \xi_a(p_b) = m \, \delta_{ab} \, ( w_a v_a - v_a w_a ) = 0 .
```

The flow of ``\xi_a`` is ``\dot{w}_a = w_a``, ``\dot{v}_a = - v_a``, that is
``(v_a, w_a) \mapsto (e^{-s} v_a, e^{s} w_a)``. With ``\lambda_a = e^{-s}`` it is the rescaling
``(v_a, w_a) \mapsto (\lambda_a v_a, w_a / \lambda_a)`` at a fixed momentum. In the
``(v_a, w_a)`` plane of a 1D-1V marker it moves the marker along the hyperbola
``w_a v_a = p_a / m``. The fibres
``\Pi^{-1}(x, p) = \{ (x, p / (m w), w) : w \in \mathbb{R}^N_{>0} \}`` are connected.
[`euler_scaling_fields`](@ref ParticleRemapping.Redistribution.euler_scaling_fields) returns the
fields as the columns of a matrix, and ``\omega_N`` annihilates them:

```@example symplectic
norm(Ω * euler_scaling_fields(y0))
```

```@example symplectic
using CairoMakie
plot_gauge_orbits()
```

The figure shows (a) the hyperbolae ``w v = p/m`` along which the gauge direction moves a marker,
with the images of the marker ``(v, w) = (1, 1)`` under the stretches ``\lambda = 1/2, 1, 2``, and
(b) the ``(x, v)`` cell of one marker under the same stretches: its area scales by ``\lambda``
while its weight scales by ``1/\lambda``.

## [The factorisation lemma](@id factorisation-lemma)

!!! info "Lemma (factorisation)"
    Let ``\mathcal{R}`` be a smooth map of the states, and let ``\mathrm{d}(\Pi \circ \mathcal{R})``
    have rank ``2Nd`` everywhere. Then ``\mathcal{R}^* \omega_N = \omega_N`` if and only if there
    is a canonical map ``\psi`` of ``(\mathcal{Z}, \omega)`` with
    ``\Pi \circ \mathcal{R} = \psi \circ \Pi``. The new weights are not restricted:
    ``\bar{w} = \bar{W}(x, v, w)`` is any positive function, and the new velocities are
    ``\bar{v}_a = \bar{p}_a / (m \bar{w}_a)`` with ``(\bar{x}, \bar{p}) = \psi(x, p)``.

!!! details "Proof"
    If ``\Pi \circ \mathcal{R} = \psi \circ \Pi`` with ``\psi^* \omega = \omega``, then
    ``\mathcal{R}^* \omega_N = \mathcal{R}^* \Pi^* \omega = (\Pi \circ \mathcal{R})^* \omega
    = (\psi \circ \Pi)^* \omega = \Pi^* \psi^* \omega = \Pi^* \omega = \omega_N``.

    Conversely, let ``\mathcal{R}^* \omega_N = \omega_N``, let ``\xi \in \ker \mathrm{d}\Pi(z)``
    and let ``\eta`` be any tangent vector at ``z``. Then
    ``\omega \big( \mathrm{d}(\Pi \circ \mathcal{R}) \xi, \mathrm{d}(\Pi \circ \mathcal{R}) \eta \big)
    = (\mathcal{R}^* \omega_N)(\xi, \eta) = \omega_N(\xi, \eta) = 0``. Since
    ``\mathrm{d}(\Pi \circ \mathcal{R}) \eta`` ranges over the whole tangent space of
    ``\mathcal{Z}`` and ``\omega`` is non-degenerate, ``\mathrm{d}(\Pi \circ \mathcal{R}) \xi = 0``.
    Hence ``\Pi \circ \mathcal{R}`` is constant on the connected fibres of ``\Pi`` and factors as
    ``\psi \circ \Pi``. The map ``\psi = \Pi \circ \mathcal{R} \circ \sigma`` with the smooth section
    ``\sigma(x, p) = (x, p / m, 1)`` of ``\Pi`` is smooth. From
    ``\Pi^* \psi^* \omega = \omega_N = \Pi^* \omega`` and the injectivity of ``\Pi^*`` on forms
    (``\Pi`` is a submersion), ``\psi^* \omega = \omega``. The weights enter only through the fibre
    coordinate, so they are free.

Every symplectic redistribution is therefore a canonical map ``\psi`` of ``(x, p)``, followed by
the choice of the new weights at the new momenta. The four constructions of the manual are four
choices of ``\psi`` and of ``\bar{W}`` ([the constructions](@ref constructions)).

## [Three forms of the symplecticity condition](@id symplecticity-conditions)

The lemma is the condition on the full state space. Two further forms are useful: one for a
redistribution that is applied at fixed input weights, and one for a redistribution generated by
a flow in a pseudo-time ``\tau``.

!!! info "Proposition (three forms of the condition)"
    1. **Full space.** ``\mathcal{R}^* \omega_N = \omega_N`` on all states if and only if
       ``(\bar{x}, \bar{p}) = \psi(x, p)`` with a canonical map ``\psi`` that depends on
       ``(x, v, w)`` only through ``(x, p)``, and ``\bar{w}_a = \bar{W}_a(x, v, w) > 0``,
       ``\bar{v}_a = \bar{p}_a / (m \bar{w}_a)`` with any positive weight rule ``\bar{W}``.
    2. **Leaf of constant weights.** For fixed input weights ``w``, the restriction of
       ``\mathcal{R}`` to the leaf ``L_w = \{(x, v, w)\}`` preserves
       ``\omega_N|_{L_w} = m \sum_a w_a \, \mathrm{d}x_a \wedge \mathrm{d}v_a`` if and only if
       ``(\bar{x}, \bar{p}) = \psi_w(x, p)`` with a canonical ``\psi_w`` that may depend on ``w``
       as a parameter, and ``\bar{v}_a = \bar{p}_a / (m \bar{w}_a)``.
    3. **Pseudo-time flow.** The flow of a vector field on the states preserves ``\omega_N`` if
       and only if there are a generator ``\mathcal{G}_\tau(x, p)`` and weight rates ``\gamma_a``
       (any functions of the state) with
       ```math
       \frac{\mathrm{d}x_a}{\mathrm{d}\tau} = \frac{\partial \mathcal{G}_\tau}{\partial p_a} , \qquad
       \frac{\mathrm{d}w_a}{\mathrm{d}\tau} = \gamma_a , \qquad
       \frac{\mathrm{d}v_a}{\mathrm{d}\tau} = -\frac{1}{m w_a} \frac{\partial \mathcal{G}_\tau}{\partial x_a}
                                              - \frac{\gamma_a}{w_a} \, v_a .
       ```
       The last term is the recoil: it changes the velocity with the weight so that
       ``\mathrm{d}p_a / \mathrm{d}\tau = -\partial \mathcal{G}_\tau / \partial x_a``.

!!! details "Proof"
    *1.* This is the factorisation lemma.

    *2.* On the leaf, ``\Pi_w = \Pi|_{L_w} : (x, v) \mapsto (x, m w v)`` is a diffeomorphism onto
    ``\mathcal{Z}`` because ``w > 0``, and ``\omega_N|_{L_w} = \Pi_w^* \omega``. Define
    ``\psi_w = (\Pi \circ \mathcal{R})|_{L_w} \circ \Pi_w^{-1}``. Then
    ``(\mathcal{R}|_{L_w})^* \omega_N = \Pi_w^* \psi_w^* \omega``, and this equals
    ``\Pi_w^* \omega`` if and only if ``\psi_w^* \omega = \omega``, since ``\Pi_w`` is a
    diffeomorphism. Nothing restricts the dependence of ``\psi_w`` on ``w``.

    *3.* Let ``X`` be the vector field and use the coordinates ``(x, p, w)``. Since ``\omega_N``
    has no ``\mathrm{d}w`` terms in these coordinates,
    ``\iota_X \omega_N = \sum_a (\dot{x}_a \cdot \mathrm{d}p_a - \dot{p}_a \cdot \mathrm{d}x_a)``,
    with ``\dot{x}_a = \mathrm{d}x_a(X)`` and ``\dot{p}_a = \mathrm{d}p_a(X)`` functions of
    ``(x, p, w)``. The flow preserves ``\omega_N`` if and only if
    ``\mathcal{L}_X \omega_N = \mathrm{d} \, \iota_X \omega_N = 0``. The parts of
    ``\mathrm{d} \, \iota_X \omega_N`` along ``\mathrm{d}w_b \wedge \mathrm{d}p_a`` and
    ``\mathrm{d}w_b \wedge \mathrm{d}x_a`` are ``\partial_{w_b} \dot{x}_a`` and
    ``-\partial_{w_b} \dot{p}_a``, so ``\dot{x}``, ``\dot{p}`` do not depend on ``w``. The
    remaining one-form on ``\mathcal{Z}`` is closed, hence exact,
    ``\iota_X \omega_N = \mathrm{d}\mathcal{G}_\tau``, which is
    ``\dot{x}_a = \partial \mathcal{G}_\tau / \partial p_a``,
    ``\dot{p}_a = -\partial \mathcal{G}_\tau / \partial x_a``. The rates
    ``\gamma_a = \mathrm{d}w_a(X)`` do not enter. Finally ``v_a = p_a / (m w_a)`` gives
    ``\dot{v}_a = \dot{p}_a / (m w_a) - p_a \gamma_a / (m w_a^2)``, the stated equation. The
    converse is the same computation read backwards.

The second form is the one that the constructions use when they compute ``\psi`` from the
markers at fixed weights; the [state-dependent parameters](@ref state-dependent-parameters) of the
next chapter say what this means for the map on the full space. The suites 1 and 2 of the
[sliced experiments](@ref sliced-experiments) check all three forms with ForwardDiff Jacobians and
controls that drop the recoil or read ``w`` at fixed ``(x, p)``.

## [The gauge step](@id gauge-step)

Every symplectic redistribution therefore splits into two steps: a canonical map ``\psi`` of the
ensemble in ``(x, p)``, and a *gauge step*

```math
G_\lambda : (x_a, v_a, w_a) \mapsto (x_a, \lambda_a v_a, w_a / \lambda_a) ,
\qquad \lambda_a > 0 ,
```

which moves along the kernel of ``\omega_N``. The gauge step is the identity on ``(x, p)``:
``\bar{p}_a = m (w_a / \lambda_a) (\lambda_a v_a) = p_a`` as functions on the states, whatever
``\lambda`` depends on. Hence ``\Pi \circ G_\lambda = \Pi`` and ``G_\lambda^* \omega_N = \omega_N``
for every positive ``\lambda``, *also when ``\lambda`` depends on the state of the whole
ensemble*. It preserves the total momentum ``P``. The gauge step is also the only part that is
not volume-preserving in ``(x, v)``: the Jacobian of ``(x_a, v_a) \mapsto (x_a, \lambda_a v_a)``
at fixed ``\lambda_a`` is ``\lambda_a^d``, so a marker's ``(x, v)``-volume scales by
``\lambda_a^d`` while its weight scales by ``1/\lambda_a``. This is the freedom that the weights
add.

Read in the canonical coordinates, the gauge step is *a change of weight at fixed* ``(x, p)``:
any new weights ``\bar{w}_a`` with the velocities

```math
\bar{v}_a = \frac{p_a}{m \bar{w}_a} = \frac{w_a}{\bar{w}_a} \, v_a ,
```

which changes the kinetic energy from ``\sum_a |p_a|^2 / (2 m w_a)`` to
``\sum_a |p_a|^2 / (2 m \bar{w}_a)``. All the symplectic content of a redistribution is
therefore in the motion of ``(x_a, p_a)``; the weights act as labels of the markers, as the
charge-to-mass ratio does for several species, and the physical Hamiltonian depends on these
labels. A change of weight at fixed *velocity* is different: it maps ``p_a`` to
``(\bar{w}_a / w_a) p_a``, and if ``\bar{w}_a`` depends on the state its Jacobian is not
symplectic in general.

The common layer implements the gauge step as [`gauge_step`](@ref). Its pull-back residual
([`pullback_residual`](@ref), central differences with step ``10^{-5}``) vanishes to the
truncation error also for a stretch that depends on the state of the marker:

```@example symplectic
λ(z) = 1 + z[1]^2 + sin(z[2]) * z[3]^2 + 0.3z[3]
gauge(z) = [z[1], λ(z) * z[2], z[3] / λ(z)]
pullback_residual(gauge, [0.3, -1.2, 0.8]; m = 1.7)
```

## [Kicks, drifts and cotangent lifts](@id elementary-canonical-maps)

The following elementary canonical maps of the ensemble are exactly symplectic even when their
data depend on the ensemble, because each is the time-one flow of a Hamiltonian on
``\mathcal{Z}``: a *kick* ``p_a \mapsto p_a + \nabla_{x_a} \chi`` with
``\chi = \chi(x_1, \dots, x_N)``; a *drift* ``x_a \mapsto x_a + \nabla_{p_a} \eta`` with
``\eta = \eta(p_1, \dots, p_N)``; and their compositions. A *cotangent lift*
``\bar{x}_a = \varphi(x_a)``, ``\bar{p}_a = D\varphi(x_a)^{-\mathsf{T}} p_a`` of a diffeomorphism
``\varphi`` of the spatial domain is canonical for a fixed ``\varphi``; in one dimension,

```math
\mathrm{d}\bar{x}_a \wedge \mathrm{d}\bar{p}_a
= \varphi' \, \mathrm{d}x_a \wedge \Big( \frac{\mathrm{d}p_a}{\varphi'}
- \frac{p_a \varphi''}{\varphi'^2} \mathrm{d}x_a \Big) = \mathrm{d}x_a \wedge \mathrm{d}p_a .
```

The cotangent lift of the gradient of a convex potential in the momenta, with the roles of ``x``
and ``p`` exchanged, is canonical in the same way; the [sliced redistribution](@ref sliced-overview)
uses both. A kick, then a drift, followed by an arbitrary new weight, preserves ``\omega_N``, as
the lemma states; so does a cotangent lift with the weight ``\bar{w} = w/\varphi'`` and
``\bar{v} = v``:

```@example symplectic
m = 1.7
function kick_drift_reweight(z)
    x, v, w = z
    P = m * w * v - sin(x)                 # kick with χ = cos x
    X = x + P^3 / 5                        # drift with η = p⁴/20
    W = 0.5 + x^2 + w * cos(v)^2           # any positive weight
    [X, P / (m * W), W]
end
φ(x) = x + 0.3sin(x)
lift(z) = [φ(z[1]), z[2], z[3] / (1 + 0.3cos(z[1]))]
(pullback_residual(kick_drift_reweight, [1.1, 0.4, 2.5]; m),
    pullback_residual(lift, [1.1, 0.4, 2.5]; m))
```
