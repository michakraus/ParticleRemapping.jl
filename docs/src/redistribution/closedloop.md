# [Closed loops](@id rd-closed-loop)

```@meta
CurrentModule = ParticleRemapping.Redistribution
```

Once ``\theta^*`` is computed from the current ensemble ``z = (x, p)``, the open-loop map
``z' \mapsto \Psi_\tau^{\theta^*}(z')`` with ``\theta^*`` held fixed is symplectic by the
[linear parameter solver](@ref rd-linear-parameter-solver) (iii). The closed-loop map
``\Phi(z) := \Psi_\tau^{\Theta(z)}(z)`` also moves the parameters ``\theta = \Theta(z)`` with the
incoming ensemble. The [closed-loop identity](@ref closed-loop-identity) measures its defect: with
``g_k(z) := \tfrac{\partial H}{\partial \theta_k}(z_{\mathrm{stage}}(z); \Theta(z))`` at the stage
point of the implicit midpoint or the symplectic Euler step,

```math
(D\Phi)^{\mathsf{T}} \mathbb{J} (D\Phi) - \mathbb{J}
= -\tau \Big( (Dg)^{\mathsf{T}} (D\Theta) - (D\Theta)^{\mathsf{T}} (Dg) \Big) .
```

Its corollary states that a rule that resets the positions onto a fixed grid is never symplectic,
and that the closed loop is symplectic when ``\Theta`` satisfies the
[variational closure](@ref variational-closure)
``\partial_\theta H(z_{\mathrm{stage}}; \theta) = \partial_\theta \mathcal{V}(\theta)``. This
page adds the tools and the cases of this chapter: [`closed_loop_defect`](@ref) evaluates both
sides of the identity from the Jacobians ``(D\Phi, D\Theta, Dg)``, and
[`fixed_grid_closed_loop`](@ref) is the rule of the fixed-grid no-go: its Jacobian has zero
position rows, and ``(D\Phi)^{\mathsf{T}} \mathbb{J} (D\Phi) = 0``. The remark below constructs a
variational closure that spreads spatial clusters, and the last section evaluates the identity on
two closed loops.

### [Remark (spatial de-clustering by closed-loop cotangent lifts)](@id rd-cotangent-lift)

!!! note "Remark (spatial de-clustering by closed-loop cotangent lifts)"
    An autonomous potential Hamiltonian ``H_{\mathrm{eff}}(z) = \operatorname{stat}_\theta
    [H(z; \theta) - \mathcal{V}(\theta)]`` conserves ``H_{\mathrm{eff}}`` along its flow
    (``\{H_{\mathrm{eff}}, H_{\mathrm{eff}}\} = 0``), so pure functions of ``x`` alone cannot drive
    a gradient descent in ``x``. To spread spatial clusters while satisfying the variational
    closure, couple ``\theta`` to translation-invariant cotangent-lift generators
    ```math
    H_{\mathrm{var}}(x, p; \theta) = -\frac{1}{N} \sum_{m=1}^{M} \theta_m \sum_{a,b=1}^{N}
    (p_a - p_b) \cdot \nabla\psi_m(x_a - x_b) , \qquad
    \mathcal{V}(\theta) = \frac{\kappa}{2} \| \theta - \theta_0 \|_2^2 ,
    ```
    where ``\psi_m`` are repulsive kernels and ``\theta_{0,m} > 0``. At leading order
    ``\dot{x}_a \approx -\tfrac{2}{N} \sum_{m,b} \theta_{0,m} \nabla\psi_m(x_a - x_b)`` drives a
    repulsive gradient flow on ``x``, while the conjugate update
    ``\dot{p}_a = -\nabla_{x_a} H_{\mathrm{var}}`` and the stationarity closure
    ``\theta_m = \theta_{0,m} + \kappa^{-1} \partial_{\theta_m} H_{\mathrm{var}}(z_{\mathrm{stage}})``
    preserve both ``\sum_a p_a`` and the closed-loop symplecticity ``\Phi^*\omega = \omega``. In
    contrast, contracting an open ``2 N d``-dimensional neighbourhood of ``(x, v)`` states
    towards a lower-dimensional manifold of uniform layouts changes the ``(x, p)``-volume when
    ``\bar{w}_a`` is eliminated as a function of ``(\bar{x}_a, \bar{v}_a)``, so full
    ``(x, v)`` target matching necessarily operates in the open-loop sense of the [linear parameter
    solver](@ref rd-linear-parameter-solver) (iii).

With the stationarity closure, the closed-loop step is the implicit midpoint step of the
autonomous Hamiltonian ``H_{\mathrm{eff}}(z) = \theta_0 \cdot \Psi(z) + |\Psi(z)|^2 / (2\kappa)``,
since ``\nabla_z H_{\mathrm{eff}} = \sum_m (\theta_{0,m} + \Psi_m / \kappa) \nabla_z \Psi_m``;
this is why the closure is symplectic.

## The two closed loops of the verification

The verification uses four markers, ``z = (0.5, 0.9, 2.8, 4.7, 0.04, -0.06, 0.03, -0.01)``, and
a pairwise spline Hamiltonian with ``p_{\mathrm{scale}} = 0.25`` and quartic bases of
``5 \times 4`` cells ([`closed_loop_example`](@ref)). SimpleSplines 0.4 rejects a periodic quartic
basis of ``4`` cells in ``x`` (it needs more cells than the degree); ``5`` is the smallest
admissible number.

1. *Control E, target matching* ([`target_matching_closed_loop`](@ref)): the targets repel the
   positions, ``\bar{x}_a = x_a + 0.04 \sum_{b \ne a} \Delta_{ab} e^{-\Delta_{ab}^2/2}``, keep
   ``\bar{p} = p``, and ``\Theta(z)`` is the Tikhonov-regularised (``\lambda = 10^{-2}``)
   solution of the midpoint parameter system ([`target_matching_parameters`](@ref)). The map is
   one midpoint step with ``\Theta(z)``.
2. *Variational closure* ([`variational_closed_loop`](@ref)): one midpoint step of
   ``H_{\mathrm{var}}`` with the three modes ``\nabla\psi_k(\Delta x) = \sin(k \Delta x)``
   ([`cotangent_lift_modes`](@ref)), ``\theta_0 = (0.06, 0.03, 0.015)`` and ``\kappa = 2``,
   solved by Newton's method jointly with the closure.

Both return exact Jacobians by implicit differentiation of the step equations
([`implicit_step_jacobian`](@ref)); fourth-order differences give a second check.

```@example rd-closed-loop
using ParticleRemapping, CairoMakie, LinearAlgebra
ex = closed_loop_example()
cl = target_matching_closed_loop(ex.H, ex.z)
lhs, rhs = closed_loop_defect(cl.jacobians...)
vl = variational_closed_loop(ex.z)
DΦ, DΘ = vl.jacobians
J = canonical_J(4)
(control_E_defect = norm(lhs), identity_error = norm(lhs - rhs) / norm(lhs),
 variational_defect = norm(transpose(DΦ) * J * DΦ - J), parameter_jacobian = norm(DΘ))
```

```@example rd-closed-loop
plot_closed_loop_defect(lhs, rhs; variational = transpose(DΦ) * J * DΦ - J)
```

The figure shows the ``8 \times 8`` matrices of the closed-loop defect of Control E (a) and of
the right-hand side of the identity (b), the defect of the variational closure on the same signed
colour scale (c, zero to round-off), and the difference of the first two on its own scale (d, of
order ``10^{-16}``).

The script `scripts/redistribution/verify_symplectic_redistribution.jl` gives a defect of
``0.396`` for Control E (the control requires more than ``10^{-3}``), a relative error of the
identity of ``3.6 \cdot 10^{-15}`` with the exact Jacobians and ``6.0 \cdot 10^{-10}`` with
fourth-order differences (the bound is ``10^{-8}``), and for the symplectic Euler version of
Control E a defect of ``0.395`` and a relative error of ``1.7 \cdot 10^{-15}``. The variational
closure has ``\|D\Theta\| = 1.04`` and a closed-loop defect of ``6.5 \cdot 10^{-16}`` with the
exact Jacobian; with fourth-order differences the defect is ``5.3 \cdot 10^{-10}``, the round-off
of the differences. [The verification](@ref rd-verification) compares these numbers with the
bounds of the sources.
