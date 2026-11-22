# [Symplectic discretisation and parameter solver](@id rd-solver)

```@meta
CurrentModule = ParticleRemapping.Redistribution
```

## Linearity of the inverse parameter problem

Write the linear parametrisation of the ensemble Hamiltonian compactly as

```math
H_\theta(z) = \sum_{k=1}^{N_\theta} \theta_k \Psi_k(z) , \qquad
z = (x, p) \in \mathbb{R}^{2 N d} , \qquad
\mathbb{J} = \begin{pmatrix} 0 & I_{N d} \\ -I_{N d} & 0 \end{pmatrix} .
```

The flow of ``H_\theta`` in artificial time over one fixed step ``\tau > 0`` (normalised to
``\tau = 1``) is approximated by the second-order symmetric implicit midpoint method or by the
first-order Type-2 symplectic Euler method [HairerLubichWanner:2006](@cite):

```math
\begin{aligned}
&\text{implicit midpoint:} &
\bar{x}_a - x_a &= \tau \nabla_{p_a} H_\theta\Big(\frac{x + \bar{x}}{2}, \frac{p + \bar{p}}{2}\Big) , &
\bar{p}_a - p_a &= -\tau \nabla_{x_a} H_\theta\Big(\frac{x + \bar{x}}{2}, \frac{p + \bar{p}}{2}\Big) , \\
&\text{symplectic Euler:} &
\bar{x}_a - x_a &= \tau \nabla_{\bar{p}_a} H_\theta(x, \bar{p}) , &
\bar{p}_a - p_a &= -\tau \nabla_{x_a} H_\theta(x, \bar{p}) .
\end{aligned}
```

The two schemes are [`ImplicitMidpointStep`](@ref) and [`SymplecticEulerStep`](@ref).

### [Proposition (linear parameter solver and forward branch invertibility)](@id rd-linear-parameter-solver)

!!! info "Proposition (linear parameter solver and forward branch invertibility)"
    Let ``z = (x, p)`` and ``\bar{z} = (\bar{x}, \bar{p})`` satisfy the discrete momentum
    compatibility ``\sum_a \bar{p}_a = \sum_a p_a`` that the [KKT
    projection](@ref rd-kkt-projection) enforces.

    1. *Linearity in ``\theta``.* In both schemes the stage point —
       ``z_{1/2} = \tfrac12 (z + \bar{z})`` for implicit midpoint, ``z_{\mathrm{SE}} = (x, \bar{p})``
       for symplectic Euler — is known from ``(z, \bar{z})``. The parameters
       ``\theta \in \mathbb{R}^{N_\theta}`` therefore solve the linear system
       ```math
       G(z_{\mathrm{stage}}) \, \theta = \bar{z} - z , \qquad
       G_{:,k}(z_{\mathrm{stage}}) := \tau \mathbb{J} \nabla_z \Psi_k(z_{\mathrm{stage}}) \in \mathcal{M}_0 ,
       ```
       where ``\mathcal{M}_0 := \{ (\delta x, \delta p) \in \mathbb{R}^{2 N d} :
       \sum_a \delta p_a = 0 \}`` has dimension ``2 N d - d``.
    2. *Dimensionless conditioning (``N_\theta \ge 2 N d - d``).* Set
       ``S = \operatorname{diag}(I_{N d}, p_{\mathrm{scale}}^{-1} I_{N d})``. If the stage
       points ``z_{\mathrm{stage},a}`` are pairwise distinct and
       ``\operatorname{rank} G(z_{\mathrm{stage}}) = 2 N d - d``, the minimum-norm solution
       ```math
       \theta^* = (S G)^+ S (\bar{z} - z)
       ```
       satisfies the linear system identically and makes ``\bar{z}`` an exact fixed point of the
       implicit midpoint or the symplectic Euler step.
    3. *Contraction condition for the forward map.* For implicit midpoint, if
       ```math
       \sup_{z' \in \mathcal{B}} \frac{\tau}{2}
       \big\| S \mathbb{J} \nabla_z^2 H_{\theta^*}(z') S^{-1} \big\| < 1
       ```
       on a ball ``\mathcal{B}`` that contains ``z`` and ``\bar{z}``, then the forward nonlinear
       map ``\Psi_\tau^{\theta^*}(z)`` is a strict contraction around ``z`` whose unique root in
       ``\mathcal{B}`` is ``\bar{z}``, with the nonsingular symplectic Jacobian
       ```math
       D_z \Psi_\tau^{\theta^*}(z) = \Big( I - \tfrac{\tau}{2} \mathbb{J} \nabla_z^2 H_{\theta^*}(z_{1/2})
       \Big)^{-1} \Big( I + \tfrac{\tau}{2} \mathbb{J} \nabla_z^2 H_{\theta^*}(z_{1/2}) \Big) .
       ```
       For symplectic Euler, the step is locally invertible on the identity-connected branch when
       ``\| \tau \nabla_x \nabla_{\bar{p}} H_{\theta^*}(x, \bar{p}) \| < 1``.

!!! details "Proof"
    Substituting the linear parametrisation into the two schemes gives the linear system at
    once, and Noether's identity ``\sum_a \nabla_{x_a} H_\theta = 0`` gives
    ``\sum_a G_{p,a,k} = 0``, so ``\operatorname{span} G \subseteq \mathcal{M}_0``. Since
    ``\bar{z} - z \in \mathcal{M}_0`` by the discrete momentum compatibility, full rank
    ``2 N d - d`` implies ``\bar{z} - z \in \operatorname{span} G``.

    For (iii), write implicit midpoint in the scaled coordinates ``\hat{z} = S z`` as
    ``\hat{z}_{\mathrm{out}} = T(\hat{z}_{\mathrm{out}}) := \hat{z}_{\mathrm{in}} + \tau S \mathbb{J}
    \nabla_z H_{\theta^*}\big(S^{-1} \tfrac{\hat{z}_{\mathrm{in}} + \hat{z}_{\mathrm{out}}}{2}\big)``.
    The Fréchet derivative of ``T`` with respect to ``\hat{z}_{\mathrm{out}}`` is
    ``\tfrac{\tau}{2} S \mathbb{J} \nabla_z^2 H_{\theta^*} S^{-1}``. When the contraction condition
    holds, Banach's fixed-point theorem gives the uniqueness of the root ``\bar{z}`` in
    ``\mathcal{B}``, and implicit differentiation gives the Jacobian, a Cayley transform of the
    Hamiltonian matrix ``\tfrac{\tau}{2} \mathbb{J} \nabla_z^2 H_{\theta^*}``, which satisfies
    ``(D_z \Psi)^{\mathsf{T}} \mathbb{J} (D_z \Psi) = \mathbb{J}``.

The parameter solve and the forward step of the package:

1. Evaluate the stage point of the scheme from ``(z, \bar{z})``
   ([`redistribution_stage_point`](@ref)).
2. Assemble ``G`` from the gradients of the basis Hamiltonians at the stage point
   ([`redistribution_parameter_matrix`](@ref)).
3. Solve ``(S G) \theta = S (\bar{z} - z)`` for the minimum-norm ``\theta^*``
   ([`solve_redistribution_parameters`](@ref), with [`redistribution_scaling`](@ref) for ``S``).
4. Advance ``z`` with ``\theta^*`` fixed by Newton's method, following the identity-connected
   branch by the homotopy ``\theta \to s \theta^*``, ``s = 0.25, 0.5, 0.75, 1`` (midpoint) or
   ``s = 0.2, 0.4, \dots, 1`` (symplectic Euler), each to a scaled residual below ``10^{-14}``
   ([`forward_redistribution_step`](@ref)).
5. Return the Jacobian of the step by the implicit function theorem: the Cayley transform for the
   midpoint, and for symplectic Euler
   ``\partial \bar{p} / \partial (x, p) = (I + H_{x\bar{p}})^{-1} (-H_{xx}, I)`` and
   ``\partial \bar{x} / \partial (x, p) = (I + H_{\bar{p}x}, 0) + H_{\bar{p}\bar{p}} \,
   \partial \bar{p} / \partial (x, p)``.

[`redistribution_contraction_number`](@ref) evaluates the left-hand side of the contraction
condition at the stage point. On the 16 markers of the verification:

```@example rd-solver
using ParticleRemapping, CairoMakie, LinearAlgebra
ex = redistribution_example()
mk = ex.markers
q, p = canonical(mk)
N = length(q)
res = remap(SplineRedistribution(), mk)              # projection, targets, KKT weights, one step
q̄, v̄ = res.info.targets_x, res.info.targets_v
p̄ = res.info.target_weights .* mk.m .* v̄
H = res.info.hamiltonian
S = redistribution_scaling(N, H.p_scale)
G = redistribution_parameter_matrix(H, ImplicitMidpointStep(), q, p, q̄, p̄)
θ = solve_redistribution_parameters(H, ImplicitMidpointStep(), q, p, q̄, p̄)
mid = forward_redistribution_step(H, ImplicitMidpointStep(), q, p, θ)
J = canonical_J(N)
(rank = rank(Diagonal(S) * G), parameters = length(θ),
 linear_residual = norm(S .* (G * θ - vcat(q̄ .- q, p̄ .- p))),
 target_error = norm(S .* (mid.z - vcat(q̄, p̄))),
 symplecticity = norm(transpose(mid.jacobian) * J * mid.jacobian - J),
 momentum_change = sum(mid.z[(N + 1):end]) - sum(p),
 contraction = redistribution_contraction_number(H, ImplicitMidpointStep(), q, p, q̄, p̄, θ))
```

The rank is ``2 N d - d = 31``, and the midpoint step reaches the targets. In parentheses below
are the numbers for a second relaxation of eight steps of size ``0.02``, which
[the verification](@ref rd-verification) compares with the numbers of the sources. The script
`scripts/redistribution/verify_symplectic_redistribution.jl` gives a target error
``9.9 \cdot 10^{-15}`` (``8.04 \cdot 10^{-15}``), a symplecticity defect ``5.5 \cdot 10^{-15}``
by the Cayley transform and ``1.5 \cdot 10^{-8}`` by fourth-order differences, and a relative
difference ``5.8 \cdot 10^{-10}`` between the two Jacobians.

The contraction number at the stage point is ``2.27`` in the verification (``2.01``): the
condition of (iii) is violated, and the midpoint homotopy still reaches ``\bar{z}``. The
condition is sufficient, not necessary. For symplectic Euler the quantity
``\| \nabla_x \nabla_{\bar{p}} H_{\theta^*}(x, \bar{p}) \|`` is ``1.39`` (``1.47``), and the
identity-connected branch ends on another root of the same equations: the target error is
``0.184`` (``0.150``), a caustic. ``\bar{z}`` itself solves the symplectic Euler equations to
round-off; the homotopy does not reach it, also with 1000 homotopy stages (the script checks
both).

[`redistribution_substeps`](@ref) composes ``K`` steps towards the intermediate targets
``z + (k/K)(\bar{z} - z)``, each momentum-compatible, with one parameter solve per substep. Each
substep is symplectic at its fixed parameters, and so is the composition. With two substeps the
symplectic Euler bound holds (``0.78``) and the composition reaches ``\bar{z}``. With four
substeps the midpoint bound holds too (``0.80``, and ``0.99`` for the second relaxation);
`test/redistribution/solver.jl` checks that three suffice in the verification set-up:

```@example rd-solver
se1 = redistribution_substeps(H, SymplecticEulerStep(), q, p, q̄, p̄)
se2 = redistribution_substeps(H, SymplecticEulerStep(), q, p, q̄, p̄; substeps = 2)
z̄ = vcat(q̄, p̄)
(one_step = (contraction = se1.contraction, target_error = norm(S .* (se1.z - z̄))),
 two_substeps = (contraction = se2.contraction, target_error = norm(S .* (se2.z - z̄))))
```

```@example rd-solver
plot_newton_convergence(mid.history, only(se1.history), vcat(se2.history...);
    labels = ["midpoint, 1 step", "symplectic Euler, 1 step", "symplectic Euler, 2 substeps"])
```

The figure shows the scaled residual of every Newton iterate of the forward solves; the dashed
lines mark the first iterate of each homotopy stage. Every stage converges quadratically to
round-off. The single symplectic Euler step converges as fast as the others, to the wrong root.

The animation shows the midpoint step in pseudo-time: the end points
``\Psi_{s\tau}^{\theta^*}(z)`` of the step with the artificial step ``s \tau``, ``s \in [0, 1]``
([`pseudotime_path`](@ref)), drawn in ``(x, v)`` at the old weights, ``v = p / (m w)``, and
coloured and sized by weight. At ``s = 1`` the markers reach the targets (open circles) in
``(x, p)``; the gauge step then assigns the KKT weights at fixed ``(x, p)`` and moves them onto
the targets in ``(x, v)`` ([`remap_frames`](@ref ParticleRemapping.remap_frames)).

```@example rd-solver
animate_stages("redistribution-step.mp4", remap_frames(SplineRedistribution(), mk, res);
    targets = (q̄, v̄), period = (0, ex.Lx))
nothing # hide
```

```@raw html
<video src="./redistribution-step.mp4" controls loop muted playsinline style="width:100%"></video>
```

### [Remark (caustics and the coarse-basis regime)](@id rd-coarse-basis)

!!! note "Remark (caustics and the coarse-basis regime ``N_\theta \ll 2 N d - d``)"
    Two practical regimes arise when two clustered markers sit within
    ``|z_a - z_b| = \varepsilon \ll 1``:

    1. Since ``\nabla_{z_a} H_\theta = \nabla_{z_b} H_\theta`` when ``z_a = z_b``, separating an
       ``\varepsilon``-close pair by ``O(1)`` in a single step forces
       ``\| \nabla^2 H_{\theta^*} \| \gtrsim \varepsilon^{-1}``, which violates the contraction
       condition and creates spurious root branches (caustics). Keeping the relaxation
       displacement of a single step proportional to ``\varepsilon`` (or composing several
       symplectic redistribution substeps) maintains the contraction condition.
    2. For large ``N`` with a coarse collective basis (``N_\theta \ll 2 N d - d``), solve the
       Tikhonov-regularised least-squares problem
       ``\min_\theta \| S (G \theta - (\bar{z} - z)) \|_2^2 + \lambda \| \theta \|_2^2``, advance
       ``(\hat{x}, \hat{p}) = \Psi_\tau^{\theta^*}(x, p)`` forward by implicit midpoint or
       symplectic Euler, and recover ``\hat{v}_a`` near ``\bar{v}_a`` along the continuous branch
       of ``m \Delta\hat{\Omega}_a f_h(\hat{x}_a, \hat{v}_a) \hat{v}_a = \hat{p}_a``, with
       ``\hat{w}_a = \Delta\hat{\Omega}_a f_h(\hat{x}_a, \hat{v}_a)``.

The Tikhonov solve is [`solve_redistribution_parameters`](@ref) with `λ > 0`; the branch recovery
is [`recover_velocity_branch`](@ref), Newton's method in ``v`` from ``\bar{v}_a``. The
recovery follows the branch of the [velocity–momentum fold](@ref velocity-momentum-fold) on which
``\bar{v}_a`` lies. The verification runs the regime (a) — the single symplectic Euler step and the
substeps above — and checks the branch recovery from perturbed starts in
`test/redistribution/projection.jl`.

```@example rd-solver
θλ = solve_redistribution_parameters(H, ImplicitMidpointStep(), q, p, q̄, p̄; λ = 1e-6)
coarse = forward_redistribution_step(H, ImplicitMidpointStep(), q, p, θλ)
q̂, p̂ = coarse.z[1:N], coarse.z[(N + 1):end]
ps = PhaseSpaceSplines()
# the volumes of the targets stand in for those of the end points
rec = recover_velocity_branch(ps, res.info.coefficients, q̂, p̂, v̄, res.info.volumes; m = mk.m)
(target_error = norm(S .* (coarse.z - z̄)), velocity_offset = maximum(abs, rec.v - v̄),
 momentum_change = sum(rec.w .* mk.m .* rec.v) - sum(p))
```
