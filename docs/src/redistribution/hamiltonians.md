# [Translation- and permutation-invariant ensemble Hamiltonians](@id rd-hamiltonians)

```@meta
CurrentModule = ParticleRemapping.Redistribution
```

## Failure of single-particle Hamiltonians and ad-hoc force projection

Before ``H_\theta(x, p)`` is specified, two natural simplifications are ruled out. An uncoupled
map moves a hole but does not shrink it ([uncoupled maps](@ref uncoupled-maps)). Part (i) of the
lemma adds a second obstruction for uncoupled Hamiltonians: if they keep the total momentum, they
change no momentum at all.

### [Lemma (obstructions for uncoupled Hamiltonians and mean-force subtraction)](@id rd-obstructions)

!!! info "Lemma (obstructions for uncoupled Hamiltonians and mean-force subtraction)"
    Let ``N \ge 2``.

    1. Let ``H(x, p) = \sum_{a=1}^{N} h(x_a, p_a)`` be an uncoupled sum of ``C^1``
       single-particle Hamiltonians. Then ``H(x_1 + \xi, \dots, x_{N} + \xi, p) = H(x, p)``
       for all ``(x, p)`` and all ``\xi \in \mathbb{R}^d`` if and only if ``\nabla_x h \equiv 0``,
       and then ``\dot{p}_a = 0`` for every marker.
    2. Let ``h \in C^2(\Omega_x \times \mathbb{R}^d)`` with ``\nabla_p \nabla_x h(x_0, p_0) \ne 0``
       at some point ``(x_0, p_0)``. Then the zero-mean-force vector field
       ```math
       \dot{x}_a = \nabla_p h(x_a, p_a) , \qquad
       \dot{p}_a = -\nabla_x h(x_a, p_a) + \frac{1}{N} \sum_{c=1}^{N} \nabla_x h(x_c, p_c)
       ```
       conserves ``\sum_a p_a``, but ``\mathcal{L}_X \omega \ne 0``: it is neither Hamiltonian
       nor symplectic.

!!! details "Proof"
    (i) Differentiating ``H(x + \xi, p) = H(x, p)`` at ``\xi = 0`` gives
    ``\sum_a \nabla_x h(x_a, p_a) = 0`` for all ``(x, p)``. Setting ``(x_a, p_a) = (y, u)`` for all
    ``a`` gives ``N \nabla_x h(y, u) = 0``, hence ``\nabla_x h \equiv 0`` (Control A).

    (ii) A vector field ``(\dot{x}, \dot{p})`` on ``\mathbb{R}^{2 N d}`` satisfies
    ``\mathcal{L}_X \omega = 0`` only if ``\partial \dot{x}_b / \partial x_a =
    -(\partial \dot{p}_a / \partial p_b)^{\mathsf{T}}`` for all ``a, b``. In the field above, for
    ``a \ne b``, ``\partial \dot{x}_b / \partial x_a = 0``, whereas
    ``-\partial \dot{p}_a / \partial p_b = -\tfrac{1}{N} \nabla_p \nabla_x h(x_b,
    p_b)^{\mathsf{T}} \ne 0`` at ``(x_b, p_b) = (x_0, p_0)`` (Control B).

The two controls of the verification are single steps on the packed state ``z = (x, p)`` of the
16 clustered markers: [`uncoupled_control_step`](@ref) (Control A) changes the total momentum, and
[`mean_force_control_step`](@ref) (Control B) keeps it but is not symplectic.

```@example rd-hamiltonians
using ParticleRemapping, CairoMakie, ForwardDiff, LinearAlgebra
ex = redistribution_example()
q, p = canonical(ex.markers)
N = length(q)
p_scale = 0.6986828002504696           # the momentum scale of the verification (solver page)
z = vcat(q, p)
D = ForwardDiff.jacobian(y -> mean_force_control_step(y, p_scale), z)
J = canonical_J(N)
(control_A_momentum_change = sum(uncoupled_control_step(z, p_scale)[(N + 1):end]) - sum(p),
 control_B_momentum_change = sum(mean_force_control_step(z, p_scale)[(N + 1):end]) - sum(p),
 control_B_symplecticity_defect = norm(transpose(D) * J * D - J))
```

The script `scripts/redistribution/verify_symplectic_redistribution.jl` gives ``0.157`` for the
momentum change of Control A and ``0.0914`` for the symplecticity defect of Control B; both
controls require more than ``10^{-2}``.

## Pairwise and collective symmetric Hamiltonians

The lemma shows that a momentum-conserving ensemble Hamiltonian that can change individual
momenta ``p_a`` must couple the markers through the relative positions ``x_a - x_b`` and treat all
markers symmetrically. Since ``\sum_a w_a = \mathcal{M}`` is fixed independently of ``N``, the
individual weights and momenta scale as ``w_a = O(N^{-1})`` and ``p_a = O(N^{-1})``, whereas
``x_a = O(1)``. A characteristic momentum scale ``p_{\mathrm{scale}} \sim \max_a |p_a| =
O(N^{-1})`` makes all basis functions and derivatives dimensionless and of order one.

### [Proposition (symmetric pairwise and collective Hamiltonians)](@id rd-symmetric-hamiltonians)

!!! info "Proposition (symmetric pairwise and collective Hamiltonians)"
    Let ``\{\chi_r : [-1, 1]^d \to \mathbb{R}\}_{r=1}^{M_p}`` be a tensor-product B-spline basis in
    the scaled momentum ``\hat{p} = p / p_{\mathrm{scale}}``, and
    ``\{\phi_m : \Omega_x \to \mathbb{R}\}_{m=1}^{M_x}`` a smooth periodic basis on
    ``\Omega_x = \mathbb{T}^d`` (such as periodic B-splines). Define the parametrised ensemble
    Hamiltonian
    ```math
    \begin{aligned}
    H_\theta(x, p) ={}& p_{\mathrm{scale}} \sum_{a=1}^{N} \sum_{r=1}^{M_p} \theta^{(0)}_r
    \chi_r\Big(\frac{p_a}{p_{\mathrm{scale}}}\Big) \\
    &+ \frac{p_{\mathrm{scale}}}{2 N} \sum_{a,b=1}^{N} \sum_{m=1}^{M_x} \sum_{l,r=1}^{M_p}
    \theta^{(2)}_{m,l,r} \, \phi_m(x_a - x_b) \, \chi_l\Big(\frac{p_a}{p_{\mathrm{scale}}}\Big)
    \chi_r\Big(\frac{p_b}{p_{\mathrm{scale}}}\Big) ,
    \end{aligned}
    ```
    with the parameter vector ``\theta = (\theta^{(0)}, \theta^{(2)}) \in \mathbb{R}^{N_\theta}``,
    ``N_\theta = M_p + M_x M_p^2``. Then:

    1. *Permutation invariance:* ``H_\theta(\sigma x, \sigma p) = H_\theta(x, p)`` for every
       permutation ``\sigma \in S_{N}`` and every ``\theta``.
    2. *Continuous translation invariance:* ``H_\theta(x_1 + \xi, \dots, x_{N} + \xi, p) =
       H_\theta(x, p)`` for every ``\xi \in \mathbb{R}^d``, so Noether's identity
       ```math
       \sum_{a=1}^{N} \nabla_{x_a} H_\theta(x, p) = 0
       ```
       holds identically on ``(\mathbb{T}^d)^{N} \times \mathbb{R}^{N d}``.
    3. *Collective factorisation in ``O(N)``:* when ``\{\phi_m(\Delta x)\}`` is the Fourier basis
       ``\{\cos(k \cdot \Delta x), \sin(k \cdot \Delta x)\}_{k \in \mathcal{K}}``, the double sum
       factors through the collective moments
       ```math
       C_{k,r}(x, p) = \frac{1}{N} \sum_{b=1}^{N} \cos(k \cdot x_b) \,
       \chi_r\Big(\frac{p_b}{p_{\mathrm{scale}}}\Big) , \qquad
       S_{k,r}(x, p) = \frac{1}{N} \sum_{b=1}^{N} \sin(k \cdot x_b) \,
       \chi_r\Big(\frac{p_b}{p_{\mathrm{scale}}}\Big) ,
       ```
       and ``H_\theta``, ``\nabla_x H_\theta`` and ``\nabla_p H_\theta`` cost
       ``O(N |\mathcal{K}| M_p^2)`` operations instead of ``O(N^2)``.
    4. *Centroid-frame B-spline alternative on ``\mathbb{R}^d``:* on ``\Omega_x = \mathbb{R}^d``,
       let ``x_{\mathrm{c}}(x) = \tfrac{1}{N} \sum_b x_b`` and ``h_\theta(\tilde{x}, p) =
       p_{\mathrm{scale}} \sum_I \theta_I B_I(\tilde{x}, p / p_{\mathrm{scale}})``. Then
       ``H^{\mathbb{R}^d}_\theta(x, p) := \sum_a h_\theta(x_a - x_{\mathrm{c}}(x), p_a)`` is
       invariant under permutations and translations, costs ``O(N)`` operations, and generates the
       recoil-coupled equations
       ```math
       \nabla_{p_a} H^{\mathbb{R}^d}_\theta = \nabla_p h_\theta(x_a - x_{\mathrm{c}}, p_a) , \qquad
       \nabla_{x_a} H^{\mathbb{R}^d}_\theta = \nabla_{\tilde{x}} h_\theta(x_a - x_{\mathrm{c}}, p_a)
       - \frac{1}{N} \sum_{c=1}^{N} \nabla_{\tilde{x}} h_\theta(x_c - x_{\mathrm{c}}, p_c) .
       ```

!!! details "Proof"
    Relabelling the summation indices ``(a, b) \mapsto (\sigma(a), \sigma(b))`` leaves
    ``H_\theta`` unchanged, which proves (i). Since ``(x_a + \xi) - (x_b + \xi) = x_a - x_b``,
    ``H_\theta(x + \xi, p) = H_\theta(x, p)``, and differentiating at ``\xi = 0`` gives Noether's
    identity (ii). The identities
    ``\cos(k \cdot (x_a - x_b)) = \cos(k \cdot x_a) \cos(k \cdot x_b) + \sin(k \cdot x_a)
    \sin(k \cdot x_b)`` and ``\sin(k \cdot (x_a - x_b)) = \sin(k \cdot x_a) \cos(k \cdot x_b) -
    \cos(k \cdot x_a) \sin(k \cdot x_b)`` reduce the pair sum to quadratic forms in
    ``(C_{k,r}, S_{k,r})``, which proves (iii). For (iv), differentiating
    ``\tilde{x}_c = x_c - \tfrac{1}{N} \sum_b x_b`` with respect to ``x_a`` gives
    ``(\delta_{ca} - N^{-1}) I_d``, which produces the recoil-coupled equations. Unlike the
    ad-hoc field of the lemma, ``\nabla_{p_a} H^{\mathbb{R}^d}_\theta`` depends on all ``x_b``
    through ``x_{\mathrm{c}}(x)``, which restores the exact symmetry of
    ``\nabla^2 H^{\mathbb{R}^d}_\theta``.

In the package, [`PairwiseSplineHamiltonian`](@ref) is the pairwise form in 1D-1V with periodic
quartic B-splines ``\phi_m`` on ``[0, 2\pi)`` and clamped quartic B-splines ``\chi_r`` on
``[-1, 1]``. The verification uses ``10`` cells in ``x`` and ``8`` in ``\hat{p}``, so
``M_x = 10``, ``M_p = 12`` and ``N_\theta = 12 + 10 \cdot 12^2 = 1452``.
[`basis_hamiltonians`](@ref) evaluates ``\Psi_k``, [`basis_hamiltonian_gradients`](@ref) their
gradients, and [`pairwise_hamiltonian_derivatives`](@ref) the gradient and the Hessian of
``H_\theta``. In the pair sum, the factor ``\tfrac12`` and the sum over ordered pairs count each
unordered pair twice with ``\phi_m(x_a - x_b)`` and ``\phi_m(x_b - x_a)``; the self pairs
``a = b`` contribute ``\phi_m(0)``, a function of the momenta alone.

```@example rd-hamiltonians
H = PairwiseSplineHamiltonian(ex.Lx, p_scale; nx_cells = 10, np_cells = 8)
θ = randn(hamiltonian_parameter_count(H)) ./ 100
Hq, Hp, _ = pairwise_hamiltonian_derivatives(H, q, p, θ)
perm = sortperm([sin(17.3 * k) for k in 1:N])
Hq_p, Hp_p, _ = pairwise_hamiltonian_derivatives(H, q[perm], p[perm], θ)
Hq_s, _, _ = pairwise_hamiltonian_derivatives(H, q .+ 1.2345, p, θ)
(noether = sum(Hq), permutation = norm(Hq_p - Hq[perm]) + norm(Hp_p - Hp[perm]),
 translation = norm(Hq_s - Hq))
```

At the parameters of the verification step, the script gives ``|\sum_a \nabla_{x_a} H_\theta| =
7.6 \cdot 10^{-17}``, a translation error ``6.2 \cdot 10^{-16}`` and a permutation error
``1.8 \cdot 10^{-16}``; with eight relaxation steps of size ``0.02`` it gives
``4.86 \cdot 10^{-17}`` and ``2.10 \cdot 10^{-16}``.

The Fourier-moment form is [`FourierMomentHamiltonian`](@ref), with
[`collective_moments`](@ref); its [`ensemble_hamiltonian`](@ref) and
[`ensemble_hamiltonian_gradient`](@ref) go through the moments, and
[`fourier_pairwise_hamiltonian`](@ref) is the direct double sum, the reference:

```@example rd-hamiltonians
F = FourierMomentHamiltonian([1, 2, 3], p_scale; np_cells = 4)
θF = randn(hamiltonian_parameter_count(F))
(moments = ensemble_hamiltonian(F, q, p, θF), pairs = fourier_pairwise_hamiltonian(F, q, p, θF),
 noether = sum(ensemble_hamiltonian_gradient(F, q, p, θF)[1]))
```

The centroid-frame form is [`CentroidFrameHamiltonian`](@ref), on ``\mathbb{R}`` with clamped
B-splines in ``\tilde{x} \in [-X_{\max}, X_{\max}]``. The Jacobian of its recoil-coupled
gradient is symmetric, so ``\mathbb{J} \nabla H`` is a Hamiltonian field:

```@example rd-hamiltonians
C = CentroidFrameHamiltonian(4.0, p_scale; nq_cells = 6, np_cells = 4)
θC = randn(hamiltonian_parameter_count(C))
qc = q .- 3.0
D = ForwardDiff.jacobian(y -> vcat(ensemble_hamiltonian_gradient(C, y[1:N], y[(N + 1):end],
    θC)...), vcat(qc, p))
(noether = sum(ensemble_hamiltonian_gradient(C, qc, p, θC)[1]),
 asymmetry = norm(D - transpose(D)) / norm(D))
```

### [Remark (no circular-mean centering on the torus)](@id rd-circular-mean)

!!! note "Remark (no circular-mean centering on the torus)"
    On ``\mathbb{T}^1 = \mathbb{R} / L\mathbb{Z}``, replacing ``x_{\mathrm{c}}`` in the centroid
    frame by the circular mean
    ``\bar{\vartheta}(x) = \tfrac{L}{2\pi} \operatorname{atan2}\big(\sum_b \sin\tfrac{2\pi x_b}{L},
    \sum_b \cos\tfrac{2\pi x_b}{L}\big)`` is translation-covariant away from
    ``\hat{\rho}_1(x) := \tfrac{1}{N} \sum_b e^{\mathrm{i} 2\pi x_b / L} = 0``, but it becomes
    singular precisely when the markers are uniformly spaced in ``x`` (``\hat{\rho}_1 = 0``). The
    pairwise periodic-spline and Fourier-moment forms stay globally smooth for all configurations
    on ``\mathbb{T}^d``.

The gradient of the circular mean grows like ``1 / |\hat{\rho}_1|`` as the markers approach a
uniform spacing:

```@example rd-hamiltonians
ϑ(q) = atan(sum(sin, q), sum(cos, q))                     # L = 2π
uniform = collect(range(0, 2π; length = 9)[1:8])
bump = sin.(uniform .+ 0.4)
[(ε = ε, ρ̂₁ = abs(sum(cis, uniform .+ ε .* bump)) / 8,
  gradient = norm(ForwardDiff.gradient(ϑ, uniform .+ ε .* bump)))
 for ε in (1e-1, 1e-2, 1e-3)]
```
