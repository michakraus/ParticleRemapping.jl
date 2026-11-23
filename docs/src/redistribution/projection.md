# [Positive, moment-compatible spline projection](@id rd-projection)

```@meta
CurrentModule = ParticleRemapping.Redistribution
```

## The ``L^2`` projection and its positivity defect

Let ``\{\Lambda_I(x, v)\}_{I \in \mathcal{I}}`` be a tensor-product B-spline basis of degree
``(k_x, k_v)`` on the domain ``\Omega = \Omega_x \times \Omega_v``, periodic on
``\Omega_x = \mathbb{T}^d`` and on an open (clamped) knot vector on
``\Omega_v = [-V_{\max}, V_{\max}]^d`` that contains all marker velocities. With
``f_h(x, v) = \sum_I c_I \Lambda_I(x, v)``, the Galerkin ``L^2`` projection of the
Klimontovich distribution ``f_N`` solves

```math
\sum_{J \in \mathcal{I}} M_{IJ} \, c_J^{L^2} = b_I := \sum_{a=1}^{N} w_a \Lambda_I(x_a, v_a) ,
\qquad M_{IJ} = \int_\Omega \Lambda_I \Lambda_J \, \mathrm{d}x \, \mathrm{d}v .
```

The basis is [`PhaseSpaceSplines`](@ref) (cubic, ``8 \times 10`` cells on
``[0, 2\pi) \times [-3.5, 3.5]`` by default), the load vector ``b`` is
[`deposit_phase_space`](@ref), and the projection [`l2_spline_projection`](@ref).

### [Lemma (continuous moment conservation versus pointwise positivity)](@id rd-l2-lemma)

!!! info "Lemma (continuous moment conservation versus pointwise positivity)"
    Let ``k_v \ge 1`` with an open knot vector on ``\Omega_v``, and ``(x_a, v_a) \in \Omega`` for
    all ``a``. Then the ``L^2`` projection conserves the total mass and the total momentum
    continuously:
    ```math
    \int_\Omega f_h^{L^2} \, \mathrm{d}x \, \mathrm{d}v = \sum_{a=1}^{N} w_a , \qquad
    \int_\Omega m v \, f_h^{L^2} \, \mathrm{d}x \, \mathrm{d}v = \sum_{a=1}^{N} p_a .
    ```
    For spline degrees ``\max(k_x, k_v) \ge 1``, however, ``M^{-1}`` has negative off-diagonal
    entries, so ``f_h^{L^2}`` does not preserve positivity and can take negative values in
    depleted holes of phase space.

!!! details "Proof"
    By the partition of unity, ``\sum_I \Lambda_I(x, v) = 1`` on ``\Omega``, so the sum of the
    projection equations over ``I`` gives ``\int_\Omega f_h^{L^2} = \sum_a w_a``. On an open knot
    vector of degree ``k_v \ge 1``, Marsden's identity expresses each coordinate
    ``v_i = \sum_I \gamma^v_{I,i} \Lambda_I(x, v)`` through the velocity Greville abscissae
    ``\gamma^v_I`` (independent of ``x``). Contracting the projection equations with
    ``m \gamma^v_{I,i}`` gives ``\int_\Omega m v_i f_h^{L^2} = \sum_a m w_a v_{a,i} =
    \sum_a p_{a,i}``. The test file `test/redistribution/projection.jl` checks both moments to
    ``10^{-14}`` and the negative minimum (Control D).

[`spline_moments`](@ref) evaluates the mass, the momentum and the smallest value at the nodes of
the Gauss quadrature of the basis. On the clustered markers of the verification:

```@example rd-projection
using ParticleRemapping, CairoMakie
ex = redistribution_example()
mk = ex.markers
x, v, w, m = mk.x, mk.v, mk.w, mk.m
p = momentum(mk)
ps = PhaseSpaceSplines()
c_l2 = l2_spline_projection(ps, x, v, w)
l2 = spline_moments(ps, c_l2; m = m)
(mass_error = l2.mass - sum(w), momentum_error = l2.momentum - sum(p), minimum = l2.minimum)
```

The script `scripts/redistribution/verify_symplectic_redistribution.jl` gives the minimum
``\min f_h^{L^2} = -0.1327`` for these 16 markers, and ``-0.0551`` for 36 markers on a
``6 \times 6`` lattice.

## De-clustering relaxation and discrete moment compatibility

To redistribute the clustered markers smoothly, without combinatorial scrambling of the labels or
crossing trajectories, the points ``(x_a, v_a)`` move by a few damped steps of repulsive gradient
relaxation on a smoothed particle-count energy,

```math
\mathcal{E}_{\mathrm{rep}}(x, v) = \frac{1}{2} \sum_{\substack{a, b = 1 \\ a \ne b}}^{N}
K_\sigma(x_a - x_b, v_a - v_b) + \sum_{a=1}^{N} U_{\mathrm{wall}}(v_a) ,
```

where ``K_\sigma`` is a periodic-in-``x`` Gaussian kernel of widths ``(\sigma_x, \sigma_v)``
matched to the nominal spacing of the markers and ``U_{\mathrm{wall}}`` confines the velocities
softly inside ``\Omega_v``. Let ``(\bar{x}_a, \bar{v}_a)`` denote the relaxed, approximately
equidistributed coordinates.

The relaxation of the package ([`declustering_relaxation`](@ref)):

1. Start from ``(\bar{x}, \bar{v}) = (x, v)``.
2. For every marker, sum the repulsive forces
   ``\sum_{b \ne a} (\Delta x / \sigma_x^2, \Delta v / \sigma_v^2) \, K_\sigma`` with the periodic
   difference ``\Delta x \in [-L/2, L/2)`` and the kernel
   ``K_\sigma = \exp(-\tfrac12 (\Delta x/\sigma_x)^2 - \tfrac12 (\Delta v/\sigma_v)^2)``.
3. Add the wall force ``-k_{\mathrm{wall}} (v_a \mp v_{\mathrm{wall}})`` for
   ``|v_a| > v_{\mathrm{wall}}`` (``U_{\mathrm{wall}} = \tfrac{k_{\mathrm{wall}}}{2}
   (|v| - v_{\mathrm{wall}})^2``, ``v_{\mathrm{wall}} = 2``, ``k_{\mathrm{wall}} = 2.5``).
4. Move every point by the step size times its force (``0.006``), and repeat from 2 (ten steps).

The widths are ``\sigma_x = 2\pi / 4`` and ``\sigma_v = 4.8 / 4``, the spacing of the initial
lattice. [`declustering_energy`](@ref) is ``\mathcal{E}_{\mathrm{rep}}``, and
[`kernel_density_cv`](@ref) measures the clustering by the coefficient of variation of the
kernel density on a ``20 \times 20`` grid. The script gives ``0.2199 \to 0.1939`` for the ten
steps of size ``0.006``, and ``0.2199 \to 0.1580`` for eight steps of size ``0.02``.

```@example rd-projection
relaxed = declustering_relaxation(x, v; Lx = ex.Lx, σx = ex.σx, σv = ex.σv)
cv(h) = kernel_density_cv(mod.(h[1], ex.Lx), h[2]; Lx = ex.Lx, σx = ex.σx, σv = ex.σv)
plot_target_relaxation(relaxed.history; Lx = ex.Lx, cv = cv.(relaxed.history))
```

The figure shows (a) the paths of the 16 points during the ten relaxation steps, from the start
(dark grey dots) to the targets (open circles), and (b) the decrease of the coefficient of
variation of the kernel density. The relaxation moves points, not weighted markers; the weights
of the targets come from the KKT projection below. The animation repeats the relaxation with
``40`` steps of size ``0.02``, where the motion is visible:

```@example rd-projection
long = declustering_relaxation(x, v; Lx = ex.Lx, σx = ex.σx, σv = ex.σv, step = 0.02,
    iterations = 40)
animate_target_relaxation("redistribution-relaxation.mp4", long.history; Lx = ex.Lx)
nothing # hide
```

```@raw html
<video src="./redistribution-relaxation.mp4" controls loop muted playsinline style="width:100%"></video>
```

Since ``(\bar{x}_a, \bar{v}_a)`` is only *approximately* uniform, a rigid constant
``\Delta\Omega = |\Omega| / N`` would bias the local mass wherever residual variations of the
density persist. Instead, the deposition of unit weights at ``(\bar{x}_a, \bar{v}_a)`` into the
positive lumped-mass spline

```math
n_h(x, v) = \sum_{I} d_I \Lambda_I(x, v) , \qquad
d_I = \max\!\Big( \frac{1}{\mu_I} \sum_{a=1}^{N} \Lambda_I(\bar{x}_a, \bar{v}_a) ,
\; \varepsilon_n \Big) , \qquad \mu_I = \int_\Omega \Lambda_I \, \mathrm{d}x \, \mathrm{d}v > 0 ,
```

defines the local importance-sampling volume
``\Delta\bar{\Omega}_a := 1 / n_h(\bar{x}_a, \bar{v}_a) > 0``, which reduces to
``|\Omega| / N`` when the markers are uniformly distributed
([`number_density_volumes`](@ref), ``\varepsilon_n = 10^{-6}``).

### [Proposition (discrete momentum compatibility and weighted KKT projection)](@id rd-kkt-projection)

!!! info "Proposition (discrete momentum compatibility and weighted KKT spline projection)"
    Let ``(x, p) \mapsto (\bar{x}, \bar{p})`` be generated by symplectic Euler or implicit
    midpoint applied to a translation-invariant Hamiltonian ``H_\theta(x, p)``.

    1. *Necessary compatibility condition.* Every parameter vector ``\theta`` gives
       ``\sum_a \bar{p}_a = \sum_a p_a``. Hence, if ``\bar{w}_a = \Delta\bar{\Omega}_a
       f_h(\bar{x}_a, \bar{v}_a)`` and ``\bar{p}_a = m \bar{w}_a \bar{v}_a``, the discrete
       momentum identity
       ```math
       \sum_{a=1}^{N} \Delta\bar{\Omega}_a f_h(\bar{x}_a, \bar{v}_a) \, m \bar{v}_a
       = \sum_{a=1}^{N} p_a
       ```
       is necessary for a ``\theta`` that generates ``(\bar{x}, \bar{p})`` to exist.
    2. *Degeneracy of uniform velocity shifts.* An attempt to enforce the identity at a fixed
       ``f_h`` (with ``\operatorname{supp} f_h \subset \Omega_x \times
       \operatorname{int}\Omega_v``) and a constant ``\Delta\Omega`` by shifting
       ``\bar{v}_a \mapsto \bar{v}_a + \delta v`` has a nearly singular Jacobian, because the
       continuous integral of its derivative vanishes by integration by parts:
       ```math
       \int_\Omega m \big( f_h(x, v) I_d + v \otimes \nabla_v f_h(x, v) \big)
       \, \mathrm{d}x \, \mathrm{d}v = 0 .
       ```
    3. *Weighted KKT spline projection.* Start from the strictly positive lumped-mass
       coefficients ``c^{\mathrm{pos}}_I = \max(b_I / \mu_I, c_{\min}) > 0``, let
       ``E_{a,I} = \Lambda_I(\bar{x}_a, \bar{v}_a)`` be the ``N \times |\mathcal{I}|``
       evaluation matrix, and assemble the ``(d + 1) \times |\mathcal{I}|`` discrete moment
       matrix
       ```math
       A = \begin{pmatrix} \Delta\bar{\Omega}^{\mathsf{T}} E \\
       m \, (\operatorname{diag}(\Delta\bar{\Omega}) \bar{V})^{\mathsf{T}} E \end{pmatrix} ,
       \qquad
       \mu_{\mathrm{mom}} = \begin{pmatrix} \sum_a w_a \\ \sum_a p_a \end{pmatrix}
       \in \mathbb{R}^{d+1} ,
       ```
       where ``\bar{V} \in \mathbb{R}^{N \times d}`` has the rows ``\bar{v}_a^{\mathsf{T}}``. If
       the velocities ``\{\bar{v}_a\}`` span ``\mathbb{R}^d`` affinely and ``k_v \ge 1``, then
       ``A`` has full row rank ``d + 1``. With ``W = \operatorname{diag}((c^{\mathrm{pos}})^{-2})``,
       the unique minimiser of ``\tfrac12 (c - c^{\mathrm{pos}})^{\mathsf{T}} W
       (c - c^{\mathrm{pos}})`` subject to ``A c = \mu_{\mathrm{mom}}`` is
       ```math
       c^{\mathrm{KKT}} = c^{\mathrm{pos}} - W^{-1} A^{\mathsf{T}}
       \big( A W^{-1} A^{\mathsf{T}} \big)^{-1} \big( A c^{\mathrm{pos}} - \mu_{\mathrm{mom}} \big) ,
       ```
       and ``\bar{w}_a := \Delta\bar{\Omega}_a f_h^{\mathrm{KKT}}(\bar{x}_a, \bar{v}_a)`` conserves
       the discrete mass ``\sum_a \bar{w}_a = \sum_a w_a`` and the discrete momentum exactly, with
       ``\bar{w}_a > 0`` whenever the relative quadrature discrepancy of the moments is small
       enough for ``c^{\mathrm{KKT}} > 0``.

!!! details "Proof"
    (i) ``\sum_a \nabla_{x_a} H_\theta = 0`` gives ``\sum_a (\bar{p}_a - p_a) = -\tau \sum_a
    \nabla_{x_a} H_\theta = 0`` for both symplectic Euler and implicit midpoint. If ``\bar{z} - z``
    has ``\sum_a (\bar{p}_a - p_a) \ne 0``, it lies outside the column space of the discrete
    Hamiltonian operator (Control C, checked in `test/redistribution/projection.jl` and
    `test/redistribution/hamiltonians.jl`).

    (ii) Integrating ``v_i \partial_{v_j} f_h`` by parts on ``\Omega_v`` with vanishing boundary
    values gives ``-\delta_{ij} \int_{\Omega_v} f_h \, \mathrm{d}v``, which cancels
    ``f_h \delta_{ij}``.

    (iii) Let ``\Gamma \in \mathbb{R}^{|\mathcal{I}| \times (d+1)}`` have the columns
    ``\mathbf{1}_{|\mathcal{I}|}`` and ``\gamma^v_{\cdot,1}, \dots, \gamma^v_{\cdot,d}``. By the
    partition of unity and Marsden's identity, ``E \Gamma = [\mathbf{1}_{N}, \bar{V}]``, so
    ```math
    A \Gamma = \begin{pmatrix} 1 & 0 \\ 0 & m I_d \end{pmatrix}
    \begin{pmatrix} \mathbf{1}_{N}^{\mathsf{T}} \\ \bar{V}^{\mathsf{T}} \end{pmatrix}
    \operatorname{diag}(\Delta\bar{\Omega})
    \begin{pmatrix} \mathbf{1}_{N} & \bar{V} \end{pmatrix} ,
    ```
    which is congruent to the weighted empirical moment matrix of ``(1, \bar{v}_a)`` and strictly
    nonsingular whenever the ``\bar{v}_a`` span affinely and ``\Delta\bar{\Omega}_a > 0``. Hence
    ``\operatorname{rank} A = d + 1``, ``A W^{-1} A^{\mathsf{T}}`` is symmetric positive definite,
    and ``c^{\mathrm{KKT}}`` is the unique KKT stationary point. The metric
    ``W^{-1} = \operatorname{diag}((c^{\mathrm{pos}})^2)`` scales every component of the update
    ``c^{\mathrm{KKT}}_I - c^{\mathrm{pos}}_I`` with ``(c^{\mathrm{pos}}_I)^2``, which protects
    small tail coefficients from relative sign flips.

The weights of the targets ([`kkt_redistribution_weights`](@ref)):

1. Compute the positive prior ``c^{\mathrm{pos}}`` from the clustered markers
   ([`lumped_spline_projection`](@ref), ``c_{\min} = 10^{-4}``).
2. Compute the volumes ``\Delta\bar{\Omega}_a`` of the targets ([`number_density_volumes`](@ref)).
3. Assemble the evaluation matrix ``E`` ([`spline_evaluation_matrix`](@ref)) and the moment
   matrix ``A`` ([`kkt_moment_matrix`](@ref)).
4. Project: ``c^{\mathrm{KKT}}`` ([`kkt_spline_projection`](@ref)).
5. Weigh: ``\bar{w}_a = \Delta\bar{\Omega}_a f_h^{\mathrm{KKT}}(\bar{x}_a, \bar{v}_a)`` and
   ``\bar{p}_a = m \bar{w}_a \bar{v}_a``.

```@example rd-projection
x̄, v̄ = relaxed.x, relaxed.v
ΔΩ = number_density_volumes(ps, x̄, v̄)
c_pos = lumped_spline_projection(ps, x, v, w)
kkt = kkt_redistribution_weights(ps, c_pos, x̄, v̄, ΔΩ, sum(w), sum(p); m = m)
w̄ = kkt.weights
(min_weight = minimum(w̄), mass_error = sum(w̄) - sum(w),
 momentum_error = sum(w̄ .* m .* v̄) - sum(p),
 raw_momentum_defect = sum(kkt.raw_weights .* m .* v̄) - sum(p))
```

The last entry is Control C: the raw evaluation ``\Delta\bar{\Omega}_a f_h^{\mathrm{pos}}`` of
the prior misses the discrete momentum. The script gives ``\min \bar{w}_a = 0.01975``, mass and
momentum errors ``1.1 \cdot 10^{-16}`` and ``5.6 \cdot 10^{-17}``, and a raw momentum defect of
``0.01417``. With eight relaxation steps of size ``0.02`` it gives ``0.0168`` and ``0.01117``.
[The verification](@ref rd-verification) compares these numbers with those of the sources.

```@example rd-projection
plot_spline_projections(phase_space_spline_grid(ps, c_l2),
    phase_space_spline_grid(ps, kkt.coefficients); markers = mk, targets = (x̄, v̄))
```

The figure shows (a) the ``L^2`` projection of the 16 clustered markers (dark grey dots, the
area proportional to the weight) with its zero contour; the blue regions between the markers are
its negative overshoot (Control D); (b) the KKT projection on the same signed colour scale, with
the targets (open circles); (c) the smallest value over ``x`` of each projection as a function of
``v``: the ``L^2`` projection is negative at almost every ``v``, the KKT projection is positive.
The colour shows the sign of the projected density here, so the markers carry no weight colour.

The script `scripts/redistribution/make_figures.jl` writes these figures and the relaxation to
`results/redistribution/` (`redistribution-projections.pdf`, `redistribution-relaxation.pdf`,
`redistribution-relaxation.mp4`).
