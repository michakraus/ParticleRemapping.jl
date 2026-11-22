# [Properties, limitations and higher dimensions](@id rd-limitations-page)

```@meta
CurrentModule = ParticleRemapping.Redistribution
```

## [Properties](@id redistribution-properties)

| quantity | status |
|:--|:--|
| symplectic form ``\omega`` | exact for the frozen parameters: each step of the implicit midpoint or the symplectic Euler method with fixed ``\theta`` is symplectic, and so is the composition of substeps, followed by the gauge step ([parameter solver](@ref rd-linear-parameter-solver), [factorisation lemma](@ref factorisation-lemma)); the closed loop ``\theta = \Theta(z)`` of target matching is not symplectic (defect ``0.396`` in [Control E](@ref rd-closed-loop)) |
| total weight | exact: the KKT projection constrains ``\sum_a \bar{w}_a`` ([KKT projection](@ref rd-kkt-projection), [example](@ref rd-experiments-standard)) |
| total momentum ``P`` | exact when the step converges: Noether's identity holds for every ``\theta``, and the KKT projection makes the targets compatible ([KKT projection](@ref rd-kkt-projection), [schemes](@ref rd-experiments-schemes)); not kept when the forward Newton solve fails (``0.165`` in [K3](@ref rd-holes)) |
| kinetic energy | not kept: no constraint holds it; ``0.1323`` of ``1.0812`` in the [example](@ref rd-experiments-standard) |
| grid charge | not kept: no constraint holds it |
| grid current | not kept: no constraint holds it |
| represented density (consistency) | the new weights are ``\bar{w}_a = \Delta\bar{\Omega}_a f_h^{\mathrm{KKT}}(\bar{x}_a, \bar{v}_a)``: the positive spline projection of the markers, corrected to the discrete mass and momentum, sampled at the targets with the volumes of their number density ([KKT projection](@ref rd-kkt-projection)); the markers sit at the targets only when the step reaches them |
| equidistribution and hole filling | the relaxation lowers the clustering (kernel-density CV ``0.2199 \to 0.1939``) on 16 markers ([projection](@ref rd-projection)); the step does not fill the hole of the 244-marker [hole-filling case](@ref rd-holes) (K3) |
| parameters computed from the state | yes: the targets, the weights, ``p_{\mathrm{scale}}`` and ``\theta`` come from the markers; the map is symplectic for their frozen values (open loop, [closed loops](@ref rd-closed-loop)) |
| cost of one remap | relaxation ``O(I_r N^2 d)``; projection and KKT weights ``O(N (k+1)^{2d} + \lvert\mathcal{I}\rvert)``; assembly of ``G`` ``O(N^2 (k_x+1)^d (k_p+1)^{2d})`` into a dense ``2Nd \times N_\theta`` matrix; the minimum-norm solve ``O((2Nd)^2 N_\theta)``; each Newton iterate ``O(N^2 (k+1)^{3d} + (2Nd)^3)``; 0.04 s for ``N = 16``, 44 s for ``N = 244`` ([higher dimensions](@ref redistribution-higher-dimensions)) |
| dimensions implemented in the package | ``d = 1`` (1D-1V); `remap` rejects other markers (`test/redistribution/method.jl`) |

## [Limitations](@id redistribution-limitations)

- **The step reaches its targets only on the identity-connected branch.** The parameters make
  ``\bar{z}`` a root of the step equations, but the forward Newton solve follows the branch from
  ``z``; when the contraction condition fails by much, it ends on another root (the single
  symplectic Euler step of the [experiments](@ref rd-experiments-schemes)) or not at all (K3).
  Two close markers separated by ``O(1)`` in one step force
  ``\| \nabla^2 H_{\theta^*} \| \gtrsim \varepsilon^{-1}`` ([remark](@ref rd-coarse-basis) (a));
  substeps reduce the displacement of each step.
- **The relaxation does not scale with ``N``.** The repulsion is a sum over all pairs, so the
  force on a point grows with the number of neighbours within the kernel widths, and a fixed step
  moves many points out of the velocity box ``\Omega_v`` (K2). The kernel widths match the
  spacing of the 16-marker lattice.
- **The targets must lie in the spline box.** The projection, the number density and the KKT
  weights live on ``\Omega_v = [-V_{\max}, V_{\max}]^d``; no basis function is nonzero at a
  target outside it, and the KKT solve fails (K2).
- **The weights are positive only when ``c^{\mathrm{KKT}} > 0``.** The KKT correction keeps the
  discrete moments exactly, and its metric protects small coefficients, but positivity holds only
  when the quadrature discrepancy of the moments is small enough
  ([KKT projection](@ref rd-kkt-projection) (iii)).
- **Symplectic in the open-loop sense only.** The parameters are computed from the state; target
  matching in ``(x, v)`` with weights tied to the density contracts an open set of states towards
  a lower-dimensional set of layouts, which no closed-loop symplectic map does
  ([cotangent lifts](@ref rd-cotangent-lift)). The variational closure is symplectic as a closed
  loop, but it spreads clusters by a repulsion and does not match targets.
- **Targets in ``(\bar{x}, \bar{v})``, not in ``\bar{p}``.** Weights tied to the density fold the
  velocity–momentum map ([fold](@ref velocity-momentum-fold)); the construction therefore
  evaluates ``\bar{w}_a`` and ``\bar{p}_a`` forward from ``(\bar{x}_a, \bar{v}_a)``, and the
  coarse-basis regime recovers velocities along the branch of the fold
  ([remark](@ref rd-coarse-basis) (b)).
- **Only the weight and the momentum are constrained.** The kinetic energy, the grid charge and
  the grid current change.
- **The pairwise basis grows with the square of the momentum basis.**
  ``N_\theta = M_p + M_x M_p^2``; the pair sums cost ``O(N^2)``, and the Newton solve works with
  the dense ``2Nd \times 2Nd`` Hessian ([higher dimensions](@ref redistribution-higher-dimensions)).

## [Known issues](@id redistribution-known-issues)

The entries of `KNOWN_ISSUES.md` of the package that concern this chapter; the notes have no
entry for it.

- **K1** (upstream): SimpleSplines 0.4 returns NaN derivatives of a periodic basis at the knot
  ``x = 0`` for a dual number, so the pairwise Hamiltonian evaluates the self pair at the plain
  number ``0``.
- **K2**: the relaxation step of `SplineRedistribution` does not scale with the number of markers;
  with 900 markers the targets leave the velocity box and `remap` fails in the KKT solve
  ([experiment](@ref rd-experiments-failures)).
- **K3**: `remap` with `SplineRedistribution` does not report a step that misses its targets; on
  the 244 markers of the hole-filling problem the target error is ``1.98 \cdot 10^{8}``
  ([experiment](@ref rd-holes)).
- **K4**: the verification case violates the contraction condition (``2.266`` for the midpoint
  step, ``1.391`` for one symplectic Euler step), which is sufficient, not necessary
  ([schemes and substeps](@ref rd-experiments-schemes)).

## [Open questions](@id redistribution-open-questions)

- **The coarse-basis regime at large ``N``.** With ``N_\theta \ll 2 N d - d`` the regularised
  parameter solve reaches the targets only approximately, and the velocities are recovered along
  the branch of the fold ([remark](@ref rd-coarse-basis) (b)). The package tests the recovery on
  16 markers; it has no experiment of this regime at large ``N``.
- **A relaxation that scales.** Kernel widths that follow the spacing ``N^{-1/(2d)}`` of an
  equidistribution, a step size that bounds the largest displacement, and a bound on the
  displacement of one step that keeps the contraction condition: the remark (a) gives the
  direction, ``\| \nabla^2 H_{\theta^*} \| \gtrsim \varepsilon^{-1}`` for a pair at distance
  ``\varepsilon``, but no a priori bound of the contraction number by the displacements.
- **A variational closure that matches targets.** The [variational closure](@ref variational-closure)
  makes the closed loop symplectic; whether a potential ``\mathcal{V}(\theta)`` exists whose
  stationary points move the markers towards the relaxed targets is open.
- **More moments in the KKT projection.** The kinetic energy at the targets,
  ``\sum_a \tfrac{m}{2} \bar{w}_a |\bar{v}_a|^2``, and the grid charge
  ``\sum_a \bar{w}_a S(x_j - \bar{x}_a)`` are linear in the spline coefficients ``c``, like the
  mass and the momentum, so they can enter ``A`` as further rows with the same KKT formula. Open
  are the rank of the larger ``A`` and the positivity of ``c^{\mathrm{KKT}}`` under more
  constraints; the kinetic energy of the result equals the constrained value only when the step
  reaches the targets.
- **The collective form in the solver.** The Fourier-moment Hamiltonian costs ``O(N)`` per
  evaluation ([symmetric Hamiltonians](@ref rd-symmetric-hamiltonians) (iii)); the package
  evaluates it, but the parameter solve and the Newton step use the pairwise spline form only.
- **Repeated remaps.** The long-time behaviour of a simulation that redistributes with
  state-dependent parameters ([open questions](@ref open-questions) of the comparison).

## [Higher dimensions](@id redistribution-higher-dimensions)

### Statements in ``d`` dimensions

Every statement of the chapter is written for ``x \in \mathbb{T}^d`` and
``v \in \Omega_v = [-V_{\max}, V_{\max}]^d``: the pull-back of ``\omega`` and the factorisation,
the fold determinant ``(m \Delta\Omega)^d f_h^{d-1} (f_h + \bar{v} \cdot \nabla_{\bar{v}} f_h)``,
the moment conservation of the ``L^2`` projection (Marsden's identity in each ``v_i``), the KKT
projection with ``d + 1`` constraints and ``\operatorname{rank} A = d + 1`` when the ``\bar{v}_a``
span ``\mathbb{R}^d`` affinely, the obstructions for uncoupled Hamiltonians, the symmetry and
Noether's identity of the pairwise, Fourier-moment and centroid-frame Hamiltonians with tensor
bases, the linear system for ``\theta`` on the subspace ``\mathcal{M}_0`` of dimension
``2 N d - d``, the contraction condition, and the closed-loop identity, which does not depend on
``d``.

### Implementation in the package

``d = 1``: [`PhaseSpaceSplines`](@ref) is a 1D-1V basis, [`PairwiseSplineHamiltonian`](@ref)
takes scalar positions and momenta, [`declustering_relaxation`](@ref) works in ``(x, v)`` with one
component each, and `remap` with [`SplineRedistribution`](@ref) throws an `ArgumentError` for
markers with matrix positions. The fold determinant is checked in two dimensions
(`test/redistribution/canonical.jl`).

### The cost in 3d3v

Let ``n_x`` and ``n_v`` be the cells per direction of the phase-space spline of degree ``k``
(periodic in ``x``, clamped in ``v``), ``n_{x,c}`` and ``n_{p,c}`` those of the pairwise
Hamiltonian of degree ``k_H`` (periodic in ``x``, clamped in ``\hat{p}``), and ``I_r`` the number
of relaxation steps. The package defaults are ``n_x = 8``, ``n_v = 10``, ``k = 3``,
``n_{x,c} = 10``, ``n_{p,c} = 8``, ``k_H = 4``, ``I_r = 10``. The script
`scripts/redistribution/failure_cases.jl` prints the sizes below for ``d = 1, 2, 3``.

| step | size | cost | ``d = 1`` | ``d = 3`` (3d3v) |
|:--|:--|:--|:--|:--|
| phase-space spline | ``\lvert\mathcal{I}\rvert = n_x^d (n_v + k)^d`` | | ``104`` | ``1\,124\,864`` |
| evaluation matrix ``E``, number density, lumped projection | ``N \times \lvert\mathcal{I}\rvert``, ``(k+1)^{2d}`` nonzeros per row | ``O(N (k+1)^{2d} + \lvert\mathcal{I}\rvert)`` | ``16`` per row | ``4096`` per row |
| KKT projection | ``A`` is ``(d+1) \times \lvert\mathcal{I}\rvert`` | ``O(N (k+1)^{2d} + d^2 \lvert\mathcal{I}\rvert + d^3)`` | ``2`` rows | ``4`` rows |
| relaxation | ``N`` points in ``2d`` dimensions | ``O(I_r N^2 d)``, pair sums | | |
| Hamiltonian basis | ``M_x = n_{x,c}^d``, ``M_p = (n_{p,c} + k_H)^d``, ``N_\theta = M_p + M_x M_p^2`` | | ``N_\theta = 1452`` | ``N_\theta \approx 2.99 \cdot 10^{9}`` |
| pair sums of ``\Psi_k``, ``\nabla \Psi_k`` | ``(k_H + 1)^{3d}`` nonzero ``(m, l, r)`` per ordered pair | ``O(N^2 (k_H+1)^{3d})`` | ``125`` per pair | ``1\,953\,125`` per pair |
| linear system for ``\theta`` | ``G`` is ``2Nd \times N_\theta``, rank ``2Nd - d`` | dense QR ``O((2Nd)^2 N_\theta)``; minimum norm through ``G G^{\mathsf{T}}``: ``O(N^3 d^2 (k_H+1)^{3d})`` | ``31`` at ``N = 16`` | ``100\,663\,293`` at ``N = 16^6`` |
| each Newton iterate | Hessian ``2Nd \times 2Nd``, dense | ``O(N^2 d^2 (k_H+1)^{3d} + (2Nd)^3)`` | ``32 \times 32`` | |

The rank ``2 N d - d`` of ``G`` counts the unknowns: ``2d`` coordinates per marker, less the
``d`` components of the total momentum, which Noether's identity fixes. The minimum-norm solution
needs ``N_\theta \ge 2 N d - d``. The Gram matrix ``G G^{\mathsf{T}}`` has
``2Nd \times 2Nd`` entries, each a sum over the parameters that two rows share; a row of
``\nabla_{x_a}`` touches the pairs ``(a, b)`` for all ``b``, so it has
``O(N (k_H+1)^{3d})`` nonzeros.

**The coarse-basis regime.** With ``N_\theta \ll 2 N d - d`` ([remark](@ref rd-coarse-basis)
(b)) the Tikhonov solve forms the normal equations, ``O(2Nd \, N_\theta^2 + N_\theta^3)``, the
targets are reached only approximately, and each marker recovers its velocity by a Newton solve
in ``d`` unknowns on the branch of the fold, ``O(N ((k+1)^{2d} + d^3))`` per iteration. The
Fourier-moment form with the wave vectors ``\mathcal{K}`` has
``N_\theta = M_p + 2 \lvert\mathcal{K}\rvert M_p^2`` and evaluates ``H_\theta`` and its gradient in
``O(N \lvert\mathcal{K}\rvert M_p^2)``, without pair sums.

### Obstacles and a route to 3d3v

The obstacles, in the order in which they bind:

1. **The parameter count.** ``M_p = (n_{p,c} + k_H)^d`` enters squared: ``N_\theta \approx
   3 \cdot 10^9`` for the defaults in 3d3v, which no dense matrix holds. The full-rank regime
   ``N_\theta \ge 2Nd - d`` is the regime of exact targets, and it is out of reach in 3d3v.
2. **The pair sums and the dense Hessian.** ``O(N^2)`` pairs and a ``2Nd \times 2Nd`` Newton
   matrix limit the step to small ``N``; the 244-marker case already takes 44 s in ``d = 1``.
3. **The relaxation** (K2) and **the caustics** (K3, K4): both grow with ``N`` at fixed
   parameters.

A route that the statements of the chapter support: the coarse collective basis — the
Fourier-moment form, with ``O(N)`` evaluations — with the Tikhonov solve and the velocity
recovery of remark (b); substeps, each with its contraction number below one; a relaxation with
kernel widths and step size tied to the spacing, and a kernel cutoff that makes it ``O(N)``; the
KKT projection unchanged, since its cost is linear in ``N`` and in
``\lvert\mathcal{I}\rvert`` and its solve is ``(d+1) \times (d+1)``. Each piece exists in the
package for ``d = 1`` except the coarse-basis solver with the Fourier-moment form; none is tested
at large ``N``.
