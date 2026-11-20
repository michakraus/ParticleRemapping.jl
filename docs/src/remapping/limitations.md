# [Properties, limitations and higher dimensions](@id remap-limitations-page)

```@meta
CurrentModule = ParticleRemapping.Remapping
```

## [Properties](@id remapping-properties)

The numbers are those of the standard examples of the [overview page](@ref remap-overview):
``6400`` markers for the gauge rezoning and ``10^4`` markers for the canonical transport with the
admissible cotangent lift.

| quantity | gauge rezoning | canonical transport |
|:--|:--|:--|
| symplectic form ``\omega`` | exact for every stretch, also one that depends on the whole ensemble: the map is a gauge step ([exact properties](@ref remap-gauge), `test/remapping/canonicity.jl`) | exact for a fixed canonical ``\psi``; the matching is a gauge step at the mapped ``(X, P)`` ([discrete algorithm](@ref remap-matching)) |
| total weight | exact: the charge constraints sum to it by the partition of unity ([example](@ref remap-overview)) | up to the quantile error: it changes by ``0.005`` of ``2\pi`` ([example](@ref remap-overview)) |
| total momentum ``P`` | exact: every ``p_a`` is kept ([example](@ref remap-overview)) | every ``P_a`` of ``\psi`` is kept, so ``P`` is that of ``\psi``; the cotangent lift changes it by ``0.0034`` ([example](@ref remap-overview)) |
| kinetic energy | exact: the constraints ``c_i^{(2)} = 0`` sum to it ([constraints](@ref remap-gauge), [example](@ref remap-overview)) | up to the quantile error: ``-0.0063`` ([example](@ref remap-overview)) |
| grid charge | exact by constraint, ``2.6 \times 10^{-15}`` ([experiments](@ref remap-gauge-example)) | up to the quantile error: ``0.078`` relative at the worst node ([example](@ref remap-overview)) |
| grid current | exact for every stretch, ``1.2 \times 10^{-16}`` ([exact properties](@ref remap-gauge), [experiments](@ref remap-gauge-example)) | set by ``\psi``: exact for ``\psi = \mathrm{id}``, ``0.13`` relative at the worst node after the lift ([example](@ref remap-overview)) |
| represented density (consistency) | kept in the constrained moments only; the third moment changes by ``0.22`` ([experiments](@ref remap-gauge-example)) | consistent with ``\hat{f}`` up to the quantile error: binned ``L^1`` error ``0.0117`` against the sampling noise ``0.0273`` ([experiments](@ref remap-transport-example)) |
| equidistribution and hole filling | along the velocity rays of each cell and sign: the gap ratio falls from ``77.5`` to ``27.5`` before the fold at ``\tau = 0.594`` ([experiments](@ref remap-gauge-example)); no move in ``x`` and no crossing of ``v = 0`` | not an aim: the admissible lift moves markers in ``x`` and keeps the area of every level set of the canonical marker density ([limit of equidistribution](@ref remap-limit)) |
| parameters computed from the state | the targets and the stretches; the remap is symplectic for every state dependence ([parameters that depend on the state](@ref state-dependent-parameters)) | the weights, by a gauge step, symplectic for every state dependence; ``\psi`` is given, and a ``\psi`` computed from the state is symplectic for its frozen value only |
| cost of one remap | per Newton iteration ``O(\lvert\mathcal{K}\rvert N)`` for ``c``, ``\mathsf{G}``, ``\nabla J``, ``H`` and ``O(\lvert\mathcal{K}\rvert^2 N)`` for the Schur complement; its solve ``O((\lvert\mathcal{K}\rvert N_g)^3)`` as a dense matrix in the package (K9); times the Newton iterations of each homotopy step and at least ``1/\Delta\tau_0`` steps ([higher dimensions](@ref remapping-higher-dimensions)) | ``\psi`` once per marker (the cotangent lift: ``t/\Delta t`` Runge–Kutta steps); the sort per cell and branch, ``O(N \log (N / N_g))``; one inverse of ``\hat{F}_c`` per marker |
| dimensions implemented in the package | 1D-1V; [`GaugeRezoningProblem`](@ref) rejects other markers | 1D-1V; `remap` rejects other markers |

## [Limitations](@id remapping-limitations)

**Gauge rezoning.**

- **Rays.** A marker moves only along its velocity ray. It cannot cross ``v = 0``, and it does
  not move in ``x``. A hole that lies across ``v = 0``, or a hole in ``x``, needs a canonical map
  in addition. A per-marker momentum translation ``p_a \mapsto p_a + \kappa_a`` is canonical for
  fixed ``\kappa_a``; it keeps the grid current if ``\sum_a \kappa_a S_i(x_a) = 0`` for every
  node, and the kinetic-energy constraint then includes the translation.
- **The fold.** The homotopy stops where some curvature ``H_a`` reaches zero: no nearby
  minimiser continues the branch ([second-order condition](@ref remap-gauge),
  [one constraint in closed form](@ref remap-charge-example)). The clustered example stops at
  ``\tau = 0.594``, and each further constraint lowers the ``\tau`` reached: ``0.402`` with
  ``\mathcal{K} = \{0, 2, 3\}`` ([experiments](@ref remap-gauge-example)).
- **Representation of ``f``.** The weight follows from the move: ``W_a = p_a / (m V_a)``. The
  construction does not ensure that ``W_a`` is close to ``f(x_a, V_a)`` times the phase-space area
  of the marker; the moments outside ``\mathcal{K}`` change. Targets taken from an estimate of
  ``f``, for example those of the [current-quantile matching](@ref remap-matching), combine the
  two aims ([combining the two constructions](@ref remap-comparison)).
- **Targets.** The default targets keep the order and the range of the speeds of each group of
  cell and sign; a group of fewer than three markers keeps its speeds. With few markers a single
  constraint need not keep the order of the speeds
  ([one constraint in closed form](@ref remap-charge-example)).
- **Local energy constraints.** The package constrains the kinetic energy at every node. The
  single global constraint ``\sum_i c_i^{(2)} = 0`` suffices for exact energy conservation and
  leaves more freedom; the package does not implement it.

**Canonical transport.**

- **Compatibility.** Consistent weights exist only if ``\psi`` preserves the sign-resolved
  current profile ([compatibility](@ref remap-transport)); after an inadmissible lift the
  compatibility residual reaches ``5.4`` and the representation of ``f`` fails
  ([experiments](@ref remap-transport-example)).
- **All freedom is in ``\psi``.** For a given ``\psi`` the weights are unique, and for
  ``\psi = \mathrm{id}`` they are the old ones ([uniqueness](@ref remap-transport)).
- **No equidistribution.** A smooth canonical map with consistent weights cannot make the marker
  distribution uniform, and it cannot fill a region in which the canonical marker density is
  small; it can only move it ([the limit](@ref remap-limit)).
- **One sheet.** The matching treats one layer of the velocity–momentum sheet
  ([the fold](@ref velocity-momentum-fold)); a marker with ``P = 0`` lies on no branch, and
  [`current_quantile_match`](@ref) rejects it.
- **The quantile rule.** The total weight, the kinetic energy and the grid charge are reproduced
  only to the accuracy of the quantile rule, which decreases with the number of markers per cell
  and branch ([five markers in one cell](@ref remap-five-example),
  [example](@ref remap-overview)). The slowest marker of a branch gets the largest change of
  weight, since both its ``P`` and its ``V`` are small ([example](@ref remap-overview)).
- **The target density.** The package implements the cumulative current of a separable
  density ``\hat{f} = J(x) M(v)`` only ([`SeparableMaxwellianTarget`](@ref)); an estimate of
  ``\hat{f}`` from the markers before the remap fits the algorithm but is not implemented.

## [Known issues](@id remapping-known-issues)

- **K9** — the Newton step of [`GaugeRezoning`](@ref) with `solver = :schur` factorises the
  Schur complement as a dense matrix ([`gauge_newton_step`](@ref)), at the cost
  ``O((\lvert\mathcal{K}\rvert N_g)^3)``, although the matrix is sparse and periodic
  block-tridiagonal ([the Newton step](@ref remap-gauge)).

No entry of the notes' `KNOWN_ISSUES.md` concerns this chapter.

## [Open questions](@id remapping-open-questions)

- **Admissible generators.** Which admissible generators ``h`` reduce clustering most, and how
  much such a remap gains within the [limit of equidistribution](@ref remap-limit).
- **Folded sheets.** How the matching splits ``f`` between the layers of a folded sheet
  ([the fold](@ref velocity-momentum-fold)).
- **Several velocity dimensions.** The matching along rays holds in the continuum
  ([higher dimensions](@ref remapping-higher-dimensions)); a discrete matching with direction
  bins, its error in the angular resolution, and the admissible canonical maps in ``d``
  dimensions are open.
- **Repeated remaps.** The long-time behaviour of a simulation with repeated remaps, given that a
  state-dependent remap breaks the symplecticity of the composite map
  ([two readings](@ref two-readings)).
- **The solver in 3d3v.** The number of conjugate-gradient iterations, or the fill-in of a sparse
  factorisation, of the Schur complement on ``N_g^3`` nodes; the package measures neither.
- **A larger reach.** Whether the global energy constraint, or targets from an estimate of ``f``,
  move the fold of the homotopy to a larger ``\tau``.
- **A dynamical test.** Holes that a simulation produces, rather than the clustered and sampled
  sets of the examples.

## [Higher dimensions](@id remapping-higher-dimensions)

### What holds in d dimensions

Let ``x_a, v_a \in \mathbb{R}^d`` and ``p_a = m w_a v_a``. The form
``\omega_N = \sum_a m (w_a \, \mathrm{d} x_a \wedge \mathrm{d} v_a + \mathrm{d} x_a \wedge v_a
\, \mathrm{d} w_a)``, with the sum over the components, has the kernel spanned by
``w_a \partial_{w_a} - v_a \cdot \partial_{v_a}``, and the proof of the
[factorisation lemma](@ref factorisation-lemma) uses only that ``(x, v, w) \mapsto (x, p)`` is a
submersion with connected fibres. Both constructions therefore keep their form.

**Gauge rezoning.**

1. The gauge step stretches each velocity vector along its ray, ``V_a = \lambda_a v_a``,
   ``W_a = w_a / \lambda_a``: one scalar ``\mu_a = \ln \lambda_a`` per marker, as in 1D. It keeps
   every ``x_a`` and ``p_a``, and so the grid current, a vector, for every ``\mu``.
2. The identity of the [exact properties](@ref remap-gauge) holds with
   ``\varphi(x, e^{\mu_a} v_a)``, and the moment of order ``k`` is the tensor
   ``c_i^{(k)} = \sum_a w_a S_i(x_a) \, v_a^{\otimes k} (e^{(k-1)\mu_a} - 1)``. The charge is one
   scalar per node, the kinetic energy (the trace of ``k = 2``) one scalar, the stress tensor
   ``d(d+1)/2`` scalars and the heat flux ``d``. With ``n_K`` scalar moments per node and
   tensor-product hat functions on ``N_g^d`` nodes there are ``n_K N_g^d`` constraints.
3. Each constraint is still a sum of single-marker terms in the scalar ``\mu_a``, so the
   stationarity decouples marker by marker and the Hessian of the Lagrangian stays diagonal,
   ``H_a = w_a [1/s^2 + \sum_k (k-1)^2 \, A_a^{(k)} \cdot v_a^{\otimes k} \, e^{(k-1)\mu_a}]``,
   with the multipliers interpolated from the ``2^d`` nodes of the cell of ``x_a``. The
   second-order condition and the argument for the fold use only the diagonal Hessian and hold
   unchanged.
4. The Schur complement ``\mathsf{G} \operatorname{diag}(H)^{-1} \mathsf{G}^T`` couples two nodes
   when they share a cell: ``3^d`` nodes around each, in blocks of size ``n_K``; it is symmetric
   and positive definite under the same conditions.
5. A gauge step changes the speed only. The targets become equispaced speeds in each group of
   cell and direction, the direction bins taking the place of the sign of ``v``.

**Canonical transport.** A canonical map in ``2d`` dimensions, such as the cotangent lift
``X = \varphi(x)``, ``P = D\varphi(x)^{-T} p``, moves the markers. The gauge step keeps the
direction ``\theta = P / |P|``, so with ``V = r \theta`` and ``P = \rho \theta``,
``\rho = m W r``, the Jacobian of ``V \mapsto P`` has the determinant
``m^d W^{d-1} \partial_r (W r)``, and the counting argument of the
[matching proposition](@ref remap-transport) gives, along every ray,

```math
\int_0^r f'(X, r' \theta) \, r'^{\,d} \, \mathrm{d} r'
= \frac{1}{m} \int_0^{\rho(r)} \tilde{n}(X, \rho' \theta) \, \rho'^{\,d} \, \mathrm{d} \rho' .
```

For ``d = 1`` this is the cumulative current of the proposition. The matching therefore holds
along rays: ``f' = f`` if and only if the radial cumulative current
``\int_0^r f \, r'^{\,d} \, \mathrm{d} r'`` matches the radial cumulative momentum at
``\rho = m W r`` for every ``X`` and ``\theta``. Consistent weights exist if and only if ``\psi``
preserves the *direction-resolved* current profile
``\int_0^\infty f(X, r\theta) \, r^d \, \mathrm{d} r`` for every ``X`` and ``\theta``, the
weights are unique for a given ``\psi``, and the limit of equidistribution holds along each ray
with ``\partial_r (f r) > 0`` in place of ``\partial_V (f V) > 0``. The test
`test/remapping/transport.jl` checks the identity along three rays in ``d = 2``, with the
determinant from automatic differentiation, to ``10^{-12}``; the frozen-weight determinant
``m^d W^d`` fails it by more than ``10^{-3}``.

### What the package implements

Both constructions in 1D-1V only: [`PeriodicGrid`](@ref ParticleRemapping.PeriodicGrid) is
one-dimensional with linear hat functions, [`GaugeRezoningProblem`](@ref) and
`remap(::CurrentQuantileMatching, …)` reject markers with matrix positions, and the target
density is the separable [`SeparableMaxwellianTarget`](@ref). The identity along rays in
``d = 2`` is checked as a continuum statement; no discrete matching in ``d > 1`` exists.

### The cost in 3d3v

One Newton step of the gauge rezoning with ``N`` markers, ``N_g^d`` nodes and ``n_K`` scalar
moments per node:

| step | count |
|:--|:--|
| unknowns of the KKT system | ``N + n_K N_g^d`` |
| constraints, multipliers, size of the Schur complement | ``n_K N_g^d`` |
| ``c``, ``\nabla J``, ``H`` and the non-zeros of ``\mathsf{G}`` | ``2^d n_K N`` |
| assembly of the Schur complement | ``4^d n_K^2 N`` products |
| non-zeros of the Schur complement | ``3^d n_K^2 N_g^d`` |
| dense factorisation, as in the package | ``\tfrac23 (n_K N_g^d)^3`` operations, ``8 (n_K N_g^d)^2`` bytes |
| one product with the sparse Schur complement | ``2 \cdot 3^d n_K^2 N_g^d`` operations |
| the step ``\Delta\mu`` | ``O(2^d n_K N)``, explicit per marker |

A remap takes a few Newton iterations per homotopy step (four or five in the
[eight-marker example](@ref remap-eight-example)) and at least ``1/\Delta\tau_0`` homotopy steps.
The script `scripts/remapping/cost_3d3v.jl` checks the 1D counts against the matrices of the
package and evaluates them for ``N_g = 32`` nodes per dimension and ``N = 10^7`` markers. For
``d = 3`` and the charge and the kinetic energy, ``n_K = 2``: ``65536`` constraints,
``1.6 \times 10^8`` non-zeros of ``\mathsf{G}``, ``2.6 \times 10^9`` products for the Schur
complement, ``3.5 \times 10^6`` non-zeros in it, and ``7.1 \times 10^6`` operations per product
with it; the dense factorisation would take ``1.9 \times 10^{14}`` operations and
``3.4 \times 10^{10}`` bytes. With the stress tensor, ``n_K = 7``, the dense factorisation would
take ``8.1 \times 10^{15}`` operations.

The canonical transport costs the map ``\psi`` once per marker: a cotangent lift in ``d``
dimensions integrates ``\varphi`` and the ``d \times d`` matrix ``D\varphi``, ``O(d^3)`` per
Runge–Kutta step. The matching sorts the markers of each cell and direction bin,
``O(N \log (N / (N_g^d n_\theta)))`` with ``n_\theta`` direction bins, and inverts one radial
cumulative current per marker.

### Obstacles and a route to 3d3v

- **The Schur solve.** The dense factorisation of the package (K9) does not reach 3d3v. The
  Schur complement is a sparse, symmetric, positive definite matrix with a ``27``-node stencil in
  blocks of size ``n_K``, so a sparse factorisation or the conjugate-gradient method applies; how
  many iterations the latter needs is open.
- **Markers per group.** The targets need at least three markers per group of cell and
  direction. With ``N = 10^7``, ``N_g = 32`` and the ``2^3 = 8`` orthants as directions, a group
  holds ``38`` markers on average.
- **Rays in three dimensions.** The gauge step cannot redistribute in angle; a hole in the
  direction of ``v`` needs a canonical map, as a hole across ``v = 0`` does in 1D.
- **The fold.** The argument is the same in every dimension, and every additional constrained
  moment lowers the ``\tau`` reached in 1D ([experiments](@ref remap-gauge-example)).
- **A discrete matching along rays.** The continuum identity holds; the discrete algorithm needs
  direction bins, the radial cumulative current of ``\hat{f}`` averaged over each bin, and
  canonical maps that keep the direction-resolved current. None of the three is implemented, and
  the admissible generators in ``d`` dimensions are not derived.

A route to 3d3v: tensor-product hat functions on ``N_g^3`` nodes, the scalar stretch per marker
and the diagonal Hessian as in 1D, and a sparse solve of the Schur complement for the gauge
rezoning; for the transport, cotangent lifts with ``D\varphi^{-T}``, direction bins, and the
radial cumulative current of the target density.
