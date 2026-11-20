# [Charge-conserving reweighting](@id projection-charge)

```@meta
CurrentModule = ParticleRemapping.Projection
```

The map moves the markers, and the [reweighting](@ref projection-final-reweighting) gives them
new weights. Both change the charge density that the field solver sees. Let ``\{\Lambda^g_i\}``
be the shape functions of the field solver in ``x``, periodic, with ``\Lambda^g_i \geq 0`` and
``\sum_i \Lambda^g_i = 1``, and let ``\rho_i = q \sum_a w_a \Lambda^g_i(x_a)`` be the grid charge
before the map, with the charge ``q`` of a marker of unit weight. Every choice of positive weights
after the map is a [gauge step](@ref gauge-step), so it keeps ``\omega`` and ``P``. The
charge-conserving reweighting uses this freedom: it replaces ``w^*_a`` in the final reweighting by

```math
w'_a = w^*_a \exp\Big( \sum_i \lambda_i \Lambda^g_i(x'_a) \Big) , \qquad
q \sum_a w'_a \Lambda^g_i(x'_a) = \rho_i \quad \text{for all } i ,
```

with one multiplier ``\lambda_i`` per shape function. Let ``A_{ia} = \Lambda^g_i(x'_a)``.

!!! info "Proposition (charge-conserving reweighting)"
    Let the positions ``x'_a`` after the map and the target weights ``w^*_a > 0`` be given.

    1. If ``\lambda`` solves the constraint, then ``w'_a > 0``, the grid charge after the map is
       ``\rho``, ``\sum_a w'_a = \sum_a w_a``, and the reweighting with ``w'`` keeps ``\omega``
       and ``P``.
    2. The weights ``w'`` are the unique minimiser of the relative entropy
       ``\sum_a ( w_a \ln (w_a / w^*_a) - w_a + w^*_a )`` among positive weights with the grid
       charge ``\rho``. The multipliers minimise the convex function

       ```math
       \Phi(\lambda) = \sum_a w^*_a \exp\big( (A^\top \lambda)_a \big) - \frac{1}{q} \lambda^\top \rho .
       ```

    3. A solution exists if and only if some positive weights ``u_a > 0`` have the grid charge
       ``\rho``. In particular, there is none if a shape function with ``\rho_i > 0`` contains no
       marker after the map in its support.

!!! details "Proof"
    (1) The exponential is positive. The constraint gives the charge. The partition of unity
    gives ``\sum_a w'_a = \sum_i \sum_a w'_a \Lambda^g_i(x'_a) = \sum_i \rho_i / q = \sum_a w_a``.
    The reweighting leaves ``(x, p)`` unchanged; it is a [gauge step](@ref gauge-step).

    (2) The relative entropy is strictly convex on positive weights and the constraints are
    linear. Its stationarity condition under the constraints is
    ``\ln (w_a / w^*_a) = (A^\top \lambda)_a``, which is the exponential form, so a solution is
    the unique minimiser. The gradient of ``\Phi`` is ``A w'(\lambda) - \rho / q`` and its Hessian
    ``A \, \mathrm{diag}(w') A^\top`` is positive semidefinite.

    (3) A solution is such a set of weights. Conversely, let ``A u = \rho / q`` with ``u_a > 0``,
    and let ``s = A^\top \lambda``. Then ``\lambda^\top \rho / q = s^\top u`` and
    ``\Phi = \sum_a (w^*_a e^{s_a} - u_a s_a)``. Every term is bounded below and grows without
    bound as ``|s_a| \to \infty``, so ``\Phi`` attains its minimum on the range of ``A^\top``, and
    the minimiser solves the constraint. If row ``i`` of ``A`` vanishes and ``\rho_i > 0``, no
    weights give the charge ``\rho_i``.

Newton's method on ``\Phi`` with a line search solves for ``\lambda``. The Hessian is banded with
the width of the shape functions, so one iteration costs ``O(N)`` operations. The correction acts
on the weights only: the momenta, and therefore ``P`` and the grid current
``(q/m) \sum_a p'_a \Lambda^g_i(x'_a)``, stay as the map made them. The velocities change by the
factor ``w^*_a / w'_a``. The additive correction
``w' = w^* + A^\top (A A^\top)^{-1} (\rho / q - A w^*)`` also gives the charge ``\rho``, but it can
give negative weights where small target weights meet an excess of charge.

## The code

[`charge_reweight`](@ref) solves for the multipliers on periodic cubic B-splines (``q = 1``) by
Newton's method with a backtracking line search on ``\Phi`` (Armijo constant ``10^{-4}``); it stops
when the charge residual falls below ``10^{-13} \max_i \rho_i`` and reports a Hessian that is not
positive definite, the case of item 3. [`additive_charge_correction`](@ref) is the additive
correction. The method type [`ChargeConservingReweighting`](@ref) applies the reweighting to
markers that carry the positions after the map and the target weights:

1. Compute the charge matrix ``A_{ia} = \Lambda^g_i(x'_a)`` ([`spline_charge_matrix`](@ref)).
2. Starting from ``\lambda = 0``, repeat: compute ``w'(\lambda) = w^* \exp(A^\top \lambda)`` and the
   residual ``F = A w' - \rho``; stop if ``|F| < 10^{-13} \max \rho``; solve
   ``A \operatorname{diag}(w') A^\top \delta = -F``; halve the step ``s`` from 1 until
   ``\Phi(\lambda + s\delta) \le \Phi(\lambda) + 10^{-4} s F^\top \delta``; set
   ``\lambda \leftarrow \lambda + s \delta``.
3. Keep the positions and the momenta, replace the weights by ``w'`` and the velocities by
   ``p / (m w')``.

The test `test/projection/charge.jl` uses 400 random markers on ``[0, 2\pi)``, 16 shape functions
and a kick and a drift with ``\kappa`` and ``\eta`` chosen by hand: Newton's method converges in 5
iterations; the grid charge is exact to ``7.1 \cdot 10^{-16}``; the total weight and ``P`` are
exact and ``\min w' > 0``; the gradient of the relative entropy is orthogonal to the null space of
``A`` to ``9.1 \cdot 10^{-16}``, and 50 feasible perturbations increase the relative entropy. The
cases that must fail: the target weights alone have the charge error ``0.29``; the additive
correction gives the exact charge and ``\min w' = -0.015`` where the exponential form gives
``2.4 \cdot 10^{-6}``; markers compressed to ``[0, \pi)`` have no solution, and the code reports
it.

```@example charge
using ParticleRemapping, Random
Random.seed!(20261121)
L, N = 2π, 400
x = L .* rand(N)
w = exp.(0.5 .* randn(N));
w .*= L / sum(w)
ρ = spline_grid_charge(UniformBasis(0.0, L / 16, 16, true), x, w)   # charge before the map
xs = mod.(x .+ 0.3 .* sin.(x), L)                                   # positions after a map
wstar = w .* exp.(0.3 .* randn(N));
wstar .*= sum(w) / sum(wstar)                                       # target weights
markers = Markers(xs, randn(N), wstar)
result = remap(ChargeConservingReweighting(L, 16), markers; reference_charge = ρ)
result.info
```
