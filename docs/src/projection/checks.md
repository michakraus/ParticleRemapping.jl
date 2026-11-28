# [Structural statements](@id projection-checks)

```@meta
CurrentModule = ParticleRemapping.Projection
```

This page gives the checks of the structural statements. Every statement of this chapter has a
test under `test/projection/`, and each test runs the case that must fail. Derivatives are
computed by forward-mode automatic differentiation, so the tolerances are set by rounding error.
The script `scripts/projection/checks.jl` runs the tests of the structural statements and prints
the measured numbers of the table.

## The structural statements

The first three rows check statements of
[Symplectic redistribution](@ref obstructions) with the functions of this package;
the other rows check the statements of this chapter.

| statement | result | case that must fail |
|:--|:--|:--|
| [compressibility](@ref weight-dependent-hamiltonians), ``d = 1, 2, 3`` | ``\lvert\operatorname{div} u - \{W, \partial_w \Psi\}\rvert \leq 1.1 \cdot 10^{-15}`` at random points, for random ``\Psi`` and ``W`` | opposite sign: ``\geq 7.8 \cdot 10^{-2}``; frozen weight: ``\geq 3.9 \cdot 10^{-2}`` |
| [fold](@ref velocity-momentum-fold), ``d = 1, 2, 3`` | relative error of the determinant ``\leq 1.7 \cdot 10^{-15}`` for random positive ``f``; Maxwellian: sign change at ``\lvert v\rvert = v_{\mathrm{th}}``, two roots at ``p = p_{\max}/2``, none at ``p = 3 p_{\max}/2`` | ``f^d`` in place of ``f^{d-1}``: ``\geq 2.4 \cdot 10^{-3}`` (``d \geq 2``) |
| [mean-field Hamiltonian](@ref one-map-for-all-markers), [holes persist](@ref holes-persist) | symplectic Euler and implicit midpoint on 274 markers: ``\lvert\psi(z_a) - z'_a\rvert \leq 6.8 \cdot 10^{-16}`` while markers move up to ``0.39``; ``\lvert\det D\psi - 1\rvert \leq 2.2 \cdot 10^{-16}``; hole area ``0.384845`` unchanged, no marker in the image of the hole | two coincident markers with weights ``0.5`` and ``2`` under a weight-dependent Hamiltonian: images differ by ``0.14`` |
| [kick–drift map](@ref projection-hamiltonian), ``d = 1, 2`` | ``\max \lvert DR_\theta^\top \mathbb{J} DR_\theta - \mathbb{J}\rvert \leq 2.9 \cdot 10^{-15}``; ``\lvert\Delta P\rvert \leq 8.9 \cdot 10^{-16}``; boost spread ``2.8 \cdot 10^{-17}``; kick and drift agree with symplectic Euler of ``H^K + H^D`` to ``3.1 \cdot 10^{-15}`` | naive kick: defect ``0.86`` and ``1.6`` with ``\lvert\Delta P\rvert \leq 8.9 \cdot 10^{-16}``; ``\theta = \theta(z)`` ([parameters that depend on the state](@ref state-dependent-parameters)): defect ``0.17`` and ``0.16``, against ``8.6 \cdot 10^{-15}`` with ``\theta`` frozen |
| [stage Jacobians](@ref projection-stage-jacobians) | error ``\leq 2.2 \cdot 10^{-16}``; at ``\pm v``: ``\det DK - 1 = \mp 0.047`` | none; without ``\ell``-dependence ``\det DK = 1`` exactly |
| [non-separable stages](@ref projection-nonseparable), ``d = 1, 2`` | the velocity form agrees with the implicit midpoint rule for ``H^N_s`` to ``4.4 \cdot 10^{-16}``; ``\max \lvert DR^\top \mathbb{J} DR - \mathbb{J}\rvert \leq 1.1 \cdot 10^{-15}`` for two stages; ``\lvert\Delta P\rvert \leq 4.4 \cdot 10^{-16}``; error of the determinant ``\leq 4.4 \cdot 10^{-15}``, of the divergence ``\leq 8.9 \cdot 10^{-16}``; at ``x = \pi``, ``\pm v``: ``\det DN - 1 = \mp 0.120`` for a separable, ``-0.144`` at both for a non-separable ``\phi_s`` | trapezoidal rule: defect ``7.5 \cdot 10^{-4}`` and ``3.6 \cdot 10^{-4}`` with ``\Delta P = 0``; no centroid shift: ``\lvert\Delta P\rvert = 0.29`` and ``1.7``; first-order determinant: ``\geq 0.315``; sign of ``\tfrac{1}{2} \mathbb{J} K`` reversed: ``\geq 0.475`` |
| [charge-conserving reweighting](@ref projection-charge), ``d = 1``, 400 markers, 16 shape functions | after a kick and a drift, Newton converges in 5 iterations; grid charge to ``7.1 \cdot 10^{-16}``; total weight and ``P`` exact; ``\min w' > 0``; gradient of the relative entropy orthogonal to the null space of ``A`` to ``9.1 \cdot 10^{-16}``, and 50 feasible perturbations increase it | ``w^*`` alone: charge error ``0.29``; additive correction: exact charge, ``\min w' = -0.015``; markers compressed to ``[0, \pi)``: no solution, reported |
| [BUG integrator](@ref projection-bug), ``64 \times 256`` grid, ``T = 3`` | free streaming from rank four: numerical rank 3, error against the full grid ``7.6 \cdot 10^{-7}``, ``1.9 \cdot 10^{-7}`` at half the step; full grid against characteristics for ``E = -0.5 \sin x``: ``2.7 \cdot 10^{-5}`` (``3.6 \cdot 10^{-4}`` on half the velocity points); ``\vartheta = 10^{-4}, 10^{-6}``: rank ``13, 31``, error ``3.6 \cdot 10^{-3}, 2.8 \cdot 10^{-5}`` below ``n_{\mathrm{steps}} \vartheta``; mass error ``6.8 \cdot 10^{-5}, 3.2 \cdot 10^{-7}`` | free streaming from rank one: stays at rank one, error ``0.28``; fixed rank one with the field: error ``0.16`` |

The [verification](@ref projection-verification) lists the test file and the script of every
row. All checks pass. The figures of three of the statements — the image of a hole under one
weight-independent step, the fold of the implicit weight for a Maxwellian, and the ``\pm v``
coupling at ``x = \pi`` — are in [one map for all markers](@ref one-map-for-all-markers),
[the velocity–momentum fold](@ref velocity-momentum-fold) and
[the non-separable stages](@ref projection-nonseparable).
