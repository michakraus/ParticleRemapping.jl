# [Target weights from a dynamical low-rank density](@id projection-lowrank)

```@meta
CurrentModule = ParticleRemapping.Projection
```

The target weights need ``f`` at the targets with a small relative error
([reweighting](@ref projection-final-reweighting)). In a simulation, ``f`` is known only through
the markers or through a second representation that runs alongside them. One such representation
is a dynamical low-rank approximation [KochLubich:2007](@cite) of the Vlasov equation, as in
[EinkemmerLubich:2018](@citet):

```math
f(t, x, v) \approx f_r(t, x, v) = \sum_{i, j = 1}^{r} X_i(t, x) \, S_{ij}(t) \, V_j(t, v) ,
```

with orthonormal ``\{X_i\}`` in ``L^2(\mathbb{T}^d)`` and ``\{V_j\}`` in ``L^2`` of the velocity
space. The target weights become ``w^*_a = c \, f_r(y^*_a)``. For ``d = 1`` and
``\partial_t f + v \, \partial_x f + E(x) \, \partial_v f = 0``, with
``C^1_{jl} = \langle V_j, v V_l \rangle``, ``C^2_{jl} = \langle V_j, \partial_v V_l \rangle``,
``D^1_{ik} = \langle X_i, \partial_x X_k \rangle`` and ``D^2_{ik} = \langle X_i, E X_k \rangle``,
one step of the rank-adaptive basis-update and Galerkin integrator [CerutiKuschLubich:2022](@cite)
solves

```math
\partial_t K = -\partial_x K \, C^{1\top} - E \, K \, C^{2\top} , \qquad
\partial_t L = -v \, L \, D^{1\top} - \partial_v L \, D^{2\top} , \qquad
\partial_t S = -D^1 S \, C^{1\top} - D^2 S \, C^{2\top} .
```

The ``K``-step starts from ``K = X S`` with ``V`` fixed, and the ``L``-step from
``L = V S^\top`` with ``X`` fixed. The new bases span ``[K, X]`` and ``[L, V]``, the ``S``-step
is the Galerkin system in these doubled bases, and a singular value decomposition of ``S``
truncates the rank to the smallest one whose discarded singular values have the norm
``\vartheta \|S\|_F`` or less.

Three properties of the low-rank approximation bear on the weights. First, the integrator
conserves neither the mass nor the sign: ``f_r`` can be negative in the tails of ``f``, and a
target weight ``w^*_a \leq 0`` makes the final reweighting meaningless. Second, the truncation
controls an absolute error in ``L^2``, not the relative error of the
[reweighting](@ref projection-final-reweighting); where ``f`` is small the relative error of
``w^*_a`` is large. Third, the projected dynamics cannot leave the rank at which it starts. For
free streaming, ``f_0 = a(x) \, b(v)`` with an even ``b`` is a fixed point of the rank-one
dynamics, because ``C^1 = \langle b, v b \rangle / \|b\|^2 = 0`` and ``D^1 = 0``, while the exact
solution has rank three for ``t > 0``. The solver therefore starts with at least four basis
functions, the extra ones with zero coefficients. The [charge-conserving
reweighting](@ref projection-charge) then removes the mass error of ``f_r`` from the grid charge.

## [The BUG integrator](@id projection-bug)

The code discretises ``x`` on ``n_x`` equispaced points of ``[0, L)`` with spectral
differentiation and ``v`` on ``n_v`` points of ``[-V_d, V_d]`` with fourth-order central
differences and ``f = 0`` outside ([`VlasovGrid`](@ref)). ``X`` is orthonormal in
``h_x \sum``, ``V`` in ``h_v \sum``, so that ``\|f\|_{L^2} = \|S\|_F``. One step of
[`bug_step`](@ref) with the step ``\Delta t`` from ``(X_0, S_0, V_0)``:

1. Compute ``C^1 = h_v V_0^\top \operatorname{diag}(v) V_0`` and
   ``C^2 = h_v V_0^\top D_v V_0``, and integrate the ``K``-step from ``K = X_0 S_0`` over
   ``\Delta t`` with one step of the classical Runge–Kutta method.
2. Compute ``D^1 = h_x X_0^\top D_x X_0`` and ``D^2 = h_x X_0^\top \operatorname{diag}(E) X_0``,
   and integrate the ``L``-step from ``L = V_0 S_0^\top`` in the same way.
3. Orthonormalise ``[K, X_0]`` to ``\hat{X}`` and ``[L, V_0]`` to ``\hat{V}``, at most twice the
   rank of columns each.
4. Project the coefficients, ``\hat{S} = (h_x \hat X^\top X_0) \, S_0 \, (h_v V_0^\top \hat V)``,
   recompute ``C^1, C^2`` in ``\hat V`` and ``D^1, D^2`` in ``\hat X``, and integrate the
   ``S``-step from ``\hat S`` with one Runge–Kutta step.
5. Compute the singular value decomposition ``\hat S = P \Sigma Q^\top`` and the numerical rank
   ``r_n``, the smallest rank whose discarded singular values have the norm
   ``\vartheta \|\Sigma\|_F`` or less ([`numerical_rank`](@ref)); keep ``r`` between the minimum
   and the maximum rank; return ``(\hat X P_{:,1:r}, \Sigma_{1:r}, \hat V Q_{:,1:r})``.

[`dlra_solve`](@ref) starts from the truncated singular value decomposition of the initial grid
function, padded with random orthonormal columns with zero coefficients up to the minimum rank
([`lowrank_factorise`](@ref)), and repeats the step. [`lowrank_evaluate`](@ref) evaluates
``f_r`` at the target points by trigonometric interpolation in ``x`` and cubic Lagrange
interpolation in ``v``; [`vlasov_full_solve`](@ref) is the full-grid reference with the same
discretisation, and [`vlasov_characteristics`](@ref) follows the characteristics.

```@example lowrank
using ParticleRemapping, CairoMakie, Random
Random.seed!(20261122)
f0(x, v) = (1 + 0.3 * cos(x)) * exp(-v^2 / 2) / sqrt(2π)
g = VlasovGrid(64, 256, 2π, 6.0, x -> -0.5 * sin(x))
F0 = f0.(g.x, g.v')
fr, ranks = dlra_solve(g, F0, 3.0, 0.01; tol = 1e-4, rmin = 4)
F = vlasov_full_solve(g, F0, 3.0, 0.01)
l2(G) = sqrt(g.hx * g.hv * sum(abs2, G))
(; largest_rank = maximum(ranks), error = l2(lowrank_full(fr) .- F) / l2(F0))
```

```@example lowrank
fig = Figure(; size = (900, 330))
ax = Axis(fig[1, 1]; title = "f_r at T = 3, ϑ = 1e-4", xlabel = "x", ylabel = "v",
    limits = (0, 2π, -4, 4))
heatmap!(ax, g.x, g.v, lowrank_full(fr); colormap = :Blues)
ax = Axis(fig[1, 2]; title = "numerical rank", xlabel = "t", ylabel = "r")
stairs!(ax, (1:length(ranks)) .* 0.01, ranks)
fig
```

The low-rank solution for ``E = -0.5 \sin x`` at ``\vartheta = 10^{-4}`` and its numerical rank
over time.

The test `test/projection/lowrank.jl` checks the solver on ``64 \times 256`` points of
``[0, 2\pi) \times [-6, 6]`` up to ``T = 3`` with the step ``0.01``. Free streaming from rank four
keeps the numerical rank 3, with the error ``7.6 \cdot 10^{-7}`` against the full grid and
``1.9 \cdot 10^{-7}`` at half the step. The full grid agrees with the characteristics for
``E = -0.5 \sin x`` to ``2.7 \cdot 10^{-5}`` (``3.6 \cdot 10^{-4}`` on half the velocity points:
fourth order in ``v``). The BUG integrator is not exact at full rank: its error has a term in the
step [CerutiKuschLubich:2022](@cite), of second order here (``2.3 \cdot 10^{-8}`` at the step
``0.01`` and full rank, ``3.8 \cdot 10^{-9}`` at half the step). At ``\vartheta = 10^{-4}`` and
``10^{-6}`` the largest ranks are 13 and 31, the errors ``3.6 \cdot 10^{-3}`` and
``2.8 \cdot 10^{-5}``, below ``n_{\mathrm{steps}} \vartheta``, and the mass errors
``6.8 \cdot 10^{-5}`` and ``3.2 \cdot 10^{-7}``. The cases that must fail: free streaming from
rank one stays at rank one with the error ``0.28``; the fixed rank one with the field has the
error ``0.16``.
