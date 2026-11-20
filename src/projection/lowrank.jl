# Dynamical low-rank approximation f(t, x, v) ≈ Σ_ij X_i(t, x) S_ij(t) V_j(t, v) of the Vlasov
# equation with a fixed field, ∂_t f + v ∂_x f + E(x) ∂_v f = 0, d = 1, integrated with the
# rank-adaptive basis-update and Galerkin (BUG) integrator. Spectral differentiation in x on
# [0, L), fourth-order central differences in v on [−Vd, Vd] with f = 0 outside; X is orthonormal
# in hx Σ, V in hv Σ, so that ‖f‖_{L²} = ‖S‖_F.

"""
    VlasovGrid(nx, nv, L, Vd, E)

The grid of `nx` equispaced points of `[0, L)` and `nv` points of `[−Vd, Vd]` (including both
ends), with the spectral differentiation matrix `Dx`, the fourth-order central differences `Dv`
(zero values outside) and the field `E` at the points of `x`.
"""
struct VlasovGrid{T <: Real}
    x::Vector{T}
    v::Vector{T}
    hx::T
    hv::T
    Dx::Matrix{T}
    Dv::Matrix{T}
    E::Vector{T}
end

# The spectral differentiation matrix on n (even) equispaced points of a period L.
function _fourier_diff(n, L)
    h = 2π / n
    [i == j ? 0.0 : 0.5 * (-1)^(i - j) * cot((i - j) * h / 2) * 2π / L
     for i in 1:n, j in 1:n]
end

# Fourth-order central differences on n points of spacing h, with zero values outside.
function _central_diff4(n, h)
    D = zeros(n, n)
    for i in 1:n, (o, c) in ((-2, 1), (-1, -8), (1, 8), (2, -1))

        1 <= i + o <= n && (D[i, i + o] = c / (12h))
    end
    D
end

function VlasovGrid(nx::Integer, nv::Integer, L::Real, Vd::Real, E)
    hx, hv = L / nx, 2Vd / (nv - 1)
    x = collect((0:(nx - 1)) .* hx)
    v = collect(range(-Vd, Vd; length = nv))
    VlasovGrid(x, v, hx, hv, _fourier_diff(nx, L), _central_diff4(nv, hv), float.(E.(x)))
end

# The right-hand side of the Vlasov equation for the full grid function F[i, j] = f(x_i, v_j).
_vlasov_rhs(g::VlasovGrid, F) = -(g.Dx * F) .* g.v' .- g.E .* (F * g.Dv')

function _rk4(rhs, y, dt)
    k1 = rhs(y)
    k2 = rhs(y .+ dt / 2 .* k1)
    k3 = rhs(y .+ dt / 2 .* k2)
    k4 = rhs(y .+ dt .* k3)
    y .+ dt / 6 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
end

"""
    LowRankDensity(X, S, V)

The low-rank grid function `X S Vᵀ` with orthonormal columns of `X` (in `hx Σ`) and of `V` (in
`hv Σ`).
"""
struct LowRankDensity{T <: Real}
    X::Matrix{T}
    S::Matrix{T}
    V::Matrix{T}
end

"""
    lowrank_rank(f::LowRankDensity)

The rank `r` of the factorisation: the size of `S`.
"""
lowrank_rank(f::LowRankDensity) = size(f.S, 1)

"""
    lowrank_full(f::LowRankDensity)

The full grid function `X S Vᵀ`.
"""
lowrank_full(f::LowRankDensity) = f.X * f.S * f.V'

"""
    lowrank_mass(g::VlasovGrid, f::LowRankDensity)

The mass `hx hv Σ_ij f(x_i, v_j)` of the low-rank density.
"""
function lowrank_mass(g::VlasovGrid, f::LowRankDensity)
    g.hx * g.hv * sum(f.X * f.S * (f.V' * ones(length(g.v))))
end

# An orthonormal basis of the columns of M in the inner product h Σ, at most size(M, 1) columns.
_orth(h, M) = Matrix(qr(sqrt(h) .* M).Q)[:, 1:min(size(M)...)] ./ sqrt(h)

"""
    numerical_rank(σ, tol)

The smallest rank `r` whose discarded singular values `σ[r+1:end]` have the norm `tol` or less.
"""
function numerical_rank(σ, tol)
    tail = sqrt.(reverse(cumsum(reverse(σ .^ 2))))     # tail[r] = ‖σ[r:end]‖
    something(
        findfirst(r -> (r == length(σ) ? zero(tol) : tail[r + 1]) <= tol, eachindex(σ)),
        length(σ))
end

"""
    lowrank_factorise(g::VlasovGrid, F; tol, rmin = 1, rmax = typemax(Int),
        rng = Random.default_rng())

The low-rank factorisation of the grid function `F`, truncated at `tol ‖F‖_{L²}`, with at least
`rmin` and at most `rmax` columns. Columns beyond the numerical rank come from random vectors with
zero coefficients: the projected dynamics cannot leave the rank it starts at.
"""
function lowrank_factorise(g::VlasovGrid, F; tol, rmin = 1, rmax = typemax(Int),
        rng::AbstractRNG = Random.default_rng())
    U, σ, W = svd(sqrt(g.hx * g.hv) .* F)
    r = clamp(numerical_rank(σ, tol * norm(σ)), 1, rmax)
    X, V = U[:, 1:r] ./ sqrt(g.hx), W[:, 1:r] ./ sqrt(g.hv)
    S = Matrix(Diagonal(σ[1:r]))
    if r < rmin
        X = _orth(g.hx, hcat(X, randn(rng, length(g.x), rmin - r)))
        V = _orth(g.hv, hcat(V, randn(rng, length(g.v), rmin - r)))
        S = [S zeros(r, rmin - r); zeros(rmin - r, rmin)]
    end
    LowRankDensity(X, S, V)
end

# The coefficient matrices ⟨V_j, v V_l⟩, ⟨V_j, ∂_v V_l⟩ and ⟨X_i, ∂_x X_k⟩, ⟨X_i, E X_k⟩.
_coefficients_v(g, V) = (g.hv .* (V' * (g.v .* V)), g.hv .* (V' * (g.Dv * V)))
_coefficients_x(g, X) = (g.hx .* (X' * (g.Dx * X)), g.hx .* (X' * (g.E .* X)))

"""
    bug_step(g::VlasovGrid, f::LowRankDensity, dt; tol, rmin = 1, rmax = typemax(Int))
        -> (f′, numerical rank)

One step of the rank-adaptive BUG integrator, each substep by the classical Runge–Kutta method:
the K-step `∂_t K = −∂_x K C¹ᵀ − E K C²ᵀ` from `K = X S` with `V` fixed, the L-step
`∂_t L = −v L D¹ᵀ − ∂_v L D²ᵀ` from `L = V Sᵀ` with `X` fixed, the new bases spanning `[K, X]` and
`[L, V]`, the Galerkin S-step `∂_t S = −D¹ S C¹ᵀ − D² S C²ᵀ` in the doubled bases, and a truncated
singular value decomposition of `S` to the smallest rank whose discarded singular values have the
norm `tol ‖S‖_F` or less, kept between `rmin` and `rmax`.
"""
function bug_step(g::VlasovGrid, f::LowRankDensity, dt; tol, rmin = 1, rmax = typemax(Int))
    X0, S0, V0 = f.X, f.S, f.V
    C1, C2 = _coefficients_v(g, V0)
    K = _rk4(K -> -(g.Dx * K) * C1' .- (g.E .* K) * C2', X0 * S0, dt)
    D1, D2 = _coefficients_x(g, X0)
    Lm = _rk4(L -> -(g.v .* L) * D1' .- (g.Dv * L) * D2', V0 * S0', dt)
    Xh, Vh = _orth(g.hx, hcat(K, X0)), _orth(g.hv, hcat(Lm, V0))
    S = (g.hx .* (Xh' * X0)) * S0 * (g.hv .* (V0' * Vh))
    C1, C2 = _coefficients_v(g, Vh)
    D1, D2 = _coefficients_x(g, Xh)
    S = _rk4(S -> -D1 * S * C1' .- D2 * S * C2', S, dt)
    P, σ, Q = svd(S)
    rn = numerical_rank(σ, tol * norm(σ))
    r = clamp(rn, rmin, min(rmax, length(σ)))
    LowRankDensity(Xh * P[:, 1:r], Matrix(Diagonal(σ[1:r])), Vh * Q[:, 1:r]), rn
end

"""
    dlra_solve(g::VlasovGrid, F0, T, dt; tol, rmin = 1, rmax = typemax(Int),
        rng = Random.default_rng()) -> (f, ranks)

The low-rank solution at time `T` from the grid function `F0` ([`lowrank_factorise`](@ref), then
[`bug_step`](@ref) with the step `dt`), and the numerical rank after each step.
"""
function dlra_solve(g::VlasovGrid, F0, T, dt; tol, rmin = 1, rmax = typemax(Int),
        rng::AbstractRNG = Random.default_rng())
    f = lowrank_factorise(g, F0; tol, rmin, rmax, rng)
    ranks = Int[]
    for _ in 1:round(Int, T / dt)
        f, rn = bug_step(g, f, dt; tol, rmin, rmax)
        push!(ranks, rn)
    end
    f, ranks
end

"""
    vlasov_full_solve(g::VlasovGrid, F0, T, dt)

The full-grid solution at time `T` with the same discretisation and the classical Runge–Kutta
method with the step `dt`, the reference of the low-rank solver.
"""
function vlasov_full_solve(g::VlasovGrid, F0, T, dt)
    F = copy(F0)
    for _ in 1:round(Int, T / dt)
        F = _rk4(F -> _vlasov_rhs(g, F), F, dt)
    end
    F
end

# Trigonometric interpolation weights of the periodic grid g.x at the point x.
function _trig_weights(g::VlasovGrid, x)
    n, L = length(g.x), length(g.x) * g.hx
    map(g.x) do xj
        y = 2π * (x - xj) / L
        s = sin(y / 2)
        abs(s) < 1e-14 ? 1.0 : sin(n * y / 2) / (n * tan(y / 2))
    end
end

# Cubic Lagrange weights of the velocity grid at the point v, zero outside [−Vd, Vd].
function _cubic_weights(g::VlasovGrid, v)
    a = zeros(length(g.v))
    t = (v - g.v[1]) / g.hv
    (0 <= t <= length(g.v) - 1) || return a
    j0 = clamp(floor(Int, t) - 1, 0, length(g.v) - 4)
    for j in j0:(j0 + 3)
        a[j + 1] = prod((t - i) / (j - i) for i in j0:(j0 + 3) if i != j)
    end
    a
end

"""
    lowrank_evaluate(g::VlasovGrid, f::LowRankDensity, x, v)

The values of the low-rank density at the points `(x_k, v_k)`: trigonometric interpolation in `x`,
cubic Lagrange interpolation in `v`, zero outside `[−Vd, Vd]`.
"""
function lowrank_evaluate(g::VlasovGrid, f::LowRankDensity, x, v)
    [dot(_trig_weights(g, x[k])' * f.X, f.S * (f.V' * _cubic_weights(g, v[k])))
     for k in eachindex(x)]
end

"""
    vlasov_characteristics(E, x, v, T; h = 1e-3)

The flow of `ẋ = v`, `v̇ = E(x)` over the time `T` from the points `(x, v)` (vectors), by the
classical Runge–Kutta method with a step of about `h`; `T < 0` runs backwards.
"""
function vlasov_characteristics(E, x, v, T; h = 1e-3)
    n = max(1, round(Int, abs(T) / h))
    dt = T / n
    x, v = copy(x), copy(v)
    for _ in 1:n
        k1x, k1v = v, E.(x)
        k2x, k2v = v .+ dt / 2 .* k1v, E.(x .+ dt / 2 .* k1x)
        k3x, k3v = v .+ dt / 2 .* k2v, E.(x .+ dt / 2 .* k2x)
        k4x, k4v = v .+ dt .* k3v, E.(x .+ dt .* k3x)
        x = x .+ dt / 6 .* (k1x .+ 2 .* k2x .+ 2 .* k3x .+ k4x)
        v = v .+ dt / 6 .* (k1v .+ 2 .* k2v .+ 2 .* k3v .+ k4v)
    end
    x, v
end
