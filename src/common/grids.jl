"""
    PeriodicGrid(L, Ng)

A periodic grid on `[0, L)` with `Ng` nodes `x_i = (i - 1) Δx`, `Δx = L / Ng`, and the linear
B-spline (hat) shape functions `S_i`. On the cell `[x_i, x_{i+1})` the only non-zero shape
functions are `S_i(x) = 1 − (x − x_i)/Δx` and `S_{i+1}(x) = (x − x_i)/Δx`, with `i + 1` taken
modulo `Ng`; they form a partition of unity.
"""
struct PeriodicGrid{T <: Real}
    L::T
    Ng::Int
    Δx::T
end

PeriodicGrid(L::Real, Ng::Integer) = PeriodicGrid(float(L), Int(Ng), float(L) / Ng)

"""
    cellof(g::PeriodicGrid, x)

The 1-based index of the cell `[x_i, x_{i+1})` that holds the position `x` (taken modulo `L`).
"""
cellof(g::PeriodicGrid, x::Real) = min(floor(Int, mod(x, g.L) / g.Δx), g.Ng - 1) + 1

"""
    shape(g::PeriodicGrid, x) -> (i, j, Si, Sj)

The two nodes `i` and `j = i + 1 (mod Ng)` of the cell of `x` and the values `S_i(x)`,
`S_j(x)` of their hat functions; all other shape functions vanish at `x`.
"""
function shape(g::PeriodicGrid, x::Real)
    i = cellof(g, x)
    t = mod(x, g.L) / g.Δx - (i - 1)
    i, mod1(i + 1, g.Ng), 1 - t, t
end

"""
    shape_matrix(g::PeriodicGrid, x)

The `N × Ng` matrix with entries `S_i(x_a)`.
"""
function shape_matrix(g::PeriodicGrid, x::AbstractVector)
    S = zeros(float(eltype(x)), length(x), g.Ng)
    for (a, xa) in pairs(x)
        i, j, Si, Sj = shape(g, xa)
        S[a, i] += Si
        S[a, j] += Sj
    end
    S
end

"""
    grid_moment(g::PeriodicGrid, x, v, w, k)

The grid moments `M_i^(k) = Σ_a w_a S_i(x_a) v_a^k`, `i = 1, …, Ng`, of 1D-1V markers. The grid
charge is `q M^(0)`, the grid current `q M^(1)`, and `(m/2) M^(2)` the kinetic energy that each
shape function collects.
"""
function grid_moment(
        g::PeriodicGrid, x::AbstractVector, v::AbstractVector, w::AbstractVector,
        k::Integer)
    M = zeros(promote_type(eltype(v), eltype(w)), g.Ng)
    for a in eachindex(x, v, w)
        i, j, Si, Sj = shape(g, x[a])
        c = w[a] * v[a]^k
        M[i] += c * Si
        M[j] += c * Sj
    end
    M
end

grid_moment(g::PeriodicGrid, mk::Markers, k::Integer) = grid_moment(g, mk.x, mk.v, mk.w, k)

export PeriodicGrid, cellof, shape, shape_matrix, grid_moment
