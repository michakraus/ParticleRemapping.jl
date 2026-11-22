# Construction 2: canonical transport with consistent weights (1D-1V).
#
# A canonical map ψ of (x, p), the same for every marker, moves the markers; the weights then
# follow from a current-quantile matching, a gauge step at fixed (X, P) that makes the markers
# represent the density f again.

"""
    simpson_rule(g, a, b; n = 4000)

The composite Simpson rule for `∫_a^b g` with `n` (even) subintervals.
"""
function simpson_rule(g, a::Real, b::Real; n::Integer = 4000)
    iseven(n) || throw(ArgumentError("Simpson's rule needs an even number of subintervals"))
    h = (b - a) / n
    s = g(a) + g(b)
    for k in 1:(n - 1)
        s += (isodd(k) ? 4 : 2) * g(a + k * h)
    end
    s * h / 3
end

"""
    represented_density(ñ, W, dWV; m = 1, frozen = false)

The density `f'(X, V) = ñ(X, m W V) m ∂_V(W V) W` that markers with the canonical density
`ñ(X, P)` and the weight field `W(X, V)` represent; `dWV(X, V) = ∂_V(W V)`. With
`frozen = true` it is the frozen-weight density `ñ(X, m W V) m W²`, which drops `V ∂_V W` and is
wrong as soon as `W` depends on `V`.
"""
function represented_density(ñ, W, dWV; m::Real = 1, frozen::Bool = false)
    frozen && return (X, V) -> ñ(X, m * W(X, V) * V) * m * W(X, V)^2
    (X, V) -> ñ(X, m * W(X, V) * V) * m * dWV(X, V) * W(X, V)
end

"""
    cumulative_current_error(f′, ñ, W, X, V; m = 1, n = 4000)

The relative difference between the cumulative current `∫_0^V f′(X, V′) V′ dV′` of a represented
density `f′` and the cumulative momentum `(1/m) ∫_0^{m W V} ñ(X, P′) P′ dP′`, both by
[`simpson_rule`](@ref). It vanishes for the exact [`represented_density`](@ref), the identity in
the proof of the matching proposition.
"""
function cumulative_current_error(
        fp, ñ, W, X::Real, V::Real; m::Real = 1, n::Integer = 4000)
    lhs = simpson_rule(u -> fp(X, u) * u, 0.0, V; n)
    rhs = simpson_rule(P -> ñ(X, P) * P, 0.0, m * W(X, V) * V; n) / m
    abs(lhs - rhs) / abs(rhs)
end

# one classical Runge–Kutta step of ẋ = a(x), ξ̇ = a'(x) ξ
function _rk4_lift(a, da, x, ξ, dt)
    F(x, ξ) = (a(x), da(x) * ξ)
    k1 = F(x, ξ)
    k2 = F(x + dt / 2 * k1[1], ξ + dt / 2 * k1[2])
    k3 = F(x + dt / 2 * k2[1], ξ + dt / 2 * k2[2])
    k4 = F(x + dt * k3[1], ξ + dt * k3[2])
    x + dt / 6 * (k1[1] + 2k2[1] + 2k3[1] + k4[1]),
    ξ + dt / 6 * (k1[2] + 2k2[2] + 2k3[2] + k4[2])
end

"""
    cotangent_lift_flow(a, da, x0, t; dt = 1e-3) -> (φ_t(x0), φ_t'(x0))

The flow `φ_t` of `ẋ = a(x)` and its derivative, by the classical Runge–Kutta method on
`(x, ξ)` with `ξ̇ = a'(x) ξ`, `ξ(0) = 1`, in `round(t/dt)` steps; `da = a'`.
"""
function cotangent_lift_flow(a, da, x0::Real, t::Real; dt::Real = 1e-3)
    x, ξ = float(x0), one(float(x0))
    for _ in 1:round(Int, t / dt)
        x, ξ = _rk4_lift(a, da, x, ξ, dt)
    end
    x, ξ
end

"""
    cotangent_lift(a, da, x, p, t; L = 2π, dt = 1e-2) -> (X, P)

The cotangent lift `X = φ_t(x)`, `P = p / φ_t'(x)` of the flow of `ẋ = a(x)` on the periodic
interval `[0, L)`, applied to every marker ([`cotangent_lift_flow`](@ref)). It is canonical for a
fixed `a`, the time-`t` flow of the Hamiltonian `h = p a(x)`.
"""
function cotangent_lift(a, da, x::AbstractVector, p::AbstractVector, t::Real; L::Real = 2π,
        dt::Real = 1e-2)
    X, P = similar(x, float(eltype(x))), similar(p, float(eltype(p)))
    for i in eachindex(x, p)
        φ, ξ = cotangent_lift_flow(a, da, x[i], t; dt)
        X[i] = mod(φ, L)
        P[i] = p[i] / ξ
    end
    X, P
end

"""
    lift_compatibility_residual(J, x, X, dφ)

The relative compatibility residual `R(x) = J(x) / (φ'(x)² J(φ(x))) − 1` of a cotangent lift with
`X = φ(x)` and `dφ = φ'(x)`, for a density with current profile proportional to `J`: zero exactly
when the lift keeps the sign-resolved current profile.
"""
lift_compatibility_residual(J, x::Real, X::Real, dφ::Real) = J(x) / (dφ^2 * J(X)) - 1

"""
    maxwellian_current(V)

The cumulative current `G(V) = ∫_0^{|V|} v M(v) dv = (1 − e^{−V²/2}) / √(2π)` of the Maxwellian
`M(v) = e^{−v²/2}/√(2π)`, with `G(∞) = 1/√(2π)`.
"""
maxwellian_current(V::Real) = (1 - exp(-V^2 / 2)) / sqrt(2π)

"""
    maxwellian_current_inverse(y)

The inverse `G⁻¹(y) = √(−2 ln(1 − √(2π) y))` of [`maxwellian_current`](@ref) on `[0, ∞)`.
"""
maxwellian_current_inverse(y::Real) = sqrt(-2log(1 - sqrt(2π) * y))

"""
    SeparableMaxwellianTarget(Jbar)
    SeparableMaxwellianTarget(grid::PeriodicGrid, α)

The density `f̂(x, v) = J(x) M(v)` to be preserved by the current-quantile matching, through the
cell means `J̄_c` of `J`; for `J(x) = 1 + α cos x` on `[0, 2π)` they are
`J̄_c = 1 + α (sin(c Δx) − sin((c − 1) Δx)) / Δx`. Its cell-averaged cumulative current is
`F̂_c(V) = J̄_c G(V)` on each branch ([`maxwellian_current`](@ref)).
"""
struct SeparableMaxwellianTarget{T <: Real}
    Jbar::Vector{T}
end

function SeparableMaxwellianTarget(grid::PeriodicGrid, α::Real)
    Δx = grid.Δx
    SeparableMaxwellianTarget([1 + α * (sin(c * Δx) - sin((c - 1) * Δx)) / Δx
                               for c in 1:(grid.Ng)])
end

"""
    cell_cumulative_current(target, c, V)

The cell-averaged cumulative current `F̂_c(V)` of the density `target` in cell `c`.
"""
function cell_cumulative_current(t::SeparableMaxwellianTarget, c::Integer, V::Real)
    t.Jbar[c] *
    maxwellian_current(V)
end

"""
    cell_current_total(target, c, s)

The total `F̂_c(s ∞)` of the cumulative current of the density `target` in cell `c` on the
branch of sign `s`.
"""
function cell_current_total(t::SeparableMaxwellianTarget, c::Integer, s::Integer)
    t.Jbar[c] / sqrt(2π)
end

"""
    cell_current_inverse(target, c, s, y)

The speed `|V|` on the branch of sign `s` of cell `c` at which the cumulative current of the
density `target` takes the value `y`, `F̂_c^{−1}(y)`.
"""
function cell_current_inverse(t::SeparableMaxwellianTarget, c::Integer, s::Integer, y::Real)
    maxwellian_current_inverse(y / t.Jbar[c])
end

"""
    DensityTarget(grid::PeriodicGrid, f; vmax, nx = 8, nv = 400)

The density `f(x, v)` to be preserved by the current-quantile matching, for any non-negative
function `f`, periodic in `x` with the period `grid.L`. The cell-averaged cumulative currents
`F̂_c(sV) = (1/Δx) ∫_{C_c} ∫_0^{V} f(x, s v') v' dv' dx` are tabulated on `nv + 1` equidistant
speeds in `[0, vmax]` (the trapezoidal rule in `v`, the midpoint rule on `nx` points per cell in
`x`); `f` is taken to vanish for `|v| > vmax`. The inverse interpolates the table linearly and
returns the smallest speed at which the current reaches its value.
"""
struct DensityTarget{T <: Real}
    V::Vector{T}
    F::Array{T, 3}
end

function DensityTarget(
        grid::PeriodicGrid, f; vmax::Real, nx::Integer = 8, nv::Integer = 400)
    V = collect(range(0, float(vmax); length = nv + 1))
    h = V[2] - V[1]
    F = zeros(eltype(V), nv + 1, grid.Ng, 2)
    for c in 1:(grid.Ng), s in (1, -1)

        xs = ((c - 1) .+ ((1:nx) .- 0.5) ./ nx) .* grid.Δx
        j = [sum(f(x, s * v) for x in xs) / nx * v for v in V]
        F[2:end, c, _branch(s)] = cumsum((j[1:(end - 1)] .+ j[2:end]) ./ 2 .* h)
    end
    DensityTarget(V, F)
end

_branch(s::Integer) = s > 0 ? 1 : 2

function cell_cumulative_current(t::DensityTarget, c::Integer, V::Real)
    s, aV = V >= 0 ? 1 : -1, abs(V)
    aV >= t.V[end] && return t.F[end, c, _branch(s)]
    i = searchsortedlast(t.V, aV)
    θ = (aV - t.V[i]) / (t.V[i + 1] - t.V[i])
    (1 - θ) * t.F[i, c, _branch(s)] + θ * t.F[i + 1, c, _branch(s)]
end

cell_current_total(t::DensityTarget, c::Integer, s::Integer) = t.F[end, c, _branch(s)]

function cell_current_inverse(t::DensityTarget, c::Integer, s::Integer, y::Real)
    F = view(t.F, :, c, _branch(s))
    i = clamp(searchsortedfirst(F, y), 2, length(F))
    θ = (y - F[i - 1]) / (F[i] - F[i - 1])
    (1 - θ) * t.V[i - 1] + θ * t.V[i]
end

"""
    branch_quantile_matching(aP, total, inverse; m = 1, Δx = 1)

Steps 3–5 of the algorithm *current-quantile matching* on one branch of one cell, for the sorted
speeds of momentum `aP = |P_(1)| ≤ … ≤ |P_(n)|`, the total `total = F̂_c(s∞)` of the cumulative
current and its inverse `inverse(y) = F̂_c^{−1}(y)`. Returns `(; Π, R, Πc, V, W)`: the midpoint
cumulative momenta `Π_(r) = (Σ_{j<r} |P_(j)| + |P_(r)|/2) / (m Δx)`, the compatibility residual
`R = Σ_j |P_(j)| / (m Δx F̂_c(s∞)) − 1`, the closed quantiles `Πc = Π/(1 + R)`, the speeds
`|V_(r)| = F̂_c^{−1}(Πc_(r))` and the weights `W_(r) = |P_(r)| / (m |V_(r)|)`.
"""
function branch_quantile_matching(aP::AbstractVector, total::Real, inverse; m::Real = 1,
        Δx::Real = 1)
    issorted(aP) || throw(ArgumentError("the momenta |P| must be sorted"))
    Π = (cumsum(aP) .- aP ./ 2) ./ (m * Δx)
    sumP = sum(aP) / (m * Δx)
    Πc = Π .* (total / sumP)
    V = inverse.(Πc)
    (; Π, R = sumP / total - 1, Πc, V, W = aP ./ (m .* V))
end

"""
    current_quantile_match(grid, X, P, target; m = 1) -> (V, W, R)

Steps 2–5 of the algorithm *current-quantile matching*: for each cell `c` of `grid` and each sign
`s`, sort the markers with `X_a ∈ C_c` and `sgn P_a = s` by `|P_a|` and apply
[`branch_quantile_matching`](@ref) with the cumulative current of `target`. Returns the new
velocities `V`, the new weights `W` (with `m W V = P`) and the residuals `R[c, 1]` (`s = +1`) and
`R[c, 2]` (`s = −1`), `NaN` for an empty branch.
"""
function current_quantile_match(grid::PeriodicGrid, X::AbstractVector, P::AbstractVector,
        target; m::Real = 1)
    any(iszero, P) && throw(ArgumentError("a marker with P = 0 lies on no branch"))
    V, W = similar(P, float(eltype(P))), similar(P, float(eltype(P)))
    R = fill(float(eltype(P))(NaN), grid.Ng, 2)
    groups = Dict{Tuple{Int, Int}, Vector{Int}}()
    for a in eachindex(X, P)
        push!(get!(groups, (cellof(grid, X[a]), Int(sign(P[a]))), Int[]), a)
    end
    for ((c, s), idx) in groups
        o = idx[sortperm(abs.(P[idx]))]
        b = branch_quantile_matching(abs.(P[o]), cell_current_total(target, c, s),
            y -> cell_current_inverse(target, c, s, y); m, Δx = grid.Δx)
        V[o] = s .* b.V
        W[o] = b.W
        R[c, s > 0 ? 1 : 2] = b.R
    end
    V, W, R
end

"""
    CurrentQuantileMatching(grid, canonical_map, target)

The canonical transport with consistent weights: the canonical map
`canonical_map(x, p) -> (X, P)`, applied to every marker, followed by the current-quantile
matching ([`current_quantile_match`](@ref)) on the cells of `grid`, which sets the weights and
velocities at fixed `(X, P)` so that the markers represent the density `target`
(e.g. a [`SeparableMaxwellianTarget`](@ref)).
"""
struct CurrentQuantileMatching{G <: PeriodicGrid, F, T} <: AbstractRemap
    grid::G
    canonical_map::F
    target::T
end

"""
    remap(method::CurrentQuantileMatching, markers::Markers) -> RemapResult

Run the algorithm *current-quantile matching*: `p_a = m w_a v_a`, `(X, P) = ψ(x, p)`, then the
matching. The `info` holds the canonical coordinates `X`, `P` after the map, the compatibility
residuals `R` per cell and sign, and their largest modulus `max_residual`.
"""
function remap(method::CurrentQuantileMatching, mk::Markers)
    mk.x isa AbstractVector || throw(ArgumentError("the matching acts on 1D-1V markers"))
    X, P = method.canonical_map(mk.x, momentum(mk))
    V, W, R = current_quantile_match(method.grid, X, P, method.target; m = mk.m)
    RemapResult(Markers(X, V, W; m = mk.m),
        (; X, P, R, max_residual = maximum(abs, filter(!isnan, R))))
end

"""
    binned_weights(grid, X, V, W, vedges)

The sums of the weights `W` of the markers in each cell of `grid` and each velocity bin
`[vedges[j], vedges[j+1])`, an `N_g × (length(vedges) − 1)` matrix; markers outside the bins are
left out.
"""
function binned_weights(grid::PeriodicGrid, X::AbstractVector, V::AbstractVector,
        W::AbstractVector, vedges::AbstractVector)
    nb = length(vedges) - 1
    H = zeros(eltype(W), grid.Ng, nb)
    for a in eachindex(X, V, W)
        j = searchsortedlast(vedges, V[a])
        1 ≤ j ≤ nb || continue
        H[cellof(grid, X[a]), j] += W[a]
    end
    H
end

"""
    binned_weight_error(target::SeparableMaxwellianTarget, grid, X, V, W;
                        vedges = range(-4, 4; length = 25))

The relative `L¹` difference between the [`binned_weights`](@ref) and the exact bin integrals
`J̄_c Δx ∫_bin M` of the density `target` (Simpson's rule with 200 subintervals per bin).
"""
function binned_weight_error(target::SeparableMaxwellianTarget, grid::PeriodicGrid,
        X::AbstractVector, V::AbstractVector, W::AbstractVector;
        vedges::AbstractVector = range(-4, 4; length = 25))
    H = binned_weights(grid, X, V, W, vedges)
    Mbin = [simpson_rule(u -> exp(-u^2 / 2) / sqrt(2π), vedges[j], vedges[j + 1]; n = 200)
            for j in 1:(length(vedges) - 1)]
    E = [target.Jbar[c] * grid.Δx * Mbin[j] for c in 1:(grid.Ng), j in eachindex(Mbin)]
    sum(abs.(H .- E)) / sum(E)
end
