# The de-clustering relaxation of the target points and the kernel-density measure of clustering.

# The periodic difference of two positions, in [−L/2, L/2).
_periodic_difference(x, y, L) = mod(x - y + L / 2, L) - L / 2

"""
    declustering_energy(x, v; Lx, σx, σv, v_wall = 2, k_wall = 2.5)

The smoothed particle-count energy
`E_rep = ½ Σ_{a≠b} K_σ(x_a − x_b, v_a − v_b) + Σ_a U_wall(v_a)` with the periodic-in-`x`
Gaussian kernel `K_σ = exp(−½ (Δx/σx)² − ½ (Δv/σv)²)` and the soft wall
`U_wall(v) = (k_wall/2) (|v| − v_wall)²` for `|v| > v_wall`.
"""
function declustering_energy(
        x::AbstractVector, v::AbstractVector; Lx::Real, σx::Real, σv::Real,
        v_wall::Real = 2, k_wall::Real = 2.5)
    E = zero(promote_type(eltype(x), eltype(v)))
    for a in eachindex(x, v)
        for b in eachindex(x, v)
            a == b && continue
            dx = _periodic_difference(x[a], x[b], Lx)
            E += exp(-(dx / σx)^2 / 2 - ((v[a] - v[b]) / σv)^2 / 2) / 2
        end
        abs(v[a]) > v_wall && (E += k_wall / 2 * (abs(v[a]) - v_wall)^2)
    end
    E
end

"""
    declustering_relaxation(x, v; Lx, σx, σv, step = 0.006, iterations = 10, v_wall = 2, k_wall = 2.5)
        -> (x, v, history)

The target points `(x̄, v̄)` after `iterations` damped steps `(x, v) ← (x, v) − step ∇E_rep` of
gradient relaxation on the [`declustering_energy`](@ref), and the `history` of the points
`(x, v)` after each step (the first entry is the start). Positions are not wrapped.
"""
function declustering_relaxation(
        x::AbstractVector, v::AbstractVector; Lx::Real, σx::Real, σv::Real,
        step::Real = 0.006, iterations::Integer = 10, v_wall::Real = 2, k_wall::Real = 2.5)
    N = length(x)
    x̄, v̄ = copy(x), copy(v)
    history = [(copy(x̄), copy(v̄))]
    for _ in 1:iterations
        dx_r, dv_r = zeros(eltype(x̄), N), zeros(eltype(v̄), N)
        for a in 1:N
            for b in 1:N
                a == b && continue
                dx = _periodic_difference(x̄[a], x̄[b], Lx)
                dv = v̄[a] - v̄[b]
                ker = exp(-0.5 * (dx / σx)^2 - 0.5 * (dv / σv)^2)
                dx_r[a] += (dx / σx^2) * ker
                dv_r[a] += (dv / σv^2) * ker
            end
            v̄[a] > v_wall && (dv_r[a] -= k_wall * (v̄[a] - v_wall))
            v̄[a] < -v_wall && (dv_r[a] -= k_wall * (v̄[a] + v_wall))
        end
        x̄ .+= step .* dx_r
        v̄ .+= step .* dv_r
        push!(history, (copy(x̄), copy(v̄)))
    end
    (x = x̄, v = v̄, history = history)
end

"""
    kernel_density_cv(x, v; Lx, σx, σv, nx = 20, nv = 20, vlims = (-2.4, 2.4))

The coefficient of variation (standard deviation over mean) of the Gaussian kernel density
`Σ_a exp(−½ (Δx/σx)² − ½ (Δv/σv)²)` of the points on an `nx × nv` grid of
`[0, Lx) × [vlims[1], vlims[2]]`: a measure of the clustering of the points.
"""
function kernel_density_cv(
        x::AbstractVector, v::AbstractVector; Lx::Real, σx::Real, σv::Real,
        nx::Integer = 20, nv::Integer = 20, vlims = (-2.4, 2.4))
    gx = range(zero(Lx), Lx; length = nx + 1)[1:nx]
    gv = range(vlims[1], vlims[2]; length = nv)
    vals = [begin
                s = zero(promote_type(eltype(x), eltype(v)))
                for a in eachindex(x, v)
                    dx = mod(xi - x[a] + 0.5 * Lx, Lx) - 0.5 * Lx
                    s += exp(-0.5 * (dx / σx)^2 - 0.5 * ((vj - v[a]) / σv)^2)
                end
                s
            end
            for xi in gx for vj in gv]
    μ = sum(vals) / length(vals)
    sqrt(sum(abs2, vals .- μ) / length(vals)) / μ
end
