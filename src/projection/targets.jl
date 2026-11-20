# The lumped projection onto hat functions and the density-interpolation flow to the target points,
# d = 1. Both work in the scaled variables s = x / L ∈ [0, 1) (periodic) and r = v / V ∈ [−1, 1].

"""
    HatProjection(nX = 32, nV = 33)

The tensor-product basis of `nX × nV` bilinear hat functions in the scaled variables
`s = x / L ∈ [0, 1)`, periodic, and `r = v / V ∈ [−1, 1]`, clamped, with the spacings `hs = 1/nX`
and `hr = 2/(nV − 1)`. The hats are nonnegative and form a partition of unity.
"""
struct HatProjection{T <: Real}
    nX::Int
    nV::Int
    hs::T
    hr::T
end

function HatProjection(nX::Integer = 32, nV::Integer = 33)
    HatProjection(Int(nX), Int(nV), 1 / nX,
        2 / (nV - 1))
end

# The bilinear hat weights of the point (s, r): four (node i, node j, weight) triples.
function _hat_weights(hp::HatProjection, s, r)
    t = mod(s, 1) / hp.hs
    i = floor(Int, t)
    a = t - i
    u = clamp((r + 1) / hp.hr, 0, hp.nV - 1 - 1e-12)
    j = floor(Int, u)
    b = u - j
    i1, i2 = mod(i, hp.nX) + 1, mod(i + 1, hp.nX) + 1
    ((i1, j + 1, (1 - a) * (1 - b)), (i2, j + 1, a * (1 - b)),
        (i1, j + 2, (1 - a) * b), (i2, j + 2, a * b))
end

# The integral of the hat of node j in r (times the cell width in s): half at the clamped ends.
_node_volume(hp::HatProjection, j) = hp.hs * hp.hr * ((j == 1 || j == hp.nV) ? 0.5 : 1.0)

"""
    lumped_projection(hp::HatProjection, s, r, q)

The nodal values `Σ_k q_k Λ_α(s_k, r_k) / ∫Λ_α` of the lumped projection of the weighted points
`(s_k, r_k)`, an `nX × nV` matrix. With `q = w` it is the projection `f_h` of the represented
density, with `q = 1` the projection `g_h` of the marker density, both in scaled variables.
"""
function lumped_projection(hp::HatProjection, s, r, q)
    g = zeros(promote_type(eltype(q), Float64), hp.nX, hp.nV)
    for k in eachindex(s), (i, j, λ) in _hat_weights(hp, s[k], r[k])

        g[i, j] += q[k] * λ
    end
    for j in 1:(hp.nV)
        g[:, j] ./= _node_volume(hp, j)
    end
    g
end

"""
    hat_interpolate(hp::HatProjection, g, s, r)

The value at `(s, r)` of the function with the nodal values `g` on the hats `hp`.
"""
function hat_interpolate(hp::HatProjection, g, s, r)
    sum(λ * g[i, j] for (i, j, λ) in _hat_weights(hp, s, r))
end

"""
    moser_laplacian(hp::HatProjection)

The five-point Laplacian on the nodes of `hp`, periodic in `s` and with the Neumann condition in `r`
through a ghost node, as a dense matrix on the nodal values (`s` fastest).
"""
function moser_laplacian(hp::HatProjection)
    nX, nV, hs, hr = hp.nX, hp.nV, hp.hs, hp.hr
    n = nX * nV
    idx(i, j) = (j - 1) * nX + i
    A = zeros(n, n)
    for j in 1:nV, i in 1:nX

        k = idx(i, j)
        A[k, idx(mod1(i - 1, nX), j)] += 1 / hs^2
        A[k, idx(mod1(i + 1, nX), j)] += 1 / hs^2
        A[k, k] -= 2 / hs^2 + 2 / hr^2
        jm = j == 1 ? 2 : j - 1
        jp = j == nV ? nV - 1 : j + 1
        A[k, idx(i, jm)] += 1 / hr^2
        A[k, idx(i, jp)] += 1 / hr^2
    end
    A
end

"""
    moser_targets(hp::HatProjection, s, r; steps = 40, history = nothing) -> (s*, r*)

The target points `χ_1(s_k, r_k)` of the density-interpolation flow, in scaled variables: the lumped
marker density `g_h` is interpolated linearly in `t` to its mean `ḡ = N/2`, the potential solves
`Δφ = g_h − ḡ` (periodic in `s`, Neumann in `r`, by the pseudo-inverse of
[`moser_laplacian`](@ref)), and the points follow `u_t = ∇φ / ((1 − t) g_h + t ḡ)` with `steps`
steps of the classical Runge–Kutta method. A vector `history` receives the points after each step.
"""
function moser_targets(hp::HatProjection, s, r; steps = 40, history = nothing)
    nX, nV, hs, hr = hp.nX, hp.nV, hp.hs, hp.hr
    g = lumped_projection(hp, s, r, ones(length(s)))
    gbar = length(s) / 2                          # N over the scaled area 1 × 2
    φ = reshape(pinv(moser_laplacian(hp)) * vec(g .- gbar), nX, nV)
    φs = (circshift(φ, (-1, 0)) - circshift(φ, (1, 0))) / (2hs)
    φr = zeros(nX, nV)
    φr[:, 2:(end - 1)] = (φ[:, 3:end] - φ[:, 1:(end - 2)]) / (2hr)
    function u(t, s, r)
        gt = (1 - t) * hat_interpolate(hp, g, s, r) + t * gbar
        hat_interpolate(hp, φs, s, r) / gt, hat_interpolate(hp, φr, s, r) / gt
    end
    S, R = copy(s), copy(r)
    dt = 1 / steps
    history === nothing || push!(history, (copy(S), copy(R)))
    for n in 0:(steps - 1)
        t = n * dt
        for k in eachindex(S)
            a1 = u(t, S[k], R[k])
            a2 = u(t + dt / 2, S[k] + dt / 2 * a1[1], R[k] + dt / 2 * a1[2])
            a3 = u(t + dt / 2, S[k] + dt / 2 * a2[1], R[k] + dt / 2 * a2[2])
            a4 = u(t + dt, S[k] + dt * a3[1], R[k] + dt * a3[2])
            S[k] += dt / 6 * (a1[1] + 2a2[1] + 2a3[1] + a4[1])
            R[k] += dt / 6 * (a1[2] + 2a2[2] + 2a3[2] + a4[2])
        end
        history === nothing || push!(history, (copy(S), copy(R)))
    end
    S, R
end

"""
    projection_targets(hp::HatProjection, x, v, w; L, V, steps = 40, density = nothing)
        -> (; x, v, w)

The target points and target weights of the markers `(x_k, v_k, w_k)` on `[0, L) × [−V, V]`: the
lifted target positions and the target velocities of [`moser_targets`](@ref), and the target
weights `w*_k = c f(x*_k, v*_k)` with `c = Σ_j w_j / Σ_j f(x*_j, v*_j)`, so that the total weight
is exact. With `density = nothing`, `f` is the lumped projection `f_h` of the represented density;
otherwise `density(xs, vs)` returns the values of `f` at the vectors of target points.
"""
function projection_targets(hp::HatProjection, x, v, w; L, V, steps = 40, density = nothing)
    s_t, r_t = moser_targets(hp, x ./ L, v ./ V; steps)
    xt, vt = L .* s_t, V .* r_t
    (; x = xt, v = vt, w = projection_target_weights(hp, x, v, w, xt, vt; L, V, density))
end

"""
    projection_target_weights(hp::HatProjection, x, v, w, xt, vt; L, V, density = nothing)

The target weights `w*_k = c f(x*_k, v*_k)`, `c = Σ_j w_j / Σ_j f(x*_j, v*_j)`, at the target
points `(xt, vt)` of the markers `(x, v, w)`: `f` is the lumped projection `f_h` of the markers
(`density = nothing`, on the hats `hp`) or `density(xt, vt)`.
"""
function projection_target_weights(
        hp::HatProjection, x, v, w, xt, vt; L, V, density = nothing)
    if density === nothing
        fh = lumped_projection(hp, x ./ L, v ./ V, w)
        fval = [hat_interpolate(hp, fh, xt[k] / L, vt[k] / V) / (L * V)
                for k in eachindex(xt)]
    else
        fval = density(xt, vt)
    end
    (sum(w) / sum(fval)) .* fval
end
