# The structural facts of the Hamiltonian particle projection: mean-field Hamiltonians move all
# markers with one map, the compressibility of a weight-dependent Hamiltonian, and the fold of the
# implicit weight.

# ---- mean-field Hamiltonians ----------------------------------------------------------------------

"""
    mean_field_example(; τ = 0.2, center = [0.1, 0.2], radius = 0.35, n = 17, half_width = 1.2,
        jitter = 0.01, ncircle = 4000, rng = Random.default_rng())

The mean-field example of the documentation section on holes: a jittered `n × n` lattice of
markers `ζ = (ξ, η)` in `[−half_width, half_width]²` with the disk of `radius` around `center` cut
out, the moments `φ(ζ) = (sin ξ η, cos ξ + 0.3 η²)`, the single-particle Hamiltonian
`h(ζ; Λ) = η²/2 + 0.4 cos ξ + Λ₁ sin(ξ + η/2) + Λ₁ Λ₂ cos ξ`, which is nonlinear in the moments
`Λ = (1/N) Σ_k φ(z_k)`, and its weight-dependent control `h_w(ζ, w; Λ) = h(ζ; Λ) + ½ ln w sin 2ξ η`.

Returns a `NamedTuple` with the step `τ`, the functions `φ`, `h`, `h_w`, the markers `Z0` (a
vector of 2-vectors), the boundary `circle` of the hole (`ncircle` points), its `center` and
`radius`.
"""
function mean_field_example(; τ = 0.2, center = [0.1, 0.2], radius = 0.35, n = 17,
        half_width = 1.2, jitter = 0.01, ncircle = 4000, rng::AbstractRNG = Random.default_rng())
    φ(ζ) = [sin(ζ[1]) * ζ[2], cos(ζ[1]) + 0.3 * ζ[2]^2]
    h(ζ, Λ) = 0.5 * ζ[2]^2 + 0.4 * cos(ζ[1]) + Λ[1] * sin(ζ[1] + 0.5ζ[2]) +
              Λ[1] * Λ[2] * cos(ζ[1])
    h_w(ζ, w, Λ) = h(ζ, Λ) + 0.5 * log(w) * sin(2ζ[1]) * ζ[2]
    lattice = [[ξ, η] for ξ in range(-half_width, half_width; length = n)
               for η in range(-half_width, half_width; length = n)]
    Z0 = filter(ζ -> norm(ζ - center) > radius, [ζ + jitter * randn(rng, 2)
                                                 for ζ in lattice])
    circle = [center + radius * [cos(t), sin(t)]
              for t in range(0, 2π; length = ncircle + 1)[1:(end - 1)]]
    (; τ, φ, h, h_w, Z0, circle, center, radius)
end

"""
    mean_field_hamiltonian(h, φ)

The ensemble Hamiltonian `H(Z) = Σ_k h(z_k; Λ)` with the moments `Λ = (1/N) Σ_k φ(z_k)` of the
markers `Z` (a vector of points).
"""
function mean_field_hamiltonian(h, φ)
    Z -> begin
        Λ = sum(φ, Z) / length(Z)
        sum(ζ -> h(ζ, Λ), Z)
    end
end

"""
    mean_field_single_hamiltonian(h, φ, Z)

The single-particle Hamiltonian `Ψ_μ(ζ) = h(ζ; Λ) + (1/N) Σ_j ∂_Λ h(z_j; Λ) · φ(ζ)` of
[`mean_field_hamiltonian`](@ref)`(h, φ)` for the empirical measure of the points `Z`. Its gradient
at `z_k` is the gradient of the ensemble Hamiltonian with respect to `z_k`.
"""
function mean_field_single_hamiltonian(h, φ, Z)
    Λ = sum(φ, Z) / length(Z)
    G = sum(ζ -> ForwardDiff.gradient(L -> h(ζ, L), Λ), Z)
    N = length(Z)
    ζ -> h(ζ, Λ) + dot(G, φ(ζ)) / N
end

_to_points(v) = [v[(2k - 1):(2k)] for k in 1:(length(v) ÷ 2)]
_jvec(g) = [g[2], -g[1]]

"""
    ensemble_gradient(H, Z)

The gradients `∇_{z_k} H` of an ensemble Hamiltonian `H` of the points `Z` (2-vectors), one per
point, by forward-mode differentiation.
"""
function ensemble_gradient(H, Z)
    _to_points(ForwardDiff.gradient(v -> H(_to_points(v)), reduce(vcat, Z)))
end

"""
    ensemble_symplectic_euler(H, Z, τ; iters = 200) -> (Z′, Ẑ)

One step of symplectic Euler, `η′ = η − τ ∂_ξ H(ξ, η′)`, `ξ′ = ξ + τ ∂_η H(ξ, η′)`, for the
ensemble Hamiltonian `H` of the points `Z = (ξ_k, η_k)`, solved by `iters` fixed-point iterations.
Returns the new points and the points `(ξ_k, η′_k)` at the mixed arguments.
"""
function ensemble_symplectic_euler(H, Z, τ; iters = 200)
    P = [ζ[2] for ζ in Z]
    for _ in 1:iters
        G = ensemble_gradient(H, [[Z[k][1], P[k]] for k in eachindex(Z)])
        P = [Z[k][2] - τ * G[k][1] for k in eachindex(Z)]
    end
    G = ensemble_gradient(H, [[Z[k][1], P[k]] for k in eachindex(Z)])
    [[Z[k][1] + τ * G[k][2], P[k]] for k in eachindex(Z)],
    [[Z[k][1], P[k]] for k in eachindex(Z)]
end

"""
    ensemble_implicit_midpoint(H, Z, τ; iters = 200) -> (Z′, Z̄)

One step of the implicit midpoint rule `z′ = z + τ J ∇H((z + z′)/2)` for the ensemble Hamiltonian
`H` of the points `Z`, solved by `iters` fixed-point iterations. Returns the new points and the
midpoints.
"""
function ensemble_implicit_midpoint(H, Z, τ; iters = 200)
    Zn = copy(Z)
    for _ in 1:iters
        M = [(Z[k] + Zn[k]) / 2 for k in eachindex(Z)]
        G = ensemble_gradient(H, M)
        Zn = [Z[k] + τ * _jvec(G[k]) for k in eachindex(Z)]
    end
    Zn, [(Z[k] + Zn[k]) / 2 for k in eachindex(Z)]
end

"""
    particle_symplectic_euler(Ψ, ζ, τ; iters = 200)

One step of symplectic Euler of the single-particle Hamiltonian `Ψ` from the point `ζ = (ξ, η)`.
"""
function particle_symplectic_euler(Ψ, ζ, τ; iters = 200)
    η = ζ[2]
    for _ in 1:iters
        η = ζ[2] - τ * ForwardDiff.gradient(Ψ, [ζ[1], η])[1]
    end
    [ζ[1] + τ * ForwardDiff.gradient(Ψ, [ζ[1], η])[2], η]
end

"""
    particle_implicit_midpoint(Ψ, ζ, τ; iters = 200)

One step of the implicit midpoint rule of the single-particle Hamiltonian `Ψ` from the point `ζ`.
"""
function particle_implicit_midpoint(Ψ, ζ, τ; iters = 200)
    ζn = ζ
    for _ in 1:iters
        ζn = ζ + τ * _jvec(ForwardDiff.gradient(Ψ, (ζ + ζn) / 2))
    end
    ζn
end

"""
    shoelace_area(P)

The area of the closed polygon with the vertices `P` (2-vectors), by the shoelace formula.
"""
function shoelace_area(P)
    0.5 * abs(sum(P[k][1] * P[mod1(k + 1, end)][2] - P[mod1(k + 1, end)][1] * P[k][2]
    for k in eachindex(P)))
end

"""
    inside_polygon(P, q)

Whether the point `q` lies inside the closed polygon with the vertices `P`, by the even–odd rule.
"""
function inside_polygon(P, q)
    c = false
    for k in eachindex(P)
        a, b = P[k], P[mod1(k + 1, length(P))]
        if (a[2] > q[2]) != (b[2] > q[2])
            xcross = a[1] + (q[2] - a[2]) * (b[1] - a[1]) / (b[2] - a[2])
            q[1] < xcross && (c = !c)
        end
    end
    c
end

"""
    mean_field_step_example(ex = mean_field_example(); every = 10)

One implicit-midpoint step of the mean-field example `ex`: the markers before and after the step,
the single-particle map `ψ` that moves all of them, and the hole boundary (every `every`-th point
of `ex.circle`) and its image under `ψ`. The data of [`plot_mean_field_step`](@ref).
"""
function mean_field_step_example(ex = mean_field_example(); every = 10)
    H = mean_field_hamiltonian(ex.h, ex.φ)
    Zn, M = ensemble_implicit_midpoint(H, ex.Z0, ex.τ)
    Ψ = mean_field_single_hamiltonian(ex.h, ex.φ, M)
    ψ(ζ) = particle_implicit_midpoint(Ψ, ζ, ex.τ)
    circle = ex.circle[1:every:end]
    (; before = ex.Z0, after = Zn, ψ, circle, image = ψ.(circle))
end

# ---- compressibility of a weight-dependent Hamiltonian ------------------------------------------

"""
    canonical_bracket(A, B, ζ, d)

The canonical Poisson bracket `{A, B} = Σ_i ∂_{x^i}A ∂_{p_i}B − ∂_{p_i}A ∂_{x^i}B` of two functions
of `ζ = (x, p) ∈ ℝ^{2d}` at `ζ`. With `ζ = (x, v)` it is the bracket `{A, B}_v` in `(x, v)`.
"""
function canonical_bracket(A, B, ζ, d)
    gA = ForwardDiff.gradient(A, ζ)
    gB = ForwardDiff.gradient(B, ζ)
    dot(gA[1:d], gB[(d + 1):(2d)]) - dot(gA[(d + 1):(2d)], gB[1:d])
end

"""
    effective_divergence(Ψ, W, ζ, d)

The divergence at `ζ = (x, p)` of the effective velocity field `u(ζ) = J ∇_ζ Ψ(x, p, W(x, p))` of
markers whose weights are the values of the field `W`; the gradient is taken at fixed weight.
"""
function effective_divergence(Ψ, W, ζ, d)
    function u(z)
        x, p = z[1:d], z[(d + 1):(2d)]
        w = W(x, p)
        g = ForwardDiff.gradient(η -> Ψ(η[1:d], η[(d + 1):(2d)], w), z)
        vcat(g[(d + 1):(2d)], -g[1:d])
    end
    tr(ForwardDiff.jacobian(u, ζ))
end

"""
    frozen_weight_divergence(Ψ, W, ζ, d)

The divergence at `ζ` of the field `J ∇_ζ Ψ(x, p, w₀)` with the weight frozen at `w₀ = W(ζ)`. It
vanishes: this is the case that must disagree with [`effective_divergence`](@ref).
"""
function frozen_weight_divergence(Ψ, W, ζ, d)
    w0 = W(ζ[1:d], ζ[(d + 1):(2d)])
    function u(z)
        g = ForwardDiff.gradient(η -> Ψ(η[1:d], η[(d + 1):(2d)], w0), z)
        vcat(g[(d + 1):(2d)], -g[1:d])
    end
    tr(ForwardDiff.jacobian(u, ζ))
end

"""
    compressibility_bracket(Ψ, W, ζ, d; total = false)

The bracket `{W, ∂_wΨ}` at `ζ = (x, p)`, evaluated at `w = W(ζ)`. With `total = false` the
derivative `∂_wΨ` is taken as a function of `(x, p)` at the fixed weight `W(ζ)`; with
`total = true` at the weight `W(x, p)` of each point, which adds a term `∂²_wΨ ∇W` that the
bracket with `W` cancels.
"""
function compressibility_bracket(Ψ, W, ζ, d; total = false)
    split(z) = (z[1:d], z[(d + 1):(2d)])
    Wζ(z) = W(split(z)...)
    w0 = Wζ(ζ)
    ∂wΨ(z) = ForwardDiff.derivative(w -> Ψ(split(z)..., w), total ? Wζ(z) : w0)
    canonical_bracket(Wζ, ∂wΨ, ζ, d)
end

# ---- the fold of the implicit weight ------------------------------------------------------------

"""
    implicit_weight_map(f, v; m = 1, c = 1)

The map `Φ_{x′}(v) = m c v f(v)` at fixed position, whose inverse solves the implicit weight
equation `w′ = c f(x′, p′/(m w′))` for `v = p′/(m w′)`.
"""
implicit_weight_map(f, v; m = 1, c = 1) = m * c * v * f(v)

"""
    fold_determinant(f, v; m = 1, c = 1)

The determinant `(m c)^d f^{d−1} (f + v · ∇_v f)` of the Jacobian of
[`implicit_weight_map`](@ref) at `v ∈ ℝ^d`.
"""
function fold_determinant(f, v::AbstractVector; m = 1, c = 1)
    d = length(v)
    fv = f(v)
    (m * c)^d * fv^(d - 1) * (fv + dot(v, ForwardDiff.gradient(f, v)))
end

_maxwellian(r, vth) = exp(-r^2 / (2vth^2)) / sqrt(2π * vth^2)

"""
    maxwellian_fold_peak(; m = 1, c = 1, vth = 1)

The largest value `p_max = m c v_th f_M(v_th)` of `|Φ_{x′}|` on a ray for the one-dimensional
Maxwellian `f_M` of thermal speed `vth`.
"""
maxwellian_fold_peak(; m = 1, c = 1, vth = 1) = m * c * vth * _maxwellian(vth, vth)

"""
    maxwellian_fold_roots(p; m = 1, c = 1, vth = 1, n = 200_001)

The number of solutions `r > 0` of `m c r f_M(r) = p` on a ray, for the Maxwellian `f_M` of
thermal speed `vth`, counted by the sign changes on `n` points of `(0, 12 vth]`.
"""
function maxwellian_fold_roots(p; m = 1, c = 1, vth = 1, n = 200_001)
    r = range(1e-6, 12vth; length = n)
    g = m * c .* r .* _maxwellian.(r, vth) .- p
    count(i -> sign(g[i]) != sign(g[i + 1]), 1:(length(g) - 1))
end
