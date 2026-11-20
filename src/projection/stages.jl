# The stages of the remap map: kicks with the centroid shift, drifts, the naive kick of the remark
# on the centroid shift, and the non-separable implicit-midpoint stage in its velocity form.
#
# An ensemble of N markers in d dimensions holds the lifted positions x (d × N), the canonical
# momenta p (d × N) and the weights w (N).

"""
    centroid_kick(x, p, w, gκ)

The exact time-one flow of the kick Hamiltonian `H^K(x) = Σ_k κ(x_k − X_c(x), w_k)`, with the
weighted centroid `X_c = Σ_j α_j x_j`, `α_j = w_j / Σ_i w_i`:
`p_k ↦ p_k − ∇_ξκ(x_k − X_c, w_k) + α_k Σ_j ∇_ξκ(x_j − X_c, w_j)`. The function `gκ(ξ, w)` returns
the gradient `∇_ξ κ` as a vector; `x` and `p` are `d × N` matrices.
"""
function centroid_kick(x, p, w, gκ)
    α = w ./ sum(w)
    Xc = x * α
    G = reduce(hcat, [gκ(x[:, k] .- Xc, w[k]) for k in axes(x, 2)])
    S = sum(G; dims = 2)
    p .- G .+ S .* α'
end

"""
    naive_kick(x, p, w, gκ)

The kick `p_k ↦ p_k − ∇κ(x_k, w_k) + α_k Σ_j ∇κ(x_j, w_j)` with the mean force removed at
unshifted positions. It conserves the total momentum but is not symplectic: the case that must
fail against [`centroid_kick`](@ref).
"""
function naive_kick(x, p, w, gκ)
    α = w ./ sum(w)
    G = reduce(hcat, [gκ(x[:, k], w[k]) for k in axes(x, 2)])
    S = sum(G; dims = 2)
    p .- G .+ S .* α'
end

"""
    velocity_drift(x, p, w, m, gη)

The exact time-one flow of the drift Hamiltonian `H^D(p) = Σ_k η(p_k/(m w_k), w_k)`:
`x_k ↦ x_k + ∇_vη(p_k/(m w_k), w_k) / (m w_k)`. The function `gη(v, w)` returns `∇_v η`.
"""
function velocity_drift(x, p, w, m, gη)
    Δ = reduce(hcat, [gη(p[:, k] ./ (m * w[k]), w[k]) ./ (m * w[k]) for k in axes(x, 2)])
    x .+ Δ
end

# The size of a number for a convergence test: for a dual number, the value and all partials count,
# so that an iteration differentiated by ForwardDiff converges in its derivatives too.
convergence_size(a::Real) = abs(a)
function convergence_size(a::ForwardDiff.Dual)
    abs(ForwardDiff.value(a)) + sum(convergence_size, ForwardDiff.partials(a))
end

"""
    nonseparable_stage(x, p, w, m, gφ; tol = 1e-13, maxit = 200) -> (x′, p′, iterations, converged)

One step of the implicit midpoint rule for `H(x, p) = Σ_k m w_k φ(x_k − X_c(x), p_k/(m w_k), w_k)`,
solved in the velocity form: each marker takes an implicit-midpoint step of `φ` in `(ξ, v)`, and
all markers share the centroid shift `c = Σ_j α_j ∇_vφ_j` and the boost `b = Σ_j α_j ∇_ξφ_j` at the
midpoints. The function `gφ(ξ, v, w)` returns the pair `(∇_ξφ, ∇_vφ)`. The fixed-point iteration
stops when the update, partials of dual numbers included, is below `tol` relative to the size of
the iterate; without convergence in `maxit` iterations the input is returned with
`converged = false`.
"""
function nonseparable_stage(x, p, w, m, gφ; tol = 1e-13, maxit = 200)
    α = w ./ sum(w)
    ξ = x .- x * α
    v = p ./ (m .* w')
    function slopes(ξ1, v1)
        G = [gφ((ξ[:, k] .+ ξ1[:, k]) ./ 2, (v[:, k] .+ v1[:, k]) ./ 2, w[k])
             for k in axes(x, 2)]
        reduce(hcat, first.(G)), reduce(hcat, last.(G))
    end
    ξ1, v1 = ξ, v
    for it in 1:maxit
        Gξ, Gv = slopes(ξ1, v1)
        ξn = ξ .+ Gv .- Gv * α
        vn = v .- Gξ .+ Gξ * α
        δ = max(maximum(convergence_size, ξn .- ξ1), maximum(convergence_size, vn .- v1))
        scale = max(maximum(convergence_size, ξn), maximum(convergence_size, vn))
        ξ1, v1 = ξn, vn
        if δ < tol * (1 + scale)
            Gξ, Gv = slopes(ξ1, v1)
            return x .+ Gv, p .- m .* w' .* (Gξ .- Gξ * α), it, true
        end
    end
    x, p, maxit, false
end

"""
    nonseparable_hamiltonian(φ, d, N, w; m = 1, centroid = true)

The ensemble Hamiltonian `H(z) = Σ_k m w_k φ(x_k − X_c(x), p_k/(m w_k), w_k)` of `N` markers in
`d` dimensions as a function of the packed state `z` ([`pack`](@ref ParticleRemapping.pack)). With `centroid = false`
the positions are not shifted by the centroid: the case that must fail to conserve the momentum.
"""
function nonseparable_hamiltonian(φ, d, N, w; m = 1, centroid = true)
    z -> begin
        x, p = unpack(z, d, N)
        Xc = centroid ? x * (w ./ sum(w)) : zeros(d)
        sum(m * w[k] * φ(x[:, k] .- Xc, p[:, k] ./ (m * w[k]), w[k]) for k in 1:N)
    end
end

"""
    hamiltonian_one_step(H, z; rule = :midpoint, tol = 1e-14, maxit = 400)

One step of the implicit midpoint rule `z′ = z + J ∇H((z + z′)/2)` (`rule = :midpoint`) or of the
trapezoidal rule `z′ = z + ½ J (∇H(z) + ∇H(z′))` (`rule = :trapezoidal`) for the Hamiltonian `H`
of the packed state `z`, by fixed-point iteration with gradients by forward-mode differentiation.
"""
function hamiltonian_one_step(H, z; rule = :midpoint, tol = 1e-14, maxit = 400)
    rule in (:midpoint, :trapezoidal) ||
        throw(ArgumentError("rule must be :midpoint or :trapezoidal, got $rule"))
    J = canonical_J(length(z) ÷ 2)
    gH(y) = ForwardDiff.gradient(H, y)
    g0 = rule === :trapezoidal ? gH(z) : nothing
    z1 = z
    for _ in 1:maxit
        zn = rule === :midpoint ? z + J * gH((z + z1) / 2) : z + J * (g0 + gH(z1)) / 2
        δ = maximum(convergence_size, zn - z1)
        z1 = zn
        δ < tol && return z1
    end
    error("the $rule rule did not converge in $maxit iterations")
end

"""
    marker_midpoint_step(φ, W, y; tol = 1e-14, maxit = 400)

The implicit-midpoint step `y′ = y + J ∇φ((y + y′)/2, W(y))` of one marker in `y = (x, v)` with the
centroid and the boost held fixed; the weight is the label `W(y)` at the start point.
"""
function marker_midpoint_step(φ, W, y; tol = 1e-14, maxit = 400)
    n = length(y) ÷ 2
    J = canonical_J(n)
    w = W(y)
    gφ(u) = ForwardDiff.gradient(s -> φ(s[1:n], s[(n + 1):end], w), u)
    y1 = y
    for _ in 1:maxit
        yn = y + J * gφ((y + y1) / 2)
        δ = maximum(convergence_size, yn - y1)
        y1 = yn
        δ < tol && return y1
    end
    error("marker_midpoint_step did not converge in $maxit iterations")
end

"""
    nonseparable_jacobian_determinant(φ, W, y, y′) -> (; exact, first, wrong)

The Jacobian determinant `1 + ∇W(y)ᵀ (I + ½ J K)⁻¹ J ∇∂_wφ` of the marker step
[`marker_midpoint_step`](@ref) from `y` to `y′`, with `K = ∇²φ` and `∇∂_wφ` at the midpoint and the
weight `W(y)` (`exact`); its first-order value `1 + ∇Wᵀ J ∇∂_wφ` (`first`); and the variant with
`(I − ½ J K)⁻¹` (`wrong`). The last two are the cases that must fail.
"""
function nonseparable_jacobian_determinant(φ, W, y, y1)
    n = length(y) ÷ 2
    J = canonical_J(n)
    w = W(y)
    yb = (y + y1) / 2
    Φ(u, ω) = φ(u[1:n], u[(n + 1):end], ω)
    K = ForwardDiff.hessian(u -> Φ(u, w), yb)
    g = ForwardDiff.gradient(u -> ForwardDiff.derivative(ω -> Φ(u, ω), w), yb)
    ∇W = ForwardDiff.gradient(W, y)
    (exact = 1 + dot(∇W, (I + J * K / 2) \ (J * g)),
        first = 1 + dot(∇W, J * g),
        wrong = 1 + dot(∇W, (I - J * K / 2) \ (J * g)))
end

"""
    nonseparable_divergence(φ, W, y) -> (div, bracket)

The divergence of `u(y) = J ∇φ(y, W(y))` in `y = (x, v)` and the bracket `{W, ∂_wφ}_v` at `y`.
"""
function nonseparable_divergence(φ, W, y)
    n = length(y) ÷ 2
    J = canonical_J(n)
    Φ(u, ω) = φ(u[1:n], u[(n + 1):end], ω)
    u(s) = J * ForwardDiff.gradient(r -> Φ(r, W(s)), s)
    div = tr(ForwardDiff.jacobian(u, y))
    g = ForwardDiff.gradient(r -> ForwardDiff.derivative(ω -> Φ(r, ω), W(y)), y)
    div, dot(ForwardDiff.gradient(W, y), J * g)
end

"""
    kick_drift_jacobian_determinants(Δv, Δx, W, z) -> (kick, drift)

The Jacobian determinants `1 + ∂_ℓΔv ∂_v ln W` and `1 + ∂_ℓΔx ∂_x ln W` of a kick
`v ↦ v + Δv(x, ℓ)` and a drift `x ↦ x + Δx(v, ℓ)` in `z = (x, v)`, `d = 1`, for markers with
`ℓ = ln W(x, v)`, with the centroid and the boost held fixed.
"""
function kick_drift_jacobian_determinants(Δv, Δx, W, z)
    ℓ = log(W(z...))
    gW = ForwardDiff.gradient(q -> log(W(q...)), z)
    (1 + ForwardDiff.derivative(l -> Δv(z[1], l), ℓ) * gW[2],
        1 + ForwardDiff.derivative(l -> Δx(z[2], l), ℓ) * gW[1])
end

"""
    pm_v_coupling(; vs = range(-3, 3; length = 121), x = π)

`det DN − 1` of one marker step [`marker_midpoint_step`](@ref) at the position `x`, for the
velocities `vs`, the weight field `W = (1 + 0.3 cos x) exp(−v²/2)`, the separable
`φ = 0.1 sin ξ ln w + 0.05 v² ln w` and the non-separable `φ = 0.1 v sin ξ ln w`. The data of
[`plot_pm_v_coupling`](@ref).
"""
function pm_v_coupling(; vs = range(-3, 3; length = 121), x = π)
    W(y) = (1 + 0.3 * cos(y[1])) * exp(-y[2]^2 / 2)
    separable(ξ, v, w) = 0.1 * sin(ξ[1]) * log(w) + 0.05 * v[1]^2 * log(w)
    nonseparable(ξ, v, w) = 0.1 * v[1] * sin(ξ[1]) * log(w)
    detm1(φ, v) = det(ForwardDiff.jacobian(s -> marker_midpoint_step(φ, W, s), [
        Float64(x), v])) - 1
    (; v = collect(vs), separable = [detm1(separable, v) for v in vs],
        nonseparable = [detm1(nonseparable, v) for v in vs], x)
end
