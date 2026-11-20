# The canonical structure of weighted markers: the pull-back of ω to (x, v, w) and the
# velocity–momentum fold of a density-tied weight.

"""
    kinematic_projection(y; m = 1)

The kinematic projection `π(x, v, w) = (x, p = m w v)` of the 1D-1V marker state
`y = (x₁, v₁, w₁, x₂, v₂, w₂, …)` to the packed canonical state `z = (q₁, …, q_N, p₁, …, p_N)`
of [`pack`](@ref).
"""
function kinematic_projection(y::AbstractVector; m::Real = 1)
    length(y) % 3 == 0 || throw(DimensionMismatch("y must hold (x, v, w) triples"))
    x, v, w = y[1:3:end], y[2:3:end], y[3:3:end]
    pack(x, m .* w .* v)
end

"""
    kinematic_pullback(y; m = 1, jacobian = ForwardDiff.jacobian)

The matrix `Dπᵀ J Dπ` of the pull-back `π*ω` of the canonical form `ω = Σ dq_a ∧ dp_a` by the
[`kinematic_projection`](@ref) at the state `y = (x₁, v₁, w₁, …)`. It equals the presymplectic
matrix [`presymplectic_matrix`](@ref ParticleRemapping.presymplectic_matrix) of `Σ m w_a dx_a ∧ dv_a + Σ m v_a dx_a ∧ dw_a`. Pass
`jacobian = (F, y) -> fd_jacobian(F, y; order = 4)` for a finite-difference Jacobian.
"""
function kinematic_pullback(y::AbstractVector; m::Real = 1, jacobian = ForwardDiff.jacobian)
    Dπ = jacobian(z -> kinematic_projection(z; m = m), y)
    transpose(Dπ) * canonical_J(size(Dπ, 1) ÷ 2) * Dπ
end

"""
    euler_scaling_fields(y)

The `3N × N` matrix whose column `a` is the Euler scaling field `K_a = w_a ∂_{w_a} − v_a ∂_{v_a}`
at the 1D-1V state `y = (x₁, v₁, w₁, …)`. The fields span the kernel of the pull-back
[`kinematic_pullback`](@ref): they keep every `w_a v_a` at fixed `x_a`.
"""
function euler_scaling_fields(y::AbstractVector)
    N = length(y) ÷ 3
    K = zeros(eltype(y), length(y), N)
    for a in 1:N
        K[3a - 1, a] = -y[3a - 1]
        K[3a, a] = y[3a]
    end
    K
end

"""
    velocity_fold_momentum(fh, v; m = 1, ΔΩ = 1)

The canonical momentum `P(v) = m ΔΩ v f_h(v)` of a marker at velocity `v` (a number or a vector
of length `d`) whose weight is tied to the density `fh` at fixed position by `w = ΔΩ f_h(v)`.
"""
velocity_fold_momentum(fh, v; m::Real = 1, ΔΩ::Real = 1) = m * ΔΩ .* v .* fh(v)

"""
    velocity_fold_determinant(fh, v; m = 1, ΔΩ = 1)

The determinant `(m ΔΩ)^d f_h^{d−1} (f_h + v ⋅ ∇_v f_h)` of the Jacobian of
[`velocity_fold_momentum`](@ref) at `v`, with `∇_v f_h` by automatic differentiation. It vanishes on the
fold `f_h + v ⋅ ∇_v f_h = 0`.
"""
function velocity_fold_determinant(fh, v::Real; m::Real = 1, ΔΩ::Real = 1)
    m * ΔΩ * (fh(v) + v * ForwardDiff.derivative(fh, v))
end

function velocity_fold_determinant(fh, v::AbstractVector; m::Real = 1, ΔΩ::Real = 1)
    f = fh(v)
    (m * ΔΩ)^length(v) * f^(length(v) - 1) * (f + dot(v, ForwardDiff.gradient(fh, v)))
end

"""
    isotropic_maxwellian(v; vth = 1, C = 1)

The isotropic Maxwellian `C exp(−|v|² / (2 v_th²))` at the velocity `v` (a number or a vector).
"""
isotropic_maxwellian(v; vth::Real = 1, C::Real = 1) = C * exp(-sum(abs2, v) / (2 * vth^2))

"""
    maxwellian_fold_momentum_bound(; vth = 1, C = 1, m = 1, ΔΩ = 1)

The largest canonical momentum `p_max = m ΔΩ v_th C e^{−1/2}` that the
[`velocity_fold_momentum`](@ref) of the [`isotropic_maxwellian`](@ref) reaches, on the thermal sphere
`|v| = v_th`.
"""
function maxwellian_fold_momentum_bound(; vth::Real = 1, C::Real = 1, m::Real = 1,
        ΔΩ::Real = 1)
    m * ΔΩ * vth * C * exp(-one(float(vth)) / 2)
end
