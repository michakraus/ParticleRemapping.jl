"""
    pack(x, p)

The packed state `z = [vec(x); vec(p)]` of an ensemble with positions `x` and momenta `p`. Its
canonical symplectic matrix is [`canonical_J`](@ref)`(length(x))`.
"""
pack(x, p) = vcat(vec(x), vec(p))

"""
    unpack(z, d, N) -> (x, p)

The inverse of [`pack`](@ref) for `N` markers in `d` dimensions; `x` and `p` are `d × N` matrices.
"""
function unpack(z, d::Integer, N::Integer)
    n = d * N
    reshape(z[1:n], d, N), reshape(z[(n + 1):(2n)], d, N)
end

"""
    canonical_J(n)

The canonical symplectic matrix `[0 I; -I 0]` of size `2n`, the matrix of `Σ dq ∧ dp` in the
packed coordinates `z = (q, p)`.
"""
function canonical_J(n::Integer)
    [zeros(n, n) Matrix{Float64}(I, n, n); -Matrix{Float64}(I, n, n) zeros(n, n)]
end

"""
    symplecticity_defect(R, z; jacobian = ForwardDiff.jacobian) -> (defect, ‖DR‖)

The largest entry of `DRᵀ J DR − J`, with `DR` the Jacobian of the map `R` at the packed state
`z`, and the operator norm of `DR`. The defect is zero (to round-off) for a symplectic map. Pass
`jacobian = fd_jacobian` for a map that does not accept dual numbers.
"""
function symplecticity_defect(R, z; jacobian = ForwardDiff.jacobian)
    DR = jacobian(R, z)
    J = canonical_J(length(z) ÷ 2)
    maximum(abs, DR' * J * DR - J), opnorm(DR)
end

"""
    fd_jacobian(F, z; h = 1e-5, order = 2)

The Jacobian of `F` at `z` by central differences of second (`order = 2`) or fourth
(`order = 4`) order with step `h`.
"""
function fd_jacobian(F, z::AbstractVector; h::Real = 1e-5, order::Integer = 2)
    order in (2, 4) || throw(ArgumentError("order must be 2 or 4, got $order"))
    F0 = F(z)
    Jac = zeros(promote_type(eltype(F0), Float64), length(F0), length(z))
    e = zeros(length(z))
    for j in eachindex(z)
        e[j] = h
        Jac[:, j] = order == 2 ? (F(z .+ e) .- F(z .- e)) ./ (2h) :
                    (8 .* (F(z .+ e) .- F(z .- e)) .- (F(z .+ 2e) .- F(z .- 2e))) ./ (12h)
        e[j] = 0
    end
    Jac
end

"""
    presymplectic_matrix(z; m = 1)

The matrix of the presymplectic form `ω_N = Σ_a m (w_a dx_a ∧ dv_a + v_a dx_a ∧ dw_a)` at the
state `z = (x_1, v_1, w_1, x_2, v_2, w_2, …)` of 1D-1V markers. It is the pull-back of
`Σ dx_a ∧ dp_a` by `(x, v, w) ↦ (x, p = m w v)` and is degenerate: its kernel is spanned by
`w_a ∂_{w_a} − v_a ∂_{v_a}`.
"""
function presymplectic_matrix(z::AbstractVector; m::Real = 1)
    length(z) % 3 == 0 || throw(DimensionMismatch("z must hold (x, v, w) triples"))
    Ω = zeros(eltype(z), length(z), length(z))
    for a in 0:(length(z) ÷ 3 - 1)
        ix, iv, iw = 3a + 1, 3a + 2, 3a + 3
        v, w = z[iv], z[iw]
        Ω[ix, iv], Ω[iv, ix] = m * w, -m * w
        Ω[ix, iw], Ω[iw, ix] = m * v, -m * v
    end
    Ω
end

"""
    pullback_residual(Φ, z; m = 1, jacobian = fd_jacobian)

The largest entry of `DΦᵀ Ω(Φ(z)) DΦ − Ω(z)` with `Ω` the [`presymplectic_matrix`](@ref): zero
when the map `Φ` of 1D-1V marker states `(x, v, w)` preserves the presymplectic form `ω_N`.
"""
function pullback_residual(Φ, z::AbstractVector; m::Real = 1, jacobian = fd_jacobian)
    DΦ = jacobian(Φ, z)
    maximum(abs, DΦ' * presymplectic_matrix(Φ(z); m = m) * DΦ -
                 presymplectic_matrix(z; m = m))
end

export pack, unpack, canonical_J, symplecticity_defect, fd_jacobian
export presymplectic_matrix, pullback_residual
