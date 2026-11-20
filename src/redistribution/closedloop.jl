# Open-loop versus closed-loop symplecticity: the identity
# Φ*ω − ω = −τ Σ_k d(∂_{θ_k} H)(z_stage) ∧ dΘ_k, its two set-ups (target matching and the
# variational cotangent-lift closure), and the exact Jacobians by implicit differentiation.

"""
    closed_loop_defect(DΦ, DΘ, Dg; τ = 1) -> (lhs, rhs)

Both sides of the closed-loop identity in matrix form,
`lhs = DΦᵀ J DΦ − J` and `rhs = −τ (Dgᵀ DΘ − DΘᵀ Dg)`, for the Jacobian `DΦ` of the closed-loop
map, the Jacobian `DΘ` of the parameter map and the Jacobian `Dg` of
`g_k(z) = ∂H/∂θ_k(z_stage(z); Θ(z))`.
"""
function closed_loop_defect(DΦ::AbstractMatrix, DΘ::AbstractMatrix, Dg::AbstractMatrix;
        τ::Real = 1)
    J = canonical_J(size(DΦ, 1) ÷ 2)
    transpose(DΦ) * J * DΦ - J, -τ .* (transpose(Dg) * DΘ - transpose(DΘ) * Dg)
end

"""
    implicit_step_jacobian(R, z̄, z)

The Jacobian `DΦ(z) = −(∂R/∂z̄)⁻¹ ∂R/∂z` of the map `z ↦ z̄` defined implicitly by
`R(z̄, z) = 0`, at a root `z̄`, with both partial Jacobians by automatic differentiation.
"""
function implicit_step_jacobian(R, z̄::AbstractVector, z::AbstractVector)
    Rz̄ = ForwardDiff.jacobian(y -> R(y, z), z̄)
    Rz = ForwardDiff.jacobian(y -> R(z̄, y), z)
    -(Rz̄ \ Rz)
end

"""
    closed_loop_example(; z = [0.5, 0.9, 2.8, 4.7, 0.04, -0.06, 0.03, -0.01], Lx = 2π,
        p_scale = 0.25, nx_cells = 5, np_cells = 4)

The four-marker set-up of the closed-loop checks: the packed state `z`, the
[`PairwiseSplineHamiltonian`](@ref) on `[0, Lx)` with `nx_cells × np_cells` cells, and `Lx`.
A periodic quartic basis needs more than four cells, so `nx_cells = 5` is the smallest
admissible value.
"""
function closed_loop_example(; z = [0.5, 0.9, 2.8, 4.7, 0.04, -0.06, 0.03, -0.01],
        Lx::Real = 2π, p_scale::Real = 0.25, nx_cells::Integer = 5, np_cells::Integer = 4)
    H = PairwiseSplineHamiltonian(Lx, p_scale; nx_cells = nx_cells, np_cells = np_cells)
    (z = collect(float.(z)), H = H, Lx = float(Lx))
end

# The target of the target-matching closed loop: q̄_a = q_a + gain Σ_{b≠a} Δ exp(−Δ²/2),
# Δ = q_a − q_b periodic; p̄ = p.
function _repelled_targets(q::AbstractVector, Lx::Real, gain::Real)
    q̄ = copy(q)
    for a in eachindex(q), b in eachindex(q)

        a == b && continue
        Δ = _periodic_difference(q[a], q[b], Lx)
        q̄[a] += gain * Δ * exp(-Δ^2 / 2)
    end
    q̄
end

"""
    target_matching_parameters(H, z; scheme = ImplicitMidpointStep(), Lx, gain = 0.04,
        λ = 1e-2)

Control E: the parameter map `Θ(z)` of a target-matching closed loop. The targets repel the
positions, `q̄_a = q_a + gain Σ_{b≠a} Δ_{ab} exp(−Δ_{ab}²/2)`, keep `p̄ = p`, and `Θ(z)` is the
Tikhonov-regularised solution (`λ`) of the parameter system of `scheme`
([`solve_redistribution_parameters`](@ref)).
"""
function target_matching_parameters(H::PairwiseSplineHamiltonian, z::AbstractVector;
        scheme::RedistributionScheme = ImplicitMidpointStep(), Lx::Real = H.Lx,
        gain::Real = 0.04, λ::Real = 1e-2)
    n = length(z) ÷ 2
    q, p = z[1:n], z[(n + 1):(2n)]
    q̄ = _repelled_targets(q, Lx, gain)
    solve_redistribution_parameters(H, scheme, q, p, q̄, p; λ = λ)
end

# The packed stage point of `scheme` for the step y ↦ ȳ.
function _packed_stage(scheme, y, ȳ)
    n = length(y) ÷ 2
    vcat(redistribution_stage_point(
        scheme, y[1:n], y[(n + 1):(2n)], ȳ[1:n], ȳ[(n + 1):(2n)])...)
end

"""
    target_matching_closed_loop(H, z; scheme = ImplicitMidpointStep(), τ = 1,
        jacobians = true, kwargs...) -> (z̄, θ, g, jacobians)

The closed-loop map `Φ(z) = Ψ_τ^{Θ(z)}(z)` of Control E: the step of `scheme` with the
[`target_matching_parameters`](@ref) `Θ(z)` (keywords `kwargs`). Returns the end point, the
parameters, `g = Ψ(z_stage)` (the basis Hamiltonians at the actual stage point,
`g_k = ∂H/∂θ_k`) and, for `jacobians = true`, the exact Jacobians `(DΦ, DΘ, Dg)` by implicit
differentiation of the step equations `q̄ − q = τ ∇_p H(z_stage; Θ(z))`,
`p̄ − p = −τ ∇_q H(z_stage; Θ(z))` (otherwise `nothing`).
"""
function target_matching_closed_loop(H::PairwiseSplineHamiltonian, z::AbstractVector;
        scheme::RedistributionScheme = ImplicitMidpointStep(), τ::Real = 1,
        jacobians::Bool = true, kwargs...)
    n = length(z) ÷ 2
    Θ(y) = target_matching_parameters(H, y; scheme = scheme, kwargs...)
    θ = Θ(z)
    z̄ = forward_redistribution_step(H, scheme, z[1:n], z[(n + 1):(2n)], θ; τ = τ).z
    Ψ(y) = basis_hamiltonians(H, y[1:n], y[(n + 1):(2n)])
    zs = _packed_stage(scheme, z, z̄)
    jacobians || return (z = z̄, θ = θ, g = Ψ(zs), jacobians = nothing)
    function R(ȳ, y)
        ys = _packed_stage(scheme, y, ȳ)
        Hq, Hp, _ = pairwise_hamiltonian_derivatives(H, ys[1:n], ys[(n + 1):(2n)], τ .*
                                                                                   Θ(y))
        vcat(ȳ[1:n] .- y[1:n] .- Hp, ȳ[(n + 1):(2n)] .- y[(n + 1):(2n)] .+ Hq)
    end
    DΦ = implicit_step_jacobian(R, z̄, z)
    DΘ = ForwardDiff.jacobian(Θ, z)
    Dstage = ForwardDiff.jacobian(y -> _packed_stage(scheme, y, z̄), z) +
             ForwardDiff.jacobian(ȳ -> _packed_stage(scheme, z, ȳ), z̄) * DΦ
    Dg = ForwardDiff.jacobian(Ψ, zs) * Dstage
    (z = z̄, θ = θ, g = Ψ(zs), jacobians = (DΦ, DΘ, Dg))
end

"""
    cotangent_lift_modes(q, p, M) -> (Ψ, ∂qΨ, ∂pΨ)

The `M` translation-invariant cotangent-lift generators
`Ψ_k(q, p) = −N⁻¹ Σ_{a,b} (p_a − p_b) sin(k (q_a − q_b))`, `k = 1, …, M` (repulsive kernels with
`∇ψ_k(Δq) = sin(k Δq)`), and their `N × M` gradient matrices.
"""
function cotangent_lift_modes(q::AbstractVector, p::AbstractVector, M::Integer)
    n = length(q)
    T = promote_type(eltype(q), eltype(p))
    vals = zeros(T, M)
    dq_mat = zeros(T, n, M)
    dp_mat = zeros(T, n, M)
    for k in 1:M
        for a in 1:n, b in 1:n

            s = sin(k * (q[a] - q[b]))
            c = cos(k * (q[a] - q[b]))
            vals[k] += -(1 / n) * (p[a] - p[b]) * s
            g_q = -(1 / n) * (p[a] - p[b]) * k * c
            dq_mat[a, k] += g_q
            dq_mat[b, k] -= g_q
            g_p = -(1 / n) * s
            dp_mat[a, k] += g_p
            dp_mat[b, k] -= g_p
        end
    end
    vals, dq_mat, dp_mat
end

# The residual of the variational closed loop: the midpoint equations of H_var(·; θ) with the
# stationarity closure θ = θ₀ + Ψ(z_{1/2}) / κ.
function _variational_residual(z̄, z, θ0, κ)
    n = length(z) ÷ 2
    zm = 0.5 .* (z .+ z̄)
    vals, dq, dp = cotangent_lift_modes(zm[1:n], zm[(n + 1):(2n)], length(θ0))
    θ = θ0 .+ vals ./ κ
    vcat(z̄[1:n] .- z[1:n] .- dp * θ, z̄[(n + 1):(2n)] .- z[(n + 1):(2n)] .+ dq * θ)
end

"""
    variational_closed_loop(z; θ0 = [0.06, 0.03, 0.015], κ = 2, tol = 1e-15, maxiter = 35,
        jacobians = true) -> (z̄, θ, jacobians)

The closed-loop map of the variational closure: one midpoint step of
`H_var(q, p; θ) = Σ_k θ_k Ψ_k(q, p)` ([`cotangent_lift_modes`](@ref)) whose parameters satisfy
the stationarity condition `∂H/∂θ = ∂V/∂θ` with `V(θ) = (κ/2) |θ − θ₀|²`, that is
`θ = θ₀ + Ψ(z_{1/2}) / κ`, solved by Newton's method. Returns the end point, the parameters and,
for `jacobians = true`, the exact Jacobians `(DΦ, DΘ)` by implicit differentiation (otherwise
`nothing`).
"""
function variational_closed_loop(z::AbstractVector; θ0 = [0.06, 0.03, 0.015], κ::Real = 2,
        tol::Real = 1e-15, maxiter::Integer = 35, jacobians::Bool = true)
    n = length(z) ÷ 2
    z̄ = copy(z)
    for _ in 1:maxiter
        F = _variational_residual(z̄, z, θ0, κ)
        norm(F) < tol && break
        z̄ .-= ForwardDiff.jacobian(y -> _variational_residual(y, z, θ0, κ), z̄) \ F
    end
    zm = 0.5 .* (z .+ z̄)
    Ψ(y) = cotangent_lift_modes(y[1:n], y[(n + 1):(2n)], length(θ0))[1]
    θ = θ0 .+ Ψ(zm) ./ κ
    jacobians || return (z = z̄, θ = θ, jacobians = nothing)
    DΦ = implicit_step_jacobian((ȳ, y) -> _variational_residual(ȳ, y, θ0, κ), z̄, z)
    DΘ = ForwardDiff.jacobian(Ψ, zm) * (0.5 .* (I + DΦ)) ./ κ
    (z = z̄, θ = θ, jacobians = (DΦ, DΘ))
end

"""
    fixed_grid_closed_loop(z, q_grid)

The closed-loop rule that resets the positions onto the state-independent grid `q_grid` and
keeps the momenta: `Φ(q, p) = (q_grid, p)`. Its pull-back of `ω` vanishes identically.
"""
function fixed_grid_closed_loop(z::AbstractVector, q_grid::AbstractVector)
    n = length(z) ÷ 2
    vcat(q_grid .+ zero(eltype(z)), z[(n + 1):(2n)])
end
