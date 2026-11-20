# The linear parameter solve G(z_stage) θ = z̄ − z and the forward Newton solves of one implicit
# midpoint or Type-2 symplectic Euler step of a linear ensemble Hamiltonian (1D-1V).

"""
    RedistributionScheme

The one-step symplectic scheme of a redistribution step: [`ImplicitMidpointStep`](@ref) or
[`SymplecticEulerStep`](@ref).
"""
abstract type RedistributionScheme end

"""
    ImplicitMidpointStep()

The implicit midpoint step `z̄ − z = τ J ∇H_θ((z + z̄)/2)`, with stage point `z_{1/2}`.
"""
struct ImplicitMidpointStep <: RedistributionScheme end

"""
    SymplecticEulerStep()

The Type-2 symplectic Euler step `q̄ − q = τ ∇_p̄ H_θ(q, p̄)`, `p̄ − p = −τ ∇_q H_θ(q, p̄)`, with
stage point `(q, p̄)`.
"""
struct SymplecticEulerStep <: RedistributionScheme end

"""
    redistribution_stage_point(scheme, q, p, q̄, p̄) -> (q_stage, p_stage)

The stage point of `scheme`, known from the start `(q, p)` and the target `(q̄, p̄)`:
`((q + q̄)/2, (p + p̄)/2)` for [`ImplicitMidpointStep`](@ref) and `(q, p̄)` for
[`SymplecticEulerStep`](@ref).
"""
function redistribution_stage_point(::ImplicitMidpointStep, q, p, q̄, p̄)
    (0.5 .* (q .+ q̄), 0.5 .* (p .+ p̄))
end
redistribution_stage_point(::SymplecticEulerStep, q, p, q̄, p̄) = (q, p̄)

"""
    redistribution_parameter_matrix(H, scheme, q, p, q̄, p̄; τ = 1)

The `2N × N_θ` matrix `G` with columns `G_{:,k} = τ J ∇_z Ψ_k(z_stage)`, so that one step of
`scheme` from `(q, p)` reaches `(q̄, p̄)` exactly when `G θ = z̄ − z`. Each column lies in the
zero-net-momentum subspace `M₀`.
"""
function redistribution_parameter_matrix(
        H, scheme::RedistributionScheme, q, p, q̄, p̄; τ::Real = 1)
    qs, ps = redistribution_stage_point(scheme, q, p, q̄, p̄)
    dq, dp = basis_hamiltonian_gradients(H, qs, ps)
    τ .* vcat(dp, -dq)
end

"""
    redistribution_scaling(N, p_scale)

The diagonal of the scaling `S = diag(I_N, p_s⁻¹ I_N)` that makes positions and momenta
dimensionless.
"""
function redistribution_scaling(N::Integer, p_scale::Real)
    vcat(ones(typeof(float(p_scale)), N),
        fill(1 / p_scale, N))
end

"""
    solve_redistribution_parameters(H, scheme, q, p, q̄, p̄; τ = 1, λ = 0)

The parameters `θ` of the linear system `G(z_stage) θ = z̄ − z` ([`redistribution_parameter_matrix`](@ref)) in
the scaled form `(S G) θ = S (z̄ − z)`: the minimum-norm least-squares solution for `λ = 0`, the
Tikhonov-regularised `(AᵀA + λI)⁻¹ Aᵀ S (z̄ − z)` with `A = S G` for `λ > 0`.
"""
function solve_redistribution_parameters(H, scheme::RedistributionScheme, q, p, q̄, p̄;
        τ::Real = 1, λ::Real = 0)
    S = redistribution_scaling(length(q), H.p_scale)
    A = Diagonal(S) * redistribution_parameter_matrix(H, scheme, q, p, q̄, p̄; τ = τ)
    b = S .* vcat(q̄ .- q, p̄ .- p)
    λ == 0 && return A \ b
    (transpose(A) * A + λ * I) \ (transpose(A) * b)
end

"""
    redistribution_contraction_number(H, scheme, q, p, q̄, p̄, θ; τ = 1)

The quantity whose being smaller than one makes the forward step invertible on the
identity-connected branch: `(τ/2) ‖S J ∇²H_θ(z_{1/2}) S⁻¹‖` for [`ImplicitMidpointStep`](@ref),
`‖τ ∇_q ∇_p̄ H_θ(q, p̄)‖` for [`SymplecticEulerStep`](@ref) (spectral norms).
"""
function redistribution_contraction_number(
        H, ::ImplicitMidpointStep, q, p, q̄, p̄, θ; τ::Real = 1)
    n = length(q)
    S = redistribution_scaling(n, H.p_scale)
    qs, ps = redistribution_stage_point(ImplicitMidpointStep(), q, p, q̄, p̄)
    _, _, Hzz = pairwise_hamiltonian_derivatives(H, qs, ps, τ .* θ)
    opnorm(Diagonal(S) * canonical_J(n) * Hzz * Diagonal(1 ./ S)) / 2
end

function redistribution_contraction_number(
        H, ::SymplecticEulerStep, q, p, q̄, p̄, θ; τ::Real = 1)
    n = length(q)
    _, _, Hzz = pairwise_hamiltonian_derivatives(H, q, p̄, τ .* θ)
    opnorm(Hzz[1:n, (n + 1):(2n)])
end

default_homotopy(::ImplicitMidpointStep) = (0.25, 0.5, 0.75, 1.0)
default_homotopy(::SymplecticEulerStep) = (0.2, 0.4, 0.6, 0.8, 1.0)

"""
    forward_redistribution_step(H, scheme, q, p, θ; τ = 1, tol = 1e-14, maxiter = 30,
        homotopy) -> (z, jacobian, residual, history)

One step of `scheme` with the fixed parameters `θ` from `(q, p)`, solved by Newton's method
with a homotopy `θ → s θ` over the values `s` of `homotopy` (the default is `(0.25, 0.5, 0.75,
1)` for the midpoint and `(0.2, 0.4, …, 1)` for symplectic Euler), which follows the
identity-connected branch. Returns the packed end point `z`, the Jacobian `DΨ` of the step by
the implicit function theorem (the Cayley transform for the midpoint), the scaled residual at
the end point and the `history` of `(s, residual)` at every Newton iterate.
"""
function forward_redistribution_step(H, scheme::ImplicitMidpointStep, q::AbstractVector,
        p::AbstractVector, θ::AbstractVector; τ::Real = 1, tol::Real = 1e-14,
        maxiter::Integer = 30, homotopy = default_homotopy(scheme))
    n = length(q)
    S_row = redistribution_scaling(n, H.p_scale)
    J = canonical_J(n)
    z = vcat(q, p)
    history = Tuple{Float64, Float64}[]
    for s in homotopy
        θs = (τ * s) .* θ
        for _ in 1:maxiter
            qm = 0.5 .* (q .+ z[1:n])
            pm = 0.5 .* (p .+ z[(n + 1):(2n)])
            Hq, Hp, Hzz = pairwise_hamiltonian_derivatives(H, qm, pm, θs)
            F = vcat(z[1:n] .- q .- Hp, z[(n + 1):(2n)] .- p .+ Hq)
            r = norm(S_row .* F)
            push!(history, (s, r))
            r < tol && break
            z .-= (I - 0.5 .* (J * Hzz)) \ F
        end
    end
    qm = 0.5 .* (q .+ z[1:n])
    pm = 0.5 .* (p .+ z[(n + 1):(2n)])
    Hq, Hp, Hzz = pairwise_hamiltonian_derivatives(H, qm, pm, τ .* θ)
    F = vcat(z[1:n] .- q .- Hp, z[(n + 1):(2n)] .- p .+ Hq)
    DΨ = (I - 0.5 .* (J * Hzz)) \ (I + 0.5 .* (J * Hzz))
    (z = z, jacobian = DΨ, residual = norm(S_row .* F), history = history)
end

function forward_redistribution_step(H, scheme::SymplecticEulerStep, q::AbstractVector,
        p::AbstractVector, θ::AbstractVector; τ::Real = 1, tol::Real = 1e-14,
        maxiter::Integer = 30, homotopy = default_homotopy(scheme))
    n = length(q)
    p̄ = copy(p)
    history = Tuple{Float64, Float64}[]
    for s in homotopy
        θs = (τ * s) .* θ
        for _ in 1:maxiter
            Hq, _, Hzz = pairwise_hamiltonian_derivatives(H, q, p̄, θs)
            Fp = p̄ .- p .+ Hq
            r = norm(Fp) / H.p_scale
            push!(history, (s, r))
            r < tol && break
            p̄ .-= (I + Hzz[1:n, (n + 1):(2n)]) \ Fp
        end
    end
    Hq, Hp, Hzz = pairwise_hamiltonian_derivatives(H, q, p̄, τ .* θ)
    q̄ = q .+ Hp
    Fp = p̄ .- p .+ Hq
    Hqq = Hzz[1:n, 1:n]
    Hqp = Hzz[1:n, (n + 1):(2n)]
    Hpp = Hzz[(n + 1):(2n), (n + 1):(2n)]
    M⁻¹ = inv(I + Hqp)
    dp̄_dq = -M⁻¹ * Hqq
    dp̄_dp = M⁻¹
    dq̄_dq = I + transpose(Hqp) + Hpp * dp̄_dq
    dq̄_dp = Hpp * dp̄_dp
    (z = vcat(q̄, p̄), jacobian = [dq̄_dq dq̄_dp; dp̄_dq dp̄_dp],
        residual = norm(Fp) / H.p_scale,
        history = history)
end

"""
    redistribution_substeps(H, scheme, q, p, q̄, p̄; substeps = 1, τ = 1, kwargs...)
        -> (z, θ, contraction, residual, history)

The composition of `substeps` redistribution steps towards the targets
`z_k = z + (k/K)(z̄ − z)`: each solves [`solve_redistribution_parameters`](@ref) from the end
point of the previous step and advances by [`forward_redistribution_step`](@ref). Every step
is a symplectic map at fixed parameters, and so is the composition. Returns the end point, the
parameters, the [`redistribution_contraction_number`](@ref) and the final residual of every step, and the
Newton histories.
"""
function redistribution_substeps(H, scheme::RedistributionScheme, q, p, q̄, p̄;
        substeps::Integer = 1, τ::Real = 1, kwargs...)
    n = length(q)
    z0, Δ = vcat(q, p), vcat(q̄ .- q, p̄ .- p)
    z = copy(z0)
    θs, κs, rs, hs = Vector{Vector{Float64}}(), Float64[], Float64[], Any[]
    for k in 1:substeps
        zt = z0 .+ (k / substeps) .* Δ
        qk, pk, qt, pt = z[1:n], z[(n + 1):(2n)], zt[1:n], zt[(n + 1):(2n)]
        θ = solve_redistribution_parameters(H, scheme, qk, pk, qt, pt; τ = τ)
        push!(κs, redistribution_contraction_number(H, scheme, qk, pk, qt, pt, θ; τ = τ))
        step = forward_redistribution_step(H, scheme, qk, pk, θ; τ = τ, kwargs...)
        z = step.z
        push!(θs, θ)
        push!(rs, step.residual)
        push!(hs, step.history)
    end
    (z = z, θ = θs, contraction = κs, residual = rs, history = hs)
end

"""
    pseudotime_path(H, scheme, q, p, θ; s = range(0, 1; length = 41), τ = 1)

The end points `Ψ_{s τ}^θ(q, p)` of the step of `scheme` with the fixed parameters `θ` and the
artificial step `s τ`, for each `s`: a path in pseudo-time from `(q, p)` (`s = 0`) to the step's
end point (`s = 1`). Returns a vector of packed states.
"""
function pseudotime_path(H, scheme::RedistributionScheme, q, p, θ;
        s = range(0, 1; length = 41), τ::Real = 1)
    [iszero(si) ? vcat(q, p) :
     forward_redistribution_step(H, scheme, q, p, θ; τ = si * τ).z
     for si in s]
end
