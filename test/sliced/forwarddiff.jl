# The port against the verification script it comes from. `ScriptReference` holds the script's
# hand-written forward-mode `Dual` number, its linear solve, its Jacobian, the sliced lift and the
# weight squeeze, copied verbatim. The package computes Jacobians with ForwardDiff instead; the
# values and the Jacobians must agree to round-off.

using ParticleRemapping
using ForwardDiff: ForwardDiff
using LinearAlgebra: cond, opnorm
using Random: MersenneTwister
using Test

#! format: off
module ScriptReference

using LinearAlgebra

# Exact first-order forward-mode Dual numbers (zero external dependencies)
struct Dual <: Real
    val::Float64
    der::Float64
end
Dual(x::Dual) = x
Dual(x::Real) = Dual(Float64(x), 0.0)
Base.promote_rule(::Type{Dual}, ::Type{<:Real}) = Dual
Base.convert(::Type{Dual}, x::Dual) = x
Base.convert(::Type{Dual}, x::Real) = Dual(Float64(x), 0.0)
Base.:+(a::Dual, b::Dual) = Dual(a.val + b.val, a.der + b.der)
Base.:-(a::Dual, b::Dual) = Dual(a.val - b.val, a.der - b.der)
Base.:-(a::Dual) = Dual(-a.val, -a.der)
Base.:*(a::Dual, b::Dual) = Dual(a.val * b.val, a.der * b.val + a.val * b.der)
Base.:/(a::Dual, b::Dual) = Dual(a.val / b.val, (a.der * b.val - a.val * b.der) / (b.val^2))
Base.sqrt(a::Dual) = (s = sqrt(a.val); Dual(s, 0.5 * a.der / s))
Base.exp(a::Dual) = (e = exp(a.val); Dual(e, a.der * e))
Base.log(a::Dual) = Dual(log(a.val), a.der / a.val)
Base.tanh(a::Dual) = (t = tanh(a.val); Dual(t, a.der * (1.0 - t^2)))
Base.cos(a::Dual) = Dual(cos(a.val), -a.der * sin(a.val))
Base.sin(a::Dual) = Dual(sin(a.val), a.der * cos(a.val))
Base.abs2(a::Dual) = a * a
Base.zero(::Type{Dual}) = Dual(0.0, 0.0)
Base.one(::Type{Dual}) = Dual(1.0, 0.0)
Base.conj(a::Dual) = a
Base.real(a::Dual) = a

solve_linear(A_in::Matrix{Float64}, b_in::AbstractVector{Float64}) = A_in \ b_in

function solve_linear(A_in::Matrix{T}, b_in::AbstractVector{T}) where {T}
    n = length(b_in)
    A = copy(A_in)
    b = Vector{T}(b_in)
    for k in 1:n
        piv = A[k, k]
        for i in (k + 1):n
            fac = A[i, k] / piv
            b[i] -= fac * b[k]
            for j in k:n
                A[i, j] -= fac * A[k, j]
            end
        end
    end
    x_sol = Vector{T}(undef, n)
    for i in n:-1:1
        s = b[i]
        for j in (i + 1):n
            s -= A[i, j] * x_sol[j]
        end
        x_sol[i] = s / A[i, i]
    end
    return x_sol
end

function jacobian_dual(f, z::Vector{Float64})
    y0 = f(z)
    m_out = length(y0)
    n_in = length(z)
    J = zeros(m_out, n_in)
    zd = Dual.(z, 0.0)
    for j in 1:n_in
        zd[j] = Dual(z[j], 1.0)
        yd = f(zd)
        for i in 1:m_out
            J[i, j] = yd[i].der
        end
        zd[j] = Dual(z[j], 0.0)
    end
    return y0, J
end

canonical_omega(n_dof::Int) = [zeros(n_dof, n_dof) I(n_dof); -I(n_dof) zeros(n_dof, n_dof)]

# 1D sliced strictly convex translation-equivariant quantile optimal-transport cotangent lift.
# Exact convex potential:
#   S(u) = (1 / 2N) (sum_c u_c)^2
#        + (1 - beta_mean)/2 * sum_b n_b (ubar_b - ubar)^2
#        + (1 - beta_intra)/2 * sum_b sum_{a in S_b} (u_a - ubar_b)^2
#        + sum_b (beta_intra / (2 n_b att_b)) * sum_{r,s} |tau_r - tau_s| sqrt((u_{a_r} - u_{a_s})^2 + eps^2)
function split_quantile_ot_lift(
    u::AbstractVector{T},
    v_conj::AbstractVector{T},
    bins::Vector{Vector{Int}},
    bin_taus::Vector{Vector{Float64}};
    anchors::Union{Nothing, Vector{Float64}} = nothing,
    beta_mean::Float64 = 0.60,
    beta_intra::Float64 = 0.72,
    eps_abs::Float64 = 0.4,
    with_cotangent::Bool = true,
) where {T}
    N = length(u)
    ubar = sum(u) / N
    unew = Vector{T}(undef, N)
    H = fill(T(beta_mean) / N, N, N)
    eps2 = T(eps_abs)^2
    for (b, idx) in enumerate(bins)
        nb = length(idx)
        ubar_b = sum(u[idx]) / nb
        bm_diff = T(beta_intra - beta_mean) / nb
        for ia in 1:nb, ic in 1:nb
            H[idx[ia], idx[ic]] += bm_diff
        end
        taus = bin_taus[b]
        span_tau = taus[end] - taus[1]
        r_ratio = eps_abs / max(span_tau / 2, 1e-12)
        att = sqrt(4.0 + r_ratio^2) / 2.0 - r_ratio / 2.0
        coef = T(beta_intra / att) / nb
        for r in 1:nb
            a = idx[r]
            ua = u[a]
            rep_a = zero(T)
            diag_h = zero(T)
            for s in 1:nb
                r == s && continue
                c = idx[s]
                wt = T(abs(taus[r] - taus[s]))
                diff_u = ua - u[c]
                denom = sqrt(diff_u^2 + eps2)
                rep_a += wt * diff_u / denom
                hac = -coef * wt * eps2 / (denom * denom * denom)
                H[a, c] += hac
                diag_h -= hac
            end
            unew[a] =
                T(beta_mean) * ubar +
                T(beta_intra - beta_mean) * ubar_b +
                (one(T) - T(beta_intra)) * ua +
                coef * rep_a
            H[a, a] += (one(T) - T(beta_intra)) + diag_h
        end
    end
    if !with_cotangent
        return unew, copy(v_conj)
    elseif anchors === nothing
        vnew = solve_linear(H, v_conj)
    else
        v_shift = [v_conj[a] - T(anchors[a]) for a in 1:N]
        vnew = solve_linear(H, v_shift) .+ T.(anchors)
    end
    return unew, vnew
end

# Exact O(N) Sherman-Morrison momentum-conserving canonical weight squeeze T_mu (for sum(mu) == N).
function canonical_weight_squeeze(
    qt::AbstractVector{T},
    pt::AbstractVector{T},
    mu::AbstractVector{TM},
    q_anc::Vector{Float64},
) where {T, TM}
    N = length(qt)
    TR = promote_type(T, TM)
    p_corr = zero(TR)
    q_corr = zero(TR)
    for c in 1:N
        p_corr += (one(TR) - TR(mu[c])) * TR(pt[c])
        q_corr += TR(qt[c]) - TR(q_anc[c])
    end
    p_corr /= N
    q_corr /= N
    pbar = Vector{TR}(undef, N)
    qbar = Vector{TR}(undef, N)
    for a in 1:N
        ma = TR(mu[a])
        pbar[a] = ma * (TR(pt[a]) + p_corr)
        qbar[a] = TR(q_anc[a]) + ((TR(qt[a]) - TR(q_anc[a])) - (one(TR) - ma) * q_corr) / ma
    end
    return qbar, pbar
end

end
#! format: on

const R = ScriptReference

# the plan, the multipliers and the state of the autonomous map of suite 3
ex = autonomous_map_experiment()
plan, μ = ex.plan, ex.μ_env
N = length(plan.q_anc)
rng = MersenneTwister(303)
q0 = rand(rng, N) .* 2pi .- pi
v0 = randn(rng, N)
w0 = [1.0 + 0.3 * cos(q0[a]) * exp(-v0[a]^2 / 2) for a in 1:N]
w0 ./= sum(w0)
p0 = w0 .* 1.25 .* v0
ε_p = plan.ε_p

port(q, p) = begin
    st = apply_sliced_lift(plan, q, p)
    canonical_weight_squeeze(st.q̂, st.p̂, μ, plan.q_anc)
end
function reference(q, p)
    qt, pt = R.split_quantile_ot_lift(q, p, plan.bins_v, plan.taus_x; anchors = plan.p_anc)
    pt2, qt2 = R.split_quantile_ot_lift(pt, qt, plan.bins_x, plan.taus_p;
        anchors = plan.q_anc, eps_abs = ε_p)
    R.canonical_weight_squeeze(qt2, pt2, μ, plan.q_anc)
end

@testset "values agree with the script" begin
    for (a, b) in zip(port(q0, p0), reference(q0, p0))
        @test isapprox(a, b; rtol = 8 * eps())
    end
    for kw in ((;), (; anchors = plan.p_anc))
        for (a, b) in zip(split_quantile_ot_lift(q0, p0, plan.bins_v, plan.taus_x; kw...),
            R.split_quantile_ot_lift(q0, p0, plan.bins_v, plan.taus_x; kw...))
            @test isapprox(a, b; rtol = 8 * eps())
        end
    end
end

@testset "ForwardDiff Jacobians agree with the script's Dual numbers" begin
    z0 = vcat(q0, p0)
    J_port = ForwardDiff.jacobian(z -> vcat(port(z[1:N], z[(N + 1):(2N)])...), z0)
    _, J_ref = R.jacobian_dual(z -> vcat(reference(z[1:N], z[(N + 1):(2N)])...), z0)
    # the two differ by the pivoting of the linear solves only
    @test isapprox(J_port, J_ref; atol = 8 * eps() * N * cond(J_port) * opnorm(J_port))
    lift(z) = vcat(split_quantile_ot_lift(z[1:N], z[(N + 1):(2N)], plan.bins_v, plan.taus_x)...)
    lift_ref(z) = vcat(R.split_quantile_ot_lift(z[1:N], z[(N + 1):(2N)], plan.bins_v,
        plan.taus_x)...)
    J1 = ForwardDiff.jacobian(lift, z0)
    @test isapprox(J1, R.jacobian_dual(lift_ref, z0)[2];
        atol = 8 * eps() * N * cond(J1) * opnorm(J1))
end
