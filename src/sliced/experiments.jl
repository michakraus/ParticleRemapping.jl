# The set-ups of the verification of the sliced redistribution (five suites) and the probes of
# the constructions that do not work. Each function takes the parameters of its run as keyword
# arguments, with the values of the verification as defaults, and returns a NamedTuple.

# The largest entry of DFᵀ Ω DF − Ω_ref with Ω the canonical matrix of the packed (q, p).
function _canonical_defect(DF, Ω_ref)
    maximum(abs, transpose(DF) * canonical_J(size(DF, 1) ÷ 2) * DF -
                 Ω_ref)
end

"""
    factorisation_pullback_experiment(; seed = 101, N = 8, m = 1.35)

The factorisation on the full state space `(x, v, w)` (suite 1). The state `z = (x, v, w)` of `N`
random markers is mapped to the canonical output `(q̄, p̄ = m w̄ v̄)` by three maps, and the pull-back
of `Σ dq̄ ∧ dp̄` is compared with that of `Σ dq ∧ dp` (the presymplectic form on `(x, v, w)`):

- `full`: a [`split_quantile_ot_lift`](@ref) of `(q, p) = (x, m w v)` (slices `1:4`, `5:8`,
  targets `±0.5, ±1.5`, `β_mean = 0.5`, `β_intra = 0.7`), new weights
  `w̄ ∝ w (1 + 0.3 cos q̄ e^{−v²/4})` with `Σ w̄ = Σ w`, and `v̄ = p̄ / (m w̄)`;
- `no_recoil`: the same with `v̄ = p̄ / (m w)` (a control that must fail);
- `leaf`: `q̄ = q`, `p̄ = p − w sin q`, a map that reads `w` at fixed `(q, p)`: canonical on each
  leaf `w = const`, not across leaves (a control that must fail across `dw`).

Returns the largest entries `err_full`, `err_no_recoil`, `err_on_leaf` (variations of `x`, `v`
only) and `err_across_dw` (all variations).
"""
function factorisation_pullback_experiment(; seed::Integer = 101, N::Integer = 8,
        m::Real = 1.35)
    rng = MersenneTwister(seed)
    x0 = randn(rng, N)
    v0 = randn(rng, N)
    w0 = 0.5 .+ rand(rng, N)
    z = vcat(x0, v0, w0)
    Π = ForwardDiff.jacobian(u -> vcat(u[1:N], u[(2N + 1):(3N)] .* m .* u[(N + 1):(2N)]), z)
    Ω_xvw = transpose(Π) * canonical_J(N) * Π
    bins = [collect(1:4), collect(5:8)]
    taus = [[-1.5, -0.5, 0.5, 1.5], [-1.5, -0.5, 0.5, 1.5]]
    function reweighted(u; recoil::Bool)
        T = eltype(u)
        xv, vv, wv = u[1:N], u[(N + 1):(2N)], u[(2N + 1):(3N)]
        qbar, pbar = split_quantile_ot_lift(xv, wv .* T(m) .* vv, bins, taus; β_mean = 0.5,
            β_intra = 0.7)
        raw_w = [wv[a] * (one(T) + T(0.3) * cos(qbar[a]) * exp(-vv[a]^2 / 4)) for a in 1:N]
        wbar = (sum(wv) / sum(raw_w)) .* raw_w
        vbar = [pbar[a] / (T(m) * (recoil ? wbar[a] : wv[a])) for a in 1:N]
        vcat(qbar, wbar .* T(m) .* vbar)
    end
    function leaf_only(u)
        xv, vv, wv = u[1:N], u[(N + 1):(2N)], u[(2N + 1):(3N)]
        pv = wv .* eltype(u)(m) .* vv
        vcat(copy(xv), [pv[a] - wv[a] * sin(xv[a]) for a in 1:N])
    end
    J_full = ForwardDiff.jacobian(u -> reweighted(u; recoil = true), z)
    J_bad = ForwardDiff.jacobian(u -> reweighted(u; recoil = false), z)
    J_leaf = ForwardDiff.jacobian(leaf_only, z)
    J_xv = J_leaf[:, 1:(2N)]
    (err_full = _canonical_defect(J_full, Ω_xvw),
        err_no_recoil = _canonical_defect(J_bad, Ω_xvw),
        err_on_leaf = _canonical_defect(J_xv, Ω_xvw[1:(2N), 1:(2N)]),
        err_across_dw = _canonical_defect(J_leaf, Ω_xvw))
end

"""
    pseudo_time_flow_experiment(; seed = 202, N = 6, m = 1.1, Δτ = 0.12, iterations = 25,
                                rate = 0.25)

The pseudo-time flow with a weight rate (suite 2). One implicit-midpoint step `Δτ` of the
generator `G(q, p) = (1/4) Σ_{a,b} exp(−(q_a − q_b)²/2) (p_a − p_b)²` (fixed-point iteration,
`iterations` sweeps) moves `(q, p)`; the weights change as `w̄_a = w_a exp(Δτ · rate · cos q̄_a)`.
With the recoil the velocities are `v̄ = p̄ / (m w̄)`; without it `v̄ = p̄ / (m w)`. Returns the
symplecticity defects `err_recoil` and `err_no_recoil` of `(q, p) ↦ (q̄, m w̄ v̄)`.
"""
function pseudo_time_flow_experiment(; seed::Integer = 202, N::Integer = 6, m::Real = 1.1,
        Δτ::Real = 0.12, iterations::Integer = 25, rate::Real = 0.25)
    rng = MersenneTwister(seed)
    x0 = randn(rng, N)
    v0 = randn(rng, N)
    w0 = 0.8 .+ 0.4 .* rand(rng, N)
    function flow_step(z::AbstractVector{T}; recoil::Bool) where {T}
        q, p = z[1:N], z[(N + 1):(2N)]
        dτ = T(Δτ)
        qn, pn = copy(q), copy(p)
        for _ in 1:iterations
            qm, pm = (q .+ qn) ./ 2, (p .+ pn) ./ 2
            dG_dq, dG_dp = zeros(T, N), zeros(T, N)
            for a in 1:N, b in 1:N

                a == b && continue
                dq, dp = qm[a] - qm[b], pm[a] - pm[b]
                k = exp(-dq^2 / 2)
                dG_dp[a] += k * dp
                dG_dq[a] += (-dq * k) * (dp^2) / 2
            end
            qn = q .+ dτ .* dG_dp
            pn = p .- dτ .* dG_dq
        end
        wbar = [T(w0[a]) * exp(dτ * T(rate) * cos(qn[a])) for a in 1:N]
        vbar = [pn[a] / (T(m) * (recoil ? wbar[a] : T(w0[a]))) for a in 1:N]
        vcat(qn, wbar .* T(m) .* vbar)
    end
    z0 = vcat(x0, w0 .* m .* v0)
    Ω = canonical_J(N)
    (
        err_recoil = _canonical_defect(
            ForwardDiff.jacobian(z -> flow_step(z; recoil = true), z0),
            Ω),
        err_no_recoil = _canonical_defect(
            ForwardDiff.jacobian(z -> flow_step(z; recoil = false), z0), Ω))
end

# The final reweighting of the autonomous map: pilot weights w_ref μ_env, a kernel ratio at the
# squeezed points limited to a factor exp(±limit) by a tanh, normalised to Σ w̄ = Σ w.
function _limited_ratio_weights(qbar, pbar, q, v_in, w0, pilot_w, m, limit)
    T = promote_type(eltype(qbar), eltype(pbar), eltype(q), eltype(v_in))
    N = length(qbar)
    v_pil = [pbar[a] / (T(m) * pilot_w[a]) for a in 1:N]
    raw_w = Vector{T}(undef, N)
    for a in 1:N
        fn = zero(T)
        rd = zero(T)
        for c in 1:N
            fn += T(w0[c]) * exp(-((qbar[a] - q[c])^2 + (v_pil[a] - v_in[c])^2))
            rd += exp(-((qbar[a] - qbar[c])^2 + (v_pil[a] - v_pil[c])^2))
        end
        ratio = (fn / (rd + T(1e-12))) / pilot_w[a]
        raw_w[a] = pilot_w[a] * exp(T(limit) * tanh(log(ratio + T(1e-12)) / T(limit)))
    end
    (T(sum(w0)) / sum(raw_w)) .* raw_w
end

"""
    autonomous_map_experiment(; seed = 303, N = 16, m = 1.25, Lx = π, Vmax = 2.0, B_v = 4,
                              B_x = 4, N2 = 12)

The autonomous four-stage map in `d = 1` and `d = 2` with its controls (suite 3).

In `d = 1`, `N` random markers on `[−π, π) × ℝ` define a frozen [`SlicedLiftPlan`](@ref)
(`B_v = B_x = 4`, `ε_q = 0.4`, `ε_p = 1.5 · 2Lp / n_x`), slot multipliers
`μ_env ∝ 0.25 + exp(−v_slot²/2.5)` that depend on the rank only, and the final weights of a
kernel ratio limited to `w_ref μ_env exp(±0.25)`. The map `(q, p) ↦ (q̄, m w̄ v̄)` is checked for
symplecticity (`err_symp_1d`) and momentum (`err_ptot_1d`), with two controls that must fail:
the lifts without their cotangent parts (`err_nocot`) and a squeeze with state-dependent
multipliers `μ ∝ 0.3 + exp(−q̄²/4 − (p̄/Lp)²/2)` (`err_badmu`).

In `d = 2`, `N2` markers in 2D-2V: the [`per_axis_sliced_lift`](@ref) with the fixed index slices
`1:6`, `7:12` ordered by each position component, the squeeze of each axis with one set of
multipliers from `|v|²`, and weights `∝ μ (1 + 0.15 cos(q̄₁ + q̄₂))`: `err_symp_2d` and the
momentum error `err_pvec_2d`.
"""
function autonomous_map_experiment(; seed::Integer = 303, N::Integer = 16, m::Real = 1.25,
        Lx::Real = Float64(π), Vmax::Real = 2.0, B_v::Integer = 4, B_x::Integer = 4,
        N2::Integer = 12)
    rng = MersenneTwister(seed)
    q0 = rand(rng, N) .* 2pi .- pi
    v0 = randn(rng, N)
    w0 = [1.0 + 0.3 * cos(q0[a]) * exp(-v0[a]^2 / 2) for a in 1:N]
    w0 ./= sum(w0)
    wref = sum(w0) / N
    nb = N ÷ B_x
    Lp = maximum(abs, [m * wref * vj for vj in uniform_quantile_targets(Vmax, nb)])
    plan = sliced_lift_plan(
        q0, v0, w0 .* m .* v0; m = m, w_ref = wref, Lx = Lx, Vmax = Vmax,
        B_v = B_v, B_x = B_x, ε_p = 1.5 * (2Lp / nb))
    mu_raw = [0.25 + exp(-plan.slot_v[a]^2 / 2.5) for a in 1:N]
    mu_env = (N / sum(mu_raw)) .* mu_raw
    function map_1d(z::AbstractVector{T}; with_cotangent::Bool = true,
            state_mu::Bool = false) where {T}
        q, p = z[1:N], z[(N + 1):(2N)]
        v_in = [p[a] / (T(m) * T(w0[a])) for a in 1:N]
        st = apply_sliced_lift(plan, q, p; with_cotangent = with_cotangent)
        μ = if state_mu
            mu_st = [T(0.3) + exp(-st.q̂[a]^2 / 4 - (st.p̂[a] / T(Lp))^2 / 2) for a in 1:N]
            (T(N) / sum(mu_st)) .* mu_st
        else
            mu_env
        end
        qbar, pbar = canonical_weight_squeeze(st.q̂, st.p̂, μ, plan.q_anc)
        pilot_w = [T(wref * mu_env[a]) for a in 1:N]
        wbar = _limited_ratio_weights(qbar, pbar, q, v_in, w0, pilot_w, m, 0.25)
        vbar = [pbar[a] / (T(m) * wbar[a]) for a in 1:N]
        vcat(qbar, wbar .* T(m) .* vbar)
    end
    z0 = vcat(q0, w0 .* m .* v0)
    Ω = canonical_J(N)
    zbar = map_1d(z0)
    J_1d = ForwardDiff.jacobian(map_1d, z0)
    J_nocot = ForwardDiff.jacobian(z -> map_1d(z; with_cotangent = false), z0)
    J_badmu = ForwardDiff.jacobian(z -> map_1d(z; state_mu = true), z0)

    X0 = randn(rng, 2, N2)
    V0 = randn(rng, 2, N2)
    W0 = 0.5 .+ rand(rng, N2)
    W0 ./= sum(W0)
    wref2 = 1.0 / N2
    halves = [collect(1:(N2 ÷ 2)), collect((N2 ÷ 2 + 1):N2)]
    bins2 = [[idx[sortperm(X0[ℓ, idx])] for idx in halves] for ℓ in 1:2]
    taus2 = [[-1.5, -0.9, -0.3, 0.3, 0.9, 1.5] for _ in 1:2]
    env2 = 0.3 .+ exp.(-vec(sum(abs2, V0; dims = 1)) ./ 3)
    mu2 = (N2 / sum(env2)) .* env2
    anc2 = zeros(N2)
    function map_2d(z::AbstractVector{T}) where {T}
        X = permutedims(reshape(z[1:(2N2)], N2, 2))
        P = permutedims(reshape(z[(2N2 + 1):(4N2)], N2, 2))
        Xt, Pt = per_axis_sliced_lift(X, P, bins2, taus2)
        x1b, p1b = canonical_weight_squeeze(Xt[1, :], Pt[1, :], mu2, anc2)
        x2b, p2b = canonical_weight_squeeze(Xt[2, :], Pt[2, :], mu2, anc2)
        raw_w = [T(wref2 * mu2[a]) * (one(T) + T(0.15) * cos(x1b[a] + x2b[a]))
                 for a in 1:N2]
        wbar = (one(T) / sum(raw_w)) .* raw_w
        v1b = [p1b[a] / (T(m) * wbar[a]) for a in 1:N2]
        v2b = [p2b[a] / (T(m) * wbar[a]) for a in 1:N2]
        vcat(x1b, x2b, wbar .* T(m) .* v1b, wbar .* T(m) .* v2b)
    end
    z2 = vcat(X0[1, :], X0[2, :], W0 .* m .* V0[1, :], W0 .* m .* V0[2, :])
    zbar2 = map_2d(z2)
    J_2d = ForwardDiff.jacobian(map_2d, z2)
    psum(z, r) = sum(z[r])
    (err_symp_1d = _canonical_defect(J_1d, Ω),
        err_ptot_1d = abs(psum(zbar, (N + 1):(2N)) - psum(z0, (N + 1):(2N))),
        err_nocot = _canonical_defect(J_nocot, Ω), err_badmu = _canonical_defect(J_badmu, Ω),
        err_symp_2d = _canonical_defect(J_2d, canonical_J(2N2)),
        err_pvec_2d = max(abs(psum(zbar2, (2N2 + 1):(3N2)) - psum(z2, (2N2 + 1):(3N2))),
            abs(psum(zbar2, (3N2 + 1):(4N2)) - psum(z2, (3N2 + 1):(4N2)))),
        plan, μ_env = mu_env)
end

"""
    composition_defect_experiment(; m = 1.4, Δt = 0.7, z0 = [0.6, -0.3, 0.9, -0.5])

The composition of a redistribution with a later flow step (suite 4), for two markers in 1D. A
canonical shear `q̄ = q`, `p̄_a = p_a − cos q_a` with state-dependent weights
`w̄ = (1 + 0.4 q̄₁² + 0.2 p̄₂, 1.2 + 0.3 q̄₂ p̄₁)` is followed by the free flow of
`H = Σ p²/(2 m w̄)` for the time `Δt` at frozen `w̄`. Returns

- `defect`: the largest entry of the symplecticity defect of the composition,
- `err_formula`: its distance to `−Σ_a d(Δt ∂H/∂w̄_a) ∧ dw̄_a`,
- `err_const_w`: the defect with constant weights `w̄ = 1` (zero),
- `err_extended` and `err_extended_wrong_sign`: the defects of the flow on the extended space
  `(q, p, w̄, θ)` with the form `Σ dq ∧ dp + Σ dw̄ ∧ dθ` and `θ̇ = −∂H/∂w̄` (zero) or
  `θ̇ = +∂H/∂w̄` (a control that must fail).
"""
function composition_defect_experiment(; m::Real = 1.4, Δt::Real = 0.7,
        z0::AbstractVector = [0.6, -0.3, 0.9, -0.5])
    wbar_rule(q, p) = [1.0 + 0.4 * q[1]^2 + 0.2 * p[2], 1.2 + 0.3 * q[2] * p[1]]
    shear(z) = [z[3] - cos(z[1]), z[4] - cos(z[2])]
    function composed(z; constant_weights::Bool)
        T = eltype(z)
        qbar, pbar = z[1:2], shear(z)
        wb = constant_weights ? [one(T), one(T)] : wbar_rule(qbar, pbar)
        vb = pbar ./ (T(m) .* wb)
        vcat(qbar .+ T(Δt) .* vb, pbar)
    end
    Ω = canonical_J(2)
    J_var = ForwardDiff.jacobian(z -> composed(z; constant_weights = false), z0)
    J_cst = ForwardDiff.jacobian(z -> composed(z; constant_weights = true), z0)
    defect = transpose(J_var) * Ω * J_var - Ω
    J_wb = ForwardDiff.jacobian(z -> wbar_rule(z[1:2], shear(z)), z0)
    J_A = ForwardDiff.jacobian(
        z -> begin
            pbar = shear(z)
            wb = wbar_rule(z[1:2], pbar)
            [-eltype(z)(Δt) * pbar[a]^2 / (2 * eltype(z)(m) * wb[a]^2) for a in 1:2]
        end,
        z0)
    formula = zeros(4, 4)
    for a in 1:2
        dA, dw = J_A[a, :], J_wb[a, :]
        formula .-= (dA * transpose(dw) .- dw * transpose(dA))
    end
    # the extended flow on (q₁, q₂, w̄₁, w̄₂, p₁, p₂, θ₁, θ₂): (q, w̄) and (p, θ) are conjugate
    function extended(y; sign::Int)
        T = eltype(y)
        q, w, p, θ = y[1:2], y[3:4], y[5:6], y[7:8]
        vcat(q .+ T(Δt) .* p ./ (T(m) .* w), w, p,
            θ .+ sign .* (T(Δt) .* (-p .^ 2 ./ (2 .* T(m) .* w .^ 2))))
    end
    y0 = [z0[1], z0[2], 1.3, 0.8, z0[3], z0[4], 0.2, -0.4]
    Ω4 = canonical_J(4)
    (defect = maximum(abs, defect), err_formula = maximum(abs, defect - formula),
        err_const_w = _canonical_defect(J_cst, Ω),
        err_extended = _canonical_defect(
            ForwardDiff.jacobian(y -> extended(y; sign = -1), y0),
            Ω4),
        err_extended_wrong_sign = _canonical_defect(
            ForwardDiff.jacobian(y -> extended(y; sign = 1), y0), Ω4))
end

"""
    checkerboard_experiment(; seed = 505, N = 400, m = 1.0, B_v = 10, B_x = 20, β_mean = 0.8,
                            β_intra = 0.85, ε_x = 0.15, ε_v = 0.15, ε_control = 0.05)

Hole filling on a checkerboard (suite 5a). `N` equal-weight markers (`w = 1/N`) fill the quadrants
`Q1 = [0, 1)²` and `Q3 = [−1, 0)²` of `[−1, 1]²`; `Q2` and `Q4` are empty, and both marginals are
uniform. Two maps, both made of the lifts `Ψ_p ∘ Ψ_q` only (the weights do not change):

- the sliced Knothe–Rosenblatt lifts of a [`sliced_lift_plan`](@ref) with `B_v` velocity and
  `B_x` position slices on `[−1, 1]²` (`ε_q = ε_x`, `ε_p = ε_v m w`);
- the control: one slice for each lift (marginal quantile lifts; `ε_q = ε_control`,
  `ε_p = ε_control m w`), which must leave `Q2` and `Q4` nearly empty.

Returns the markers `before`, `sliced`, `marginal`, their quadrant counts and the `plan`.
"""
function checkerboard_experiment(; seed::Integer = 505, N::Integer = 400, m::Real = 1.0,
        B_v::Integer = 10, B_x::Integer = 20, β_mean::Real = 0.8, β_intra::Real = 0.85,
        ε_x::Real = 0.15, ε_v::Real = 0.15, ε_control::Real = 0.05)
    rng = MersenneTwister(seed)
    qc = vcat(rand(rng, N ÷ 2), -rand(rng, N ÷ 2))
    vc = vcat(rand(rng, N ÷ 2), -rand(rng, N ÷ 2))
    wc = fill(1.0 / N, N)
    pc = wc .* m .* vc
    w_ref = 1.0 / N
    # the control: one slice in each lift, targets on all N markers
    qm1, pm1 = split_quantile_ot_lift(qc, pc, [sortperm(qc)],
        [uniform_quantile_targets(1.0, N)]; β_mean = β_mean, β_intra = β_intra, ε = ε_control)
    qm2, pm2 = dual_quantile_ot_lift(qm1, pm1, [sortperm(pm1)],
        [[m * w_ref * vj for vj in uniform_quantile_targets(1.0, N)]]; β_mean = β_mean,
        β_intra = β_intra, ε = ε_control / N)
    plan = sliced_lift_plan(
        qc, vc, pc; m = m, w_ref = w_ref, Lx = 1.0, Vmax = 1.0, B_v = B_v,
        B_x = B_x, β_mean_q = β_mean, β_intra_q = β_intra, ε_q = ε_x, β_mean_p = β_mean,
        β_intra_p = β_intra, ε_p = ε_v / N)
    st = apply_sliced_lift(plan, qc, pc)
    before = Markers(qc, vc, wc; m = m)
    sliced = Markers(st.q̂, st.p̂ ./ (m .* wc), wc; m = m)
    marginal = Markers(qm2, pm2 ./ (m .* wc), wc; m = m)
    (; before, sliced, marginal, quadrants_before = phase_quadrant_counts(qc, vc),
        quadrants_sliced = phase_quadrant_counts(sliced.x, sliced.v),
        quadrants_marginal = phase_quadrant_counts(marginal.x, marginal.v), plan, stages = st)
end

"""
    electric_wave_markers(; Nx = 30, Nv = 30, m = 1.2, Lx = π, V0 = 1.4)

The initial markers of the long-time run: an `Nx × Nv` lattice of cell centres of
`[−Lx, Lx) × [−V0, V0]` with weights `∝ f₀(x, v) Δx Δv`,
`f₀ = (1 + 0.35 cos x)(0.40 + v²) exp(−v²/1.4)`, normalised to `Σ w = 1`.
"""
function electric_wave_markers(; Nx::Integer = 30, Nv::Integer = 30, m::Real = 1.2,
        Lx::Real = Float64(π), V0::Real = 1.4)
    dx0 = 2Lx / Nx
    dv0 = 2V0 / Nv
    f0(xx, vv) = (1.0 + 0.35 * cos(xx)) * (0.40 + vv^2) * exp(-vv^2 / 1.4)
    x, v, w = Float64[], Float64[], Float64[]
    for ix in 1:Nx, iv in 1:Nv

        xa = -Lx + (ix - 0.5) * dx0
        va = -V0 + (iv - 0.5) * dv0
        push!(x, xa)
        push!(v, va)
        push!(w, f0(xa, va) * dx0 * dv0)
    end
    w ./= sum(w)
    Markers(x, v, w; m = m)
end

"""
    electric_wave_simulation(markers; Lx = π, Δt = 0.05, steps = 500, record_every = 0,
                             remap_method = nothing, remap_every = 0)

The markers pushed through the two-wave field `E(t, x) = −0.85 sin x − 0.75 sin(2x − 1.3 t)` by
`steps` drift–kick–drift steps of size `Δt` (positions wrapped into `[−Lx, Lx)`, weights
unchanged). With `remap_method` and `remap_every > 0`, the markers are remapped by
`remap(remap_method, ·)` after every `remap_every` steps. With `record_every > 0` the states after
every `record_every` steps are recorded. Returns `(; markers, times, snapshots, remap_times)`, with
`snapshots` a vector of `Markers` (the initial state first).
"""
function electric_wave_simulation(markers::Markers; Lx::Real = Float64(π), Δt::Real = 0.05,
        steps::Integer = 500, record_every::Integer = 0, remap_method = nothing,
        remap_every::Integer = 0)
    E(t, xx) = -0.85 * sin(xx) - 0.75 * sin(2xx - 1.3 * t)
    m = markers.m
    x, v, w = copy(markers.x), copy(markers.v), copy(markers.w)
    times, snapshots, remap_times = [0.0], [Markers(copy(x), copy(v), copy(w); m = m)],
    Float64[]
    for step in 1:steps
        tm = (step - 0.5) * Δt
        for a in eachindex(x)
            xh = x[a] + 0.5 * Δt * v[a]
            vn = v[a] + Δt * E(tm, xh) / m
            x[a] = periodic_wrap(xh + 0.5 * Δt * vn, Lx)
            v[a] = vn
        end
        if remap_method !== nothing && remap_every > 0 && step % remap_every == 0 &&
           step < steps
            r = remap(remap_method, Markers(x, v, w; m = m)).markers
            x, v, w = periodic_wrap.(r.x, Lx), copy(r.v), copy(r.w)
            push!(remap_times, step * Δt)
        end
        if record_every > 0 && step % record_every == 0
            push!(times, step * Δt)
            push!(snapshots, Markers(copy(x), copy(v), copy(w); m = m))
        end
    end
    (; markers = Markers(x, v, w; m = m), times, snapshots, remap_times)
end

"""
    long_time_remap_experiment(; Nx = 30, Nv = 30, m = 1.2, Lx = π, V0 = 1.4, Δt = 0.05,
                               steps = 500, method = SlicedQuantileRemap(), squeeze = true)

The long-time run (suite 5b): the [`electric_wave_markers`](@ref) pushed by
[`electric_wave_simulation`](@ref) to `T = steps Δt = 25`, then one [`remap`](@ref) with
`method`. Returns the markers `before` and `after` the remap, the cell statistics of
[`phase_cell_statistics`](@ref) (`stats_before`, `stats_after`), the
[`smoothed_density_difference`](@ref) `l2_difference`, the errors of total momentum and weight,
the remap `info`, and `min_hessian_q`, `min_hessian_p`, the smallest eigenvalues of the Hessians
of the two lifts.

With `squeeze = false` (a control) the squeeze is skipped: the weights `w_ref μ` of the density
ratio are set after the lifts at fixed momenta `p̂` (`v̄ = p̂ / (m w̄)`).
"""
function long_time_remap_experiment(; Nx::Integer = 30, Nv::Integer = 30, m::Real = 1.2,
        Lx::Real = Float64(π), V0::Real = 1.4, Δt::Real = 0.05, steps::Integer = 500,
        method::SlicedQuantileRemap = SlicedQuantileRemap(), squeeze::Bool = true)
    initial = electric_wave_markers(; Nx = Nx, Nv = Nv, m = m, Lx = Lx, V0 = V0)
    before = electric_wave_simulation(initial; Lx = Lx, Δt = Δt, steps = steps).markers
    result = remap(method, before)
    info = result.info
    after = if squeeze
        result.markers
    else
        w̄ = info.w_ref .* info.μ
        Markers(info.q̂, info.p̂ ./ (m .* w̄), w̄; m = m)
    end
    stats(mk) = phase_cell_statistics(mk.x, mk.v; Lx = method.Lx, Vmax = method.Vmax)
    plan = info.plan
    p = momentum(before)
    hq = sliced_quantile_hessian(
        before.x, plan.bins_v, plan.taus_x; β_mean = plan.β_mean_q,
        β_intra = plan.β_intra_q, ε = plan.ε_q)
    hp = sliced_quantile_hessian(info.p̃, plan.bins_x, plan.taus_p; β_mean = plan.β_mean_p,
        β_intra = plan.β_intra_p, ε = plan.ε_p)
    (; initial, before, after, stats_before = stats(before), stats_after = stats(after),
        l2_difference = smoothed_density_difference(before.x, before.v, before.w,
            periodic_wrap.(after.x, method.Lx), after.v, after.w; hx = method.hx,
            hv = method.hv, Lx = method.Lx, Vmax = method.Vmax),
        momentum_error = abs(sum(momentum(after)) - sum(p)),
        weight_error = abs(sum(after.w) - sum(before.w)), info,
        min_hessian_q = eigmin(Symmetric(hq)), min_hessian_p = eigmin(Symmetric(hp)))
end

"""
    sliced_stages_example(; Nx = 12, Nv = 10, steps = 500, B_v = 5, B_x = 12)

A small example of the stages of the sliced map: `Nx × Nv` [`electric_wave_markers`](@ref)
pushed to `T = 25` by [`electric_wave_simulation`](@ref) and remapped by a
[`SlicedQuantileRemap`](@ref) with `B_v` velocity and `B_x` position slices (other parameters
at their defaults). Returns `(; before, result)`.
"""
function sliced_stages_example(; Nx::Integer = 12, Nv::Integer = 10, steps::Integer = 500,
        B_v::Integer = 5, B_x::Integer = 12)
    before = electric_wave_simulation(electric_wave_markers(; Nx = Nx, Nv = Nv);
        steps = steps).markers
    (; before, result = remap(SlicedQuantileRemap(; B_v = B_v, B_x = B_x), before))
end

# --- probes of constructions that do not work ----------------------------------------------------

# The markers of the probes: a checkerboard in [−1, 1]² with weights ∝ f(x, v) / N.
function _checkerboard_probe_markers(seed, N, amplitude)
    rng = MersenneTwister(seed)
    q0 = vcat(rand(rng, N ÷ 2), -rand(rng, N ÷ 2))
    v0 = vcat(rand(rng, N ÷ 2), -rand(rng, N ÷ 2))
    w0 = [(1.0 + amplitude * cos(pi * q0[a]) * exp(-v0[a]^2)) / N for a in 1:N]
    q0, v0, w0
end

# One lift in each of the rotated frames (Q, P) = (c q s + s p / s_q, −s q s_q + c p / s_q) with
# the scale s_q = √(P0 / Lx); `step!(Q, P)` returns the lifted (Q, P) or `nothing` to stop.
function _rotated_frame_passes(step, q0, p0, scale_q, passes, angles)
    q, p = copy(q0), copy(p0)
    for pass in 1:passes, (k, θ) in pairs(angles)

        c, s = cos(θ), sin(θ)
        qs, ps = q .* scale_q, p ./ scale_q
        out = step(c .* qs .+ s .* ps, -s .* qs .+ c .* ps, pass, k)
        out === nothing && return q, p, (pass, k)
        Qr2, Pr2 = out
        q = (c .* Qr2 .- s .* Pr2) ./ scale_q
        p = (s .* Qr2 .+ c .* Pr2) .* scale_q
    end
    q, p, nothing
end

"""
    repulsion_wall_probe(; seed = 42, N = 400, passes = 12, α = 0.12, σ = 0.35, k_w = 0.4)

The cotangent lift of `Φ(Q) = Q − α ∇U(Q)` with a Gaussian pair potential of width `σ` and the
quartic wall `k_w Q⁴/4` (in units of the action `u₀ = √(m/N)`), applied in six rotated frames per
pass to the checkerboard markers of the probe (weights `∝ 1 + 0.5 cos(πx) e^{−v²}`). The
Jacobian `DΦ = I − α ∇²U` is symmetric; the wall adds `−3 α k_w Q_a²` to its diagonal. Returns the
first sub-step `(pass, frame)` and the eigenvalue at which `DΦ` stops being positive definite (so
that `Φ` folds), and the sub-step at which `DΦ` stops being finite (`nothing` if never).
"""
function repulsion_wall_probe(; seed::Integer = 42, N::Integer = 400, passes::Integer = 12,
        α::Real = 0.12, σ::Real = 0.35, k_w::Real = 0.4)
    q0, v0, w0 = _checkerboard_probe_markers(seed, N, 0.5)
    m = 1.0
    P0 = m * 1.0 / N
    u0 = sqrt(P0 * 1.0)
    alpha, sigma, kw = α * u0^2, σ * u0, k_w / u0^4
    first_fold = Ref{Any}(nothing)
    function step(Q, P, pass, k)
        n = length(Q)
        dPhi = Matrix{Float64}(I, n, n)
        Qnew = similar(Q)
        inv_s2 = 1.0 / sigma^2
        for a in 1:n
            grad_a = zero(eltype(Q))
            diag_h = zero(eltype(Q))
            qa = Q[a]
            for b in 1:n
                a == b && continue
                r = qa - Q[b]
                kr = exp(-r^2 * inv_s2 / 2)
                kp = -r * inv_s2 * kr
                kpp = (r^2 * inv_s2^2 - inv_s2) * kr
                grad_a += kp / n
                hess_ab = -kpp / n
                dPhi[a, b] = -alpha * hess_ab
                diag_h -= hess_ab
            end
            Qnew[a] = qa - alpha * (grad_a + kw * qa^3)
            dPhi[a, a] = 1.0 - alpha * (diag_h + 3 * kw * qa^2)
        end
        all(isfinite, dPhi) || return nothing
        λ = eigmin(Symmetric(dPhi))
        λ <= 0 && first_fold[] === nothing &&
            (first_fold[] = (pass = pass, frame = k, λ = λ))
        Qnew, transpose(dPhi) \ P
    end
    angles = [0.0, pi / 4, pi / 2, 3pi / 4, pi / 6, 2pi / 3]
    _, _, stop = _rotated_frame_passes(
        step, q0, w0 .* m .* v0, sqrt(P0 / 1.0), passes, angles)
    (first_fold = first_fold[], first_nonfinite = stop)
end

"""
    rotated_frame_lift_probe(; seed = 42, N = 400, passes = 4, β = 0.25, ε = 0.08)

The convex quantile lift of all markers in one slice, `Q̃_a = Q̄ + (1 − β)(Q_a − Q̄) +
(β u₀ / N) Σ_b (Q_a − Q_b)/√((Q_a − Q_b)² + (ε u₀)²)` with its cotangent lift, applied in six
rotated frames per pass to the checkerboard markers of the probe. Returns the quadrant counts of
`(q, p / (m w))` before and after: the lift is symplectic but moves almost no marker into the
empty quadrants.
"""
function rotated_frame_lift_probe(; seed::Integer = 42, N::Integer = 400,
        passes::Integer = 4, β::Real = 0.25, ε::Real = 0.08)
    q0, v0, w0 = _checkerboard_probe_markers(seed, N, 0.5)
    m = 1.0
    P0 = m * 1.0 / N
    u0 = sqrt(P0 * 1.0)
    L, eps = u0, ε * u0
    function step(Q, P, _, _)
        n = length(Q)
        qbar = sum(Q) / n
        Qnew = similar(Q)
        dPhi = fill(β / n, n, n)
        eps2 = eps^2
        coef = β * L / n
        for a in 1:n
            qa = Q[a]
            rep_a = zero(eltype(Q))
            diag_h = zero(eltype(Q))
            for b in 1:n
                a == b && continue
                r = qa - Q[b]
                denom = sqrt(r^2 + eps2)
                rep_a += r / denom
                hab = -coef * eps2 / (denom * denom * denom)
                dPhi[a, b] += hab
                diag_h -= hab
            end
            Qnew[a] = qbar + (1 - β) * (qa - qbar) + coef * rep_a
            dPhi[a, a] += (1 - β) + diag_h
        end
        Qnew, dPhi \ P
    end
    angles = [0.0, pi / 4, pi / 2, 3pi / 4, pi / 6, 2pi / 3]
    q, p, _ = _rotated_frame_passes(step, q0, w0 .* m .* v0, sqrt(P0 / 1.0), passes, angles)
    (quadrants_before = phase_quadrant_counts(q0, v0),
        quadrants_after = phase_quadrant_counts(q, p ./ (m .* w0)))
end

"""
    checkerboard_reweighting_control(; h = 0.25, kwargs...)

The checkerboard of [`checkerboard_experiment`](@ref) after its sliced lifts, reweighted at fixed
momenta: `w̄ ∝ r` with `r` the [`kernel_density_ratio`](@ref) of the input markers at the lifted
points (kernel width `h` in `x` and `v`, periodic on `[−1, 1)`), `Σ w̄ = Σ w`, and
`v̄ = p̂ / (m w̄)`. The input density vanishes in the empty quadrants, so the weights there are
tiny and the velocities large. Returns the smallest and largest `w̄ / w`, the largest `|v̄|` and
the markers.
"""
function checkerboard_reweighting_control(; h::Real = 0.25, kwargs...)
    ex = checkerboard_experiment(; kwargs...)
    before, lifted = ex.before, ex.sliced
    r = kernel_density_ratio(
        lifted.x, lifted.v, before.x, before.v, before.w; hx = h, hv = h,
        Lx = 1.0)
    w̄ = (sum(before.w) / sum(r)) .* r
    p̂ = momentum(lifted)
    v̄ = p̂ ./ (before.m .* w̄)
    (; weight_ratio = extrema(w̄ ./ before.w), max_speed = maximum(abs, v̄),
        markers = Markers(lifted.x, v̄, w̄; m = before.m))
end
