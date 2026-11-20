# The composed sliced redistribution R = F ∘ T_μ ∘ Ψ_p ∘ Ψ_q in 1D-1V.

"""
    SlicedLiftPlan

The data of the two lifts `Ψ_q` and `Ψ_p` of a sliced redistribution, read from one input state
and frozen: the velocity slices `bins_v` with their position targets `taus_x` and momentum
anchors `p_anc`, the position slices `bins_x` with their momentum targets `taus_p` and position
anchors `q_anc`, the slot velocity `slot_v` of each marker (the velocity target of its rank
within its position slice), and the parameters `β_mean`, `β_intra`, `ε` of each lift. For a
frozen plan, [`apply_sliced_lift`](@ref) is a smooth symplectic map of `(q, p)`.
Build it with [`sliced_lift_plan`](@ref).
"""
struct SlicedLiftPlan{T <: Real}
    bins_v::Vector{Vector{Int}}
    taus_x::Vector{Vector{T}}
    p_anc::Vector{T}
    bins_x::Vector{Vector{Int}}
    taus_p::Vector{Vector{T}}
    q_anc::Vector{T}
    slot_v::Vector{T}
    β_mean_q::T
    β_intra_q::T
    ε_q::T
    β_mean_p::T
    β_intra_p::T
    ε_p::T
end

"""
    sliced_lift_plan(x, v, p; m, w_ref, Lx, Vmax, B_v, B_x, β_mean_q = 0.60, β_intra_q = 0.72,
                     ε_q = 0.4, β_mean_p = 0.60, β_intra_p = 0.72, ε_p) -> SlicedLiftPlan

The plan of the sliced lifts for markers with positions `x`, velocities `v` and momenta `p`:

1. `B_v` velocity slices of `n_v = N / B_v` markers, each ordered by `x`; the position targets are
   the [`uniform_quantile_targets`](@ref) of `[−Lx, Lx]` with `n_v` points; the momentum anchor
   of a marker is `m w_ref` times the centre of its velocity slice in `[−Vmax, Vmax]`;
2. the lift `Ψ_q` with these data, run once to find `q̃`, `p̃`;
3. `B_x` position slices of `n_x = N / B_x` markers by `q̃`, each ordered by `p̃`; the position
   anchor is the centre of the slice in `[−Lx, Lx]`; the momentum targets are `m w_ref` times the
   `n_x` uniform quantile targets of `[−Vmax, Vmax]`.

The slices and anchors are piecewise constant functions of the state: the composed map is smooth
and symplectic where they do not change, and they change where two markers swap their order.
"""
function sliced_lift_plan(x::AbstractVector, v::AbstractVector, p::AbstractVector; m::Real,
        w_ref::Real, Lx::Real, Vmax::Real, B_v::Integer, B_x::Integer, β_mean_q::Real = 0.60,
        β_intra_q::Real = 0.72, ε_q::Real = 0.4, β_mean_p::Real = 0.60,
        β_intra_p::Real = 0.72, ε_p::Real)
    N = length(x)
    n_v = N ÷ B_v
    n_x = N ÷ B_x
    bins_v = sorted_slices(v, x, B_v)
    taus_x = [uniform_quantile_targets(Lx, n_v) for _ in 1:B_v]
    p_anc = slice_anchors(bins_v, [m * w_ref * vb
                                   for vb in uniform_quantile_targets(Vmax, B_v)])
    q̃, p̃ = split_quantile_ot_lift(x, p, bins_v, taus_x; anchors = p_anc, β_mean = β_mean_q,
        β_intra = β_intra_q, ε = ε_q)
    bins_x = sorted_slices(q̃, p̃, B_x)
    q_anc = slice_anchors(bins_x, uniform_quantile_targets(Lx, B_x))
    v_grid = uniform_quantile_targets(Vmax, n_x)
    taus_p = [[m * w_ref * vj for vj in v_grid] for _ in 1:B_x]
    slot_v = zeros(N)
    for idx in bins_x, (r, a) in enumerate(idx)

        slot_v[a] = v_grid[r]
    end
    SlicedLiftPlan(bins_v, taus_x, p_anc, bins_x, taus_p, q_anc, slot_v,
        promote(β_mean_q, β_intra_q, ε_q, β_mean_p, β_intra_p, ε_p)...)
end

"""
    apply_sliced_lift(plan::SlicedLiftPlan, q, p; with_cotangent = true) -> (; q̃, p̃, q̂, p̂)

The two lifts of a frozen plan: `(q̃, p̃) = Ψ_q(q, p)` equalises the positions within each
velocity slice ([`split_quantile_ot_lift`](@ref)), `(q̂, p̂) = Ψ_p(q̃, p̃)` the momenta within each
position slice ([`dual_quantile_ot_lift`](@ref)). Both keep `Σ dq ∧ dp` and `Σ p`.
`with_cotangent = false` drops the cotangent parts (a control).
"""
function apply_sliced_lift(plan::SlicedLiftPlan, q::AbstractVector, p::AbstractVector;
        with_cotangent::Bool = true)
    q̃, p̃ = split_quantile_ot_lift(q, p, plan.bins_v, plan.taus_x; anchors = plan.p_anc,
        β_mean = plan.β_mean_q, β_intra = plan.β_intra_q, ε = plan.ε_q,
        with_cotangent = with_cotangent)
    q̂, p̂ = dual_quantile_ot_lift(q̃, p̃, plan.bins_x, plan.taus_p; anchors = plan.q_anc,
        β_mean = plan.β_mean_p, β_intra = plan.β_intra_p, ε = plan.ε_p,
        with_cotangent = with_cotangent)
    (; q̃, p̃, q̂, p̂)
end

"""
    SlicedQuantileRemap(; Lx = π, Vmax = 2.25, B_v = 15, B_x = 90, β_mean_q = 0.75,
                        β_intra_q = 0.78, ε_q_spacings = 2.0, β_mean_p = 0.78, β_intra_p = 0.80,
                        ε_p_spacings = 2.0, hx = 0.45, hv = 0.35, ratio_floor = 0.08,
                        iterations = 3)

The sliced symplectic redistribution `R = F ∘ T_μ ∘ Ψ_p ∘ Ψ_q` on the box
`[−Lx, Lx) × [−Vmax, Vmax]` (periodic in `x`), with the parameters of the long-time experiment:

- `Ψ_q`, `Ψ_p`: the lifts of [`sliced_lift_plan`](@ref) with `B_v` velocity and `B_x` position
  slices; the mollifier widths are `ε_q_spacings` target spacings `2Lx / n_v` and `ε_p_spacings`
  spacings `2Lp / n_x`, with `Lp` the largest momentum target;
- `T_μ`: the [`canonical_weight_squeeze`](@ref) with the [`squeeze_multipliers`](@ref) (kernel
  widths `hx`, `hv`, floor `ratio_floor`, `iterations` predictor–corrector passes);
- `F`: the weights `w̄ = w_ref μ` and velocities `v̄ = p̄ / (m w̄)`.

`N` must be a multiple of `B_v` and `B_x`. For 1D-1V markers the map keeps `Σ w`, `Σ p` and, for
the frozen plan and multipliers, `Σ dx ∧ dp`. For `d × N` markers stages 1–3 act on each axis
with one common `μ` (see [`remap`](@ref)).
"""
struct SlicedQuantileRemap{T <: Real} <: AbstractRemap
    Lx::T
    Vmax::T
    B_v::Int
    B_x::Int
    β_mean_q::T
    β_intra_q::T
    ε_q_spacings::T
    β_mean_p::T
    β_intra_p::T
    ε_p_spacings::T
    hx::T
    hv::T
    ratio_floor::T
    iterations::Int
end

function SlicedQuantileRemap(; Lx::Real = π, Vmax::Real = 2.25, B_v::Integer = 15,
        B_x::Integer = 90, β_mean_q::Real = 0.75, β_intra_q::Real = 0.78,
        ε_q_spacings::Real = 2.0, β_mean_p::Real = 0.78, β_intra_p::Real = 0.80,
        ε_p_spacings::Real = 2.0, hx::Real = 0.45, hv::Real = 0.35, ratio_floor::Real = 0.08,
        iterations::Integer = 3)
    Lx, Vmax, β_mean_q, β_intra_q, ε_q_spacings, β_mean_p,
    β_intra_p, ε_p_spacings, hx, hv,
    ratio_floor = promote(float(Lx), Vmax, β_mean_q, β_intra_q, ε_q_spacings, β_mean_p,
        β_intra_p, ε_p_spacings, hx, hv, ratio_floor)
    SlicedQuantileRemap(Lx, Vmax, Int(B_v), Int(B_x), β_mean_q, β_intra_q, ε_q_spacings,
        β_mean_p, β_intra_p, ε_p_spacings, hx, hv, ratio_floor, Int(iterations))
end

# The plan of one axis: positions x, velocities v, momenta p of that axis.
function _plan(method::SlicedQuantileRemap, x, v, p, m, w_ref)
    N = length(x)
    n_v = N ÷ method.B_v
    n_x = N ÷ method.B_x
    Lp = maximum(abs, [m * w_ref * vj for vj in uniform_quantile_targets(method.Vmax, n_x)])
    sliced_lift_plan(x, v, p; m = m, w_ref = w_ref, Lx = method.Lx, Vmax = method.Vmax,
        B_v = method.B_v, B_x = method.B_x, β_mean_q = method.β_mean_q,
        β_intra_q = method.β_intra_q, ε_q = method.ε_q_spacings * (2method.Lx / n_v),
        β_mean_p = method.β_mean_p, β_intra_p = method.β_intra_p,
        ε_p = method.ε_p_spacings * (2Lp / n_x))
end

"""
    remap(method::SlicedQuantileRemap, mk::Markers) -> RemapResult

Redistribute the markers by the sliced map `R = F ∘ T_μ ∘ Ψ_p ∘ Ψ_q` (see
[`SlicedQuantileRemap`](@ref)). The `info` holds the `plan` (a [`SlicedLiftPlan`](@ref), one per
axis for `d × N` markers), the intermediate states `q̃`, `p̃`, `q̂`, `p̂`, the squeezed state `q̄`,
`p̄`, the multipliers `μ` and the reference weight `w_ref = Σ w / N`. The positions are not wrapped
into `[−Lx, Lx)`.
"""
function remap(method::SlicedQuantileRemap, mk::Markers{T, <:AbstractVector}) where {T}
    m = mk.m
    w_ref = sum(mk.w) / length(mk)
    p = momentum(mk)
    plan = _plan(method, mk.x, mk.v, p, m, w_ref)
    stages = apply_sliced_lift(plan, mk.x, p)
    μ = squeeze_multipliers(stages.q̂, stages.p̂, plan.q_anc, mk.x, mk.v, mk.w; m = m,
        w_ref = w_ref, hx = method.hx, hv = method.hv, ratio_floor = method.ratio_floor,
        iterations = method.iterations, Lx = method.Lx)
    q̄, p̄ = canonical_weight_squeeze(stages.q̂, stages.p̂, μ, plan.q_anc)
    w̄ = w_ref .* μ
    v̄ = p̄ ./ (m .* w̄)
    RemapResult(Markers(q̄, v̄, w̄; m = m), merge(stages, (; q̄, p̄, μ, w_ref, plan)))
end
