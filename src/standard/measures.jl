# The measures of the standard test set: one marker set, or a marker set before and after a
# remap, measured on the box of a case.

# The cells per axis of the evenness and the nodes of the grid of the charge.
_case_cells(case::StandardCase) = case.experiment === :checkerboard ? 8 : 10
_case_nodes(case::StandardCase) = case.experiment === :checkerboard ? 8 : 16

"""
    case_evenness(case::StandardCase, mk::Markers)
        -> (; empty, rstd, support_cells, support_empty, support_rstd, window)

The evenness of the markers `mk` on the box of `case`: `empty` and `rstd` of
[`phase_cell_statistics`](@ref) on the `g × g` cells of the box (`g = 8` for the experiment
`:checkerboard`, `10` otherwise; every marker in the box counted), and the same two over the
`support_cells` cells whose centre has a positive represented density (`support_rstd` with the
mean over those cells). For family `:W`, `window` holds `(; empty, rstd)` of
`phase_cell_statistics(mk.x, mk.v)`, the old window `[−π, π) × [−2.25, 2.25]`; otherwise
`nothing`.
"""
function case_evenness(case::StandardCase, mk::Markers)
    g = _case_cells(case)
    hL, V = case.L / 2, case.V
    st = phase_cell_statistics(mk.x, mk.v; gx = g, gv = g, Lx = hL, Vmax = V)
    support = [case.represented(-hL + (i - 0.5) * 2hL / g, -V + (j - 0.5) * 2V / g) > 0
               for i in 1:g, j in 1:g]
    c = st.counts[support]
    μ = sum(c) / length(c)
    window = case.family === :W ?
             (w = phase_cell_statistics(mk.x, mk.v); (; w.empty, w.rstd)) : nothing
    (; st.empty, st.rstd, support_cells = count(support), support_empty = count(==(0), c),
        support_rstd = sqrt(sum((c .- μ) .^ 2) / length(c)) / μ, window)
end

"""
    case_outside(case::StandardCase, mk::Markers) -> (; count, weight_share, max_v, nonfinite)

The markers of `mk` outside the velocity box of `case`: their `count` (the markers with
`!(abs(v) <= V)`, a velocity `NaN` included), the `weight_share` of the finite markers among
them in the total of the finite weights, the largest finite `|v|` of `mk` (`0.0` when there is
none), and the number `nonfinite` of markers with a non-finite position, velocity or weight.
"""
function case_outside(case::StandardCase, mk::Markers)
    outside = .!(abs.(mk.v) .<= case.V)
    finite = isfinite.(mk.x) .& isfinite.(mk.v) .& isfinite.(mk.w)
    finite_w = isfinite.(mk.w)
    (; count = count(outside),
        weight_share = sum(mk.w[outside .& finite]; init = 0.0) / sum(mk.w[finite_w]),
        max_v = maximum(abs, filter(isfinite, mk.v); init = 0.0),
        nonfinite = count(!, finite))
end

# The midpoints of n equal cells of [a, b].
_midpoints(a, b, n) = [a + (i - 0.5) * (b - a) / n for i in 1:n]

# The kernel K * f of the representation error on the ngrid × ngrid cell centres of the box,
# by the midpoint rule on nq × nq cells.
function _smoothed_density(case::StandardCase, variant::Symbol, ngrid::Integer, nq::Integer)
    f = case_density(case, variant)
    L, V = case.L, case.V
    s = L / 2π
    hx, hv = 0.45 * s, 0.35 * V / 2.25
    xg, vg = _midpoints(-L / 2, L / 2, ngrid), _midpoints(-V, V, ngrid)
    xq, vq = _midpoints(-L / 2, L / 2, nq), _midpoints(-V, V, nq)
    F = [f(a, b) for a in xq, b in vq] .* (L / nq) .* (2V / nq)
    Kx = [exp(-SlicedRedistribution._periodic_distance2(a - b, s) / hx^2 / 2)
          for a in xg, b in xq]
    Kv = [exp(-(a - b)^2 / hv^2 / 2) for a in vg, b in vq]
    Kx * F * transpose(Kv)
end

# The exact grid charge E_i = ∫ S_i(x) ρ(x) dx, ρ = ∫_{−V}^{V} f dv, by the midpoint rule on nq
# points in each variable.
function _exact_charge(case::StandardCase, variant::Symbol, nq::Integer)
    f = case_density(case, variant)
    L, V = case.L, case.V
    g = PeriodicGrid(L, _case_nodes(case))
    vq = _midpoints(-V, V, nq)
    E = zeros(g.Ng)
    for xq in _midpoints(-L / 2, L / 2, nq)
        ρ = sum(f(xq, b) for b in vq) * 2V / nq
        i, j, Si, Sj = shape(g, xq)
        E[i] += Si * ρ * L / nq
        E[j] += Sj * ρ * L / nq
    end
    E
end

function _check_reference(reference, variant::Symbol)
    reference.variant === variant || throw(ArgumentError(
        "the reference is of the variant :$(reference.variant), the measure asks for :$variant"))
end

"""
    case_reference(case::StandardCase, variant; ngrid = 30, nq = 240, nq_charge = 960)
        -> (; variant, smoothed, charge)

The integrals of the density of `case` of `variant` that the measures compare the markers with,
computed once: the `ngrid × ngrid` matrix `smoothed` of `K * f` of
[`case_representation`](@ref) (with `nq`) and the vector `charge` of the exact grid charge of
[`case_charge_noise`](@ref) (with `nq_charge`). Pass it to both measures, and to
[`run_case`](@ref), as the keyword `reference`.
"""
function case_reference(case::StandardCase, variant::Symbol; ngrid::Integer = 30,
        nq::Integer = 240, nq_charge::Integer = 960)
    _check_variant(variant)
    (; variant, smoothed = _smoothed_density(case, variant, ngrid, nq),
        charge = _exact_charge(case, variant, nq_charge))
end

"""
    case_representation(case::StandardCase, mk::Markers, variant; ngrid = 30, nq = 240,
                        reference = nothing)

The representation error `‖Σ_a w_a K(· − x_a, · − v_a) − K * f‖₂ / ‖K * f‖₂` of the markers
`mk` on the `ngrid × ngrid` cell centres of the box of `case`, with `f` the density of
`variant` ([`case_density`](@ref)), the kernel of [`smoothed_density_difference`](@ref)
(periodic in `x` with the period `L`) at the widths `hx = 0.45 L/2π` and `hv = 0.35 V/2.25`,
and `K * f` the midpoint rule on the `nq × nq` cells of the box. With a
[`case_reference`](@ref) as `reference`, `K * f` is its `smoothed`, which must be of the same
`variant` (an `ArgumentError`) and `ngrid × ngrid` (a `DimensionMismatch`).
"""
function case_representation(case::StandardCase, mk::Markers, variant::Symbol;
        ngrid::Integer = 30, nq::Integer = 240, reference = nothing)
    _check_variant(variant)
    Kf = if reference === nothing
        _smoothed_density(case, variant, ngrid, nq)
    else
        _check_reference(reference, variant)
        size(reference.smoothed) == (ngrid, ngrid) || throw(DimensionMismatch(
            "the reference holds a $(size(reference.smoothed)) grid, the measure asks for " *
            "$ngrid × $ngrid"))
        reference.smoothed
    end
    L, V = case.L, case.V
    s = L / 2π
    inv_hx2, inv_hv2 = 1 / (0.45 * s)^2, 1 / (0.35 * V / 2.25)^2
    xg, vg = _midpoints(-L / 2, L / 2, ngrid), _midpoints(-V, V, ngrid)
    fm = [sum(mk.w[c] *
              exp(-(SlicedRedistribution._periodic_distance2(a - mk.x[c], s) * inv_hx2 +
                    (b - mk.v[c])^2 * inv_hv2) / 2) for c in eachindex(mk.w))
          for a in xg, b in vg]
    norm(fm - Kf) / norm(Kf)
end

"""
    case_charge_noise(case::StandardCase, mk::Markers, variant; nq = 960, reference = nothing)

The charge noise `‖M − E‖₂ / ‖E‖₂` of the markers `mk`: `M` the grid charge
`grid_moment(PeriodicGrid(L, Ng), mk.x, mk.v, mk.w, 0)` (`Ng = 8` for the experiment
`:checkerboard`, `16` otherwise) and `E_i = ∫ S_i(x) ρ(x) dx` the exact grid charge of the
density `f` of `variant`, `ρ = ∫_{−V}^{V} f dv`, both integrals by the midpoint rule on `nq`
points. With a [`case_reference`](@ref) as `reference`, `E` is its `charge`, which must be of
the same `variant` (an `ArgumentError`) and of length `Ng` (a `DimensionMismatch`).
"""
function case_charge_noise(case::StandardCase, mk::Markers, variant::Symbol;
        nq::Integer = 960, reference = nothing)
    _check_variant(variant)
    Ng = _case_nodes(case)
    E = if reference === nothing
        _exact_charge(case, variant, nq)
    else
        _check_reference(reference, variant)
        length(reference.charge) == Ng || throw(DimensionMismatch(
            "the reference holds $(length(reference.charge)) grid charges, the case has " *
            "$Ng nodes"))
        reference.charge
    end
    M = grid_moment(PeriodicGrid(case.L, Ng), mk.x, mk.v, mk.w, 0)
    norm(M - E) / norm(E)
end

"""
    case_weights(mk::Markers) -> (; min, max, ratio, small_share)

The spread of the weights of `mk`: the smallest and the largest weight, their `ratio`
(`max/min`, `Inf` when `min ≤ 0`), and the `small_share` of the total weight in the markers
whose weight lies below `1e-2` of the mean weight.
"""
function case_weights(mk::Markers)
    w = mk.w
    lo, hi = extrema(w)
    small = w .< 1e-2 * sum(w) / length(w)
    (; min = lo, max = hi, ratio = lo > 0 ? hi / lo : Inf,
        small_share = sum(w[small]; init = 0.0) / sum(w))
end

"""
    case_conservation(case::StandardCase, before::Markers, after::Markers)
        -> (; weight, momentum, kinetic_energy, charge, current)

The relative changes of a remap from `before` to `after`, each with the scale of `before`:
`|ΔΣ w| / Σ |w|`, `|ΔP| / Σ m |w| |v|`, `|ΔKE| / KE`, and the largest change of the grid
moments `M^(k)`, `k = 0, 1` (the charge and the current), on `PeriodicGrid(L, Ng)` divided by
the largest moment `M^(k)(|v|)` of `before`, so that it is finite also where a node holds no
marker before.
"""
function case_conservation(case::StandardCase, before::Markers, after::Markers)
    g = PeriodicGrid(case.L, _case_nodes(case))
    moment(k) = maximum(abs, grid_moment(g, after, k) .- grid_moment(g, before, k)) /
                maximum(grid_moment(g, before.x, abs.(before.v), before.w, k))
    (; weight = abs(total_weight(after) - total_weight(before)) / sum(abs, before.w),
        momentum = abs(total_momentum(after) - total_momentum(before)) /
                   (before.m * sum(abs.(before.w) .* abs.(before.v))),
        kinetic_energy = abs(kinetic_energy(after) - kinetic_energy(before)) /
                         kinetic_energy(before),
        charge = moment(0), current = moment(1))
end

"""
    case_displacement(case::StandardCase, before::Markers, after::Markers)
        -> (; max_dx, max_dv, max_weight_change, kinetic_energy_change)

The changes of a remap marker by marker: the largest periodic `|Δx|` (wrapped into
`[−L/2, L/2)`), the largest `|Δv|`, the largest `|Δw| / |w|`, and the signed relative change
`(KE_after − KE_before) / KE_before` of the kinetic energy. A `DimensionMismatch` refuses two
marker sets of different sizes.
"""
function case_displacement(case::StandardCase, before::Markers, after::Markers)
    length(before) == length(after) || throw(DimensionMismatch(
        "$(length(before)) markers before the remap but $(length(after)) after"))
    (; max_dx = maximum(abs, periodic_wrap.(after.x .- before.x, case.L / 2)),
        max_dv = maximum(abs, after.v .- before.v),
        max_weight_change = maximum(abs.(after.w .- before.w) ./ abs.(before.w)),
        kinetic_energy_change = (kinetic_energy(after) - kinetic_energy(before)) /
                                kinetic_energy(before))
end

export case_evenness, case_outside, case_representation, case_charge_noise, case_weights,
       case_conservation, case_displacement, case_reference
