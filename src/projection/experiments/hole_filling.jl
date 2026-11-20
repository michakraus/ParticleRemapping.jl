# The hole-filling experiment of the Hamiltonian particle projection, d = 1: markers on a lattice
# with an elliptic hole, remapped by the kick–drift map, the non-separable map and their controls.

"""
    hole_filling_problem(; hole_x = π, hole_v = 1.2, hole_a = 0.45, hole_b = 0.40, nL = 48,
        L = 2π, V = 4, m = 1)

The markers of the hole-filling experiment: an `nL × nL` lattice of cell centres in
`[0, L) × [−V, V]`, without the markers in the elliptic hole with the semi-axes `hole_a`,
`hole_b` around `(hole_x, hole_v)`, with the weights `w = c f` for
`f = (1 + 0.3 cos x) (2π)^{−1/2} e^{−v²/2}` and `c = 2LV / nL²`.

Returns a `NamedTuple` with the `markers`, their fields `x0`, `v0`, `w0` and momenta `p0`, the
number `N`, the box `L`, `V`, the mass `m`, the lattice spacings `dx`, `dv`, the `hole`
([`EllipticHole`](@ref)), the phase-space `volume`, the constant `c0`, the `density` `f` and
`nL`.
"""
function hole_filling_problem(; hole_x = π, hole_v = 1.2, hole_a = 0.45, hole_b = 0.40,
        nL = 48, L = 2π, V = 4.0, m = 1.0)
    density(x, v) = (1 + 0.3 * cos(x)) * exp(-v^2 / 2) / sqrt(2π)
    dx, dv = L / nL, 2V / nL
    hole = EllipticHole(hole_x, hole_v, hole_a, hole_b; L)
    lattice = [((i - 0.5) * dx, -V + (j - 0.5) * dv) for i in 1:nL for j in 1:nL]
    volume = L * 2V
    c0 = volume / length(lattice)
    kept = filter(z -> !inside_hole(hole, z...), lattice)
    x0 = [z[1] for z in kept]
    v0 = [z[2] for z in kept]
    w0 = c0 .* density.(x0, v0)
    markers = Markers(x0, v0, w0; m)
    (; markers, x0, v0, w0, p0 = momentum(markers),
        N = length(kept), L, V, m, dx, dv, hole,
        volume, c0, density, nL)
end

# The diagnostics of markers (x, v) of the hole problem.
function _hole_diagnostics(prob, hp, x, v)
    (; empty = local_empty_disk(prob.hole, x, v; prob.dx, prob.dv),
        fill = hole_fill(prob.hole, x, v; prob.volume),
        largest_empty = largest_empty_disk(x, v; prob.L, prob.dx, prob.dv),
        discrepancy = density_discrepancy(hp, x, v; prob.L, prob.V))
end

"""
    hole_remap_method(variant; kwargs...)

The method of a variant of the hole-filling experiment: [`KickDriftRemap`](@ref) for
`:weight_dependent`, `:control_A`, `:control_B`, [`NonSeparableRemap`](@ref) for `:non_separable`
and, with `weighted = false`, for `:control_C`. The keyword arguments go to the constructor.
"""
function hole_remap_method(variant::Symbol; kwargs...)
    variant in (:weight_dependent, :control_A, :control_B) &&
        return KickDriftRemap(; variant, kwargs...)
    variant === :non_separable && return NonSeparableRemap(; kwargs...)
    variant === :control_C && return NonSeparableRemap(; weighted = false, kwargs...)
    throw(ArgumentError("unknown variant $variant"))
end

"""
    hole_filling_experiment(; hole_x = π, weights = :exact,
        variants = (:weight_dependent, :control_A, :control_B, :non_separable, :control_C),
        problem = hole_filling_problem(; hole_x), projection = HatProjection(32, 33),
        moser_steps = 40, ncharge = 32, β = 1e-3, maxit = 40, verbose = false)

The hole-filling experiment: the targets of the markers of `problem` ([`projection_targets`](@ref),
with the target weights `c f` for `weights = :exact` or from the lumped projection for
`weights = :projected`), then for each variant ([`hole_remap_method`](@ref)) the fit, the map, the
reweighting with the target weights and the charge-conserving reweighting on `ncharge` periodic
cubic B-splines against the grid charge `ρ0` of the markers before the map. With `verbose = true`
one line per variant is printed while the experiment runs.

Returns a `NamedTuple` with the `problem`, the `projection`, the `targets`, the reference grid
charge `ρ0`, the charge method `charge`, the diagnostics `before` and `target_diagnostics`
(`empty`, `fill`, `largest_empty`, `discrepancy`), and a dictionary `results` with one entry per
variant and one for `:θ0`, the reweighting alone. Each entry holds the lifted positions `X`, the
momenta `P`, the velocities `v1 = P/(m w*)` and `vc = P/(m w′)`, the charge-conserving weights
`wc`, the diagnostics with `w*` (`diagnostics`) and with `w′` (`diagnostics_charge`), the grid
charge residuals `charge_res` of `w*` and `charge_res_c` of `w′`, `charge_ok`, `charge_its`,
`dlnw = max |ln(w′/w*)|`, and for the maps `θ`, `ratio` (objective after the fit relative to
`θ = 0`), `iterations`, `ΔP`, `Δw`, `fit_time`, `map_time`, `stage_iterations`, the
`parametrisation` and the `method` ([`experiment_remap`](@ref)).
"""
function hole_filling_experiment(; hole_x = π, weights::Symbol = :exact,
        variants = (:weight_dependent, :control_A, :control_B, :non_separable, :control_C),
        problem = hole_filling_problem(; hole_x), projection = HatProjection(32, 33),
        moser_steps = 40, ncharge = 32, β = 1e-3, maxit = 40, verbose = false)
    prob = problem
    L, V, m = prob.L, prob.V, prob.m
    density = weights === :exact ? (xs, vs) -> prob.density.(xs, vs) : nothing
    targets = projection_targets(projection, prob.x0, prob.v0, prob.w0; L, V,
        steps = moser_steps, density)
    wt = targets.w
    charge = ChargeConservingReweighting(L, ncharge)
    ρ0 = spline_grid_charge(charge.basis, prob.x0, prob.w0)
    residual(x, w) = maximum(abs, spline_grid_charge(charge.basis, x, w) .- ρ0) /
                     maximum(ρ0)
    function charge_conserving(X, P)
        r = charge_reweight(charge.basis, X, wt, ρ0)
        v1, vc = P ./ (m .* wt), P ./ (m .* r.w)
        (; X, P, v1, wc = r.w, vc, charge_ok = r.converged, charge_its = r.iterations,
            charge_res = residual(X, wt), charge_res_c = residual(X, r.w),
            dlnw = maximum(abs, log.(r.w ./ wt)),
            diagnostics = _hole_diagnostics(prob, projection, X, v1),
            diagnostics_charge = _hole_diagnostics(prob, projection, X, vc))
    end
    results = Dict{Symbol, Any}(:θ0 => charge_conserving(prob.x0, prob.p0))
    for variant in variants
        method = hole_remap_method(variant; L, V, projection, moser_steps, β, maxit)
        res = remap(method, prob.markers; targets)
        i = res.info
        r = (;
            charge_conserving(i.X, i.P)..., i.θ, ratio = i.objective_ratio, i.iterations,
            ΔP = i.momentum_change, Δw = i.weight_change, i.fit_time, i.map_time,
            i.stage_iterations, i.parametrisation, method)
        verbose &&
            @printf("  %-17s empty disk %5.2f, fill %4.2f, objective %.3f, %d its, %.0f s\n",
                variant, r.diagnostics.empty, r.diagnostics.fill, r.ratio, r.iterations,
                r.fit_time)
        results[variant] = r
    end
    (; problem = prob, projection, targets, ρ0, charge,
        before = _hole_diagnostics(prob, projection, prob.x0, prob.v0),
        target_diagnostics = _hole_diagnostics(prob, projection, targets.x, targets.v), results)
end

"""
    projected_weight_errors(problem = hole_filling_problem(); projection = HatProjection(32, 33),
        moser_steps = 40, bands = ((0, 1), (1, 2), (2, 3), (3, 4)))

The relative error of the target weights from the lumped projection against `c f` at the target
points, by velocity band ([`weight_errors_by_band`](@ref)).
"""
function projected_weight_errors(problem = hole_filling_problem();
        projection = HatProjection(32, 33), moser_steps = 40,
        bands = ((0, 1), (1, 2), (2, 3), (3, 4)))
    prob = problem
    tg = projection_targets(projection, prob.x0, prob.v0, prob.w0; prob.L, prob.V,
        steps = moser_steps)
    we = projection_target_weights(
        projection, prob.x0, prob.v0, prob.w0, tg.x, tg.v; prob.L,
        prob.V, density = (xs, vs) -> prob.density.(xs, vs))
    weight_errors_by_band(tg.v, tg.w, we; bands)
end
