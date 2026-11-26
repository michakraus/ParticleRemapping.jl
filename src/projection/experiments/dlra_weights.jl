# The experiment of the target weights after a push, d = 1: the markers of the hole-filling
# experiment follow the characteristics of the Vlasov equation with a fixed field; the target weights
# come from the exact density, from dynamical low-rank solutions and from the lumped projection.

_tolerance_name(ϑ) = "ϑ = 1e" * string(round(Int, log10(ϑ)))

"""
    dlra_weights_experiment(; hole_x = π, field_amplitude = 0.15, T = 3.0, nx = 64, nv = 256,
        Vd = 6.0, dt = 0.01, tolerances = (1e-3, 1e-4, 1e-5, 1e-6), fit = true,
        problem = hole_filling_problem(; hole_x), projection = HatProjection(32, 33),
        moser_steps = 40, ncharge = 32, β = 1e-3, maxit = 40, rng = Random.default_rng(),
        verbose = false)

The experiment of the target weights after a push: the markers of `problem` follow the
characteristics of `∂_t f + v ∂_x f + E(x) ∂_v f = 0` with `E = −field_amplitude sin x` up to `T`
([`vlasov_characteristics`](@ref)); `f(T) = f_0 ∘ Φ_{−T}` comes from backward characteristics. The
BUG integrator ([`dlra_solve`](@ref)) solves the same equation on `nx × nv` points of
`[0, L) × [−Vd, Vd]` with the step `dt`, from rank four at each tolerance of `tolerances` and at
the fixed rank one; `rng` draws the padding columns of the start. The target weights come from
`f(T)`, from each low-rank solution ([`lowrank_evaluate`](@ref)) and from the lumped projection.
For each source without nonpositive target weights, the non-separable map is fitted, and the
charge-conserving reweighting follows against the grid charge `ρT` of the pushed markers; `fit`
is `true` (every source), `false` (none) or a collection of source names.

Returns a `NamedTuple` with the `problem`, the `grid`, the field `E`, the time `T`, the step `dt`,
the pushed markers `xT`, `vT`, `pT`, `ρT`, the diagnostics of the pushed markers `pushed`
(`largest_empty`, `fill` of the image of the hole, `discrepancy`), the `targets` (`x`, `v` and the
exact weights `w`), the low-rank solutions `solutions` (a dictionary from the source name to
`(; f, ranks, mass_error, time)`), the grid values `Fe` of `f(T)`, the image of the hole outline
`hole_image`, and the dictionary `results` from the source name (`"f(T)"`, `"ϑ = 1e-3"`, …,
`"rank 1"`, `"projection"`) to the target weights `wt`, the band errors `errs`, the number `nneg`
of nonpositive weights and, after a fit, `X`, `P`, `v1`, `vc`, `wc`, `θ`, `ratio`, `iterations`,
`status` (the `SolverStatus` of the fit), `ΔP`, `fit_time`, the diagnostics with `w*`
(`diagnostics`) and with `w′` (`diagnostics_charge`), the residuals `charge_res`, `charge_res_c`,
`charge_ok` and `dlnw`.
"""
function dlra_weights_experiment(; hole_x = π, field_amplitude = 0.15, T = 3.0, nx = 64,
        nv = 256, Vd = 6.0, dt = 0.01, tolerances = (1e-3, 1e-4, 1e-5, 1e-6), fit = true,
        problem = hole_filling_problem(; hole_x), projection = HatProjection(32, 33),
        moser_steps = 40, ncharge = 32, β = 1e-3, maxit = 40,
        rng::AbstractRNG = Random.default_rng(), verbose = false)
    prob = problem
    L, V, m, hole = prob.L, prob.V, prob.m, prob.hole
    E(x) = -field_amplitude * sin(x)
    grid = VlasovGrid(nx, nv, L, Vd, E)
    xT, vT = vlasov_characteristics(E, prob.x0, prob.v0, T)
    maximum(abs, vT) < V || error("a pushed marker leaves |v| < V = $V")
    pT = m .* prob.w0 .* vT
    charge = ChargeConservingReweighting(L, ncharge)
    ρT = spline_grid_charge(charge.basis, xT, prob.w0)
    f_pushed(xs, vs) = prob.density.(vlasov_characteristics(E, xs, vs, -T)...)
    function pushed_fill(x, v)
        xb, vb = vlasov_characteristics(E, x, v, -T)
        count(k -> inside_hole(hole, xb[k], vb[k]), eachindex(xb)) /
        (π * hole.a * hole.b * length(x) / prob.volume)
    end
    diagnostics(x, v) = (; largest_empty = largest_empty_disk(x, v; L, prob.dx, prob.dv),
        fill = pushed_fill(x, v), discrepancy = density_discrepancy(projection, x, v; L, V))

    F0 = prob.density.(grid.x, grid.v')
    mass0 = grid.hx * grid.hv * sum(F0)
    solutions = Dict{String, Any}()
    sources = Pair{String, Any}["f(T)" => f_pushed]
    cases = [[(_tolerance_name(ϑ), (; tol = ϑ, rmin = 4)) for ϑ in tolerances];
             ("rank 1", (; tol = 0.0, rmin = 1, rmax = 1))]
    for (name, kw) in cases
        t = @elapsed f, ranks = dlra_solve(grid, F0, T, dt; kw..., rng)
        mass_error = abs(lowrank_mass(grid, f) - mass0) / mass0
        solutions[name] = (; f, ranks, mass_error, time = t)
        verbose &&
            @printf("  BUG %-9s  largest rank %2d, final rank %2d, mass error %.1e\n",
                name, maximum(ranks), lowrank_rank(f), mass_error)
        push!(sources, name => (xs, vs) -> lowrank_evaluate(grid, f, xs, vs))
    end
    push!(sources, "projection" => nothing)

    tg = projection_targets(projection, xT, vT, prob.w0; L, V, steps = moser_steps,
        density = f_pushed)
    results = Dict{String, Any}()
    for (name, src) in sources
        wt = name == "f(T)" ? tg.w :
             projection_target_weights(projection, xT, vT, prob.w0, tg.x, tg.v; L, V,
            density = src)
        errs = weight_errors_by_band(tg.v, wt, tg.w)
        nneg = count(<=(0), wt)
        r = (; wt, errs, nneg)
        if nneg == 0 && (fit === true || (fit !== false && name in fit))
            method = NonSeparableRemap(; L, V, projection, moser_steps, β, maxit, charge)
            res = remap(
                method, Markers(xT, vT, prob.w0; m); targets = (; tg.x, tg.v, w = wt),
                reference_charge = ρT)
            i = res.info
            v1, vc = i.P ./ (m .* wt), i.P ./ (m .* i.wc)
            r = (; r..., i.X, i.P, v1, vc, i.wc, i.θ, ratio = i.objective_ratio,
                i.iterations, i.status, ΔP = abs(sum(i.P) - sum(pT)) / sum(abs, pT), i.fit_time,
                charge_ok = i.charge_converged, charge_res = i.charge_residual_target,
                charge_res_c = i.charge_residual, dlnw = i.max_log_change,
                diagnostics = diagnostics(i.X, v1), diagnostics_charge = diagnostics(i.X, vc))
            verbose &&
                @printf("  %-11s empty disk %5.2f, fill %4.2f, density %.3f, %.0f s\n",
                    name, r.diagnostics.largest_empty, r.diagnostics.fill,
                    r.diagnostics.discrepancy, r.fit_time)
        end
        results[name] = r
    end
    xs, vs = vlasov_characteristics(E, repeat(grid.x, length(grid.v)),
        repeat(grid.v; inner = length(grid.x)), -T)
    Fe = reshape(prob.density.(xs, vs), length(grid.x), length(grid.v))
    hx, hv = hole_outline(hole)
    (; problem = prob, grid, E, T, dt, xT, vT, pT, ρT, pushed = diagnostics(xT, vT),
        targets = tg, solutions, Fe, hole_image = vlasov_characteristics(E, hx, hv, T), results)
end
