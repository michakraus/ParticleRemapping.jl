# The driver of the standard test set: the constructions on the box of a case, and one row of
# every measure for a remap of a case.

# The divisor of N nearest to t, ties to the smaller.
_nearest_divisor(N, t) = argmin(d -> (abs(d - t), d), [d for d in 1:N if N % d == 0])

"""
    case_method(key, case::StandardCase; variant = :positive)

The method of the construction `key`, a key of [`STANDARD_METHODS`](@ref) or of
[`REFERENCE_METHODS`](@ref), for `case`: the parameters of
[`standard_method`](@ref)`(key, case.experiment; N, density)` with the box of the case in place
of the box of the experiment (`L = case.L` everywhere, `V = case.V` for the projection and the
non-separable stages, `Vmax = case.V` for the spline-projection redistribution and the sliced
map, `Lx = case.L/2` for the sliced map, `vmax = case.V` for the target of the transport) and
the density [`case_density`](@ref)`(case, variant)`, shifted by `L/2` as `standard_method`
shifts it. The projection, the non-separable stages and the spline-projection redistribution
act inside a [`ShiftedRemap`](@ref); the spline-projection redistribution scales its
relaxation step to `0.006 · 16 / N` at every `N` other than the full size (400 markers for the
experiment `:checkerboard`, 900 for `:wave`). The sliced map takes for `B_v` and `B_x` the
divisors of `N` nearest to `B_v0 √(N/N0)` and `B_x0 √(N/N0)`, ties to the smaller, with
`(B_v0, B_x0, N0) = (10, 20, 400)` for `:checkerboard` and `(15, 90, 900)` for `:wave`.
`:grid_lumped` gives the [`GridResample`](@ref) of the box with the hats of `standard_method`,
`:grid_exact` the same with the density of `variant`. An `ArgumentError` refuses any other key.
"""
function case_method(key::Symbol, case::StandardCase; variant::Symbol = :positive)
    cb = case.experiment === :checkerboard
    L, V = case.L, case.V
    shift = L / 2
    N = length(case.markers)
    hats = cb ? HatProjection(12, 13) : HatProjection(16, 17)
    f = case_density(case, variant)
    shifted = (xs, vs) -> f.(periodic_wrap.(xs .- shift, L / 2), vs)
    if key === :projection
        ShiftedRemap(
            KickDriftRemap(; L, V, projection = hats, maxit = 20, density = shifted),
            shift, L)
    elseif key === :nonseparable
        ShiftedRemap(
            NonSeparableRemap(; L, V, nξ = 6, nv = 7, nℓ = 3, projection = hats, maxit = 5,
                density = shifted),
            shift, L)
    elseif key === :spline
        Nr = _reduced(N, cb ? 400 : 900)
        relax = Nr === nothing ? (;) : (; relax_step = 0.006 * 16 / Nr)
        ShiftedRemap(
            SplineRedistribution(; Lx = L, Vmax = V, σx = L / 4, σv = cb ? 0.5 : 1.2,
                relax...),
            shift, L)
    elseif key === :gauge
        GaugeRezoning(PeriodicGrid(L, _case_nodes(case)))
    elseif key === :transport
        grid = PeriodicGrid(L, _case_nodes(case))
        CurrentQuantileMatching(grid, (x, p) -> (copy(x), copy(p)),
            DensityTarget(grid, (x, v) -> f(periodic_wrap(x, L / 2), v); vmax = V))
    elseif key === :sliced
        B_v0, B_x0, N0 = cb ? (10, 20, 400) : (15, 90, 900)
        B_v = _nearest_divisor(N, B_v0 * sqrt(N / N0))
        B_x = _nearest_divisor(N, B_x0 * sqrt(N / N0))
        cb ? SlicedQuantileRemap(; Lx = L / 2, Vmax = V, B_v, B_x, hx = 0.15, hv = 0.15) :
        SlicedQuantileRemap(; Lx = L / 2, Vmax = V, B_v, B_x)
    elseif key === :grid_lumped
        GridResample(; L, V, projection = hats)
    elseif key === :grid_exact
        GridResample(; L, V, projection = hats, density = f)
    else
        throw(ArgumentError("unknown construction :$key; the keys are " *
                            "$((STANDARD_METHODS..., REFERENCE_METHODS...))"))
    end
end

# The measures of one marker set of a row.
function _case_measures(case, mk, variant, reference)
    (; evenness = case_evenness(case, mk), outside = case_outside(case, mk),
        representation = case_representation(case, mk, variant; reference),
        charge_noise = case_charge_noise(case, mk, variant; reference),
        weights = case_weights(mk))
end

# The state of the solver of a remap: the fit of the projection and of the non-separable stages,
# the homotopy parameter of the gauge rezoning, nothing of the other constructions.
_solver_row(method::ShiftedRemap, info) = _solver_row(method.method, info)
function _solver_row(::Union{KickDriftRemap, NonSeparableRemap}, info)
    (; code = string(info.status.code), info.iterations,
        info.objective_ratio, tau = missing)
end
function _solver_row(::GaugeRezoning, info)
    (; code = missing, iterations = missing, objective_ratio = missing, tau = info.τ)
end
function _solver_row(::AbstractRemap, info)
    (; code = missing, iterations = missing, objective_ratio = missing, tau = missing)
end

"""
    run_case(method::AbstractRemap, case::StandardCase; variant = :positive,
             label = string(nameof(typeof(method))), warmup = false,
             reference = case_reference(case, variant))
    run_case(key::Symbol, case::StandardCase; variant = :positive, warmup = false,
             reference = case_reference(case, variant))

Run `remap(method, case.markers)` once (twice with `warmup = true`, the first untimed) and
return one row of every measure, a `NamedTuple` with the fields `label`, `family`,
`parameters`, `N`, `variant`, `status` (`:success` or `:failure`), `error_class` and
`error_message` (`nothing` on success), `time` (the seconds of the timed remap), `julia`,
`package` (the versions), `before` and `after` (each `(; evenness, outside, representation,
charge_noise, weights)` of [`case_evenness`](@ref), [`case_outside`](@ref),
[`case_representation`](@ref), [`case_charge_noise`](@ref) and [`case_weights`](@ref)),
`conservation` ([`case_conservation`](@ref)), `report` ([`conservation_report`](@ref) without
a grid), `displacement` ([`case_displacement`](@ref)) and `solver`, `(; code, iterations,
objective_ratio, tau)`: the name of the return code, the iterations and `F/F0` of the fit of
the projection and of the non-separable stages, the homotopy parameter `τ` of the gauge
rezoning, `missing` where a construction has none. `after` and the last four are `nothing` on
failure.

An exception of a remap that is a `DomainError`, an `ArgumentError`, a `SingularException` or a
`LAPACKException` is a failure, recorded with its class and its message; any other exception is
rethrown. A result with a non-finite position, velocity or weight is a failure of the class
`NonFiniteResult`. A remap that leaves markers outside the box is a success, and `after.outside`
counts them.

`reference` is a [`case_reference`](@ref) of `case`, which the measures before and after use;
pass one to run several methods on one case with one computation of the integrals of its
density. A reference of another variant than `variant` throws an `ArgumentError` before the
remap runs. The method of `key` is [`case_method`](@ref)`(key, case; variant)`, with the label
`string(key)`.
"""
function run_case(method::AbstractRemap, case::StandardCase; variant::Symbol = :positive,
        label = string(nameof(typeof(method))), warmup::Bool = false,
        reference = case_reference(case, variant))
    before = case.markers
    common = (; label, case.family, case.parameters, N = length(before), variant)
    versions = (; julia = string(VERSION), package = string(pkgversion(ParticleRemapping)))
    # the measures before the remap check the variant and the reference, so a call with a
    # wrong one throws before the remap runs
    measures_before = _case_measures(case, before, variant, reference)
    failure(class, message, time) = (; common..., status = :failure, error_class = class,
        error_message = message, time, versions..., before = measures_before,
        after = nothing, conservation = nothing, report = nothing, displacement = nothing,
        solver = nothing)
    start = time_ns()
    result = try
        warmup && remap(method, before)
        start = time_ns()
        remap(method, before)
    catch e
        e isa _REMAP_FAILURES || rethrow()
        return failure(string(nameof(typeof(e))), sprint(showerror, e),
            (time_ns() - start) / 1e9)
    end
    elapsed = (time_ns() - start) / 1e9
    after = result.markers
    nonfinite = count(.!(isfinite.(after.x) .& isfinite.(after.v) .& isfinite.(after.w)))
    nonfinite == 0 || return failure("NonFiniteResult",
        "$nonfinite of $(length(after)) markers have a non-finite position, velocity or " *
        "weight after the remap", elapsed)
    (; common..., status = :success, error_class = nothing, error_message = nothing,
        time = elapsed, versions..., before = measures_before,
        after = _case_measures(case, after, variant, reference),
        conservation = case_conservation(case, before, after),
        report = conservation_report(before, after),
        displacement = case_displacement(case, before, after),
        solver = _solver_row(method, result.info))
end

function run_case(key::Symbol, case::StandardCase; variant::Symbol = :positive,
        warmup::Bool = false, reference = case_reference(case, variant))
    run_case(case_method(key, case; variant), case; variant, label = string(key), warmup,
        reference)
end

export case_method, run_case
