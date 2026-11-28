# The runner of the family P: VlasovMethods' leapfrog march of a case, the field-energy history and
# its growth rate, and a remap every R steps against the march without remaps.

function _check_pic_case(case::StandardCase)
    case.family in (:P1, :P2) || throw(ArgumentError(
        "run_pic needs a case of the family :P1 or :P2, got :$(case.family)"))
    case
end

# The default number of steps; it checks the family first, because only a case of family P holds
# the parameters of a run.
_pic_steps(case::StandardCase) = _check_pic_case(case).parameters.steps

function _check_steps(steps::Integer)
    steps >= 1 || throw(ArgumentError("run_pic needs steps >= 1, got steps = $steps"))
end

# Half the least-squares slope of log W over the times in `window`, since the field energy grows
# as e^{2γt}; `nothing` when the history ends before the window does.
function _growth_rate(t, W, window)
    t[end] >= window[2] || return nothing
    sel = findall(s -> window[1] <= s <= window[2], t)
    τ, y = t[sel], log.(W[sel])
    d = τ .- sum(τ) / length(τ)
    sum(d .* (y .- sum(y) / length(y))) / sum(abs2, d) / 2
end

# One march of `case` over `steps` steps of VlasovMethods' leapfrog, with a remap by `method`
# after every step k with k % R == 0 (none for `method === nothing`). The positions of the march
# are those of the case shifted by L/2 onto [0, L), left unwrapped; the periodic basis of the
# potential reduces them. `Mref` holds the grid charges of the reference march, against which
# the charge deviation is measured; `nothing` for the reference itself.
function _march(case::StandardCase, steps::Integer, method, R::Integer, Mref)
    p = case.parameters
    L, dt = case.L, p.dt
    hL = L / 2
    x = case.markers.x .+ hL
    v = copy(case.markers.v)
    w = copy(case.markers.w)
    m = case.markers.m
    a = zeros(length(x))
    field = PoissonField(Potential(PeriodicBasisSpline((0.0, L), p.degree + 1, p.cells)))
    history = (;
        t = Float64[], W = Float64[], M = Vector{Float64}[], deviation = Float64[],
        grid = PeriodicGrid(L, 16), dt, hL, Mref)
    _record!(history, 0, _fresh_energy(field, x, w, 0.0), x, v, w)
    max_v = 0.0
    remaps = NamedTuple[]
    time_push = time_remap = 0.0
    failure = nothing
    for k in 1:steps
        time_push += @elapsed full_leapfrog_step!(x, v, a, w, field, dt, (k - 0.5) * dt)
        max_v = max(max_v, maximum(abs, v))
        if method !== nothing && k % R == 0
            # the constructor converts its arrays, so `before` holds copies of the state
            before = Markers(periodic_wrap.(x .- hL, hL), v, w; m)
            W_before = _fresh_energy(field, x, w, k * dt)
            start = time_ns()
            result = try
                remap(method, before)
            catch e
                e isa _REMAP_FAILURES || rethrow()
                failure = (; step = k, class = string(nameof(typeof(e))),
                    message = sprint(showerror, e))
                nothing
            end
            elapsed = (time_ns() - start) / 1e9
            time_remap += elapsed
            if result !== nothing
                after = result.markers
                n = count(.!(isfinite.(after.x) .& isfinite.(after.v) .& isfinite.(after.w)))
                n == 0 || (failure = (; step = k, class = "NonFiniteResult",
                    message = "$n of $(length(after)) markers have a non-finite position, " *
                              "velocity or weight after the remap"))
            end
            if failure !== nothing
                _record!(history, k, W_before, x, v, w)  # the state before the failed remap
                break
            end
            x = after.x .+ hL
            v = copy(after.v)
            w = copy(after.w)
            length(a) == length(x) || (a = zeros(length(x)))
            max_v = max(max_v, maximum(abs, v))
            Wk = _fresh_energy(field, x, w, k * dt)
            push!(remaps,
                (; step = k, time = elapsed, W_before, W_after = Wk,
                    outside_push = case_outside(case, before),
                    outside = case_outside(case, after),
                    conservation = case_conservation(case, before, after),
                    weights = case_weights(after), result.info))
        else
            Wk = _fresh_energy(field, x, w, k * dt)
        end
        _record!(history, k, Wk, x, v, w)
        if !isfinite(Wk)
            failure = (; step = k, class = "NonFiniteField",
                message = "the field energy is not finite after step $k")
            break
        end
    end
    (; history.t, history.W, history.M, history.deviation, max_v, remaps, time_push,
        time_remap, failure)
end

# The field energy of a fresh deposit of the state. The field after a leapfrog step is that of
# the half-drifted positions, so the energy of a state needs its own deposit.
function _fresh_energy(field, x, w, time)
    update!(field, x, w, time)
    energy(field)
end

# Record the state after step k: its time, its field energy `Wk`, its grid charge and, against a
# reference, the charge deviation (NaN past the end of a reference that failed before step k).
function _record!(h, k, Wk, x, v, w)
    push!(h.t, k * h.dt)
    push!(h.W, Wk)
    Mk = grid_moment(h.grid, x .- h.hL, v, w, 0)
    push!(h.M, Mk)
    if h.Mref !== nothing
        push!(h.deviation,
            k + 1 <= length(h.Mref) ? norm(Mk - h.Mref[k + 1]) / norm(h.Mref[k + 1]) : NaN)
    end
    nothing
end

# The row of a march, with the march without remaps of the same case as `reference`.
function _pic_row(case, label, R, steps, variant, charge_noise, march, reference)
    window = case.parameters.window
    f = march.failure
    (; label, case.family, case.parameters, N = length(case.markers), R, steps, variant,
        status = f === nothing ? :success : :failure,
        error_class = f === nothing ? nothing : f.class,
        error_message = f === nothing ? nothing : f.message,
        failed_at = f === nothing ? nothing : f.step,
        march.t, march.W, rate = _growth_rate(march.t, march.W, window),
        W_reference = reference === nothing ? nothing : reference.W,
        rate_reference = reference === nothing ? nothing :
                         _growth_rate(reference.t, reference.W, window),
        charge_noise, charge_deviation = reference === nothing ? nothing : march.deviation,
        march.max_v, march.remaps, march.time_push, march.time_remap,
        julia = string(VERSION), package = string(pkgversion(ParticleRemapping)),
        vlasovmethods = string(pkgversion(VlasovMethods)))
end

function ParticleRemapping.run_pic(case::StandardCase; steps::Integer = _pic_steps(case),
        variant::Symbol = :positive)
    _check_pic_case(case)
    _check_steps(steps)
    charge_noise = case_charge_noise(case, case.markers, variant)  # checks the variant
    march = _march(case, steps, nothing, 0, nothing)
    _pic_row(case, "none", nothing, steps, variant, charge_noise, march, nothing)
end

function ParticleRemapping.run_pic(method::AbstractRemap, case::StandardCase; R::Integer,
        steps::Integer = _pic_steps(case), variant::Symbol = :positive,
        label = string(nameof(typeof(method))))
    _check_pic_case(case)
    _check_steps(steps)
    1 <= R <= steps || throw(ArgumentError(
        "run_pic needs 1 <= R <= steps, got R = $R, steps = $steps"))
    charge_noise = case_charge_noise(case, case.markers, variant)
    reference = _march(case, steps, nothing, 0, nothing)
    march = _march(case, steps, method, R, reference.M)
    _pic_row(case, label, R, steps, variant, charge_noise, march, reference)
end

function ParticleRemapping.run_pic(key::Symbol, case::StandardCase; R::Integer,
        steps::Integer = _pic_steps(case), variant::Symbol = :positive)
    ParticleRemapping.run_pic(pic_method(key, case), case; R, steps, variant,
        label = string(key))
end
