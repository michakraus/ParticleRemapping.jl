# The downstream family P of the standard test set: self-consistent Vlasov–Poisson runs of
# VlasovMethods' particle model with a remap every R steps. P1 and the constructions of a run need
# no VlasovMethods; P2 and the runner are methods of the extension ParticleRemappingVlasovMethodsExt.

"""
    two_stream_case(; N = 900, loader = :lattice, seed = 0, V = 6.0)

The case P1 of family `:P1`, the two-stream instability of VlasovMethods' reference: the density
`f₀(x, v) = (1 + α cos(k x)) v² e^{−v²/2}/√(2π)`, `α = 10⁻³`, `k = 0.5`, on `[0, L)`,
`L = 4π`, in case coordinates, so that the case's density at `x` is `f₀(x + L/2, v)`, on the box
`V`, with mass `1` and the experiment `:wave`. The markers are truncated at `|v| ≤ VL = 5`: the
represented density is `f₀(x + L/2, v)` for `|v| ≤ VL` and zero beyond, the positive density
`f₀(x + L/2, v)` on the whole box.

With `loader = :lattice`, the markers are the `n × n` cell centres
`(−L/2 + (i − ½) dx, −VL + (j − ½) dv)`, `n = √N`, `dx = L/n`, `dv = 2VL/n`, `i` outer, `j`
inner, with the weights `c f₀(x + L/2, v) dx dv` and `c` such that `Σ w = L`, the normalisation of
the model, whose plasma frequency is `1`. With `loader = :sample`, `N` markers of weight `L/N`
are drawn with `Xoshiro(seed)`: per attempt `x = −L/2 + L rand(rng)`, rejected when
`(1 + α) rand(rng) > 1 + α cos(k (x + L/2))`, then `s = √(r₁² + r₂² + r₃²)` from three
`randn(rng)` (`v² e^{−v²/2}/√(2π)` is the density of `±χ₃`), rejected when `s > VL` (the attempt
starts again), then `v = −s` when `rand(rng) < 0.5` and `v = s` otherwise. The sample is the noisy
loader of the family: its charge noise, of the order of `1/√N`, lies far above `α`.

`parameters` holds the run of [`run_pic`](@ref): `cells = 16` and `degree = 3` of the potential,
the step `dt = 0.1`, `steps = 400`, the fit `window = (10.0, 25.0)` and the reference `rate`,
VlasovMethods' `TWO_STREAM_RATE = 0.2592499317`.

An `ArgumentError` refuses a loader other than `:lattice` and `:sample`, a lattice `N` that is
no square `n²` with `n ≥ 1` and a lattice seed other than `0`, a sample `N < 1` and a negative
sample seed; [`StandardCase`](@ref) refuses a box `V` that does not hold every marker.
"""
function two_stream_case(; N::Integer = 900, loader::Symbol = :lattice, seed::Integer = 0,
        V::Real = 6.0)
    L, VL, α, k = 4π, 5.0, 1e-3, 0.5
    f0(x, v) = (1 + α * cos(k * (x + L / 2))) * v^2 * exp(-v^2 / 2) / sqrt(2π)
    markers = if loader === :lattice
        n = _lattice_side(N)
        n === nothing && throw(ArgumentError(
            "two_stream_case :lattice needs N = n^2 with n >= 1, got N = $N"))
        seed == 0 || throw(ArgumentError(
            "two_stream_case :lattice is deterministic; its only seed is 0, got seed = $seed"))
        x, v = _cell_centres(L, VL, n)
        raw = f0.(x, v) .* (L / n * (2VL / n))
        Markers(x, v, raw .* (L / sum(raw)); m = 1.0)
    elseif loader === :sample
        N >= 1 && seed >= 0 || throw(ArgumentError(
            "two_stream_case :sample needs N >= 1 and a seed >= 0, got N = $N, seed = $seed"))
        _two_stream_sample(N, seed, L, VL, α, k)
    else
        throw(ArgumentError("two_stream_case needs loader :lattice or :sample, got :$loader"))
    end
    represented(x, v) = abs(v) <= VL ? f0(x, v) : 0.0
    parameters = (; N, loader, seed, V = Float64(V), VL, α, k, cells = 16, degree = 3,
        dt = 0.1, steps = 400, window = (10.0, 25.0), rate = 0.2592499317)
    StandardCase(:P1, parameters, markers, L, V, represented, f0, :wave)
end

# N markers of weight L/N of the truncated two-stream density, drawn with Xoshiro(seed).
function _two_stream_sample(N, seed, L, VL, α, k)
    rng = Xoshiro(seed)
    x, v = Float64[], Float64[]
    while length(x) < N
        xa = -L / 2 + L * rand(rng)
        (1 + α) * rand(rng) > 1 + α * cos(k * (xa + L / 2)) && continue
        s = sqrt(randn(rng)^2 + randn(rng)^2 + randn(rng)^2)
        s > VL && continue
        push!(x, xa)
        push!(v, rand(rng) < 0.5 ? -s : s)
    end
    Markers(x, v, fill(L / N, N); m = 1.0)
end

"""
    bump_on_tail_case(; N = 10_000, seed = 1234, V = 8.0)

The case P2 of family `:P2`, the bump-on-tail run of VlasovMethods' scripts: the markers of
VlasovMethods' own draw, `Random.seed!(seed)` and
`initialize!(ParticleDistribution(1, 1, N), BumpOnTail(κ = 0.3))`, with their positions moved
into the case's box by `mod(x, L) − L/2`, on `L = 2π/0.3` and the box `V`, with mass `1` and the
experiment `:wave`. The state of the default generator is saved before the draw and restored
after it, so the caller's stream does not change. The one density in both fields is `f(x + L/2,
v)`, `f = (1 − ε cos κx) ((1 − α) e^{−v²/2} + (α/σ) e^{−(v − v₀)²/(2σ²)})/√(2π)` of
`BumpOnTail`, `ε = 0.03`, `κ = 0.3`, `α = 0.1`, `σ = 0.5`, `v₀ = 4.5`.

`parameters` holds the run of [`run_pic`](@ref): `cells = 16` and `degree = 3` of the potential,
the step `dt = 0.1`, `steps = 250`, the fit `window = (5.0, 17.5)`, and `rate = nothing`: P2 has
no reference rate. An `ArgumentError` refuses `N < 1` and a negative seed; [`StandardCase`](@ref)
refuses a box `V` that does not hold every marker.

Needs `using VlasovMethods`, which loads the extension `ParticleRemappingVlasovMethodsExt`.
"""
function bump_on_tail_case end

"""
    run_pic(case::StandardCase; steps = case.parameters.steps, variant = :positive)
    run_pic(method::AbstractRemap, case::StandardCase; R, steps = case.parameters.steps,
            variant = :positive, label = string(nameof(typeof(method))))
    run_pic(key::Symbol, case::StandardCase; R, steps = case.parameters.steps,
            variant = :positive)

March the markers of a case of the family `:P1` or `:P2` by VlasovMethods' leapfrog
(`full_leapfrog_step!` over the `PoissonField` of `PeriodicBasisSpline((0, L), degree + 1,
cells)`, positions `x + L/2` on `[0, L)`, the step `dt` and the kick time `(k − ½) dt` of step
`k`), and return one row, a `NamedTuple`. The field energy `W[k + 1]` is that of a fresh deposit
of the state after step `k` (`W[1]` of the initial markers), and `rate` is half the
least-squares slope of `log W` over the times in `parameters.window`, since the field energy
grows as `e^{2γt}`; `nothing` when the history ends before the window does.

With a `method`, the run without remaps of the same case runs first (the reference), then the
march again, and after every step `k` with `k % R == 0` the state is remapped: the markers
`before` (positions wrapped into the case's box) go to `remap(method, before)`, and the
positions, velocities and weights of the result go back into the march, which carries nothing
else across the remap; `W[k + 1]` is then the field energy of the remapped state. The method of
`key` is [`pic_method`](@ref)`(key, case)`, with the label `string(key)`.

The row has the fields `label` (`"none"` without remaps), `family`, `parameters`, `N` (the
markers of the case), `R`, `steps`, `variant`, `status` (`:success` or `:failure`),
`error_class` and `error_message` (the class and the message of a failure), `failed_at` (its step), `t` (`t[k + 1] = k dt`), `W`, `rate`, `W_reference` and
`rate_reference` (of the reference), `charge_noise` ([`case_charge_noise`](@ref) of the initial
markers and `variant`), `charge_deviation` (`‖M_k − M_k^ref‖₂ / ‖M_k^ref‖₂` of the grid charge
`M` on `PeriodicGrid(L, 16)` of the state after step `k` and `M^ref` of the reference, `NaN`
past the last step of a reference that failed), `max_v` (the largest `|v|` of the initial
markers and after every push and every remap), `remaps`, `time_push` and `time_remap`
(the seconds of all steps and of all remaps), `julia`, `package` and `vlasovmethods` (the
versions). Without remaps `R`, `W_reference`, `rate_reference` and `charge_deviation` are
`nothing` and `remaps` is empty. `remaps` holds one `(; step, time, W_before, W_after,
outside_push, outside, conservation, weights, info)` per remap: the field energies before and
after it, [`case_outside`](@ref) of the markers before the remap and after it,
[`case_conservation`](@ref), [`case_weights`](@ref) of the markers after it, the seconds `time`
of the remap and its own `info`.

The row classifies nothing beyond the status: markers outside the box after a push or a remap
are counted and the run goes on. A remap that throws a `DomainError`, an `ArgumentError`, a
`SingularException` or a `LAPACKException`, a remap result with a non-finite position, velocity
or weight (`"NonFiniteResult"`) and a field energy that is not finite after a step
(`"NonFiniteField"`) are failures: the run stops, and `t`, `W` and `charge_deviation` hold the
states up to that step (the state before a failed remap), `remaps` the remaps that succeeded.
Every other exception is rethrown. An `ArgumentError` refuses, before any step, a case of
another family, `steps < 1` and an `R` outside `1 ≤ R ≤ steps`.

Needs `using VlasovMethods`, which loads the extension `ParticleRemappingVlasovMethodsExt`.
"""
function run_pic end

"""
    PIC_METHODS

The keys of the constructions that [`pic_method`](@ref) builds for a remap in the middle of a
particle-in-cell run: `:projection`, `:nonseparable`, `:spline`, `:gauge`, `:sliced` and
`:grid_lumped`. The transport and `:grid_exact` need the exact density, which a run knows only
at `t = 0`.
"""
const PIC_METHODS = (:projection, :nonseparable, :spline, :gauge, :sliced, :grid_lumped)

"""
    pic_method(key, case::StandardCase)

The method of the construction `key` (one of [`PIC_METHODS`](@ref)) for a remap in the middle of
a run of `case`, where the exact density is no longer known. For `:spline`, `:gauge`, `:sliced`
and `:grid_lumped`, [`case_method`](@ref)`(key, case)`, which reads no density for them. For
`:projection` and `:nonseparable`, the parameters of [`standard_method`](@ref)`(key, :wave)`
with the box of the case, `L = case.L` and `V = case.V`, and `density = nothing`, the projection
of the markers' own density, inside a [`ShiftedRemap`](@ref) by `L/2`. An `ArgumentError` refuses
`:transport` and `:grid_exact` and every other key.
"""
function pic_method(key::Symbol, case::StandardCase)
    L = case.L
    if key === :projection
        ShiftedRemap(
            KickDriftRemap(; L, V = case.V, projection = HatProjection(16, 17), maxit = 20),
            L / 2, L)
    elseif key === :nonseparable
        ShiftedRemap(
            NonSeparableRemap(; L, V = case.V, nξ = 6, nv = 7, nℓ = 3,
                projection = HatProjection(16, 17), maxit = 5),
            L / 2, L)
    elseif key in PIC_METHODS
        case_method(key, case)
    elseif key in (:transport, :grid_exact)
        throw(ArgumentError(
            "pic_method(:$key) needs the exact density, which a particle-in-cell run knows " *
            "only at t = 0; the keys are $PIC_METHODS"))
    else
        throw(ArgumentError("unknown construction :$key; the keys are $PIC_METHODS"))
    end
end

export two_stream_case, bump_on_tail_case, run_pic, PIC_METHODS, pic_method
