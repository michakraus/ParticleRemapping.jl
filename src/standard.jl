# The standard experiments: two marker sets that every construction remaps with the same measures,
# the checkerboard and the long-time run in a two-wave electric field.

using .Remapping: CurrentQuantileMatching, DensityTarget, GaugeRezoning
using .Projection: HatProjection, KickDriftRemap, NonSeparableRemap
using .Redistribution: SplineRedistribution
using .SlicedRedistribution: SlicedQuantileRemap, cell_count_history, checkerboard_markers,
                             electric_wave_markers, electric_wave_simulation, periodic_wrap,
                             phase_cell_statistics, phase_quadrant_counts,
                             smoothed_density_difference

"""
    ShiftedRemap(method, shift, L)

The remap `method` of a construction on `[0, L)`, applied to markers on `[−L/2, L/2)`: the
positions are shifted by `shift` into `[0, L)` before the remap and back into `[−L/2, L/2)`
after it. The velocities and the weights are those of `method`. `remap` returns the `info` of
`method`.
"""
struct ShiftedRemap{M <: AbstractRemap, T <: Real} <: AbstractRemap
    method::M
    shift::T
    L::T
end

function ShiftedRemap(method::AbstractRemap, shift::Real, L::Real)
    ShiftedRemap(
        method, promote(float(shift), float(L))...)
end

function remap(s::ShiftedRemap, mk::Markers)
    r = remap(s.method, Markers(mod.(mk.x .+ s.shift, s.L), mk.v, mk.w; m = mk.m))
    after = r.markers
    RemapResult(
        Markers(periodic_wrap.(after.x .- s.shift, s.L / 2), after.v, after.w;
            m = after.m), r.info)
end

"""
    STANDARD_METHODS

The keys of the constructions that run the standard experiments, in the order of the manual:
`:projection` (kick–drift map), `:nonseparable` (non-separable stages), `:spline`
(spline-projection redistribution), `:gauge` (gauge rezoning), `:transport` (canonical transport
with the identity map) and `:sliced` (sliced redistribution).
"""
const STANDARD_METHODS = (:projection, :nonseparable, :spline, :gauge, :transport, :sliced)

function _check_variant(variant)
    variant in (:represented, :positive) ||
        throw(ArgumentError("variant must be :represented or :positive, got :$variant"))
end

"""
    checkerboard_density(x, v; variant = :represented)

A density of the checkerboard, with total weight 1 on `[−1, 1]²` and zero outside. With
`variant = :represented`, the density that the [`checkerboard_markers`](@ref) represent: `1/2`
in the quadrants `Q1 = [0, 1)²` and `Q3 = [−1, 0)²`, zero in `Q2` and `Q4`. With
`variant = :positive`, `1/4` on `[−1, 1]²`: the markers are a sample of it with a hole, since
their weights `1/N` are those of a sample of `1/4` with the uniform marker density `N/4` on
`[−1, 1]²`, and `Q2` and `Q4` hold no marker.
"""
function checkerboard_density(x::Real, v::Real; variant::Symbol = :represented)
    _check_variant(variant)
    abs(x) <= 1 && abs(v) <= 1 || return 0.0
    variant === :positive ? 0.25 : ((x >= 0) == (v >= 0) ? 0.5 : 0.0)
end

"""
    electric_wave_density(x, v; variant = :represented, steps = 500, Δt = 0.05, m = 1.2,
                          Lx = π, V0 = 1.4, Nx = 30, Nv = 30)

The density of the long-time run after `steps` steps of [`electric_wave_simulation`](@ref):
`f₀ = (1 + 0.35 cos x)(0.40 + v²) e^{−v²/1.4}` at the foot of the discrete characteristic
through `(x, v)`, found by the exact inverse of the drift–kick–drift steps, and normalised as
the weights of the [`electric_wave_markers`](@ref) (`Σ w = 1` on the `Nx × Nv` lattice of
`[−Lx, Lx) × [−V0, V0]`). With `variant = :represented`, `f₀` vanishes for `|v| > V0`, where
the lattice holds no marker: the density that the markers represent. With `variant = :positive`,
`f₀` is positive everywhere, and the markers sample it only where the characteristics start in
`|v| ≤ V0`.
"""
function electric_wave_density(x::Real, v::Real; variant::Symbol = :represented,
        steps::Integer = 500, Δt::Real = 0.05, m::Real = 1.2, Lx::Real = Float64(π),
        V0::Real = 1.4, Nx::Integer = 30, Nv::Integer = 30)
    _check_variant(variant)
    cut = variant === :represented
    f0(xx, vv) = cut && abs(vv) > V0 ? 0.0 :
                 (1.0 + 0.35 * cos(xx)) * (0.40 + vv^2) * exp(-vv^2 / 1.4)
    dx0, dv0 = 2Lx / Nx, 2V0 / Nv
    Z = sum(f0(-Lx + (i - 0.5) * dx0, -V0 + (j - 0.5) * dv0) for i in 1:Nx, j in 1:Nv) *
        dx0 *
        dv0
    E = SlicedRedistribution._electric_wave_field
    for step in steps:-1:1
        xh = x - 0.5 * Δt * v
        v -= Δt * E((step - 0.5) * Δt, xh) / m
        x = xh - 0.5 * Δt * v
    end
    f0(periodic_wrap(x, Lx), v) / Z
end

const _REMAP_FAILURES = Union{DomainError, ArgumentError, LinearAlgebra.SingularException,
    LinearAlgebra.LAPACKException}

# The number of markers of a reduced run, `nothing` at the full size `full`: the spline-projection
# redistribution keeps its defaults at the full size.
_reduced(N, full) = N == full ? nothing : N

# The value of `f()`, or the message of a failure of a construction on its marker set.
function _try_remap(f)
    try
        f(), nothing
    catch e
        e isa _REMAP_FAILURES || rethrow()
        nothing, sprint(showerror, e)
    end
end

"""
    standard_method(key, experiment; N = nothing, density = nothing)

The method of the construction `key` (one of [`STANDARD_METHODS`](@ref)) for the standard
experiment `experiment`, `:checkerboard` or `:wave`, with the parameters of the manual. `N` is
the number of markers of a reduced run, for which the spline-projection redistribution scales its
relaxation step to `0.006 · 16 / N`; `nothing`, the full size of the standard experiments, keeps
its defaults. `density(x, v)` is the density on the natural box of the
experiment that the constructions with a target density receive, the projection and the
transport; the other constructions ignore it. Without it, the projection projects the density of
the markers and the transport takes the represented density of the experiment.
"""
function standard_method(key::Symbol, experiment::Symbol; N = nothing, density = nothing)
    experiment in (:checkerboard, :wave) ||
        throw(ArgumentError("experiment must be :checkerboard or :wave, got :$experiment"))
    cb = experiment === :checkerboard
    L = cb ? 2.0 : 2π
    shift = L / 2
    hats = cb ? HatProjection(12, 13) : HatProjection(16, 17)
    shifted = density === nothing ? nothing :
              (xs, vs) -> density.(periodic_wrap.(xs .- shift, L / 2), vs)
    if key === :projection
        ShiftedRemap(
            KickDriftRemap(; L, V = cb ? 1.0 : 3.0, projection = hats, maxit = 20,
                density = shifted),
            shift,
            L)
    elseif key === :nonseparable
        ShiftedRemap(
            NonSeparableRemap(; L, V = cb ? 1.0 : 3.0, nξ = 6, nv = 7, nℓ = 3,
                projection = hats, maxit = 5, density = shifted),
            shift,
            L)
    elseif key === :spline
        relax = N === nothing ? (;) : (; relax_step = 0.006 * 16 / N)
        ShiftedRemap(
            SplineRedistribution(; Lx = L, Vmax = cb ? 1.5 : 3.5, σx = L / 4,
                σv = cb ? 0.5 : 1.2, relax...),
            shift,
            L)
    elseif key === :gauge
        GaugeRezoning(PeriodicGrid(L, cb ? 8 : 16))
    elseif key === :transport
        grid = PeriodicGrid(L, cb ? 8 : 16)
        f = density === nothing ? (cb ? checkerboard_density : electric_wave_density) :
            density
        CurrentQuantileMatching(grid, (x, p) -> (copy(x), copy(p)),
            DensityTarget(grid, (x, v) -> f(periodic_wrap(x, L / 2), v); vmax = cb ? 1.0 :
                                                                                3.5))
    elseif key === :sliced
        cb ?
        SlicedQuantileRemap(; Lx = 1, Vmax = 1, B_v = 10, B_x = 20, hx = 0.15,
            hv = 0.15) : SlicedQuantileRemap()
    else
        throw(ArgumentError("unknown construction :$key; the keys are $STANDARD_METHODS"))
    end
end

"""
    standard_checkerboard(key; variant = :represented, N = 400, seed = 505)

The standard checkerboard experiment for the construction `key`: the
[`checkerboard_markers`](@ref) remapped by `remap(standard_method(key, :checkerboard), ·)`; the
projection and the transport receive the [`checkerboard_density`](@ref) of `variant`
(`:represented`, which vanishes in `Q2` and `Q4`, or `:positive`). The box is periodic in `x` on
`[−1, 1)`, and `after` holds the positions wrapped into it.

Returns `before`, `after` (`nothing` if the construction fails), the `error` message
(`nothing` if it does not), the `result` of `remap`, the quadrant counts before and after, the
weight `hole_weight` in `Q2 ∪ Q4` after, the largest periodic displacement `max_shift` of a
position and the largest `|v|` after, the [`conservation_report`](@ref) and the cell
statistics on `8 × 8` cells of `[−1, 1]²` before and after.
"""
function standard_checkerboard(
        key::Symbol; variant::Symbol = :represented, N::Integer = 400,
        seed::Integer = 505)
    _check_variant(variant)
    before = checkerboard_markers(; seed, N)
    method = standard_method(key, :checkerboard; N = _reduced(N, 400),
        density = (x, v) -> checkerboard_density(x, v; variant))
    result, error = _try_remap(() -> remap(method, before))
    stats(mk) = phase_cell_statistics(mk.x, mk.v; gx = 8, gv = 8, Lx = 1, Vmax = 1)
    common = (; key, variant, N, before, error, result,
        quadrants_before = phase_quadrant_counts(before.x, before.v),
        stats_before = stats(before))
    result === nothing && return (; common..., after = nothing)
    after = Markers(
        periodic_wrap.(result.markers.x, 1), result.markers.v, result.markers.w;
        m = result.markers.m)
    hole = (after.x .< 0) .!= (after.v .< 0)
    (; common..., after, quadrants_after = phase_quadrant_counts(after.x, after.v),
        hole_weight = sum(after.w[hole]; init = 0.0),
        max_shift = maximum(abs, periodic_wrap.(after.x .- before.x, 1)),
        max_v = maximum(abs, after.v), report = conservation_report(before, after),
        stats_after = stats(after))
end

"""
    standard_wave_markers(; Nx = 30, Nv = 30, steps = 500)

The markers of the standard long-time experiment: the `Nx × Nv` [`electric_wave_markers`](@ref)
after `steps` steps of [`electric_wave_simulation`](@ref) (`T = 25` for the defaults).
"""
function standard_wave_markers(; Nx::Integer = 30, Nv::Integer = 30, steps::Integer = 500)
    electric_wave_simulation(electric_wave_markers(; Nx, Nv); steps).markers
end

"""
    standard_long_time(key; variant = :represented, Nx = 30, Nv = 30, steps = 500)

The standard long-time experiment for the construction `key`: the
[`standard_wave_markers`](@ref) remapped once by `remap(standard_method(key, :wave), ·)`; the
projection and the transport receive the [`electric_wave_density`](@ref) of `variant`.

Returns `before`, `after` (`nothing` if the construction fails; the positions wrapped into
`[−π, π)`), the `error` message, the `result` of `remap`, the cell statistics of
[`phase_cell_statistics`](@ref) before and after, the [`smoothed_density_difference`](@ref) of
the two marker sets, the largest `|v|` after and the [`conservation_report`](@ref).
"""
function standard_long_time(key::Symbol; variant::Symbol = :represented, Nx::Integer = 30,
        Nv::Integer = 30, steps::Integer = 500)
    _check_variant(variant)
    before = standard_wave_markers(; Nx, Nv, steps)
    method = standard_method(key, :wave; N = _reduced(Nx * Nv, 900),
        density = (x, v) -> electric_wave_density(x, v; variant, steps, Nx, Nv))
    result, error = _try_remap(() -> remap(method, before))
    stats(mk) = phase_cell_statistics(mk.x, mk.v)
    common = (; key, variant, N = Nx * Nv, before, error, result,
        stats_before = stats(before))
    result === nothing && return (; common..., after = nothing)
    after = Markers(
        periodic_wrap.(result.markers.x, π), result.markers.v, result.markers.w;
        m = result.markers.m)
    (; common..., after, stats_after = stats(after),
        l2_difference = smoothed_density_difference(before.x, before.v, before.w, after.x,
            after.v, after.w),
        max_v = maximum(abs, after.v), report = conservation_report(before, after))
end

"""
    standard_remap_history(key; variant = :represented, Nx = 30, Nv = 30, steps = 500,
                           remap_every = 100, record_every = 10)

The long-time run with a remap by `standard_method(key, :wave)` after every `remap_every` steps,
beside the run without remaps; the projection and the transport receive the
[`electric_wave_density`](@ref) of `variant` at the time of each remap. Returns both runs of
[`electric_wave_simulation`](@ref) (`plain`, `remapped`), their [`cell_count_history`](@ref)
(`h_plain`, `h_remap`), the mean relative standard deviation of the cell counts for
`t ≥ remap_every Δt` of each (`rstd_plain`, `rstd_remap`), and the `error` message if a remap of
the construction fails (then `remapped`, `h_remap` and `rstd_remap` are `nothing`).
"""
function standard_remap_history(key::Symbol; variant::Symbol = :represented,
        Nx::Integer = 30, Nv::Integer = 30, steps::Integer = 500, remap_every::Integer = 100,
        record_every::Integer = 10)
    _check_variant(variant)
    initial = electric_wave_markers(; Nx, Nv)
    Δt = 0.05
    plain = electric_wave_simulation(initial; steps, record_every)
    h_plain = cell_count_history(plain.snapshots, plain.times)
    late(h) = sum(h.rstd[h.times .>= remap_every * Δt]) /
              count(>=(remap_every * Δt), h.times)
    method(t) = standard_method(key, :wave; N = _reduced(Nx * Nv, 900),
        density = (x, v) -> electric_wave_density(x, v; variant,
            steps = round(Int, t / Δt), Nx, Nv))
    remapped, error = _try_remap(() -> electric_wave_simulation(
        initial; steps, record_every,
        remap_method = method, remap_every))
    common = (; key, variant, plain, h_plain, rstd_plain = late(h_plain), error)
    remapped === nothing &&
        return (; common..., remapped, h_remap = nothing, rstd_remap = nothing)
    h_remap = cell_count_history(remapped.snapshots, remapped.times)
    (; common..., remapped, h_remap, rstd_remap = late(h_remap))
end

export ShiftedRemap, STANDARD_METHODS, checkerboard_density, electric_wave_density,
       standard_method, standard_checkerboard, standard_wave_markers, standard_long_time,
       standard_remap_history

# the standard test set: the cases, the measures, the reference and the driver
include("standard/cases.jl")
include("standard/measures.jl")
include("standard/reference.jl")
include("standard/driver.jl")
include("standard/pic.jl")
