# The cases of the standard test set: marker sets on one box each, with their exact densities.
# Every case lies on Ω = [−L/2, L/2) × [−V, V]; the constructions on [0, L) see it through a
# `ShiftedRemap`.

using Random: MersenneTwister, Xoshiro

const _CASE_FAMILIES = (:W, :C, :E, :S1, :S2, :S3, :S4)

# The box of a case must be finite and positive: with V = NaN no marker passes |v| <= V, with
# V = Inf every one does, and with V <= 0 the rejection of S1 at |v| > V accepts no draw.
function _check_box(L::Float64, V::Float64)
    isfinite(L) && L > 0 && isfinite(V) && V > 0 || throw(ArgumentError(
        "StandardCase needs a finite L > 0 and V > 0, got L = $L, V = $V"))
end

"""
    StandardCase(family, parameters, markers, L, V, represented, positive, experiment)

A case of the standard test set: the `markers` of the family `family` (`:W`, `:C`, `:E`, `:S1`,
`:S2`, `:S3` or `:S4`) built with `parameters`, on the box `Ω = [−L/2, L/2) × [−V, V]`, with
the two scalar densities `represented(x, v)` and `positive(x, v)` (see [`case_density`](@ref);
one function in both fields where the family has one density) and the parameter set
`experiment` (`:wave` or `:checkerboard`) of [`standard_method`](@ref) that
[`case_method`](@ref) takes. The mass is `markers.m`.

`L` and `V` are stored as `Float64`. An `ArgumentError` refuses, in this order, a box that is
not finite and positive, a family or an experiment outside the lists, an empty marker set, a
marker with `|v| > V` (or a velocity `NaN`), and a position outside `[−L/2, L/2)`.
"""
struct StandardCase{P <: NamedTuple, M <: Markers, R, Q}
    family::Symbol
    parameters::P
    markers::M
    L::Float64
    V::Float64
    represented::R
    positive::Q
    experiment::Symbol

    function StandardCase(family::Symbol, parameters::P, markers::M, L::Real, V::Real,
            represented::R, positive::Q, experiment::Symbol) where {P <: NamedTuple,
            M <: Markers, R, Q}
        L, V = Float64(L), Float64(V)
        _check_box(L, V)  # first, because the checks of the markers cannot see a bad box
        family in _CASE_FAMILIES || throw(ArgumentError(
            "unknown family :$family; the families are $_CASE_FAMILIES"))
        experiment in (:wave, :checkerboard) || throw(ArgumentError(
            "unknown experiment :$experiment; the experiments are (:wave, :checkerboard)"))
        length(markers) >= 1 ||
            throw(ArgumentError("StandardCase needs at least one marker"))
        all(v -> abs(v) <= V, markers.v) || throw(ArgumentError(
            "StandardCase needs every |v| <= V = $V, got max |v| = $(maximum(abs, markers.v))"))
        all(x -> -L / 2 <= x < L / 2, markers.x) || throw(ArgumentError(
            "StandardCase needs every x in [-L/2, L/2) with L = $L"))
        new{P, M, R, Q}(family, parameters, markers, L, V, represented, positive,
            experiment)
    end
end

"""
    case_density(case::StandardCase, variant)

The density of `case` of `variant`: `case.represented` for `:represented`, the density that the
markers represent, and `case.positive` for `:positive`, a density positive on the box that the
markers sample with holes where the two differ. Both are the same function for the families E
and S.
"""
function case_density(case::StandardCase, variant::Symbol)
    _check_variant(variant)
    variant === :represented ? case.represented : case.positive
end

"""
    STANDARD_LADDERS

The rungs, sizes and seeds of the families of the standard test set: `W` the wave ladder
([`wave_case`](@ref)) over the times `T`, the sizes `N` and the seeds; `C` the hole ladder
([`hole_case`](@ref)) over the orders `k` of the checkerboard; `E` the plain and the deformed
lattice ([`lattice_case`](@ref)); `S` the stress cases ([`stress_case`](@ref)), whose seeds
apply to `:S1` and `:S3`, while `:S2` and `:S4` have the one seed `0`.
"""
const STANDARD_LADDERS = (
    W = (N = (16, 144, 576, 900), T = (0.0, 1.25, 2.5, 5.0, 10.0, 15.0, 25.0), seed = 0:4),
    C = (k = (1, 2, 4, 8), N = 400, seed = 505:509),
    E = (N = (144, 576), deformed = (false, true)),
    S = (name = (:S1, :S2, :S3, :S4), N = (144, 576), seed = 0:4))

# The side n of an n × n lattice of N markers, or `nothing` when N is no square n² with n ≥ 1.
function _lattice_side(N::Integer)
    N >= 1 || return nothing
    n = isqrt(N)
    n^2 == N ? n : nothing
end

# The positions and the velocities of the n × n cell centres of [−L/2, L/2) × [−V, V], i outer,
# j inner.
function _cell_centres(L, V, n)
    dx, dv = L / n, 2V / n
    ([-L / 2 + (i - 0.5) * dx for i in 1:n for j in 1:n],
        [-V + (j - 0.5) * dv for i in 1:n for j in 1:n])
end

"""
    wave_case(; N = 144, T = 5.0, seed = 0, V = 3.0)

The case `W(N, T)` of family `:W`: `N = n²` markers of the long run in the two-wave field
([`electric_wave_simulation`](@ref)) after `steps = T / 0.05` steps, on the box `L = 2π` and
`V`, with mass `1.2` and the experiment `:wave`. Every rung, size and seed takes the one box
`V`; a `V` that does not hold every marker is refused by [`StandardCase`](@ref). At `T = 0` the
markers fill only `|v| ≤ 1.4` of the box, so the rung `T = 0` is the first rung of the ladder,
not an equidistribution of its box.

Seed `0` pushes the `n × n` [`electric_wave_markers`](@ref). A seed `s ≥ 1` draws with
`Xoshiro(s)`, for each lattice point in the order of `electric_wave_markers`, the offsets `ux`,
then `uv`, moves the point to `(−π + (ix − 1 + ux) dx0, −1.4 + (iv − 1 + uv) dv0)`, weights it
by `f₀ dx0 dv0` at the moved point, normalises the weights to `Σ w = 1`, and pushes it. The
densities are [`electric_wave_density`](@ref) of the same `steps` and lattice, for the two
variants; they differ wherever the characteristic through `(x, v)` starts at `|v| > 1.4`.

An `ArgumentError` refuses an `N` that is no square `n²` with `n ≥ 1`, a `T` that is not
finite, is negative, lies off the step grid by more than `1e-9` steps or needs more steps than
an `Int` holds, and a negative seed.
"""
function wave_case(; N::Integer = 144, T::Real = 5.0, seed::Integer = 0, V::Real = 3.0)
    n = _lattice_side(N)
    n === nothing &&
        throw(ArgumentError("wave_case needs N = n^2 with n >= 1, got N = $N"))
    k = T / 0.05
    # k < typemax(Int): a step count that `round(Int, k)` can hold
    isfinite(T) && T >= 0 && k < typemax(Int) && abs(k - round(k)) <= 1e-9 ||
        throw(ArgumentError(
            "wave_case needs T >= 0 on the step grid T = 0.05 k, got T = $T"))
    seed >= 0 || throw(ArgumentError("wave_case needs a seed >= 0, got seed = $seed"))
    steps = round(Int, k)
    initial = seed == 0 ? electric_wave_markers(; Nx = n, Nv = n) :
              _stratified_wave_markers(n, seed)
    markers = electric_wave_simulation(initial; steps).markers
    represented(x, v) = electric_wave_density(x, v; variant = :represented, steps,
        Nx = n, Nv = n)
    positive(x, v) = electric_wave_density(x, v; variant = :positive, steps, Nx = n, Nv = n)
    StandardCase(:W, (; N, T = Float64(T), seed, steps, V = Float64(V)), markers, 2π, V,
        represented, positive, :wave)
end

# The n × n lattice of `electric_wave_markers`, each point moved inside its cell by the offsets
# that `Xoshiro(seed)` draws, weighted by f₀ at the moved point.
function _stratified_wave_markers(n::Integer, seed::Integer)
    rng = Xoshiro(seed)
    Lx, V0 = Float64(π), 1.4
    dx0, dv0 = 2Lx / n, 2V0 / n
    f0(x, v) = (1.0 + 0.35 * cos(x)) * (0.40 + v^2) * exp(-v^2 / 1.4)
    x, v, w = Float64[], Float64[], Float64[]
    for ix in 1:n, iv in 1:n

        ux = rand(rng)
        uv = rand(rng)
        xa = -Lx + (ix - 1 + ux) * dx0
        va = -V0 + (iv - 1 + uv) * dv0
        push!(x, xa)
        push!(v, va)
        push!(w, f0(xa, va) * dx0 * dv0)
    end
    w ./= sum(w)
    Markers(x, v, w; m = 1.2)
end

"""
    hole_case(; k = 2, seed = 505, N = 400)

The case `C(k)` of family `:C`: `N` markers of weight `1/N` on the occupied cells of the
`k × k` checkerboard of `[−1, 1]²` (the cells with `i + j` even), on the box `L = 2`, `V = 1`,
with mass `1` and the experiment `:checkerboard`. The cell `(i, j)` is
`[−1 + 2i/k, −1 + 2(i + 1)/k) × [−1 + 2j/k, −1 + 2(j + 1)/k)`, with the last cell of each axis
closed at its right or top edge.

For `k = 2` the markers are the [`checkerboard_markers`](@ref). For every other `k` they are
drawn with `MersenneTwister(seed)`: the occupied cells, `i` outer and `j` inner, hold
`N ÷ n_occ` markers each and the first `N mod n_occ` of them one more, each drawn as
`x = x_lo + (2/k) rand(rng)`, then `v = v_lo + (2/k) rand(rng)`. The represented density is
`1/A` on the occupied cells, `A = (2/k)² n_occ` their area, and zero elsewhere; the positive
density is `1/4` on `[−1, 1]²`; both vanish outside `[−1, 1]²`. For `k = 2` both equal
[`checkerboard_density`](@ref).

An `ArgumentError` refuses `k < 1` or `N < 1`, an odd `N` at `k = 2`, and a negative seed.
"""
function hole_case(; k::Integer = 2, seed::Integer = 505, N::Integer = 400)
    k >= 1 && N >= 1 || throw(ArgumentError(
        "hole_case needs k >= 1 and N >= 1, got k = $k, N = $N"))
    k == 2 && isodd(N) &&
        throw(ArgumentError("hole_case needs an even N at k = 2, got N = $N"))
    seed >= 0 || throw(ArgumentError("hole_case needs a seed >= 0, got seed = $seed"))
    occupied = [(i, j) for i in 0:(k - 1) for j in 0:(k - 1) if iseven(i + j)]
    markers = k == 2 ? checkerboard_markers(; seed, N) : _hole_markers(k, occupied, seed, N)
    h = 2 / k
    A = h^2 * length(occupied)
    cell(z) = min(floor(Int, (z + 1) * k / 2), k - 1)
    function represented(x, v)
        abs(x) <= 1 && abs(v) <= 1 || return 0.0
        iseven(cell(x) + cell(v)) ? 1 / A : 0.0
    end
    positive(x, v) = abs(x) <= 1 && abs(v) <= 1 ? 0.25 : 0.0
    StandardCase(
        :C, (; k, N, seed), markers, 2.0, 1.0, represented, positive, :checkerboard)
end

# N markers of weight 1/N spread over the `occupied` cells of the k × k checkerboard.
function _hole_markers(k, occupied, seed, N)
    rng = MersenneTwister(seed)
    h = 2 / k
    q, r = divrem(N, length(occupied))
    x, v = Float64[], Float64[]
    for (c, (i, j)) in enumerate(occupied), _ in 1:(q + (c <= r))

        push!(x, -1 + 2i / k + h * rand(rng))
        push!(v, -1 + 2j / k + h * rand(rng))
    end
    Markers(x, v, fill(1.0 / N, N); m = 1.0)
end

"""
    lattice_case(; N = 144, deformed = false)

The case of family `:E`, do no harm: the markers of [`lattice_example`](@ref) with
`nL = √N` and the maps `a = 0.2`, `b = 0.3` when `deformed` (`a = b = 0`, the plain lattice,
otherwise), shifted by `−π` into `[−π, π)`, with the density of the example shifted with them,
`(1 − 0.3 cos x) e^{−v²/2}/√(2π)`, in both fields, on the box `L = 2π`, `V = 4`, with mass `1`
and the experiment `:wave`. An `ArgumentError` refuses an `N` that is no square `n²` with
`n ≥ 1`.
"""
function lattice_case(; N::Integer = 144, deformed::Bool = false)
    n = _lattice_side(N)
    n === nothing &&
        throw(ArgumentError("lattice_case needs N = n^2 with n >= 1, got N = $N"))
    ex = lattice_example(; nL = n, a = deformed ? 0.2 : 0.0, b = deformed ? 0.3 : 0.0)
    markers = Markers(ex.markers.x .- π, ex.markers.v, ex.markers.w; m = ex.m)
    density(x, v) = ex.density(x + π, v)
    StandardCase(:E, (; N, deformed), markers, 2π, 4.0, density, density, :wave)
end

"""
    stress_case(name; N = 144, seed = 0, V = nothing)

The stress case `name` of family `:S1`, `:S2`, `:S3` or `:S4`, with the one density in both
fields, the experiment `:wave` and `L = 2π`:

- `:S1`, sparse tails: `f = (1 + 0.1 cos x) e^{−v²/2}/√(2π)` on the box `V` (`4` by default),
  mass `1`, `N` markers of weight `2π g(V)/N`, `g(V)` the mass of `e^{−v²/2}/√(2π)` on
  `[−V, V]` (the midpoint rule on 10 000 cells), drawn with `Xoshiro(seed)` by rejection: per
  attempt `x = −π + 2π rand(rng)`, rejected when `1.1 rand(rng) > 1 + 0.1 cos x`, then
  `v = randn(rng)`, rejected when `|v| > V`. A box with `g(V) < 1e-3` is refused before the
  draw.
- `:S2`, weights over decades: the `n × n` cell centres of `[−π, π) × [−6, 6]`, `i` outer,
  `j` inner, with the weights `f dx dv` of the bump-on-tail density
  `f = (1 + 0.04 cos x)(0.9 M(v; 0, 1) + 0.1 M(v; 4, 0.5))`, `M(v; μ, σ)` the normal density;
  `V = 6`, mass `1`.
- `:S3`, one dense cluster: `N` markers of weight `1/N` drawn with `Xoshiro(seed)`, per attempt
  `x = periodic_wrap(σx randn(rng), π)`, then `v = σv randn(rng)`, rejected when `|v| > V`, with
  `σx = 0.1 L`, `σv = 0.6`, `V = 3`, mass `1`; the density is the Gaussian with its two periodic
  images in `x`, `Σ_{k=−1}^{1} e^{−(x + kL)²/(2σx²)} e^{−v²/(2σv²)} / (2π σx σv)`.
- `:S4`, the clustered example of the projection: the markers of
  [`clustered_example`](@ref)`(; nx0 = n, nv0 = n, V)` shifted by `−π`, with its density
  shifted with them, on the box `V` (`3.5` by default), mass `1.5`; `clustered_example` refuses
  a `V` below the largest `|v|` of its markers.

`V` is a keyword of `:S1` and `:S4` only; `:S2` and `:S3` keep the box of their density. `:S2`
and `:S4` need `N = n²` with `n ≥ 1` and take only the seed `0`; `:S1` and `:S3` need `N ≥ 1`
and a seed `≥ 0`. An `ArgumentError` refuses every other input.
"""
function stress_case(name::Symbol; N::Integer = 144, seed::Integer = 0, V = nothing)
    name in (:S1, :S2, :S3, :S4) || throw(ArgumentError(
        "unknown stress case :$name; the names are (:S1, :S2, :S3, :S4)"))
    name in (:S2, :S3) && V !== nothing &&
        throw(ArgumentError(
            "stress_case :$name keeps the box of its density; V is a keyword of :S1 and :S4 only"))
    if name in (:S2, :S4)
        n = _lattice_side(N)
        n === nothing && throw(ArgumentError(
            "stress_case :$name needs N = n^2 with n >= 1, got N = $N"))
        seed == 0 || throw(ArgumentError(
            "stress_case :$name is deterministic; its only seed is 0, got seed = $seed"))
    else
        N >= 1 || throw(ArgumentError("stress_case :$name needs N >= 1, got N = $N"))
        seed >= 0 ||
            throw(ArgumentError("stress_case :$name needs a seed >= 0, got seed = $seed"))
    end
    if name in (:S1, :S4)
        box = V === nothing ? (name === :S1 ? 4.0 : 3.5) : Float64(V)
        _check_box(2π, box)  # before S1 draws or S4 builds its markers
    end
    if name === :S1
        _stress_s1(N, seed, box)
    elseif name === :S2
        _stress_s2(n)
    elseif name === :S3
        _stress_s3(N, seed)
    else
        _stress_s4(n, box)
    end
end

function _stress_s1(N, seed, V)
    nq = 10_000
    g = sum(exp(-(-V + (j - 0.5) * 2V / nq)^2 / 2) for j in 1:nq) / sqrt(2π) * 2V / nq
    # a draw is accepted with probability g/1.1, so at the floor a marker takes about 1100 attempts on average
    g >= 1e-3 || throw(ArgumentError(
        "stress_case :S1 needs g(V) >= 0.001, the mass of e^(-v^2/2)/sqrt(2π) on [-V, V], " *
        "got g(V) = $g at V = $V"))
    rng = Xoshiro(seed)
    x, v = Float64[], Float64[]
    while length(x) < N
        xa = -π + 2π * rand(rng)
        1.1 * rand(rng) > 1 + 0.1 * cos(xa) && continue
        va = randn(rng)
        abs(va) > V && continue
        push!(x, xa)
        push!(v, va)
    end
    density(x, v) = (1 + 0.1 * cos(x)) * exp(-v^2 / 2) / sqrt(2π)
    StandardCase(:S1, (; N, seed, V), Markers(x, v, fill(2π * g / N, N); m = 1.0), 2π, V,
        density, density, :wave)
end

function _stress_s2(n)
    L, V = 2π, 6.0
    dx, dv = L / n, 2V / n
    M(v, μ, σ) = exp(-(v - μ)^2 / (2σ^2)) / (σ * sqrt(2π))
    density(x, v) = (1 + 0.04 * cos(x)) * (0.9 * M(v, 0.0, 1.0) + 0.1 * M(v, 4.0, 0.5))
    x, v = _cell_centres(L, V, n)
    StandardCase(:S2, (; N = n^2, seed = 0, V), Markers(x, v, dx * dv .* density.(x, v)),
        L, V, density, density, :wave)
end

function _stress_s3(N, seed)
    L, V = 2π, 3.0
    σx, σv = 0.1 * L, 0.1 * 2V
    rng = Xoshiro(seed)
    x, v = Float64[], Float64[]
    while length(x) < N
        xa = periodic_wrap(σx * randn(rng), π)
        va = σv * randn(rng)
        abs(va) > V && continue
        push!(x, xa)
        push!(v, va)
    end
    density(x, v) = sum(exp(-(x + k * L)^2 / (2σx^2)) for k in -1:1) *
                    exp(-v^2 / (2σv^2)) / (2π * σx * σv)
    StandardCase(:S3, (; N, seed, V), Markers(x, v, fill(1.0 / N, N)), L, V, density,
        density, :wave)
end

function _stress_s4(n, V)
    ex = clustered_example(; nx0 = n, nv0 = n, V)
    markers = Markers(ex.markers.x .- π, ex.markers.v, ex.markers.w; m = ex.m)
    density(x, v) = ex.density(x + π, v)
    StandardCase(:S4, (; N = n^2, seed = 0, V), markers, 2π, V, density, density, :wave)
end

export StandardCase, STANDARD_LADDERS, wave_case, hole_case, lattice_case, stress_case,
       case_density
