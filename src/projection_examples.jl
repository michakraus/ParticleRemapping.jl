# The two examples of the projection chapter: a mildly deformed lattice with its exact density,
# and the clustered markers of the spline chapter with the exact inverse of their clustering as
# the density. The second needs `redistribution_example` and serves the projection, so it lives
# here, after the four modules.

"""
    lattice_example(; nL = 24, a = 0.2, b = 0.3) -> (; initial, markers, density, L, V, m, N)

The main example of the projection: the lattice of the hole problem without a hole, its cells
moved by a mild smooth map. `initial` holds the `nL × nL` cell centres `(ξ, η)` of
`[0, L) × [−V, V]`, in the order of [`hole_filling_problem`](@ref), with the weights `c f(ξ, η)`,
`c = 2LV / nL²` and `f(x, v) = (1 + 0.3 cos x) e^{−v²/2} / √(2π)`. `markers` holds these cells
moved by `x = (ξ + a sin(ξ − 0.5)) mod L`, `v = η + b sin(π (η + V) / 2V)`, with the weights
`c f(x, v) J`, where
`J = (1 + a cos(ξ − 0.5)) (1 + (π b / 2V) cos(π (η + V) / 2V))` is the Jacobian of the map at the
cell centre, so that `c J` is the area of the moved cell to second order. The scalar
`density(x, v) = f(x, v)` is the exact density of both marker sets. The box is `L = 2π`,
`V = 4`, the mass `m = 1` and `N = nL²`. With `a = b = 0`, `markers` is the plain lattice.

The map of `x` is monotone for `|a| < 1`, the map of `v` for `|b| < 2V/π = 8/π`; it keeps
`±V` fixed, so every marker stays inside `(−V, V)`. An `ArgumentError` refuses `nL < 1` and an
`a` or a `b` outside these bounds, `NaN` included.
"""
function lattice_example(; nL::Integer = 24, a::Real = 0.2, b::Real = 0.3)
    nL >= 1 || throw(ArgumentError("lattice_example needs nL ≥ 1, got nL = $nL"))
    abs(a) < 1 || throw(ArgumentError(
        "lattice_example needs |a| < 1 for a monotone map in x, got a = $a"))
    abs(b) < 8 / π || throw(ArgumentError(
        "lattice_example needs |b| < 8/π for a monotone map in v, got b = $b"))
    L, V, m = 2π, 4.0, 1.0
    density(x, v) = (1 + 0.3 * cos(x)) * exp(-v^2 / 2) / sqrt(2π)
    dx, dv = L / nL, 2V / nL
    c = L * 2V / Int(nL)^2  # Int: nL^2 overflows in a narrow integer type
    ξ = [(i - 0.5) * dx for i in 1:nL for j in 1:nL]
    η = [-V + (j - 0.5) * dv for i in 1:nL for j in 1:nL]
    x = mod.(ξ .+ a .* sin.(ξ .- 0.5), L)
    v = η .+ b .* sin.(π .* (η .+ V) ./ (2V))
    J = (1 .+ a .* cos.(ξ .- 0.5)) .* (1 .+ (π * b / (2V)) .* cos.(π .* (η .+ V) ./ (2V)))
    initial = Markers(ξ, η, c .* density.(ξ, η); m)
    markers = Markers(x, v, c .* density.(x, v) .* J; m)
    (; initial, markers, density, L, V, m, N = length(markers))
end

"""
    clustered_example(; nx0 = 24, nv0 = 24, V = 3.5, Lx = 2π, vlims = (-2.4, 2.4), m = 1.5,
        amplitude = 0.45, dt = 0.05, nsteps = 60,
        force = x -> 1.4 sin(x − 0.5) + 0.8 sin(2x)) -> (; initial, markers, density, L, V, m, N)

The second example of the projection: the clustered markers of the spline chapter,
[`redistribution_example`](@ref) with the same keywords, in the velocity box `[−V, V]` of the
projection. `initial` and `markers` are the markers before and after the `nsteps` leapfrog steps,
`L = float(Lx)`, `m` the mass and `N` the number of markers.

The scalar `density(x, v)` is the density that the clustered markers carry,
`f₀(Φ⁻¹(x, v))` with `f₀ = (1 + amplitude cos x) e^{−v²/2} / (Lx √(2π))` and `Φ⁻¹` the exact
inverse of the leapfrog steps, taken in reverse order:
`v −= ½ dt F(x)/m; x = (x − dt v) mod Lx; v −= ½ dt F(x)/m`. It is `f₀` continued beyond the
lattice, so it does not vanish outside the velocity range of the lattice. `density(x, v) Δx Δv`
is the weight of the marker at `(x, v)`, with the spacings `Δx = Lx / nx0` and
`Δv = (vlims[2] − vlims[1]) / nv0` of the lattice.

An `ArgumentError` refuses `nx0 < 1` or `nv0 < 1`, and a box `V` below the largest `|v|` of the
clustered markers.
"""
function clustered_example(; nx0::Integer = 24, nv0::Integer = 24, V::Real = 3.5,
        Lx::Real = 2π, vlims = (-2.4, 2.4), m::Real = 1.5, amplitude::Real = 0.45,
        dt::Real = 0.05, nsteps::Integer = 60,
        force = x -> 1.4 * sin(x - 0.5) + 0.8 * sin(2x))
    nx0 >= 1 && nv0 >= 1 || throw(ArgumentError(
        "clustered_example needs nx0 ≥ 1 and nv0 ≥ 1, got nx0 = $nx0, nv0 = $nv0"))
    ex = redistribution_example(;
        nx0, nv0, Lx, vlims, m, amplitude, dt, nsteps, force)
    vmax = maximum(abs, ex.markers.v)
    vmax <= V || throw(ArgumentError("clustered_example needs V ≥ $vmax, the largest |v| " *
                        "of the clustered markers, got V = $V"))
    function density(x, v)
        for _ in 1:nsteps
            v -= 0.5 * dt * force(x) / m
            x = mod(x - dt * v, Lx)
            v -= 0.5 * dt * force(x) / m
        end
        (1 + amplitude * cos(x)) * exp(-0.5 * v^2) / (Lx * sqrt(2π))
    end
    (; ex.initial, ex.markers, density, L = ex.Lx, V = float(V), m = ex.markers.m,
        N = length(ex.markers))
end

export lattice_example, clustered_example
