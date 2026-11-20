# Plotting functions of the remapping. The methods are in the Makie extension, which loads with
# any Makie backend (`using CairoMakie`).

"""
    plot_charge_example(ex = charge_constraint_example())

The figure of the worked example with one charge constraint: (a) the speeds before, at the
targets and after the constrained remap, with the marker area proportional to the weight; (b) the
constraint `c⁽⁰⁾(A)` for the target scales `β = 1`, `β_f` and `1.6`, with `A_lim = −1/e` dashed;
(c) the curvatures `H_a/w_a` as functions of `β` up to the fold `β_f`. Needs a Makie backend.
"""
function plot_charge_example end

"""
    animate_charge_fold(path, ex = charge_constraint_example(); framerate = 30)

An animation, written to `path`, of the stationary point of the charge-constraint example as the
target scale `β` grows from `0` to the fold: the targets and the remapped markers above, the
curvatures `H_a/w_a` below. Needs a Makie backend.
"""
function animate_charge_fold end

"""
    plot_eight_marker_example(data = eight_marker_homotopy())

The figure of the eight-marker example of the gauge rezoning: (a) the markers in `(x, v)` before
(open) and after (filled), the area proportional to the weight, the targets as grey bars and the
nodes dashed; (b) the scaled constraint residual during the Newton iterations of each homotopy
step; (c) the smallest curvature `H_a/w_a` along the homotopy, with a cross at the rejected step
without homotopy. Needs a Makie backend.
"""
function plot_eight_marker_example end

"""
    animate_eight_marker_homotopy(path, data = eight_marker_homotopy(); framerate = 25)

An animation, written to `path`, of the eight markers along the homotopy from `τ = 0` to `τ = 1`,
with the targets `τ μ*` and the smallest curvature. Needs a Makie backend.
"""
function animate_eight_marker_homotopy end

"""
    plot_gauge_experiment(markers, result)

The figure of the numerical example of the gauge rezoning, for the `markers` before and the
`RemapResult` of a [`GaugeRezoning`](@ref): (a), (b) the markers in `(x, v)` before and after,
coloured by `log₁₀ w`; (c) the distribution of the stretches `λ_a`; (d) the ratio of the largest
to the mean gap per group of cell and sign, before against after, with the identity dashed. Needs
a Makie backend.
"""
function plot_gauge_experiment end

"""
    animate_gauge_homotopy(path, markers, result; framerate = 15)

An animation, written to `path`, of the accepted homotopy steps in the history of the
`RemapResult` of a [`GaugeRezoning`](@ref): the first two cells enlarged, and all cells. Needs a
Makie backend.
"""
function animate_gauge_homotopy end

"""
    plot_lift_example(setup, experiment; t = setup.t, cell = 5)

The figure of the cotangent-lift example and the numerical example of the canonical transport,
for a [`quantile_matching_setup`](@ref) and its [`quantile_matching_experiment`](@ref): (a) the
characteristics of the admissible lift `a = J^{−1/2}`; (b) at time `t` the represented density
`J(X)`, the marker density `(J(x) J(X))^{1/2}` and the weight ratio `1/φ'` (lines), with the
sampled marker counts (circles) and charges (triangles) per cell; (c) the compatibility residual
`R_c^s` per cell for the admissible and the inadmissible lift; (d) the velocity marginal in
`cell`: exact, after the matching, with the old weights kept, and after the inadmissible lift.
Needs a Makie backend.
"""
function plot_lift_example end

"""
    animate_lift_matching(path, frames, setup; framerate = 20)

An animation, written to `path`, of the [`lift_matching_frames`](@ref): the markers in `(X, V)`
coloured by `W/w`, and the marker counts and charges per cell against `J`. Needs a Makie backend.
"""
function animate_lift_matching end

"""
    plot_matching_example(ex = quantile_matching_example())

The figure of the five-marker example of the current-quantile matching: (a) the cumulative
current `F̂(V)`, each marker's closed quantile (dashed) mapped to its new velocity (filled), and
the old velocities (open); (b) the weights before (grey) and after the matching. Needs a Makie
backend.
"""
function plot_matching_example end

"""
    animate_matching_example(path, ex = quantile_matching_example(); framerate = 15)

An animation, written to `path`, of the current-quantile matching of the five-marker example,
one marker at a time. Needs a Makie backend.
"""
function animate_matching_example end

"""
    plot_kkt_structure(problem::GaugeRezoningProblem, μ, α)

The structure of the Newton step of the gauge rezoning at `(μ, α)`: the sparsity pattern of the
KKT matrix `[diag(H) Gᵀ; G 0]`, and the Schur complement `G diag(H)⁻¹ Gᵀ` with its periodic
block-tridiagonal pattern. Needs a Makie backend.
"""
function plot_kkt_structure end

"""
    plot_gauge_orbits(; m = 1, momenta = (0.25, 0.5, 1, 2), λs = (0.5, 1, 2))

The gauge step in pictures: (a) the hyperbolae `w v = p/m` in the `(v, w)` plane, along which the
gauge step moves a marker at fixed momentum, with the images of one marker under the stretches
`λs`; (b) the `(x, v)` cell of one marker, whose area scales by `λ` while its weight scales by
`1/λ`. Needs a Makie backend.
"""
function plot_gauge_orbits end

"""
    plot_velocity_momentum_sheet(; vth = 1, m = 1)

The sheet map `v ↦ p = m w(v) v` at fixed `x`, for equal weights (a straight line) and for
weights proportional to a Maxwellian of thermal speed `vth`, which folds at `|v| = vth` where
`∂_v(w v) = 0`; below, `∂_v(w v)` itself. Needs a Makie backend.
"""
function plot_velocity_momentum_sheet end
