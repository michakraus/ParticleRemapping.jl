# Plotting functions of the Hamiltonian particle projection. The methods are in the Makie extension,
# which loads with any Makie backend (`using CairoMakie`).

"""
    plot_mean_field_step(data = mean_field_step_example())

Two panels: the markers and the hole (shaded) of [`mean_field_step_example`](@ref) before one
implicit-midpoint step of the mean-field Hamiltonian, and the markers and the image of the hole
under the single-particle map after the step; the image contains no marker. The markers of this
example carry no weights (the Hamiltonian does not depend on them), so they are drawn in one
colour. Needs a Makie backend.
"""
function plot_mean_field_step end

"""
    plot_implicit_weight_fold(; vmax = 4, levels = (0.5, 1.5))

The fold of the implicit weight for a Maxwellian: `p / p_max = |v| f_M(v) / (v_th f_M(v_th))`
against `|v| / v_th`, with the levels `p = levels .* p_max` and the number of roots of each
([`maxwellian_fold_roots`](@ref)); the Jacobian changes sign at `|v| = v_th` (dotted). Needs a
Makie backend.
"""
function plot_implicit_weight_fold end

"""
    plot_pm_v_coupling(data = pm_v_coupling())

`det DN − 1` of one marker step at `x = π` against `v` for the separable and the non-separable
`φ` of [`pm_v_coupling`](@ref): opposite signs at `±v` for the separable one, the same sign for the
non-separable one. Needs a Makie backend.
"""
function plot_pm_v_coupling end

"""
    plot_hole_filling(experiment; window = 1.2)

The markers near the hole (dashed outline) of a [`hole_filling_experiment`](@ref) run with the
variants `:weight_dependent`, `:non_separable` and `:control_C`, coloured and sized by weight, in
`(x, v)` with `v = p/(m w)` for the weights of each panel: (a) before, (b) the targets, (c) the
reweighting alone, (d) control C, (e) the kicks and drifts, (f) the non-separable map, (g) and (h)
the maps of (e) and (f) after the charge-conserving reweighting ([`plot_hole_stages`](@ref)).
Needs a Makie backend.
"""
function plot_hole_filling end

"""
    plot_charge_correction(experiments...; variant = :non_separable, labels)

The charge-conserving reweighting after the map `variant` of each [`hole_filling_experiment`](@ref)
run: (a) the grid charge of the target weights (lines) and of the corrected weights (dots) at the
positions after the map, against the grid charge before the map, relative to its maximum; (b) the
correction `ln(w′/w*)` per marker. The colour tells the runs apart, so the markers of (b) are not
coloured by weight. Needs a Makie backend.
"""
function plot_charge_correction end

"""
    plot_dlra_weights(experiment; tolerance = 1e-6, error_tolerance = 1e-4)

Target weights after a push, from a [`dlra_weights_experiment`](@ref) run: (a) the pushed markers
and the image of the hole; (b) the non-separable map with the weights of `f_r` at `tolerance`,
after the charge-conserving reweighting, both coloured and sized by weight on one scale; (c)
`f_r − f(T)` at `error_tolerance`; (d) the numerical rank of `f_r` over time; (e) the relative
error of the target weights against `|v*|` for the tolerances and the lumped projection, coloured
by the source of the weights. Needs a Makie backend.
"""
function plot_dlra_weights end

"""
    plot_lumped_projection(problem = hole_filling_problem(); projection = HatProjection(32, 33))

The first step of the algorithm: (a) the markers of `problem`, coloured and sized by weight,
(b) the lumped projection `f_h` of the represented density and (c) the lumped projection `g_h` of
the marker density ([`lumped_projection`](@ref)), in `(x, v)`, with the hole outline where
`problem` has a `hole` field. Needs a Makie backend.
"""
function plot_lumped_projection end

"""
    plot_target_flow(problem, targets; window = nothing)

The second step of the algorithm: the markers of `problem`, coloured and sized by weight, their
target points from [`projection_targets`](@ref) (open circles) and the displacement of each marker
(segments), with the hole outline where `problem` has a `hole` field; `window = (x₀, x₁, v₀, v₁)`
restricts the view. Needs a Makie backend.
"""
function plot_target_flow end

"""
    plot_hole_outline!(ax, hole::EllipticHole)

Draws the outline of `hole` into the axis `ax`, dashed, in the outline colour of the package; the
`decorate!` of [`plot_stages`](@ref ParticleRemapping.plot_stages) and
[`animate_stages`](@ref ParticleRemapping.animate_stages) for the remaps of this chapter. Needs a
Makie backend.
"""
function plot_hole_outline! end

"""
    plot_hole_stages(stages, hole::EllipticHole; window = 1.2, columns = 4, kwargs...)

[`plot_stages`](@ref ParticleRemapping.plot_stages) of the markers near `hole`: the view is the
square of half-width `window` around the centre of the hole, with the hole outline in every panel.
Only the markers that lie near the view in some stage are drawn, so that the marker sizes suit the
markers in view. The keywords `kwargs` go to `plot_stages`. Needs a Makie backend.
"""
function plot_hole_stages end

"""
    animate_dlra(path, frames, grid; framerate = 20)

Writes the animation of [`dlra_push_frames`](@ref) to `path`: the low-rank density on `grid` as a
heat map, the markers and the image of the hole along the push. The point of the animation is the
density, so the markers are drawn in one light colour and the image of the hole in vermillion,
both visible over the heat map. Needs a Makie backend.
"""
function animate_dlra end
