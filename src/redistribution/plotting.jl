# Plotting functions of the spline-projection redistribution. The methods are in the Makie
# extension, which loads with any Makie backend (`using CairoMakie`). The standard pictures of the
# remap are the common `plot_remap`, `plot_stages` and `animate_stages` with `remap_states` and
# `remap_frames`.

"""
    plot_velocity_momentum_fold(v, P, D; vth, pmax)

Two panels over the velocities `v`: (a) the canonical momentum `P(v) = m ΔΩ v f_h(v)`
([`velocity_fold_momentum`](@ref)) with the bound `±pmax`, and (b) the Jacobian determinant
`D(v)` ([`velocity_fold_determinant`](@ref)), which changes sign on the fold `|v| = vth`. The
suprathermal velocities `|v| > vth` are shaded. Curves, not markers: the figure has no weights.
Needs a Makie backend.
"""
function plot_velocity_momentum_fold end

"""
    plot_spline_projections(l2, kkt; markers = nothing, targets = nothing)

The `L²` projection `f_h^{L²}` and the KKT projection `f_h^KKT` of the markers as heat maps on
one symmetric signed colour scale, each a `(xs, vs, F)` tuple of
[`phase_space_spline_grid`](@ref), with the zero contour (Control D), the clustered `markers` on
the first panel and the target points `targets = (x̄, v̄)` as open circles on the second. A third
panel shows the smallest value of each projection over `x` as a function of `v`. The point of the
figure is the sign of the projected density, so the colour belongs to the heat map: the markers
are dark grey, with the area proportional to the weight. Needs a Makie backend.
"""
function plot_spline_projections end

"""
    plot_target_relaxation(history; Lx, cv = nothing)

(a) The paths of the points during the de-clustering relaxation
([`declustering_relaxation`](@ref)), from the start (dark grey dots) to the targets (open
circles), positions modulo `Lx`; with a vector `cv` of the [`kernel_density_cv`](@ref) after each
step, (b) its decrease. The relaxation moves points, not weighted markers, so the points carry no
weight colour. Needs a Makie backend.
"""
function plot_target_relaxation end

"""
    plot_newton_convergence(histories...; labels)

The scaled residual of the forward Newton solves against the Newton iterate, on a logarithmic
scale, one line per `history` of [`forward_redistribution_step`](@ref), with the first iterate of
each homotopy stage marked by a dashed line. Needs a Makie backend.
"""
function plot_newton_convergence end

"""
    plot_closed_loop_defect(lhs, rhs; variational = nothing)

Heat maps of the closed-loop defect `DΦᵀ J DΦ − J` (`lhs`), of the right-hand side
`−τ (Dgᵀ DΘ − DΘᵀ Dg)` of the closed-loop identity (`rhs`) and of their difference, from
[`closed_loop_defect`](@ref); with the defect of the variational closure (`variational`), a
further panel on the scale of the first two. Signed colour scales, symmetric about zero. Needs a
Makie backend.
"""
function plot_closed_loop_defect end

"""
    animate_target_relaxation(path, history; Lx, framerate = 20)

Record to `path` the de-clustering relaxation of the target points, one frame per step of the
`history` of [`declustering_relaxation`](@ref), positions modulo `Lx`: the start points as dark
grey dots, the moving points as open circles, the style of the targets. Needs a Makie backend.
"""
function animate_target_relaxation end
