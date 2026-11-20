# Plotting functions of the spline-projection redistribution. The methods are in the Makie
# extension, which loads with any Makie backend (`using CairoMakie`).

"""
    plot_velocity_momentum_fold(v, P, D; vth, pmax)

Two panels over the velocities `v`: the canonical momentum `P(v) = m ΔΩ v f_h(v)`
([`velocity_fold_momentum`](@ref)) with the bound `±pmax`, and the Jacobian determinant `D(v)`
([`velocity_fold_determinant`](@ref)), which changes sign on the fold `|v| = vth`. Needs a Makie backend.
"""
function plot_velocity_momentum_fold end

"""
    plot_spline_projections(l2, kkt; markers = nothing, targets = nothing)

The `L²` projection `f_h^{L²}` and the KKT projection `f_h^KKT` of the markers as heat maps, each
a `(xs, vs, F)` tuple of [`phase_space_spline_grid`](@ref), with the region `f_h^{L²} < 0`
outlined (Control D), the clustered `markers` on the first panel and the target points `targets
= (x̄, v̄)` on the second. A third panel shows the smallest value of each projection over `x` as a
function of `v`. Needs a Makie backend.
"""
function plot_spline_projections end

"""
    plot_target_relaxation(history; Lx, cv = nothing)

The paths of the points during the de-clustering relaxation ([`declustering_relaxation`](@ref)),
from the start (grey) to the targets (black), positions modulo `Lx`; with a vector `cv` of the
[`kernel_density_cv`](@ref) after each step, a second panel shows its decrease. Needs a Makie
backend.
"""
function plot_target_relaxation end

"""
    plot_redistribution_markers(before, after; targets = nothing)

The markers before and after a redistribution in `(x, v)`, coloured and sized by weight on a
shared colour scale, with the target points `targets = (x̄, v̄)` as open circles on the second
panel. Needs a Makie backend.
"""
function plot_redistribution_markers end

"""
    plot_newton_convergence(histories...; labels)

The scaled residual of the forward Newton solves against the Newton iterate, on a logarithmic
scale, one line per `history` of [`forward_redistribution_step`](@ref), with the homotopy
stages marked. Needs a Makie backend.
"""
function plot_newton_convergence end

"""
    plot_closed_loop_defect(lhs, rhs; variational = nothing)

Heat maps of the closed-loop defect `DΦᵀ J DΦ − J` (`lhs`), of the right-hand side
`−τ (Dgᵀ DΘ − DΘᵀ Dg)` of the closed-loop identity (`rhs`) and of their difference, from
[`closed_loop_defect`](@ref); with the defect of the variational closure (`variational`), a
fourth panel on the same colour scale. Needs a Makie backend.
"""
function plot_closed_loop_defect end

"""
    animate_redistribution(path, states, targets; p_scale, framerate = 12)

Record to `path` the symplectic redistribution step in pseudo-time: the markers move along the
packed states `states` of [`pseudotime_path`](@ref) in the plane `(q, p/p_s)` towards the
packed `targets` (open circles). Needs a Makie backend.
"""
function animate_redistribution end

"""
    animate_target_relaxation(path, history; Lx, framerate = 6)

Record to `path` the de-clustering relaxation of the target points, one frame per step of the
`history` of [`declustering_relaxation`](@ref), positions modulo `Lx`. Needs a Makie backend.
"""
function animate_target_relaxation end
