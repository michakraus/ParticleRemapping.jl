"""
    SlicedRedistribution

Sliced symplectic redistribution: quantile optimal-transport cotangent lifts within velocity and
position slices, and the momentum-conserving canonical weight squeeze.
"""
module SlicedRedistribution

using ..ParticleRemapping: ParticleRemapping, AbstractRemap, Markers, RemapResult,
                           canonical_J,
                           momentum
import ..ParticleRemapping: remap

using ForwardDiff: ForwardDiff
using LinearAlgebra: I, Symmetric, eigmin
using Random: MersenneTwister

include("binning.jl")
include("lift.jl")
include("squeeze.jl")
include("reweighting.jl")
include("map.jl")
include("peraxis.jl")
include("diagnostics.jl")
include("experiments.jl")
include("plotting.jl")

export uniform_quantile_targets, sorted_slices, slice_anchors
export sliced_quantile_potential, sliced_quantile_hessian, split_quantile_ot_lift,
       dual_quantile_ot_lift
export canonical_weight_squeeze
export kernel_density_ratio, squeeze_multipliers
export SlicedLiftPlan, sliced_lift_plan, apply_sliced_lift, SlicedQuantileRemap
export per_axis_sliced_lift
export periodic_wrap, phase_quadrant_counts, phase_cell_statistics, cell_count_history,
       smoothed_density_difference
export factorisation_pullback_experiment, pseudo_time_flow_experiment,
       autonomous_map_experiment, composition_defect_experiment, checkerboard_experiment,
       electric_wave_markers, electric_wave_simulation, long_time_remap_experiment,
       sliced_stages_example
export repulsion_wall_probe, rotated_frame_lift_probe, checkerboard_reweighting_control
export plot_sliced_stages, plot_checkerboard_filling, plot_cell_count_history,
       plot_wave_remap, animate_sliced_stages, animate_long_time_run

end
