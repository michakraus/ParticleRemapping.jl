"""
    Remapping

Hamiltonian remapping of weighted markers in 1D-1V: the factorisation of every symplectic
redistribution into a canonical map and a gauge step, velocity rezoning by the gauge step, and
canonical transport with current-quantile matching.
"""
module Remapping

using ..ParticleRemapping: ParticleRemapping, AbstractRemap, Markers, PeriodicGrid,
                           RemapResult, RemapStage,
                           cellof, fd_jacobian, from_canonical, grid_moment,
                           interpolate_stages, momentum, pullback_residual, shape
import ..ParticleRemapping: remap, remap_frames, remap_states
using LinearAlgebra: Diagonal
using Random: Xoshiro
using SparseArrays: dropzeros!, sparse, spdiagm, spzeros

include("gauge.jl")
include("gauge_examples.jl")
include("transport.jl")
include("transport_examples.jl")
include("experiments.jl")
include("canonicity.jl")
include("frames.jl")
include("plotting.jl")

# Construction 1: velocity rezoning by the gauge step
export GaugeRezoning, GaugeRezoningProblem
export velocity_groups, gauge_targets, speed_gap_ratio, mean_gap_ratio
export gauge_constraints, gauge_objective, gauge_objective_gradient, gauge_stationarity
export gauge_hessian, gauge_kkt_matrix, gauge_schur_matrix, gauge_newton_step, gauge_merit
export gauge_newton, gauge_homotopy
export lambert_w0, charge_constraint_example, charge_stretches, charge_constraint_residual
export charge_multiplier_limit, charge_stationary_point, charge_fold
export eight_marker_example, eight_marker_homotopy

# Construction 2: canonical transport with consistent weights
export CurrentQuantileMatching, SeparableMaxwellianTarget
export simpson_rule, represented_density, cumulative_current_error
export cotangent_lift_flow, cotangent_lift, lift_compatibility_residual
export maxwellian_current, maxwellian_current_inverse
export cell_cumulative_current, cell_current_total, cell_current_inverse
export branch_quantile_matching, current_quantile_match, binned_weights, binned_weight_error
export cosine_profile, lift_theta, lift_theta_inverse, admissible_lift,
       admissible_lift_matching
export quantile_matching_example

# the numerical examples
export two_stream_density, gauge_rezoning_setup, quantile_matching_setup
export quantile_matching_experiment, lift_matching_frames

# the symplecticity condition
export frozen_weight_residual, symplectic_condition_checks

# plotting (Makie extension)
export plot_charge_example, animate_charge_fold, plot_eight_marker_example
export animate_eight_marker_homotopy, plot_gauge_experiment
export plot_lift_example, animate_lift_matching, plot_matching_example
export animate_matching_example, plot_kkt_structure, plot_gauge_orbits
export plot_velocity_momentum_sheet

end
