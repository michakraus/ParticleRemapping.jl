"""
    Projection

Hamiltonian particle projection: remap Hamiltonians built from kicks, drifts and non-separable
implicit-midpoint stages with a centroid shift, fitted to target points from a projected density;
charge-conserving reweighting and target weights from a dynamical low-rank density.
"""
module Projection

using ..ParticleRemapping: ParticleRemapping, AbstractRemap, Markers, RemapResult,
                           RemapStage, UniformBasis, interpolate_stages, remap_states,
                           canonical_J, eachnonzero, momentum, unpack, velocity
import ..ParticleRemapping: remap

using ForwardDiff: ForwardDiff
using GeometricSolvers: DogLeg, GaussNewton, LeastSquaresProblem, solve
using LinearAlgebra: Diagonal, I, Symmetric, det, dot, isposdef, norm, pinv, qr, svd, tr
using Printf: @printf
using Random: Random, AbstractRNG, Xoshiro

include("structure.jl")
include("stages.jl")
include("targets.jl")
include("parametrisation.jl")
include("fit.jl")
include("charge.jl")
include("lowrank.jl")
include("diagnostics.jl")
include("remap.jl")
include("experiments/hole_filling.jl")
include("experiments/dlra_weights.jl")
include("frames.jl")
include("plotting.jl")

# structural facts
export mean_field_example, mean_field_hamiltonian, mean_field_single_hamiltonian,
       ensemble_gradient, ensemble_symplectic_euler, ensemble_implicit_midpoint,
       particle_symplectic_euler, particle_implicit_midpoint, shoelace_area, inside_polygon,
       mean_field_step_example
export canonical_bracket, effective_divergence, frozen_weight_divergence,
       compressibility_bracket
export implicit_weight_map, fold_determinant, maxwellian_fold_peak, maxwellian_fold_roots

# stages
export centroid_kick, naive_kick, velocity_drift, nonseparable_stage,
       nonseparable_hamiltonian,
       hamiltonian_one_step, marker_midpoint_step, nonseparable_jacobian_determinant,
       nonseparable_divergence, kick_drift_jacobian_determinants, pm_v_coupling

# projection and targets
export HatProjection, lumped_projection, hat_interpolate, moser_laplacian, moser_targets,
       projection_targets, projection_target_weights

# parametrisation and fit
export AbstractRemapParametrisation, KickDriftParametrisation, NonSeparableParametrisation,
       remap_parameter_count, kick_drift_map, nonseparable_map, remap_stages,
       nonseparable_sensitivity_error, fit_remap_parameters

# charge-conserving reweighting
export spline_grid_charge, spline_charge_matrix, charge_reweight,
       additive_charge_correction,
       ChargeConservingReweighting

# dynamical low-rank approximation
export VlasovGrid, LowRankDensity, lowrank_rank, lowrank_full, lowrank_mass, numerical_rank,
       lowrank_factorise, bug_step, dlra_solve, vlasov_full_solve, lowrank_evaluate,
       vlasov_characteristics

# diagnostics
export EllipticHole, inside_hole, hole_outline, hole_fill, largest_empty_disk,
       local_empty_disk,
       density_discrepancy, weight_errors_by_band

# the full algorithm and the experiments
export KickDriftRemap, NonSeparableRemap, remap_parametrisation
export hole_filling_problem, hole_remap_method, hole_filling_experiment,
       projected_weight_errors, dlra_weights_experiment
export experiment_remap, moser_flow_history, dlra_push_frames

# plotting
export plot_mean_field_step, plot_implicit_weight_fold, plot_pm_v_coupling,
       plot_hole_filling,
       plot_charge_correction, plot_dlra_weights, plot_lumped_projection, plot_target_flow,
       plot_hole_outline!, plot_hole_stages, animate_dlra

end
