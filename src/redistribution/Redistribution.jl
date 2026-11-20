"""
    Redistribution

Symplectic redistribution in phase space: a positive, moment-compatible B-spline projection and
translation- and permutation-invariant pairwise spline Hamiltonians whose coefficients solve a
linear inverse problem.
"""
module Redistribution

using ..ParticleRemapping: AbstractRemap, Markers, RemapResult, RemapStage, canonical_J,
                           interpolate_stages, pack
import ..ParticleRemapping: remap, remap_states, remap_frames
using ForwardDiff: ForwardDiff
using LinearAlgebra: Diagonal, I, dot, norm, opnorm
using SimpleSplines: BSplineBasis, Free, Periodic, Spline, TensorProductQuadrature,
                     UniformMesh,
                     bases, basis_index, evaluate, evaluate_all!, local_width,
                     mass_operator,
                     nbasis, quadrature_nodes, quadrature_weights, ⊗

include("canonical.jl")
include("projection.jl")
include("relaxation.jl")
include("hamiltonians.jl")
include("solver.jl")
include("closedloop.jl")
include("method.jl")
include("frames.jl")
include("plotting.jl")

# canonical structure
export kinematic_projection, kinematic_pullback, euler_scaling_fields
export velocity_fold_momentum, velocity_fold_determinant, isotropic_maxwellian
export maxwellian_fold_momentum_bound
# spline projection and targets
export PhaseSpaceSplines, deposit_phase_space, l2_spline_projection,
       lumped_spline_projection
export spline_moments, evaluate_phase_space_spline, phase_space_spline_grid
export number_density_volumes, spline_evaluation_matrix, kkt_moment_matrix
export kkt_spline_projection, kkt_redistribution_weights, recover_velocity_branch
export declustering_energy, declustering_relaxation, kernel_density_cv
# ensemble Hamiltonians
export PairwiseSplineHamiltonian, FourierMomentHamiltonian, CentroidFrameHamiltonian
export hamiltonian_parameter_count, basis_hamiltonians, basis_hamiltonian_gradients
export pairwise_hamiltonian_derivatives, ensemble_hamiltonian, ensemble_hamiltonian_gradient
export collective_moments, fourier_pairwise_hamiltonian
export uncoupled_control_step, mean_force_control_step
# parameter solve and forward step
export RedistributionScheme, ImplicitMidpointStep, SymplecticEulerStep
export redistribution_stage_point, redistribution_parameter_matrix, redistribution_scaling
export solve_redistribution_parameters, redistribution_contraction_number
export forward_redistribution_step, redistribution_substeps
export pseudotime_path
# closed loop
export closed_loop_defect, implicit_step_jacobian, closed_loop_example
export target_matching_parameters, target_matching_closed_loop
export cotangent_lift_modes, variational_closed_loop, fixed_grid_closed_loop
# the method
export redistribution_example, SplineRedistribution
# plotting
export plot_velocity_momentum_fold, plot_spline_projections, plot_target_relaxation
export plot_newton_convergence, plot_closed_loop_defect, animate_target_relaxation

end
