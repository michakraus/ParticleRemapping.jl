# [Library: Redistribution](@id rd-library)

```@meta
CurrentModule = ParticleRemapping.Redistribution
```

```@docs
Redistribution
```

## The method

```@docs
SplineRedistribution
remap(::SplineRedistribution, ::Markers)
remap_states(::SplineRedistribution, ::Markers, ::RemapResult)
remap_frames(::SplineRedistribution, ::Markers, ::RemapResult)
redistribution_example
```

## Canonical structure

```@docs
kinematic_projection
kinematic_pullback
euler_scaling_fields
velocity_fold_momentum
velocity_fold_determinant
isotropic_maxwellian
maxwellian_fold_momentum_bound
```

## Spline projection and targets

```@docs
PhaseSpaceSplines
deposit_phase_space
l2_spline_projection
lumped_spline_projection
spline_moments
evaluate_phase_space_spline
phase_space_spline_grid
number_density_volumes
spline_evaluation_matrix
kkt_moment_matrix
kkt_spline_projection
kkt_redistribution_weights
recover_velocity_branch
declustering_energy
declustering_relaxation
kernel_density_cv
```

## Ensemble Hamiltonians

```@docs
PairwiseSplineHamiltonian
FourierMomentHamiltonian
CentroidFrameHamiltonian
hamiltonian_parameter_count
basis_hamiltonians
basis_hamiltonian_gradients
pairwise_hamiltonian_derivatives
ensemble_hamiltonian
ensemble_hamiltonian_gradient
collective_moments
fourier_pairwise_hamiltonian
uncoupled_control_step
mean_force_control_step
```

## Parameter solve and forward step

```@docs
RedistributionScheme
ImplicitMidpointStep
SymplecticEulerStep
redistribution_stage_point
redistribution_parameter_matrix
redistribution_scaling
solve_redistribution_parameters
redistribution_contraction_number
forward_redistribution_step
redistribution_substeps
pseudotime_path
```

## Closed loop

```@docs
closed_loop_defect
implicit_step_jacobian
closed_loop_example
target_matching_parameters
target_matching_closed_loop
cotangent_lift_modes
variational_closed_loop
fixed_grid_closed_loop
```

## Plotting

```@docs
plot_velocity_momentum_fold
plot_spline_projections
plot_target_relaxation
plot_newton_convergence
plot_closed_loop_defect
animate_target_relaxation
```
