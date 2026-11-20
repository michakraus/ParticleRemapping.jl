# [Library: Projection](@id projection-library)

```@meta
CurrentModule = ParticleRemapping.Projection
```

```@docs
Projection
```

## The method types

```@docs
KickDriftRemap
NonSeparableRemap
ChargeConservingReweighting
remap_parametrisation
```

## Structural facts

```@docs
mean_field_example
mean_field_hamiltonian
mean_field_single_hamiltonian
ensemble_gradient
ensemble_symplectic_euler
ensemble_implicit_midpoint
particle_symplectic_euler
particle_implicit_midpoint
mean_field_step_example
shoelace_area
inside_polygon
canonical_bracket
effective_divergence
frozen_weight_divergence
compressibility_bracket
implicit_weight_map
fold_determinant
maxwellian_fold_peak
maxwellian_fold_roots
```

## Stages

```@docs
centroid_kick
naive_kick
velocity_drift
nonseparable_stage
nonseparable_hamiltonian
hamiltonian_one_step
marker_midpoint_step
nonseparable_jacobian_determinant
nonseparable_divergence
kick_drift_jacobian_determinants
pm_v_coupling
```

## Projection and target points

```@docs
HatProjection
lumped_projection
hat_interpolate
moser_laplacian
moser_targets
projection_targets
projection_target_weights
```

## Parametrisation and fit

```@docs
AbstractRemapParametrisation
KickDriftParametrisation
NonSeparableParametrisation
remap_parameter_count
kick_drift_map
nonseparable_map
remap_stages
nonseparable_sensitivity_error
fit_remap_parameters
```

## Charge-conserving reweighting

```@docs
spline_grid_charge
spline_charge_matrix
charge_reweight
additive_charge_correction
```

## Dynamical low-rank approximation

```@docs
VlasovGrid
LowRankDensity
lowrank_rank
lowrank_full
lowrank_mass
numerical_rank
lowrank_factorise
bug_step
dlra_solve
vlasov_full_solve
lowrank_evaluate
vlasov_characteristics
```

## Diagnostics

```@docs
EllipticHole
inside_hole
hole_outline
hole_fill
largest_empty_disk
local_empty_disk
density_discrepancy
weight_errors_by_band
```

## Experiments

```@docs
hole_filling_problem
hole_remap_method
hole_filling_experiment
projected_weight_errors
dlra_weights_experiment
experiment_remap
moser_flow_history
dlra_push_frames
```

## Plotting

These functions need a Makie backend, for example `using CairoMakie`.

```@docs
plot_mean_field_step
plot_implicit_weight_fold
plot_pm_v_coupling
plot_hole_filling
plot_charge_correction
plot_dlra_weights
plot_lumped_projection
plot_target_flow
plot_hole_outline!
plot_hole_stages
animate_dlra
```
