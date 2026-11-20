# Library: Remapping

```@meta
CurrentModule = ParticleRemapping.Remapping
```

```@docs
Remapping
```

## Velocity rezoning by the gauge step

```@docs
GaugeRezoning
remap(::GaugeRezoning, ::ParticleRemapping.Markers)
GaugeRezoningProblem
velocity_groups
gauge_targets
speed_gap_ratio
mean_gap_ratio
gauge_constraints
gauge_objective
gauge_objective_gradient
gauge_stationarity
gauge_hessian
gauge_kkt_matrix
gauge_schur_matrix
gauge_newton_step
gauge_merit
gauge_newton
gauge_homotopy
```

## Worked examples of the gauge rezoning

```@docs
lambert_w0
charge_constraint_example
charge_stretches
charge_constraint_residual
charge_multiplier_limit
charge_stationary_point
charge_fold
eight_marker_example
eight_marker_homotopy
```

## Canonical transport with consistent weights

```@docs
CurrentQuantileMatching
remap(::CurrentQuantileMatching, ::ParticleRemapping.Markers)
SeparableMaxwellianTarget
cell_cumulative_current
cell_current_total
cell_current_inverse
branch_quantile_matching
current_quantile_match
maxwellian_current
maxwellian_current_inverse
represented_density
cumulative_current_error
simpson_rule
cotangent_lift_flow
cotangent_lift
lift_compatibility_residual
binned_weights
binned_weight_error
```

## Worked examples of the canonical transport

```@docs
cosine_profile
lift_theta
lift_theta_inverse
admissible_lift
admissible_lift_matching
quantile_matching_example
```

## The numerical examples

```@docs
two_stream_density
gauge_rezoning_setup
quantile_matching_setup
quantile_matching_experiment
lift_matching_frames
```

## The symplecticity condition

```@docs
frozen_weight_residual
symplectic_condition_checks
```

## The states of a remap

```@docs
ParticleRemapping.remap_states(::GaugeRezoning, ::ParticleRemapping.Markers, ::ParticleRemapping.RemapResult)
ParticleRemapping.remap_frames(::GaugeRezoning, ::ParticleRemapping.Markers, ::ParticleRemapping.RemapResult)
ParticleRemapping.remap_states(::CurrentQuantileMatching, ::ParticleRemapping.Markers, ::ParticleRemapping.RemapResult)
ParticleRemapping.remap_frames(::CurrentQuantileMatching, ::ParticleRemapping.Markers, ::ParticleRemapping.RemapResult)
```

## Plotting

The plotting functions need a Makie backend, for example `using CairoMakie`. The standard
pictures of the two remaps are [`plot_remap`](@ref ParticleRemapping.plot_remap),
[`plot_stages`](@ref ParticleRemapping.plot_stages) and
[`animate_stages`](@ref ParticleRemapping.animate_stages) of the common layer.

```@docs
plot_charge_example
animate_charge_fold
plot_eight_marker_example
animate_eight_marker_homotopy
plot_kkt_structure
plot_gauge_experiment
plot_lift_example
animate_lift_matching
plot_matching_example
animate_matching_example
plot_gauge_orbits
plot_velocity_momentum_sheet
```
