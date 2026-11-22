# [Library: SlicedRedistribution](@id sliced-library)

```@meta
CurrentModule = ParticleRemapping.SlicedRedistribution
```

```@docs
SlicedRedistribution
```

## Slices and targets

```@docs
uniform_quantile_targets
sorted_slices
slice_anchors
```

## The lifts

```@docs
sliced_quantile_potential
sliced_quantile_hessian
split_quantile_ot_lift
dual_quantile_ot_lift
per_axis_sliced_lift
```

## The weight squeeze and the weights

```@docs
canonical_weight_squeeze
kernel_density_ratio
squeeze_multipliers
```

## The composed map

```@docs
SlicedLiftPlan
sliced_lift_plan
apply_sliced_lift
SlicedQuantileRemap
remap(::SlicedQuantileRemap, ::Markers{T, <:AbstractVector}) where {T}
remap(::SlicedQuantileRemap, ::Markers{T, <:AbstractMatrix}) where {T}
```

## Stages and frames

```@docs
remap_states(::SlicedQuantileRemap, ::Markers, ::RemapResult)
ParticleRemapping.remap_frames(::SlicedQuantileRemap, ::Markers, ::RemapResult)
```

## Diagnostics

```@docs
periodic_wrap
phase_quadrant_counts
phase_cell_statistics
cell_count_history
smoothed_density_difference
```

## Experiments

```@docs
factorisation_pullback_experiment
pseudo_time_flow_experiment
autonomous_map_experiment
composition_defect_experiment
checkerboard_markers
checkerboard_experiment
electric_wave_markers
electric_wave_simulation
long_time_remap_experiment
sliced_stages_example
```

## Probes of constructions that do not work

```@docs
repulsion_wall_probe
rotated_frame_lift_probe
checkerboard_reweighting_control
```

## Plotting

These functions need a Makie backend, for example `using CairoMakie`.

```@docs
plot_sliced_stages
plot_checkerboard_filling
plot_cell_count_history
plot_wave_remap
animate_sliced_stages
animate_long_time_run
```
