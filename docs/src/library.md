# Common library

```@meta
CurrentModule = ParticleRemapping
```

```@docs
ParticleRemapping
```

The common layer holds what the four constructions share: the marker sets, the canonical
structure, the grids and bases, the remap interface and the plots of marker sets. The names of
each construction are on its own library page; `using ParticleRemapping` makes all of them
available.

## Markers

```@docs
Markers
spatial_dimension
momentum
velocity
canonical
from_canonical
gauge_step
total_weight
total_momentum
kinetic_energy
```

## Canonical structure

```@docs
pack
unpack
canonical_J
symplecticity_defect
fd_jacobian
presymplectic_matrix
pullback_residual
```

## Grids and bases

```@docs
PeriodicGrid
cellof
shape
shape_matrix
grid_moment
bspline3
dbspline3
UniformBasis
eachnonzero
basis_values
```

## The remap interface

```@docs
AbstractRemap
RemapResult
remap
conservation_report
```

## The states of a remap

Every construction gives the states of its remap, from the markers before it to the markers after
it, as a vector of [`RemapStage`](@ref): [`remap_states`](@ref) for a figure with one panel per
stage, and [`remap_frames`](@ref) for an animation.

```@docs
RemapStage
remap_states
remap_frames
interpolate_stages
```

## Plots

The plotting functions need a Makie backend, for example `using CairoMakie`, which loads the
extension `ParticleRemappingMakieExt`. Every figure and animation of the manual uses the style of
[`remap_theme`](@ref); the [figures and animations](@ref figures) page states its rules.

```@docs
remap_theme
plot_markers
plot_markers!
plot_remap
plot_stages
animate_stages
```

## SVG pictures of the steps

The overview page of each construction introduces its algorithm with SVG pictures, which need no
plotting package: one panel per step and an animation through the steps (SMIL, no script).

```@docs
SVGFigure
svg_stage_panels
svg_stage_animation
svg_line
svg_rect
svg_text
svg_arrow
svg_category_colors
algorithm_svg
```

## Standard experiments

The [standard experiments](@ref standard-experiments), which every construction runs.

```@docs
STANDARD_METHODS
standard_method
ShiftedRemap
checkerboard_density
electric_wave_density
standard_wave_markers
standard_checkerboard
standard_long_time
standard_remap_history
```
