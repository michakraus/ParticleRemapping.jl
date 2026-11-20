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

## Plots

The plotting functions need a Makie backend, for example `using CairoMakie`, which loads the
extension `ParticleRemappingMakieExt`.

```@docs
plot_markers
plot_markers!
plot_remap
```
