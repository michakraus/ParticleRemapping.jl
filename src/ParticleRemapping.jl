"""
    ParticleRemapping

Symplectic remapping, projection and redistribution of weighted particle-in-cell markers.

A marker carries a position `x`, a velocity `v` and a weight `w`; its canonical momentum is
`p = m w v`. The package collects the constructions of four notes on maps that move, reweight or
redistribute such markers while they keep the canonical form `Σ dx ∧ dp`:

- [`Remapping`](@ref): the factorisation of every symplectic redistribution into a canonical map
  and a gauge step, velocity rezoning by the gauge step under charge and energy constraints, and
  canonical transport with current-quantile matching (1D-1V).
- [`Projection`](@ref): kick–drift and non-separable remap Hamiltonians with a centroid shift that
  conserve the total momentum, fitted to target points from a projected density; charge-conserving
  reweighting and target weights from a dynamical low-rank density.
- [`Redistribution`](@ref): a positive, moment-compatible B-spline projection and translation- and
  permutation-invariant pairwise spline Hamiltonians, solved for their coefficients through a
  linear inverse problem.
- [`SlicedRedistribution`](@ref): sliced quantile optimal-transport cotangent lifts and the
  canonical weight squeeze.

The common layer — [`Markers`](@ref), the canonical structure, grids and bases, the
[`remap`](@ref) interface — is shared by all four. Plotting lives in the Makie extension.
"""
module ParticleRemapping

using ForwardDiff: ForwardDiff
using LinearAlgebra: LinearAlgebra, I, norm, opnorm

# the common layer
include("common/markers.jl")
include("common/symplectic.jl")
include("common/grids.jl")
include("common/bsplines.jl")
include("common/interface.jl")
include("common/stages.jl")
include("common/svg.jl")
include("common/plotting.jl")

# the four constructions, each in its own module
include("remapping/Remapping.jl")
include("projection/Projection.jl")
include("redistribution/Redistribution.jl")
include("sliced/SlicedRedistribution.jl")

using .Remapping
using .Projection
using .Redistribution
using .SlicedRedistribution

# Re-export what each construction exports, so that `using ParticleRemapping` gives the whole API.
for mod in (Remapping, Projection, Redistribution, SlicedRedistribution)
    for name in names(mod)
        name === nameof(mod) || @eval export $name
    end
end

export Remapping, Projection, Redistribution, SlicedRedistribution

# the standard experiments, which every construction runs
include("standard.jl")
# the examples of the projection chapter, and the graphical introductions of the overview pages
include("projection_examples.jl")
include("algorithm_figures.jl")

end
