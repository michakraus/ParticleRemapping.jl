module ParticleRemappingVlasovMethodsExt

# The downstream family P of the standard test set on VlasovMethods' particle model: the case P2
# of VlasovMethods' own bump-on-tail draw, and the runner of a case of family P with a remap every
# R steps. `full_leapfrog_step!`, `energy` and `update!` are internal names of VlasovMethods;
# `Potential` and `PeriodicBasisSpline` are PoissonSolvers' names, which VlasovMethods uses.

using LinearAlgebra: norm
using ParticleRemapping: ParticleRemapping, AbstractRemap, Markers, PeriodicGrid,
                         StandardCase, case_charge_noise,
                         case_conservation, case_outside, case_weights, grid_moment,
                         periodic_wrap, pic_method, remap, _REMAP_FAILURES
using Random: Random
using VlasovMethods: VlasovMethods, BumpOnTail, ParticleDistribution, PeriodicBasisSpline,
                     PoissonField, Potential, energy, full_leapfrog_step!, initialize!,
                     update!

include("bump_on_tail.jl")
include("run.jl")

end
