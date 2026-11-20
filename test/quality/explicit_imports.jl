using ExplicitImports
using ParticleRemapping
using Test

@test check_no_stale_explicit_imports(ParticleRemapping) === nothing
@test check_no_implicit_imports(ParticleRemapping;
    allow_unanalyzable = (ParticleRemapping.Remapping,
        ParticleRemapping.Projection, ParticleRemapping.Redistribution,
        ParticleRemapping.SlicedRedistribution)) === nothing
