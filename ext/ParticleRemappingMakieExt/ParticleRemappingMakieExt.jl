module ParticleRemappingMakieExt

using Makie
using ParticleRemapping
using ParticleRemapping: Markers

include("common.jl")
include("remapping.jl")
include("projection.jl")
include("redistribution.jl")
include("sliced.jl")

end
