# The sizes of the Hamiltonian particle projection in d dimensions: the parameters of the two
# parametrisations, the nodes of the hat projection and of the Poisson problem of the
# density-interpolation flow, the shape functions of the charge-conserving reweighting, and the
# basis functions that one marker touches. The bases are those of the hole-filling experiment in
# every direction. The counts for d = 1 are checked against the code.
#
#   julia --startup-file=no --project=scripts scripts/projection/higher_dimensions.jl

using ParticleRemapping
using Printf

# the bases of the hole-filling experiment
const S = 2                                   # stages
const KD = (nξ = 32, nv = 40, nw = 16)        # kick–drift splines in ξ, v, ln w
const NS = (nξ = 16, nv = 17, nℓ = 6)         # non-separable splines in ξ, v, ln w
const HATS = (nX = 32, nV = 33)               # hats in x and v
const NG = 32                                 # charge shape functions in x

function sizes(d)
    (;
        kick_drift = S * (KD.nξ^d + KD.nv^d) * KD.nw,       # (n_ξ^d + n_v^d) n_w per stage
        non_separable = S * NS.nξ^d * NS.nv^d * NS.nℓ,      # n_ξ^d n_v^d n_ℓ per stage
        hat_nodes = (HATS.nX * HATS.nV)^d,                  # nodes of f_h, g_h and φ
        charge_shapes = NG^d,                               # multipliers λ_i
        hats_per_marker = 2^(2d),                           # bilinear hats in 2d dimensions
        kick_terms = 4^(d + 1),                             # cubic splines in ξ (or v) and ln w
        nonseparable_terms = 4^(2d + 1),                    # cubic splines in ξ, v and ln w
        charge_terms = 4^d,                                 # cubic shape functions in x
        hessian_row = 7^d)                                  # nonzeros of a row of A diag(w) Aᵀ
end

# the counts of d = 1 are those of the code
w = ParticleRemapping.Projection.hole_filling_problem().w0
s1 = sizes(1)
@assert s1.kick_drift ==
        remap_parameter_count(KickDriftParametrisation(w; KD..., stages = S))
@assert s1.non_separable == remap_parameter_count(NonSeparableParametrisation(w; NS...,
    stages = S))
@assert s1.hat_nodes == size(moser_laplacian(HatProjection(HATS.nX, HATS.nV)), 1)

println("quantity                                  d = 1          d = 2          d = 3")
for (k, name) in (:kick_drift => "kick–drift parameters, two stages",
    :non_separable => "non-separable parameters, two stages",
    :hat_nodes => "nodes of the hat projection",
    :charge_shapes => "charge shape functions",
    :hats_per_marker => "hats at one marker",
    :kick_terms => "spline terms of one kick at one marker",
    :nonseparable_terms => "spline terms of one midpoint iteration",
    :charge_terms => "shape functions at one marker",
    :hessian_row => "nonzeros of a row of the charge Hessian")
    @printf("%-40s %14d %14d %14d\n", name, (getfield(sizes(d), k) for d in 1:3)...)
end
s3 = sizes(3)
@printf("\nd = 3: the normal matrix of the fit has %.2e entries (kick–drift), %.2e (non-separable)\n",
    Float64(s3.kick_drift)^2, Float64(s3.non_separable)^2)
@printf("d = 3: the dense Laplacian of the density-interpolation flow has %.2e entries\n",
    Float64(s3.hat_nodes)^2)
@printf("d = 3: one field on the nodes of the hats takes %.1f GB in Float64\n",
    8 * s3.hat_nodes / 1e9)

# the grid of the low-rank density of the test after a push, 64 points in x and 256 in v
const NX, NV = 64, 256
for d in 1:3
    @printf("d = %d: full grid of the low-rank test %.3e points; low rank: %d values per unit of rank\n",
        d, Float64(NX * NV)^d, NX^d + NV^d)
end
