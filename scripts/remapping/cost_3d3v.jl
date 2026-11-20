# The size and the cost of one Newton step of the gauge rezoning in d dimensions.
#
# Checks the counts of the 1D-1V implementation against the matrices that the package assembles
# (non-zeros of the Jacobian G, size and stencil of the Schur complement), then evaluates the same
# counts in d dimensions for tensor-product linear hat functions on N_g^d nodes, each marker on
# the 2^d nodes of its cell, and n_K scalar moments per node. Prints the numbers that the section
# "Higher dimensions" of the chapter quotes, with the average number of markers per group of cell
# and orthant of the targets in 3d3v.
#
# Run, from the package root:
#   julia --startup-file=no --project=scripts scripts/remapping/cost_3d3v.jl

using ParticleRemapping
using Printf
using SparseArrays: nnz, sparse

results = Pair{String, Bool}[]
function check(name, ok, detail = "")
    push!(results, name => ok)
    println(rpad(name, 64), detail, ok ? "  ok" : "  UNEXPECTED")
end

# the counts of one Newton step in d dimensions
function counts(; d, Ng, N, nK)
    nodes = Ng^d
    nc = nK * nodes
    (; unknowns = N + nc, constraints = nc,
        nnz_G = 2^d * nK * N,                 # each column: the 2^d nodes of the cell, nK moments
        assembly = 4^d * nK^2 * N,            # products added into G diag(H)⁻¹ Gᵀ, per marker
        nnz_schur = 3^d * nK^2 * nodes,       # each node couples with the 3^d nodes around it
        dense_flops = 2 / 3 * float(nc)^3,    # the LU factorisation of the package's dense solve
        dense_bytes = 8 * float(nc)^2,
        cg_flops = 2 * 3^d * nK^2 * nodes)    # one product with the sparse Schur complement
end

# 1D-1V: the package's own matrices against the counts
println("1D-1V: the package's matrices against the counts")
for (Ng, mpn, K) in ((8, 50, (0, 2)), (16, 40, (0, 2, 3)))
    s = gauge_rezoning_setup(; Ng, markers_per_node = mpn, K)
    mk = s.markers
    prob = GaugeRezoningProblem(s.method, mk)
    μ = 0.1 .* sin.(1:length(mk))
    _, G = gauge_constraints(prob, μ)
    H = gauge_hessian(prob, μ, zeros(length(K) * Ng))
    S = sparse(Matrix(gauge_schur_matrix(H, G)))
    c = counts(; d = 1, Ng, N = length(mk), nK = length(K))
    check("  N_g = $Ng, |K| = $(length(K)): nnz(G) = 2|K|N",
        nnz(G) == c.nnz_G, "$(nnz(G))")
    check("  N_g = $Ng, |K| = $(length(K)): Schur complement |K| N_g square",
        size(S) == (c.constraints, c.constraints), "$(size(S, 1))")
    check("  N_g = $Ng, |K| = $(length(K)): nnz(Schur) = 3 |K|² N_g",
        nnz(S) == c.nnz_schur, "$(nnz(S))")
end

# control: the Schur complement of the package is not banded without the periodic wrap
s = gauge_rezoning_setup(; Ng = 8, markers_per_node = 50)
prob = GaugeRezoningProblem(s.method, s.markers)
_, G = gauge_constraints(prob, zeros(length(s.markers)))
S = gauge_schur_matrix(prob.w, G)
check("Control: node 1 couples with node N_g (periodic wrap)", !iszero(S[1, 8]),
    @sprintf("%.2e", S[1, 8]))

println("\nOne Newton step in d dimensions, N_g = 32 nodes per dimension, N = 10^7 markers;")
println("n_K = 2: charge and kinetic energy, n_K = 1 + d(d+1)/2: charge and the stress tensor")
@printf("%-3s %-4s %12s %12s %12s %12s %12s %12s\n", "d", "n_K", "constraints", "nnz(G)",
    "assembly", "nnz(Schur)", "dense flops", "dense bytes")
for d in (1, 2, 3), nK in unique((2, 1 + d * (d + 1) ÷ 2))

    c = counts(; d, Ng = 32, N = 10^7, nK)
    @printf("%-3d %-4d %12d %12.2e %12.2e %12.2e %12.2e %12.2e\n", d, nK, c.constraints,
        c.nnz_G, c.assembly, c.nnz_schur, c.dense_flops, c.dense_bytes)
end
c3 = counts(; d = 3, Ng = 32, N = 10^7, nK = 2)
@printf("\n3d3v, K = {charge, energy}: one product with the sparse Schur complement %.2e flops\n",
    c3.cg_flops)
@printf("3d3v: markers per group of cell and orthant, N / (N_g^3 2^3) = %.1f\n",
    10^7 / (32^3 * 2^3))

nfail = count(!last, results)
println("\n", length(results) - nfail, "/", length(results), " checks as expected")
nfail == 0 || error("$nfail check(s) not as expected")
