# R2: the de-clustering relaxation of the target points lowers the clustering.

using ParticleRemapping
using Test

ex = redistribution_example()
x, v = ex.markers.x, ex.markers.v
cv(xs, vs) = kernel_density_cv(mod.(xs, ex.Lx), vs; Lx = ex.Lx, σx = ex.σx, σv = ex.σv)

@testset "Relaxation lowers the clustering" begin
    relaxed = declustering_relaxation(x, v; Lx = ex.Lx, σx = ex.σx, σv = ex.σv)
    @test cv(relaxed.x, relaxed.v) < cv(x, v)
    @test length(relaxed.history) == 11
    @test relaxed.history[1] == (x, v)
    @test relaxed.history[end] == (relaxed.x, relaxed.v)
    # damped gradient steps lower the repulsion energy at every step
    E = [declustering_energy(h...; Lx = ex.Lx, σx = ex.σx, σv = ex.σv)
         for h in relaxed.history]
    @test all(diff(E) .< 0)
end

@testset "The research log's relaxation: 0.220 → 0.158" begin
    # the numbers of the log come from a relaxation step 0.02 with eight steps
    log_case = declustering_relaxation(
        x, v; Lx = ex.Lx, σx = ex.σx, σv = ex.σv, step = 0.02,
        iterations = 8)
    @test cv(x, v) ≈ 0.220 atol = 5e-4
    @test cv(log_case.x, log_case.v) ≈ 0.158 atol = 5e-4
end

@testset "Soft wall confines the velocities" begin
    # a point beyond the wall is pushed back, and the wall adds energy
    xs, vs = [0.0, 3.0], [2.6, -0.1]
    r = declustering_relaxation(
        xs, vs; Lx = 2π, σx = 0.1, σv = 0.1, step = 0.1, iterations = 1)
    @test r.v[1] ≈ 2.6 - 0.1 * 2.5 * 0.6 rtol = 1e-12
    @test declustering_energy(xs, vs; Lx = 2π, σx = 0.1, σv = 0.1) >= 2.5 / 2 * 0.6^2
end
