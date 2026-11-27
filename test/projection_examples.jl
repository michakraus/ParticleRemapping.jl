using ParticleRemapping
using LinearAlgebra: norm
using Random: Xoshiro
using Test

# the density of the lattice example and of the hole problem
f_lattice(x, v) = (1 + 0.3 * cos(x)) * exp(-v^2 / 2) / sqrt(2π)

# the markers of `lattice_example`, computed from the cell centres (ξ, η) of its plain lattice
function moved_lattice(ξ, η, nL; a, b, L = 2π, V = 4.0)
    x = mod.(ξ .+ a .* sin.(ξ .- 0.5), L)
    v = η .+ b .* sin.(π .* (η .+ V) ./ (2V))
    J = (1 .+ a .* cos.(ξ .- 0.5)) .* (1 .+ (π * b / (2V)) .* cos.(π .* (η .+ V) ./ (2V)))
    (; x, v, w = (2L * V / nL^2) .* f_lattice.(x, v) .* J)
end

@testset "lattice_example: the deformed markers" begin
    for nL in (24, 10)
        ex = lattice_example(; nL)
        @test ex.N == nL^2 == length(ex.markers) == length(ex.initial)
        @test ex.L ≈ 2π
        @test ex.V == 4.0
        @test ex.m == 1.0 == ex.markers.m == ex.initial.m
        ref = moved_lattice(ex.initial.x, ex.initial.v, nL; a = 0.2, b = 0.3)
        @test isapprox(ex.markers.x, ref.x; atol = 1e-12, rtol = 1e-12)
        @test isapprox(ex.markers.v, ref.v; atol = 1e-12, rtol = 1e-12)
        @test isapprox(ex.markers.w, ref.w; atol = 1e-12, rtol = 1e-12)
        plain = lattice_example(; nL, a = 0, b = 0).markers
        @test isapprox(ex.initial.x, plain.x; rtol = 1e-14)
        @test isapprox(ex.initial.v, plain.v; rtol = 1e-14)
        @test isapprox(ex.initial.w, plain.w; rtol = 1e-14)
        @test all(x -> 0 <= x < ex.L, ex.markers.x)
        @test maximum(abs, ex.markers.v) < ex.V
        @test abs(sum(ex.markers.w) - sum(ex.initial.w)) / sum(ex.initial.w) < 1e-6
    end
end

@testset "lattice_example: the plain lattice is the hole problem's lattice" begin
    for (nL, kw) in ((10, (; hole_x = π / 2, hole_a = 0.9, hole_b = 0.8)), (48, (;)))
        ex = lattice_example(; nL, a = 0, b = 0)
        prob = hole_filling_problem(; nL, kw...)
        keep = [!inside_hole(prob.hole, x, v) for (x, v) in zip(ex.markers.x, ex.markers.v)]
        @test count(keep) == prob.N
        @test ex.markers.x[keep] == prob.markers.x
        @test ex.markers.v[keep] == prob.markers.v
        @test ex.markers.w[keep] == prob.markers.w
        @test isapprox(ex.markers.x, ex.initial.x; rtol = 1e-14)
        @test isapprox(ex.markers.v, ex.initial.v; rtol = 1e-14)
        @test isapprox(ex.markers.w, ex.initial.w; rtol = 1e-14)
        for (x, v) in ((0.1, 0.2), (3.0, -1.0), (6.0, 3.9), (π, 5.0), (1.0, -7.0))
            @test isapprox(ex.density(x, v), prob.density(x, v); rtol = 1e-14)
        end
    end
end

@testset "lattice_example: the wrap of x" begin
    ex = lattice_example(; a = 0.9)
    @test all(x -> 0 <= x < ex.L, ex.markers.x)
    ξ = ex.initial.x
    @test isapprox(ex.markers.x, mod.(ξ .+ 0.9 .* sin.(ξ .- 0.5), ex.L); atol = 1e-12)
end

@testset "lattice_example: a 1 × 1 lattice" begin
    ex = lattice_example(; nL = 1)
    @test ex.N == 1
    x, v = π + 0.2 * sin(π - 0.5), 0.3
    @test isapprox(ex.markers.x[1], x; rtol = 1e-12)
    @test isapprox(ex.markers.v[1], v; rtol = 1e-12)
    J = (1 + 0.2 * cos(π - 0.5)) * (1 + (π * 0.3 / 8) * cos(π / 2))
    @test isapprox(ex.markers.w[1], 16π * f_lattice(x, v) * J; rtol = 1e-12)
end

@testset "lattice_example: an empty lattice" begin
    @test_throws r"^ArgumentError: lattice_example needs nL ≥ 1, got nL = 0$" lattice_example(;
        nL = 0)
    @test_throws r"^ArgumentError: lattice_example needs nL ≥ 1, got nL = -1$" lattice_example(;
        nL = -1)
end

@testset "lattice_example: the bounds of a and b" begin
    ma(s) = Regex("^ArgumentError: lattice_example needs \\|a\\| < 1 for a monotone map in x, " *
                  "got a = $(s)\$")
    mb(s) = Regex("^ArgumentError: lattice_example needs \\|b\\| < 8/π for a monotone map in v, " *
                  "got b = $(s)\$")
    @test_throws ma("1\\.0") lattice_example(; a = 1.0)
    @test_throws ma("-1\\.0") lattice_example(; a = -1.0)
    @test_throws ma("NaN") lattice_example(; a = NaN)
    @test_throws mb("2\\.6") lattice_example(; b = 2.6)
    @test_throws mb("-2\\.6") lattice_example(; b = -2.6)
    @test_throws mb("NaN") lattice_example(; b = NaN)
    for kw in ((; a = 0.9), (; a = -0.9), (; b = 2.5))
        ex = lattice_example(; kw...)
        @test all(>(0), ex.markers.w)
        @test all(x -> 0 <= x < ex.L, ex.markers.x)
        @test maximum(abs, ex.markers.v) < ex.V
    end
end

# the density f0 of the spline chapter's example, and the leapfrog steps in its force
function f0_clustered(x, v; Lx = 2π, amplitude = 0.45)
    (1 + amplitude * cos(x)) * exp(-0.5 * v^2) /
    (Lx * sqrt(2π))
end
force_clustered(x) = 1.4 * sin(x - 0.5) + 0.8 * sin(2x)
function leapfrog(x, v; nsteps = 60, dt = 0.05, m = 1.5, Lx = 2π)
    for _ in 1:nsteps
        v += 0.5 * dt * force_clustered(x) / m
        x = mod(x + dt * v, Lx)
        v += 0.5 * dt * force_clustered(x) / m
    end
    x, v
end

@testset "clustered_example: the clustered 24 × 24 markers" begin
    ex = clustered_example()
    ref = redistribution_example(; nx0 = 24, nv0 = 24)
    @test ex.N == 576 == length(ex.markers)
    @test ex.L ≈ 2π
    @test ex.V == 3.5
    @test ex.m == 1.5 == ex.markers.m
    for (a, b) in ((ex.markers, ref.markers), (ex.initial, ref.initial))
        @test a.x == b.x
        @test a.v == b.v
        @test a.w == b.w
    end
    @test ex.markers.w == ex.initial.w
    @test norm(ex.markers.x - ex.initial.x) > 1
    @test maximum(abs, ex.markers.v) < ex.V
end

@testset "clustered_example: the density is the exact inverse of the clustering" begin
    for (kw, Lx, nx0, nv0) in (((;), 2π, 24, 24), ((; nx0 = 10, nv0 = 10), 2π, 10, 10),
        ((; nx0 = 6, nv0 = 4, Lx = 4.0), 4.0, 6, 4))
        ex = clustered_example(; kw...)
        Δx, Δv = Lx / nx0, 4.8 / nv0
        @test isapprox(ex.density.(ex.markers.x, ex.markers.v) .* Δx .* Δv, ex.markers.w;
            rtol = 1e-12)
    end
    density = clustered_example().density
    rng = Xoshiro(20261123)
    for _ in 1:50
        x0, v0 = 2π * rand(rng), -5 + 10 * rand(rng)
        @test isapprox(density(leapfrog(x0, v0)...), f0_clustered(x0, v0); rtol = 1e-12)
    end
end

@testset "clustered_example: a 1 × 1 and a non-square lattice" begin
    ex = clustered_example(; nx0 = 1, nv0 = 1)
    ref = redistribution_example(; nx0 = 1, nv0 = 1)
    @test ex.N == 1
    @test ex.markers.x == ref.markers.x && ex.markers.v == ref.markers.v
    @test isapprox(
        ex.density(ex.markers.x[1], ex.markers.v[1]) * 2π * 4.8, ex.markers.w[1];
        rtol = 1e-12)
    ex = clustered_example(; nx0 = 6, nv0 = 4)
    @test ex.N == 24
    @test isapprox(ex.density.(ex.markers.x, ex.markers.v) .* (2π / 6) .* (4.8 / 4),
        ex.markers.w; rtol = 1e-12)
end

@testset "clustered_example: an empty clustered lattice" begin
    @test_throws r"^ArgumentError: clustered_example needs nx0 ≥ 1 and nv0 ≥ 1, got nx0 = 0, nv0 = 24$" clustered_example(;
        nx0 = 0)
    @test_throws r"^ArgumentError: clustered_example needs nx0 ≥ 1 and nv0 ≥ 1, got nx0 = 24, nv0 = 0$" clustered_example(;
        nv0 = 0)
end

@testset "clustered_example: no force, and no steps" begin
    ex = clustered_example(; nx0 = 10, nv0 = 10, force = x -> 0.0)
    @test ex.markers.v == ex.initial.v
    @test isapprox(ex.markers.x, mod.(ex.initial.x .+ 3.0 .* ex.initial.v, 2π); atol = 1e-12)
    Δ = (2π / 10) * (4.8 / 10)
    @test isapprox(ex.density.(ex.markers.x, ex.markers.v) .* Δ, ex.markers.w; rtol = 1e-12)
    ex = clustered_example(; nx0 = 10, nv0 = 10, nsteps = 0)
    @test ex.markers.x == ex.initial.x
    @test ex.markers.v == ex.initial.v
    @test isapprox(ex.density.(ex.markers.x, ex.markers.v),
        f0_clustered.(ex.markers.x, ex.markers.v); rtol = 1e-14)
end

@testset "clustered_example: the density outside the velocity box" begin
    density = clustered_example().density
    for v in (3.5, -3.5, 4.0, 6.0)
        d = density(π, v)
        @test isfinite(d) && d > 0
    end
end

@testset "clustered_example: a box that does not hold the markers" begin
    @test_throws r"^ArgumentError: clustered_example needs V ≥ 2\.75356\d*, the largest \|v\| of the clustered markers, got V = 2\.0$" clustered_example(;
        nx0 = 10, nv0 = 10, V = 2.0)
    # a box whose edge is the largest |v| holds the markers
    vmax = maximum(abs, clustered_example(; nx0 = 10, nv0 = 10).markers.v)
    @test clustered_example(; nx0 = 10, nv0 = 10, V = vmax).V == vmax
end

# the method of the SVG steps of the projection
function steps_method(ex; kw...)
    KickDriftRemap(; L = ex.L, V = ex.V,
        projection = HatProjection(8, 9), nξ = 16, nv = 20, nw = 8, maxit = 20,
        solver = :gauss_newton, density = (xs, vs) -> ex.density.(xs, vs), kw...)
end

@testset "the remap of the 10 × 10 lattice example" begin
    ex = lattice_example(; nL = 10)
    r = remap(steps_method(ex), ex.markers)
    info = r.info
    @test all(>(0), info.targets.w)
    @test string(info.status.code) == "MAXITERS"
    @test info.iterations == 20
    @test info.F < info.F0
    @test info.objective_ratio < 0.02
    @test info.momentum_change < 1e-14
    @test info.weight_change < 1e-14
    @test length(r.markers) == 100
    @test maximum(abs, r.markers.v) < 4.0
end

@testset "the remap of the plain 10 × 10 lattice: the fitted map is the identity" begin
    ex = lattice_example(; nL = 10, a = 0, b = 0)
    r = remap(steps_method(ex), ex.markers)
    info = r.info
    @test string(info.status.code) == "SUCCESS"
    @test info.status.stationary
    @test info.iterations == 1
    @test maximum(abs, info.θ) < 1e-12
    @test abs(info.objective_ratio - 1) < 1e-12
    @test isapprox(r.markers.x, ex.markers.x; atol = 1e-12)
end

@testset "the remap of the clustered example on a basis finer than the marker spacing" begin
    ex = clustered_example(; nx0 = 10, nv0 = 10)
    fine(density) = KickDriftRemap(;
        L = ex.L, V = ex.V, projection = HatProjection(24, 25),
        nξ = 16, nv = 20, nw = 8, maxit = 2, solver = :gauss_newton, density)
    r = remap(fine((xs, vs) -> ex.density.(xs, vs)), ex.markers)
    @test all(>(0), r.info.targets.w)
    @test_throws r"^DomainError with .*the target weights must be positive"s remap(
        fine(nothing), ex.markers)
end
