# The diagnostics of the hole-filling experiments against closed forms.

using ParticleRemapping
using Test

@testset "Diagnostics" begin
    L, V, n = 2π, 4.0, 48
    dx, dv = L / n, 2V / n
    lattice = [((i - 0.5) * dx, -V + (j - 0.5) * dv) for i in 1:n for j in 1:n]
    x, v = first.(lattice), last.(lattice)

    @testset "largest empty disk" begin
        # a perfect lattice: the centre of a cell is √2/2 lattice units from its corners
        @test largest_empty_disk(x, v; L, dx, dv) ≈ sqrt(0.5) rtol = 1e-12
        # positions are taken modulo L
        @test largest_empty_disk(x .+ L, v; L, dx, dv) ≈ sqrt(0.5) rtol = 1e-12
        hole = EllipticHole(π, 1.2, 0.45, 0.40)
        @test local_empty_disk(hole, x, v; dx, dv) ≈ sqrt(0.5) rtol = 1e-12
        keep = [!inside_hole(hole, x[k], v[k]) for k in eachindex(x)]
        @test count(!, keep) == 24
        # the hole of the experiment leaves a disk of 2.75 lattice units
        @test round(local_empty_disk(hole, x[keep], v[keep]; dx, dv); digits = 2) == 2.75
    end

    @testset "hole and fill" begin
        hole = EllipticHole(π / 2, 1.2, 0.45, 0.40)
        @test inside_hole(hole, π / 2, 1.2) && inside_hole(hole, π / 2 + L, 1.2)
        @test !inside_hole(hole, π / 2 + 0.46, 1.2) && !inside_hole(hole, π / 2, 1.61)
        hx, hv = hole_outline(hole; n = 9)
        @test all(((hx .- hole.x) ./ hole.a) .^ 2 .+ ((hv .- hole.v) ./ hole.b) .^ 2 .≈ 1)
        # the fill of the lattice is its count in the hole relative to π a b N / |Ω|
        expected = count(k -> inside_hole(hole, x[k], v[k]), eachindex(x)) /
                   (π * 0.45 * 0.40 * n^2 / (2L * V))
        @test hole_fill(hole, x, v; volume = 2L * V) ≈ expected rtol = 8eps()
        @test hole_fill(hole, x, v; volume = 2L * V) ≈ 1 atol = 0.1
    end

    @testset "density discrepancy" begin
        hp = HatProjection(32, 33)
        # a constant density on a lattice that the hats resolve exactly
        xl = [(i - 1) * L / 32 for i in 1:32 for j in 1:33]
        vl = [-V + (j - 1) * 2V / 32 for i in 1:32 for j in 1:33]
        # every inner node collects one marker: g_h = 1/(hs hr) = 512 against the mean N/2 = 528
        @test density_discrepancy(hp, xl, vl; L, V) ≈ 16 / 528 rtol = 1e-9
        @test density_discrepancy(hp, x, v; L, V) > 0
    end

    @testset "weight errors by band" begin
        vb = [0.5, 0.5, 1.5, 2.5, 3.5]
        w = [1.1, 0.9, 1.2, 1.0, 0.5]
        e = weight_errors_by_band(vb, w, ones(5))
        @test e[1].rms ≈ 0.1 rtol = 8eps()
        @test e[1].max ≈ 0.1 rtol = 8eps()
        @test e[2].rms ≈ 0.2 rtol = 8eps()
        @test e[3].rms == 0
        @test e[4].max ≈ 0.5 rtol = 8eps()
    end
end
