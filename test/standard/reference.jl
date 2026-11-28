using ParticleRemapping
using Test

# The message that `f()` throws as an exception of type `T`, as `showerror` prints it.
function shown(T, f)
    try
        f()
    catch e
        e isa T && return sprint(showerror, e)
        rethrow()
    end
    "no exception"
end

# The n × n cell centres of [−L/2, L/2) × [−V, V], i outer, j inner.
function lattice(n, L, V)
    dx, dv = L / n, 2V / n
    ([-L / 2 + (i - 0.5) * dx for i in 1:n for j in 1:n],
        [-V + (j - 0.5) * dv for i in 1:n for j in 1:n])
end

@testset "GridResample: the lumped variant" begin
    method = GridResample(; L = 2π, V = 3.0)
    @test method.projection.nX == 16 && method.projection.nV == 17
    @test method.density === nothing
    for (N, T) in ((900, 25), (144, 5))
        mk = wave_case(; N, T).markers
        r = remap(method, mk)
        a = r.markers
        @test abs(sum(a.w) - sum(mk.w)) <= 1e-14 * sum(mk.w)
        x, v = lattice(isqrt(N), 2π, 3.0)
        @test a.x == x && a.v == v
        @test r.info.variant === :lumped
        @test r.info.zero_weights == count(iszero, a.w)
        @test a.m == 1.2
        @test r.info.scale ≈ sum(mk.w) / r.info.raw_total
        N == 900 && @test r.info.zero_weights > 0
        N == 900 && @test isapprox(r.info.raw_total / sum(mk.w), 1.001076; rtol = 1e-6)
    end
    one = Markers(zeros(9), zeros(9), [1.0; zeros(8)])
    a = remap(method, one).markers
    x, v = lattice(3, 2π, 3.0)
    @test a.x == x && a.v == v
    @test abs(a.w[5] - 1) < 1e-14
    @test all(iszero, a.w[[1:4; 6:9]])
end

@testset "GridResample: the exact variant" begin
    case = wave_case(; N = 900, T = 25)
    mk = case.markers
    f = case.positive
    r = remap(GridResample(; L = 2π, V = 3.0, density = f), mk)
    x, v = lattice(30, 2π, 3.0)
    @test r.markers.x == x && r.markers.v == v
    @test isapprox(r.markers.w, f.(x, v) .* (sum(mk.w) / sum(f.(x, v))); rtol = 1e-12)
    @test r.info.variant === :exact
    @test r.info.zero_weights == 0
    hole = hole_case(; k = 4)
    r = remap(GridResample(; L = 2.0, V = 1.0, density = hole.positive), hole.markers)
    @test all(isapprox.(r.markers.w, 1 / 400; rtol = 1e-14))
    @test r.markers.m == 1.0
    # zero weights where the density vanishes
    r = remap(GridResample(; L = 2.0, V = 1.0, density = hole.represented), hole.markers)
    @test r.info.zero_weights == 200
    @test isapprox(sum(r.markers.w), 1; rtol = 1e-14)
end

@testset "GridResample: the edges" begin
    method = GridResample(; L = 2π, V = 2.5)
    for N in (150, 0)
        mk = Markers(zeros(N), zeros(N), ones(N))
        @test shown(ArgumentError, () -> remap(method, mk)) ==
              "ArgumentError: GridResample places an n × n lattice, so N must be a square " *
              "n^2 >= 1, got N = $N"
    end
    @test shown(ArgumentError, () -> remap(method, Markers(zeros(2, 16), zeros(2, 16),
        ones(16)))) == "ArgumentError: GridResample acts on 1D-1V markers"
    zero = GridResample(; L = 2π, V = 2.5, density = (x, v) -> 0.0)
    @test startswith(shown(DomainError, () -> remap(zero, wave_case(; N = 16).markers)),
        "DomainError with 0.0:\nGridResample needs a density that is positive at some " *
        "lattice point")
    # a marker beyond the box lands on the edge nodes of the hats and raises the top row: the
    # 4 × 4 lattice reads the nodes 3, 7, 11 and 15 of 17 in r, so the mass on the nodes 16 and
    # 17 leaves the raw total, and the scale c raises every other weight
    mk = wave_case(; N = 16, T = 1.25).markers
    k = argmin(mk.v)
    v = copy(mk.v)
    v[k] = 3.7
    out = Markers(mk.x, v, mk.w; m = mk.m)
    hp = method.projection
    g = ParticleRemapping.Projection.lumped_projection(hp, out.x ./ 2π, out.v ./ 2.5, out.w)
    g0 = ParticleRemapping.Projection.lumped_projection(hp, mk.x ./ 2π, mk.v ./ 2.5, mk.w)
    @test sum(g[:, 16:17]) > sum(g0[:, 16:17])
    a = remap(method, out).markers
    a0 = remap(method, mk).markers
    @test abs(sum(a.w) - sum(mk.w)) <= 1e-14 * sum(mk.w)
    top = a.v .== maximum(a.v)
    @test count(top) == 4
    @test sum(a.w[top]) > sum(a0.w[top])
end
