# The two examples of the projection overview page (`docs/src/projection/index.md`): the mildly
# deformed 24 × 24 lattice (`lattice_example()`) and the 24 × 24 clustered markers
# (`clustered_example()`), each remapped by the same 16/20/8 kick–drift map on 12 × 13 hats. Checks
# every number the overview quotes with its digits or inside its range, verifies that each quoted
# string appears on `docs/src/projection/index.md`, checks conservation of total weight and total
# momentum, and verifies that the lumped projection on 12 × 13 hats refuses the clustered markers
# with a `DomainError`.

using ParticleRemapping
using Test

# One unit of the last digit of the number `quoted`.
function quoted_unit(quoted::AbstractString)
    m = match(r"^-?(\d+)(?:\.(\d+))?(?:e(-?\d+))?$", quoted)
    decimals = m.captures[2] === nothing ? 0 : length(m.captures[2])
    exponent = m.captures[3] === nothing ? 0 : parse(Int, m.captures[3])
    10.0^(exponent - decimals)
end

# `x` agrees with `quoted` to the digits that `quoted` shows: within half a unit of its last digit.
# A range `lo–hi` holds `x` between its rounded ends.
function agrees(x, quoted::AbstractString)
    if occursin('–', quoted)
        lo, hi = parse.(Float64, split(quoted, '–'))
        return lo <= x <= hi
    end
    abs(x - parse(Float64, quoted)) <= 0.5 * quoted_unit(quoted) * (1 + 1e-9)
end

# The LaTeX inline-math token for a quoted number or range on `docs/src/projection/index.md`.
function doc_math(s::AbstractString)
    m = match(r"^(-?\d+(?:\.\d+)?)e(-?\d+)$", s)
    m === nothing ? "``$s``" : "``$(m.captures[1]) \\cdot 10^{$(m.captures[2])}``"
end

function doc_token(quoted::AbstractString)
    occursin('–', quoted) ? join(doc_math.(split(quoted, '–')), '–') : doc_math(quoted)
end

function median_value(v::AbstractVector)
    s = sort(v)
    n = length(s)
    isodd(n) ? s[cld(n, 2)] : 0.5 * (s[n ÷ 2] + s[n ÷ 2 + 1])
end

# The section of `page` from the heading carrying `anchor` to the next `## ` heading.
function page_section(page::AbstractString, anchor::AbstractString)
    lines = split(page, '\n')
    i = findfirst(l -> occursin("(@id $anchor)", l), lines)
    i === nothing && error("anchor $anchor not found")
    j = findnext(l -> startswith(l, "## "), lines, i + 1)
    join(lines[i:(j === nothing ? end : j - 1)], '\n')
end

const QUOTED_MAIN = (
    code = "MAXITERS",
    iterations = "20",
    K0 = "3.138",
    objective_ratio = "0.220",
    dK = "0.431",
    outside = "7",
    max_v = "4.5"
)

const QUOTED_CLUSTERED = (
    outside = "210–224",
    max_v = "130–180",
    weight_share = "0.094–0.11",
    min_target_w = "2.9e-6",
    median_target_w = "1.2e-4–1.4e-4",
    mean_w = "1.7e-3"
)

@testset "Overview examples" begin
    ex = lattice_example()
    cl = clustered_example()

    @testset "sizes and boxes" begin
        @test ex.N == 576
        @test cl.N == 576
        @test ex.V == 4.0
        @test cl.V == 3.5
    end

    @testset "main example: deformed lattice" begin
        method = KickDriftRemap(; L = ex.L, V = ex.V, projection = HatProjection(12, 13),
            nξ = 16, nv = 20, nw = 8, maxit = 20,
            density = (xs, vs) -> ex.density.(xs, vs))
        r = remap(method, ex.markers)
        i = r.info
        K0 = kinetic_energy(ex.markers)
        dK = kinetic_energy(r.markers) - K0
        out = abs.(r.markers.v) .> ex.V

        @test string(i.status.code) == QUOTED_MAIN.code
        @test !i.status.stationary
        @test i.iterations == parse(Int, QUOTED_MAIN.iterations)
        @test agrees(K0, QUOTED_MAIN.K0)
        @test agrees(i.objective_ratio, QUOTED_MAIN.objective_ratio)
        @test agrees(dK, QUOTED_MAIN.dK)
        @test count(out) == parse(Int, QUOTED_MAIN.outside)
        @test agrees(maximum(abs, r.markers.v), QUOTED_MAIN.max_v)
        @test i.weight_change <= 1e-13
        @test i.momentum_change <= 1e-13
        @test abs(total_weight(r.markers) - total_weight(ex.markers)) /
              total_weight(ex.markers) <= 1e-13
        @test abs(total_momentum(r.markers) - total_momentum(ex.markers)) /
              sum(abs, momentum(ex.markers)) <= 1e-13
        rep = conservation_report(ex.markers, r.markers; grid = PeriodicGrid(2π, 12))
        @test rep.kinetic_energy ≈ dK rtol = 8eps()
    end

    @testset "second example: clustered markers" begin
        method = KickDriftRemap(; L = cl.L, V = cl.V, projection = HatProjection(12, 13),
            nξ = 16, nv = 20, nw = 8, maxit = 20,
            density = (xs, vs) -> cl.density.(xs, vs))
        r = remap(method, cl.markers)
        out = abs.(r.markers.v) .> cl.V
        w_out = r.markers.w[out]

        @test string(r.info.status.code) == "MAXITERS"
        @test r.info.iterations == 20
        @test agrees(count(out), QUOTED_CLUSTERED.outside)
        @test agrees(maximum(abs, r.markers.v), QUOTED_CLUSTERED.max_v)
        @test agrees(sum(w_out) / sum(r.markers.w), QUOTED_CLUSTERED.weight_share)
        @test agrees(minimum(w_out), QUOTED_CLUSTERED.min_target_w)
        @test agrees(median_value(w_out), QUOTED_CLUSTERED.median_target_w)
        @test agrees(sum(cl.markers.w) / cl.N, QUOTED_CLUSTERED.mean_w)

        # The lumped projection on 12 × 13 hats vanishes at a target of the clustered markers.
        @test_throws r"^DomainError(?s:.)*the target weights must be positive" remap(
            KickDriftRemap(; L = cl.L, V = cl.V, projection = HatProjection(12, 13),
                nξ = 16, nv = 20, nw = 8, maxit = 2),
            cl.markers)
    end

    @testset "overview page quotes every number" begin
        page = read(
            joinpath(@__DIR__, "..", "..", "docs", "src", "projection", "index.md"),
            String)
        main_sec = page_section(page, "projection-example")
        @test occursin("`$(QUOTED_MAIN.code)`", main_sec)
        @test occursin("scripts/projection/example.jl", main_sec)
        for q in (QUOTED_MAIN.iterations, QUOTED_MAIN.K0, QUOTED_MAIN.objective_ratio,
            QUOTED_MAIN.dK, QUOTED_MAIN.outside, QUOTED_MAIN.max_v)
            @test occursin(doc_token(q), main_sec)
        end

        cl_sec = page_section(page, "projection-clustered-example")
        @test occursin("scripts/projection/example.jl", cl_sec)
        for q in values(QUOTED_CLUSTERED)
            @test occursin(doc_token(q), cl_sec)
        end
    end
end
