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

@testset "the evenness of the wave rungs" begin
    case = wave_case(; N = 144, T = 0)
    e = case_evenness(case, case.markers)
    @test e.empty == 40
    @test isapprox(e.rstd, 1.0773; rtol = 5e-3)
    @test e.support_cells == 40
    @test e.support_empty == 0
    @test isapprox(e.support_rstd, 0.3944; rtol = 5e-3)
    @test e.window.empty == 40
    @test isapprox(e.window.rstd, 0.9230; rtol = 5e-3)
    case = wave_case(; N = 900, T = 25)
    e = case_evenness(case, case.markers)
    @test e.empty in 30:32
    @test isapprox(e.rstd, 0.8498; rtol = 5e-3)
    @test e.support_cells == 44
    @test e.support_empty == 0
    @test isapprox(e.support_rstd, 0.2644; rtol = 5e-3)
    @test e.window.empty in 13:15
    @test isapprox(e.window.rstd, 0.5976; rtol = 5e-3)
    case = hole_case()
    e = case_evenness(case, case.markers)
    @test e.window === nothing
    @test e.support_cells == 32
    @test e.support_empty == 0
    @test e.empty == 32
end

@testset "markers outside the box and non-finite markers" begin
    case = wave_case(; N = 16, T = 0)
    o = case_outside(case, case.markers)
    @test o.count == 0 && o.weight_share == 0 && o.nonfinite == 0
    @test o.max_v == maximum(abs, case.markers.v)
    v, w = copy(case.markers.v), copy(case.markers.w)
    v[1] = 3.5
    v[2] = 3.0
    v[3] = NaN
    w[4] = Inf
    o = case_outside(case, Markers(case.markers.x, v, w; m = 1.2))
    @test o.count == 2
    @test o.nonfinite == 2
    @test o.max_v == 3.5
    @test isapprox(o.weight_share, w[1] / sum(w[isfinite.(w)]); rtol = 1e-14)
    x = copy(case.markers.x)
    x[5] = NaN
    @test case_outside(case, Markers(x, case.markers.v, case.markers.w; m = 1.2)).nonfinite ==
          1
    @test case_outside(case, Markers([0.0], [NaN], [1.0])).max_v == 0.0
end

@testset "the representation error and the charge noise" begin
    bounds = Dict(
        (900, 0, :represented) => ((0.0, 0.01), (0.0, 0.02)),
        (144, 0, :represented) => ((0.0, 0.01), (0.0, Inf)),
        (900, 0, :positive) => ((0.30, 0.36), (0.27, 0.33)),
        (900, 25, :represented) => ((0.05, 0.10), (0.06, 0.11)),
        (900, 25, :positive) => ((0.25, 0.31), (0.27, 0.33)))
    for ((N, T, variant), (rb, cb)) in bounds
        case = wave_case(; N, T)
        mk = case.markers
        ref = case_reference(case, variant)
        r = case_representation(case, mk, variant; reference = ref)
        c = case_charge_noise(case, mk, variant; reference = ref)
        @test rb[1] < r < rb[2]
        @test cb[1] < c < cb[2]
        @test isfinite(r) && isfinite(c)
        cp = copy(mk)
        @test case_representation(case, cp, variant; reference = ref) == r
        @test case_charge_noise(case, cp, variant; reference = ref) == c
    end
end

@testset "the precomputed reference" begin
    case = wave_case(; N = 144, T = 5)
    ref = case_reference(case, :represented)
    @test ref.variant === :represented
    @test size(ref.smoothed) == (30, 30)
    @test length(ref.charge) == 16
    for mk in (case.markers, wave_case(; N = 144, T = 5, seed = 1).markers)
        @test case_representation(case, mk, :represented; reference = ref) ==
              case_representation(case, mk, :represented)
        @test case_charge_noise(case, mk, :represented; reference = ref) ==
              case_charge_noise(case, mk, :represented)
    end
    msg = "ArgumentError: the reference is of the variant :represented, the measure asks for :positive"
    @test shown(ArgumentError,
        () -> case_representation(case, case.markers, :positive; reference = ref)) == msg
    @test shown(ArgumentError,
        () -> case_charge_noise(case, case.markers, :positive; reference = ref)) == msg
    coarse = case_reference(case, :represented; ngrid = 20, nq = 60, nq_charge = 60)
    @test_throws DimensionMismatch case_representation(case, case.markers, :represented;
        reference = coarse)
    short = (; ref.variant, ref.smoothed, charge = ref.charge[1:8])
    @test_throws DimensionMismatch case_charge_noise(case, case.markers, :represented;
        reference = short)
    @test_throws ArgumentError case_reference(case, :uniform)
end

@testset "the weights" begin
    s = case_weights(stress_case(:S2; N = 144).markers)
    @test isapprox(s.min, 4.8790e-08; rtol = 1e-3)
    @test isapprox(s.max, 1.7232e-01; rtol = 1e-3)
    @test isapprox(s.ratio, 3.5318e+06; rtol = 1e-3)
    @test isapprox(s.small_share, 8.010e-04; rtol = 1e-3)
    c = case_weights(hole_case().markers)
    @test c.min == c.max == 1 / 400
    @test c.ratio == 1
    @test c.small_share == 0
    @test case_weights(Markers([0.0, 1.0], [0.0, 0.0], [0.0, 1.0])).ratio == Inf
end

@testset "the conservation and the displacement" begin
    case = stress_case(:S3; N = 576)
    before = case.markers
    after = stress_case(:S3; N = 576, seed = 1).markers
    # the per-node report divides by empty nodes on this pair
    @test isnan(conservation_report(before, after; grid = PeriodicGrid(2π, 16)).charge)
    c = case_conservation(case, before, after)
    @test keys(c) == (:weight, :momentum, :kinetic_energy, :charge, :current)
    @test all(isfinite, values(c))
    @test c.charge > 0 && c.current > 0 && c.kinetic_energy > 0
    d = case_displacement(case, before, after)
    @test keys(d) == (:max_dx, :max_dv, :max_weight_change, :kinetic_energy_change)
    @test all(isfinite, values(d))
    @test 0 < d.max_dx <= π
    @test all(iszero, values(case_conservation(case, before, before)))
    @test all(iszero, values(case_displacement(case, before, before)))
    short = Markers(after.x[1:575], after.v[1:575], after.w[1:575]; m = after.m)
    @test_throws DimensionMismatch case_displacement(case, before, short)
    # the relative changes against values computed here
    shifted = Markers(before.x, before.v .+ 0.1, 2 .* before.w; m = before.m)
    c = case_conservation(case, before, shifted)
    @test isapprox(c.weight, 1.0; rtol = 1e-14)
    ke(mk) = mk.m / 2 * sum(mk.w .* mk.v .^ 2)
    @test isapprox(c.kinetic_energy, abs(ke(shifted) - ke(before)) / ke(before); rtol = 1e-12)
    p(mk) = mk.m * sum(mk.w .* mk.v)
    @test isapprox(c.momentum,
        abs(p(shifted) - p(before)) / (before.m * sum(abs.(before.w .* before.v)));
        rtol = 1e-12)
    d = case_displacement(case, before, shifted)
    @test isapprox(d.max_dv, 0.1; rtol = 1e-12)
    @test isapprox(d.max_weight_change, 1.0; rtol = 1e-12)
    @test isapprox(d.kinetic_energy_change, (ke(shifted) - ke(before)) / ke(before);
        rtol = 1e-12)
    # a position moved across the periodic edge counts its short way round
    x = copy(before.x)
    x[1] = before.x[1] + 2π - 0.05
    @test isapprox(
        case_displacement(case, before, Markers(x, before.v, before.w; m = before.m)).max_dx,
        0.05; rtol = 1e-10)
end

@testset "the exported measures" begin
    for name in (:case_evenness, :case_outside, :case_representation, :case_charge_noise,
        :case_weights, :case_conservation, :case_displacement, :case_reference)
        @test Base.isexported(ParticleRemapping, name)
    end
end
