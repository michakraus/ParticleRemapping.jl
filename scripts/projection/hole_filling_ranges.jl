# The 14 cells of the hole-filling tables that the documentation quotes as a range `lo–hi`, rebuilt:
# each cell is the range of its value over the experiment and the experiment with the targets moved
# by 1e-13 relative with Xoshiro(1:10), with each end rounded outward to the digits of the quote.
# The script prints each run's value and the rebuilt range, and checks that the rebuilt range lies
# within the quoted one. A quoted range also holds the values of the other Julia and of the CI
# runners, so the range rebuilt on one Julia can be narrower than the quote. Run it on Julia 1.13
# and on Julia 1.14:
#
#   julia --startup-file=no --project=scripts scripts/projection/hole_filling_ranges.jl
#   julia +1.14 --startup-file=no --project=scripts scripts/projection/hole_filling_ranges.jl
#
# Runtime: about 85 minutes on each Julia.

using ParticleRemapping
using Printf
using Random
using Test

# The cells quoted as a range: the hole position, the variant, the column of the tables (1 empty
# disk, 2 fill, 3 objective, 4 charge, 5 disk, 6 fill and 7 max |ln| with the charge-conserving
# weights) and the quote.
const RANGES = (
    (π / 2, :weight_dependent, 7, "0.175–0.176"),
    (π / 2, :control_A, 3, "0.834–0.835"),
    (π / 2, :control_B, 7, "0.061–0.062"),
    (π / 2, :control_C, 2, "0.54–0.55"),
    (π / 2, :control_C, 3, "0.667–0.668"),
    (π / 2, :control_C, 5, "1.85–1.86"),
    (π / 2, :control_C, 6, "0.54–0.55"),
    (Float64(π), :non_separable, 1, "1.51–1.52"),
    (Float64(π), :non_separable, 2, "0.93–0.94"),
    (Float64(π), :non_separable, 5, "1.55–1.56"),
    (Float64(π), :control_C, 2, "0.54–0.55"),
    (Float64(π), :control_C, 3, "0.713–0.714"),
    (Float64(π), :control_C, 5, "2.05–2.06"),
    (Float64(π), :control_C, 6, "0.54–0.55"))

const MOVES = 1:10

# The values of the documentation's tables of a variant's entry `r` of `results`.
function row_values(r)
    (r.diagnostics.empty, r.diagnostics.fill, r.ratio, r.charge_res,
        r.diagnostics_charge.empty, r.diagnostics_charge.fill, r.dlnw)
end

# The ends of the quote `lo–hi` in units of its last digit, and that unit's number of decimals.
function quoted_units(quoted)
    lo, hi = split(quoted, '–')
    d = length(split(hi, '.')[2])
    round(Int, parse(Float64, lo) * 10^d), round(Int, parse(Float64, hi) * 10^d), d
end

# The values of the cells of `hole_x` in the experiment with the targets `targets`.
function cell_values(hole_x, targets)
    cells = filter(c -> c[1] == hole_x, RANGES)
    variants = Tuple(unique(c[2] for c in cells))
    Random.seed!(20261119)
    e = hole_filling_experiment(; hole_x, variants, targets)
    e.targets, [row_values(e.results[v])[col] for (_, v, col, _) in cells]
end

@printf("Julia %s, check-bounds %d\n", VERSION, Base.JLOptions().check_bounds)
@testset "the cells quoted as a range" begin
    for hole_x in (π / 2, Float64(π))
        cells = filter(c -> c[1] == hole_x, RANGES)
        tg, values = cell_values(hole_x, nothing)
        runs = [values]
        for k in MOVES
            x = tg.x .* (1 .+ 1e-13 .* randn(Xoshiro(k), length(tg.x)))
            push!(runs, last(cell_values(hole_x, merge(tg, (; x)))))
        end
        @printf("\nx_h = %.4f\n", hole_x)
        for (i, (_, v, col, quoted)) in enumerate(cells)
            xs = [run[i] for run in runs]
            qlo, qhi, d = quoted_units(quoted)
            lo, hi = floor(Int, minimum(xs) * 10^d), ceil(Int, maximum(xs) * 10^d)
            @printf("%-17s column %d: quoted %-12s rebuilt %.*f–%.*f; moves 0:%d %s\n", v,
                col, quoted, d, lo / 10^d, d, hi / 10^d,
                last(MOVES),
                join([@sprintf("%.6f", x) for x in xs], " "))
            @test qlo <= lo && hi <= qhi
        end
    end
end
