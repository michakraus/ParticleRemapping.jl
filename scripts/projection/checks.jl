# The numerical checks of the structural statements of the Hamiltonian particle projection: runs
# the test files of the compressibility, the fold and the mean-field map, of the kick–drift and
# non-separable stages, of the charge-conserving reweighting and of the low-rank solver, each in a
# module of its own, and prints the measured numbers of the table of checks in the documentation.
#
#   julia --startup-file=no --project=scripts scripts/projection/checks.jl
#
# Runtime: about one minute.

const TESTS = joinpath(@__DIR__, "..", "..", "test", "projection")

for file in ("structure.jl", "stages.jl", "charge.jl", "lowrank.jl")
    println("\n==== ", file)
    m = Module(Symbol(first(splitext(file))))
    Core.eval(m, :(include(x) = Base.include($m, x)))
    Base.include(m, joinpath(TESTS, file))
end
