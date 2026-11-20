# The checks of the symplecticity condition for maps of one weighted marker (x, v, w): the gauge
# step with a state-dependent stretch, a canonical map followed by a free weight, free streaming,
# and the family X = φ(x), V = c(x) v + b(x), W = w/(φ'c), with two controls that must fail.
#
# Run, from the package root:
#   julia --startup-file=no --project=scripts scripts/remapping/symplectic_condition.jl

using ParticleRemapping
using Printf

checks = symplectic_condition_checks()
ok = map(checks) do c
    # central differences with h = 1e-5 leave residuals below 1e-8 for a symplectic map
    good = c.holds ? c.residual < 1e-8 : c.residual > 1e-2
    println(rpad(c.name, 54), rpad(c.holds ? "holds" : "fails", 7),
        @sprintf("max residual %.2e", c.residual), good ? "  ok" : "  UNEXPECTED")
    good
end
println("\n", count(ok), "/", length(ok), " checks as expected")
all(ok) || error("$(count(!, ok)) check(s) not as expected")
