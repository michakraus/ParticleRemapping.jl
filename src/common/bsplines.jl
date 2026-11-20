"""
    bspline3(t)

The centred uniform cubic B-spline with support `(−2, 2)`, normalised so that its integer
translates form a partition of unity.
"""
function bspline3(t)
    a = abs(t)
    a < 1 && return 2 / 3 - a^2 + a^3 / 2
    a < 2 && return (2 - a)^3 / 6
    zero(t)
end

"""
    dbspline3(t)

The derivative of [`bspline3`](@ref).
"""
function dbspline3(t)
    a = abs(t)
    s = sign(t)
    a < 1 && return s * (-2a + 3a^2 / 2)
    a < 2 && return -s * (2 - a)^2 / 2
    zero(t)
end

"""
    UniformBasis(s0, h, n, periodic)

`n` cubic B-splines `B_j(s) = bspline3((s − s0)/h − (j − 1))`, centred at `s0 + (j − 1) h`,
`j = 1, …, n`. A periodic basis has period `n h` and is a partition of unity; a non-periodic one
keeps only the `n` splines whose centres lie in `[s0, s0 + (n − 1) h]`.
"""
struct UniformBasis{T <: Real}
    s0::T
    h::T
    n::Int
    periodic::Bool
end

function UniformBasis(s0::Real, h::Real, n::Integer, periodic::Bool)
    UniformBasis(promote(float(s0), float(h))..., Int(n), periodic)
end

Base.length(b::UniformBasis) = b.n

"""
    eachnonzero(f, b::UniformBasis, s)

Call `f(j, B_j(s), B_j'(s))` for the (at most four) basis functions that do not vanish at `s`.
The argument `s` may be a dual number; only its value selects the splines.
"""
function eachnonzero(f, b::UniformBasis, s)
    t = (s - b.s0) / b.h
    j0 = floor(Int, ForwardDiff.value(t))
    for j in (j0 - 1):(j0 + 2)
        τ = t - j
        if b.periodic
            jj = mod(j, b.n) + 1
        else
            (0 <= j < b.n) || continue
            jj = j + 1
        end
        f(jj, bspline3(τ), dbspline3(τ) / b.h)
    end
    nothing
end

"""
    basis_values(b::UniformBasis, s)

The vector `(B_1(s), …, B_n(s))`.
"""
function basis_values(b::UniformBasis, s)
    B = zeros(typeof(float(s)), b.n)
    eachnonzero((j, Bj, _) -> (B[j] += Bj), b, s)
    B
end

export bspline3, dbspline3, UniformBasis, eachnonzero, basis_values
