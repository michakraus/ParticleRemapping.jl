# The phase-space B-spline projection: the L² projection, the positive lumped projection, the
# local sampling volumes from the number density, and the weighted KKT projection.

"""
    PhaseSpaceSplines(; Lx = 2π, Vmax = 3.5, nx = 8, nv = 10, degree = 3)

The tensor-product B-spline basis `Λ_I(x, v)` of degree `degree` on `[0, Lx) × [−Vmax, Vmax]`,
periodic in `x` with `nx` cells and clamped (open knot vector) in `v` with `nv` cells, with its
Gauss quadrature, its mass operator and the integrals `μ_I = ∫ Λ_I dx dv`.
"""
struct PhaseSpaceSplines{B, Q, M, A <: AbstractMatrix, T <: Real}
    basis::B
    quadrature::Q
    mass::M
    integrals::A
    Lx::T
    Vmax::T
end

function PhaseSpaceSplines(;
        Lx::Real = 2π, Vmax::Real = 3.5, nx::Integer = 8, nv::Integer = 10,
        degree::Integer = 3)
    Lx, Vmax = promote(float(Lx), float(Vmax))
    bx = BSplineBasis(UniformMesh(nx, (zero(Lx), Lx)), degree, Periodic())
    bv = BSplineBasis(UniformMesh(nv, (-Vmax, Vmax)), degree, Free())
    B = bx ⊗ bv
    quad = TensorProductQuadrature(B)
    M = mass_operator(quad)
    # partition of unity: μ_I = ∫ Λ_I = Σ_J M_IJ
    μ = M * ones(typeof(Lx), size(B)...)
    PhaseSpaceSplines(B, quad, M, μ, Lx, Vmax)
end

"""
    deposit_phase_space(ps::PhaseSpaceSplines, x, v, w)

The load vector `b_I = Σ_a w_a Λ_I(x_a, v_a)`, as an array of the size of the basis. Positions
are taken modulo `Lx`.
"""
function deposit_phase_space(ps::PhaseSpaceSplines, x::AbstractVector, v::AbstractVector,
        w::AbstractVector)
    B = ps.basis
    sz = size(B)
    c = zeros(promote_type(eltype(x), eltype(v), eltype(w)), sz...)
    bufs = ntuple(k -> zeros(eltype(c), local_width(bases(B)[k])), 2)
    for a in eachindex(x, v, w)
        j0 = evaluate_all!(bufs, B, (mod(x[a], ps.Lx), v[a]), (0, 0))
        for t in CartesianIndices(map(length, bufs))
            Ix = basis_index(bases(B)[1], j0[1] + t[1] - 1)
            Iv = basis_index(bases(B)[2], j0[2] + t[2] - 1)
            (1 <= Ix <= sz[1] && 1 <= Iv <= sz[2]) || continue
            c[Ix, Iv] += w[a] * bufs[1][t[1]] * bufs[2][t[2]]
        end
    end
    c
end

"""
    l2_spline_projection(ps::PhaseSpaceSplines, x, v, w)

The coefficients `c^{L²} = M⁻¹ b` of the Galerkin `L²` projection of the markers' Klimontovich
distribution onto the basis, with `b` from [`deposit_phase_space`](@ref). The projection conserves
mass and momentum but is not positive.
"""
function l2_spline_projection(ps::PhaseSpaceSplines, x, v, w)
    ps.mass \ deposit_phase_space(ps, x, v, w)
end

"""
    lumped_spline_projection(ps::PhaseSpaceSplines, x, v, w; c_min = 1e-4)

The strictly positive lumped-mass coefficients `c^pos_I = max(b_I / μ_I, c_min)`, as a vector.
"""
function lumped_spline_projection(ps::PhaseSpaceSplines, x, v, w; c_min::Real = 1e-4)
    max.(vec(deposit_phase_space(ps, x, v, w)) ./ vec(ps.integrals), c_min)
end

"""
    spline_moments(ps::PhaseSpaceSplines, c; m = 1) -> (mass, momentum, minimum)

The mass `∫ f_h`, the momentum `∫ m v f_h` by the Gauss quadrature of the basis, and the smallest
value of `f_h` at the quadrature nodes, for the spline `f_h` with coefficients `c`.
"""
function spline_moments(ps::PhaseSpaceSplines, c::AbstractArray; m::Real = 1)
    fh = Spline(ps.basis, reshape(c, size(ps.basis)))
    nodes, weights = quadrature_nodes(ps.quadrature), quadrature_weights(ps.quadrature)
    T = eltype(c)
    mass, mom, fmin = zero(T), zero(T), T(Inf)
    for ix in eachindex(nodes[1]), iv in eachindex(nodes[2])

        val = fh((nodes[1][ix], nodes[2][iv]))
        wq = weights[1][ix] * weights[2][iv]
        mass += wq * val
        mom += wq * m * nodes[2][iv] * val
        fmin = min(fmin, val)
    end
    (mass = mass, momentum = mom, minimum = fmin)
end

"""
    evaluate_phase_space_spline(ps::PhaseSpaceSplines, c, x, v)

The values `f_h(x_a, v_a)` of the spline with coefficients `c` (positions modulo `Lx`).
"""
function evaluate_phase_space_spline(
        ps::PhaseSpaceSplines, c::AbstractArray, x::AbstractVector,
        v::AbstractVector)
    fh = Spline(ps.basis, reshape(c, size(ps.basis)))
    [fh((mod(x[a], ps.Lx), v[a])) for a in eachindex(x, v)]
end

"""
    phase_space_spline_grid(ps::PhaseSpaceSplines, c; nx = 160, nv = 160)
        -> (xs, vs, F)

The spline with coefficients `c` on a uniform `nx × nv` grid of `[0, Lx] × [−Vmax, Vmax]`, with
`F[i, j] = f_h(xs[i], vs[j])`.
"""
function phase_space_spline_grid(
        ps::PhaseSpaceSplines, c::AbstractArray; nx::Integer = 160,
        nv::Integer = 160)
    fh = Spline(ps.basis, reshape(c, size(ps.basis)))
    xs = range(zero(ps.Lx), ps.Lx; length = nx)
    vs = range(-ps.Vmax, ps.Vmax; length = nv)
    xs, vs, [fh((mod(x, ps.Lx), v)) for x in xs, v in vs]
end

"""
    number_density_volumes(ps::PhaseSpaceSplines, x̄, v̄; ε_n = 1e-6)

The local sampling volumes `ΔΩ̄_a = 1 / n_h(x̄_a, v̄_a)` of the target points, with the positive
lumped-mass number density `n_h = Σ d_I Λ_I`, `d_I = max(μ_I⁻¹ Σ_a Λ_I(x̄_a, v̄_a), ε_n)`.
"""
function number_density_volumes(
        ps::PhaseSpaceSplines, x̄::AbstractVector, v̄::AbstractVector;
        ε_n::Real = 1e-6)
    n_dep = deposit_phase_space(ps, x̄, v̄, ones(eltype(x̄), length(x̄)))
    d = reshape(max.(vec(n_dep) ./ vec(ps.integrals), ε_n), size(ps.basis))
    1 ./ evaluate_phase_space_spline(ps, d, x̄, v̄)
end

"""
    spline_evaluation_matrix(ps::PhaseSpaceSplines, x̄, v̄)

The `N × |I|` matrix `E_{a,I} = Λ_I(x̄_a, v̄_a)`, columns in the linear order of the basis.
"""
function spline_evaluation_matrix(ps::PhaseSpaceSplines, x̄::AbstractVector, v̄::AbstractVector)
    [evaluate(ps.basis, j, (mod(x̄[a], ps.Lx), v̄[a]))
     for a in eachindex(x̄, v̄),
    j in 1:nbasis(ps.basis)]
end

"""
    kkt_moment_matrix(E, ΔΩ, v̄; m = 1)

The `2 × |I|` discrete moment matrix `A = [ΔΩᵀ E; m (ΔΩ ∘ v̄)ᵀ E]` of the mass and momentum
of the weights `w̄_a = ΔΩ_a f_h(x̄_a, v̄_a)` (1D-1V).
"""
function kkt_moment_matrix(E::AbstractMatrix, ΔΩ::AbstractVector, v̄::AbstractVector; m::Real = 1)
    vcat(transpose(ΔΩ) * E, transpose(ΔΩ .* m .* v̄) * E)
end

"""
    kkt_spline_projection(c_pos, A, moments)

The weighted KKT projection
`c^KKT = c^pos − W⁻¹ Aᵀ (A W⁻¹ Aᵀ)⁻¹ (A c^pos − moments)` with `W = diag(c^pos)⁻²`: the
minimiser of `½ (c − c^pos)ᵀ W (c − c^pos)` subject to `A c = moments`. With the
[`kkt_moment_matrix`](@ref) `A` and `moments = [Σ w_a, Σ p_a]` the weights
`ΔΩ̄_a f_h^KKT(x̄_a, v̄_a)` conserve the discrete mass and momentum.
"""
function kkt_spline_projection(c_pos::AbstractVector, A::AbstractMatrix, moments::AbstractVector)
    W⁻¹ = Diagonal(c_pos .^ 2)
    c_pos - W⁻¹ * transpose(A) * ((A * W⁻¹ * transpose(A)) \ (A * c_pos - moments))
end

"""
    recover_velocity_branch(ps::PhaseSpaceSplines, c, q̂, p̂, v_start, ΔΩ; m = 1, tol = 1e-14,
        maxiter = 50) -> (v, w, slope)

The velocities `v̂_a` on the branch of `m ΔΩ_a f_h(q̂_a, v) v = p̂_a` through `v_start` (the
coarse-basis regime of the redistribution, where the step reaches `(q̂, p̂)` near the targets
only), by Newton's method in `v` from `v_start_a`, the weights `ŵ_a = ΔΩ_a f_h(q̂_a, v̂_a)`, and
the slopes `m ΔΩ_a (f_h + v̂_a ∂_v f_h)` at the solutions; a slope of the sign it has at
`v_start` keeps the solution on the branch of the start, off the fold.
"""
function recover_velocity_branch(
        ps::PhaseSpaceSplines, c::AbstractArray, q̂::AbstractVector,
        p̂::AbstractVector, v_start::AbstractVector, ΔΩ::AbstractVector; m::Real = 1,
        tol::Real = 1e-14, maxiter::Integer = 50)
    fh = Spline(ps.basis, reshape(c, size(ps.basis)))
    v̂ = float.(copy(v_start))
    slope = similar(v̂)
    for a in eachindex(v̂)
        x = mod(q̂[a], ps.Lx)
        for _ in 1:maxiter
            f, fv = fh((x, v̂[a])), fh((x, v̂[a]), (0, 1))
            g = m * ΔΩ[a] * f * v̂[a] - p̂[a]
            slope[a] = m * ΔΩ[a] * (f + v̂[a] * fv)
            abs(g) <= tol * max(abs(p̂[a]), eps()) && break
            v̂[a] -= g / slope[a]
        end
    end
    (v = v̂, w = ΔΩ .* evaluate_phase_space_spline(ps, c, q̂, v̂), slope = slope)
end

"""
    kkt_redistribution_weights(ps, c_pos, x̄, v̄, ΔΩ, mass, momentum; m = 1)
        -> (weights, coefficients, raw_weights)

The target weights `w̄_a = ΔΩ_a f_h^KKT(x̄_a, v̄_a)` of the [`kkt_spline_projection`](@ref) of the
prior `c_pos` at the target points `(x̄, v̄)` with volumes `ΔΩ`, its coefficients, and the raw
weights `ΔΩ_a f_h^pos(x̄_a, v̄_a)` of the prior, which do not conserve the moments.
"""
function kkt_redistribution_weights(ps::PhaseSpaceSplines, c_pos::AbstractVector, x̄, v̄, ΔΩ,
        mass::Real, momentum::Real; m::Real = 1)
    E = spline_evaluation_matrix(ps, x̄, v̄)
    c_kkt = kkt_spline_projection(c_pos, kkt_moment_matrix(E, ΔΩ, v̄; m = m), [
        mass, momentum])
    w̄ = ΔΩ .* evaluate_phase_space_spline(ps, c_kkt, x̄, v̄)
    (weights = w̄, coefficients = c_kkt, raw_weights = ΔΩ .* (E * c_pos))
end
