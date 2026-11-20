# The whole redistribution as a method of the common interface, and the 1D-1V set-up of the
# verification: a perturbed Maxwellian on a lattice, clustered by a leapfrog run.

"""
    redistribution_example(; nx0 = 4, nv0 = 4, Lx = 2π, vlims = (-2.4, 2.4), m = 1.5,
        amplitude = 0.45, dt = 0.05, nsteps = 60,
        force = x -> 1.4 sin(x − 0.5) + 0.8 sin(2x))

The verification set-up: `nx0 × nv0` markers at the cell centres of a lattice on
`[0, Lx) × [vlims[1], vlims[2]]`, weighted by `w = Δx Δv f₀(x, v)` with
`f₀ = (1 + amplitude cos x) exp(−v²/2) / (Lx √(2π))`, then advanced by `nsteps` leapfrog steps
of length `dt` in the static force field `force` so that they cluster. Returns the initial and
the clustered [`Markers`](@ref) and the kernel widths `σx = Lx/nx0`, `σv = Δv`, which match the
initial spacing.
"""
function redistribution_example(; nx0::Integer = 4, nv0::Integer = 4, Lx::Real = 2π,
        vlims = (-2.4, 2.4), m::Real = 1.5, amplitude::Real = 0.45, dt::Real = 0.05,
        nsteps::Integer = 60, force = x -> 1.4 * sin(x - 0.5) + 0.8 * sin(2x))
    Δx = Lx / nx0
    Δv = (vlims[2] - vlims[1]) / nv0
    f0(x, v) = (1 + amplitude * cos(x)) * exp(-0.5 * v^2) / (Lx * sqrt(2π))
    x0, v0, w0 = Float64[], Float64[], Float64[]
    for i in 1:nx0, j in 1:nv0

        xi = (i - 0.5) * Δx
        vj = vlims[1] + (j - 0.5) * Δv
        push!(x0, xi)
        push!(v0, vj)
        push!(w0, Δx * Δv * f0(xi, vj))
    end
    x, v = copy(x0), copy(v0)
    for _ in 1:nsteps
        v .+= 0.5 * dt .* force.(x) ./ m
        x .= mod.(x .+ dt .* v, Lx)
        v .+= 0.5 * dt .* force.(x) ./ m
    end
    (initial = Markers(x0, v0, copy(w0); m = m), markers = Markers(x, v, w0; m = m),
        Lx = float(Lx), σx = Δx, σv = Δv)
end

"""
    SplineRedistribution(; Lx = 2π, Vmax = 3.5, nx = 8, nv = 10, degree = 3, c_min = 1e-4,
        ε_n = 1e-6, σx = π/2, σv = 1.2, relax_step = 0.006, relax_iterations = 10, v_wall = 2,
        k_wall = 2.5, nx_cells = 10, np_cells = 8, p_scale_factor = 2.5,
        scheme = ImplicitMidpointStep(), substeps = 1, τ = 1)

The symplectic redistribution of 1D-1V markers by a positive, moment-compatible spline
projection and a pairwise spline Hamiltonian:

1. project the markers onto the [`PhaseSpaceSplines`](@ref) (`nx × nv` cells of degree
   `degree`) by the positive lumped projection (`c_min`);
2. relax the points to targets `(x̄, v̄)` ([`declustering_relaxation`](@ref) with `σx`, `σv`,
   `relax_step`, `relax_iterations`, `v_wall`, `k_wall`);
3. weigh the targets by the [`kkt_redistribution_weights`](@ref) with the
   [`number_density_volumes`](@ref) (`ε_n`), so that `p̄ = m w̄ v̄` conserves the discrete mass and
   momentum;
4. solve for the parameters of a [`PairwiseSplineHamiltonian`](@ref) (`nx_cells`, `np_cells`,
   `p_s = p_scale_factor max |p|, |p̄|`) and advance by the symplectic `scheme` in `substeps`
   steps of artificial step `τ` ([`redistribution_substeps`](@ref));
5. assign the weights `w̄` and the velocities `v = p / (m w̄)` to the end point.
"""
struct SplineRedistribution{T <: Real, S <: RedistributionScheme} <: AbstractRemap
    Lx::T
    Vmax::T
    nx::Int
    nv::Int
    degree::Int
    c_min::T
    ε_n::T
    σx::T
    σv::T
    relax_step::T
    relax_iterations::Int
    v_wall::T
    k_wall::T
    nx_cells::Int
    np_cells::Int
    p_scale_factor::T
    scheme::S
    substeps::Int
    τ::T
end

function SplineRedistribution(; Lx::Real = 2π, Vmax::Real = 3.5, nx::Integer = 8,
        nv::Integer = 10, degree::Integer = 3, c_min::Real = 1e-4, ε_n::Real = 1e-6,
        σx::Real = Lx / 4, σv::Real = 1.2, relax_step::Real = 0.006,
        relax_iterations::Integer = 10, v_wall::Real = 2, k_wall::Real = 2.5,
        nx_cells::Integer = 10, np_cells::Integer = 8, p_scale_factor::Real = 2.5,
        scheme::RedistributionScheme = ImplicitMidpointStep(), substeps::Integer = 1,
        τ::Real = 1)
    T = float(promote_type(typeof.((
        Lx, Vmax, c_min, ε_n, σx, σv, relax_step, v_wall, k_wall,
        p_scale_factor, τ))...))
    SplineRedistribution(T(Lx), T(Vmax), Int(nx), Int(nv), Int(degree), T(c_min), T(ε_n),
        T(σx), T(σv), T(relax_step), Int(relax_iterations), T(v_wall), T(k_wall),
        Int(nx_cells), Int(np_cells), T(p_scale_factor), scheme, Int(substeps), T(τ))
end

"""
    remap(method::SplineRedistribution, markers::Markers) -> RemapResult

Redistribute 1D-1V `markers` by the steps of [`SplineRedistribution`](@ref). The `info` holds
the targets `x̄`, `v̄`, the target weights `w̄`, the volumes `ΔΩ̄`, the raw weights of the prior,
the KKT coefficients, the momentum scale, the parameters, contraction numbers, residuals and
Newton histories of the steps, the scaled distance `target_error` of the end point from the
targets, and the kernel-density coefficients of variation before and after.
"""
function remap(method::SplineRedistribution, markers::Markers)
    markers.x isa AbstractVector ||
        throw(ArgumentError("SplineRedistribution redistributes 1D-1V markers"))
    x, v, w, m = markers.x, markers.v, markers.w, markers.m
    p = w .* m .* v
    ps = PhaseSpaceSplines(; Lx = method.Lx, Vmax = method.Vmax, nx = method.nx,
        nv = method.nv, degree = method.degree)
    c_pos = lumped_spline_projection(ps, x, v, w; c_min = method.c_min)
    relaxed = declustering_relaxation(x, v; Lx = method.Lx, σx = method.σx, σv = method.σv,
        step = method.relax_step, iterations = method.relax_iterations,
        v_wall = method.v_wall, k_wall = method.k_wall)
    x̄, v̄ = relaxed.x, relaxed.v
    ΔΩ = number_density_volumes(ps, x̄, v̄; ε_n = method.ε_n)
    kkt = kkt_redistribution_weights(ps, c_pos, x̄, v̄, ΔΩ, sum(w), sum(p); m = m)
    w̄ = kkt.weights
    p̄ = w̄ .* m .* v̄
    p_scale = method.p_scale_factor * maximum(max.(abs.(p), abs.(p̄)))
    H = PairwiseSplineHamiltonian(method.Lx, p_scale; nx_cells = method.nx_cells,
        np_cells = method.np_cells)
    steps = redistribution_substeps(
        H, method.scheme, x, p, x̄, p̄; substeps = method.substeps,
        τ = method.τ)
    N = length(x)
    q̂, p̂ = steps.z[1:N], steps.z[(N + 1):(2N)]
    S = redistribution_scaling(N, p_scale)
    cv(xs, vs) = kernel_density_cv(mod.(xs, method.Lx), vs; Lx = method.Lx, σx = method.σx,
        σv = method.σv)
    after = Markers(mod.(q̂, method.Lx), p̂ ./ (m .* w̄), w̄; m = m)
    info = (targets_x = x̄, targets_v = v̄, target_weights = w̄, volumes = ΔΩ,
        raw_weights = kkt.raw_weights, coefficients = kkt.coefficients, p_scale = p_scale,
        hamiltonian = H, θ = steps.θ, contraction = steps.contraction,
        residual = steps.residual, history = steps.history,
        target_error = norm(S .* (steps.z .- vcat(x̄, p̄))),
        cv_before = cv(x, v), cv_after = cv(x̄, v̄), relaxation = relaxed.history)
    RemapResult(after, info)
end
