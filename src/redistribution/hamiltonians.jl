# Translation- and permutation-invariant ensemble Hamiltonians H_θ(q, p), linear in θ, and the
# controls that show why uncoupled Hamiltonians and mean-force subtraction fail (1D-1V).

"""
    PairwiseSplineHamiltonian(Lx, p_scale; nx_cells = 8, np_cells = 6, degree = 4)

The pairwise ensemble Hamiltonian

    H_θ(q, p) = p_s Σ_a Σ_r θ⁰_r χ_r(p_a/p_s)
              + p_s/(2N) Σ_{a,b} Σ_{m,l,r} θ²_{m,l,r} φ_m(q_a − q_b) χ_l(p_a/p_s) χ_r(p_b/p_s)

with periodic B-splines `φ_m` of degree `degree` on `[0, Lx)` (`nx_cells` cells, so `M_x =
nx_cells`) and clamped B-splines `χ_r` of the same degree on `[−1, 1]` (`np_cells` cells,
`M_p = np_cells + degree`), in the scaled momentum `p/p_s`. The parameter vector holds `θ⁰`
followed by `θ²` in the column-major order of `(m, l, r)`; its length is
`M_p + M_x M_p²` ([`hamiltonian_parameter_count`](@ref)).
"""
struct PairwiseSplineHamiltonian{BX, BP, T <: Real}
    bx::BX
    bp::BP
    Lx::T
    p_scale::T
    Mx::Int
    Mp::Int
end

function PairwiseSplineHamiltonian(Lx::Real, p_scale::Real; nx_cells::Integer = 8,
        np_cells::Integer = 6, degree::Integer = 4)
    Lx, p_scale = promote(float(Lx), float(p_scale))
    bx = BSplineBasis(UniformMesh(nx_cells, (zero(Lx), Lx)), degree, Periodic())
    bp = BSplineBasis(UniformMesh(np_cells, (-one(Lx), one(Lx))), degree, Free())
    PairwiseSplineHamiltonian(bx, bp, Lx, p_scale, nbasis(bx), nbasis(bp))
end

"""
    hamiltonian_parameter_count(H)

The number `N_θ` of parameters of the linear ensemble Hamiltonian `H`.
"""
hamiltonian_parameter_count(H::PairwiseSplineHamiltonian) = H.Mp + H.Mx * H.Mp^2

# The relative position (q_a − q_b) mod Lx. The self pair is the constant 0, whose derivative
# vanishes identically; SimpleSplines 0.4 returns NaN partials for a dual number at the knot 0
# of a periodic basis, so the self pair must not carry the partials of q_a − q_a.
function _pair_difference(H::PairwiseSplineHamiltonian, q::AbstractVector, a::Integer,
        b::Integer)
    a == b ? zero(H.Lx) : mod(q[a] - q[b], H.Lx)
end

# The momentum basis χ_r^{(d)}(p_a / p_s) / p_s^d as an N × M_p matrix.
function _momentum_basis(H::PairwiseSplineHamiltonian, p::AbstractVector, d::Integer)
    [evaluate(H.bp, r, p[a] / H.p_scale, d) / H.p_scale^d
     for a in eachindex(p), r in 1:(H.Mp)]
end

"""
    basis_hamiltonians(H::PairwiseSplineHamiltonian, q, p)

The values `Ψ_k(q, p)` of the basis Hamiltonians, `H_θ = Σ_k θ_k Ψ_k`.
"""
function basis_hamiltonians(H::PairwiseSplineHamiltonian, q::AbstractVector, p::AbstractVector)
    n = length(q)
    chi0 = _momentum_basis(H, p, 0)
    vals = zeros(eltype(chi0), hamiltonian_parameter_count(H))
    for r in 1:(H.Mp)
        vals[r] = H.p_scale * sum(@view chi0[:, r])
    end
    lin_idx = LinearIndices((H.Mx, H.Mp, H.Mp))
    c_fac = H.p_scale / (2 * n)
    for a in 1:n, b in 1:n

        dq_ab = _pair_difference(H, q, a, b)
        for m in 1:(H.Mx)
            phi0 = evaluate(H.bx, m, dq_ab, 0)
            phi0 == 0 && continue
            for l in 1:(H.Mp)
                cl_a = chi0[a, l]
                cl_a == 0 && continue
                for r in 1:(H.Mp)
                    cr_b = chi0[b, r]
                    cr_b == 0 && continue
                    vals[H.Mp + lin_idx[m, l, r]] += c_fac * phi0 * cl_a * cr_b
                end
            end
        end
    end
    vals
end

"""
    ensemble_hamiltonian(H, q, p, θ)

The value `H_θ(q, p) = Σ_k θ_k Ψ_k(q, p)` of a linear ensemble Hamiltonian.
"""
function ensemble_hamiltonian(H::PairwiseSplineHamiltonian, q, p, θ::AbstractVector)
    dot(θ, basis_hamiltonians(H, q, p))
end

"""
    basis_hamiltonian_gradients(H::PairwiseSplineHamiltonian, q, p) -> (∂qΨ, ∂pΨ)

The `N × N_θ` matrices of the gradients `∂Ψ_k/∂q_a` and `∂Ψ_k/∂p_a` of the basis Hamiltonians.
By translation invariance every column of `∂qΨ` sums to zero.
"""
function basis_hamiltonian_gradients(H::PairwiseSplineHamiltonian, q::AbstractVector,
        p::AbstractVector)
    n = length(q)
    chi0 = _momentum_basis(H, p, 0)
    chi1 = _momentum_basis(H, p, 1)
    T = promote_type(eltype(chi0), eltype(q))
    dH_dq = zeros(T, n, hamiltonian_parameter_count(H))
    dH_dp = zeros(T, n, hamiltonian_parameter_count(H))
    for r in 1:(H.Mp)
        dH_dp[:, r] .= H.p_scale .* (@view chi1[:, r])
    end
    lin_idx = LinearIndices((H.Mx, H.Mp, H.Mp))
    c_fac = H.p_scale / (2 * n)
    for a in 1:n, b in 1:n

        dq_ab = _pair_difference(H, q, a, b)
        for m in 1:(H.Mx)
            phi0 = evaluate(H.bx, m, dq_ab, 0)
            phi1 = evaluate(H.bx, m, dq_ab, 1)
            (phi0 == 0 && phi1 == 0) && continue
            for l in 1:(H.Mp)
                cl_a, dcl_a = chi0[a, l], chi1[a, l]
                (cl_a == 0 && dcl_a == 0) && continue
                for r in 1:(H.Mp)
                    cr_b, dcr_b = chi0[b, r], chi1[b, r]
                    (cr_b == 0 && dcr_b == 0) && continue
                    col = H.Mp + lin_idx[m, l, r]
                    if phi1 != 0 && cl_a != 0 && cr_b != 0
                        g = c_fac * phi1 * cl_a * cr_b
                        dH_dq[a, col] += g
                        dH_dq[b, col] -= g
                    end
                    if phi0 != 0 && dcl_a != 0 && cr_b != 0
                        dH_dp[a, col] += c_fac * phi0 * dcl_a * cr_b
                    end
                    if phi0 != 0 && cl_a != 0 && dcr_b != 0
                        dH_dp[b, col] += c_fac * phi0 * cl_a * dcr_b
                    end
                end
            end
        end
    end
    dH_dq, dH_dp
end

"""
    pairwise_hamiltonian_derivatives(H::PairwiseSplineHamiltonian, q, p, θ) -> (∇qH, ∇pH, ∇²H)

The gradients `∇_q H_θ`, `∇_p H_θ` (vectors of length `N`) and the Hessian `∇²_z H_θ`
(`2N × 2N`, in the packed order `z = (q, p)`) of the pairwise Hamiltonian.
"""
function pairwise_hamiltonian_derivatives(H::PairwiseSplineHamiltonian, q::AbstractVector,
        p::AbstractVector, θ::AbstractVector)
    n = length(q)
    chi0 = _momentum_basis(H, p, 0)
    chi1 = _momentum_basis(H, p, 1)
    chi2 = _momentum_basis(H, p, 2)
    T = promote_type(eltype(chi0), eltype(q), eltype(θ))
    Hq, Hp = zeros(T, n), zeros(T, n)
    Hqq, Hqp, Hpp = zeros(T, n, n), zeros(T, n, n), zeros(T, n, n)
    for a in 1:n, r in 1:(H.Mp)

        Hp[a] += θ[r] * H.p_scale * chi1[a, r]
        Hpp[a, a] += θ[r] * H.p_scale * chi2[a, r]
    end
    Θ2 = reshape(@view(θ[(H.Mp + 1):end]), H.Mx, H.Mp, H.Mp)
    c_fac = H.p_scale / (2 * n)
    for a in 1:n, b in 1:n

        dq_ab = _pair_difference(H, q, a, b)
        for m in 1:(H.Mx)
            phi0 = evaluate(H.bx, m, dq_ab, 0)
            phi1 = evaluate(H.bx, m, dq_ab, 1)
            phi2 = evaluate(H.bx, m, dq_ab, 2)
            (phi0 == 0 && phi1 == 0 && phi2 == 0) && continue
            F00, F10, F01 = zero(T), zero(T), zero(T)
            F20, F02, F11 = zero(T), zero(T), zero(T)
            for l in 1:(H.Mp)
                c0l, c1l, c2l = chi0[a, l], chi1[a, l], chi2[a, l]
                (c0l == 0 && c1l == 0 && c2l == 0) && continue
                for r in 1:(H.Mp)
                    t_mlr = Θ2[m, l, r]
                    t_mlr == 0 && continue
                    c0r, c1r, c2r = chi0[b, r], chi1[b, r], chi2[b, r]
                    F00 += t_mlr * c0l * c0r
                    F10 += t_mlr * c1l * c0r
                    F01 += t_mlr * c0l * c1r
                    F20 += t_mlr * c2l * c0r
                    F02 += t_mlr * c0l * c2r
                    F11 += t_mlr * c1l * c1r
                end
            end
            Hq[a] += c_fac * phi1 * F00
            Hq[b] -= c_fac * phi1 * F00
            Hp[a] += c_fac * phi0 * F10
            Hp[b] += c_fac * phi0 * F01

            g_qq = c_fac * phi2 * F00
            Hqq[a, a] += g_qq
            Hqq[b, b] += g_qq
            Hqq[a, b] -= g_qq
            Hqq[b, a] -= g_qq

            Hqp[a, a] += c_fac * phi1 * F10
            Hqp[a, b] += c_fac * phi1 * F01
            Hqp[b, a] -= c_fac * phi1 * F10
            Hqp[b, b] -= c_fac * phi1 * F01

            Hpp[a, a] += c_fac * phi0 * F20
            Hpp[b, b] += c_fac * phi0 * F02
            Hpp[a, b] += c_fac * phi0 * F11
            Hpp[b, a] += c_fac * phi0 * F11
        end
    end
    Hq, Hp, [Hqq Hqp; transpose(Hqp) Hpp]
end

"""
    FourierMomentHamiltonian(wavenumbers, p_scale; np_cells = 6, degree = 4)

The pairwise ensemble Hamiltonian of [`PairwiseSplineHamiltonian`](@ref) with the Fourier basis
`φ = {cos(k Δq), sin(k Δq)}_{k ∈ K}` on the torus in place of the periodic splines. Its double
sum factors through the collective moments ([`collective_moments`](@ref)), so that
[`ensemble_hamiltonian`](@ref) and [`ensemble_hamiltonian_gradient`](@ref) cost
`O(N |K| M_p²)`. The parameter vector holds `θ⁰` (length `M_p`), then the cosine coefficients
`θᶜ_{k,l,r}` and the sine coefficients `θˢ_{k,l,r}`, each `|K| × M_p × M_p` in column-major
order.
"""
struct FourierMomentHamiltonian{BP, V <: AbstractVector, T <: Real}
    bp::BP
    wavenumbers::V
    p_scale::T
    Mp::Int
end

function FourierMomentHamiltonian(wavenumbers::AbstractVector, p_scale::Real;
        np_cells::Integer = 6, degree::Integer = 4)
    p_scale = float(p_scale)
    bp = BSplineBasis(UniformMesh(np_cells, (-one(p_scale), one(p_scale))), degree, Free())
    FourierMomentHamiltonian(bp, float.(wavenumbers), p_scale, nbasis(bp))
end

function hamiltonian_parameter_count(H::FourierMomentHamiltonian)
    H.Mp + 2 * length(H.wavenumbers) * H.Mp^2
end

function _momentum_basis(H::FourierMomentHamiltonian, p::AbstractVector, d::Integer)
    [evaluate(H.bp, r, p[a] / H.p_scale, d) / H.p_scale^d
     for a in eachindex(p), r in 1:(H.Mp)]
end

function _fourier_parameters(H::FourierMomentHamiltonian, θ::AbstractVector)
    K, Mp = length(H.wavenumbers), H.Mp
    n2 = K * Mp^2
    θ0 = @view θ[1:Mp]
    θc = reshape(@view(θ[(Mp + 1):(Mp + n2)]), K, Mp, Mp)
    θs = reshape(@view(θ[(Mp + n2 + 1):(Mp + 2n2)]), K, Mp, Mp)
    θ0, θc, θs
end

"""
    collective_moments(H::FourierMomentHamiltonian, q, p) -> (C, S)

The `|K| × M_p` collective moments `C_{k,r} = N⁻¹ Σ_b cos(k q_b) χ_r(p_b/p_s)` and
`S_{k,r} = N⁻¹ Σ_b sin(k q_b) χ_r(p_b/p_s)`.
"""
function collective_moments(H::FourierMomentHamiltonian, q::AbstractVector, p::AbstractVector)
    χ = _momentum_basis(H, p, 0)
    n = length(q)
    C = [sum(cos(k * q[b]) * χ[b, r] for b in 1:n) / n
         for k in H.wavenumbers, r in 1:(H.Mp)]
    S = [sum(sin(k * q[b]) * χ[b, r] for b in 1:n) / n
         for k in H.wavenumbers, r in 1:(H.Mp)]
    C, S
end

# The pair sum (p_s / 2N) Σ_ab Σ_lr [θᶜ cos(kΔq) + θˢ sin(kΔq)] χ_l χ_r of one wave number as a
# quadratic form in the moments, (p_s N / 2) [Cᵀθᶜ C + Sᵀθᶜ S + Sᵀθˢ C − Cᵀθˢ S].
function _fourier_pair_energy(H, N, C, S, θc, θs)
    E = zero(promote_type(eltype(C), eltype(θc)))
    for k in axes(C, 1)
        Ck, Sk = view(C, k, :), view(S, k, :)
        Θc, Θs = view(θc, k, :, :), view(θs, k, :, :)
        E += dot(Ck, Θc, Ck) + dot(Sk, Θc, Sk) + dot(Sk, Θs, Ck) - dot(Ck, Θs, Sk)
    end
    H.p_scale * N / 2 * E
end

function ensemble_hamiltonian(H::FourierMomentHamiltonian, q, p, θ::AbstractVector)
    θ0, θc, θs = _fourier_parameters(H, θ)
    χ = _momentum_basis(H, p, 0)
    C, S = collective_moments(H, q, p)
    H.p_scale * sum(χ * θ0) + _fourier_pair_energy(H, length(q), C, S, θc, θs)
end

"""
    fourier_pairwise_hamiltonian(H::FourierMomentHamiltonian, q, p, θ)

The value of the Fourier-moment Hamiltonian by the direct `O(N²)` double sum over pairs, the
reference for the collective factorisation of [`ensemble_hamiltonian`](@ref).
"""
function fourier_pairwise_hamiltonian(H::FourierMomentHamiltonian, q::AbstractVector,
        p::AbstractVector, θ::AbstractVector)
    θ0, θc, θs = _fourier_parameters(H, θ)
    χ = _momentum_basis(H, p, 0)
    n = length(q)
    E = H.p_scale * sum(χ * θ0)
    for a in 1:n, b in 1:n, (ik, k) in enumerate(H.wavenumbers)
        c, s = cos(k * (q[a] - q[b])), sin(k * (q[a] - q[b]))
        E += H.p_scale / (2n) *
             (c * dot(view(χ, a, :), view(θc, ik, :, :), view(χ, b, :)) +
              s * dot(view(χ, a, :), view(θs, ik, :, :), view(χ, b, :)))
    end
    E
end

"""
    ensemble_hamiltonian_gradient(H, q, p, θ) -> (∇qH, ∇pH)

The gradients of the ensemble Hamiltonian. For the [`FourierMomentHamiltonian`](@ref) they are
assembled through the collective moments in `O(N |K| M_p²)` operations; for the
[`CentroidFrameHamiltonian`](@ref) they are the recoil-coupled gradients.
"""
function ensemble_hamiltonian_gradient(H::FourierMomentHamiltonian, q::AbstractVector,
        p::AbstractVector, θ::AbstractVector)
    θ0, θc, θs = _fourier_parameters(H, θ)
    n = length(q)
    χ0, χ1 = _momentum_basis(H, p, 0), _momentum_basis(H, p, 1)
    C, S = collective_moments(H, q, p)
    # ∂E/∂C and ∂E/∂S of the quadratic form, one row per wave number
    dC = similar(C, promote_type(eltype(C), eltype(θ)))
    dS = similar(dC)
    for k in axes(C, 1)
        Θc, Θs = view(θc, k, :, :), view(θs, k, :, :)
        Ck, Sk = view(C, k, :), view(S, k, :)
        dC[k, :] = (Θc + transpose(Θc)) * Ck + transpose(Θs) * Sk - Θs * Sk
        dS[k, :] = (Θc + transpose(Θc)) * Sk + Θs * Ck - transpose(Θs) * Ck
    end
    dC .*= H.p_scale * n / 2
    dS .*= H.p_scale * n / 2
    Hq = zeros(eltype(dC), n)
    Hp = H.p_scale .* (χ1 * θ0)
    for a in 1:n, (ik, k) in enumerate(H.wavenumbers)

        ca, sa = cos(k * q[a]), sin(k * q[a])
        # ∂C_{k,r}/∂q_a = −k sin(k q_a) χ_r(a)/N, ∂S_{k,r}/∂q_a = k cos(k q_a) χ_r(a)/N
        Hq[a] += k / n * (-sa * dot(view(dC, ik, :), view(χ0, a, :)) +
                  ca * dot(view(dS, ik, :), view(χ0, a, :)))
        Hp[a] += (ca * dot(view(dC, ik, :), view(χ1, a, :)) +
                  sa * dot(view(dS, ik, :), view(χ1, a, :))) / n
    end
    Hq, Hp
end

"""
    CentroidFrameHamiltonian(Qmax, p_scale; nq_cells = 8, np_cells = 6, degree = 4)

The centroid-frame ensemble Hamiltonian on `ℝ`, `H_θ(q, p) = Σ_a h_θ(q_a − q̄(q), p_a)` with the
centroid `q̄ = N⁻¹ Σ_b q_b` and `h_θ(q̃, p) = p_s Σ_{i,r} θ_{i,r} B_i(q̃) χ_r(p/p_s)`, clamped
B-splines `B_i` on `[−Qmax, Qmax]` (`nq_cells` cells) and `χ_r` on `[−1, 1]` (`np_cells` cells).
The parameters `θ_{i,r}` are in column-major order.
"""
struct CentroidFrameHamiltonian{BQ, BP, T <: Real}
    bq::BQ
    bp::BP
    p_scale::T
    Mq::Int
    Mp::Int
end

function CentroidFrameHamiltonian(Qmax::Real, p_scale::Real; nq_cells::Integer = 8,
        np_cells::Integer = 6, degree::Integer = 4)
    Qmax, p_scale = promote(float(Qmax), float(p_scale))
    bq = BSplineBasis(UniformMesh(nq_cells, (-Qmax, Qmax)), degree, Free())
    bp = BSplineBasis(UniformMesh(np_cells, (-one(Qmax), one(Qmax))), degree, Free())
    CentroidFrameHamiltonian(bq, bp, p_scale, nbasis(bq), nbasis(bp))
end

hamiltonian_parameter_count(H::CentroidFrameHamiltonian) = H.Mq * H.Mp

# h_θ(q̃, p) and its partial derivatives ∂_q̃ h, ∂_p h.
function _centroid_frame_single(H::CentroidFrameHamiltonian, q̃, p, θ::AbstractVector)
    Θ = reshape(θ, H.Mq, H.Mp)
    B0 = [evaluate(H.bq, i, q̃, 0) for i in 1:(H.Mq)]
    B1 = [evaluate(H.bq, i, q̃, 1) for i in 1:(H.Mq)]
    χ0 = [evaluate(H.bp, r, p / H.p_scale, 0) for r in 1:(H.Mp)]
    χ1 = [evaluate(H.bp, r, p / H.p_scale, 1) / H.p_scale for r in 1:(H.Mp)]
    H.p_scale .* (dot(B0, Θ, χ0), dot(B1, Θ, χ0), dot(B0, Θ, χ1))
end

function ensemble_hamiltonian(H::CentroidFrameHamiltonian, q, p, θ::AbstractVector)
    q̄ = sum(q) / length(q)
    sum(_centroid_frame_single(H, q[a] - q̄, p[a], θ)[1] for a in eachindex(q, p))
end

function ensemble_hamiltonian_gradient(H::CentroidFrameHamiltonian, q::AbstractVector,
        p::AbstractVector, θ::AbstractVector)
    q̄ = sum(q) / length(q)
    parts = [_centroid_frame_single(H, q[a] - q̄, p[a], θ) for a in eachindex(q, p)]
    hq = [h[2] for h in parts]
    # ∇_{q_a} H = ∂_q̃ h(q_a − q̄, p_a) − N⁻¹ Σ_c ∂_q̃ h(q_c − q̄, p_c): the recoil term
    hq .- sum(hq) / length(hq), [h[3] for h in parts]
end

function ensemble_hamiltonian_gradient(H::PairwiseSplineHamiltonian, q::AbstractVector,
        p::AbstractVector, θ::AbstractVector)
    Hq, Hp, _ = pairwise_hamiltonian_derivatives(H, q, p, θ)
    Hq, Hp
end

"""
    uncoupled_control_step(z, p_scale)

Control A: the step `p̄_a = (p_a − 0.1 p_s cos q_a) / (1 + 0.05 cos q_a)`,
`q̄_a = q_a + 0.05 sin q_a` of an uncoupled, position-dependent single-particle map on the packed
state `z = (q, p)`. It changes the total momentum.
"""
function uncoupled_control_step(z::AbstractVector, p_scale::Real)
    n = length(z) ÷ 2
    q, p = z[1:n], z[(n + 1):(2n)]
    vcat(q .+ 0.05 .* sin.(q), (p .- 0.1 .* p_scale .* cos.(q)) ./ (1 .+ 0.05 .* cos.(q)))
end

"""
    mean_force_control_step(z, p_scale)

Control B: the step `p̄_a = p_a + F_a − N⁻¹ Σ_c F_c` with the force
`F_a = −0.1 p_s cos q_a (1 + p_a / (2 p_s))`, and `q̄_a = q_a + 0.05 sin q_a`, on the packed state
`z = (q, p)`. The subtraction of the mean force keeps the total momentum but breaks
symplecticity.
"""
function mean_force_control_step(z::AbstractVector, p_scale::Real)
    n = length(z) ÷ 2
    q, p = z[1:n], z[(n + 1):(2n)]
    f_raw = -0.1 .* p_scale .* cos.(q) .* (1 .+ 0.5 .* p ./ p_scale)
    vcat(q .+ 0.05 .* sin.(q), p .+ (f_raw .- sum(f_raw) / n))
end
