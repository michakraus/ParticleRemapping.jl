# The worked examples of the canonical transport: the admissible cotangent lift of a separable
# density in closed form, and the current-quantile matching of five markers in one cell.

"""
    cosine_profile(x, α = 0.4)

The spatial profile `J(x) = 1 + α cos x` of the separable density `f = J(x) M(v)`.
"""
cosine_profile(x::Real, α::Real = 0.4) = 1 + α * cos(x)

"""
    lift_theta(x; α = 0.4, n = 2000)

`Θ(x) = ∫_0^x J(y)^{1/2} dy` for `J` = [`cosine_profile`](@ref), by Simpson's rule with `n`
subintervals. The flow of `ẋ = J(x)^{−1/2}` is `Θ(φ_t(x)) = Θ(x) + t`.
"""
function lift_theta(x::Real; α::Real = 0.4, n::Integer = 2000)
    h = x / n
    s = sqrt(cosine_profile(0.0, α)) + sqrt(cosine_profile(x, α))
    for k in 1:(n - 1)
        s += (isodd(k) ? 4 : 2) * sqrt(cosine_profile(k * h, α))
    end
    s * h / 3
end

"""
    lift_theta_inverse(y; α = 0.4)

The inverse of [`lift_theta`](@ref), by Newton's method from the linear interpolant.
"""
function lift_theta_inverse(y::Real; α::Real = 0.4)
    x = y / lift_theta(2π; α) * 2π
    for _ in 1:50
        δ = (lift_theta(x; α) - y) / sqrt(cosine_profile(x, α))
        x -= δ
        abs(δ) < 1e-14 && break
    end
    x
end

"""
    admissible_lift(x, t; α = 0.4) -> (X, φ_t'(x))

The admissible cotangent lift `a = J^{−1/2}` in closed form: `X = φ_t(x) = Θ⁻¹(Θ(x) + t)` (not
reduced modulo `2π`) and `φ_t'(x) = (J(x)/J(X))^{1/2}` ([`lift_theta`](@ref)).
"""
function admissible_lift(x::Real, t::Real; α::Real = 0.4)
    X = lift_theta_inverse(lift_theta(x; α) + t; α)
    X, sqrt(cosine_profile(x, α) / cosine_profile(X, α))
end

"""
    admissible_lift_matching(x, v, t, w0; α = 0.4, m = 1)

The continuum current-quantile matching after the admissible lift ([`admissible_lift`](@ref)) of
a marker at `(x, v)` with weight `w0`, for `f = J(x) M(v)` and equal weights: the cumulative
momentum `Φ_X(P) = (J(x)/φ'²) G(|P| φ'/(m w0))` matched to `F_X(V) = J(X) G(V)`. Returns
`(; X, P, V, W, dφ)`; the matching gives `V = v` and `W = w0/φ'`.
"""
function admissible_lift_matching(x::Real, v::Real, t::Real, w0::Real; α::Real = 0.4,
        m::Real = 1)
    X, dφ = admissible_lift(x, t; α)
    P = m * w0 * v / dφ
    Φ = cosine_profile(x, α) / dφ^2 * maxwellian_current(abs(P) * dφ / (m * w0))
    V = sign(P) * maxwellian_current_inverse(Φ / cosine_profile(X, α))
    (; X, P, V, W = P / (m * V), dφ)
end

"""
    quantile_matching_example(; v = [0.2, 0.5, 0.75, 1.1, 1.45], w0 = 0.1, Δx = 1, Jbar = 1,
                              m = 1)

The worked example of the current-quantile matching: one cell of width `Δx`, `f̂ = J̄ M(v)`, and
markers on the positive branch with equal weights `w0` and the velocities `v`, as a canonical map
has left them. Returns `(; v, w0, Δx, Jbar, m, aP)` and the fields `Π, R, Πc, V, W` of
[`branch_quantile_matching`](@ref).
"""
function quantile_matching_example(; v::AbstractVector = [0.2, 0.5, 0.75, 1.1, 1.45],
        w0::Real = 0.1, Δx::Real = 1.0, Jbar::Real = 1.0, m::Real = 1.0)
    all(>(0), v) || throw(ArgumentError("the example lies on the positive branch"))
    vs = sort(v)
    aP = m * w0 .* vs
    b = branch_quantile_matching(aP, Jbar * maxwellian_current(Inf),
        y -> maxwellian_current_inverse(y / Jbar); m, Δx)
    (; v = vs, w0, Δx, Jbar, m, aP, b...)
end
