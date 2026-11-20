# Checks of the symplecticity condition for maps of one weighted marker (x, v, w) in 1D-1V.

"""
    frozen_weight_residual(T, Wfun, x, v, w; m = 1)

The canonicity residual of a map `(x, v) ↦ (X, V) = T(x, v)` with the new weight
`W = Wfun(x, v, w)` for a frozen weight `w`: `|det ∂(X, P)/∂(x, v) − m w|` with `P = m W V`, the
difference from `det ∂(x, p)/∂(x, v) = m w`. It vanishes when the map is canonical in `(x, p)`.
The Jacobian is taken by central differences ([`fd_jacobian`](@ref)).
"""
function frozen_weight_residual(T, Wfun, x::Real, v::Real, w::Real; m::Real = 1)
    function XP(y)
        X, V = T(y[1], y[2])
        [X, Wfun(y[1], y[2], w) * m * V]
    end
    J = fd_jacobian(XP, [x, v])
    abs((J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1]) - w * m)
end

"""
    symplectic_condition_checks(; m = 1.7, samples = …)

The checks of the symplecticity condition for maps of one marker `(x, v, w)`, each the largest
residual over the sample states `samples`, as `(; name, residual, holds)` with `holds` the
expected outcome:

1. the gauge step `(x, λ v, w/λ)` with a stretch `λ(x, v, w)` that depends on the whole state
   preserves `ω` ([`pullback_residual`](@ref));
2. a canonical map of `(x, p)` (a kick, then a drift) followed by an arbitrary new weight `W`,
   `V = P/(m W)`, preserves `ω`;
3. control: `X = x + v`, which does not factor through `(x, p)`, does not;
4. free streaming `X = x + 0.9 v` does not preserve `ω`;
5. free streaming is canonical in `(x, p)` for a frozen weight ([`frozen_weight_residual`](@ref));
6. `X = φ(x)`, `V = c(x) v + b(x)`, `W = w/(φ' c)` is canonical for a frozen weight;
7. control: `V = v + 0.4 v³` with `W = w / det DT` is not.
"""
function symplectic_condition_checks(; m::Real = 1.7,
        samples = [[0.3, -1.2, 0.8], [1.1, 0.4, 2.5], [-0.7, 2.0, 0.15], [2.2, -0.3, 1.3]])
    pull(Φ) = maximum(pullback_residual(Φ, z; m) for z in samples)
    frozen(T, Wfun) = maximum(frozen_weight_residual(T, Wfun, z...; m) for z in samples)

    λgauge(z) = 1 + z[1]^2 + sin(z[2]) * z[3]^2 + 0.3z[3]
    gauge(z) = [z[1], λgauge(z) * z[2], z[3] / λgauge(z)]
    function kick_drift_reweight(z)
        x, v, w = z
        p = w * m * v
        P = p - sin(x)
        X = x + P^3 / 5
        W = 0.5 + x^2 + w * cos(v)^2
        [X, P / (m * W), W]
    end
    φ(x) = x + 0.3sin(x)
    dφ(x) = 1 + 0.3cos(x)
    c(x) = 1.5 + 0.5cos(2x)
    b(x) = 0.4sin(x)
    [(; name = "gauge step with λ(x, v, w)", residual = pull(gauge), holds = true),
        (; name = "canonical ψ(x, p) and arbitrary W, V = P/(mW)",
            residual = pull(kick_drift_reweight), holds = true),
        (; name = "control: X = x + v", residual = pull(z -> [z[1] + z[2], z[2], z[3]]),
            holds = false),
        (; name = "free streaming on (x, v, w)",
            residual = pull(z -> [z[1] + 0.9z[2], z[2], z[3]]), holds = false),
        (; name = "free streaming, frozen w, in (x, p)",
            residual = frozen((x, v) -> (x + 0.9v, v), (x, v, w) -> w), holds = true),
        (; name = "X = φ(x), V = c(x) v + b(x), W = w/(φ'c), frozen w",
            residual = frozen((x, v) -> (φ(x), c(x) * v + b(x)),
                (x, v, w) -> w / (dφ(x) * c(x))), holds = true),
        (; name = "control: V = v + 0.4 v³, W = w/det DT",
            residual = frozen((x, v) -> (φ(x), v + 0.4v^3),
                (x, v, w) -> w / (dφ(x) * (1 + 1.2v^2))), holds = false)]
end
