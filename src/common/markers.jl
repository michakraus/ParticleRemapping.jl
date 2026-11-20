"""
    Markers(x, v, w; m = 1)

An ensemble of `N` weighted markers of mass `m`.

In one dimension `x` and `v` are vectors of length `N`; in `d` dimensions they are `d × N`
matrices with one column per marker. The weights `w` are a vector of length `N`. The canonical
momentum of marker `a` is `p_a = m w_a v_a` ([`momentum`](@ref)).
"""
struct Markers{T <: Real, A <: AbstractVecOrMat{T}, W <: AbstractVector{T}}
    x::A
    v::A
    w::W
    m::T

    function Markers(x::A, v::A, w::W, m::T) where {T <: Real, A <: AbstractVecOrMat{T},
            W <: AbstractVector{T}}
        size(x) == size(v) ||
            throw(DimensionMismatch("x has size $(size(x)), v has size $(size(v))"))
        size(x, ndims(x)) == length(w) ||
            throw(DimensionMismatch("$(size(x, ndims(x))) markers but $(length(w)) weights"))
        new{T, A, W}(x, v, w, m)
    end
end

function Markers(x::AbstractVecOrMat, v::AbstractVecOrMat, w::AbstractVector; m::Real = 1)
    T = promote_type(eltype(x), eltype(v), eltype(w), typeof(m))
    Markers(T.(x), T.(v), T.(w), T(m))
end

Base.length(mk::Markers) = length(mk.w)
Base.copy(mk::Markers) = Markers(copy(mk.x), copy(mk.v), copy(mk.w), mk.m)

"""
    spatial_dimension(mk::Markers)

The number `d` of position (and velocity) components of each marker.
"""
spatial_dimension(mk::Markers) = mk.x isa AbstractVector ? 1 : size(mk.x, 1)

# A weight vector broadcast against `x` or `v`: itself for vectors, a row for `d × N` matrices.
_rowwise(w::AbstractVector, ::AbstractVector) = w
_rowwise(w::AbstractVector, ::AbstractMatrix) = reshape(w, 1, :)

"""
    momentum(mk::Markers)
    momentum(v, w, m)

The canonical momenta `p_a = m w_a v_a`, with the shape of `v`.
"""
momentum(v::AbstractVecOrMat, w::AbstractVector, m::Real) = m .* _rowwise(w, v) .* v
momentum(mk::Markers) = momentum(mk.v, mk.w, mk.m)

"""
    velocity(p, w, m)

The velocities `v_a = p_a / (m w_a)` of canonical momenta `p` at weights `w`.
"""
velocity(p::AbstractVecOrMat, w::AbstractVector, m::Real) = p ./ (m .* _rowwise(w, p))

"""
    canonical(mk::Markers) -> (x, p)

The canonical coordinates of the markers: positions and momenta `p = m w v`.
"""
canonical(mk::Markers) = (mk.x, momentum(mk))

"""
    from_canonical(x, p, w; m = 1)

The markers with positions `x`, weights `w` and the velocities `v = p / (m w)` that give them the
canonical momenta `p`.
"""
function from_canonical(x::AbstractVecOrMat, p::AbstractVecOrMat, w::AbstractVector; m::Real = 1)
    Markers(x, velocity(p, w, m), w; m = m)
end

"""
    gauge_step(mk::Markers, λ)

The gauge step `(x, v, w) ↦ (x, λ v, w / λ)` with positive `λ` (a number or one value per marker).
It keeps every position and every canonical momentum, so it is symplectic for any `λ`, also one
that depends on the whole ensemble.
"""
function gauge_step(mk::Markers, λ)
    all(>(0), λ) || throw(DomainError(λ, "the gauge stretch must be positive"))
    Markers(
        copy(mk.x), _rowwise(λ isa AbstractVector ? λ : fill(λ, length(mk)), mk.v) .* mk.v,
        mk.w ./ λ; m = mk.m)
end

"""
    total_weight(mk::Markers)

The sum of the weights, `Σ_a w_a`.
"""
total_weight(mk::Markers) = sum(mk.w)

"""
    total_momentum(mk::Markers)

The total canonical momentum `P = Σ_a m w_a v_a` (a number in 1D, a vector of length `d` otherwise).
"""
total_momentum(mk::Markers) = _sum_markers(momentum(mk))

"""
    kinetic_energy(mk::Markers)

The kinetic energy `(m/2) Σ_a w_a |v_a|²`.
"""
kinetic_energy(mk::Markers) = mk.m / 2 * sum(_rowwise(mk.w, mk.v) .* abs2.(mk.v))

_sum_markers(p::AbstractVector) = sum(p)
_sum_markers(p::AbstractMatrix) = vec(sum(p; dims = 2))

export Markers, spatial_dimension, momentum, velocity, canonical, from_canonical, gauge_step
export total_weight, total_momentum, kinetic_energy
