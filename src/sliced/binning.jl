# Slices and quantile targets of the sliced lifts.
#
# A slice (bin) is a set of markers with neighbouring values of one coordinate, ordered by another
# coordinate. The lifts move the ordered markers of a slice towards equispaced targets.

"""
    uniform_quantile_targets(L, n)

The `n` midpoints `τ_r = −L + (2r − 1) L / n`, `r = 1, …, n`, of the `n` equal cells of
`[−L, L]`: the quantile targets of a uniform distribution on `[−L, L]`. They also serve as the
bin centres of `n` equal bins.
"""
uniform_quantile_targets(L::Real, n::Integer) = [-L + (2r - 1) * L / n for r in 1:n]

"""
    sorted_slices(slice_key, order_key, B)

Split the markers into `B` slices of `N / B` markers each by the values of `slice_key` (the
markers with the smallest `N / B` values form slice 1, and so on) and order each slice by
`order_key`. Returns a vector of `B` index vectors. `N` must be a multiple of `B`.
"""
function sorted_slices(slice_key::AbstractVector, order_key::AbstractVector, B::Integer)
    N = length(slice_key)
    N % B == 0 || throw(ArgumentError("$N markers do not split into $B equal slices"))
    length(order_key) == N ||
        throw(DimensionMismatch("$N slice keys but $(length(order_key)) order keys"))
    n = N ÷ B
    perm = sortperm(slice_key)
    bins = [perm[((b - 1) * n + 1):(b * n)] for b in 1:B]
    [idx[sortperm(order_key[idx])] for idx in bins]
end

"""
    slice_anchors(bins, centres)

The vector that holds `centres[b]` at the indices of slice `b`: the anchor of each marker.
"""
function slice_anchors(bins::AbstractVector{<:AbstractVector{<:Integer}},
        centres::AbstractVector{T}) where {T <: Real}
    length(bins) == length(centres) ||
        throw(DimensionMismatch("$(length(bins)) slices but $(length(centres)) centres"))
    anchors = zeros(T, sum(length, bins))
    for (b, idx) in pairs(bins), a in idx

        anchors[a] = centres[b]
    end
    anchors
end
