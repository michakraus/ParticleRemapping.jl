# Changelog

## [Unreleased]

### Added

- The initial package. It consolidates the notes and scripts on the symplectic remapping,
  projection and redistribution of weighted particle-in-cell markers behind one interface:
  `remap(method, markers)` returns a `RemapResult` with the new markers and the diagnostics of the
  method, and `conservation_report` gives the changes of the conserved quantities.
- A common layer: the `Markers` type with canonical momenta and conserved quantities, the
  canonical structure (`canonical_J`, `symplecticity_defect`, `presymplectic_matrix`,
  `pullback_residual`, `gauge_step`), periodic grids, cubic B-splines and spline bases.
- Four modules, one for each construction: `Remapping` (gauge rezoning and canonical transport
  with consistent weights), `Projection` (Hamiltonian particle projection with kick–drift and
  non-separable stages, charge-conserving reweighting and dynamical low-rank target weights),
  `Redistribution` (spline-projection redistribution) and `SlicedRedistribution` (sliced quantile
  cotangent lifts and the canonical weight squeeze). Each module also holds the worked examples and
  experiments of its note.
- The Makie extension `ParticleRemappingMakieExt`, with the figures and animations of every note.
- One manual, built by `docs/make.jl`: foundations, a chapter for each construction in one
  notation, a comparison with the open questions, a verification chapter with the findings against
  the source notes, and a library page for each module.
- A test suite of 41 test files in the groups `core` and `slow`, and scripts in
  `scripts/<construction>/` that reproduce the numbers and figures of each note.
- `KNOWN_ISSUES.md` with the entries K1–K8.
