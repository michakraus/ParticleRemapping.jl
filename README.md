# ParticleRemapping.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://michakraus.github.io/ParticleRemapping.jl/dev/)

Symplectic remapping, projection and redistribution of weighted particle-in-cell markers.

A particle-in-cell method samples a phase-space density with weighted markers. Over a long run the
markers cluster, and phase space develops regions with few or no markers. A redistribution step
moves the markers back towards an even cover. This package asks which redistributions of weighted
markers preserve the canonical symplectic structure of the marker system, and what they can
achieve. It collects four constructions behind one interface, `remap(method, markers)`, on a
common layer of marker sets, canonical structure, grids and spline bases.

| module | construction |
|:--|:--|
| `Projection` | Hamiltonian particle projection: kicks and drifts, or non-separable stages, of a weight-dependent Hamiltonian fitted to target points, with a charge-conserving reweighting and target weights from a dynamical low-rank density |
| `Redistribution` | spline-projection redistribution: one symplectic step of a pairwise spline Hamiltonian, with weights from a positive, moment-compatible spline projection |
| `Remapping` | gauge rezoning (positions and momenta fixed, velocities and weights from a constrained optimisation) and canonical transport with consistent weights (current-quantile matching), in 1D-1V |
| `SlicedRedistribution` | sliced quantile cotangent lifts followed by a canonical weight squeeze |

`using ParticleRemapping` makes the names of all four modules available. The plotting functions
are in the extension `ParticleRemappingMakieExt` and need a Makie backend, for example
`using CairoMakie`.

## Quick start

```julia
using ParticleRemapping

mk = Markers([0.3, 1.2, 2.0], [0.5, -1.0, 2.0], [0.2, 0.5, 0.3]; m = 2.0)
momentum(mk), total_momentum(mk), kinetic_energy(mk)

# the Hamiltonian particle projection with non-separable stages on a lattice with an elliptic hole
prob = hole_filling_problem(; nL = 16, hole_a = 0.9, hole_b = 0.8)
method = NonSeparableRemap(; nξ = 6, nv = 7, nℓ = 3, projection = HatProjection(12, 13), maxit = 5)
r = remap(method, prob.markers)
conservation_report(prob.markers, r.markers)

# the sliced redistribution of 900 markers after 500 steps in a two-wave electric field
using CairoMakie
mk = electric_wave_simulation(electric_wave_markers()).markers
r = remap(SlicedQuantileRemap(), mk)
plot_remap(mk, r.markers)
```

The home page of the manual has one example for each construction.

## Tests, scripts and documentation

The package root holds three environments besides the package itself. From the package root:

```bash
julia --startup-file=no --project=test test/runtests.jl
```

```bash
julia --startup-file=no --project=test test/runtests.jl core
```

```bash
julia --startup-file=no --project=docs docs/make.jl
```

The test group `core` runs in about 5 minutes. The group `slow` holds the two experiments of the
Hamiltonian particle projection and runs in about 7 minutes. Each construction has scripts in
`scripts/<construction>/` that print the numbers of its chapter next to the values of its source
note and check them, and a `make_figures.jl` that writes figures and animations to
`results/<construction>/`. Run them with `--project=scripts`.

## The manual

`docs/make.jl` builds one manual into `docs/build/`. It starts with the foundations (the canonical
structure of weighted markers, the factorisation lemma and the obstructions), has one chapter for
each construction, and closes with a comparison of the four, the open questions, a page on the
figures and animations, and a verification chapter. Every chapter of a construction has the same
pages: an overview with one example shown by the same pictures in every chapter, the theory, the
numerical experiments, and the properties, limitations, known issues, open questions and the cost
in higher dimensions. The verification chapter lists every test and script and the findings
against the source notes.

## Status

The package is not registered. Its repository is the private
[michakraus/ParticleRemapping.jl](https://github.com/michakraus/ParticleRemapping.jl). `KNOWN_ISSUES.md` lists the known defects and gaps.
`CHANGELOG.md` records the changes.
