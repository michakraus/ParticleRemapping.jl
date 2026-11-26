# ParticleRemapping.jl

```@meta
CurrentModule = ParticleRemapping
```

*Symplectic remapping, projection and redistribution of weighted particle-in-cell markers.*

A particle-in-cell method samples a phase-space density with weighted markers. Over a long run
the markers cluster, and phase space develops regions with few or no markers. A redistribution
step moves the markers back towards an even cover. This package and its manual answer one
question: **which redistributions of weighted markers preserve the canonical symplectic structure
of the marker system, and what can they achieve?**

The answer has one common core. The canonical momentum of a marker is ``p_a = m w_a v_a``, so a
change of weight at a fixed momentum moves nothing in the canonical coordinates. Every symplectic
redistribution is a canonical map of positions and momenta followed by such a change of weights
(the [factorisation lemma](@ref factorisation-lemma)). A canonical map that is the same for
every marker cannot fill a hole ([Liouville](@ref liouville-holes)); a redistribution that does
fill holes changes the weights, or moves the markers with a map that depends on their weights or
couples them. The manual builds four constructions on this core, each with its own choice of the
canonical map and of the weights.

## The manual

**Foundations.**
[Background](background.md) states the problem and the requirements on a redistribution.
[Canonical structure of weighted markers](symplectic.md) derives the canonical coordinates, the
presymplectic form on positions, velocities and weights, the factorisation lemma and the gauge
step. [Symplectic redistribution](@ref obstructions) derives the obstructions that
shape every construction: Liouville's theorem for holes, weight-dependent Hamiltonians, the
velocity–momentum fold, and parameters that depend on the state.

**Constructions.**

1. [Hamiltonian particle projection](@ref projection-index) ([`Projection`](@ref)): kicks and
   drifts, or non-separable stages, of a weight-dependent Hamiltonian, fitted to target points
   of a density-interpolation flow; with a charge-conserving reweighting and target weights from
   a dynamical low-rank density.
2. [Spline-projection redistribution](@ref redistribution-overview)
   ([`Redistribution`](@ref)): one symplectic step of a translation- and permutation-invariant
   pairwise spline Hamiltonian, whose parameters solve a linear system, with weights from a
   positive, moment-compatible spline projection.
3. [Gauge rezoning and canonical transport](@ref remap-overview) ([`Remapping`](@ref)): the two
   extreme cases of the factorisation, in 1D-1V. The gauge rezoning keeps positions and momenta
   and sets the velocities by a constrained optimisation that keeps the grid charge and the
   kinetic energy. The canonical transport moves the markers by a given canonical map and sets
   weights that reproduce the density exactly.
4. [Sliced redistribution](@ref sliced-overview) ([`SlicedRedistribution`](@ref)): explicit
   sliced quantile cotangent lifts and a canonical weight squeeze.

Every chapter of a construction has the same shape:

- an **overview**: the algorithm step by step, with one SVG figure for each step and an animation
  of the steps on a small example, and one example shown by the same pictures in every chapter
  (the markers before and after the remap, the states of the remap stage by stage, an animation
  of the states, and the changes of the conserved quantities);
- the pages of the construction;
- the **algorithm** in full;
- the **numerical experiments**: the [standard experiments](@ref standard-experiments), which
  every construction runs with the same markers and the same measures;
- the additional experiments of the construction, with the hole-filling experiments on a page of
  their own;
- a page on its **properties, limitations and higher dimensions**, with the known issues and the
  open questions.

**Closing chapters.** [Comparison and open questions](comparison.md) sets the four side by side:
the results of the standard experiments, what each keeps and needs, and what binds first in
3d3v. [Verification](@ref verification) lists every test and script, what each checks, and which
numbers of the source notes the package reproduces. The [library](library.md) documents every
function, one page per module. [Figures and animations](@ref figures) states the common style
and lists every figure and animation of the manual and of the notes.

## Quick start

A marker set is a [`Markers`](@ref) object: positions, velocities, weights and the mass. The
common layer gives the canonical momenta and the conserved quantities:

```@example home
using ParticleRemapping
mk = Markers([0.3, 1.2, 2.0], [0.5, -1.0, 2.0], [0.2, 0.5, 0.3]; m = 2.0)
(p = momentum(mk), P = total_momentum(mk), kinetic = kinetic_energy(mk))
```

Every construction is a method type that [`remap`](@ref) applies to a marker set. The result is
a [`RemapResult`](@ref) with the new markers and the diagnostics of the method, and
[`conservation_report`](@ref) gives the changes of the conserved quantities (with a grid, also of
the grid charge, current and energy). Each call below uses the example of its chapter.

The Hamiltonian particle projection with non-separable stages, on a lattice with an elliptic
hole; `objective_ratio` is the fit objective after the fit relative to the identity, and `code`
the return code of the fit, `MAXITERS` after its budget of 5 iterations:

```@example home
prob = hole_filling_problem(; nL = 16, hole_a = 0.9, hole_b = 0.8)
method = NonSeparableRemap(; nξ = 6, nv = 7, nℓ = 3, projection = HatProjection(12, 13),
    maxit = 5)
r = remap(method, prob.markers)
(objective_ratio = r.info.objective_ratio, code = r.info.status.code,
    report = conservation_report(prob.markers, r.markers))
```

The spline-projection redistribution of 16 clustered markers; `target_error` is the scaled
distance of the end point from the relaxed targets:

```@example home
ex = redistribution_example()
r = remap(SplineRedistribution(), ex.markers)
(target_error = r.info.target_error, report = conservation_report(ex.markers, r.markers))
```

The gauge rezoning of 6400 clustered markers, with constraints on the grid charge and the kinetic
energy:

```@example home
s = gauge_rezoning_setup()
r = remap(s.method, s.markers)
conservation_report(s.markers, r.markers; grid = s.grid)
```

The current-quantile matching after the identity map, for 20000 markers sampled from a perturbed
Maxwellian; the weights match the target density up to the quantile error:

```@example home
q = quantile_matching_setup(; N = 20_000)
r = remap(q.identity, q.markers)
conservation_report(q.markers, r.markers; grid = q.grid)
```

The sliced redistribution of 900 markers after 500 time steps in a two-wave electric field:

```@example home
using CairoMakie
mk = electric_wave_simulation(electric_wave_markers()).markers
r = remap(SlicedQuantileRemap(), mk)
plot_remap(mk, r.markers; period = (-π, π),
    titles = ("(a) after 500 steps in the field", "(b) after the sliced remap"))
```

The plotting functions, such as [`plot_remap`](@ref), need a Makie backend: `using CairoMakie`
loads the extension `ParticleRemappingMakieExt`.

## Running the tests, the scripts and the documentation

The package root holds three environments besides the package itself. From the package root:

```
julia --startup-file=no --project=test test/runtests.jl            # all test files
julia --startup-file=no --project=test test/runtests.jl core       # without the three slow files
julia --startup-file=no --project=test test/remapping/gauge.jl     # one test file
julia --startup-file=no --project=scripts scripts/sliced/verify_sliced_redistribution.jl
julia --startup-file=no --project=scripts scripts/projection/make_figures.jl
julia --startup-file=no --project=docs docs/make.jl                # this manual
```

The group `core` runs in about 5 minutes; the group `slow` holds the two experiments of the
Hamiltonian particle projection and the standard experiments, and runs in about 21 minutes. Each
construction has scripts in `scripts/<construction>/` that print every number of its chapter next
to the value of its source and check them, and a `make_figures.jl` that writes the figures and
animations to `results/<construction>/`. The [verification](@ref verification) chapter lists them
all.
