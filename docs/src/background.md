# Background

This chapter sets out the problem that the whole manual addresses. The weighted markers of a
particle-in-cell method cluster in phase space, and a redistribution that moves them back should
respect the Hamiltonian structure of the marker system. The chapter states the requirements on
such a redistribution and gives the road map of the constructions. The canonical structure itself
is the subject of the next chapter, [Canonical structure of weighted markers](symplectic.md).

## The Vlasov equation and its particle discretisation

Consider one species of charged particles with mass ``m`` and charge ``q`` in a self-consistent
or external electric field ``E``. The phase-space density ``f(x, v, t)`` solves the Vlasov
equation

```math
\frac{\partial f}{\partial t} + v \cdot \nabla_x f + \frac{q}{m} E \cdot \nabla_v f = 0 .
```

A particle method samples ``f`` with ``N`` markers. Marker ``a`` has a position ``x_a``, a
velocity ``v_a`` and a weight ``w_a > 0``, and the density is represented by the Klimontovich
sum

```math
f(x, v) \approx f_N(x, v) = \sum_{a=1}^{N} w_a \, \delta(x - x_a) \, \delta(v - v_a) .
```

The markers move along the characteristics ``\dot{x}_a = v_a``,
``\dot{v}_a = (q/m) E(x_a)``. The density is constant along a characteristic, so the weights do
not change in time. The total weight ``\sum_a w_a`` is the represented mass (in units of the
particle number), and ``m \sum_a w_a v_a`` the represented momentum.

### Marker density and weights

Let ``g(x, v)`` be the density of the markers in phase space, normalised to ``\int g = N``.
Comparing ``f_N`` with ``f = w \, g`` gives the consistency condition

```math
w_a = \frac{f(x_a, v_a)}{g(x_a, v_a)} .
```

Two loadings are common:

- **uniform loading:** the markers are equidistributed in a box ``\Omega`` of phase space,
  ``g = \bar{g} = N / |\Omega|``, and the weights are ``w_a = c \, f(x_a, v_a)`` with
  ``c = |\Omega| / N``;
- **importance sampling:** the markers are drawn from ``g \propto f`` and all weights are equal.

### The canonical momentum

The weight enters the mechanics of the markers. In the Lagrangian of [Low:1958](@citet),
discretised with weighted markers, marker ``a`` carries the kinetic energy
``\frac{m}{2} w_a |\dot{x}_a|^2``, so its canonical momentum is

```math
p_a = m \, w_a \, v_a .
```

The positions and these momenta are canonical coordinates of the marker system, and the total
momentum ``P = \sum_a p_a`` is the physical momentum of the represented density. The
[next chapter](@ref canonical-coordinates) derives this and the symplectic form that it carries.

## Clustering and holes

The characteristic flow preserves phase-space volume. In the continuum limit a uniform marker
density therefore stays uniform. For a finite number of markers this is not true at the
resolution scale. Phase mixing stretches the initial arrangement into filaments, the local marker
spacing becomes anisotropic, and regions appear that contain few or no markers. Such regions are
called **holes**. A hole increases the noise of every moment that the markers deposit on a grid,
and it leaves a part of phase space unresolved. With importance sampling, a hole in ``f`` is also
a hole in the marker distribution.

The following example loads ``40 \times 40`` markers on a lattice in ``[0, 2\pi) \times [-4, 4]``
with weights ``w \propto f``, for ``f(x, v) = (1 + 0.3 \cos x) \, e^{-v^2/2} / \sqrt{2\pi}``, and
pushes them in the fixed field ``E(x) = -0.5 \sin x`` (with ``q/m = 1``) with the leapfrog scheme.
The panels count the cells of a ``16 \times 16`` grid that hold no marker.

```@example background
using ParticleRemapping, CairoMakie

f(x, v) = (1 + 0.3cos(x)) * exp(-v^2 / 2) / sqrt(2π)
E(x) = -0.5sin(x)

n = 40
x0 = vec([(i - 0.5) * 2π / n for i in 1:n, j in 1:n])
v0 = vec([-4 + (j - 0.5) * 8 / n for i in 1:n, j in 1:n])
w0 = f.(x0, v0) .* (2π * 8 / n^2)

function push(x, v, T; dt = 0.01)
    x, v = copy(x), copy(v)
    for _ in 1:round(Int, T / dt)
        v .+= dt / 2 .* E.(x)
        x .= mod.(x .+ dt .* v, 2π)
        v .+= dt / 2 .* E.(x)
    end
    x, v
end

empty_cells(x, v) = 16^2 - length(unique(zip(floor.(Int, x ./ (2π / 16)),
    floor.(Int, (v .+ 4) ./ (8 / 16)))))

fig = Figure(size = (1000, 330))
for (k, T) in enumerate((0.0, 10.0, 40.0))
    x, v = push(x0, v0, T)
    ax = Axis(fig[1, k]; xlabel = "x", ylabel = "v", limits = (0, 2π, -4, 4),
        title = "t = $(T), empty cells: $(empty_cells(x, v))")
    plot_markers!(ax, Markers(x, v, w0); colorrange = extrema(w0))
end
fig
```

## Remapping

A **redistribution** (remapping, resampling) step moves the markers back towards a regular
arrangement and assigns new weights so that the represented density does not change. In its usual
form it deposits the markers on a grid and samples new ones from the grid density. This is a
coarse-graining step: it adds numerical diffusion and it does not respect the Hamiltonian structure
of the marker system. This manual asks which redistributions of weighted markers are compatible
with the canonical structure, and constructs four of them.

### [Requirements](@id requirements)

A *redistribution map* takes the markers ``(x_a, v_a, w_a)_a`` to new markers
``(\bar{x}_a, \bar{v}_a, \bar{w}_a)_a``, with the new momenta
``\bar{p}_a = m \bar{w}_a \bar{v}_a``. The requirements are:

- **(R1) Symplecticity.** The map ``(x, p) \mapsto (\bar{x}, \bar{p})`` preserves the canonical
  form ``\omega = \sum_a \mathrm{d}x_a \wedge \mathrm{d}p_a``.
- **(R2) Momentum.** ``\bar{P} = P``.
- **(R3) Consistency.** The new weights represent the same density:
  ``\bar{w}_a = c' f_h(\bar{x}_a, \bar{v}_a)``, where ``f_h`` is an approximation of the
  represented density and ``c'`` a normalisation.
- **(R4) Equidistribution.** The markers ``(\bar{x}_a, \bar{v}_a)`` are approximately
  equidistributed afterwards.

Each construction meets (R1) exactly. Which of (R2)–(R4) it meets exactly, which approximately,
and which further quantities (grid charge, kinetic energy, total weight) it keeps, is stated in its
chapter and collected in the [comparison](comparison.md).

## [The constructions](@id constructions)

The weights add freedom. A change of weight at a fixed canonical momentum,
``(x, v, w) \mapsto (x, \lambda v, w / \lambda)``, is the identity on ``(x, p)``, so it is
symplectic for any ``\lambda > 0``, even one that depends on the state of the whole ensemble. The
[factorisation lemma](@ref factorisation-lemma) of the next chapter shows that *every* symplectic
redistribution is a canonical map ``\psi`` of ``(x, p)`` followed by such a
[gauge step](@ref gauge-step). The chapter after it,
[Symplectic redistribution](@ref obstructions), shows what this leaves possible: a
smooth canonical map applied to every marker moves a hole but cannot shrink it
([Liouville](@ref liouville-holes)), so a redistribution that fills holes either changes the
weights, or moves the markers with a map that depends on their weights or couples them.

The four constructions of the manual differ in how they choose ``\psi`` and the new weights:

| chapter | module | the canonical map ``\psi`` | the new weights | kept exactly |
|:--|:--|:--|:--|:--|
| [Hamiltonian particle projection](@ref projection-index) | [`Projection`](@ref) | kicks and drifts, or non-separable stages, of a weight-dependent Hamiltonian, fitted to target points | a projected density, or a charge-conserving reweighting | ``\omega`` for fixed parameters, ``P``, total weight; grid charge in the charge-conserving variant |
| [Spline-projection redistribution](@ref redistribution-overview) | [`Redistribution`](@ref) | one symplectic step of a pairwise spline Hamiltonian, solved to reach relaxed targets | a positive spline projection with exact mass and momentum | ``\omega`` for fixed parameters, ``P``, total weight |
| [Gauge rezoning and canonical transport](@ref remap-overview) | [`Remapping`](@ref) | the identity (gauge rezoning), or a given canonical map (canonical transport) | stretches from a constrained optimisation; or current-quantile matching | ``\omega`` (the gauge rezoning for any state dependence), ``x``, ``p``, grid current; grid charge and kinetic energy by constraints |
| [Sliced redistribution](@ref sliced-overview) | [`SlicedRedistribution`](@ref) | sliced quantile cotangent lifts and a canonical weight squeeze, explicit | the squeeze multipliers | ``\omega`` for fixed slices and multipliers, ``P``, total weight |

The first two chapters move the markers to target points with a map whose parameters are fitted
or solved for. The third keeps the canonical map fixed, the identity or a given map, and chooses
the new weights. The fourth builds an explicit map that couples the markers. Every chapter runs
the same [standard experiments](@ref standard-experiments). Where the
parameters, slices or multipliers are computed from the markers, the map is symplectic for their
frozen values; [Parameters that depend on the state](@ref state-dependent-parameters) states what
this means and when it can be closed. Each
construction implements the common interface [`remap`](@ref)`(method, markers)`, which returns a
[`RemapResult`](@ref) with the new markers and the diagnostics of the method.
[`conservation_report`](@ref) compares two marker sets.

The constructions come from four sets of research notes: *Hamiltonian remapping of weighted
particles in 1D-1V*, *Hamiltonian particle projection: a symplectic redistribution of weighted
markers*, *Symplectic redistribution of weighted particles in phase space*, and a derivation of
the sliced redistribution, whose note is written from the [sliced chapter](@ref sliced-overview).
The
[verification](@ref verification) chapter lists which number of each source the package
reproduces.

## [Notation](@id notation)

| symbol | meaning |
|:--|:--|
| ``N``, ``a`` | number of markers, marker index |
| ``d`` | number of position (and velocity) components |
| ``m``, ``q`` | mass and charge of the species |
| ``x_a``, ``v_a``, ``w_a`` | position, velocity and weight of marker ``a`` |
| ``p_a = m w_a v_a`` | canonical momentum |
| ``P = \sum_a p_a`` | total momentum |
| ``\omega = \sum_a \mathrm{d}x_a \wedge \mathrm{d}p_a`` | canonical symplectic form on ``(x, p)`` |
| ``\Pi(x, v, w) = (x, p)``, ``\omega_N = \Pi^* \omega`` | kinematic projection, presymplectic form on ``(x, v, w)`` |
| ``\mathcal{Z}``, ``\mathcal{M}``, ``\zeta = (\xi, \eta)`` | ensemble phase space ``(x, p)``, single-marker phase space and its points |
| ``\mathbb{J}`` | canonical matrix ``\begin{pmatrix} 0 & I \\ -I & 0 \end{pmatrix}`` |
| ``f``, ``f_h``, ``g``, ``n`` | represented density, its approximation, marker density in ``(x, v)``, marker density in ``(x, p)`` |
| ``L``, ``N_g``, ``\Delta x`` | length of the periodic domain, number of grid nodes, grid spacing |
| ``S_i`` | linear B-spline (hat) shape function of node ``i`` |
| ``\mathcal{M}_i^{(k)}`` | grid moment ``\sum_a w_a S_i(x_a) v_a^k`` |
| ``\psi``, ``\lambda_a`` | canonical map and stretches of a redistribution |

Each chapter states how it marks the image of its map (a bar, a hat, a prime or a capital
letter) and its own symbols.
