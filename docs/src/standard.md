# [Standard experiments](@id standard-experiments)

```@meta
CurrentModule = ParticleRemapping
```

Every construction of the manual remaps the same two marker sets and reports the same measures.
The page *Numerical experiments* of each chapter shows the results of its construction; the
[comparison](@ref comparison-standard) sets them side by side. The two marker sets come from the
sliced chapter, where they are suites 5a and 5b: a checkerboard with two empty quadrants, and
the markers of a long run in a two-wave electric field. They show the behaviour that the
[obstructions](@ref obstructions) predict: a map that is the same for every marker keeps the area
of a hole, a construction that keeps the positions cannot fill a hole in ``x``, and a construction
that couples the markers can.

## [The checkerboard](@id standard-checkerboard)

[`checkerboard_markers`](@ref): ``N = 400`` markers of weight ``1/400`` (``m = 1``), uniform in
the quadrants ``Q_1 = [0, 1)^2`` and ``Q_3 = [-1, 0)^2`` of ``[-1, 1]^2``; ``Q_2`` and ``Q_4`` hold
no marker, and both marginals are uniform. The box is periodic in ``x``. Two densities go with the
markers, as the keyword `variant` of [`checkerboard_density`](@ref):

- `:represented`, the density that the markers represent: ``1/2`` in ``Q_1`` and ``Q_3``, zero in
  ``Q_2`` and ``Q_4``. A remap that keeps the represented density cannot put weight into the empty
  quadrants; it can only move markers of small weight there.
- `:positive`, the uniform density ``1/4`` on ``[-1, 1]^2``. The markers are a sample of it with a
  hole: their weights ``1/N`` are those of a sample of ``1/4`` with the uniform marker density
  ``N/4`` on ``[-1, 1]^2``. A remap that keeps this density fills the empty quadrants with the weight
  ``1/2``.

The density matters only for the constructions that take one: the Hamiltonian particle
projection (its target weights) and the canonical transport (its target density). The other
constructions compute everything from the markers, so both variants give the same result.

[`standard_checkerboard`](@ref) runs the experiment for one construction and returns the quadrant
counts before and after, the weight in ``Q_2 \cup Q_4`` after, the largest displacement of a
position, the largest speed, the [`conservation_report`](@ref), and the counts of the markers in
``8 \times 8`` cells of ``[-1, 1]^2`` (the number of empty cells and the relative standard
deviation of the counts, [`phase_cell_statistics`](@ref)).

## [The long run in an electric field](@id standard-long-time)

[`standard_wave_markers`](@ref): ``30 \times 30`` markers at the cell centres of
``[-\pi, \pi) \times [-1.4, 1.4]`` with weights ``\propto f_0(x, v)``,
``f_0 = (1 + 0.35 \cos x)(0.40 + v^2) e^{-v^2/1.4}`` (``\sum w = 1``, ``m = 1.2``), pushed by 500
drift–kick–drift steps of ``\Delta t = 0.05`` through the field
``E(t, x) = -0.85 \sin x - 0.75 \sin(2x - 1.3 t)`` to ``T = 25``
([`electric_wave_markers`](@ref), [`electric_wave_simulation`](@ref)). The markers cluster along
the filaments of the two waves and leave 14 of the ``10 \times 10`` cells of
``[-\pi, \pi) \times [-2.25, 2.25]`` empty. Two densities go with them, as the keyword `variant`
of [`electric_wave_density`](@ref), both computed exactly by the inverse of the discrete steps:

- `:represented`: ``f_0`` at the foot of the characteristic, cut to ``|v| \le 1.4`` at ``t = 0``,
  where the lattice holds markers. This is the density that the markers represent.
- `:positive`: ``f_0`` at the foot of the characteristic without the cut. It is positive
  everywhere; the markers do not sample its tails.

[`standard_long_time`](@ref) remaps the markers once at ``T = 25`` and returns the cell
statistics before and after, the relative ``L^2`` difference of the kernel-smoothed densities
([`smoothed_density_difference`](@ref)), the largest speed and the
[`conservation_report`](@ref). [`standard_remap_history`](@ref) runs the simulation again with a
remap after every 100 steps (``t = 5, 10, 15, 20``) and records the number of empty cells and
the relative standard deviation of the cell counts every 10 steps, beside the run without
remaps.

```@example standard
using ParticleRemapping, CairoMakie
set_theme!(remap_theme())
plot_remap(checkerboard_markers(), standard_wave_markers();
    titles = ("(a) the checkerboard", "(b) the long run at T = 25"))
```

## [The constructions and their parameters](@id standard-parameters)

[`standard_method`](@ref) gives the method of each construction, keyed by
[`STANDARD_METHODS`](@ref). The constructions on ``[0, L)`` see the markers through a
[`ShiftedRemap`](@ref), which shifts the positions into ``[0, L)`` and back.

| key | method | checkerboard | long run |
|:--|:--|:--|:--|
| `:projection` | [`KickDriftRemap`](@ref ParticleRemapping.Projection.KickDriftRemap), 20 fit iterations | ``L = 2``, ``V = 1``, ``12 \times 13`` hats | ``L = 2\pi``, ``V = 3``, ``16 \times 17`` hats |
| `:nonseparable` | [`NonSeparableRemap`](@ref ParticleRemapping.Projection.NonSeparableRemap), ``6 \times 7 \times 3`` splines, 5 fit iterations | as above | as above |
| `:spline` | [`SplineRedistribution`](@ref ParticleRemapping.Redistribution.SplineRedistribution) | ``L = 2``, ``V_{\max} = 1.5``, ``\sigma_x = 0.5``, ``\sigma_v = 0.5`` | ``L = 2\pi``, ``V_{\max} = 3.5``, defaults |
| `:gauge` | [`GaugeRezoning`](@ref ParticleRemapping.Remapping.GaugeRezoning), ``\mathcal{K} = \{0, 2\}`` | 8 nodes | 16 nodes |
| `:transport` | [`CurrentQuantileMatching`](@ref ParticleRemapping.Remapping.CurrentQuantileMatching), ``\psi = \mathrm{id}``, a [`DensityTarget`](@ref ParticleRemapping.Remapping.DensityTarget) | 8 cells, ``|v| \le 1`` | 16 cells, ``|v| \le 3.5`` |
| `:sliced` | [`SlicedQuantileRemap`](@ref ParticleRemapping.SlicedRedistribution.SlicedQuantileRemap) | ``L_x = V_{\max} = 1``, ``B_v = 10``, ``B_x = 20``, ``h_x = h_v = 0.15`` | defaults |

The canonical transport takes the identity as its canonical map: the experiments ask what a
remap does to these markers, and no admissible map is known for them. The spline-projection
redistribution keeps its defaults at the full sizes, where its relaxation fails
([K2](@ref rd-experiments-failures)). Its chapter adds reduced runs with 16 markers, the size of
its verification set-up, with the relaxation step scaled to ``0.006 \cdot 16 / N``.

A construction that fails on a marker set reports the failure: `after` is `nothing` and `error`
holds the message. Such a failure is a result of the experiment.

## [The standard test set](@id standard-test-set)

The standard test set grades the two experiments above into ladders of cases, so that every
construction meets a weak and a strong form of each difficulty with the same measures. Each case
is a [`StandardCase`](@ref): a marker set on one box ``\Omega = [-L/2, L/2) \times [-V, V]``,
periodic in ``x``, with the exact density that goes with it. Where the markers represent a density
with holes, the case holds both variants of [`case_density`](@ref), `:represented` and
`:positive`, as the two experiments above do; otherwise it holds one density in both.

**The box.** Every marker of a case lies in its box: [`StandardCase`](@ref) refuses a marker with
``|v| > V``. Each case has one box for all its markers. For the families W, S1 and S4 the box
``V`` is a keyword, which a benchmark sweeps; the families C, E, S2 and S3 keep the box of their
density. A marker outside the box after a remap is a structure failure: the hat projection clamps
it onto the edge nodes of the box (`_hat_weights` in `src/projection/targets.jl`), so the next
deposit, remap or push sees it on the edge of the box.

**The families.** [`STANDARD_LADDERS`](@ref) lists their rungs, sizes and seeds.

- W, the wave ladder, [`wave_case`](@ref): the markers of the long run after the time ``T``, for
  ``T \in \{0, 1.25, 2.5, 5, 10, 15, 25\}`` and ``N \in \{16, 144, 576, 900\}``, on the box
  ``L = 2\pi`` and ``V = 3`` by default, for every rung, size and seed. Seed 0 pushes the lattice
  of the long run; the seeds 1 to 4 move each lattice point to a random point of its cell and
  weight it by ``f_0`` there. The densities are those of [`electric_wave_density`](@ref) at the
  same time; they differ wherever the characteristic starts at ``|v| > 1.4``, which the lattice
  does not fill. The rung ``T = 0`` is the first rung of the ladder, not an even filling of its
  box.
- C, the hole ladder, [`hole_case`](@ref): 400 markers of weight ``1/400`` on the occupied cells
  of the ``k \times k`` checkerboard of ``[-1, 1]^2``, ``k \in \{1, 2, 4, 8\}``, seeds 505 to
  509, on the box ``L = 2``, ``V = 1``. ``k = 2`` is the checkerboard above. The represented
  density is constant on the occupied cells and zero in the holes; the positive density is ``1/4``
  on the box.
- E, do no harm, [`lattice_case`](@ref): the plain and the deformed lattice of
  [`lattice_example`](@ref) with 144 and 576 markers, on the box ``L = 2\pi``, ``V = 4``, with
  their exact density. The markers already fill the box with no hole and represent their density,
  so a remap should not make them worse: E is the do-no-harm family.
- S, the stress cases, [`stress_case`](@ref), with 144 and 576 markers: S1, sparse tails, a
  sample of ``(1 + 0.1 \cos x)\, e^{-v^2/2}/\sqrt{2\pi}`` on the box ``V``, ``4`` by default; S2,
  weights over decades, a lattice weighted by a Maxwellian with a narrow beam at ``v = 4``, on
  ``V = 6``; S3, one dense cluster, a Gaussian sample of the widths ``0.1 L`` and ``0.6``, on
  ``V = 3``; S4, the clustered markers of [`clustered_example`](@ref), on the box ``V``, ``3.5``
  by default.

**The reference.** [`GridResample`](@ref) puts the ``N = n^2`` markers on the cell centres of an
``n \times n`` lattice of the box and weights them by a density: the lumped projection of the
markers on the hats of the box, or the exact density of the case
([`REFERENCE_METHODS`](@ref)). It keeps the total weight and nothing else; it is not symplectic.
It shows what a resample that ignores the structure gives on each case.

**The measures** compare the markers with the box and the density of the case:
[`case_evenness`](@ref), the cell counts on the box and on the support of the represented density
(and on the old window ``|v| \le 2.25`` for W); [`case_outside`](@ref), the markers outside the box
and the non-finite ones; [`case_representation`](@ref), the smoothed markers against the smoothed
density; [`case_charge_noise`](@ref), the grid charge against the exact one;
[`case_weights`](@ref), the spread of the weights; and, for a remap,
[`case_conservation`](@ref) and [`case_displacement`](@ref). [`case_reference`](@ref) computes
the integrals of the density once for both measures of a case.

**The driver.** [`case_method`](@ref) gives the method of a construction for a case: the
parameters of [`standard_method`](@ref) with the box of the case, and for the sliced map the slice
counts nearest to those of the experiment scaled with ``\sqrt{N}`` among the divisors of ``N``.
[`run_case`](@ref) runs a method on a case and returns one row of every measure before and after,
the conservation, the displacement and the state of the solver. A `DomainError`, an
`ArgumentError` or a failure of a factorisation in the remap is recorded with its class and
message, any other exception is thrown again, and a non-finite result is a failure; the row
records the markers outside the box after a remap and classifies nothing more.
