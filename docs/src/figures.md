# [Figures and animations](@id figures)

```@meta
CurrentModule = ParticleRemapping
```

Every chapter of a construction shows its algorithm, its remap and its standard experiments with
the same pictures, and uses one style for every figure and animation. This page states the style
and lists every figure and animation of the manual and of the notes.

## The standard pictures of a remap

| id | picture | function |
|:--|:--|:--|
| S1 | the markers before and after the remap, side by side, with the target points where the construction has them | [`plot_remap`](@ref) |
| S2 | one panel for each stage of the remap: the markers before it, after each stage of the canonical map ``\psi`` at the old weights, and after the gauge step | [`plot_stages`](@ref) of [`remap_states`](@ref) |
| A1 | the animation of S2: the markers move through the stages | [`animate_stages`](@ref) of [`remap_frames`](@ref) |
| T1 | the changes of the conserved quantities, as a table | [`conservation_report`](@ref) |
| G1 | the steps of the algorithm on a small example: one SVG panel per step, with the overlays of the step (grids, slices, rays, targets, arrows) | `algorithm_svg(key).panels` ([`algorithm_svg`](@ref)) |
| G2 | the interactive Vue/SVG player (`AlgorithmPlayer`) through the steps of G1, fed from JSON | `algorithm_svg(key).player`, [`write_algorithm_player`](@ref) |
| X1 | the markers before and after the remap of each [standard experiment](@ref standard-experiments), and the cell counts of the run with repeated remaps | [`plot_remap`](@ref), `plot_cell_count_history` |

Between two states, the animation changes the positions, the momenta ``p = m w v`` and the
weights linearly ([`interpolate_stages`](@ref)), unless the construction follows a path of its own,
such as a homotopy or the flow of a Hamiltonian. A canonical stage at fixed weights then moves the
markers on straight lines in ``(x, v)``, and a gauge step at fixed ``(x, p)`` moves each marker
along its line ``v \propto 1/w``.

## Style

[`remap_theme`](@ref) sets the fonts, Wong's colour cycle and the defaults; the documentation and
the figure scripts apply it with `set_theme!(remap_theme())`. The plots of every construction
follow these rules:

- **Markers in phase space** are coloured and sized by weight: the colour map `viridis`, the
  marker area proportional to the weight, and one colour range and one size reference for all
  panels or frames of a figure. For 200 markers or more, the colour range and the size reference
  run from the 0.5 % to the 99.5 % quantile of the weights, so that a few extreme weights do not
  wash out the others; the colour scale is linear unless a figure asks for `colorscale = log10`.
  A figure whose point is another property of the markers, such as the slice of each
  marker, colours by that property and says so in its caption.
- **Target points** are open circles in vermillion; **displacements** are thin grey segments;
  **outlines** of holes and regions are dark grey.
- **Curves and series** take the colours of Wong's palette in order.
- **Fields** use `cividis` when they are non-negative (a density, a count, a magnitude) and
  `balance` with a symmetric colour range when they are signed.
- **Categories**, such as slices, cycle through Wong's palette.
- **Panels** of a multi-panel figure are titled "(a) …", "(b) …". Axis labels of quantities are
  LaTeX strings, such as ``x`` and ``v``; labels in words are plain text.
- **Animations** are `.mp4` files at 20 frames per second, with the stage or the parameter in the
  title, and they hold the first and the last frame for a moment.

## Inventory

The kinds of pictures: the standard pictures S1, S2, A1, T1, G1, G2 and X1 above; C, a step of the
construction; W, a worked example in closed form; E, a numerical experiment; N, a solver and its
convergence; L, a long run with repeated remaps; D, a demonstration of a structural fact. A
figure of a note is labelled by its `\label`; an animation is named in the note by its file. The
figure scripts `scripts/<construction>/make_figures.jl` write every file to
`results/<construction>/`, from where the notes include them.

**Coverage.** Every construction has G1 and G2 on its overview page, S1, S2, A1 and T1 on the
example of its overview page, and X1 on its page *Numerical experiments*; each note shows S1 and
S2. The kinds C, E, N and D differ between the constructions, because each shows the steps and
facts of its own construction. The long run with repeated remaps (L) of the standard experiments
runs for the gauge rezoning, the canonical transport and the sliced redistribution; the
projection fails at its second or third remap, and the spline-projection redistribution does not
reach the marker counts of the run with its defaults (K2).

### Hamiltonian particle projection

| picture | kind | function | page | note | file | weight-coded |
|:--|:--|:--|:--|:--|:--|:--|
| before and after, with targets | S1 | [`plot_remap`](@ref) | [overview](@ref projection-example) | `fig:standard-remap` | `projection-remap.pdf` | yes |
| before, kicks, drifts, gauge step | S2 | [`plot_stages`](@ref) | [overview](@ref projection-example) | `fig:standard-stages` | `projection-stages.pdf` | yes |
| animation with the stage parameters scaled from zero | A1 | [`animate_stages`](@ref) | [overview](@ref projection-example) | named | `projection-stages.mp4` | yes |
| conserved quantities | T1 | [`conservation_report`](@ref) | [overview](@ref projection-example) | text | — | — |
| clustered markers before and after, with targets | S1 | [`plot_remap`](@ref) | [clustered markers](@ref projection-clustered-example) | — | — | yes |
| clustered markers, animated | A1 | [`animate_stages`](@ref) | [clustered markers](@ref projection-clustered-example) | — | `projection-clustered.mp4` (manual only) | yes |
| markers, ``f_h`` and ``g_h`` | C | `plot_lumped_projection` | [construction](@ref projection-construction) | — | `projection-projection.pdf` | yes in (a); (b), (c) are fields |
| markers, targets and displacements | C | `plot_target_flow` | [construction](@ref projection-construction) | — | `projection-target-flow.pdf` | yes |
| the density-interpolation flow | C | [`animate_stages`](@ref) | [construction](@ref projection-construction) | named | `projection-moser-targets.mp4` | yes |
| deformed lattice with the lumped projection | C | [`plot_remap`](@ref) | [algorithm](@ref projection-algorithm) | — | — | yes |
| kick–drift stages on the deformed lattice, charge-conserving | C | [`plot_stages`](@ref) | [algorithm](@ref projection-algorithm) | — | — | yes |
| non-separable stages on the deformed lattice, charge-conserving | C | [`plot_stages`](@ref) | [non-separable stages](@ref projection-nonseparable) | — | — | yes |
| hole filling, ``x_h = \pi/2`` and ``x_h = \pi`` | E | `plot_hole_filling` | [hole filling](@ref projection-hole-filling) | `fig:hole-pi2`, `fig:hole-pi` | `projection-hole-pi2.pdf`, `projection-hole-pi.pdf` | yes |
| the charge correction | E | `plot_charge_correction` | [hole filling](@ref projection-hole-filling) | `fig:charge` | `projection-charge.pdf` | no: the colour marks the run |
| target flow near the hole, ``x_h = \pi`` | E | `plot_target_flow` | [hole filling](@ref projection-hole-filling) | — | — | yes |
| kick–drift stages near the hole, ``x_h = \pi/2`` | E | `plot_hole_stages` | [hole filling](@ref projection-hole-filling) | — | `projection-stages-kick-drift-pi2.pdf` | yes |
| non-separable stages near the hole, ``x_h = \pi`` | E | `plot_hole_stages` | [hole filling](@ref projection-hole-filling) | — | `projection-stages-non-separable-pi.pdf` | yes |
| the hole-filling remaps, animated | E | [`animate_stages`](@ref) | [hole filling](@ref projection-hole-filling) | named | `projection-remap-kick-drift-pi2.mp4`, `projection-remap-non-separable-pi.mp4` | yes |
| target weights after a push | E | `plot_dlra_weights` | [hole filling](@ref projection-push-test) | `fig:dlra` | `projection-dlra.pdf` | yes in (a), (b); (c) is a signed field, (d), (e) are series |
| the low-rank density along the push | E | `animate_dlra` | [hole filling](@ref projection-push-test) | named | `projection-dlra.mp4` | no: the density is the point |
| a mean-field step and the image of a hole | D | `plot_mean_field_step` | [obstructions](@ref obstructions) | `fig:structure` | `projection-mean-field.pdf` | no: the markers carry no weights |
| the fold of the implicit weight | D | `plot_implicit_weight_fold` | [non-separable stages](@ref projection-nonseparable) | `fig:structure` | `projection-fold.pdf` | — (a curve) |
| the coupling of ``+v`` and ``-v`` | D | `plot_pm_v_coupling` | [non-separable stages](@ref projection-nonseparable) | `fig:structure` | `projection-pm-v.pdf` | — (curves) |

### Spline-projection redistribution

| picture | kind | function | page | note | file | weight-coded |
|:--|:--|:--|:--|:--|:--|:--|
| before and after, with the relaxed targets | S1 | [`plot_remap`](@ref) | [overview](@ref redistribution-overview) | Fig. 1 | `redistribution-remap.pdf` | yes |
| before, targets, step at the old weights, gauge step | S2 | [`plot_stages`](@ref) | [overview](@ref redistribution-overview) | Fig. 2 | `redistribution-stages.pdf` | yes |
| the step in pseudo-time, then the gauge step | A1 | [`animate_stages`](@ref) | [overview](@ref redistribution-overview) | named | `redistribution-stages.mp4` | yes |
| conserved quantities | T1 | [`conservation_report`](@ref) | [overview](@ref redistribution-overview) | properties table | — | — |
| the velocity–momentum fold | D | `plot_velocity_momentum_fold` | [obstructions](@ref velocity-momentum-fold) | Fig. 3 | `redistribution-fold.pdf` | — (curves) |
| ``L^2`` and KKT projections | C | `plot_spline_projections` | [projection](@ref rd-projection) | Fig. 4 | `redistribution-projections.pdf` | size; the colour shows the signed density |
| relaxation paths and clustering | C | `plot_target_relaxation` | [projection](@ref rd-projection) | Fig. 5 | `redistribution-relaxation.pdf` | no: the relaxation moves points without weights |
| the relaxation, animated | C | `animate_target_relaxation` | [projection](@ref rd-projection) | named | `redistribution-relaxation.mp4` | no: as above |
| Newton residuals | N | `plot_newton_convergence` | [solver](@ref rd-solver) | Fig. 6 | `redistribution-newton.pdf` | — (curves) |
| closed-loop defect matrices | D | `plot_closed_loop_defect` | [closed loops](@ref rd-closed-loop) | Fig. 7 | `redistribution-closed-loop.pdf` | no: signed matrices |
| 900 targets leaving the box (K2) | E | `plot_target_relaxation` | [verification set-up](@ref rd-setup) | Fig. 8 | `redistribution-relaxation-900.pdf` | no: points without weights |
| the lattice and the clustered set-up | E | [`plot_remap`](@ref) | [verification set-up](@ref rd-setup) | — | — | yes |
| the stages of a single symplectic Euler step | E | [`plot_stages`](@ref) | [verification set-up](@ref rd-setup) | — | — | yes |

### Gauge rezoning and canonical transport

| picture | kind | function | page | note | file | weight-coded |
|:--|:--|:--|:--|:--|:--|:--|
| gauge rezoning, before and after | S1 | [`plot_remap`](@ref) | [overview](@ref remap-overview) | `fig:gauge-remap` | `remapping-gauge-remap.pdf` | yes; no targets, since 6400 open circles would hide the markers |
| gauge rezoning, states along the homotopy | S2 | [`plot_stages`](@ref) | [overview](@ref remap-overview) | `fig:gauge-stages` | `remapping-gauge-stages.pdf` | yes |
| gauge rezoning, animation along the homotopy | A1 | [`animate_stages`](@ref) | [overview](@ref remap-overview) | named | `remapping-gauge-stages.mp4` | yes |
| gauge rezoning, conserved quantities | T1 | [`conservation_report`](@ref) | [overview](@ref remap-overview) | properties table | — | — |
| canonical transport, before and after | S1 | [`plot_remap`](@ref) | [overview](@ref remap-overview) | `fig:transport-remap` | `remapping-transport-remap.pdf` | yes |
| canonical transport, before, ``\psi`` at the old weights, matched | S2 | [`plot_stages`](@ref) | [overview](@ref remap-overview) | `fig:transport-stages` | `remapping-transport-stages.pdf` | yes |
| canonical transport, animation | A1 | [`animate_stages`](@ref) | [overview](@ref remap-overview) | named | `remapping-transport-stages.mp4` | yes |
| canonical transport, conserved quantities | T1 | [`conservation_report`](@ref) | [overview](@ref remap-overview) | properties table | — | — |
| one charge constraint and its fold | W | `plot_charge_example` | [gauge](@ref remap-charge-example) | `fig:ex1` | `remapping-ex1.pdf` | size; the colour identifies the marker across the panels |
| the fold, animated | W | `animate_charge_fold` | [gauge](@ref remap-charge-example) | named | `remapping-ex1-fold.mp4` | size; the colour identifies the marker |
| eight markers on two nodes | W | `plot_eight_marker_example` | [algorithm](@ref remap-eight-example) | `fig:ex2` | `remapping-ex2.pdf` | yes |
| eight markers, homotopy | W | `animate_eight_marker_homotopy` | [algorithm](@ref remap-eight-example) | named | `remapping-ex2-homotopy.mp4` | yes |
| pattern of the KKT matrix and its Schur complement | C | `plot_kkt_structure` | [gauge](@ref remap-gauge) | — | `remapping-kkt.pdf` | no: a matrix (field colour map) |
| gauge rezoning of 6400 markers, stretches and gaps | E | `plot_gauge_experiment` | [examples](@ref remap-gauge-example) | `fig:gauge` | `remapping-gauge.pdf` | yes in (a), (b); (c), (d) are distributions |
| the homotopy in cells 1 and 2 | E | [`animate_stages`](@ref) | [examples](@ref remap-gauge-example) | named | `remapping-gauge-homotopy.mp4` | yes |
| the cotangent lift and the matching, densities and residuals | W, E | `plot_lift_example` | [examples](@ref remap-transport-example) | `fig:ex3` | `remapping-ex3.pdf` | no: curves and counts per cell |
| the lift with the matching, animated | E | `animate_lift_matching` | [examples](@ref remap-transport-example) | named | `remapping-ex3-lift.mp4` | yes |
| five markers in one cell | W | `plot_matching_example` | [algorithm](@ref remap-five-example) | `fig:ex4` | `remapping-ex4.pdf` | size; the colour links quantile, velocity and weight of a marker |
| five markers, animated | W | `animate_matching_example` | [algorithm](@ref remap-five-example) | named | `remapping-ex4-matching.mp4` | size; the colour identifies the marker |
| gauge orbits | D | `plot_gauge_orbits` | [canonical structure](@ref gauge-step) | — | `remapping-gauge-orbits.pdf` | yes |
| the velocity–momentum sheet | D | `plot_velocity_momentum_sheet` | library | — | `remapping-sheet.pdf` | no: curves |

### Sliced redistribution

| picture | kind | function | page | note | file | weight-coded |
|:--|:--|:--|:--|:--|:--|:--|
| before and after | S1 | [`plot_remap`](@ref) | [overview](@ref sliced-overview) | `fig:standard-remap` | `sliced-remap.pdf` | yes |
| before, ``\Psi_q``, ``\Psi_p``, ``T_\mu``, gauge step | S2 | [`plot_stages`](@ref) | [overview](@ref sliced-overview) | `fig:standard-stages` | `sliced-stages.pdf` | yes |
| animation along the canonical path | A1 | [`animate_stages`](@ref) | [overview](@ref sliced-overview) | named | `sliced-stages.mp4` | yes |
| conserved quantities | T1 | [`conservation_report`](@ref) | [overview](@ref sliced-overview) | text | — | — |
| the stages coloured by slice | C | `plot_sliced_stages` | [algorithm](@ref sliced-algorithm) | `fig:stages` | `sliced-slices.pdf` | partly: the lift panels colour by slice, which is their point |
| the stages coloured by slice, animated | C | `animate_sliced_stages` | [algorithm](@ref sliced-algorithm) | named | `sliced-slices.mp4` | no: the slices are the point |
| the checkerboard | E | `plot_checkerboard_filling` | [hole filling](@ref sliced-holes) | `fig:checkerboard` | `checkerboard.pdf` | no: equal weights; the colour shows the quadrant of origin |
| a remap in the long run, markers and cell counts | E | `plot_wave_remap` | [long-time run](@ref sliced-long-time-page) | `fig:wave-remap` | `wave-remap.pdf` | yes; the counts are a field |
| empty cells and spread of the counts over time | L | `plot_cell_count_history` | [long-time run](@ref sliced-long-time-page) | `fig:history` | `cell-count-history.pdf` | — (curves) |
| the long runs with and without remaps | L | `animate_long_time_run` | [long-time run](@ref sliced-long-time-page) | named | `long-time-run.mp4`, `long-time-run-no-remap.mp4` | yes |
| per-axis lifts on a hole between the axes | D | table of an `@example` | [hole filling](@ref sliced-hole-between-axes) | numbers | — | — |
