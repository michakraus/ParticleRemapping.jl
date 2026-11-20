# [Figures and animations](@id figures)

```@meta
CurrentModule = ParticleRemapping
```

Every chapter of a construction shows its remap with the same three pictures, on the example of
its overview page, and uses one style for every figure and animation. This page states the style
and lists every figure and animation of the manual and of the notes.

## The standard pictures of a remap

| id | picture | function |
|:--|:--|:--|
| S1 | the markers before and after the remap, side by side, with the target points where the construction has them | [`plot_remap`](@ref) |
| S2 | one panel for each stage of the remap: the markers before it, after each stage of the canonical map ``\psi`` at the old weights, and after the gauge step | [`plot_stages`](@ref) of [`remap_states`](@ref) |
| A1 | the animation of S2: the markers move through the stages | [`animate_stages`](@ref) of [`remap_frames`](@ref) |
| T1 | the changes of the conserved quantities, as a table | [`conservation_report`](@ref) |

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
  panels or frames of a figure. When the weights span more than two decades, the colour scale is
  logarithmic. A figure whose point is another property of the markers, such as the slice of each
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
