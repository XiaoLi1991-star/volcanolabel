# volcanolabel 0.1.1

- Print resolved automatic plotting parameters from `volcano_plot()` with
  `message()` to make later agent-assisted layout and color tuning easier to
  reproduce.
- Refined README/release example figure with anonymized gene identifiers,
  wider outside label columns, and custom palette auto rings.
- Added optional `label_anchor_x_left` and `label_anchor_x_right` controls to
  `volcano_plot()` for publication figures that need fixed outside label
  columns.
- Increased the default gap between outside label anchor dots and label text so
  text does not visually sit on top of the anchor dot.
- Redraw labelled source points above their rings so the ring highlights the
  point without covering or clipping the original filled circle.

# volcanolabel 0.1.0

Initial release.

- Added `volcano_plot()` for plotting-ready volcano plot data.
- Added adaptive outside label layout with side-aware label anchoring, label
  wrapping, x-range expansion, and compact elbow connectors.
- Added group-aware labelled-point rings. `label_point_ring_color = "auto"`
  derives visible ring colors from each point group, including up, down, and
  normal groups; `"group"` keeps exact group colors.
- Added safeguards against white label backgrounds covering points.
- Added `save_volcano()`, `theme_volcano()`, `volcano_palette()`,
  `wrap_volcano_labels()`, and `compute_outside_label_layout()`.
