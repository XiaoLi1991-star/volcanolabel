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
