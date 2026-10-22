# Research W P0 figure map

Wall cache used: `results/qualification_source_average` (source-cell averaged).

All plots read saved results only. PNG and SVG share a stem; underlying CSVs are in `figures/raw/`.

- `figures/01_1d_bandwidth.png` / `.svg`: Independent planar analytic regression; not a 2D witness. Non-nested frequency grids. Source: `results/p0a/results.json`.
- `figures/02_geometry.png` / `.svg`: Calibration material; exact-air cells include cavity and exterior; coefficient boundary cells are area averaged. Source: `results/qualification_source_average`.
- `figures/03_first_source_field.png` / `.svg`: Cached total field, source index 0. Calibration and task use the same source; per-seed shared magnitude scale. Source: `results/qualification_source_average`.
- `figures/04_convergence.png` / `.svg`: Full vector norms, not per-channel RMS. Adjacent mesh consistency is not a continuum certificate. Source: `results/qualification; results/qualification_source_average`.
- `figures/05_reference_pair_bandwidth.png` / `.svg`: Only two cached bandwidths; connecting lines guide the eye. Mean response distances, not sampled noisy data. No searched witness. Source: `results/qualification_source_average`.
- `figures/06_partial_cost.png` / `.svg`: No end-to-end inference cost or speedup measured. Wall setup timings are absent; separate benchmark is explicitly labeled. Source: `results/qualification_source_average`.
- `figures/07_fourblock_diagnostic.png` / `.svg`: Task-restricted 21-cell ports; discrete identity only on an unqualified continuum grid. Not a wall witness/ROM; no near-resonance claim. Source: `results/four_blocks/RESULTS.json`.

## Missing / not run

- Candidate search/witness: NOT_RUN; no candidate figure.
- Posterior/NN/3D/held-out probe validation/inference speedup: NOT_RUN; no fabricated figures.

## Render verification

All seven PNGs visually inspected. Geometry titles/colorbars, field panels, convergence legends, bandwidth annotations, partial-cost labeling and four-block footer are readable after layout-only repairs. SVGs were saved by the same figure objects; raw CSVs preserve plotted values. Billed jobs `figures-001`, `figures-002` and `figures-003` completed, totaling 20.7183 charged seconds. The latter two were layout-only rerenders. No physics calls, installation or source-cache changes occurred.

Interpretation limits: original wall nx64 is a baseline without an adjacent difference, so norms begin at nx96. The corrected wall data supply only the nx128-to192 adjacent norm; squares at nx192 show that measurement. The same-target total-data curve is not joint target ambiguity. Four-block lines connect two measured frequencies only; the diagnostic does not establish near-resonance behavior.
