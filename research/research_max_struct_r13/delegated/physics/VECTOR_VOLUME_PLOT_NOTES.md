# Authorized exact-grid vector volume plotting

Prepared `src/plot_vector_volumes.py`; it has not been run against the pending physics campaign, GT, or reconstruction arrays. Only syntax/static source inspection was performed. A preparation-time syntax error was corrected before delivery; no numerical or rendering attempt occurred. Existing runtime/solver code was not modified.

The parent may run it only after authorizing the complete terminal physics seal:

```text
python src/plot_vector_volumes.py --root R13_ROOT --seal AUTHORIZED_PHYSICS_SEAL.json --out NEW_PLOT_DIRECTORY
```

NumPy and the existing standard Matplotlib installation are sufficient. There is no download/install step, physical model import, forward solve, interpolation or refined isosurface construction. The reviewed `seal.manifest_snapshot` is the only case/method/truth-location manifest; no alternate shortlist or method filter is accepted. The authorization/terminal-queue guard runs before every NPZ/GT load and before numerical imports. Output must be outside sealed result/queue directories.

All four vector conditions and all seven methods appear in original manifest order. They remain three independent material objects: the crossgrid frame is a repeated geometry. Each method's saved initial/raw/primary/corrected material is shown. Missing/invalid arrays remain explicit UNAVAILABLE panels; out-of-box values remain in raw arrays and are annotated rather than omitted. There is no secondary substitution or new physical candidate selection.

The default is32 pages:4 conditions ×4 views ×2 method pages, with at most4 methods plus GT per page. Each page contains a real3D voxel row and three exact orthogonal central-slice rows. At most one figure is resident; it is closed after PNG/PDF export. Only the current case's small6³ arrays are retained. Camera(elevation25°,azimuth−55°), physical axes, equal aspect and color range[0,8] are fixed across every panel.

All cells are the exact saved coarse6³ material values. The displayed support is strictly chi>1, matching the predeclared IoU segmentation. Color is physical chi, without GT-dependent scaling or interpolation. Subthreshold values are preserved wherever the exact slice planes intersect them; complete216-cell arrays remain referenced in the panel table. Three central sections do not display every interior cell. On the even6³ grid the central positive plane is exact index3, at+.02lambda0, and is labeled as such; there is no averaging of the two middle planes. Surface voxels can occlude closed/off-center pores, so the exact YZ/XZ/XY sections are retained on the same page; a single camera does not certify the absence of an interior pore. These sections complement the scorer's existing slice montage and use the same saved arrays, with explicit plane labels.

Every title states synthetic full-vector CM+RR DDA, matched6³ or crossgrid8³→6³, and the missing general continuum-accuracy qualification. No Mie control is presented as a nonsphere accuracy certificate. Contrast saturation outside[0,8] is explicitly counted/labeled; values are not silently clipped in source arrays. Initial chi=.25 correctly has empty chi>1 voxel support, while exact slices preserve its nonzero material.

Outputs after execution are:

- `FIGURES/V3D__<full_case_ID>__<view>__pageNN.{png,pdf}`.
- `FIGURE_MANIFEST.{json,csv}`: stable figure IDs, case/object/method/view IDs, PNG/PDF relative download paths and absolute paths, display parameters and per-figure CPU/wall/export status.
- `FIGURE_PANEL_TABLE.{json,csv}`: source-array paths, saved selection metadata, availability, raw material ranges, support/illegal-cell counts, and the fixed plotting conventions.
- `PRIVATE_ACCESS_LOG.json`, `ISSUES.json`, `PLOT_COST_RECEIPT.json`: authorization provenance, all failures/missing panels and actual plotting costs. No material-recovery statistics or scientific acceptance are invented.

The cost receipt records actual current-run render attempts, successes/failures, import time, per-figure time and total input/plot/export CPU/wall, separately from online and physical-validation cost. Previously failed plotting attempts are unknown and explicitly reported as unknown, never inferred to be zero. New forward/LU/RHS counts remain zero because the script only displays saved arrays. Only the parent may make final scientific/visual acceptance and publication decisions after reviewing actual rendered output.

## Authorized execution update

The parent subsequently authorized offline access and the fixed renderer ran successfully:32 pages,112 method/view panels plus32 repeated GT references, no missing or out-of-box panels, no additional forward solve. Actual external plotting cost was18.863529CPU seconds/19.164855wall seconds. Preparation-time statements above remain historical; the completed figure manifest, access log, cost receipt and basic two-page visual QA are in `analysis/vector_volumes/`. Full scientific/visual acceptance remains with the parent.
