# Fixed-case material and field display handoff

Code: `tools/g0_material_fields.py` — **CODE_READY_NOT_EXECUTED**. New script execution/render inspection by this worker: **NOT_RUN**.

Default generated paths, after the parent serial wrapper runs the script:

- `figures/g0_material_fields/C4_1_ACTUAL_MATERIALS.png`
- `figures/g0_material_fields/C4_1_ACTUAL_MATERIALS.pdf`
- `figures/g0_material_fields/C4_1_ACTUAL_MATERIALS_DATA.npz`
- `figures/g0_material_fields/STRONG_CENTER_SOURCE0_ACTUAL_TOTAL_FIELD.png`
- `figures/g0_material_fields/STRONG_CENTER_SOURCE0_ACTUAL_TOTAL_FIELD.pdf`
- `figures/g0_material_fields/STRONG_CENTER_SOURCE0_ACTUAL_TOTAL_FIELD_DATA.npz`
- `research/delegated/results_support/MATERIAL_FIELD_FIGURE_SUMMARY.json`

The C4_1 material input is its saved `endpoint_pose`, with the two fixed known-support worlds specified by the parent. Both complex epsilon arrays come from frozen `vendor/wslam/world.py:World.coefficient(f=1,n=384)` on a 401 × 276 display grid over `[-2,2] × [.25,3] lambda0`. Re/Im base and endpoint panels share their respective scales. Signed differences share a diverging palette; an identically zero difference uses one zero color bin with a zero-only colorbar. The NPZ retains the true complex coefficients and differences. This is the narrow TX0/RX8 known-support mechanism fixture; there is no inferred unknown-wall material image.

The currently recorded C4_1 endpoint pose is unchanged at `(.2,2)`. The script checks the saved NPZ pose exactly at runtime and adds the caption **本pair仅wall材料变化；未建立新目标位移多解** (English equivalent in the figure). An unchanged-target pair is not labeled as imaging recovery or target-compensation success. No ambiguity claim for complete 2TX/16RX current data is made by this display.

The field input is the separately fixed `strong_center.npz:fixed_field` on `FIXED_MESH.npz`. Indexing follows the original `Basis(MeshTri(points,triangles),ElementTriP3(),intorder=8)` and `free=setdiff(arange(N),basis.get_dofs().all())`; mesh-order and saved-array-size mismatches fail before display. Saved free scattered values are mapped back to global P3 DOFs, with the recorded homogeneous outer-boundary convention. At P3 nodes inside the material box only, source-0 total field is saved scattered field plus `0.25j*hankel1(0,2*pi*f*distance)`. Since material nodes have `y>=.25`, the `y=0` source singularity is excluded. Magnitude and principal phase have colorbars and coordinates in lambda0. Zero-magnitude phase, if present, is undefined and its display triangles are masked.

The field NPZ preserves original mesh arrays, global P3 DOF locations and free/boundary indices, saved/free and mapped/global scattered arrays, displayed node indices and triangulation, analytic incident and complex total field, magnitude, phase and its undefined mask. Matplotlib triangulation/Gouraud colors are rendering only. The script imports no FEM object, assembly form, `asm`, LU, or physical evaluator and executes no PDE solve. Missing required input files/keys produce `NOT_RUN` and no substitute image. The two cases remain separately labeled and were not chosen from unknown truth or posterior scores.

Wrapper accounting and rendered-image review remain with the parent. No main report or gate file is changed.
