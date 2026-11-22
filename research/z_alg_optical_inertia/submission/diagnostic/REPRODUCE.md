# Diagnostic reproduction

Run the following from the Inv_SLAM root using Python with NumPy, SciPy, python-flint and Matplotlib. The replay disables multi-threaded BLAS. All programs write their results to this diagnostic directory. When replaying a sealed run, preserve the original receipts first if the original machine timing is needed.

```sh
python -B research/z_alg_optical_inertia/submission/diagnostic/run.py
python -B research/z_alg_optical_inertia/submission/diagnostic/audit.py
```

- The plot is already included as [PDF](figure_gn_diagnostic.pdf) and [PNG](figure_gn_diagnostic.png). An authorized workspace can recalculate the R10 material scores by supplying the original reference using `audit.py --r10-reference FILE`, and redraw the two reference panels using `figure.py --r10-reference FILE`.
- [Frozen design](FROZEN_DESIGN.md) fixes the cases, objective and stopping rule. The online program reads the two saved noisy port arrays and paired operators. Its outputs contain every trial material, full complex port prediction, current, acceptance, and stop. It does not read the reference materials, the cap, or the certificates.
- R10 inputs are `research/z_alg_secant_tasks/highdim/public_vie_inputs.npz`, `research/z_alg_secant_tasks/highdim/results/r10_twod_0001.npz` and `research/z_alg_secant_tasks/highdim/sector_inertia/certificate_u5.json`. The readout checks the noisy measurement against the original MAX-STRUCT R10 `data_public/twod/r10_twod_0001.npz`. The cap applies to all nonnegative materials, including the zero entries permitted in this diagnostic. R10 originally used a positive minimum in its optimization contract. The original ring reference remains sealed. The archive includes the previously recorded scalar scores and figure while the two pointwise forward audits use the assembled public operators.
- Native inputs are `synthetic/public/fine_geometry.npz`, the frozen observation `synthetic/public/optfresh_01__pw180_10x20__n02.npz`, and the interval cap and feasibility receipt `results/theory/fine_feasibility.json`. After auditing the reconstruction, the readout opens `synthetic/evaluation/optfresh_01.npz` for scoring. There is no coarse-grid move or model-mismatch inflation in either diagnostic. The additional native radius guard is the existing allowance for generating the exact discrete observations with a finite-tolerance numerical solve.
- The two terminal materials and the first three nontrivial native updates are read back as exact binary64 arrays for peak checks. The R10 residual audit forms an Arb enclosure of an active dense inverse defect and uses a corrected forward. Native replays use the same exact Fourier, active-LU, inverse-error and corrected-port implementation as the previously published native feasibility audit. The chronology includes the two early native updates whose peaks lie below the cap. The terminal native candidate is the sole native endpoint.
- [Interval and offline receipt](interval_and_offline.json) contains bounds, time, operations, and scores. The ring interior is the set of original material values above its 0.1 host. The slotted-block interior is its nonzero support. Ground truth was exposed in previous work. Any image display is an offline explanation of these diagnostic observations.

To recompute the cap-wide proof bounds from the frozen public inputs, run the original validators from the ZIP root.

```sh
python -B research/z_alg_secant_tasks/highdim/sector_inertia/certify.py --case r10_twod_0001 --cap 5 --mu 0.001 --output r10_cap5_replay.json
python -B research/z_alg_optical_inertia/fine/verify_fourier.py --geometry 180
```

Existing study and public source entry point [Z-ALG optical-inertia project](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/z_alg_optical_inertia/).

## Measured costs and material readout

The run used one CPU thread on an Apple M4. Both predictions were sealed before the offline metrics were evaluated.

| Case | Optimization wall / CPU | Full forward evaluations | Full Jacobians | Point interval checks wall / CPU | Offline global / target relative L2 | Background RMSE |
|---|---:|---:|---:|---:|---:|---:|
| Archived ring | 0.343 / 0.342 s | 8 | 7 | 0.020 / 0.015 s | 1.0061 / 1.0000 | 0.6243 |
| Native slotted block | 3.876 / 3.844 s | 13 | 12 | 74.474 / 73.912 s | 0.9478 / 0.9459 | 0.09572 |

The optimizer timings include loading public inputs, accepted and attempted full-wave steps, and writing the paid-path arrays. Ring iterations used 96 source and 336 receiver-transpose dense right-hand sides. The native run used 120 nontrivial source and 220 receiver-transpose GMRES right-hand sides, 4,684 Krylov iterations, and 5,154 Green FFT applications. The interval timing includes its separate geometry load and independently encloses one ring point and four native points. The native terminal-only interval replay took 16.76 s wall, after shared geometry setup. The ring generating-material offline audit and scoring took a further 1.977 s wall. The native reference uses the existing independently verified feasibility interval. Existing cap certificate costs are reported with their original proof receipts and are outside the reconstruction timings.
