# Fresnel measured-data audit and actual inversion, 2026-11-04

**Status: COMPLETED_EXPLORATORY_MEASURED_INVERSION, 6/6 pipelines. No measured A advantage.** Four original dielectric TM files were freshly obtained from official IOP supplementary links. Two physical targets were inverted at 2 GHz with mature Dual-CIE, Global-Dual and original A-material, followed by the same original physical corrector. Their saved material endpoints are bitwise identical within each target. There is no exact pixel GT, no material NRMSE/IoU, and no claim of calibrated absolute permittivity recovery.

## Access and provenance

The [official Fresnel landing page](https://www.fresnel.fr/3Ddatabase/) and canonical [IOP 2001 supplementary page](https://iopscience.iop.org/article/10.1088/0266-5611/17/6/301/data) returned HTTP 200 to Python requests on 2026-11-04. The latter exposed public signed links for original files. All four dielectric TM files downloaded successfully, with original ASCII headers and data. Earlier web-wrapper navigation to the article reported a restricted URL; this did not reflect the actual Python download route. No challenge was solved or access bypassed. The already retained September originals for the two 1–8 GHz files were byte-identical to the fresh downloads; no checksum operations were performed.

Files in `data_public/fresnel`: `dielTM_dec4f.exp` (550707 bytes), `dielTM_dec8f.exp` (1101084), `twodielTM_4f.exp` (550689), `twodielTM_8f.exp` (1101067). Receipts: `delegated/public_literature/FRESNEL_FRESH_DOWNLOAD.json`, `data_public/fresnel/PROVENANCE.json`. `src/fresnel_fetch.py` repeats the official download; expiring signed URLs are regenerated from the canonical page.

The official landing permits scientific use. The supplementary page reserves underlying authors' rights; an unrestricted redistribution license was not found. Keep raw measurements, publisher/author PDFs and downloaded full HTML out of the project's permissively licensed public source release; distribute the acquisition script, source links, audit, derived numerical outputs and figures instead. In particular, per-run `endpoint.npz` contains `measured_Y` and is local-only; `RESULTS/fresnel_public_derived.npz` is the small publishable array bundle with only reconstructed m, model predictions, masks and centers.

## Primary measurement contract

Source: Belkebir and Saillard (2001), [primary descriptor](https://www.fresnel.fr/perso/belkebir/Articles/Ip01Introduction_Belkebir.pdf), DOI [10.1088/0266-5611/17/6/301](https://doi.org/10.1088/0266-5611/17/6/301), §3–6, read alongside original file headers.

| Item | Resolved meaning |
|---|---|
| Geometry | Source radius 0.720 m, receiver radius 0.760 m, each reported ±0.003 m. Source/view angle `(index−1)×10°`, receiver angle `(index−1)×5°`. |
| Acquisition | 36 views; 49 measured receivers per view. Receiver/source relative angular sector 60–300° inclusive. Missing 23/72 receiver indices per view remain absent, never zero measurements. |
| Polarization | TM: electric field parallel to long-cylinder axis. The finite scalar 2D E model is compatible with this idealization. TE/metal files are excluded. |
| File | Seven columns: source index, receiver index, physical frequency in GHz, Re/Im total E, Re/Im target-absent incident E. Scattered E = complex total − complex incident. |
| Convention | Published time dependence exp(+iωt). Conjugate both complex scattered and incident fields to reuse inherited outgoing H1 / exp(−iωt) operators. |
| Frequencies | `*8f`: 1–8 GHz in unit steps; `*4f`: 4,8,12,16 GHz. Only 2 GHz from each 8f file is inverted here; other frequencies are not extra independent objects. |
| Nominal target | Dielectric circular cylinders, radius 15 mm, reported real εr≈3±0.3; single off-center and two-cylinder configurations. Placement and exact material/loss are uncertain. |
| Units and calibration gaps | Absolute E units, detailed horn antenna factors and interior incident field are not supplied. No noise covariance is supplied. The descriptor explicitly calls for estimating or optimizing the incident field inside the target region. |

`RESULTS/fresnel_ingestion.json`: all four files have seven finite columns, zero duplicate source/receiver/frequency keys. The 8f files contain 14112 rows each, the 4f files 7056. At 2 GHz each measured target has 1764 valid complex channels and 828 absent channels on the 72×36 index grid.

## Explicit approximation and frozen split

`src/fresnel_adapter.py` reuses `research/research_z_alg_r4_diversity/src/physics.py::build_operators` and `research/research_fbg_r3_experiments/src/physics.py::MaxwellModel`, including the retained dense LU and all sources as RHS. It does not replace nonlinear Maxwell scattering with Born fitting.

The interior horn illumination is approximated by a 2D cylindrical source. For each source independently, its complex amplitude is determined only from **target-absent** incident E at the opposite, on-axis receiver. The cylindrical Green ratio transports the measured source-to-receiver field to the target cells; no object-scattered residual or material GT sets this scale. This approximation does not reproduce the full horn beam/antenna pattern and is not independently calibrated. It is part of the scientific model mismatch, not a claim that the supplied data give exact interior E.

Use 24×24 disk-pulse cells over a 0.16×0.16 m square, radius `0.45×spacing` exactly as inherited. The pulse disks cover about 63.6% of a cell's square area; reconstructed contrast is therefore the inherited disk-pulse parameter, not a certified continuous pixel εr. No mesh robustness or material accuracy is inferred from this pilot.

Native solver receiver S is shared over sources. Rather than rewriting its missing-data algebra or filling absent values, train on sources 0–80° (nine sources) and their common receiver indices 29–61 (angles 140–300°). Within that common sector reserve every fifth receiver index starting at 29. Training has 26×9=234 actual channels; 1530 actual channels are withheld, including 63 receiver channels at trained sources and all 1323 channels from the other 27 sources. Each held-out source uses only its target-absent incident calibration. Held-out object scattering never sets gain, selects endpoints, stops iterations, or tunes parameters. Training aperture is deliberately limited; extrapolation is correspondingly harder.

Every pipeline freezes its actual config before inversion in `RESULTS/fresnel_<case>_<arm>/FROZEN_CONFIG.json`: common host .1, contrast box [.02,20], DPCS threshold .001, inherited prior, 500/500/500/1500 auxiliary stages with 3/6/10/12 Fourier masks, full original material derivatives, then original L-BFGS-B physical corrector capped at 1000 accepted iterations. This is a new measured pilot contract, not the frozen 3000-corrector R6 synthetic benchmark. All six use exactly the same measured acquisition and objective. CPU is single-thread BLAS; no GPU.

## Actual results

Relative errors are norms of complex physical prediction minus measured scattered E, divided by the measured E norm in that split. They are not noise-normalized statistical discrepancies.

| Target | Pipeline | Train | All withheld | Held-out sources | Held-out receivers | CPU s | Corrector stop |
|---|---|---:|---:|---:|---:|---:|---|
| single | Dual-CIE | .141931 | .267555 | .279765 | .164767 | 11.118 | projected gradient, 140 accepted |
| single | Global-Dual | .141931 | .267555 | .279765 | .164767 | 13.907 | same |
| single | A-material | .141931 | .267555 | .279765 | .164767 | 34.440 | same |
| twin | Dual-CIE | .100624 | .234442 | .245000 | .133555 | 35.013 | 1000-iteration cap |
| twin | Global-Dual | .100624 | .234442 | .245000 | .133555 | 37.919 | same |
| twin | A-material | .100624 | .234442 | .245000 | .133555 | 97.039 | same |

Total inversion CPU = 229.436 s, including data/operator setup, auxiliary inversion, corrector and held-out prediction, excluding acquisition/retrieval/plotting and final artifact serialization. Per-run wall time, auxiliary operations, material evaluations, LU/RHS and training/prediction costs are in each `RECEIPT.json`; aggregate `RESULTS/fresnel_summary.json`. Additional audit evaluations and retrieval costs are outside this 229.436 s total and are not represented as zero.

Single/twin discrete constitutive state residuals are 2.14e−15 / 1.86e−15. This checks the algebraic full-wave solve, **not** discretization, horn fidelity or fit acceptance. Twin correctors reached their declared cap; preserve `success=false` and do not call them converged.

All three saved m arrays are bitwise equal within each target. Their first corrector objective is exactly the same as the common uniform host physical objective (single .3374385351; twin .3069630488), documented in `RESULTS/fresnel_auxiliary_selection_audit.json`. The native auxiliary solver selects the lowest original-physical-score candidate including its initial incumbent. These observations are consistent with rejecting the auxiliary stage endpoints and retaining the host; they provide **no evidence that any auxiliary method improved the measured inverse trajectory**. Original per-stage candidate material arrays were not retained by this first pilot; this recording limit is explicit. The common physical corrector accounts for the observed endpoint, while A adds overhead. Do not interpret identical final fits as three independent successes or evidence that A helped.

Figures: `FIGURES/fresnel_material_endpoints.png`, `FIGURES/fresnel_complex_predictions.png`. Both reconstructions attain the contrast upper bound 20 and concentrate contrast into a few bright cells; this compensatory, limited-aperture result does not establish nominal dielectric shape or permittivity recovery. No nominal geometric image is shown as GT. No numerical material-accuracy threshold is invented.

## Reproduce

Run `src/fresnel_fetch.py` with the project's Python, then run `experiments/fresnel_invert.py --case dielTM_dec8f --arm cie_dual` and analogous case/arm combinations (`twodielTM_8f`; `global_dual`, `a_material`). Output directories must not exist, preserving prior attempts. Use the existing project environment and one BLAS thread; no new framework or installation is required.

The present conclusion is limited: real measurements were obtained, adapted, inverted and predicted on withheld channels. Residuals remain 23–27% on the broad withheld split under an uncalibrated horn approximation. No claim of true material recovery, unknown array calibration, multi-frequency LSL, mature independent NIE, vector 3D or algorithm superiority is supported by this pilot.
