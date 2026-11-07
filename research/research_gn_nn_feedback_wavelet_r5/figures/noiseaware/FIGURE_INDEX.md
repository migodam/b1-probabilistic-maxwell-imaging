# Phase II figure index

**Scope: DEVELOPMENT VAL used for checkpoint/model selection; no independent TEST.**

Only the five primary models are included. All pictured image/mother identities are absent TRAIN; classes/fonts may be seen. The extra noise-ball mechanism arm is excluded. Inputs: common noise-aware GN anchors, the same declared TRAIN128 / VAL78 mother cohort, separate from-scratch fits.

Raw full cohort: **780 arm-observation outputs**, **351/780 physically qualified**. Separately corrected matched subset: **360 outputs (72 per arm)**, **360/360 qualified**. These denominators and stages must remain separate. Physical qualification does not establish material recovery.

## Report figures

- [01_actual_architecture.png](01_actual_architecture.png): Actual architecture; common noise-aware GN anchor; fixed prior / free control; all seven bands.
- [01_actual_architecture.pdf](01_actual_architecture.pdf): Architecture PDF.
- [02_train_val_curves.png](02_train_val_curves.png): Separate TRAIN / development VAL histories; all five primary arms completed 12 epochs.
- [03_prespecified_val_reconstructions.png](03_prespecified_val_reconstructions.png): Main raw reconstruction panel; GT + common GN + all five models; prespecified clean/noisy A=.2 examples.
- [04_absolute_haar_error_energies.png](04_absolute_haar_error_energies.png): Absolute full-Haar error-band contributions summing to voxel MSE.
- [05_quality_residual_actual_cost.png](05_quality_residual_actual_cost.png): Raw material quality, actual physical residual/qualification and charged full-LU cost.
- [06_depth_diagnostics.png](06_depth_diagnostics.png): Combined depth subset; six observations per arm at 1,2,4,8,16.
- [06_depth_diagnostics_clean.png](06_depth_diagnostics_clean.png): Clean depth subset; three observations per arm.
- [06_depth_diagnostics_noisy.png](06_depth_diagnostics_noisy.png): Noisy depth subset; three observations per arm.
- [07_raw_vs_corrected_matched_subset.png](07_raw_vs_corrected_matched_subset.png): Separate raw / corrected outputs on exactly matched360-row subset.
- [08_gn_input_mode_raw_pipeline_comparison.png](08_gn_input_mode_raw_pipeline_comparison.png): Phase I / II observed raw pipeline comparison; matched identities; no isolated causal attribution.
- [09_no_field_material_vs_physical_gap.png](09_no_field_material_vs_physical_gap.png): No-field ablation: existing VAL material error versus actual nonlinear physical residual.

## Complete raw gallery

All 156 ordered VAL observations appear on 26 pages, with all five models and failures retained. Fixed display [.02,.55]; pixel saturation counts and exact manifest exposure flags are printed.

- [Page01](gallery_01.png)
- [Page02](gallery_02.png)
- [Page03](gallery_03.png)
- [Page04](gallery_04.png)
- [Page05](gallery_05.png)
- [Page06](gallery_06.png)
- [Page07](gallery_07.png)
- [Page08](gallery_08.png)
- [Page09](gallery_09.png)
- [Page10](gallery_10.png)
- [Page11](gallery_11.png)
- [Page12](gallery_12.png)
- [Page13](gallery_13.png)
- [Page14](gallery_14.png)
- [Page15](gallery_15.png)
- [Page16](gallery_16.png)
- [Page17](gallery_17.png)
- [Page18](gallery_18.png)
- [Page19](gallery_19.png)
- [Page20](gallery_20.png)
- [Page21](gallery_21.png)
- [Page22](gallery_22.png)
- [Page23](gallery_23.png)
- [Page24](gallery_24.png)
- [Page25](gallery_25.png)
- [Page26](gallery_26.png)

## Auditable tables and scope

- [RAW_NN_FULL_COHORT_SUMMARY.csv](RAW_NN_FULL_COHORT_SUMMARY.csv)
- [CORRECTED_MATCHED_SUBSET_SUMMARY.csv](CORRECTED_MATCHED_SUBSET_SUMMARY.csv)
- [NO_FIELD_MATERIAL_PHYSICS_GAP.csv](NO_FIELD_MATERIAL_PHYSICS_GAP.csv)
- [INPUT_MODE_RAW_PIPELINE_COMPARISON.csv](INPUT_MODE_RAW_PIPELINE_COMPARISON.csv)
- [HAAR_ERROR_BAND_ENERGIES.csv](HAAR_ERROR_BAND_ENERGIES.csv)
- [DEPTH_PLOT_SCOPE.json](DEPTH_PLOT_SCOPE.json)
- [FIGURES.json](FIGURES.json)

## Raw qualification counts

| Primary model | Clean | Noisy |
|---|---:|---:|
| Pixel closed | 78/78 | 77/78 |
| Haar closed | 0/78 | 58/78 |
| db2 closed | 4/78 | 70/78 |
| Haar no field loss | 0/78 | 1/78 |
| Haar free control | 0/78 | 63/78 |

## Rendering review

One complete Phase II plotting invocation was run with the existing Gaussian Python and one BLAS thread. Viewed the main reconstruction panel, no-field physical-gap plot, combined/clean/noisy depth figures, architecture, raw/corrected comparison and GN-input-mode comparison. No display defect required a rerender. Remaining gallery pages were generated but not individually viewed.

The absolute Haar band-sum / voxel-MSE maximum absolute difference is 1.862645149230957e-09. All 576 coefficients are retained.

Depth plots cover only the declared six-observation subset per arm, not all 156 VAL observations. Free-cascade energy curves do not have a shared fixed-energy guarantee. The Phase I / II figure compares full pipelines whose GN inputs and independently trained weights differ; it does not isolate wavelet benefit or a loss-only effect.

No training, inference, physical solve or new experiment was run for plotting. TRAIN truth and TEST data were not opened. Phase I figures were not rendered or edited.
