# Authorized physical validation and exact-grid vector gallery

All results below were evaluated only after the joint offline authorization and the parent-authorized physics seal. The measured evaluator completed all 56 stored views from 14 tasks; no candidate, seal, prior, physics model or online selection was changed. The conclusions concern predictive fit on the actual heldout measurements. There is no Fresnel pixel material GT, so no material NRMSE/IoU or material-recovery superiority is reported.

## Measured predictive fit and correction

Each value is a relative complex Frobenius residual in percent. Training uses 1053 actual channels; heldout uses 711 actual channels (441 on withheld sources and 270 on withheld receiver channels). The frozen original full E0/source-port rows were verified against the training package before validation. No absent channel was zero-filled, and no gain or new source normalization was fitted.

Raw is the paid proposal incumbent selected by original training Phi; primary is the final original-Phi selection. For direct, raw is explicitly the original initialization, so its raw-to-corrected change represents direct physical optimization rather than an auxiliary proposal plus corrector.

### dielTM_dec8f

| Method | Raw train | Raw heldout | Corrected train | Corrected/primary heldout | Secondary heldout | Primary Phi |
|---|---:|---:|---:|---:|---:|---:|
| cie | 90.7787% | 90.8347% | 8.7643% | 9.1500% | 14.9231% | 0.003841090 |
| dual_cie | 90.7787% | 90.8347% | 8.7643% | 9.1499% | 14.9231% | 0.003841089 |
| global_dual | 90.7787% | 90.8347% | 8.7643% | 9.1499% | 14.9231% | 0.003841086 |
| dr | 9.7378% | 10.1152% | 8.7643% | 9.1498% | 11.6400% | 0.003841078 |
| ap | 15.2303% | 15.3294% | 8.7647% | 9.1505% | 14.8150% | 0.003841383 |
| event_integer | 10.0755% | 10.3700% | 8.7643% | 9.1499% | 9.8522% | 0.003841052 |
| direct | 90.7787% | 90.8347% | 8.7643% | 9.1499% | 14.9231% | 0.003841045 |

### twodielTM_8f

| Method | Raw train | Raw heldout | Corrected train | Corrected/primary heldout | Secondary heldout | Primary Phi |
|---|---:|---:|---:|---:|---:|---:|
| cie | 78.3539% | 78.5334% | 9.1400% | 9.8575% | 15.1846% | 0.004177867 |
| dual_cie | 78.3539% | 78.5334% | 9.1400% | 9.8575% | 15.1846% | 0.004177867 |
| global_dual | 78.3539% | 78.5334% | 9.1400% | 9.8575% | 15.1846% | 0.004177867 |
| dr | 10.5043% | 11.1660% | 9.1378% | 9.8552% | 11.3927% | 0.004175912 |
| ap | 15.0723% | 15.5804% | 9.1379% | 9.8552% | 14.8743% | 0.004175996 |
| event_integer | 9.4533% | 10.1335% | 9.1380% | 9.8551% | 10.7237% | 0.004176083 |
| direct | 78.3539% | 78.5334% | 9.1400% | 9.8575% | 15.1846% | 0.004177867 |

After correction, the seven primary fits cluster around 9.150% heldout error for the single specimen and 9.855–9.858% for the twin. The largest raw differences largely disappear under original physical correction. All 14 corrected views equal their primary selected material in this saved experiment. This is a bounded two-specimen predictive result, not evidence of material accuracy or universal method superiority.

The three native measured CIE variants each completed 300 material updates with 35–36 paid original physical audits, including native stage endpoints. Nevertheless, their raw original-Phi incumbent remained the uniform initialization for both specimens. The saved histories show those audits were eligible/legal; this is not explained by missing periodic audits or an early import failure. Their final predictive fit comes from original physical correction. The poor auxiliary candidates and original raw initial states remain preserved. No algorithm was revised to rescue these outcomes.

DR, AP and integer-event proposals produced distinct raw predictive fits. All 14 combinations of method and specimen improve heldout fit under their recorded physical correction, including direct optimization from initialization. Several correctors hit their declared process CPU allowance; numerical convergence must not be inferred from the final residual alone. Exact stop records, D/R decomposition and original online costs are retained in CORRECTOR_COMPARISON.csv and the sealed receipts.

Every stored secondary candidate exists, with no fallback used. Secondary selects minimum prior among paid candidates in the fixed 15% training discrepancy ball; it is not selected by heldout error. The table retains its different predictive result. Membership of that ball is a model-discrepancy diagnostic, not a CN99 certificate: actual measurement noise covariance and receive-horn response remain unresolved.

## Actual validation costs and scope

The 56 requested view references contain 27 distinct material groups in total (15 single-specimen and 12 twin-specimen groups); sharing is restricted to exact material equality within one calibrated object. Paid shared validation used **27 LU and 972 source RHS**: one retained LU serves all 36 original source columns for each unique material. Forward-only time was 0.218255 CPU s/0.229756 wall s. Input/model setup was 0.120713 CPU s; imports 0.270692 CPU s. Whole external evaluator cost was 0.822801 CPU s/0.947316 wall s. These are offline validation costs on the local Apple CPU, separate from the online i9 method costs. CPU timing across these hosts must not be pooled as a speed comparison. Per-reference cost shares are allocations; the integer physical work is counted once in UNIQUE_VALIDATION_COSTS.json.

No new physical validation failures occurred. Calibration/mask checks, access times after authorization, measured predictions and all 56 rows are retained under analysis/fresnel_holdout. The fixed calibrated forward model remains an approximation to the laboratory setup; a roughly 9–10% heldout predictive residual does not identify dielectric material uniquely.

## Exact full-vector 3D gallery

All **32 PNG/PDF pages** completed, covering four conditions, seven methods and initial/raw/primary/corrected views: 112 method-view panels plus 32 repeated GT reference panels. Every saved panel is available and none has out-of-box cells. Four conditions represent three independent objects; the 8³→6³ frame is a repeated geometry. Rendering used the exact saved coarse 6³ arrays, fixed chi>1 support, physical color range [0,8], camera 25°/−55°, equal aspect and exact index 3 orthogonal sections at +.02 lambda0. No interpolation, surface refinement or physical solve occurred.

External plot cost was 18.863529 CPU s/19.164855 wall s, with zero additional LU/RHS and zero failed renders in this run. Prior failed plot attempts remain explicitly unknown. The plot manifest and panel table contain every figure ID, source-array reference and PNG/PDF download path.

Visual QA checked the matched-frame primary page 1 and crossgrid-frame corrected page 2. Labels/scales and coarse voxels are readable; the sampled frames visibly retain broadened support and undercontrast relative to the coarse GT. This descriptive visual observation is not a substitute for the parent’s full quantitative scoring or review of the remaining pages. Closed/off-plane pores may remain hidden from one camera; exact orthogonal slices are included alongside it. No continuum-accuracy claim follows from the gallery or from the old bounded Mie controls.

Artifacts: [measured table](../../analysis/fresnel_holdout/MEASURED_HOLDOUT.csv), [correction comparison](../../analysis/fresnel_holdout/CORRECTOR_COMPARISON.csv), [unique validation cost ledger](../../analysis/fresnel_holdout/UNIQUE_VALIDATION_COSTS.json), [3D figure manifest](../../analysis/vector_volumes/FIGURE_MANIFEST.csv), [3D panel table](../../analysis/vector_volumes/FIGURE_PANEL_TABLE.csv), [plot cost receipt](../../analysis/vector_volumes/PLOT_COST_RECEIPT.json). Final scientific and publication acceptance remain with the parent.
