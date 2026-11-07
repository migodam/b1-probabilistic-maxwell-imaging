# E0 absolute energy audit

**COMPLETE:** 3264/3264 model-observations: all 64 low-amplitude mothers (clean/noisy), all 32 existing A2.5 mothers (clean/noisy), all 17 confirmation models (three seeds plus two cascades). No training, no R4 test truth access.

- Main rows below are **actual nonlinear Maxwell-current errors**, obtained by solving the unchanged frozen Model at the historical raw material and final top1. Auxiliary-head discrepancies and material-coordinate Green energies remain separate quantities.
- For clean low-amplitude seed-7 results, final internal-strong19 and receiver-strong9 error increases in all five arms; D509 decreases. For noisy results, final internal-strong19 decreases in all five arms, while receiver-strong9 still increases. Visible48 and D509 decrease. Thus the absolute-current results do not support a blanket strong-error improvement statement.
- The raw→final pipeline is recorded directly, including update energy and signed cross term. It includes finite-candidate correction and selection and does not isolate a TRF effect. All 1278 unclipped/clipped differences are float32 lower-bound roundoff only (max Linf 4.47035e-10, max RMS 3.56345e-10; zero above 1e-7), so no material clipping confound is evidenced. Every band records ΔE = update energy + signed cross term; a large update alone does not establish improvement.
- Receiver true-error and observation-misfit have different references. data_residual_relation.csv records E_true = E_misfit + E_noise + signed cross, using frozen physical residual receipts and stored historical GT current, with zero additional LU. Clean GN can fit far below the .001 NN correction gate; stopping accuracy can partly account for increased true receiver-field error. This audit does not establish causality.
- Complex Green projections of real material error are labelled material_complex_Green_embedding; they are not material observability or real-Jacobian null spaces. The deterministic real-J diagnostic is separate (32 observations, 16 mothers), never promoted to the full test population.

## Actual-current final minus GN absolute squared error (seed 7, one pass)

|noise|arm|internal strong19|receiver strong9|visible48|D509|
|---:|---|---:|---:|---:|---:|
|0.0|v_image|2.04815e-05|4.80975e-07|-0.000106007|-0.00942719|
|0.0|v_physics|1.17933e-05|7.53327e-08|-8.05809e-05|-0.00880125|
|0.0|j_full|2.64283e-05|3.68571e-08|-1.50537e-05|-0.00463261|
|0.0|j_null|1.08822e-05|3.63108e-08|-1.39937e-06|-0.00646363|
|0.0|j_null_current|1.14858e-06|4.87323e-07|7.75773e-07|-0.000224655|
|0.001|v_image|-0.00189696|9.34892e-06|-0.0690393|-0.712985|
|0.001|v_physics|-0.00234746|1.39804e-06|-0.0753137|-0.844633|
|0.001|j_full|-0.00223286|2.36457e-06|-0.0743625|-0.833642|
|0.001|j_null|-0.00205677|6.3515e-06|-0.0517053|-0.784718|
|0.001|j_null_current|-0.00201129|5.99748e-06|-0.0265526|-0.804989|

Negative ΔE is improvement. All entries are absolute squared norms, never shares. DECISIVE_TRUE_CURRENT_SEED7.csv includes raw/final and mother-paired 95% bootstrap intervals. The full bootstrap file keeps seeds, noise and low/stress distinct; grouped_absolute_energy.csv also separates amplitude.

Fixed geometry: P_in=19/557, S=9/39/528, visible48; B67 is the union of internal strong19 and receiver visible48, D509 its orthogonal complement. Sequential I→O and O→I, true numerical intersections (weak/weak rank20 and weak/null rank509), and the numerical-intersection remainder are separate columns.

Cost: 177.4 s main audit; 6532 LU, 79920 complex RHS; one CPU worker with BLAS=1. Full-VIE max equation defect 1.73e-15. All 3264 frozen final top1 outputs were available. No missing E0 absolute-energy cases. Exact unclipped-raw observed-misfit relations are formally NOT_RUN for 1278 float32 lower-bound-roundoff cases; actual verified raw_1 residuals are retained in a separately labelled column. No material clipping effect or missing absolute-current energy follows from these formal mismatches.

Files: per_case_energy.csv, grouped_absolute_energy.csv, mother_paired_bootstrap.csv, DECISIVE_TRUE_CURRENT_SEED7.csv, data_residual_relation.csv, real_J_diagnostic.csv, real_J_diagnostic_delta_summary.csv, GEOMETRY.json, RECEIPT.json. Figures: absolute_raw_final.png, noisy_D509_paired_delta.png, true_current_strong_paired_delta.png. Final scientific interpretation remains with the parent review.

Absolute before/raw/final current magnitudes are in absolute_before_raw_final.csv and figures/e0/absolute_before_raw_final.png. Clean receiver-strong9 GN energy is 2.08393e-13, final means 3.63110e-8 to4.87323e-7; noisy GN energy 1.25854e-7, final means 1.52389e-6 to9.47477e-6. These increases are large relative to very small GN strong-error baselines.
