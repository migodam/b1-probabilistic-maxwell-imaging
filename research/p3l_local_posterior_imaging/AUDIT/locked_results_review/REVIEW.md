# Final locked-summary arithmetic and reporting audit

2026-10-22. Independent inspection of saved JSON, frozen configurations, small ledger snapshot and actual job receipts only. No GT files, locked observation/result arrays, model execution, training or physics. This is an arithmetic/reporting audit; final scientific interpretation remains with the parent.

**No blocking arithmetic discrepancy found.** Eight concrete findings follow.

1. **Complete frozen cohort and zero fallback use.** Risk, online-imaging and Full-IS summaries each contain exactly the same ordered 128 `SBI_locked_0000…0127` identities, with all128-retained flags true. All 128 risk rows are evaluated; missing-truth list is empty. All seven methods have zero recorded raw failures and zero fallback selections; both model-MC-error diagnostics are available for every parent. Imaging, IS, generation and isolated evaluator receipts are `COMPLETE`, error null. This means complete output collection, not universal optimizer stationarity.

2. **Every reported mean and paired comparison reproduces.** Recomputed all methods × full-grid/ROI × MSE/MAE cohort means and all declared baseline paired comparisons directly from the saved per-parent JSON risks. Using seed 2026102299, 2,000 paired parent bootstrap draws, the maximum absolute difference is exactly zero for means, difference CIs, relative reductions and relative-reduction CIs. Frozen/evaluation seed, repetitions and selected-density identity agree. The preregistered primary probit64 baseline is present.

3. **Primary direction and units must remain explicit.** Model×strong-Gaussian mean MSE is **0.061502116830457**, versus probit64 coordinate point **0.011030917324877**. Difference is **+0.050471199505580**, paired 95% CI **[0.042787339485809, 0.058286498161523]**. Ratio is **5.5754308567**. Reported relative reduction **−4.5754308567** means **−457.5431% reduction / +457.5431% increase**, not −4.575%. Its CI is [−5.4956949915, −3.7687812090] in fractional units. These are the actual frozen estimator-risk numbers; they are not a comparison against a qualified Full-posterior mean.

4. **Point-estimator optimization status limits MAP wording.** All 128 native96 outputs are saved, but only **79/128** have reported scaled KKT ≤0.001; range is 1.148e−7 to 0.4885. Probit64 scaled gradient is ≤0.001 for **115/128**; range is 1.308e−6 to 0.08368. Cheap16 has **0/128** below 0.001. Therefore call these bounded point estimators/coordinate points and preserve stationarity metadata. Probit64 optimizes density in x=(z,eta), so it is not native-uniform MAP or a material Bayes action. No global-MAP certificate is established even where local stationarity is small.

5. **Weak-coordinate scores remain limited generative diagnostics.** All 128 calibration parents are assessed. Proper weak8/probit2 NLL mean **14.3094417881**, parent-level SE **0.2091279815**, empirical component 90% coverage **0.859375–0.9296875**. NLL, SE, coverage and each PIT histogram exactly reproduce from saved diagnostics; every component histogram totals 128. Ten correlated components do not create 1,280 independent parent units. Scores/ranks/intervals do not certify Full posterior, material credible coverage or calibration of the strong56 reference. Model-reference image MCSE measures finite IID model draws only, not posterior error or approximation bias.

6. **Full-target weighting did not produce a qualified reference.** All 128 IS parents completed their fixed 128 draws, with explicit finite-SNIS bias and `posterior_qualified=false`. ESS is **min 1, median 1, max 1.9970604791**; **128/128 ESS<2**. Maximum normalized weight ranges **0.519426–1**. Keep “Full-target corrected” separate from “qualified posterior”: exact target weights do not remove finite-sample degeneracy or SNIS bias. IS image full-grid MSE is **0.024000911416106**; its statistical status is still unqualified.

7. **Actual D costs reconcile and remain below cap.** Per-phase JSON deltas sum to imaging **543,912 RHS / 16,684 states / 16,684 LU**, matching both the imaging closure counter and its cumulative-ledger delta. Breakdown: cheap16 **67,896 RHS**, independent native96 **263,192 RHS**, incremental probit64 **212,824 RHS**. Full IS is exactly **131,072 RHS / 16,384 states / 16,384 LU**. Locked generation is **1,024 RHS**; isolated offline evaluator adds zero RHS. Thus the complete locked sequence uses **676,008 new RHS**. Current cumulative D snapshot is **1,040,616 RHS**, leaving **459,384** under its 1,500,000 cap; GPU seconds remain zero. Imaging closure wall is **417.5351 s**, IS **406.8126 s**, evaluator **2.8064 s**, generation **1.5843 s**. Do not describe probit64's 212,824 RHS as standalone total: its initialization also uses native96, giving **476,016 RHS** combined for that point baseline. Shared setup and historical model-selection/training costs remain separate disclosures.

8. **All seven risk actions and cohort exposure remain visible.** Full-grid mean complex MSE is shown below. All are actual complete-cohort risk actions, with raw/effective estimates identical because no fallback occurred. This locked cohort has now been evaluated and must remain labeled exposed for any subsequent tuning; additional C proxy experiments cannot retroactively become fresh held-out evidence on it. No final GO/NO-GO, posterior validity or novelty judgment is supplied by this audit.

| Fixed action | Mean complex MSE |
| --- | --- |
| cheap16 | 0.012841728810002 |
| independent native96 | 0.012493812169956 |
| probit64 coordinate point | 0.011030917324877 |
| model×strong Gaussian reference mean | 0.061502116830457 |
| conditional-mode-plugin reference mean | 0.056516952258972 |
| original-prior material mean | 0.175420681054668 |
| finite Full-target IS, unqualified | 0.024000911416106 |

[Compact recomputation and receipt evidence](REVIEW.json). Original summaries, configurations, arrays and receipts were not modified.
