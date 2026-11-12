# Final manuscript factual review

Review scope: parent `RESEARCH_REPORT_ZH.md`, `PAPER_DRAFT.md`, and `CHAT_HANDOFF.md`, checked against the existing source audits, current runner source, saved visual observations, contracts, cost ledger, and existing analysis CSV/JSON. This is factual verification, not a new novelty decision or final scientific acceptance. No parent document, frozen receipt, result, or selector was edited. No GT/material NPZ, new solver, scoring pass, optimization, or external literature search was used. The compact statistics below are ordinary aggregation of already-written metadata with the standard library.

## Findings for parent action

| Priority | Location | Finding and concrete clarification |
|---|---|---|
| P2 — factual visual ambiguity | `RESEARCH_REPORT_ZH.md:71` | The grouped sentence about 0005/0008/0031 can be read as DR and event both fragmenting on all three. On 0008, **DR retains the perforated body and its cavities; event_integer collapses to the sparse arc**. The actual pages are `diagnostic_structure_methods1_page03.png` and `methods2_page03.png`, row 1. State this difference explicitly. The 0005 and 0031 comparisons are consistent with the visual notes. |
| P2 — count clarity, not an arithmetic conflict | `RESEARCH_REPORT_ZH.md:5`; `CHAT_HANDOFF.md:11`; `PAPER_DRAFT.md:358` | The canonical total 962 is correct. `audit/RUN_CONTRACTS_FINAL/RUN_CONTRACTS.csv` has **1,004 metadata rows = 962 canonical runs + 42 recovery aliases**. Add that identity near the canonical count or its contract link so readers do not infer missing runs. The aliases retain recovery provenance and must not add science samples or duplicate numeric CPU. |
| P3 — statistical wording | `PAPER_DRAFT.md:303` | “independent native-schedule M4 reference” is a separately executed **same-16-object** reference, not an independent material cohort. The surrounding text explains this correctly; replacing “independent” with “separately executed same-object” would remove the local ambiguity. |
| P3 — causal strength | `PAPER_DRAFT.md:303,352` | “supports the importance of trajectory” and “can remove favorable continuation” should retain the stated representation/initialization/cost/host qualifications. These controls show that increased capacity is not sufficient under these runs; they do not uniquely identify continuation as the cause of any individual loss. No stronger causal attribution should be added. |
| P3 — nearby baseline qualifier | `PAPER_DRAFT.md:7,277` | “mature CIE variants” is acceptable only with the explicit real-material/common-prior adaptation qualification already present at :234 and in the audit. Consider “source-audited CIE adaptations” at first use. This is a wording improvement; the draft does **not** claim original-author reproduction elsewhere. |

No blocking formula-transcription error or newly discovered baseline-wiring contradiction was found in this review. The parent owns the proofs, novelty judgment, population claims, and final qualitative categories.

## Method/theory insertion

The manuscript block from `## Physical model and lifted geometry` through `## Paid audits, correction and provenance` is an **exact text match**, after trimming outer whitespace, to `delegated/full_current/PAPER_METHOD_THEORY_EN.md`: 17,503 characters. Thus this check found no omission or alteration during insertion; it is not a new independent proof.

The inserted text correctly preserves:

- the `exp(-i omega t)` convention; shared real scalar material; the declared state scale; scalar/vector/source-specific operator separation (`PAPER_DRAFT.md:30–63`);
- legal angle competition, boundary signs, ties, dark/zero-normal cases, and the bounded bilinear rather than Maxwell-exclusive scope (:67–110);
- exact mutual-AP fixed-shadow drift versus the conditional one-step approximate-pair bound; the lack of a future multi-step radius guarantee (:120–145);
- the exact same-shadow witness with **abstract G=-1**, explicitly not a Maxwell self term (:147–164);
- the strict state-ball margin and its a posterior interpretation when the radius was measured after a run (:166–182);
- the invertibility, propagation, Slater/KKT and finite multiplier assumptions of the Maxwell resolvent correction; actual projection multiplier rather than fitted pressure (:184–218);
- the single-beta fixed-current-block scope of the quadratic profile identity, without inventing a common dual-beta G19 Lyapunov objective (:220–232);
- AP/DR’s same first step only when the initial state is already in C (:139–145);
- paid physical selection, late-candidate exclusion, raw/corrected/secondary distinctions, and the incomplete auxiliary/PR-state archive (:254–262).

The title/header nesting in this inserted block is stylistic; it does not change the mathematics.

## Baseline and selection qualifications retained

1. **G19/full spaces.** The three G19 arms retain the same material denominator/readout, PR family, beta stages and prior; their affine current families differ. mf12 already spans full APCS on 24². Early full APCS changes continuation; full current additionally releases DPCS. The draft does not call full current a new NFFT-SOM implementation. Source: `delegated/full_current/BASELINE_SOURCE_AUDIT.md:35–65`; manuscript :230–236,303.
2. **Warm versus cold.** CIE/Dual/A+ start with zero APCS coefficients above JD; G19 projects the initial physical current into its admissible family. A shared initial material does not establish shared auxiliary J/Q or optimizer memory. Manuscript :236,277 makes this distinction.
3. **Two different “dual” controls.** global-PR Dual-CIE is distinct from physical GlobalDual. The latter uses the changing total-field source Gram at **beta=0**, including its derivative, and the common uniform0.1 reference fixes its scale. beta=.2 would be feedback, not this historical GlobalDual. Manuscript :234,238–250 is consistent with `BASELINE_SOURCE_AUDIT.md:7–29` and `src/run_experiment.py:151–165`.
4. **Full-J TRF.** It uses the existing analytic physical residual/Jacobian and the same scalar prior residual. The prior half-squared norm identity is stated correctly; no exact nonlinear Hessian claim is made. TRF differs from direct L-BFGS-B and has N independent material variables; “Full-J” is a derivative matrix, not extra material capacity. Source audit :96–136; manuscript :234,279.
5. **NFFT-equations control.** Ordinary LSIE, JD plus complete correction, source data normalization 1/||Y_s||² and state normalization 1/||JD_s||², zero correction meaning J=JD, joint L-BFGS-B, box and prior are disclosed. It is explicitly **NOT original NFFT-SOM author-algorithm reproduction**. Manuscript :248–250 and :360 retain that qualification.
6. **Classical/measured/vector domains.** The scalar/common-port GlobalDual/NFFT/TRF formulas are not silently asserted as vector or measured implementations (:252). Austria and Fresnel are explicitly qualified adaptations, with model/discretization error and no pixel-truth metrics invented for Fresnel (:344–346). CIE source auditing is not an assertion of complete original-paper recovery reproduction.
7. **Paid selection.** Runner `record` computes original D+R and only updates best for legal/eligible paid candidates (`src/run_experiment.py:79–99`). Common audits use process-CPU buckets after completed updates, retain native opportunities, and may overshoot. Already-paid optimizer trials remain candidates. Material error and auxiliary objectives are not primary selectors. The text does not falsely promise identical audit counts or exact audit times.
8. **Finite completion.** “COMPLETE” process output is not native schedule completion. The reported 120-second truncation and the separate 360-second same-object reference are properly distinguished. Different hosts/libraries are explicitly not a pure budget ablation.
9. **Archive coverage.** The draft and report retain the absence of every CIE/G19 J/Q/PR snapshot. Saved physical states and material pools cannot support exact optimizer continuation or a full PR-trajectory proof.

## Counts and CPU reconciliation

The contract CSV ordinary row count gives:

| Quantity | Count |
|---|---:|
| Fixed canonical runs | 917 |
| Of those: fresh condition-method runs | 828 |
| Of those: exposed method runs | 89 |
| Registered conditional rescue runs, 9 triggers × 5 arms | 45 |
| Canonical runs | **962** |
| Additional physics recovery alias rows | 42 |
| Contract metadata rows including aliases | **1,004** |

The 17 local TRF rows contain 16 fresh references plus one exposed diagnostic. The 60 composed physics rows already contain the 42 recovered physical jobs. The alias campaign is retained for attempt provenance, not added again to the canonical total. Forty-two initial import failures and four interrupted partial attempts are attempts, not new material objects.

The stated 80,241.94 i9 online-scientific CPU seconds and 3,995.17 Apple online-scientific CPU seconds agree with `COST_LEDGER.md:7–10`. They include development/reference work, not only fresh core48. The ledger separately retains i9 controller CPU 31.828125 and initial physics-failure CPU 36.734375, plus setup/validation/scoring/packaging and unknown partial-attempt costs. These numbers must not become a fully measured all-project CPU total or be read as parallel elapsed time. Known failures are charged; unknown work remains unknown. Parent-plus-rescue costs retain the parent G19 work for all 48, including 39 bypasses. Counterfactual uniform-DR portfolios have sums of existing job CPU but no newly measured conditional deployment wall time, as the manuscript correctly says.

## Existing finite-discrepancy material distances

Source: existing `analysis/feasible_branches/PER_OBJECT_BRANCH_SUMMARY.csv`, `FEASIBLE_MATERIAL_PAIRS.csv`, and `RECEIPT.json`. No materials, observations, GT, model or forward calculations were loaded in this review.

RMS distance is ||m1−m2||₂/sqrt(N), in susceptibility units. Box-normalized distance divides this by the public box width. **Median across per-object maxima** and **median across all pairs** have different meanings; the latter weights objects according to their number of retained pairs.

| Scope and statistic | Eligible summary/pair count | RMS material distance | Public-box normalized |
|---|---:|---:|---:|
| Core48: median of per-object maxima | 43 objects with pairs /48 | **1.341554** | **0.067145** |
| Core48: largest pair | 17,397 pairs, case r13_fresh_0024 | **3.537434** | **0.177049** |
| Core48: median across all pairs | 17,397 | 0.625469 | 0.031305 |
| All recorded domains: median of per-object maxima | 50 objects with pairs /57 | 1.302160 | 0.066441 |
| All recorded domains: largest pair | 17,614 pairs, austria_matched32 | 4.294004 | 0.214700 |
| All recorded domains: median across all pairs | 17,614 | 0.623058 | 0.031231 |

The five core objects without pairs and seven all-domain objects without pairs have unavailable maxima, not zero distance. Pair counts are not independent material counts, and these describe retained numerical candidates, not an exhaustive solution set. Source `src/analyze_feasible_branches.py:33–44` compares saved raw/corrected/primary/secondary views across runs and deduplicates exactly identical materials; it does not enumerate every trajectory candidate. Austria/Fresnel diagnostic discrepancy balls and cross-grid CN mismatch must remain qualified. Different recorded prior values do **not** demonstrate equal-prior minimizers. These numbers support only finite-discrepancy candidate separation, not exact-data or continuous-PDE nonuniqueness.

The existing branch receipt reports 2.881253 CPU seconds /4.435613 wall seconds, with zero new Maxwell LU, adjoint RHS or optimization for that original analysis. Those are existing analysis costs, not this review’s timing.

## Corrector material help/harm from existing core48 metadata

Source: `analysis/confirmation_core48/CORRECTOR_DECOMPOSITION.csv`. Counts below use the exact sign of its existing corrected-minus-raw NRMSE field: negative=lower recorded error, positive=higher recorded error, equality=exact recorded equality. They are descriptive material-error changes, not a new reconstruction success threshold. No new truth scoring was performed.

| Method | Lower NRMSE | Higher NRMSE | Exact equal | Mean corrected−raw NRMSE |
|---|---:|---:|---:|---:|
| CIE | 33 | 15 | 0 | +0.003057862 |
| Dual-CIE | 35 | 13 | 0 | −0.003105868 |
| G19 | 39 | 9 | 0 | −0.041498904 |
| G19 full APCS | 36 | 12 | 0 | −0.032416346 |
| G19 full current | 31 | 17 | 0 | −0.006538900 |
| DR | 40 | 8 | 0 | −0.023947217 |
| Integer event | 41 | 7 | 0 | −0.036729493 |
| Matched random integer | 40 | 8 | 0 | −0.018643661 |

Each row retains all 48 materials. The raw-to-corrected means printed in the report/manuscript match these metadata. For CIE, error decreases on 33/48 yet increases in the mean; the sizes of changes matter. The five principal CIE/Dual/G19/DR/event rows have no recorded original-Phi increases from raw to corrected, while material harm occurs in each. This directly illustrates why physical objective decrease and material error must remain separate.

For completeness, the CSV also gives direct raw-to-corrected counts 19 lower/29 higher/0 equal and mean change +0.100450507. **Do not call that shared-corrector harm:** direct’s raw is the initial material and its corrected field is the direct physical optimization incumbent (`src/run_experiment.py:120–146`), not the proposal-to-common-corrector transition.

## Review boundary

The review confirms the insertion and source/metadata consistency stated above. It does not certify every independent literature claim, every publication link, new stationarity/local-minimum conclusions, or final scientific attribution. No new result or success criterion was introduced. Whole-review CPU/wall and human/model review cost were not measured and remain unknown; individual tool timings are not substituted for whole-review cost.

