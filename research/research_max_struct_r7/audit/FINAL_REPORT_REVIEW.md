# Independent final report evidence review — primary and targeted scopes sealed

Review date: 2026-11-04. Scope: current `RESEARCH_REPORT_ZH.md`, `PAPER_DRAFT.md`, primary180 analysis/scoring outputs, extra47 and rational12 public score records, prior implementation/scoring audits, cost and release scope. This is a read-only evidence/claim review. No algorithm, prediction, physical call, hash or score was changed. No sealed truth or GT arrays were opened. The original review excluded targeted20 while running; the final addendum below now checks only its existing post-seal public JSON/MD scores and receipts after all80 predictions sealed. Final scientific novelty judgment remains with the parent Codex.

## Outcome and matters raised immediately

The primary quantitative story is supported: A+ retains A's observed 21-object success set and reduces total inversion CPU31.4%, while remaining11.2% more costly than sourcewise Dual. This is a useful cost result relative to original A, not evidence of being cheaper than the established controls or universally better in recovery. The report generally preserves that distinction.

Two factual/scope issues were raised to the parent as soon as found:

1. **Stale cost-ledger status — corrected during review.** `COST_LEDGER.md` originally said the primary180 and initialization18 were still running. Its current header now accurately says primary sealed/scored and targeted confirmation running. It separately retains packaging and unmetered-work limits.
2. **Rank32 material-score denominator — resolved by separate scored addendum.** [History: the initial review found this pending.] The original rational scoring artifacts contain12 rows: three anchors × four endpoints (anchor continuation, rank2, rank8, local exact origin). They contain **no rank32/coherent row**. Main-report wording initially grouped rank2/8/32 then attached “12 endpoints, zero material successes,” which could imply rank32 was included. The parent acknowledged this and requested separate post-seal scoring of the three rank32 endpoints, preserving the original12. The separate `RESULTS/scoring/rational_coherent/RECEIPT.json` and `SUMMARY.md` now verify0/3 material passes with0new physics; report/draft explicitly distinguish original12 plus coherent3, all on three old anchors. Do not turn all15 corrected endpoints into15 novel finite candidates.

A cross-benchmark qualification was also requested and is now explicit in both report and draft: vector3D has its own quadratic+l2 prior, normalization and relative-material-L2 NRMSE. “Same prior/Phi” applies within a benchmark, not across scalar2D and vector3D.

## 1. Primary counts, costs and success sets verified directly

I parsed the public `RESULTS/scoring/primary/SCORES.csv` (180 rows,30 objects, six arms) and checked the arm aggregate JSON and prose. All180 rows are retained; no failed/poor object is omitted from the quoted denominator.

| Arm | Material successes /30 | CPU s | Physical LU | Complex RHS |
|---|---:|---:|---:|---:|
| CIE single |21|1543.595881|54873|1154196|
| Dual-CIE |20|1539.147756|53408|1137060|
| Global-Dual |20|1512.969297|53923|1137684|
| Original A |21|2495.142684|54384|1148820|
| A+ |21|1711.335523|52156|1103040|
| Direct physical |10|838.182328|47754|1009764|
| Total |180 pipelines|9640.373469|316498|6690564|

The report's rounded table and primary totals match. A+ total CPU reduction versus A is31.4133202%; A+ overhead versus Dual is11.1872149%. These are **ratios of total observed CPU**, not the mean of per-object percentage speedups, wall-clock concurrency speedups or free material-kernel timings.

A and A+ have identical actual pass/fail object sets. Single-CIE has the same count21 but a **different set**. Thus “same21 as A” is correct; “same success set as CIE single” would be false. No such false equivalence is needed for the main result.

Same-grid material+CN99 counts in the report also match: single15/24; Dual/global/A/A+14/24; direct9/24. Cross-grid six use a noise-only threshold with unbudgeted model mismatch and cannot be pooled into a calibrated error certificate.

## 2. Case IDs and strong-baseline comparisons

Verified from exact public score rows:

| Case | single-CIE NRMSE/pass | Dual NRMSE/pass | Global NRMSE/pass | A NRMSE/pass | A+ NRMSE/pass |
|---|---|---|---|---|---|
|obj12 curved bar|.721134 / fail|.735001 / fail|.735114 / fail|.039646 / pass|.039145 / pass|
|obj20 slit shell|.070223 / pass|.296361 / fail|.296552 / fail|.069809 / pass|.068721 / pass|
|obj09 thin winding|.097286 / pass|.098332 / pass|.095241 / pass|.098047 / fail|.097571 / fail|
|obj17 coupled capsules|.064319 / pass|.064815 / fail|.065159 / fail|.065727 / fail|.065689 / fail|
|obj24 cross-grid curved bar|.392973 / fail|.037939 / pass|.038113 / pass|.037577 / pass|.037257 / pass|

Obj20 is an improvement over Dual/global, **not an A-only recovery against every established baseline**: sourcewise single-CIE already succeeds. Obj12 provides the stronger all-three-CIE-controls failure/ A-path recovery example. Obj09 has IoU .696629 for A/A+ versus .712644 Global (Dual .704545), so it is an IoU threshold reversal rather than catastrophic continuous-error loss. Obj17 is the separate single-CIE comparison: single IoU exactly .7 passes, A/A+ .675 fails. Obj24 illustrates that single and Dual can reverse roles. The revised main report correctly describes these contrasts and does not conflate obj09 with obj17.

Consequently, do not promote a Dual-only comparison into superiority over the best available CIE variant per object. Such an oracle baseline is also not a runnable method unless its selection rule is defined independently of GT. The complete six-arm table and fixed algorithms are the appropriate comparisons.

## 3. Confidence intervals and their interpretation

The quoted report values match public paired-bootstrap ledgers. All are left minus right,10000 paired material-object resamples, with the stated fixed heterogeneous acquisition mixture and no multiplicity adjustment.

| Comparison | Metric | Mean difference |95% interval |
|---|---|---:|---|
|A−Global|NRMSE|−.04381569|[−.10087691,−.00567822]|
|A−Global|material success|+.0333333|[−.0666667,+.1333333]|
|A+−A|NRMSE|+.00060392|[−.00092329,+.00287233]|
|A+−A|CPU seconds/object|−26.12691|[−33.81986,−18.62342]|
|A+−Dual, secondary|NRMSE|−.03777867|[−.09291871,−.00308697]|
|A+−Dual, secondary|CPU seconds/object|+5.73959|[+1.91586,+9.54081]|

A−Global's median paired NRMSE interval crosses zero; its success-rate interval crosses zero. The report correctly avoids converting a negative mean-error interval, driven partly by a few large recoveries, into stable per-case or rate superiority. A+−A's observed success-difference bootstrap [0,0] is degenerate because the observed binary labels are identical; it is not future equality or formal noninferiority. No margin was prespecified. These qualifications are explicit in the report and `INFERENCE_LIMITS.md`.

Sourcewise Dual comparisons are labeled secondary rather than silently replacing the cleaner A-vs-Global material-block contrast. The public secondary file contains both frontier-summary bootstrap intervals and a separate bootstrap at each common maximum; their slightly different percentile endpoints are separate Monte Carlo calculations. The main report does not conflate those numeric intervals. Preserve an exact source when quoting either interval in a final table.

## 4. Actual-cost matching

At each pair/object, the common observed window is max(first saved cost) to min(last saved cost), on CPU, LU and RHS separately. Fixed window fractions select the last saved original-Phi incumbent at or below budget; GT grades afterward. No stopped trajectory is extrapolated, no truth-based iterate selection or interpolation is needed.

Verified common-maximum results:

- A versus Global at CPU:17/30 versus20/30; mean NRMSE difference interval [−.10855328,+.00507330] crosses zero.
- A+ versus A at CPU:21/30 versus17/30; its mean NRMSE interval also crosses zero. Report the observed count without a new rate-significance assertion.
- A+ versus Dual, secondary:22/30 versus20/30 at common maximum CPU/LU/RHS. This is checkpoint success, **not final22/30**. Endpoints remain21/30 versus20/30.
- Original A versus Dual at common maximum CPU is15/30 versus20/30. This negative comparison is in the public frontier ledger; omission from a brief prose summary is acceptable only while the full comparison remains linked and no broad equal-CPU A superiority is claimed.

CPU includes auxiliary work. LU/RHS do not include that work as equivalent full-wave factorizations, so equal LU/RHS must not be called equal compute. Every-ten-physical-record and sparse auxiliary-stage snapshots leave unused-budget gaps; float32 offline material snapshots have their own precision limit. The independent audit reports maximum float32-vs-final NRMSE difference≈1.99e−9 and final saved-budget-to-process-CPU tail≈.042407s. Full task CPU includes the serialization tail, while the last saved budget point does not; this is disclosed and small but not zero.

Rational serial observer offsets still include preceding comparator corrections, and their narrow setup timer omits anchor/per-action preparation. The revised report explicitly refuses a rational matched-CPU global-search claim, consistent with the implementation review.

## 5. Extra runs and rational denominators

Public `extra/SCORE_RECEIPT.json` and `EXTRA_ANALYSIS.md` support47 scored runs: development20 (four disclosed objects × five readouts), old-digit9 A+, initialization18 (six existing objects × three arms). They add **zero new independent objects** to the30 denominator.

The quoted development CPU284.21→197.78 and3/4 recovery for original A/final A+ match the extra summary. Projected gives2/4. Digit A+ is8/9, and original R6 Dual9/9 is clearly historical/different wiring rather than a rerun control. Repeated m=1 initialization yields A/A+5/6, Dual4/6; no better-start GT selection replaces primary endpoints.

Original rational scores contain12 **corrected endpoint records**, not12 independent materials and not12 newly proposed rational actions: four records for each of three disclosed historical anchors. All12 material-success flags are false. This supports the stated negative result for rank2/rank8/local/continuation comparisons. [History: rank32 had physical-only receipts at the initial snapshot.] Its separate completed post-seal public score addendum now reports0/3 material passes; combined0/15 is valid only as corrected endpoint count on three disclosed anchors, including controls.

## 6. Measured and vector evidence

The Fresnel statements match the independent acquisition/inversion audit: two physical targets, three methods each, target-absent calibration only,234 training and1530 actual withheld channels per target; train residuals.14193/.10062 and withheld.26755/.23444; identical saved endpoints by target; contrast saturation; twin corrector cap. No exact per-pixel truth or fully calibrated interior horn field is available. The report correctly qualifies equal host-objective observations as consistent with rejecting auxiliary candidates, without inventing retained aux arrays or certifying the exact argmin.

Vector totals16 initial +16 common long continuations,125/216 real isotropic material coordinates and375/648 vector systems match the vector contract. The125-grid A NRMSE values quoted in report/draft (.17844,.19965,.13300) match the long-correction table. Joint cap exits remain visible. “Methods approach comparable regularized solutions” is supported in this finite development setting; exact mathematical equality/full convergence is not proved.

The current report/draft now explicitly separate the vector quadratic+l2 prior, data normalization and relative-material-L2 NRMSE from scalar2D's smooth/TV prior and amplitude-normalized NRMSE. Do not pool their thresholds or statistics. Binary Mie6.58%→27.59% failure and fixed-fraction8.40%→8.20% correction remain negative/limited forward-validation evidence. Neither small algebraic residuals nor actual3D arrays establish continuum imaging superiority. These caveats are correctly retained.

One minor manuscript clarity recommendation: typeset general-q conjugate-transpose terms unambiguously (`Re tr(ZᴴT)+½ tr(TᴴATB)`). Plain Markdown `Zq*T`/`T*AqTBq` can read as multiplication or trigger emphasis even though the implementation algebra is correct. This was notation clarity, not a newly found sign defect. Final addendum check confirms both `OPᴴ` and `ZqᴴT` are now corrected by the parent in draft§5.

## 7. Costs, release rights and outstanding targeted addendum

The updated ledger correctly distinguishes total task CPU from parallel elapsed wall, overlapping nested timers/counters, LU dimensions and auxiliary small solves. Retrieval, failed access attempts, reporting/plotting/packaging and two Fresnel host-audit evaluations have incomplete telemetry. Unknown money/energy/RSS is not zero. A total of9640.37s is the **primary inversion** total, not an all-work project bill. Packaging/publication status must remain pending until actual receipts and direct remote readback exist.

The current `PUBLICATION_SCOPE.json` and `package_evidence.py` now implement the required raw Fresnel `.exp`, measured-Y endpoint NPZ, downloaded fulltext/PDF/HTML/TXT, copied archives, caches and symlink exclusions; own generated figure PDFs remain eligible. The small Fresnel public-derived bundle excludes measured_Y. I did not execute packaging or assert an actual release exclusion manifest already exists.

**Release timing caveat:** the packager's accepted parent seal globally permits synthetic truth under R7. Existing primary180 completion alone must not authorize release of the later targeted20 truth while its predictions are running. Parent release authorization should cover completion/sealing of every cohort included in the release, including all four targeted baselines; no actual leak was observed at that earlier pending-packaging snapshot. All targeted predictions are now sealed/scored; parent reports final packaging in progress, outside this numerical evidence check.

[History: the original review correctly treated `[TARGETED_RESULTS]` as pending and did not count running jobs as complete.] The placeholder is now replaced by a completed, separately denominated targeted-family section in the report and draft. The final addendum below checks the20×four public scores without GT-array access. The family choice remains post-hoc; it is not pooled into primary30 or rebranded as a universal blind test.

## Final primary-scope disposition

Primary counts, costs, confidence intervals, case IDs and the bounded A+ cost conclusion are **supported by the public sealed evidence**. Strong single-CIE versus Dual/global differences, no rate-superiority claim, no noninferiority proof, explicit physical/continuum limits and copyright exclusions are preserved. The cost-ledger and cross-benchmark wording issues were corrected during review. Rank32 score qualification and the targeted-family evidence checks are now resolved by the final addendum below; the parent retains final interpretation, release and novelty judgment.


## Final addendum — targeted80 and coherent3 complete

This addendum reads only existing public JSON/MD: `RESULTS/targeted_analysis/{SCORE_RECEIPT,ARM_SUMMARY,FAMILY_SUMMARY,SCORES,PAIRWISE_COST_FRONTIER_SUMMARY,APLUS_FAILURES}.json`, `ANALYSIS.md`, `COST_AND_FAILURES.md`, coherent3 score receipt/summary, and the updated main report/draft/chat handoff. No new experiment, physics, GT-array access, algorithm change, bootstrap recalculation or sample replacement occurred.

**Targeted denominator and exact counts:** scorer receipt is COMPLETE, POST_COMPLETE80_TASK_PUBLIC_SEAL_ONLY,20 objects/80 tasks,0new physical calls; parsed public scores contain exactly80 rows/20 case IDs. The order below is A+/Dual/Global/single:

|Scope|A+|Dual|Global|single|
|---|---:|---:|---:|---:|
|Material pass /20|12|9|9|8|
|Material+CN99 /20|9|8|8|6|
|Curved bar pass /10|9|8|8|5|
|Slit shell pass /10|3|1|1|3|
|Total CPU s|1089.689649|1061.574844|1027.612678|855.186925|
|Physical LU|43594|46876|45126|43013|
|Complex RHS|849516|909312|870528|836184|

These values match the report, new draft§6.2 and chat handoff. A+ and single shell successes are the exact same targets11,17,18; Dual/global succeed only on17. Curved-bar target03 is A+-only among the four fixed controls. Therefore shell gains over Dual/global are not A+-unique; the selected-family result remains mixed and cannot be generalized to universal superiority.

At equal CPU versus single, the public frontier summary at common maximum has10/20 versus8/20 and mean NRMSE difference−.07217114,95% CI[−.16501200,+.00548403], crossing zero. The report correctly distinguishes this from full-allowance endpoint performance. A+−Dual endpoint mean difference−.068311 has95% CI[−.142132,−.0124211]; success difference+.15 has CI[0,.30], touching zero. Unadjusted endpoint intervals are not formal noninferiority or universal-rate proof. The20 are a post-selected pair of family/acquisition distributions and remain separate from primary30.

All eight A+ failures are present: targets02,10,12,13,14,15,16,19. They comprise one curved bar/seven shells; all fail CN99. Public stop strings give exactly two numerical-convergence exits (target10 projected gradient, target12 relative function decrease) and six iteration caps. `APLUS_FAILURES.json` records a strictly better feasible offline GT-Phi witness for all eight; this is a reported post-seal diagnostic, not GT-based candidate selection, uniqueness or identifiability evidence.

**Coherent3:** its separate COMPLETE receipt reports three full576-variable endpoints,0/3 material passes, no new physics, and preserved original12 seal/results. Thus0/15 corrected endpoints on three old anchors is verified;12 includes continued-anchor/local controls and is not12 novel proposals. The revised report/draft/handoff now state this correctly.

**Resolved history:** the earlier targeted/coherent pending statements above are historical and have been explicitly superseded. All reviewed scientific-result counts and bounded claims are now supported. Packaging/publication still requires the parent's actual release authorization/receipt; scientific sealing alone is not publication.

**Final disposition:** targeted80 and coherent3 tables and updated report/draft/handoff are factually consistent within this requested scope. The two Hermitian notation edits are confirmed resolved. No unresolved numerical-count, family-denominator, stop-count or confidence-interval mismatch remains. Original predictions, algorithms and numeric sources were not modified. This review is complete.
