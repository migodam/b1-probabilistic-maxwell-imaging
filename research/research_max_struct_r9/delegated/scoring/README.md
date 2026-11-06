# Frozen four-case offline scorer

Implemented `experiments/score_development.py`; **not executed for sealing or scoring in this delivery, no GT loaded and no solver launched**.

Default manifests are `configs/baseline_development.json` (16 noise-ball controls), `configs/baseline_historical_fixed.json` (8 historical same-code fixed-rank controls), `configs/landscape_development.json` (12 fixed/continuation controls), and `configs/baseline_cost_controls.json` (7 local hardware controls): 43 new tasks on b20d,40ac,r8_test_0000,r8_test_0027. These are already exposed development objects, not a new blind denominator. The script reads the manifests rather than a result-directory glob. With custom manifests, pass the same complete repeated `--manifest` list to both modes; every task must remain within the four fixed IDs and every ID must be represented.

From the R9 repository root, after all requested runs have finished:

```sh
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" \
  research/research_max_struct_r9/experiments/score_development.py seal

"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" \
  research/research_max_struct_r9/experiments/score_development.py score
```

Do not resolve the virtualenv executable symlink to a base interpreter. The scorer defaults to actual sibling `r8_worktree` for inherited data/metrics/receipts; `--r8-repo` explicitly overrides it. Paths recorded in the frozen manifests keep their existing semantics.

`seal` checks every requested receipt is completed or explicitly failed, with matching task identity, then stamps manifests/configs, public inputs, final predictions, paid stage metadata and available inherited endpoints using bytes/mtime (no hashes). Missing/RUNNING outputs block the seal. Successful baseline output requires `SEALED_PUBLIC_PREDICTION`; landscape's existing terminal success is `COMPLETED_PUBLIC_PREDICTION`. Both have fixed final prediction arrays. Failed runs remain declared terminal rows with missing material metrics and false success; no scientific success gate is added. Preparation-only output does not count as a completed solve. Explicit `FAILED_R9_ADAPTER` preparation failures are retained when no inherited receipt exists.

`score` revalidates the complete task set and every stamp before its first private metadata/array read. No solver module is imported. Pure material/noise/prior definitions are loaded from the actual R7 scorer used by R8. Material NRMSE is RMSE divided by the declared contrast amplitude, not relative material L2; IoU is prediction-only256-bin Otsu. Material success remains NRMSE≤.15,IoU≥.7,legal[.02,20]. Raw noise discrepancy, physical residual, original Phi and numerical stop remain separate. Historical grid equivalence is not inferred when absent from mapping metadata.

Truth resolution is copied from actual R8 scoring contracts and occurs only after the seal:

- b20d: selected entry in R6 `sealed_truth/scoring/historical_mapping_PRIVATE.json`, resolving the existing R2 truth path; amplitude from R6 public cases.
- 40ac: selected entry in R6 `sealed_truth/confirmation/mapping_PRIVATE.json`; amplitude from R6 public cases.
- r8_test_0000 and r8_test_0027: exactly those two R8 `sealed_truth/warm_test/<id>.npz` arrays; amplitudes from the existing R8 score metadata. No other test material arrays are opened.

Inherited R8 uniform controls are discovered only in the actual R8 `warm_all_manifest.json` and `history_warm_manifest.json`, with real sealed receipts/arrays. Public G/E/S/Y/k/radius/coords/sigma are compared exactly to the R9 case, and uniform0.1/no-NN/no-current-anchor initialization is verified. Other versions are not renamed R8, and missing exact inherited arms are explicitly recorded `NOT_AVAILABLE`. The new eight historical fixed controls exist because old R6 dense/explicit paths are not the same implementation. The seven cost controls add local fixed-rank runs for r8_test_0027 and local direct-physics controls. Existing0027 R8 endpoints ran on Windows/XINAN: retained historical rows do not replace fresh same-Mac fixed controls for CPU comparisons. Existing0000 Mac endpoints are labeled reused historical controls. Reused configurations and paths remain visible; no equal-CPU claim follows from matching inputs and caps.

Outputs under fresh `RESULTS/development_scoring/`:

- `SCORES.json`/`.csv`: final material metrics, residual, original Phi and independent saved-array Phi check, total/inner CPU and wall, LU/RHS, auxiliary counters, corrector stop/iterations, selected stage, exact provenance.
- `BUDGET_SCORES.json`: all saved incumbent material snapshots graded offline. Baseline costs use the existing R9-preparation-inclusive budget copy; landscape adds the recorded exact `pre_solver_CPU_seconds` when present. Already completed runs without that field use `process_CPU_seconds−cost.CPU_seconds` as `UPPER_SHIFT_including_output_tail`, never as exact startup time. Total process CPU remains the endpoint comparison.
- `STAGES.json`: paid stage/corrector metadata preserved from receipts; no GT-based stage selection.
- `INHERITED_AVAILABILITY.json`: actual matched and missing R8 endpoints.
- `COST_LEDGER.json`: new inclusive task costs and separately retained inherited costs, archived failed import-attempt index and measured archived receipts. Missing baseline failure CPU remains unknown; the known archived CPU is a lower bound. Current/receiver Gram/Cholesky operations are not counted as full-wave LU.
- `RECEIPT.json`: exact four truth paths scope, scoring-only CPU/wall and zero new physical/LU/RHS calls.

The outer R9 baseline total already includes preparation and R8 execution; never add inclusive and nested costs twice. Sum of concurrent task walls is workload, not queue elapsed time. Archived first-attempt failures are engineering costs rather than additional material samples. A supplied NN initializer's creation/training costs would remain external; current frozen development manifests use uniform initialization.

## Completed development scoring repair

After root sealed all43 outputs, the first scoring attempt raised `KeyError: cases`: actual R6 public metadata uses `records`. The scorer now accepts the actual schema (and the older `cases` spelling). `FAILED_SCORING_ATTEMPT_01.json` records the error, unknown CPU/wall, private-metadata-only access before the failure and zero physical calls. No partial scoring directory existed. Authorized attempt02 completed59 rows and10727 paid budget snapshots; source predictions/selection and the original43 result directory are unchanged. Compact paired results are in `DEVELOPMENT_COMPARISON.md`; saved Phi agrees within2.776e-17.

## Prepared six-case confirmation utility (not run by this delivery)

`experiments/score_confirmation.py` and pure `scoring_common.py` support:

```sh
python -B research/research_max_struct_r9/experiments/score_confirmation.py seal
python -B research/research_max_struct_r9/experiments/score_confirmation.py score
python -B research/research_max_struct_r9/experiments/score_confirmation.py figures
```

Use the existing virtualenv executable shown above in place of `python`. Root invokes these only after the requested66 tasks are completed/failed. The public case list comes from `configs/confirmation_public_manifest.json`; task coverage is the exact48 baseline +6 profile-continuation +6 noise-single-CIE +6 weak-profile handoff manifests, with all eleven methods for each of six cases. All stamps/identities are checked before first private metadata/GT access. Only `sealed_truth/confirmation/score_metadata.json` and the six explicitly selected material arrays are used for labels; family/shape generator metadata is not needed. No confirmation labels are read by preparation or `--help`.

The frozen `configs/portfolio_rule.json` is included in the seal. Before any GT read, `seal` selects each of the two exact candidate pairs by minimum final original Phi, with task-ID tie breaking, into `audit/CONFIRMATION_PUBLIC_SEAL_PORTFOLIOS.json`. The only portfolios are noise-global-A+ +profile and noise-GlobalDual +profile. Every portfolio charges both complete runs, including rejected work and separate current Cholesky/RHS; no other arm or GT oracle enters selection. The66 actual task cost total counts each executed branch once, so portfolio totals must not be added again as new work.

Score outputs go to fresh `RESULTS/confirmation_scoring/` with66 actual rows plus12 frozen portfolio rows, paid budgets, stages and costs. `MATERIALS.npz`/index are postseal offline exports for plotting, not inverse inputs. Figure mode reads those exports and the completed score receipt, then writes `FIGURES/confirmation/all_six_materials.png` (GT +eleven actual arms for every case; common full material range within each row, including failures) and `paired_cost_matrix.png` (NRMSE,IoU,CPU,including both portfolios). No new private reads or physics occur in figure mode. Generated figures still require visual review after root runs them.

## Prepared separate three-task development supplement utility

`experiments/score_development_supplements.py seal` then `score` covers exactly `profile_restart_control.json` plus `curvature_followup.json`:40ac fixed-penalty restart, b20d ordinary continuation from the old endpoint, and b20d continuation from the selected curvature proposal. It requires all three public outputs before opening only the two already-exposed historical labels. Outputs are separate `RESULTS/development_supplement_scoring/`; original43 scoring is not rewritten.

The complete public Hessian/eigensolve/probe receipt and arrays are stamped. Its CPU/LU/RHS are charged once to the proposal path, added to that path's budget prefix; the anchor has no diagnostic prefix. Incremental comparison reuses the already-paid development endpoint, whose original43 cost remains reported upstream. The utility does not interpret a negative restricted Hessian eigenvalue as a global or continuous-space basin result. Both new utilities are code preparation only until root authorizes their execution.

Root subsequently authorized the three supplements: they were sealed and scored in the separate directory, with no original43 rewrite. `DEVELOPMENT_SUPPLEMENT_COMPARISON.md` records that the fixed-penalty restart and both curvature follow-ups fail material recovery. Confirmation tools remain unexecuted until root's separate66-task seal/score invocation.

## Separate weak-stage handoff control

`experiments/score_weak_handoff.py seal` then `score` covers only `configs/weak_handoff_control.json`. It requires the one completed/failed public output first, and stamps the exact reused first-stage physical array, original profile receipt and original paid budget. Only40ac's existing disclosed label is graded. Outputs go to `RESULTS/weak_handoff_scoring/`, leaving original43 and three-supplement scores immutable.

`cpu_seconds`,LU and RHS are **new continuation work only**, which is the actual incremental research expenditure. `standalone_*` fields additionally charge the reused first-stage cost taken at the exact `profile_0_0.01_endpoint` paid budget, including recorded process startup prefix. If startup was absent the explicitly labeled upper-shift convention applies. The source initializer is checked equal to the declared paid stage. Existing first-stage Cholesky/RHS are included in standalone cost, and must not be charged a second time to the actual research ledger.

This conditional root authorization was executed once the weak task completed: the separate score reports material success (NRMSE .058351,IoU .969231), actual new CPU16.352524s and standalone CPU19.633731s including the reused3.281207s first-stage prefix. No new scoring-time physical calls occurred. The confirmation66 score remains unexecuted by this worker pending root's full-seal invocation.

## Final confirmation amendment:66 actual tasks

Root added `configs/confirmation_weak_handoff.json` (six tasks) based only on the disclosed development weak-handoff result, before any confirmation GT. The confirmation scorer now requires all66=6x11 public task outputs, including this fourth inverse manifest. The two existing portfolios remain exactly unchanged, so scoring has66 actual rows and12 portfolio rows, not additional weak portfolios.

Each weak task declares its public initializer, reused first-stage physical state, original profile receipt and paid budget; all are stamped in the confirmation seal. At scoring, the exact `profile_0_0.01_endpoint` budget determines the necessary initialization cost. This endpoint precedes the original full-profile poststage `details=True`/normal-residual extra factorization, so standalone current counts do not charge that unneeded later diagnostic. Recorded pre-solver startup is included, or the upper-shift limitation is explicit if absent. Initial material is checked equal to the declared paid physical stage, independently of any truth.

Weak rows and their budgets retain actual new CPU/counts plus distinct `standalone_*` fields that include the reused initialization. The experiment cost ledger sums only actual66 work, with no repeated charge for profile's first stage. The six-object material figure now has72 panels (GT+11 arms); weak CPU* and the paired matrix display standalone method cost including initialization. Raw actual continuation cost remains explicit in scores. Original43 development scoring, three supplements, one weak development score and the shared helper source underlying those seals are unchanged by this amendment. No confirmation seal, score or GT read was executed during preparation.
