# Observable19 offline analyzer implementation (not executed)

Implemented only `../../src/analyze_results.py`; no runner, factory, config, old scorer or accepted history CSV changed. This is implementation delivery, not scoring/scientific acceptance. No prediction seal was generated, GT array opened, scorer executed, or physics run by this work.

## Parent execution

`--root` is the **supplement** directory, not the parent R14 directory. The default resolves correctly from the script. Use the existing scientific Python environment (no dependencies installed).

1. Once the terminal queue accounts for all144 IDs, create an unauthorized seal:

```sh
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" research/research_max_struct_r14/supplement_observable19/src/analyze_results.py --make-seal "$SUPPLEMENT/predictions/PREDICTION_SEAL.json"
```

`SUPPLEMENT` here denotes the absolute supplement directory. Optional repeatable `--queue /absolute/attempt/RECEIPT.json` chooses terminal attempt receipts explicitly, in chronological order. Without it, all terminal campaign `RECEIPT.json` paths are sorted; later same-task jobs become the selected job. The seal preserves every supplied receipt/snapshot and all attempt rows; paths/bytes for all existing result files and missing expected views are recorded, without hashes/mtime/NPZ reading. Running or incompletely accounted queues are rejected. Child failures/missing files are retained, not required to pass quality.

2. Parent reviews/copies that seal, setting **`authorized: true`** explicitly in the separate copy. `--make-seal` always sets it false and never authorizes itself. The extra compatibility `authorized_scoring` field is false; this analyzer's authorization gate is the documented `authorized` field.
3. Only after approval:

```sh
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" research/research_max_struct_r14/supplement_observable19/src/analyze_results.py --seal "$SUPPLEMENT/predictions/PREDICTION_SEAL_APPROVED.json" --out "$SUPPLEMENT/analysis/observable19_v1"
```

Outputs must be a new supplement `analysis` directory. Existing output is rejected; raw results are never edited. Authorization revalidates the exact manifest snapshot, all task identities/status/file existence/bytes, terminal queues, and five view inventories before any truth import/load. Same-byte replacement is not detectable by this deliberately no-hash contract. New result files also invalidate the inventory. Do not seal during a continuing download/export.

## Reused contracts and schema

- Exactly144 tasks =24 fixed case IDs × G/S/O0 × PGP/GP, rank19, material576, exposed-followup label, original120/72 budget. Figure ordering is G-PGP,S-PGP,O0-PGP,G-GP,S-GP,O0-GP.
- Shared helpers are imported read-only from R14 `src/score_dualmetric.py` and `score_extensions.py`. Private metrics helper is imported **only after authorization/public plan export** from R14 `runtime/scoring_history/src/score_results.py`. Original `load_truth/material_metrics` are reused; no new NRMSE/IoU or recovery selection threshold is implemented.
- Truth path is exclusively R14 `sealed_truth/independent24/<case_id>.npz`, `material_coarse`; support `>1`, amplitude `max_truth-.1`, original material box [.02,20], exactly as `score_extensions.py:350`. No generator family/regime/private labels are read.
- Original lower-case five views and `CANDIDATES.npz:m` are used. Shared planner fixes min legal/eligible completed original-Phi candidates before truth. Saved view material is compared exactly against its candidate-pool row. Metadata candidate/Phi/D/R/CPU/phase/legal/eligible are audited against HISTORY; inconsistent saved predictions stay archived with issue flags, never silently replaced or tuned.
- Raw retains its original source plus any shared public initial-incumbent fallback flag. Original C2 normally snapshots `best` (including initial) as raw. Corrected verifies actual corrector phase; when no corrector candidate exists, original C2 `phase_best.get('corrector',best)` fallback is explicitly audited. Secondary missing uses primary deployment fallback with a reason. All missing/failed objects remain in24 denominators and null metrics have reasons.
- Original paid-prefix15/30/60/90/120, native common-process-CPU prefixes, and common componentwise final physical LU/RHS candidate-pool selections are reported. Selection uses public Phi only. Native per-candidate process time is not an exact external OS per-candidate clock; external full-child measurements stay separate. All actual spent costs remain charged; these are saved-pool views, never stopped reruns. Missing physical counters are null, not zero. Auxiliary setup/SVD/Green operations are separate and not declared free.
- External Windows queue helper `external.child_process_CPU_seconds` has precedence. `COST_LEDGER` retains actual factory setup operations from DENOMINATOR_SETUP/FAILED or wrapper, separately from original nominal `auxiliary_receipt.denominator_setup_cost`; O0 receiver solve deltas are already in physical aggregate and never added again. All raw receipts/counters and null reasons remain available. Queue wall/controller costs are per receipt; parallel/overlapping wall is not summed. `ALL_QUEUE_ATTEMPTS` preserves failures; selected task CPU totals and historical costs must not be added again to these raw attempt rows. Multiple receipts can repeat the same attempt, explicitly flagged rather than automatically double-counted.
- Old R14 independent24 primary CIE/Dual/A+-Arc(pair01)/11/G19/DR rows are copied as144 **historical reference** rows, including absent rows. They are not pooled into the new strict six-arm normalizations.

## Outputs and limitations

CSV+plain JSON: ALL_VIEWS, SELECTOR_AUDIT, SUMMARY, COST_LEDGER, QUEUE_COSTS, ALL_QUEUE_ATTEMPTS, PAID_HISTORY, CORRECTOR_EFFECTS, paired continuous differences, descriptive paired bootstrap, gain/loss membership, common-budget counterparts, HISTORICAL_REFERENCES, TRUTH_SOURCES, ISSUES. Receipt captures analyzer CPU/wall (including failed analysis) separately from solver costs.

Six material pages cover all24 cases, four cases/page, GT+six arms, shared per-case scale, dedicated colorbar column. Additional figures show continuous metrics, saved paid-cost curves, primary recovery membership, paired metric differences, raw-corrected D/R/Phi/material-error changes, and measured external CPU. Static HTML links all tables/figures. Display scales and descriptive posthoc recovery are not used for selection. Parent owns final scientific interpretation.

Checks performed: standard-library AST/compile and `--help` only; no runtime scoring, synthetic numerical test, truth loading, seal creation, remote access, or physics. Runtime behavior still requires review/execution after the real terminal seal. Audit implementation CPU is not a new solver cost; uninstrumented editing/inspection is unknown.
