# Public-only prediction acceptance script

New standalone source: `src/audit_predictions.py`. It does not import a Maxwell model, read truth/labels/observed-data quality scores, factorize a Maxwell matrix or run an optimizer. It loads only selected public operator/initializer/bounds arrays and existing prediction/metadata artifacts.

```text
python -m research.research_max_struct_r13.src.audit_predictions \
  --root research/research_max_struct_r13 \
  --manifest research/research_max_struct_r13/configs/CAMPAIGN.json \
  --output research/research_max_struct_r13/audit/NEW_DIRECTORY
```

The queue must have an explicit terminal canonical receipt and no running/pending work. Output must be a new directory outside that campaign's immutable predictions/queue. No fresh campaign execution has been performed by this worker; parent authorization is still needed before a future fresh audit.

The full manifest task denominator is retained, including missing artifacts and malformed task bundles. Exact index/selector/timing issues are separate from numerical defects. Unknown counters, unavailable state views and unproved DOF fields stay null with reasons; they are not converted to zero.

Checks include:

- Contiguous history/material indices; legal/eligible flags independently recomputed using the original runner's box tolerance and each saved absolute CPU/deadline/total cap.
- Primary minimum legal eligible Phi; secondary minimum `(R,Phi)` over eligible saved noise-feasible candidates; raw and corrected index consistency; prefix minimum-Phi choices with explicit no-future/no-ineligible checks.
- Snapshot materials against their exact CANDIDATES row and public initial m_init. Methods sharing exact G/E/receiver maps/components/m_init are grouped and their initial m/j/Q/Y_prediction max differences recorded. NN groups remain separate. Cold auxiliary CIE current is not compared as if it were the shared physical initial state.
- For initial/raw/corrected/primary/secondary saved physical states: `J-D_mQ`, `Q-E-GJ`, `Y_prediction-receive(J)` Frobenius norms with declared physical units and explicit relative denominators. No arbitrary numerical pass threshold or recovery interpretation is applied.
- Common S, source-specific `S_ports`/3D S arrays, original per-source receiver-index lists, and cell-major vector material repetition are supported directly from public matrices. No absent receiver is zero filled. Actual vector/measured validation still awaits parent execution; this worker's real-case test was scalar development_v3 only.
- Physical LU/RHS canonical cost coverage, unknown external CPU, native PR/stage metadata and ambient-versus-independent DOF metadata are preserved. A reached/setup stage is explicitly not a completed native schedule claim. Nested physical deltas are not summed twice.

Every validation G/S action and source column is counted separately from the solver's paid work. Validation CPU/wall includes the script-entry clock before NumPy import, reads, comparisons and matrix products, excluding final report serialization/exit. There is no new Maxwell factorization or optimizer work. Public arrays are reused/deduplicated only within each case, including separate NN initializers; the case cache is discarded between cases. NPZ uncompressed member-size preflight and configurable `--max-array-MiB` (512 default) bound each requested load and the cached per-case public arrays. The report labels this as an array/cache cap, not a measured process peak RSS.

Outputs are `AUDIT.json` and concise `README.md`. The final real-shape check used only the already exposed terminal `development_v3_local` campaign: all six expected task rows retained, no exact wiring/index/availability issue reported. Initial physical m/j/Q/Y_prediction matched exactly across those methods. Numerical residual values remain available for parent interpretation without a pass gate.

Selector verification uses the frozen history Phi/R/noise-feasible values; it does not independently reevaluate every historical objective or forward state. Matrix consistency checks apply to the saved physical snapshot views that actually exist. A correct minimum index is not proof that every supplied history objective is correct, nor a material-truth certificate.

Validation records: `prediction_audit_development_v3_final/` is the final version; `prediction_audit_development_v3/` is the earlier successful version before adding raw/corrected selector checks. Both runs' separate matrix-action/CPU costs are preserved in `PREDICTION_AUDIT_VALIDATION_RUNS.json`; no failed attempt was hidden or cost imputed. No active runner, source bundle or frozen result was altered.
