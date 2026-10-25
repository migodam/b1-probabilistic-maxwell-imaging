# Evidence exporter: code ready, execution NOT_RUN

`tools/export_evidence.py --root ROOT` rebuilds deterministic snapshot exports only. It never imports accounting/physics, resets budgets, writes the authoritative ledger, or changes any scientific gate. Worker compiled the source **in memory only**; exporter itself was not executed. Run after completed reports and again after publication tail.

Root outputs:

- `PHYSICS_LEDGER.jsonl`: every authoritative JSON-ledger event, in original list order with unchanged IDs/fields, including failures, preparation and corrections.
- `COST_SUMMARY.json`: per-stage original caps, used/remaining/exceeded; physical counts exclusively from `kind=reservation` RHS/factorizations. Informational observed counts never charge. CPU separates measured `job_completion`, estimated preparation, accounting corrections (including180-second plot allowance) and other charged events. Wall is summed only from job completions once. RSS receipts retain units/high-water semantics and maxima per unit; no summed RSS. Running reserved jobs have no final CPU yet.
- `FAILURE_LEDGER.jsonl`: failed job completions and their wrapper-interval physical charges; preflight events; failed/unknown/nonconverged scientific-result statuses and source pointers; E0 qualification/refinement history and mesh-memory plan; original finite-bank duplicate-keyword compiler log; original/repaired A0 analysis provenance. Explicit SMC failures/unknowns are exported only when recorded. Source reports remain authoritative.
- `A0_ANALYSIS_PROVENANCE.json`: original shape/layout/missing-prediction row counts and separate repaired-analysis path/summary. Invalid original rows cannot support a zero-witness conclusion. If repaired analysis is absent, it is marked NOT_FOUND rather than inferred.
- `PUBLICATION_ARRAY_MANIFEST.csv`: every NPZ/PNG/CSV artifact path, bytes and type under the configured scientific workspace; no hashes. The manifest excludes itself.
- `EXPORT_STATUS.json`: source snapshot/event IDs, report statuses, skipped/read-failed inputs, exclusions and untouched-gate/budget flags.

Only the configured scientific tree is enumerated; no credential scans or outside-tree traversal. Runtime/credential directories (`.git,.aws,.codex,.agents,.venv,.deps,node_modules,__pycache__`) and symlinks are excluded. JSON evidence extraction is limited to results and root scientific documents.

AST parsing alone does **not** establish compilability: the earlier finite-bank repeated keyword passed AST parsing but failed before physical execution. That failed zero-PDE job is retained via receipt/events and compiler-log provenance. The parent fixed it and separately compiled sources; export does not erase this history.
