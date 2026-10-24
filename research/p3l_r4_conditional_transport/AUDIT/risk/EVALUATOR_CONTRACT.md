# Locked risk inference and independent statistics

New files only: src/p3r4/evaluation.py and tools/evaluate_r4_risk.py. Running first-competition risk.py/train_risk.py were not modified in this follow-up.

Methods specification (write before any new locked observations/truth access):

```json
{
  "fallback": "anchor",
  "parents": ["R4L1_00000", "R4L1_00001"],
  "methods": [
    {"id":"anchor", "family":"anchor", "seed":0, "loader":"anchor"},
    {"id":"r2_selected", "family":"r2_selected", "seed":33, "loader":"r2", "checkpoint":"/absolute/original/R2/BEST.pt"},
    {"id":"phase16_seed1", "family":"phase16", "seed":2026102501, "loader":"r4", "checkpoint":"/absolute/new/BEST.pt"}
  ],
  "comparisons": [
    {"id":"phase16_vs_r2", "candidate_family":"phase16", "reference_family":"r2_selected", "action":"clipped"}
  ]
}
```

Freeze command: `PYTHONPATH=src <python> tools/evaluate_r4_risk.py freeze --methods SPEC.json --out METHODS_FROZEN.json`. Freeze is exclusive-create and stores byte-size/mtime-ns checkpoint identities, exact parent roster, family/seed labels, comparisons and fixed bootstrap settings. No hashes. The parent orchestrator owns freeze-before-new-observation chronology.

Online command: `... infer --freeze METHODS_FROZEN.json --contexts CONTEXT_DIRECTORY_OR_PACKED_NPZ --out ONLINE_DIRECTORY`. A directory contains parent.npz. Only observation,q0,chi0,J,residual are materialized; packed parent is identity metadata. No GT key or offline file is read. Original R2 loader is exact raw historical CNN; passive rectangle clipping uses original fusion.passive_box and yields fixed selected R2 alpha1 action. Each new R4 predictor produces raw and separately clipped output. Every failed prediction is explicitly replaced by that parent's anchor, recorded in failure list and fallback boolean; no parent is silently dropped. A missing/invalid anchor cannot be rescued and leaves ONLINE_COMPLETE absent. All outputs plus online-only anchor/q0 are saved before ONLINE_COMPLETE. Frozen checkpoint changes fail closed before prediction.

Offline command (launch new process): `... score --freeze METHODS_FROZEN.json --online ONLINE_DIRECTORY --truth inputs/offline_only --out SCORE_DIRECTORY`. Gate checks exact freeze snapshot, completed inference and different process PID before opening any truth. Truth files require native q66 and chi576. Per-method/per-parent arrays include full576 complex MSE/MAE, ROI MSE/MAE, exact host/transition/pure strata from anchor standardized field and GT field. Empty strata are NaN with count0, retained in NPZ; JSON means show nonempty-parent count. Full576-normalized stratum MSE contributions permit additive attribution without reweighting empty parents.

Raw and clipped actions remain separate throughout. Per-seed results are saved individually. Family comparisons average each parent's losses across frozen seeds before4000 paired parent bootstrap resamples(seed20261025099),98.33333333333333% CI. They do not ensemble predictions and do not treat seeds or pixels as independent parents. Comparison output includes absolute/relative improvements and seed roster. Claims remain with parent orchestrator.

Tiny historical-validation fixture PASS: evaluator_fixture/{METHODS_FROZEN.json,online/ONLINE_COMPLETE.json,score/RESULTS.json}. Four existing old validation parents, no training/new physics. Deliberately omitted one parent's J: phase output fell back to anchor and all4 remained in bootstrap; error retained. This fixture is API evidence only, never a locked scientific risk result.
