> 执行状态更新（2026-11-14）：112预测已终态、公开验收0问题、封存授权后完成评分，结果位于analysis/stack_comparison。代码只修正了派生cohort标签，未改预测/指标/选择；完整解释见audit/STACK_COHORT_LABEL_CORRECTION.json。

# Fixed A+ version and frozen-NN stack comparison

`compare_stack.py` is implementation only. Do not run until the parent reviews all
four terminal prediction seals and creates a separate authorization with:

```json
{
  "authorized_stack_scoring": true,
  "scope": "FIXED_EXISTING_VARIANTS_AND_FROZEN_NN_FOLLOWUP",
  "root": "/absolute/local/research_max_struct_r13",
  "supplement_root": "/absolute/local/research_max_struct_r13/supplement_aplus48",
  "stack_root": "/absolute/local/research_max_struct_r13/supplement_aplus48/stack_run",
  "seals": {
    "confirmation_core48": "/absolute/local/core48_seal.json",
    "confirmation_mechanism16": "/absolute/local/mechanism16_seal.json",
    "aplus48_extension32": "/absolute/local/extension32_seal.json",
    "aplus_stack112": "/absolute/local/stack112_seal.json"
  }
}
```

All seals independently need `authorized_scoring: true` and valid original
size/mtime/directory inventories. Parent-provided path aliases must restore the
seal paths; this code never rewrites seals, timestamps, source results or configs.
Failed/missing tasks remain in expected denominators. Terminal queues with errors
are permitted by the original seal validator; partial eligible paid candidates
remain visible, with failed status and issues retained.

After authorization, use the existing Python runtime with NumPy/Matplotlib:

```sh
python /absolute/supplement_aplus48/src/compare_stack.py \
  --root /absolute/research_max_struct_r13 \
  --supplement-root /absolute/research_max_struct_r13/supplement_aplus48 \
  --stack-root /absolute/research_max_struct_r13/supplement_aplus48/stack_run \
  --authorization /absolute/parent_stack_authorization.json \
  --out /absolute/research_max_struct_r13/supplement_aplus48/analysis/stack_comparison
```

`--include-direct-reference` adds original NN→direct16 as a separate reference,
excluded from the six-solver factorial/common-six limits. `--no-figures` skips
rendering. `--baseline-analysis` changes the exact-index metric cache directory,
otherwise existing `analysis/comparison48` is used.

The public data have no family labels. For a family-grouped gallery, provide
`--family-labels /absolute/authorized_case_family_map.json`, a JSON mapping full
case IDs to family names, and name exactly that file in the authorization as
`authorized_family_labels`. Without it, every case is still shown, in full-ID
order, with family unavailable explicitly stated; no families are guessed from
IDs. Labels are descriptive and cannot affect selection.

## Frozen comparison and outputs

- Uniform same48: CIE, Dual-CIE, G19, DR, R13 A+ quintic-MM/no guard, and the exact
  R7 final **readout options** inside the R13 protocol/current wrapper. This is
  not a rerun of the entire historical R7 experimental contract.
- Same predeclared16: six solvers × uniform/NN × initial/raw/corrected/primary/
  secondary. NN G19/DR reuse original mechanism16; NN CIE/Dual/MM/arc use stack112.
- New uniform-six common cost table uses per-object
  `min(120,actual whole-child CPU of all six)`. Any missing cost makes the common
  limit null. NN-six and matched all-twelve limits are separately labeled. This
  never replaces the prior five-method table.
- Paid prefixes15/30/60/90/120 choose the original minimum Phi over already
  audited legal eligible candidates, without GT. No checkpoint interpolation.
- Public indices, limits and run index are serialized **before** opening scored
  cache metrics or truth/material arrays. Score helpers are the original
  `material_metrics`, `load_truth`, `truth_spec`, `physical_counts`.
- Existing metrics are reused only for the exact same case/method/candidate index;
  missing/mismatched caches use the existing scorer on requested indexed states.
  No full candidate-pool GT ranking or new selector is computed.
- LF CSV + JSON contain material error, support IoU, D/R/Phi, noise qualification,
  actual whole-child CPU, candidate cumulative CPU/LU/RHS, failure status, source
  paths, version counters and all expected rows. Unknown values remain null.
- Within-object NN−uniform, corrected−raw, and arc−MM paired deltas and4000-replicate
  object bootstrap intervals are descriptive and unadjusted. Corrector cumulative
  CPU/LU/RHS deltas refer to selected audited candidates; whole-job costs repeat
  across views and must not be interpreted as phase costs or summed repeatedly.
- Compact all48 GT/MM/arc gallery and six-panel same16 NN paid-prefix curves /
  paired scatter plots supplement the existing five-method galleries.

No forward, adjoint solve, optimization, training or inference is performed.
Existing R8 physics-network training plus train/validation generation cost is
402.558889 Mac CPU seconds (`R8/NEURAL_WARM_STARTS.md:55`), not zero. Other historical
training costs and already paid R13 inference/preparation remain separate. New
112 run costs and reused old run references are separated; cumulative history,
internal timers and projector counters are never added to whole-child CPU.

These controls use existing exposed objects. They are not fresh blind confirmation,
parameter tuning, proof of all-A+ superiority, or final scientific interpretation.
