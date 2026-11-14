# Diagnostic instrumentation handoff

Implementation: `supplement_aplus48/src/diagnostic_instrumentation.py`. No frozen source edits. No inversion, physics, optimization, GT access, new scoring, or case selection was executed.

## Build interface

```python
build_diagnostic_runtime(
    frozen_C2_root,       # exact extracted C2 directory OR frozen runtime ZIP
    new_diag_root,        # separate directory that must not already exist
    legacy_repo=None,     # pass exact deployed inherited repository explicitly
    validate_only=False,
)
```

CLI, with parent-supplied platform paths:

```text
python diagnostic_instrumentation.py --frozen-c2-root <original-C2-root-or-zip> --new-diag-root <new-diagnostic-root> --legacy-repo <exact-inherited-repo>
```

Add `--validate-only` to transform/compile in memory and write no runtime. Builder copies **src/runtime only**, never RESULTS/queues/GT/materials/public arrays/configs. It writes `DIAGNOSTIC_BUILD.json` in the new root and replaces only the **private copies** of:

- runner: passive configure, DR pre/post/gap observers and final export;
- full_current: its inherited backend import points to private diagnostic_cie_fft;
- diagnostic_cie_fft: exact inherited module with setup/native-history/update observers;
- diagnostic_instrumentation: this observer implementation copied as a self-contained helper.

Parent must supply a new diagnostic manifest and unchanged public inputs. Absolute `input_file` paths to already-public files are simplest; builder does not copy input data. New manifest task campaign/root must remain separate from primary RESULTS. Standard runner CLI remains:

```text
python <new-diagnostic-root>/src/run_experiment.py --root <new-diagnostic-root> --legacy-repo <exact-inherited-repo> --manifest <diagnostic-manifest> --task-id <declared-task-id> [--cpu <parent-assigned-core>]
```

The original whole-process checks,120CPU cap, audit schedule, parameters, original selector and return/exception handling remain unchanged. Observer cost changes time buckets/stopping and is fully paid; the diagnostic output is never eligible to replace primary supplement predictions. The private observer has no physics or optimization interface.

## What is retained

Every completed CIE-family update saves its existing native `history[-1]` row immediately after append, **before** original iteration_callback can raise BudgetEnd. Native stage-end history is saved too. The optional data is copied only; no objective/gradient/norm/projection is recomputed.

Full snapshots retain configured global iterations **1,8,32,128,512,1024,1500,2000,2500,3000**, plus first/last observed completed update in each native stage. An optional task field `diagnostic_target_iterations` adds up to two positive integer checkpoint indices from already-public primary raw/proposal audit references (CIE global_iteration or DR iteration), never GT. Example: `"diagnostic_target_iterations": [417, 963]`. Duplicate indices are deduplicated; an index already fixed shares the same snapshot. A replacing “last” slot ensures the latest complete observation survives an ensuing BudgetEnd. A hard20-state cap fails explicitly rather than silently evicting required checkpoints. Four stages plus10 fixed and at most2 extra checkpoints permit at most20 unique CIE snapshots; first/last/fixed/target aliases share copied arrays.

The copied runner passes the task field directly into configure/Recorder. DIAGNOSTIC_RECEIPT records `requested_target_iterations`, `reached_target_iterations`, and `unreached_target_iterations`; state-index slots include `cie:target:<iteration>` or `dr:target:<iteration>`. Only naturally reached completed updates count. No extra updates, late audit, budget extension or primary re-selection is performed. Original initial-state index0 is not a completed native update and should use the already-saved physical initial snapshot rather than a target entry. The public iteration references permit comparison at the same recorded update despite observer overhead changing CPU audit/stop times; they do not guarantee cross-run bitwise equality.

CIE snapshot fields: m_before_readout,m,J,Q,A,grad,direction,pr,alpha; scalar stage/local/global indices, betas and native counters. At that exact phase old_grad=grad and old_direction=direction. Static actual JD/Vd, metrics, bounds, onehop blocks, frequency order and stage betas are retained once/per-stage without new actions. No inner MM root/color state is recorded by this minimal implementation.

Specifically, static arrays include JD,Vd,Od,Os,Ls and bounds; each stage metadata includes rho,sourcewise,receiver_rank,mf,frequencies,beta_current,beta_material. Inherited source decoding is explicitly UTF-8. Parent subsequently reported normalized UTF-8 text agreement for eight actual remote legacy sources; this worker did not perform a new remote comparison.

DR full snapshots retain actual pre X/Q, constitutive cX/cQ/cm, reflected A output and post X/Q at the same fixed steps plus first/latest. DR rows are recorded every completed step, but gap is null unless the frozen loop itself computed an audit-time gap; the existing value is then attached without an extra contraction. `last_A_info` is copied after the corresponding native projection, using its actual `multiplier` key. No A-shadow or post-loop C-shadow is newly computed. Pre-state material cm(k−1) and post iterator z(k) are explicitly indexed.

## Export and costs

Each new task output receives:

- `DIAGNOSTIC_NATIVE_HISTORY.jsonl`
- `DIAGNOSTIC_DR_STEPS.jsonl`
- `DIAGNOSTIC_STATES.npz` (uncompressed; retained arrays only)
- `DIAGNOSTIC_STATE_INDEX.json` (snapshot IDs, aliased slots, array prefixes, fields and phases)
- `DIAGNOSTIC_RECEIPT.json`

The modified private runner flushes before its final receipt and adds diagnostic-only metadata. An atexit fallback flushes on other ordinary exits. Flush is idempotent; failures are retained/reported and do not invent missing updates. Fatal kill or interpreter failure can prevent export; no crash-proof persistence is claimed.

The diagnostic receipt records retained/peak state counts, native/DR rows, nested observer-copy CPU, export CPU/wall, output bytes, original termination context, and failure status. These costs are **nested in** whole-process CPU, not additive to it. External whole-child CPU remains authoritative; receipt-writing/final interpreter-exit costs fall after the nested export timer. No extra CPU is granted to compensate instrumentation.

## Static validation performed

`STATIC_WIRING_VALIDATION.json` records:

- exactly one known CIE completed-update, stage-end and stage-entry anchor;
- exactly one known frozen main DR loop, pre-state, post-state, native-gap, configure and export anchor;
- after stripping only marked diagnostic statements, the original R8 loop and frozen runner ASTs are exactly restored;
- full_current changes only its one exact inherited-module import;
- transformed copies and all copied Python source compile as text.

Validation used the actual `runtime_confirmation_C2_final.zip` and existing inherited source text. It did not import NumPy or any Maxwell solver, open GT/material payloads, launch any task or write a diagnostic runtime. No SHA checking was performed. Unknown whole-review timing remains unknown.

Parent owns case choice, launch authorization, scoring order and mathematical analysis.

