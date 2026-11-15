# Fixed paid-primary diagnostic — implemented, not run

Entry point: `src/diagnose_paid_endpoints.py`. Static check only: AST parse, bytecode compile and prohibited-call inspection passed; the module was not imported/executed. No public/prediction array, GT, physics or optimizer was opened/run. Receipt: `audit/PAID_ENDPOINT_DIAGNOSTIC_STATIC_CHECK.json`.

After the parent creates and approves the **new R15 main264 prediction seal**, use the actual deployed paths:

```powershell
python "$r15Root/src/diagnose_paid_endpoints.py" --root "$r15Root" --legacy-repo "$legacyRoot" --seal "$approvedR15Seal"
```

Default output is the new `audit/paid_endpoint_diagnostic_v1` directory; `--out` may choose a different new directory under R15/audit. Existing outputs cannot be overwritten. The seal is read and validated **before output creation and before importing numerical physics**. The schema is the actual `analyze_campaign.py` schema `r15_twofold_exposed24_prediction_seal_v1`, with `authorized=true`, exact root/campaign/manifest snapshot, all264 terminal jobs and unchanged sealed result inventories/byte sizes. Old R14/old-endpoint authorizations are rejected. `authorized_scoring` is retained as provenance but not required for this no-GT validation.

The fixed44 rows are 0009/0020/0021/0023 × g19/g19_s19/twofold/twofold_hop/twofold_dr/twofold_ap/twofold_full_dr/twofold_reset/cie/dual_cie/dr. Each uses **only the already-selected sealed primary.npz `m`**, with no fallback, selection or candidate generation. Missing primary files remain `MISSING_PRIMARY`; setup/material/physics failures remain `FAILED`; each stays in denominator44. Unknown original termination, CPU and noise fields stay null.

Actual inherited APIs inspected:

- `research/research_fbg_r3_experiments/src/physics.py`: `MaxwellModel.state` stores `m/j/Q/Y/lu` and reuses exact-equal material states; `data_value_gradient(m,Y)` calls that state cache, solves one `ns`-RHS adjoint against its retained LU, and returns original normalized-data-plus-prior value/gradient.
- `research/research_z_alg_r4_diversity/src/physics.py`: `Model.forward`, `operation`, real-material smooth/TV `Regularizer.evaluate` and its original defaults `(smooth,tv,epsilon)=(1e-4,1e-4,.01)`.
- Actual R15 stateful output stores real `m` and `termination_reason`; C2 stores real `m` and `termination`. Both original CPU forms are retained separately from queue external child CPU.

One case model is shared by its11 endpoint validations. Actual repeated identical endpoint materials may reuse its current LU; that genuine zero-new-LU reuse is reported, never charged as another factorization. Each endpoint first calls original `state(m)`, then repeats `state(m.copy())` to check object/LU identity and zero new counter delta, then calls original `data_value_gradient`. Any unexpected cache miss or refactorization is charged and retained as failure. No `jacobian`, `tangent`, finite-step response or optimizer call exists. Model counters and operation CPU/wall timings are saved per endpoint; case totals include setup and failures and are summed once. Endpoint/case/whole-process CPU measurements are nested, not additive.

The diagnostic saves:

- Original physical `D=.5||Y_pred−Y||²/||Y||²`, `R`, `Phi`, residual, sigma/delta, residual/delta and residual/(sigma sqrt(number of complex ports)); missing/nonpositive noise denominators produce null, not synthesized thresholds.
- Full `g_data`, `g_prior`, `g_total` and their2/infinity norms. Data gradient is the inherited total gradient minus the unchanged prior gradient; only one adjoint is paid.
- Unit-scale box gradient mapping `m−clip(m−g_total,lower,upper)`,2/infinity norms, active masks/counts with inherited1e-8 reporting tolerance, and actual violations without clipping the material used for physics.
- Per-row safe NPZ containing `m`, physical `J/Q/Y_pred/Y_observed/residual`, all gradients, box mapping/masks and bounds. It contains no truth, optimizer state or proposed material.

Outputs are `REPORT.md`, `ENDPOINTS.csv/json`, `COST_EVENTS.csv`, `CASE_COSTS.json`, `FAILURES.json`, `PROGRESS.json`, `FIXED_PLAN.json`, `PARENT_SEAL_AT_ACTIONS.json`, per-row `states/*.npz`, and final **`RECEIPT.json`**. A whole-script interruption/error retains existing arrays and writes **`FAILED_RECEIPT.json`**, including charged active-case counters without adding them twice. No scientific gate or candidate-pool write is performed. A small gradient is explicitly not a local-minimum, nullspace, identifiability or recovery proof; parent interpretation remains required.
