# Independent static review of BothBlocks source diff

The complete unified diff is `BOTHBLOCKS_SOURCE_DIFF.patch`. Inspection found changes confined to: the source scope comment; the new feedback_current_maps helper; the current transform and actual-loss receipt labels; replacing scalar current H/T/H* with the prepared matrix denominator; replacing the old scalar/onehop historical loss values by actual full moving loss values, with small-factor reuse and explicit cost counters. The full material callback itself is the same R10 matrix_material_readout alias; no material-rule edit was made in this control.

No change appears in uniform initialization, initial current, DPCS construction/rank selection, operator geometry, source/data metrics, fixed E/G/S, Fourier masks/projection, stage construction, beta values, coefficient transport, PR update/restart, exact quadratic line, current coefficient update, material callback arguments, original physical Phi validation, public-Phi incumbent selection, physical corrector or regularizer. NUMERICAL_SCHEDULE.json is exactly identical for each corresponding old G19 and new BothBlocks object. The separate auxiliary receipt/history is authoritative for actual current/material metric labels; the inherited pre-dispatch schedule dispatch flag remains onehop for the custom material callback.

This review is static and saved-artifact based, with zero new FW, adjoint or optimization. The one actual H adjoint/directional derivative check is recorded separately. The new current runner is identical to its prelaunch executed source snapshot: True. No R8 or frozen six-arm source was changed. The new wrapper's monkeypatch is local to its own child process; ordinary arms must use the ordinary runner.

Reproduction commands from repository root use the existing Python environment and:

```
python research/research_max_struct_r11/experiments/diagnose_bothblocks_run.py --manifest research/research_max_struct_r11/configs/twod_dev_bothblocks_tasks.json --config research/research_max_struct_r11/configs/twod_dev_solver.json --task-id r10_twod_0022__fullGreenBothBlocks__480CPU
```

Output paths are immutable run destinations; reproductions require a separate manifest with new output directories. The executed startup manifest/config/source copies are in bothblocks_run_logs/. No historical source path or command was invented.
