# Historical six-dev remote handoff

Current update: root froze the core before generating24confirmation observations, then authorized remote eight-method execution. The remote recovery controller is active; see CURRENT_EXECUTION_HANDOFF.md for actual paths, attempt receipts and the corrected current package. The commands below preserve the earlier six-dev preparation record. They should not be rerun against existing paid outputs. No GT scoring has occurred.

Bundle: `R12_DEV6_BASELINE_PUBLIC_20261111.zip`. Extract its contents into a fresh directory such as `D:/AI/MAX_STRUCT_R12_20261111/confirmation`. It contains solver scripts/config and only public operators/observations. Sealed truth stays local. No old remote files are overwritten.

Existing read-only code: `D:/AI/MAX_STRUCT_R11_DEV_20261107`; existing Python: `D:/python/python.exe`; existing vendor: `D:/AI/MAX_STRUCT_R11_DEV_20261107/vendor`. No installs.

On the remote host, prepare absolute paths once (does not run a solver):

```powershell
$env:R12_LEGACY_REPO = 'D:/AI/MAX_STRUCT_R11_DEV_20261107'
& 'D:/python/python.exe' 'D:/AI/MAX_STRUCT_R12_20261111/confirmation/prepare_runtime.py' --legacy-repo 'D:/AI/MAX_STRUCT_R11_DEV_20261107'
```

This writes `RUNTIME_DEV6_BASELINE_TASKS.json` with30 mature baseline tasks, all public only. Repeated preparation refuses to overwrite. Example full CIE dev task; add `--cpu <root-allocated-core>` only after root allocates it:

```powershell
& 'D:/python/python.exe' 'D:/AI/MAX_STRUCT_R12_20261111/confirmation/run_baseline.py' --manifest 'D:/AI/MAX_STRUCT_R12_20261111/confirmation/RUNTIME_DEV6_BASELINE_TASKS.json' --config 'D:/AI/MAX_STRUCT_R12_20261111/confirmation/DEV6_BASELINE_SOLVER.json' --task-id 'r12_dev6_0000__cie_single__240CPU' --vendor-path 'D:/AI/MAX_STRUCT_R11_DEV_20261107/vendor'
```

Each of six case IDs r12_dev6_0000 through0005 has cie_single/cie_dual/global_dual/a_plus/fullG19. Existing complete paths refuse overwrite. DR is separately root-owned and intentionally excluded from this mature-baseline runtime manifest; root adds its frozen core through its own runner. Direct multistart/NN controls are not implemented here.

Baseline phase caps:144CPU total process prefix for proposal/setup, then96CPU corrector, global240CPU. Original native stage allowances remain500/500/500/1500 and original corrector3000 accepted steps; only actual paid trajectory prefixes may truncate. R11_CAP_PREDICTION and R11_INCUMBENT_TRAJECTORY prefixes retain inherited filenames, with studyR12 in receipts. External whole-child Windows GetProcessTimes should be captured by root's existing controller, including atomic overrun and serialization.

Data generation charged0.503048 absolute childCPU seconds,6 physical576LU,72 complex sourceRHS,0newoperatorbuilds. No inverseGT scoring authorized by this handoff.
