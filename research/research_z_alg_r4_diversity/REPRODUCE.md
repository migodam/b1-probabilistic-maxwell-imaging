# R4 reproduction and offline analysis

Run from this repository's R4 directory. The archived observations, frozen solver/configuration and closed result receipts define the actual study. Reusing `public_data/` reproduces those exact observations. Running the generator makes new random backgrounds, complex noise directions, anonymous IDs and public order; it is a new realization, not a byte-for-byte regeneration of the archived study. Use a separate clean checkout/output location for new data. Do not overwrite the archived realization or reuse its predictions as if they solved newly generated observations.

The CPU comparison has 12 shapes × three contrasts × two noise levels × four methods = 288 planned inference tasks. CUDA has the fixed digit2/digit7/ring subset: 72 tasks. Compare algorithms within each backend. The CUDA subset does not add independent objects to the 12 CPU shapes, and repeated contrast/noise variants do not turn 12 shapes into 72 independent objects. CUDA availability and shared workload timing do not establish a hardware speedup claim.

## Environment and archived data

Provide your own existing Python environments, SSH account/key and isolated remote project directory. No private host or key path is embedded here. The Windows environment recorded in `hardware/environment.json` used Python 3.10.9, NumPy 2.2.6, SciPy 1.15.3 and Torch 2.6.0+cu124, with i9-14900HX P cores and RTX4060 Laptop complex128 CUDA. Local generation additionally uses Pillow and the documented local Arial Bold / Georgia Bold font files. Plotting uses Matplotlib. The small `threadpoolctl` dependency was isolated under the remote project's `vendor/` after the recorded import failure; copy that archived directory if reproducing the same deployment. These instructions do not install CUDA or modify other sessions.

```sh
export LOCAL_PYTHON='/absolute/path/to/local/python'
export R4_SSH_HOST='user@your-existing-host'
export R4_SSH_KEY='/absolute/path/to/existing/private-key'
export R4_REMOTE_PROJECT='D:\AI\your-isolated-r4-project'
export R4_REMOTE_SCP='D:/AI/your-isolated-r4-project'
export R4_REMOTE_PYTHON='D:\python\python.exe'
```

Use the archived `public_data/` directly for exact observation replay. The optional command below is only for a newly authorized realization in a clean location, after freezing solver/config and preserving the local-only truth boundary:

```sh
"$LOCAL_PYTHON" experiments/generate_diversity.py --go \
  --operators ../research_z_alg_r3_digits/public_data/operators.npz
```

The generator uses 36 full nonlinear 576-dimensional VIE forward factorizations and 432 source RHS. Clean/noisy variants share each forward. The private mapping retains the measured Born/full relative difference from matrix multiplication without another physical solve. No contrast level is assumed to be Born-valid. `private_truth/` must remain local and unread by inference workers until the parent closes planned results.

## Remote CPU and CUDA inference

For a single anonymous archived observation, the same frozen API can be run directly on a CPU without the Windows orchestration. Set `OBSERVATION` to a path from the public manifest; it contains no digit label. This is a replay, outside the archived benchmark timing:

```sh
"$LOCAL_PYTHON" run_one.py --observation "$OBSERVATION" \
  --method globaldual --backend cpu --output replay_globaldual.json
# Other methods: native_trf, reciprocal_trf, data_lbfgs.
```

Create the isolated remote project directory using your existing SSH transport first. Upload only public observations/operators, frozen source/configuration, runner, manifests and the archived isolated dependency directory. Do not recursively upload the whole project or `private_truth/`.

```sh
scp -i "$R4_SSH_KEY" -o BatchMode=yes -o ConnectTimeout=8 \
  remote_runner.py "$R4_SSH_HOST:$R4_REMOTE_SCP/remote_runner.py"
scp -r -i "$R4_SSH_KEY" -o BatchMode=yes -o ConnectTimeout=8 \
  src configs public_data vendor "$R4_SSH_HOST:$R4_REMOTE_SCP/"
scp -i "$R4_SSH_KEY" -o BatchMode=yes -o ConnectTimeout=8 \
  hardware/run_cpu.json hardware/run_cuda.json \
  "$R4_SSH_HOST:$R4_REMOTE_SCP/hardware/"
```

Use the frozen manifests without increasing scientific budgets: 45 seconds, LU 2000, full-system RHS 60000 per task, shared bounds/regularization and 8→12→24 material stages. The process guard in the archived manifest is operational overhead protection, not a larger scientific solve budget. The runner binds only its workers: four CPU workers on distinct P cores and two CUDA workers on separate P-core host threads, with one BLAS/Torch thread each. Never change another session's affinity or stop its jobs. Queue CUDA while an existing workload occupies the GPU; preserve that scheduling in the audit.

The following foreground SSH command keeps the Windows process attached to the live SSH session. The initial detached `Start-Process` attempt did not survive its SSH session; its unknown unrecovered costs remain in the ledger. Set `R4_BACKEND=cpu` first, then run the same block with `R4_BACKEND=cuda` after the existing GPU workload has cleared. This runs inference only from public data and never opens truth. Leave the controlling terminal/process alive until the command returns.

```sh
export R4_BACKEND=cpu  # later: cuda
"$LOCAL_PYTHON" - <<'PY'
import base64, os, subprocess, sys
def ps_quote(s): return "'" + s.replace("'", "''") + "'"
backend = os.environ['R4_BACKEND']
script = ("$ErrorActionPreference='Stop'; Set-Location "
          + ps_quote(os.environ['R4_REMOTE_PROJECT']) + '; & '
          + ps_quote(os.environ['R4_REMOTE_PYTHON'])
          + ' -u remote_runner.py run '
          + ps_quote('hardware/run_' + backend + '.json') + '; exit $LASTEXITCODE')
encoded = base64.b64encode(script.encode('utf-16le')).decode()
command = ['ssh','-i',os.environ['R4_SSH_KEY'],'-o','BatchMode=yes',
           '-o','ConnectTimeout=8',os.environ['R4_SSH_HOST'],
           'powershell -NoProfile -EncodedCommand ' + encoded]
sys.exit(subprocess.run(command).returncode)
PY
```

Download the remote `results_cpu/` and `results_cuda/` directories after the scheduled backend completes. Preserve existing receipts; do not automatically restart completed tasks or rerun missing tasks to obtain a complete-looking table.

```sh
mkdir -p downloaded
scp -r -i "$R4_SSH_KEY" -o BatchMode=yes -o ConnectTimeout=8 \
  "$R4_SSH_HOST:$R4_REMOTE_SCP/results_cpu" ./downloaded/
scp -r -i "$R4_SSH_KEY" -o BatchMode=yes -o ConnectTimeout=8 \
  "$R4_SSH_HOST:$R4_REMOTE_SCP/results_cuda" ./downloaded/
```

If CUDA is explicitly closed as NOT_RUN, skip its download. Collect the existing task receipts into the scorer's single input directory; absent tasks are left absent:

```sh
"$LOCAL_PYTHON" - <<'PY'
import json, pathlib, shutil
output=pathlib.Path('remote_results');output.mkdir(exist_ok=True)
tasks=json.loads(pathlib.Path('hardware/run_manifest.json').read_text())['tasks']
for task in tasks:
    source=pathlib.Path('downloaded')/('results_'+task['backend'])/(task['task_id']+'.json')
    if source.exists():shutil.copyfile(source,output/source.name)
PY
```

## Close, score and figures

Only after the parent declares execution closed, record all 360 planned dispositions. The seal command does not run physics, retry missing tasks, or require successful completion. Missing receipts become `NOT_RUN`; explicit failures and unknown cost counters remain visible. The parent preserves the closed raw receipts before truth disclosure.

```sh
"$LOCAL_PYTHON" score_results.py seal \
  --manifest hardware/run_manifest.json --results remote_results \
  --output audit/ALL_PREDICTIONS_SEALED.json
"$LOCAL_PYTHON" score_results.py score \
  --seal audit/ALL_PREDICTIONS_SEALED.json --truth private_truth \
  --output RESULTS_SUMMARY.json
"$LOCAL_PYTHON" make_figures.py --scores RESULTS_SUMMARY.json \
  --truth private_truth --out figures
```

Scoring independently reports material RMSE, NRMSE=`RMSE/A`, MAE, reconstruction-only Otsu IoU and separately labeled offline known-midpoint IoU. Frozen gates are NRMSE≤.15 and Otsu IoU≥.7 for material, relative data residual≤.001 clean / ≤.003 noisy for data, and their conjunction for joint success. The offline oracle-best-K covers finite correction-phase optimizer candidates; these are not independent branches and do not replace the public-physical-score top-1 endpoint.

## Cost ledger and paired shape reports

Run reporting only after scoring. It opens the scored summary and result/audit receipts, imports no solver and accesses no truth directory. Raw result materials are not used in reporting; scalar top-1 metadata and counters are extracted from the actual solver schema.

```sh
"$LOCAL_PYTHON" analyze_results.py --summary RESULTS_SUMMARY.json \
  --receipts remote_results --audit audit --out analysis
```

`--cost-notes audit/ADDITIONAL_COST_ROWS.json` is optional and repeatable. Supply measured scalar setup/DEV/import/launch records, never invented estimates. Its JSON is an array or `{"rows":[...]}`; for example a DEV row may carry `category`, `backend`, `status`, `source`, `LU`, `actual_complex_RHS`, `Gram_RHS`, `wall_seconds`, `cpu_seconds`, `costs_status` and `note`. Omit unavailable numeric fields, and retain unknown setup/launch costs as unknown. A successful DEV run is a cost outside production comparisons, not another benchmark observation. The existing generation, import-failure and failed-detached-launch audits plus `audit/DEV_CPU.json` and `hardware/dev_gpu_actual/dev_cuda_globaldual.json` are recognized automatically. Avoid duplicating them through cost notes. Missing overhead timing gets explicit unknown rows rather than zero.

Outputs are `COST_LEDGER.csv`, `PAIRED_VARIANTS.csv`, `PAIRED_BY_BACKEND_AMPLITUDE_NOISE.csv`, `SHAPE_CLUSTERS.csv`, `SHAPE_CLUSTER_SUMMARY.csv` and `ANALYSIS.json`. The ledger preserves failed/unknown/NOT_RUN records, exact full-system/source/adjoint/receiver-transpose RHS, small-Gram RHS/factorizations, FullJ counts, wall/CPU timing, top-1 stage/origin, `all_unknowns_released` and stop reasons. Small Gram RHS are not added to full-system RHS. Solver, worker-total and operation times overlap and must not be summed together.

GlobalDual is paired against each comparator within backend/contrast/noise using `abs(delta NRMSE)<=max(.02,.10*min(left NRMSE,right NRMSE))` as practical ties; every continuous delta and raw RMSE is retained alongside joint gates. Shape-cluster rows average matched contrast/noise variants within each shape; partial shapes retain missing counts. The final equal-shape summary uses only complete shape clusters, with complete/partial counts disclosed. CUDA clusters remain the same three paired objects already represented in the CPU collection. These are descriptive summaries without invented significance or a claim that variants are independent samples.
