# R4 remote execution

The isolated remote project is `D:\AI\ZALG_R4_DIVERSITY_20261030\project` on the existing XINAN SSH account. No A9/A17 directory is used or changed. `hardware/environment.json` is the actual environment probe, obtained before any physics work. Preparation does not authorize the frozen experiment; the run manifest requires `execution_authorized: true` supplied after the parent GO.

Observed 2026-10-30: Python 3.10.9, NumPy 2.2.6, SciPy 1.15.3, PyTorch 2.6.0+cu124, CUDA build 12.4, CUDA available, RTX 4060 Laptop 8188 MiB. The initial probe saw an idle GPU; A9 work resumed later and is recorded separately. Only the missing `threadpoolctl==3.6.0` was installed into this isolated project's `vendor/`; global Python, Torch and CUDA were unchanged. CPU topology discovery used Windows `GetSystemCpuSetInformation`: eight group-0 physical P cores have logical sibling groups `[0,1]` through `[14,15]`; E-core logical IDs are 16–31.

The runner binds only its own workers with `SetProcessAffinityMask`. CPU has four workers on logical 0/2/4/6 (four physical P cores); CUDA has two workers on logical 8/10 (two other physical P cores). BLAS/OpenMP/Torch thread counts are 1. CUDA workers share device 0. This does not reserve those cores or alter any existing session. CPU and CUDA tasks/results retain separate backend labels; CUDA is a paired implementation comparison, not an extra algorithm arm. Successful affinity bindings and PIDs are recorded in `worker_affinity.json` when workers initialize.

Each task loads shared `operators_path` NPZ, then overlays the per-observation `public_path` NPZ. No private truth file is loaded. It calls `src.solver.solve(public, method, config, backend=backend)`. The solver owns exact LU/RHS accounting and its 45-second scientific budget. The runner synchronizes CUDA before/after solve for elapsed time and records return values, failures and tracebacks. Actual frozen task manifests use a 180-second external process timeout; this is an emergency process bound, not extra optimization time. An interrupted task remains `HARD_TIMEOUT` with unknown unrecovered physical counters. Import/load has a 120-second process bound outside scientific solve time. No extra numerical verification loop is run.

The parent supplies the final frozen manifest: CPU 72 observations × four algorithms; CUDA 18 paired observations × four algorithms. Task order is preserved within each pool. Return JSON must contain only the agreed lightweight solver results/scalar trace; the runner does not save truth or additional fields.

```json
{
  "execution_authorized": false,
  "solver_module": "src.solver",
  "output_dir": "results",
  "tasks": [
    {
      "task_id": "cpu_obs000_native_trf",
      "backend": "cpu",
      "operators_path": "public/operators.npz",
      "public_path": "public/obs000.npz",
      "method": "native_trf",
      "config": {"example_only": "replace with exact frozen solver configuration"},
      "hard_timeout_s": 180
    }
  ]
}
```

Run these commands from the local R4 root. Keep the existing SSH known-host policy; all connections use BatchMode and never read key contents.

```sh
scp -i "$SSH_IDENTITY" -o BatchMode=yes remote_runner.py "$REMOTE_HOST:$REMOTE_PROJECT/remote_runner.py"
scp -r -i "$SSH_IDENTITY" -o BatchMode=yes src public_data "$REMOTE_HOST:$REMOTE_PROJECT/"
scp -i "$SSH_IDENTITY" -o BatchMode=yes hardware/run_cpu.json hardware/run_cuda.json "$REMOTE_HOST:$REMOTE_PROJECT/hardware/"
```

The existing transport encodes PowerShell and avoids shell quoting ambiguity. Inspect activity immediately before launch; record unrelated new jobs and coordinate rather than stopping or altering them.

```sh
python /Volumes/migodam\'s-external-brain/Research/Inv_SLAM/Gaussian/A9/code/remote.py --script 'Get-CimInstance Win32_Process | Where-Object { $_.Name -match "python|pythonw|jupyter" } | Select-Object ProcessId,Name,CommandLine | ConvertTo-Json -Compress; nvidia-smi --query-gpu=name,memory.total,memory.used,utilization.gpu --format=csv,noheader'
```

Actual execution uses A9's foreground long-lived SSH transport. The first detached `Start-Process` attempt did not survive session closure and produced no prediction receipts; its unrecovered work is UNKNOWN (`audit/CPU_LAUNCH_CORRECTION.json`). It was replaced with the command below; no completed prediction was rerun.

```sh
python "$A9_TRANSPORT" --script '$r="D:\AI\ZALG_R4_DIVERSITY_20261030\project"; Set-Location $r; & D:\python\python.exe -u remote_runner.py run hardware/run_cpu.json' --timeout 7200
# Separate foreground session, once the GPU resource schedule allows it:
python "$A9_TRANSPORT" --script '$r="D:\AI\ZALG_R4_DIVERSITY_20261030\project"; Set-Location $r; & D:\python\python.exe -u remote_runner.py run hardware/run_cuda.json' --timeout 7200
```

Download after completion or for a checkpoint; individual task receipts are atomic:

```sh
scp -r -i "$SSH_IDENTITY" -o BatchMode=yes "$REMOTE_HOST:$REMOTE_PROJECT/results_cpu" ./
scp -r -i "$SSH_IDENTITY" -o BatchMode=yes "$REMOTE_HOST:$REMOTE_PROJECT/results_cuda" ./
```
