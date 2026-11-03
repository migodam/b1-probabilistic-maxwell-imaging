"""Read-only existing i9/4060 resource probe; no workload or service launch."""
from pathlib import Path
import importlib.util,json,time
ROOT=Path(__file__).resolve().parents[1]
workspace=ROOT.parents[5]
source=workspace/'Gaussian/A9/code/remote.py'
spec=importlib.util.spec_from_file_location('existing_remote',source)
remote=importlib.util.module_from_spec(spec);spec.loader.exec_module(remote)
script=r"""
$cpu=Get-CimInstance Win32_Processor | Select-Object Name,NumberOfCores,NumberOfLogicalProcessors
$jobs=Get-Process -Name python,pythonw -ErrorAction SilentlyContinue | Select-Object Id,ProcessName,CPU,WorkingSet64
$gpu=& nvidia-smi --query-gpu=name,memory.total,memory.used,utilization.gpu --format=csv,noheader
$gpuJobs=& nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv,noheader
[ordered]@{cpu=$cpu;python_processes=$jobs;gpu=$gpu;compute_processes=$gpuJobs;
cuda_library_present=(Test-Path 'D:\AI\ZALG_CUDA_PIPELINE_20261030\src\device.py');
python_present=(Test-Path 'D:\python\python.exe')} | ConvertTo-Json -Depth 5 -Compress
"""
start=time.perf_counter()
result=remote.powershell(script,timeout=25)
receipt=dict(machine='existing XINAN i9/RTX4060',scope='read-only activity and existing operator inventory',
             returncode=result.returncode,stdout=result.stdout.decode('utf8','replace'),
             stderr=result.stderr.decode('utf8','replace'),wall_seconds=time.perf_counter()-start,
             physical_LU=0,physical_RHS=0,launched_jobs=0)
(ROOT/'audit').mkdir(exist_ok=True)
(ROOT/'audit/COMPUTE_RESOURCE_PROBE.json').write_text(json.dumps(receipt,indent=2))
print(json.dumps(receipt))
