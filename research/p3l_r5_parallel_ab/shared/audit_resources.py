"""Read-only inventory/header audit; never reads all large payloads or secrets."""
import json,os,platform,subprocess,zipfile,shutil,statistics
from pathlib import Path
from .physics_adapter import ROOT,HISTORY,SHARED,OUT,write
import numpy as np

def command(args):
    try:return subprocess.check_output(args,text=True,stderr=subprocess.STDOUT).strip()
    except Exception as e:return repr(e)
def header(path):
    result={}
    with zipfile.ZipFile(path) as archive:
        for info in archive.infolist():
            if not info.filename.endswith('.npy'):continue
            with archive.open(info) as f:
                version=np.lib.format.read_magic(f)
                shape,order,dtype=(np.lib.format.read_array_header_1_0(f) if version==(1,0) else np.lib.format.read_array_header_2_0(f))
            result[info.filename[:-4]]=dict(shape=list(shape),dtype=str(dtype),compressed_bytes=info.compress_size,array_bytes=info.file_size)
    return result

def main():
    old=HISTORY/'research/p3l_r4_conditional_transport'
    asset=json.loads((old/'AUDIT/risk/ASSET_INVENTORY.json').read_text())
    assets=[]
    for name,item in asset['assets'].items():
        path=Path(item['path']);record=dict(name=name,path=str(path),exists=path.exists(),historical_report_bytes=item['bytes'])
        if path.is_file():record.update(actual_bytes=path.stat().st_size,fields=header(path),read_scope='NPZ central directory and NPY headers only')
        assets.append(record)
    contexts=[]
    for group in ('dev','locked','confirmation','train','validation'):
        folder=old/'results/contexts'/group
        costs=[]
        if folder.exists():
            for path in folder.glob('*.json'):
                d=json.loads(path.read_text())
                if 'new_RHS' in d:costs.append(d['new_RHS'])
        contexts.append(dict(group=group,parents_with_cost_receipts=len(costs),RHS_min=min(costs) if costs else None,RHS_median=statistics.median(costs) if costs else None,RHS_max=max(costs) if costs else None,RHS_sum=sum(costs)))
    checkpoints=[]
    for directory in (old/'models',old/'I',old/'II',old/'III',old/'IV'):
        if directory.exists():
            for path in directory.rglob('*.pt'):
                st=path.stat();checkpoints.append(dict(path=str(path.relative_to(HISTORY)),bytes=st.st_size,mtime=st.st_mtime,read_scope='stat only; weights not deserialized'))
    write(OUT/'audit/HISTORICAL_ASSET_AUDIT.json',dict(historical_source_commit='85b7bcbf39f0eeab7062ecce17531538151082bf',delivery_cleanup_commit='092da74e',large_arrays=assets,checkpoints=checkpoints,checkpoint_count=len(checkpoints),checkpoint_bytes=sum(x['bytes'] for x in checkpoints),contexts=contexts,context_algorithm=dict(source='tools/run_r4_contexts.py',native_cap=96,probit_cap=64,source_rule='frozen common context budget96+64',upper_RHS_per_parent=6400,full512_upper_RHS=3276800,qualification='Actual optimizer cost is data dependent; 40RHS supplied-x forward+J is not equivalent to a strong optimized context. New lower-cost context is a separate regime.'),version_diff_stat=command(['git','diff','85b7bcbf','092da74e','--stat']),no_large_payloads_read=True,no_SHA256=True))
    disk=shutil.disk_usage(ROOT)
    try:
        import torch
        torchdata=dict(version=torch.__version__,MPS_available=torch.backends.mps.is_available(),MPS_built=torch.backends.mps.is_built())
    except ImportError:torchdata=dict(status='NOT_INSTALLED_IN_RUNTIME')
    write(OUT/'audit/HARDWARE.json',dict(platform=platform.platform(),machine=platform.machine(),CPU_model=command(['sysctl','-n','machdep.cpu.brand_string']),CPU_logical=command(['sysctl','-n','hw.logicalcpu']),CPU_physical=command(['sysctl','-n','hw.physicalcpu']),memory_bytes=command(['sysctl','-n','hw.memsize']),disk=dict(total_bytes=disk.total,used_bytes=disk.used,free_bytes=disk.free),python=platform.python_version(),numpy=np.__version__,torch=torchdata,physics_device='CPU',physics_workers=1,BLAS_threads_requested=1,BLAS_backend='Apple Accelerate in numpy and scipy local .venv_nn',thread_environment={k:os.environ.get(k) for k in ('VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS')},no_secrets_read=True))
if __name__=='__main__':main()
