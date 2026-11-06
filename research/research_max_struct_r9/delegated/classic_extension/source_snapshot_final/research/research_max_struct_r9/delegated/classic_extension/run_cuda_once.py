"""One explicitly authorized GPU-slot job and complete public-result fetch."""
import time
START,CPU_START=time.time(),time.process_time()
import argparse,json,resource,shutil,subprocess,zipfile
from pathlib import Path
REPO=Path(__file__).resolve().parents[4];ROOT=REPO/'research/research_max_struct_r9';A=ROOT/'delegated/classic_extension'


def main():
    p=argparse.ArgumentParser();p.add_argument('--host',required=True);p.add_argument('--identity',required=True);a=p.parse_args()
    remote_root='D:/AI/MAX_STRUCT_R9_CLASSIC_20261106';base=['ssh','-i',a.identity,'-o','BatchMode=yes','-o','ConnectTimeout=10',a.host,'D:/python/python.exe','-B','-'];before=resource.getrusage(resource.RUSAGE_CHILDREN)
    candidate=ROOT/'data_public/classic_extension/cuda_candidate_material.npz'
    shutil.copy2(ROOT/'RESULTS/classic_extension/austria_eps3p5_400MHz__cie_single/CHECKPOINT_080.npz',candidate)
    paths=[candidate,A/'cuda_callback_benchmark.py',REPO/'research/research_z_alg_cuda_pipeline/src/device.py']
    bundle=A/'cuda_callback_code_input.zip'
    with zipfile.ZipFile(bundle,'w',zipfile.ZIP_DEFLATED,compresslevel=1) as z:
        for f in paths:z.write(f,f.relative_to(REPO))
    def ssh(code,log=None):
        proc=subprocess.run(base,input=code,text=True,stdout=log or subprocess.PIPE,stderr=subprocess.STDOUT if log else subprocess.PIPE)
        if proc.returncode:raise RuntimeError('Own GPU-slot helper exited nonzero; see retained log')
        return proc.stdout
    def scp(s,t):
        proc=subprocess.run(['scp','-i',a.identity,'-o','BatchMode=yes',s,t],capture_output=True,text=True)
        if proc.returncode:raise RuntimeError(proc.stderr)
    result={}
    try:
        scp(str(bundle),a.host+':'+remote_root+'/cuda_callback_bundle.zip')
        ssh("from pathlib import Path\nimport zipfile\nr=Path("+repr(remote_root)+")\nwith zipfile.ZipFile(r/'cuda_callback_bundle.zip') as z:z.extractall(r)\n")
        with (A/'CUDA_CALLBACK_CONNECTED.log').open('w') as log:
            ssh("import os,runpy\nos.chdir("+repr(remote_root)+")\nrunpy.run_path('research/research_max_struct_r9/delegated/classic_extension/cuda_callback_benchmark.py',run_name='__main__')\n",log)
        result['benchmark_process_returncode']=0
    except BaseException as exc:
        result.update(status='FAILED',error=repr(exc))
    # Fetch failure/progress artifacts too, after process exit released context.
    try:
        info=ssh("from pathlib import Path\nimport zipfile,json,subprocess,time\nc=time.process_time();r=Path("+repr(remote_root)+");o=r/'research/research_max_struct_r9/RESULTS/classic_extension_cuda_callback'\nwith zipfile.ZipFile(r/'cuda_callback_results.zip','w',zipfile.ZIP_DEFLATED,compresslevel=1) as z:\n for f in o.rglob('*'):\n  if f.is_file():z.write(f,f.relative_to(r))\napps=subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,used_gpu_memory','--format=csv,noheader'],text=True).strip()\nmem=subprocess.check_output(['nvidia-smi','--query-gpu=memory.used','--format=csv,noheader,nounits'],text=True).strip()\nprint(json.dumps(dict(compute_apps_after_exit=apps,GPU_memory_MiB_after_exit=mem,remote_fetch_helper_CPU_seconds=time.process_time()-c)))\n")
        result.update(json.loads(info.lstrip('\ufeff')))
        archive=A/'cuda_callback_results.zip';scp(a.host+':'+remote_root+'/cuda_callback_results.zip',str(archive))
        with zipfile.ZipFile(archive) as z:
            for name in z.namelist():
                if not name.startswith('research/research_max_struct_r9/RESULTS/classic_extension_cuda_callback/'):
                    raise ValueError('Unexpected archive member')
            z.extractall(REPO)
        result.setdefault('status','COMPLETED')
    except BaseException as exc:result.update(fetch_error=repr(exc),status='FAILED')
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    result.update(UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),GT_access='NONE',GPU_slot='parent authorized exclusive; hostP12',
        wall_seconds=time.time()-START,local_controller_CPU_seconds=time.process_time()-CPU_START,
        local_transport_child_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
        archive_bytes=bundle.stat().st_size,source_paths=[str(f.relative_to(REPO)) for f in paths],
        CPU16_changed=False,credentials_serialized=False)
    (A/'CUDA_TRANSPORT_RECEIPT.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result),flush=True)
    if result['status']!='COMPLETED':raise SystemExit(1)


if __name__=='__main__':main()
