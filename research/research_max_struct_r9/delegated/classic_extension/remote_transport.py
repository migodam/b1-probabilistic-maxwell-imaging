"""Own finite classic cohort transport; no truth, services, GPU or other jobs.

Credentials are CLI-only and are not written into receipts or archives.
"""
import time
START,CPU_START=time.time(),time.process_time()
import argparse
import ast
import json
from pathlib import Path
import resource
import subprocess
import zipfile

REPO=Path(__file__).resolve().parents[4]
ROOT=REPO/'research/research_max_struct_r9'
AUDIT=ROOT/'delegated/classic_extension'
SOURCE_AUDIT=AUDIT


def dump(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False))


def main():
    global AUDIT
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['prepare','run','status','fetch'])
    ap.add_argument('--host',required=True);ap.add_argument('--identity',required=True)
    ap.add_argument('--remote-root',default='D:/AI/MAX_STRUCT_R9_CLASSIC_20261106')
    ap.add_argument('--local-output-root',type=Path,help='New local cohort destination; keeps historical RESULTS untouched')
    a=ap.parse_args();before=resource.getrusage(resource.RUSAGE_CHILDREN)
    local_output=(a.local_output_root if a.local_output_root.is_absolute() else REPO/a.local_output_root) if a.local_output_root else None
    if local_output is not None:
        AUDIT=local_output/'transport';AUDIT.mkdir(parents=True,exist_ok=True)
    ssh=['ssh','-i',a.identity,'-o','BatchMode=yes','-o','ConnectTimeout=10',a.host,'D:/python/python.exe','-B','-']
    def remote(code):
        # Metadata/zip helpers also remain inside this worker's allocation.
        prelude="import psutil\npsutil.Process().cpu_affinity([0,2,4,6])\n"
        p=subprocess.run(ssh,input=prelude+code,text=True,capture_output=True,timeout=900)
        if p.returncode:raise RuntimeError(p.stderr+'\n'+p.stdout)
        return json.loads(p.stdout.lstrip('\ufeff'))
    def scp(source,target):
        p=subprocess.run(['scp','-i',a.identity,'-o','BatchMode=yes',source,target],capture_output=True,text=True,timeout=900)
        if p.returncode:raise RuntimeError(p.stderr)
    root=repr(a.remote_root)
    result={}
    try:
        if a.action=='prepare':
            seeds=[ROOT/'src/classic_extension.py',ROOT/'experiments/classic_extension.py',SOURCE_AUDIT/'remote_fifo.py']
            found=set();pending=list(seeds)
            while pending:
                p=pending.pop()
                if p in found:continue
                found.add(p)
                for node in ast.walk(ast.parse(p.read_text())):
                    mods=([node.module] if isinstance(node,ast.ImportFrom) and node.module else
                          [n.name for n in node.names] if isinstance(node,ast.Import) else [])
                    for name in mods:
                        if name.startswith('research.'):
                            target=REPO/Path(*name.split('.')).with_suffix('.py')
                            if target.exists():pending.append(target)
            public=[p for p in (ROOT/'data_public/classic_extension').rglob('*') if not p.name.startswith('cuda_')]
            files=sorted(found|{p for p in public if p.is_file()}|
                {ROOT/'configs/classic_extension_main.json',ROOT/'configs/classic_extension_paid_budgets.json'})
            archive=AUDIT/'classic_xinan_code_public_only.zip'
            with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=1) as z:
                for p in files:z.write(p,p.relative_to(REPO))
            result=remote("from pathlib import Path\nimport json\nr=Path("+root+")\nr.mkdir(parents=True,exist_ok=True)\nprint(json.dumps(dict(root=str(r),exists=True)))")
            scp(str(archive),a.host+':'+a.remote_root+'/classic_bundle.zip')
            result.update(remote("from pathlib import Path\nimport json,zipfile,sys\nr=Path("+root+")\nwith zipfile.ZipFile(r/'classic_bundle.zip') as z:z.extractall(r)\nprint(json.dumps(dict(extracted=True,python=sys.executable)))"))
            result.update(archive_bytes=archive.stat().st_size,files=[str(p.relative_to(REPO)) for p in files],
                          GT_uploaded=False,allocated_Pcores=[0,2,4,6],GPU=False)
        elif a.action=='run':
            code="import os,runpy\nos.chdir("+root+")\nrunpy.run_path('research/research_max_struct_r9/delegated/classic_extension/remote_fifo.py',run_name='__main__')\n"
            with (AUDIT/'REMOTE_CONNECTED_QUEUE.log').open('w') as log:
                p=subprocess.run(ssh,input=code,text=True,stdout=log,stderr=subprocess.STDOUT)
            result=dict(returncode=p.returncode,connected_SSH_kept_until_finite_FIFO_exit=True)
            if p.returncode:raise RuntimeError('Remote connected queue exited nonzero; inspect own log')
        else:
            code="""import time
helper_cpu,helper_wall=time.process_time(),time.perf_counter()
from pathlib import Path
import json,zipfile
r=Path(REMOTE_ROOT);out=r/'research/research_max_struct_r9/RESULTS/classic_extension_xinan';q=out/'queue'
rows=json.loads((q/'TASK_PROCESS_RECEIPTS.json').read_text()) if (q/'TASK_PROCESS_RECEIPTS.json').exists() else []
done=json.loads((q/'QUEUE_COMPLETE.json').read_text()) if (q/'QUEUE_COMPLETE.json').exists() else None
answer=dict(completed=len(rows),failed=sum(x['returncode']!=0 for x in rows),queue_complete=done is not None,last=rows[-2:])
progress=[]
for p in out.glob('*/PROGRESS.json'):
 if not (p.parent/'RECEIPT.json').exists() and not (p.parent/'FAILURE.json').exists():
  d=json.loads(p.read_text());progress.append(dict(task_id=p.parent.name,phase=d['phase'],stage=d.get('stage'),iteration=d.get('iteration'),accepted=d.get('accepted'),CPU_seconds=d['process_cpu_seconds']))
answer['active_progress']=progress
if FETCH:
 archive=r/'classic_public_predictions.zip'
 with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=1) as z:
  for p in out.rglob('*'):
   if p.is_file():z.write(p,p.relative_to(r))
 answer['archive_bytes']=archive.stat().st_size
answer['remote_helper_CPU_seconds']=time.process_time()-helper_cpu
answer['remote_helper_wall_seconds']=time.perf_counter()-helper_wall
print(json.dumps(answer))
""".replace('REMOTE_ROOT',root).replace('FETCH',repr(a.action=='fetch'))
            result=remote(code)
            if a.action=='fetch':
                archive=AUDIT/'classic_public_predictions_xinan.zip'
                scp(a.host+':'+a.remote_root+'/classic_public_predictions.zip',str(archive))
                with zipfile.ZipFile(archive) as z:
                    prefix=Path('research/research_max_struct_r9/RESULTS/classic_extension_xinan')
                    target_root=local_output if local_output is not None else REPO/prefix
                    for name in z.namelist():
                        if not name.startswith('research/research_max_struct_r9/RESULTS/classic_extension_xinan/'):
                            raise ValueError('Unexpected remote archive member')
                        target=target_root/Path(name).relative_to(prefix)
                        if target.exists():raise FileExistsError('Refuse to overwrite existing prediction/artifact: '+str(target))
                    for member in z.infolist():
                        if member.is_dir():continue
                        target=target_root/Path(member.filename).relative_to(prefix)
                        target.parent.mkdir(parents=True,exist_ok=True)
                        with z.open(member) as source,target.open('wb') as output:
                            import shutil
                            shutil.copyfileobj(source,output)
                result['local_output_root']=str(target_root.relative_to(REPO)) if target_root.is_relative_to(REPO) else str(target_root)
        result.update(status='COMPLETED',action=a.action,remote_alias='authorized_XINAN',
            remote_root=a.remote_root,UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
            epoch_start=START,epoch_end=time.time(),wall_seconds=time.time()-START,
            local_controller_CPU_seconds=time.process_time()-CPU_START,
            local_transport_child_CPU_seconds=(resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime+
                resource.getrusage(resource.RUSAGE_CHILDREN).ru_stime-before.ru_utime-before.ru_stime),
            local_CPU_scope='controller+its SSH/SCP subprocesses only; remote compute counted by GetProcessTimes',
            GT_access='NONE')
        path=AUDIT/('REMOTE_'+a.action.upper()+'_RECEIPT.json')
        if path.exists():
            history=json.loads(path.read_text());history=history if isinstance(history,list) else [history]
            history.append(result);dump(path,history)
        else:dump(path,result)
        print(json.dumps({k:v for k,v in result.items() if k!='files'}),flush=True)
    except BaseException as exc:
        dump(AUDIT/('REMOTE_'+a.action.upper()+'_FAILURE.json'),dict(status='FAILED',error=repr(exc),
            action=a.action,remote_alias='authorized_XINAN',GT_access='NONE',
            local_controller_CPU_seconds=time.process_time()-CPU_START,
            wall_seconds=time.time()-START))
        raise


if __name__=='__main__':main()
