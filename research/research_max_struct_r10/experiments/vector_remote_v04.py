"""Finite foreground Windows v04 CPU deployment; no CUDA/no existing job edits."""
import time,subprocess,resource,sys,json,zipfile
from pathlib import Path
T=time.perf_counter();C=time.process_time();ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
OUT=ROOT/'delegated/vector/remote_v04_transport';OUT.mkdir(parents=True,exist_ok=False)
REMOTE='D:/AI/MAX_STRUCT_R10_VECTOR_20261107';HOST='jayzh@100.121.97.39';KEY='/Users/migodam/.ssh/heatcryst_xinan_ed25519'
SSH=['ssh','-i',KEY,'-o','BatchMode=yes','-o','ServerAliveInterval=15','-o','ServerAliveCountMax=3',HOST]
paths=[REPO/'research/research_max_struct_r9/src/vector_extension.py',*sorted((ROOT/'src').glob('vector_*.py')),ROOT/'experiments/vector_run.py',ROOT/'experiments/vector_windows_cpu_queue.py',*sorted((ROOT/'configs').glob('vector_*.json'))]
steps=[]
def run(command,label):
    start=time.perf_counter()
    with (OUT/(label+'.log')).open('w') as log:
        p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        for line in p.stdout:log.write(line);log.flush();print(line,end='',flush=True)
        status=p.wait()
    steps.append(dict(operation=label,exit=status,wall_seconds=time.perf_counter()-start,command=command))
    if status:raise RuntimeError(label+' failed '+str(status))
try:
    with zipfile.ZipFile(OUT/'CPU_SOURCE_UPDATE.zip','w',zipfile.ZIP_DEFLATED) as z:
        for path in paths:z.write(path,str(path.relative_to(REPO)))
    (OUT/'BUNDLE_MANIFEST.json').write_text(json.dumps(dict(files=[dict(path=str(p.relative_to(REPO)),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns) for p in paths],private_arrays_included=False),indent=2)+'\n')
    run(['scp','-i',KEY,str(OUT/'CPU_SOURCE_UPDATE.zip'),HOST+':'+REMOTE+'/CPU_SOURCE_UPDATE.zip'],'upload')
    code="import zipfile,time,json;from pathlib import Path;c=time.process_time();t=time.perf_counter();r=Path('"+REMOTE+"');zipfile.ZipFile(r/'CPU_SOURCE_UPDATE.zip').extractall(r);(r/'CPU_V04_EXTRACTION_COST.json').write_text(json.dumps(dict(platform='XINAN',CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t,scope='Own CPU source update only')))"
    run(SSH+['D:/python/python.exe -c "'+code+'"'],'extract')
    run(SSH+['D:/python/python.exe -u '+REMOTE+'/research/research_max_struct_r10/experiments/vector_windows_cpu_queue.py'],'connected_cpu_v04')
    status='COMPLETE'
except BaseException as exc:status='FAILED';print(repr(exc),flush=True)
finally:
    child=resource.getrusage(resource.RUSAGE_CHILDREN)
    (OUT/'TRANSPORT_COST.json').write_text(json.dumps(dict(status=status,platform='Mac',local_controller_CPU_seconds=time.process_time()-C,
        local_transport_child_CPU_seconds=child.ru_utime+child.ru_stime,wall_seconds=time.perf_counter()-T,
        steps=steps,scope='Own local connected SSH/SCP only; remote CPU scientific processes/helper CPU excluded',GPU_use=False),indent=2)+'\n')
if status!='COMPLETE':raise SystemExit(1)
