"""Finite connected transport/controller; no truth and no detached jobs."""
import subprocess,sys,time,json,resource,zipfile
from pathlib import Path
ENTRY=time.perf_counter();CPU=time.process_time();ROOT=Path(__file__).resolve().parents[3];REPO=ROOT.parents[1];HERE=Path(__file__).parent
REMOTE='D:/AI/MAX_STRUCT_R9_PHYSREF_20261106';HOST='jayzh@100.121.97.39';KEY='/Users/migodam/.ssh/heatcryst_xinan_ed25519'
SSH=['ssh','-i',KEY,'-o','BatchMode=yes','-o','ServerAliveInterval=15','-o','ServerAliveCountMax=3',HOST]
paths=['research/research_max_struct_r1/src/a9_green.py','research/research_max_struct_r1/src/vector_backend.py','research/research_max_struct_r8/src/vector_cie3d.py']
paths+=['research/research_max_struct_r9/'+p for p in ['src/vector_original_extension.py','src/vector_original_cuda.py','src/vector_physical_refinement.py','experiments/run_original_physical_refinement.py','configs/original_physical_refinement.json','delegated/3d_extension/windows_resolved.py','RESULTS/original_vector3d/cohorts/matched_grid/public/observations.npz','RESULTS/original_vector3d/cohorts/matched_grid/inversions/GlobalDual/stage2.npz','RESULTS/original_vector3d/cohorts/matched_grid/inversions/AplusGlobal/stage2.npz']]
steps=[]
def run(cmd,label):
    t=time.perf_counter();r=subprocess.run(cmd);steps.append(dict(operation=label,exit=r.returncode,wall_seconds=time.perf_counter()-t));
    if r.returncode:raise RuntimeError(label+' failed '+str(r.returncode))
try:
    with zipfile.ZipFile(HERE/'PUBLIC_BUNDLE.zip','w',zipfile.ZIP_DEFLATED) as z:
        for p in paths:z.write(REPO/p,p)
    (HERE/'BUNDLE_MANIFEST.json').write_text(json.dumps(dict(files=[dict(path=p,bytes=(REPO/p).stat().st_size) for p in paths],GT_included=False),indent=2))
    run(SSH+['D:/python/python.exe -c "from pathlib import Path; p=Path(\''+REMOTE+'\'); p.mkdir(exist_ok=False)"'],'REMOTE_ISOLATED_DIRECTORY')
    run(['scp','-i',KEY,str(HERE/'PUBLIC_BUNDLE.zip'),HOST+':'+REMOTE+'/PUBLIC_BUNDLE.zip'],'PUBLIC_BUNDLE_UPLOAD')
    code="import zipfile,time,json;from pathlib import Path;t=time.perf_counter();c=time.process_time();p=Path('"+REMOTE+"');zipfile.ZipFile(p/'PUBLIC_BUNDLE.zip').extractall(p);(p/'EXTRACTION_COST.json').write_text(json.dumps(dict(platform='XINAN',CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t,operation='PUBLIC_BUNDLE_EXTRACTION',additive=True)))"
    run(SSH+['D:/python/python.exe -c "'+code+'"'],'REMOTE_EXTRACT')
    t=time.perf_counter();cmd=SSH+['D:/python/python.exe -u '+REMOTE+'/research/research_max_struct_r9/experiments/run_original_physical_refinement.py']
    (HERE/'LAUNCH.json').write_text(json.dumps(dict(command=cmd,GT_access=False,started_unix=time.time()),indent=2))
    with (HERE/'CONNECTED_RUN.log').open('w') as log:
        p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        for line in p.stdout:log.write(line);log.flush();print(line,end='',flush=True)
        status=p.wait()
    steps.append(dict(operation='CONNECTED_GPU_PROCESS',exit=status,wall_seconds=time.perf_counter()-t))
    if status:raise RuntimeError('Connected GPU process failed '+str(status))
    status='COMPLETE'
except BaseException as exc:
    status='FAILED';print(repr(exc),flush=True)
finally:
    child=resource.getrusage(resource.RUSAGE_CHILDREN)
    (HERE/'TRANSPORT_CONTROLLER_COST.json').write_text(json.dumps(dict(status=status,platform='Mac',local_controller_CPU_seconds=time.process_time()-CPU,local_transport_child_CPU_seconds=child.ru_utime+child.ru_stime,wall_seconds=time.perf_counter()-ENTRY,steps=steps,remote_GPU_process_CPU_included=False,remote_extraction_CPU_included=False,GT_access=False),indent=2))
if status!='COMPLETE':raise SystemExit(1)
