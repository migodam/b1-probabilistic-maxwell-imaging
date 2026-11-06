"""Finite connected first vector cohort; local transport CPU is separate."""
import time,resource,subprocess,sys,json,zipfile
from pathlib import Path
T=time.perf_counter();C=time.process_time();ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];OUT=ROOT/'delegated/vector/remote_first';OUT.mkdir(parents=True,exist_ok=False)
REMOTE='D:/AI/MAX_STRUCT_R10_VECTOR_20261107';HOST='jayzh@100.121.97.39';KEY='/Users/migodam/.ssh/heatcryst_xinan_ed25519'
SSH=['ssh','-i',KEY,'-o','BatchMode=yes','-o','ServerAliveInterval=15','-o','ServerAliveCountMax=3',HOST]
paths=['research/research_max_struct_r1/src/a9_green.py','research/research_max_struct_r1/src/vector_backend.py','research/research_max_struct_r8/src/vector_cie3d.py','research/research_max_struct_r9/src/vector_original_extension.py','research/research_max_struct_r9/src/vector_original_cuda.py','research/research_max_struct_r9/delegated/3d_extension/windows_resolved.py']
paths += [str(p.relative_to(REPO)) for p in sorted((ROOT/'src').glob('vector_*.py'))]
paths += [str((ROOT/p).relative_to(REPO)) for p in ['experiments/vector_run.py','configs/vector_contract_proposed.json','configs/vector_material_design.json','configs/vector_first_queue.json']]
steps=[]
def run(command,label):
    start=time.perf_counter()
    with (OUT/(label+'.log')).open('w') as log:
        p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        for line in p.stdout:log.write(line);log.flush();print(line,end='',flush=True)
        status=p.wait()
    steps.append(dict(operation=label,exit=status,wall_seconds=time.perf_counter()-start,command=command))
    (OUT/'PROGRESS.json').write_text(json.dumps(dict(status='CONNECTED_FINITE_CONTROLLER',steps=steps),indent=2)+'\n')
    if status:raise RuntimeError(label+' failed '+str(status))
try:
    with zipfile.ZipFile(OUT/'SOURCE_BUNDLE.zip','w',zipfile.ZIP_DEFLATED) as z:
        for path in paths:z.write(REPO/path,path)
    (OUT/'BUNDLE_MANIFEST.json').write_text(json.dumps(dict(files=[dict(path=p,bytes=(REPO/p).stat().st_size) for p in paths],new_private_arrays_included=False),indent=2)+'\n')
    run(SSH+['D:/python/python.exe -c "from pathlib import Path; Path(\''+REMOTE+'\').mkdir(exist_ok=False)"'],'create_isolated')
    run(['scp','-i',KEY,str(OUT/'SOURCE_BUNDLE.zip'),HOST+':'+REMOTE+'/SOURCE_BUNDLE.zip'],'upload')
    code="import zipfile,time,json;from pathlib import Path;t=time.perf_counter();c=time.process_time();p=Path('"+REMOTE+"');zipfile.ZipFile(p/'SOURCE_BUNDLE.zip').extractall(p);(p/'EXTRACTION_COST.json').write_text(json.dumps(dict(platform='XINAN',CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t,operation='SOURCE_EXTRACTION_ONLY')))"
    run(SSH+['D:/python/python.exe -c "'+code+'"'],'extract')
    base='D:/python/python.exe -u '+REMOTE+'/research/research_max_struct_r10/experiments/vector_run.py '
    for mode,output,extra in [('probe','RESULTS/vector/queue_probe',''),('generate','RESULTS/vector/generation_initial2',' --ids v02 v09'),('invert','RESULTS/vector/first_cohort_process',' --max-tasks 3')]:
        command=SSH+[base+mode+' --output '+output+extra]
        (OUT/('LAUNCH_'+mode+'.json')).write_text(json.dumps(dict(mode=mode,command=command,started_unix=time.time(),classification='Frozen finite firstcase only; no score'),indent=2)+'\n')
        run(command,mode)
    status='COMPLETE'
except BaseException as exc:
    status='FAILED';print(repr(exc),flush=True)
finally:
    children=resource.getrusage(resource.RUSAGE_CHILDREN)
    (OUT/'TRANSPORT_CONTROLLER_COST.json').write_text(json.dumps(dict(status=status,platform='Mac',local_controller_CPU_seconds=time.process_time()-C,local_transport_child_CPU_seconds=children.ru_utime+children.ru_stime,wall_seconds=time.perf_counter()-T,steps=steps,scope='Own local controller/SSH/SCP only; remote scientific/main/helper CPU excluded',new_private_array_reads=False),indent=2)+'\n')
if status!='COMPLETE':raise SystemExit(1)
