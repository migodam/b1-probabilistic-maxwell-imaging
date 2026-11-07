"""Finite connected R11 exposed-DEV transport, no truth, no R10 writes."""
import time,subprocess,resource,sys,json,zipfile,argparse
from pathlib import Path
T=time.perf_counter();C=time.process_time();ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
p=argparse.ArgumentParser();p.add_argument('--tag',default='dev_ordinary');p.add_argument('--queue',default='configs/vector_dev_ordinary_queue.json');p.add_argument('--count',type=int,default=6);args=p.parse_args()
OUT=ROOT/('delegated/vector/'+args.tag+'_transport');OUT.mkdir(parents=True,exist_ok=False)
REMOTE='D:/AI/MAX_STRUCT_R10_VECTOR_20261107';HOST='jayzh@100.121.97.39';KEY='/Users/migodam/.ssh/heatcryst_xinan_ed25519'
SSH=['ssh','-i',KEY,'-o','BatchMode=yes','-o','ServerAliveInterval=15','-o','ServerAliveCountMax=3',HOST]
paths=[*sorted((ROOT/'src').glob('vector*.py')),ROOT/'experiments/vector_run.py',ROOT/'configs/vector_dev_contract.json',ROOT/args.queue,ROOT/'delegated/vector/source_basis/BASIS.npz',*sorted((ROOT/'delegated/vector/development_observations').rglob('*'))]
paths=[p for p in paths if p.is_file()];steps=[]
def run(command,label):
 start=time.perf_counter()
 with (OUT/(label+'.log')).open('w') as log:
  process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT);status=process.wait()
 steps.append(dict(operation=label,exit=status,wall_seconds=time.perf_counter()-start))
 if status:raise RuntimeError(label+' failed: '+str(status))
 print(json.dumps(steps[-1]),flush=True)
try:
 with zipfile.ZipFile(OUT/'SOURCE.zip','w',zipfile.ZIP_DEFLATED) as z:
  for path in paths:z.write(path,str(path.relative_to(REPO)))
 (OUT/'BUNDLE_MANIFEST.json').write_text(json.dumps(dict(paths=[dict(path=str(p.relative_to(REPO)),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns) for p in paths],GT_access=False),indent=2)+'\n')
 run(['scp','-i',KEY,str(OUT/'SOURCE.zip'),HOST+':'+REMOTE+'/R11_'+args.tag.upper()+'_SOURCE.zip'],'upload')
 code="import zipfile;zipfile.ZipFile('"+REMOTE+"/R11_"+args.tag.upper()+"_SOURCE.zip').extractall('"+REMOTE+"')"
 run(SSH+['D:/python/python.exe -c "'+code+'"'],'extract')
 cmd='D:/python/python.exe -u '+REMOTE+'/research/research_max_struct_r11/experiments/vector_run.py --config configs/vector_dev_contract.json --queue '+args.queue+' --max-tasks '+str(args.count)+' --output delegated/vector/'+args.tag+'_process'
 run(SSH+[cmd],'run')
 task_dirs=['research/research_max_struct_r11/delegated/vector/inversions/'+t['task_id'] for t in json.loads((ROOT/args.queue).read_text())['tasks'][:args.count]]
 bases=['research/research_max_struct_r11/delegated/vector/'+args.tag+'_process']+task_dirs
 code="import zipfile;from pathlib import Path;r=Path('"+REMOTE+"');z=zipfile.ZipFile(r/'R11_"+args.tag.upper()+"_RESULT.zip','w',zipfile.ZIP_DEFLATED);[z.write(p,str(p.relative_to(r))) for base in "+repr(bases)+" for p in (r/base).rglob('*') if p.is_file()];z.close()"
 run(SSH+['D:/python/python.exe -c "'+code+'"'],'pack')
 run(['scp','-i',KEY,HOST+':'+REMOTE+'/R11_'+args.tag.upper()+'_RESULT.zip',str(OUT/'RESULT.zip')],'download')
 with zipfile.ZipFile(OUT/'RESULT.zip') as z:z.extractall(REPO)
 status='COMPLETE'
except BaseException as exc:status='FAILED';print(repr(exc),flush=True)
finally:
 child=resource.getrusage(resource.RUSAGE_CHILDREN)
 (OUT/'TRANSPORT_COST.json').write_text(json.dumps(dict(status=status,platform='Mac',local_controller_CPU_seconds=time.process_time()-C,local_transport_child_CPU_seconds=child.ru_utime+child.ru_stime,wall_seconds=time.perf_counter()-T,steps=steps,scope='Local connected SSH/SCP only; remote scientific CPU excluded'),indent=2)+'\n')
if status!='COMPLETE':raise SystemExit(1)
