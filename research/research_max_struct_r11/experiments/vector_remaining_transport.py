"""Connected finite R11 remaining DEV+50inverse transport; public inputs only."""
import time,subprocess,resource,sys,json,zipfile
from pathlib import Path
T=time.perf_counter();C=time.process_time();ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
OUT=ROOT/'delegated/vector/remaining_transport';OUT.mkdir(parents=True,exist_ok=False)
REMOTE='D:/AI/MAX_STRUCT_R10_VECTOR_20261107';HOST='jayzh@100.121.97.39';KEY='/Users/migodam/.ssh/heatcryst_xinan_ed25519'
SSH=['ssh','-i',KEY,'-o','BatchMode=yes','-o','ServerAliveInterval=15','-o','ServerAliveCountMax=3',HOST]
configs=['vector_dev_contract.json','vector_dev_denominator_queue.json','vector_confirmation_contract.json','vector_confirmation_execution.json','vector_confirmation_queue.json']
paths=[*sorted((ROOT/'src').glob('vector*.py')),ROOT/'experiments/vector_run.py',ROOT/'experiments/vector_remaining_controller.py',*[ROOT/'configs'/s for s in configs],ROOT/'delegated/vector/source_basis/BASIS.npz',*sorted((ROOT/'delegated/vector/development_observations').rglob('*')),*sorted((ROOT/'delegated/vector/confirmation_observations').rglob('*'))]
paths=[p for p in paths if p.is_file()];steps=[]
def run(command,label,strict=True):
 start=time.perf_counter()
 with (OUT/(label+'.log')).open('w') as log:
  process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT);code=process.wait()
 row=dict(operation=label,exit=code,wall_seconds=time.perf_counter()-start);steps.append(row);print(json.dumps(row),flush=True)
 if code and strict:raise RuntimeError(label+' failed: '+str(code))
 return code
status='FAILED'
try:
 with zipfile.ZipFile(OUT/'SOURCE.zip','w',zipfile.ZIP_DEFLATED) as z:
  for p in paths:z.write(p,str(p.relative_to(REPO)))
 (OUT/'BUNDLE_MANIFEST.json').write_text(json.dumps(dict(files=[dict(path=str(p.relative_to(REPO)),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns) for p in paths],private_truth_included=False,material_design_file_included=False),indent=2)+'\n')
 run(['scp','-i',KEY,str(OUT/'SOURCE.zip'),HOST+':'+REMOTE+'/R11_REMAINING_SOURCE.zip'],'upload')
 code="import zipfile;zipfile.ZipFile('"+REMOTE+"/R11_REMAINING_SOURCE.zip').extractall('"+REMOTE+"')"
 run(SSH+['D:/python/python.exe -c "'+code+'"'],'extract')
 execution=run(SSH+['D:/python/python.exe -u '+REMOTE+'/research/research_max_struct_r11/experiments/vector_remaining_controller.py'],'run',strict=False)
 tasks=[]
 for name in ['vector_dev_denominator_queue.json','vector_confirmation_queue.json']:tasks.extend(json.loads((ROOT/'configs'/name).read_text())['tasks'])
 bases=['research/research_max_struct_r11/delegated/vector/'+s for s in ['remaining_controller','dev_denominator_process','confirmation_process']]+['research/research_max_struct_r11/delegated/vector/inversions/'+t['task_id'] for t in tasks]
 code="import zipfile;from pathlib import Path;r=Path('"+REMOTE+"');z=zipfile.ZipFile(r/'R11_REMAINING_RESULT.zip','w',zipfile.ZIP_DEFLATED);[z.write(p,str(p.relative_to(r))) for base in "+repr(bases)+" for p in (r/base).rglob('*') if p.is_file()];z.close()"
 run(SSH+['D:/python/python.exe -c "'+code+'"'],'pack_all_complete_or_failed_evidence')
 run(['scp','-i',KEY,HOST+':'+REMOTE+'/R11_REMAINING_RESULT.zip',str(OUT/'RESULT.zip')],'download')
 with zipfile.ZipFile(OUT/'RESULT.zip') as z:z.extractall(REPO)
 status='COMPLETE' if execution==0 else 'FAILED_SCIENCE_EVIDENCE_RETRIEVED'
except BaseException as exc:print(repr(exc),flush=True)
finally:
 child=resource.getrusage(resource.RUSAGE_CHILDREN)
 (OUT/'TRANSPORT_COST.json').write_text(json.dumps(dict(status=status,platform='Mac',local_controller_CPU_seconds=time.process_time()-C,local_transport_child_CPU_seconds=child.ru_utime+child.ru_stime,wall_seconds=time.perf_counter()-T,steps=steps,scope='Local connected SSH/SCP only; remote scientific/controller CPU excluded; elapsed wall overlaps remote work'),indent=2)+'\n')
if status!='COMPLETE':raise SystemExit(1)
