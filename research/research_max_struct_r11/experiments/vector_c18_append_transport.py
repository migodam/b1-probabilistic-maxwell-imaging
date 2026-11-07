"""Public-only c18append transport and canonical sealed collection, no GT."""
import time,json,subprocess,resource,zipfile,shutil,sys
from pathlib import Path
ENTRY=time.perf_counter();ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
OUT=ROOT/'delegated/vector/c18_append_transport';OUT.mkdir(parents=True,exist_ok=False)
REMOTE='D:/AI/MAX_STRUCT_R10_VECTOR_20261107';HOST='jayzh@100.121.97.39';KEY='/Users/migodam/.ssh/heatcryst_xinan_ed25519';SSH=['ssh','-i',KEY,'-o','BatchMode=yes','-o','ServerAliveInterval=15','-o','ServerAliveCountMax=3',HOST]
paths=[*sorted((ROOT/'src').glob('vector*.py')),ROOT/'experiments/vector_run.py',ROOT/'experiments/vector_gpu_c18_append_controller.py',ROOT/'configs/vector_confirmation_c18_gpu_execution.json',ROOT/'configs/vector_confirmation_c18_gpu_append_queue.json',ROOT/'configs/vector_confirmation_host_dispatch.json'];steps=[];status='FAILED'
def run(command,name,strict=True):
 start=time.perf_counter()
 with (OUT/(name+'.log')).open('w') as f:p=subprocess.Popen(command,stdout=f,stderr=subprocess.STDOUT);code=p.wait()
 steps.append(dict(operation=name,exit=code,wall_seconds=time.perf_counter()-start));print(json.dumps(steps[-1]),flush=True)
 if code and strict:raise RuntimeError(name+' failed')
 return code
try:
 with zipfile.ZipFile(OUT/'SOURCE.zip','w',zipfile.ZIP_DEFLATED) as z:
  for p in paths:z.write(p,str(p.relative_to(REPO)))
 (OUT/'BUNDLE_MANIFEST.json').write_text(json.dumps(dict(paths=[dict(path=str(p.relative_to(REPO)),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns) for p in paths],new_GT_or_observation_generation=False,scientific_tasks=10,already_frozen_original50_subset=True),indent=2)+'\n')
 run(['scp','-i',KEY,str(OUT/'SOURCE.zip'),HOST+':'+REMOTE+'/R11_C18_APPEND_SOURCE.zip'],'upload')
 code="import zipfile;zipfile.ZipFile('"+REMOTE+"/R11_C18_APPEND_SOURCE.zip').extractall('"+REMOTE+"')";run(SSH+['D:/python/python.exe -c "'+code+'"'],'extract')
 exitcode=run(SSH+['D:/python/python.exe -u '+REMOTE+'/research/research_max_struct_r11/experiments/vector_gpu_c18_append_controller.py'],'run',False)
 bases=['research/research_max_struct_r11/delegated/vector/'+name for name in ('c18_gpu_append_controller','c18_gpu_append_process','c18_gpu_append_inversions')]
 code="from pathlib import Path;import zipfile,time,json;r=Path('"+REMOTE+"');t=time.perf_counter();z=zipfile.ZipFile(r/'R11_C18_APPEND_RESULT.zip','w');files=[p for base in "+repr(bases)+" for p in (r/base).rglob('*') if p.is_file()];[z.write(p,str(p.relative_to(r)),compress_type=zipfile.ZIP_STORED if p.suffix=='.npz' else zipfile.ZIP_DEFLATED) for p in files];z.close();print(json.dumps(dict(CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-t,files=len(files),scope='Whole remote pack only,no physics')))"
 run(SSH+['D:/python/python.exe -c "'+code+'"'],'pack_complete_or_failed_evidence')
 run(['scp','-i',KEY,HOST+':'+REMOTE+'/R11_C18_APPEND_RESULT.zip',str(OUT/'RESULT.zip')],'download')
 with zipfile.ZipFile(OUT/'RESULT.zip') as z:z.extractall(REPO)
 mappings=[]
 if exitcode==0:
  for task in json.loads((ROOT/'configs/vector_confirmation_c18_gpu_append_queue.json').read_text())['tasks']:
   src=ROOT/'delegated/vector/c18_gpu_append_inversions'/task['task_id'];dst=ROOT/'delegated/vector/inversions'/task['task_id']
   if dst.exists():raise RuntimeError('Refuse duplicate canonical task'+task['task_id'])
   if json.loads((src/'SEALED.json').read_text())['status']!='PREDICTION_SEALED_BEFORE_ROOT_OFFLINE_SCORING':raise RuntimeError('Append not sealed')
   shutil.move(src,dst);mappings.append(dict(task_id=task['task_id'],original_archive_path=str(src.relative_to(REPO)),canonical_collection_path=str(dst.relative_to(REPO)),data_values_unchanged=True))
  (OUT/'CANONICAL_COLLECTION_MAP.json').write_text(json.dumps(dict(status='TEN_UNIQUE_GPU_APPEND_SEALS_CANONICALLY_COLLECTED',mappings=mappings,GT_access=False),indent=2)+'\n')
 status='COMPLETE' if exitcode==0 else 'FAILED_SCIENCE_EVIDENCE_RETRIEVED'
except BaseException as exc:print(repr(exc),flush=True)
finally:
 children=resource.getrusage(resource.RUSAGE_CHILDREN);(OUT/'TRANSPORT_COST.json').write_text(json.dumps(dict(status=status,local_controller_CPU_seconds=time.process_time(),local_transport_child_CPU_seconds=children.ru_utime+children.ru_stime,wall_seconds=time.perf_counter()-ENTRY,steps=steps,scope='Own local transport/controller CPU only;remote scientific/controller costs separate,elapsedwall overlaps waiting25+append10',GT_access=False),indent=2)+'\n')
raise SystemExit(0 if status=='COMPLETE' else 1)
