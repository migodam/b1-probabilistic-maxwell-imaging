"""Metered public-only continuation deployment and foreground execution."""
import time,resource,subprocess,json,zipfile
from pathlib import Path
T=time.perf_counter();C=time.process_time();usage=resource.getrusage(resource.RUSAGE_CHILDREN)
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];D=ROOT/'delegated/vector/continuation_v04_transport';D.mkdir(exist_ok=False)
import numpy as np
initials=ROOT/'data_public/vector/continuation_v04';initials.mkdir(exist_ok=False)
for arm in ('sourceCIE','GlobalDual','AplusGlobal'):
 p=ROOT/'RESULTS/vector'/('v04_'+arm)/'physical_accepted80.npz'
 with np.load(p,allow_pickle=False) as z:material=z['material'].copy()
 np.savez_compressed(initials/(arm+'.npz'),material=material)
files=list(initials.glob('*.npz'))+[ROOT/'configs/vector_continuation_v04.json',ROOT/'experiments/vector_continuation.py',ROOT/'experiments/vector_windows_continuation_queue.py']
archive=D/'PUBLIC_CONTINUATION.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_STORED) as z:
 for p in files:z.write(p,p.relative_to(REPO).as_posix())
key='/Users/migodam/.ssh/heatcryst_xinan_ed25519';host='jayzh@100.121.97.39';remote='D:/AI/MAX_STRUCT_R10_VECTOR_20261107';steps=[]
def run(args):
 t=time.perf_counter();p=subprocess.run(args,stdout=open(D/'RUN.log','a'),stderr=subprocess.STDOUT);steps.append(dict(command=args,exit=p.returncode,wall_seconds=time.perf_counter()-t));p.check_returncode()
try:
 run(['scp','-i',key,str(archive),host+':'+remote+'/PUBLIC_CONTINUATION.zip'])
 code="import zipfile,time,json;from pathlib import Path;c=time.process_time();t=time.perf_counter();r=Path('"+remote+"');zipfile.ZipFile(r/'PUBLIC_CONTINUATION.zip').extractall(r);(r/'CONTINUATION_EXTRACTION_COST.json').write_text(json.dumps(dict(platform='XINAN',CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t)))"
 run(['ssh','-i',key,'-o','BatchMode=yes',host,'D:/python/python.exe -c "'+code+'"'])
 run(['ssh','-i',key,'-o','BatchMode=yes','-o','ServerAliveInterval=15','-o','ServerAliveCountMax=3',host,'D:/python/python.exe -u '+remote+'/research/research_max_struct_r10/experiments/vector_windows_continuation_queue.py']);status='COMPLETE'
except BaseException as exc:status='FAILED';raise
finally:
 u=resource.getrusage(resource.RUSAGE_CHILDREN);d=dict(status=status,platform='Mac',local_controller_CPU_seconds=time.process_time()-C,local_transport_child_CPU_seconds=u.ru_utime+u.ru_stime-usage.ru_utime-usage.ru_stime,wall_seconds=time.perf_counter()-T,steps=steps,scope='Own public initial extraction/packaging and foreground SSH/SCP; remote science separate',GT_read=False);(D/'TRANSPORT_COST.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({k:v for k,v in d.items() if k!='steps'}),flush=True)
