"""Byte-only final own remote harvest; never interprets arrays or GT."""
import time,resource,subprocess,json,zipfile
from pathlib import Path
T=time.perf_counter();C=time.process_time();R=resource.getrusage(resource.RUSAGE_CHILDREN)
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];D=ROOT/'delegated/vector/final_harvest_transport';D.mkdir(exist_ok=True)
key='/Users/migodam/.ssh/heatcryst_xinan_ed25519';host='jayzh@100.121.97.39';remote='D:/AI/MAX_STRUCT_R10_VECTOR_20261107';steps=[]
def call(args):
 t=time.perf_counter();p=subprocess.run(args,capture_output=True,text=True);steps.append(dict(command=args,exit=p.returncode,wall_seconds=time.perf_counter()-t,stdout=p.stdout,stderr=p.stderr));p.check_returncode()
call(['ssh','-i',key,'-o','BatchMode=yes',host,'D:/python/python.exe '+remote+'/research/research_max_struct_r10/experiments/vector_status.py'])
archive=D/'VECTOR_FINAL_PUBLIC_AND_ENGINEERING.zip'
call(['scp','-i',key,host+':'+remote+'/VECTOR_FINAL_PUBLIC_AND_ENGINEERING.zip',str(archive)])
for name in ('REMOTE_STATUS_SNAPSHOT.json','STATUS_HELPERS.jsonl'):
 call(['scp','-i',key,host+':'+remote+'/research/research_max_struct_r10/delegated/vector/'+name,str(ROOT/'delegated/vector'/name)])
with zipfile.ZipFile(archive) as z:
 for a in ('sourceCIE','GlobalDual','AplusGlobal'):
  if (ROOT/'RESULTS/vector'/('v08_'+a)).exists():raise RuntimeError('Refuse duplicate primary extraction')
 for info in z.infolist():
  p=(REPO/info.filename).resolve()
  if not p.is_relative_to(REPO.resolve()):raise RuntimeError('Unsafe archive member')
 for info in z.infolist():
  if info.filename.endswith(('REMOTE_STATUS_SNAPSHOT.json','STATUS_HELPERS.jsonl')):continue
  z.extract(info,REPO)
# Verify local bytes for every original primary sealed file (mtime changed by zip extraction).
checks=[]
for a in ('sourceCIE','GlobalDual','AplusGlobal'):
 p=ROOT/'RESULTS/vector'/('v08_'+a);seal=json.loads((p/'SEALED.json').read_text())
 for row in seal['files']:
  target=REPO/row['path'].replace('\\','/')
  if target.stat().st_size!=row['bytes']:raise RuntimeError('Local byte size mismatch '+str(target))
 checks.append(p.name)
(ROOT/'delegated/vector/PRIMARY_LOCAL_BYTE_VERIFICATION.json').write_text(json.dumps(dict(status='ALL_REMOTE_SOURCE_SEALS_VERIFIED_AND_LOCAL_BYTES_MATCH',tasks=checks,GT_arrays_read=False,mtime_scope='Original Windows seal verified before zip; local extracted mtimes not used'),indent=2)+'\n')
R2=resource.getrusage(resource.RUSAGE_CHILDREN)
receipt=dict(status='COMPLETE',platform='Mac',local_controller_CPU_seconds=time.process_time()-C,local_transport_child_CPU_seconds=R2.ru_utime+R2.ru_stime-R.ru_utime-R.ru_stime,wall_seconds=time.perf_counter()-T,steps=steps,scope='Own local SSH/SCP plus byte extraction and local seal byte check; remote packaging/helper CPU separate',GT_arrays_read=False)
(D/'TRANSPORT_COST.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:v for k,v in receipt.items() if k!='steps'}))
