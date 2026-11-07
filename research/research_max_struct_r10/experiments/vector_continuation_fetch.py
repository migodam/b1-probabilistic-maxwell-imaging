"""Completed3-arm byte harvest and disjoint accounting, no GT reads."""
import time,resource,subprocess,json,zipfile
from pathlib import Path
T=time.perf_counter();C=time.process_time();usage=resource.getrusage(resource.RUSAGE_CHILDREN)
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];D=ROOT/'delegated/vector/continuation_v04_harvest';D.mkdir(exist_ok=False)
key='/Users/migodam/.ssh/heatcryst_xinan_ed25519';host='jayzh@100.121.97.39';remote='D:/AI/MAX_STRUCT_R10_VECTOR_20261107';steps=[]
def run(args):
 t=time.perf_counter();p=subprocess.run(args,capture_output=True,text=True);steps.append(dict(command=args,exit=p.returncode,stdout=p.stdout,stderr=p.stderr,wall_seconds=time.perf_counter()-t));p.check_returncode()
run(['scp','-i',key,str(ROOT/'experiments/vector_continuation_harvest.py'),host+':'+remote+'/research/research_max_struct_r10/experiments/vector_continuation_harvest.py'])
run(['ssh','-i',key,'-o','BatchMode=yes',host,'D:/python/python.exe '+remote+'/research/research_max_struct_r10/experiments/vector_continuation_harvest.py'])
archive=D/'VECTOR_CONTINUATION_FINAL.zip';run(['scp','-i',key,host+':'+remote+'/VECTOR_CONTINUATION_FINAL.zip',str(archive)])
if (ROOT/'RESULTS/vector_continuation_v04').exists():raise RuntimeError('Refuse overwrite')
with zipfile.ZipFile(archive) as z:
 for info in z.infolist():
  if not (REPO/info.filename).resolve().is_relative_to(REPO.resolve()):raise ValueError('Unsafe path')
 z.extractall(REPO)
out=ROOT/'RESULTS/vector_continuation_v04';checks=[]
for arm in ('sourceCIE','GlobalDual','AplusGlobal'):
 d=json.loads((out/arm/'SEALED.json').read_text())
 for row in d['files']:
  p=REPO/row['path'].replace('\\','/')
  if p.stat().st_size!=row['bytes']:raise RuntimeError('Local bytes mismatch')
 checks.append(arm)
(out/'LOCAL_BYTE_VERIFICATION.json').write_text(json.dumps(dict(status='ALL3_SOURCE_SEALS_AND_LOCAL_BYTES_VERIFIED',arms=checks,GT_read=False,mtime_scope='Windows ns checked at source before zip, local extracted bytes only'),indent=2)+'\n')
u=resource.getrusage(resource.RUSAGE_CHILDREN)
transport=dict(status='COMPLETE',platform='Mac',local_controller_CPU_seconds=time.process_time()-C,local_transport_child_CPU_seconds=u.ru_utime+u.ru_stime-usage.ru_utime-usage.ru_stime,wall_seconds=time.perf_counter()-T,steps=steps,scope='Own continuation byte harvest local SSH/SCP/extraction; remote packaging CPU separate',GT_read=False)
(D/'TRANSPORT_COST.json').write_text(json.dumps(transport,indent=2)+'\n')
# Only continuation incremental process rows. Original18 prefix is nonadditive.
rows=[];views=[]
def add(identity,platform,cpu,wall,p,scope,counts=None):rows.append(dict(identity=identity,platform=platform,CPU_seconds=cpu,wall_seconds=wall,source=str(p.relative_to(ROOT)),scope=scope,counts=counts,additive=True,included_in=None))
p=ROOT/'delegated/vector/windows_v04_continuation_queue/PROCESS_COST.json';q=json.loads(p.read_text());hdr=json.loads((out/'PROCESS_COST.json').read_text());add('scientific3arm-wholeprocess','XINAN',q['CPU_seconds'],q['wall_seconds'],p,q['scope'],hdr['counts'])
p=ROOT/'delegated/vector/windows_v04_continuation_queue/CONTROLLER_COST.json';q=json.loads(p.read_text());add('foreground-controller','XINAN',q['CPU_seconds'],q['wall_seconds'],p,q['scope'])
p=out/'PROCESS_COST.json';add('scientific-query-children','XINAN',hdr['GPU_gate_helper_CPU_seconds'],None,p,'Separate query child CPU; whole scientific main already charged')
for name in ('REMOTE_EXTRACTION_COST.json','HARVEST_COST.json'):
 p=out/name;q=json.loads(p.read_text());add(name,'XINAN',q['CPU_seconds'],q['wall_seconds'],p,'Own byte-only remote helper')
for name in ('continuation_v04_transport','continuation_v04_harvest'):
 p=ROOT/'delegated/vector'/name/'TRANSPORT_COST.json';q=json.loads(p.read_text());add(name+'-main','Mac',q['local_controller_CPU_seconds'],q['wall_seconds'],p,q['scope']);add(name+'-SSH-SCP','Mac',q['local_transport_child_CPU_seconds'],None,p,'Own transport subprocess CPU, remote science excluded')
for arm in checks:
 p=out/arm/'RECEIPT.json';q=json.loads(p.read_text());views.append(dict(arm=arm,source=str(p.relative_to(ROOT)),additive=False,included_in='scientific3arm-wholeprocess',CPU_seconds=q['CPU_seconds'],wall_seconds=q['wall_seconds'],counts=q['counts'],physical_Phi=q['physical_Phi'],PG_inf=q['PG_inf']))
# Compiler work after transport measurement is separately measured, no overlap.
compiler_start=time.process_time();totals={platform:sum(x['CPU_seconds'] for x in rows if x['platform']==platform and x['CPU_seconds'] is not None) for platform in ('Mac','XINAN')}
result=dict(status='FINAL_KNOWN_LOWER_BOUND',finalized_for_parent_compiler=True,classification='POSTHOC_BUDGET_DIAGNOSTIC_EXPOSED_DEVELOPMENT',rows=rows,by_platform_CPU_lower_bounds=totals,nonadditive_views=views,prefix=dict(source='RESULTS/vector/ACCOUNTING_FINAL.json',additive=False,reason='Already paid original18 prefix, not charged again'),GT_read=False,unknowns=[dict(CPU_seconds=None,reason='Initial manual GPU availability and progress queries, implementation editing and final compiler tail not metered; no invented zero')],wall_not_summed_as_makespan=True,GPU_event_scope='Stream intervals contain host/gate/IO gaps, not active GPU time')
(out/'ACCOUNTING_FINAL.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status=result['status'],sealed=len(checks),by_platform=totals)))
