"""Predeclared secondary full-operator filters, only when time permits.

Selection uses public manifest order and parent indices, never TEST truth or
scored outcomes. No weak-mode truncation is introduced.
"""
from pathlib import Path
import datetime,json,os,subprocess,sys,time
S=Path(__file__).resolve().parents[1]
ENV=os.environ.copy()
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):ENV[k]='1'
records=json.loads((S/'data/manifests/online.json').read_text())['records']
ids=[];seen=set()
for r in records:
    if r['split']!='test':continue
    i=int(r['parent_id'].split('p')[-1]);key=(r['parent_id'],r['noise'])
    # First pair per selected parent is the low-amplitude cohort by the frozen
    # generation ordering; exclude subsequent stress pair without reading GT.
    if (i//4)%4!=0 or key in seen:continue
    seen.add(key);ids.append(r['obs_id'])
plans=[dict(arm=arm,source=f'confirm1024_{arm}_seed7',run=f'confirm1024_{arm}_seed7_filter_{kind}',mode=mode)for arm in ('v_physics','j_null')for kind,mode in [('background','fixed_background_softfilter_box_qp'),('current','full_softfilter_box_qp')]]
frozen=dict(frozen_unix=time.time(),plans=plans,obs_ids=ids,selection='16 public parent indices, (index//4)%4=0; first clean/noisy pair per parent, 8ID+8OOD and all4families',time_gate='Start only before 2026-11-07T02:40:00Z; otherwise NOT_RUN. No outcome gate.',truth_access='NONE',scope='Secondary full-operator filtering, not model selection or weak truncation')
(S/'audit/SECONDARY_PLAN_FROZEN.json').write_text(json.dumps(frozen,indent=2))
while not (S/'audit/CONFIRMATION_PHYSICS_DONE.json').exists():time.sleep(30)
if time.time()>datetime.datetime(2026,11,7,2,40,tzinfo=datetime.timezone.utc).timestamp():
    (S/'audit/SECONDARY_DONE.json').write_text(json.dumps(dict(status='NOT_RUN_TIME_CAP',unix=time.time())));sys.exit(0)
start=time.time()
for plan in plans:
    dest=S/'predictions'/plan['run'];dest.mkdir(exist_ok=True)
    source=S/'predictions'/plan['source']
    (dest/'REUSED_PROPOSALS.json').write_text(json.dumps(dict(source_run=plan['source'],actual_new_NN_inference=False,scope='Same frozen public proposals; per-arm deployment preparation is conservatively re-invoiced by evaluator')))
    for obs in ids:
        for suffix in ('.npz','.json'):
            os.link(source/f'{obs}{suffix}',dest/f'{obs}{suffix}')
    os.link(source/'inference_setup.json',dest/'inference_setup.json')
    children=[];files=[]
    for index in (0,1):
        f=(S/f'evaluation/{plan["run"]}_P{index}.log').open('w');files.append(f)
        cmd=[sys.executable,str(S/'experiments/run_evaluation.py'),'evaluate','--runs',plan['run'],'--splits','test','--obs-id',*ids,'--secondary-mode',plan['mode'],'--partition-count','2','--partition-index',str(index),'--partition-block-size','8']
        children.append(subprocess.Popen(cmd,cwd=S,env=ENV,stdout=f,stderr=subprocess.STDOUT))
    codes=[p.wait()for p in children]
    for f in files:f.close()
    plan['returncodes']=codes
(S/'audit/SECONDARY_DONE.json').write_text(json.dumps(dict(status='COMPLETED',wall_seconds=time.time()-start,plans=plans,obs_ids=ids,truth_access='NONE'),indent=2))
