"""All288-terminal public-only replay seal of the frozen R10 conditional policy."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import argparse,json,time,shutil
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from twod_recovery import active_task,active_job,queue_receipt,attempt_cost_info,deployment_cost_info
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
METHODS=('cie_single','cie_dual','global_dual','fullG19','freezeA','contractedT')
def utc():return datetime.now(timezone.utc).isoformat()
def read(p):return json.loads(Path(p).read_text())
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def resolve(p):
 p=Path(str(p).replace('\\','/'));return p if p.is_absolute() else REPO/p
def method(t):return t.get('material_readout',{}).get('label',t['arm'])
def all_terminal(manifest,queue_paths):
 tasks=[active_task(t) for t in read(manifest)['tasks']];cases=sorted({t['case_id'] for t in tasks});lookup={t['task_id']:t for t in tasks};jobs={}
 if len(tasks)!=288 or len(lookup)!=288 or len(cases)!=48:raise ValueError('Original all288 manifest required')
 for p in queue_paths:
  p=queue_receipt(p);q=read(p)
  if q['status'] not in ('COMPLETED','COMPLETED_WITH_ERRORS','STOPPED_DISPATCH_CHILDREN_FINISHED','INTERRUPTED_CONTROLLER_RETAINED_TERMINAL_PREFIX') or q.get('running'):raise RuntimeError('All dispatchers must be terminal before secondary sealing')
  if q['status'] not in ('STOPPED_DISPATCH_CHILDREN_FINISHED','INTERRUPTED_CONTROLLER_RETAINED_TERMINAL_PREFIX') and (q.get('pending') or q.get('not_started_task_ids')):raise RuntimeError('Pending queue')
  for j in q['jobs']:
   if not active_job(j):continue
   if j['task_id'] in jobs or j['status'] not in ('FINISHED','FAILED','LAUNCH_FAILED'):raise ValueError('Duplicate/nonterminal job')
   jobs[j['task_id']]=dict(j,source_queue_receipt=str(p))
 if set(jobs)!=set(lookup):raise RuntimeError('Exactly288 unique terminal paid outcomes required')
 by_case={}
 for case in cases:
  ts=[t for t in tasks if t['case_id']==case]
  if sorted(method(t) for t in ts)!=sorted(METHODS):raise ValueError('Six balanced methods required')
  by_case[case]={method(t):t for t in ts}
 receipts={}
 for t in tasks:
  p=resolve(t['output_dir'])/'RECEIPT.json';r=read(p) if p.exists() else {}
  if r.get('status') in ('RUNNING','PREPARING','PREPARED'):raise RuntimeError('Nonterminal child receipt')
  j=jobs[t['task_id']]
  if j['status']=='FINISHED':
   if not p.exists() or r.get('status') not in ('SEALED_PUBLIC_PREDICTION','CPU_BUDGET_REACHED'):raise RuntimeError('Finished worker public receipt must be synced before policy seal')
   if r.get('cap_prediction_status')=='SEALED_PAID_INCUMBENT':
    if (r.get('task_id'),r.get('case_id'),r.get('arm'))!=(t['task_id'],t['case_id'],t['arm']):raise ValueError('Paid attempt identity mismatch')
    for name in ('R11_CAP_PREDICTION.npz','R11_CAP_PREDICTION.json'):
     if not (p.parent/name).exists():raise RuntimeError('Missing synchronized cap prediction before policy seal')
  receipts[t['task_id']]=r
 return by_case,jobs,receipts

def public_noise(task):
 data={}
 for filename in task['public_files']:
  path=resolve(filename)
  if {'data_private','private','truth','gt','sealed_gt','sealed_truth'} & {x.lower() for x in path.parts}:raise ValueError('Private input rejected')
  with np.load(path,allow_pickle=False) as z:
   if {'truth','gt','m_true','m_truth','y_clean','yclean','family','seed','alpha','mask'} & {k.lower() for k in z.files}:raise ValueError('Private payload rejected')
   for key in ('Y','sigma','noise_radius'):
    if key in z.files:
     v=z[key].copy()
     if key in data and not np.array_equal(data[key],v):raise ValueError('Conflicting public payload')
     data[key]=v
 if set(data)!={'Y','sigma','noise_radius'}:raise ValueError('R10 rule requires stored public Y/sigma/noise_radius')
 sigma=float(data['sigma']);radius=float(data['noise_radius'])
 if not np.isfinite(sigma) or sigma<=0 or not np.isfinite(radius) or radius<0:raise ValueError('Invalid public R10 noise contract')
 return np.asarray(data['Y'],complex),sigma,radius

def attempt(t,j,r):
 d=resolve(t['output_dir']);valid=j['status']=='FINISHED' and r.get('status') in ('SEALED_PUBLIC_PREDICTION','CPU_BUDGET_REACHED') and r.get('cap_prediction_status')=='SEALED_PAID_INCUMBENT' and (d/'R11_CAP_PREDICTION.npz').is_file()
 phi=r.get('cap_physical_objective') if valid else None
 if phi is not None and not np.isfinite(phi):valid=False;phi=None
 if valid:
  p=read(d/'R11_CAP_PREDICTION.json')
  if float(p['objective'])!=float(phi):raise ValueError('Cap prediction originalPhi mismatch')
 c=r.get('physical_model_cost') or j.get('physical_model_cost') or {}
 result=dict(task_id=t['task_id'],method=method(t),source_queue_receipt=j.get('source_queue_receipt'),cpu_logical=j.get('cpu'),path=t['output_dir'],valid=valid,Phi=phi,status=r.get('status',j['status']),full_child_CPU_seconds=j.get('child_process_CPU_seconds'),physical_LU=c.get('LU'),physical_complex_RHS=c.get('actual_complex_RHS'));result.update(attempt_cost_info(t['task_id'],j.get('child_process_CPU_seconds')));return result

def add_cost(prefix,fallback,key,trigger):
 a=prefix.get(key);b=fallback.get(key) if trigger else 0
 return a+b if a is not None and b is not None else None

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--manifest',type=Path,required=True);ap.add_argument('--queue',type=Path,action='append',required=True);ap.add_argument('--policy',type=Path,default=ROOT/'configs/twod_independent48_secondary_policy.json');ap.add_argument('--output',type=Path,default=ROOT/'delegated/twod/independent48/secondary/public_seal');a=ap.parse_args();wall=time.perf_counter();policy=read(a.policy)
 if policy['status']!='FROZEN_SECONDARY_RULE_REGISTERED_BEFORE_ANY48_GT_ACCESS':raise ValueError('Frozen secondary rule required')
 if policy['prefix']!='cie_single' or policy['fallbacks']!={'CIE_then_G19_PGP':'fullG19','CIE_then_DualCIE':'cie_dual'}:raise ValueError('Changed R10 policy')
 by_case,jobs,receipts=all_terminal(a.manifest,a.queue)
 a.output.mkdir(parents=True,exist_ok=False);started=utc();triggers=[];deployments=[]
 dump(a.output/'PUBLIC_GATE_BEFORE_SELECTION.json',dict(status='PASS_ALL288_TERMINAL_NO_GT',created_utc=started,root_registration_context=policy['registration_context'],manifest=str(a.manifest),queues=[str(p) for p in a.queue],truth_access='NONE'))
 for case,ts in sorted(by_case.items()):
  pt=ts['cie_single'];prefix=attempt(pt,jobs[pt['task_id']],receipts[pt['task_id']]);gate=dict(case_id=case,prefix_task_id=pt['task_id'],trigger=False,status='INVALID_PREFIX_NO_COMPUTABLE_TRIGGER',strict_rule='||paidY_prediction-publicY||F / public sigma > stored public noise_radius')
  if prefix['valid']:
   Y,sigma,radius=public_noise(pt)
   with np.load(resolve(pt['output_dir'])/'R11_CAP_PREDICTION.npz',allow_pickle=False) as z:prediction=z['Y_prediction']
   if prediction.shape!=Y.shape or not np.isfinite(prediction).all():raise ValueError('Invalid paid prefix observation prediction')
   statistic=float(np.linalg.norm(prediction-Y)/sigma);gate.update(status='EVALUATED',whitened_port_residual=statistic,CN99_radius=radius,public_sigma=sigma,trigger=bool(statistic>radius))
  triggers.append(gate)
  for deployment,fallback_method in policy['fallbacks'].items():
   ft=ts[fallback_method];fallback=attempt(ft,jobs[ft['task_id']],receipts[ft['task_id']]);attempts=[dict(prefix,kind='prefix')]
   if gate['trigger']:attempts.append(dict(fallback,kind='fallback'))
   eligible=[x for x in attempts if x['valid'] and x['Phi'] is not None];chosen=min(eligible,key=lambda x:x['Phi']) if eligible else None
   directory=a.output/'deployments'/(case+'__'+deployment);directory.mkdir(parents=True,exist_ok=False)
   row=dict(case_id=case,method=deployment,trigger=gate['trigger'],prefix_status=prefix['status'],fallback_status=fallback['status'] if gate['trigger'] else 'NOT_TRIGGERED',selected=chosen,selection='Minimum same originalPhi; exact ties CIE; no GT',attempts=attempts,prefix_full_child_CPU_seconds=prefix['full_child_CPU_seconds'],fallback_full_child_CPU_seconds=fallback['full_child_CPU_seconds'] if gate['trigger'] else 0,precomputed_fallback_trial_full_child_CPU_seconds=fallback['full_child_CPU_seconds'],deployment_full_child_CPU_seconds=add_cost(prefix,fallback,'full_child_CPU_seconds',gate['trigger']),deployment_physical_LU=add_cost(prefix,fallback,'physical_LU',gate['trigger']),deployment_physical_complex_RHS=add_cost(prefix,fallback,'physical_complex_RHS',gate['trigger']),deployment_cost_additive=False,shared_campaign_new_inverse_CPU_seconds=0,truth_access='NONE',execution_scope=policy['execution_scope'],budget_scope=policy['budget_scope'])
   row.update(deployment_cost_info(prefix,fallback,gate['trigger']))
   if chosen:
    source=resolve(chosen['path']);shutil.copy2(source/'R11_CAP_PREDICTION.npz',directory/'prediction.npz');shutil.copy2(source/'R11_CAP_PREDICTION.json',directory/'prediction.json')
   dump(directory/'SELECTION.json',row);deployments.append(row)
 dump(a.output/'TRIGGER_DECISIONS.json',dict(status='PUBLIC_ONLY48_DECISIONS',rows=triggers,new_forward_calls=0,truth_access='NONE'))
 dump(a.output/'DEPLOYMENT_SELECTIONS.json',dict(status='ALL96_PUBLIC_DEPLOYMENTS_SEALED',rows=deployments,new_forward_calls=0,truth_access='NONE'))
 dump(a.output/'FROZEN_POLICY_USED.json',policy)
 summary=[]
 for name in policy['deployments']:
  selected=[r for r in deployments if r['method']==name];cpu=[r['deployment_full_child_CPU_seconds'] for r in selected]
  summary.append(dict(method=name,deployments=len(selected),triggered=sum(r['trigger'] for r in selected),valid_selected=sum(r['selected'] is not None for r in selected),deployment_scientific_prefix_CPU_seconds=sum(cpu) if all(x is not None for x in cpu) else None,exact_end_to_end_cost_qualified=all(r.get('exact_end_to_end_cost_qualified',True) for r in selected),deployment_end_to_end_CPU_exact=None if any(not r.get('exact_end_to_end_cost_qualified',True) for r in selected) else sum(cpu),not_additive_with_other_view=True))
 receipt=dict(status='ALL96_PUBLIC_DEPLOYMENTS_SEALED_BEFORE48_GT',version=policy['version'],registered_utc=policy['registered_utc'],registration_context=policy['registration_context'],seal_started_utc=started,seal_finished_utc=utc(),primary_tasks_terminal=288,objects=48,deployments=96,triggers=sum(g['trigger'] for g in triggers),summary=summary,truth_access='NONE',new_inverse_runs=0,new_forward_calls=0,new_adjoint_calls=0,new_physical_LU=0,new_physical_RHS=0,primary_campaign_actual_child_CPU_seconds=sum(j['child_process_CPU_seconds'] for j in jobs.values()) if all(j.get('child_process_CPU_seconds') is not None for j in jobs.values()) else None,unknown_CPU_jobs=[j['task_id'] for j in jobs.values() if j.get('child_process_CPU_seconds') is None],campaign_cost_policy='All288 primary attempts remain charged once. Deployment views replay CIE plus only triggered corresponding fallback; never add alternative views as campaign cost.',secondary_seal_CPU_seconds=time.process_time(),CPU_scope='Absolute whole secondary seal process_time includes imports/public reads/copy/selection through receipt preparation',wall_seconds=time.perf_counter()-wall)
 dump(a.output/'RECEIPT.json',receipt);print(json.dumps({k:receipt[k] for k in ('status','objects','deployments','triggers','secondary_seal_CPU_seconds')}))
if __name__=='__main__':main()
