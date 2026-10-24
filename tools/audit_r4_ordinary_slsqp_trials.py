#!/usr/bin/env python3
"""Saved ordinary-control diagnosis only; no GT, mu, Full, optimizer or solve."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,json,csv,collections
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,write,pack,probit_derivative

def main():
 source=OUT/'results/fresh/locked/representatives';job=OUT/'IV/jobs/fresh_locked64_representatives_v1';dest=OUT/'AUDIT/ordinary_slsqp_diagnosis'
 if dest.exists():raise FileExistsError('preserve diagnosis')
 if json.loads((job/'CLOSED.json').read_text())['status']!='COMPLETE':raise ValueError('closedoriginalsourceonly')
 b=Book('II','ordinary_slsqp_saved_trial_diagnosis_v1',max_wall=600,analysis_only=True);status='COMPLETE';error=None
 try:
  dest.mkdir();sigma=json.loads((OUT/'configs/TARGET.json').read_text())['noise']['sigma'];parents=[];trials=[]
  for i in range(64):
   parent=f'R4L1_{i:05d}';case=source/parent/'ordinary_likelihood_shared';result=json.loads((case/'RESULTS.json').read_text());events=[json.loads(s) for s in (case/'TRIAL_EVENTS.jsonl').read_text().splitlines()]
   with np.load(OUT/'results/contexts/locked'/f'{parent}.npz',allow_pickle=False) as z:x0=z['x'].copy();Jq=z['J'].copy();y=z['observation'].copy();F0=z['prediction'].copy()
   J=Jq*probit_derivative(x0)[None,:];r0=pack(F0-y)*np.sqrt(2)/sigma;g=J.T@r0;phi0=float(.5*r0@r0)
   evaluated=[e for e in events if e['event']=='EVALUATED' and e['origin']!='original_qualified_anchor'] if events and 'origin' in events[0] else [e for e in events if e['event']=='EVALUATED' and e['endpoint']!=parent]
   callbacks=[e for e in events if e['event']=='ACCEPTED_ITERATE'];callback_ids={e['endpoint'] for e in callbacks};first=None;values=[]
   for k,e in enumerate(evaluated):
    x=np.asarray(e['x']);delta=x-x0;step=float(np.linalg.norm(delta));normg=float(np.linalg.norm(g));cos=float(-g@delta/(normg*step)) if step and normg else None
    path=job/'endpoints'/f"{e['endpoint']}.npz"
    with np.load(path,allow_pickle=False) as z:prediction=z['prediction'].copy()
    rw=pack(prediction-y)*np.sqrt(2)/sigma;phi=float(.5*rw@rw);linear=J@delta;model=float(.5*np.linalg.norm(r0+linear)**2)
    row=dict(parent=parent,trial=k,endpoint=e['endpoint'],anchor_Phi=phi0,trial_Phi=phi,DeltaPhi=phi-phi0,delta_x_norm=step,gradient_x_norm=normg,
     cosine_to_negative_initial_gradient=cos,initial_directional_derivative=float(g@delta),linear_model_Phi=phi0+float(g@delta),GN_model_Phi=model,
     nonlinear_whitened_secant_error=float(np.linalg.norm(rw-r0-linear)),was_callback_endpoint=e['endpoint'] in callback_ids,
     callback_feasible=e['Phi']<=result['phi_threshold'],trial_prediction_source=str(path),source_events=str(case/'TRIAL_EVENTS.jsonl'))
    trials.append(row);values.append(phi)
    if first is None:first=row
   parents.append(dict(parent=parent,paid_trials=len(evaluated),callbacks=len(callbacks),callbacks_feasible=sum(e['feasible'] for e in callbacks),callbacks_selected=sum(e['selected_best'] for e in callbacks),
    callbacks_with_Phi_below_anchor=sum(e['Phi']<phi0 for e in callbacks),anchor_fallback=result['fallback_to_original_anchor'],optimizer_iterations=result['optimizer']['iterations'],stop_reason=result['stop_reason'],
    first_trial_step_norm=first['delta_x_norm'],first_trial_cosine_to_minus_gradient=first['cosine_to_negative_initial_gradient'],first_trial_Phi=first['trial_Phi'],anchor_Phi=phi0,
    first_trial_delta_equals_minus_gradient_error=float(np.linalg.norm(np.asarray(evaluated[0]['x'])-x0+g)),minimum_paid_trial_Phi=min(values),final_Phi=result['Phi']))
  for name,data in (('PER_PARENT.csv',parents),('ALL_SAVED_ORDINARY_TRIALS.csv',trials)):
   with (dest/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
  def stats(key):
   a=np.array([p[key] for p in parents],float);return dict(minimum=float(a.min()),median=float(np.median(a)),maximum=float(a.max()))
  result=dict(status='SAVED_ONLY_ORDINARY_DIAGNOSIS_COMPLETE',parents=64,paid_trials=len(trials),all_fallback=sum(p['anchor_fallback'] for p in parents),
   zero_callback_parents=sum(p['callbacks']==0 for p in parents),total_callbacks=sum(p['callbacks'] for p in parents),feasible_callbacks=sum(p['callbacks_feasible'] for p in parents),selected_callbacks=sum(p['callbacks_selected'] for p in parents),
   callback_Phi_improvements=sum(p['callbacks_with_Phi_below_anchor'] for p in parents),first_trial_step_norm=stats('first_trial_step_norm'),first_trial_direction_cosine=stats('first_trial_cosine_to_minus_gradient'),
   first_trial_Phi=stats('first_trial_Phi'),first_minus_gradient_error=stats('first_trial_delta_equals_minus_gradient_error'),
   saved_anchor_gradient='Jq*probitPDF then Jx.Tsqrt2/sigma*pack(F0-y),allalready-paidonline',new_RHS=0,new_LU=0,new_states=0,GT_read=False,mu_read=False,old_run_modified=False,
   interpretation='SLSQPcallbackrecordsarenotproofthatobjectiveimproved; budgetstop+fallbackcannotestablishnostrengtheningpossible')
  write(dest/'RESULTS.json',result)
  (dest/'REPORT_ZH.md').write_text('# 原ordinary SLSQP：已付费试探诊断\n\n'+json.dumps(result,ensure_ascii=False,indent=2)+'\n\n只读缓存在线J/原残差、原试探预测与callback事件；不读取GT或mu，不作新求解。原320case与封账不改。callback标签沿用原ACCEPTED_ITERATE，不能把有callback等同于下降或最终可用。\n')
  print(json.dumps(result),flush=True)
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:b.close(status,error)
if __name__=='__main__':main()
