"""Append48 public CIE→Global choices; original96 seal/rule remain unchanged."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import argparse,time,shutil,json
from pathlib import Path
from datetime import datetime,timezone
from twod_secondary_seal import all_terminal,read,dump,resolve,attempt,add_cost
from twod_recovery import deployment_cost_info
ROOT=Path(__file__).resolve().parents[1]
def utc():return datetime.now(timezone.utc).isoformat()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--manifest',type=Path,required=True);ap.add_argument('--queue',type=Path,action='append',required=True);ap.add_argument('--original96',type=Path,default=ROOT/'delegated/twod/independent48/secondary/public_seal');ap.add_argument('--rule',type=Path,default=ROOT/'configs/twod_independent48_secondary_global_appendix.json');ap.add_argument('--output',type=Path,default=ROOT/'delegated/twod/independent48/secondary/global_appendix/public_seal');a=ap.parse_args();wall=time.perf_counter();rule=read(a.rule);old=read(a.original96/'RECEIPT.json')
 if rule['status']!='FROZEN_APPEND_ONLY_GLOBAL_POLICY_BEFORE48_GT' or rule['fallback']!='global_dual' or old['status']!='ALL96_PUBLIC_DEPLOYMENTS_SEALED_BEFORE48_GT':raise ValueError('Frozen appendix and unchanged original96 public seal required')
 by_case,jobs,receipts=all_terminal(a.manifest,a.queue);triggers=read(a.original96/'TRIGGER_DECISIONS.json')['rows'];oldchoices=read(a.original96/'DEPLOYMENT_SELECTIONS.json')['rows']
 if len(triggers)!=48 or len(oldchoices)!=96 or {g['case_id'] for g in triggers}!=set(by_case):raise ValueError('Unchanged48 CIE trigger inventory required')
 a.output.mkdir(parents=True,exist_ok=False);rows=[];started=utc()
 for gate in triggers:
  case=gate['case_id'];ts=by_case[case];pt,ft=ts['cie_single'],ts['global_dual'];prefix=attempt(pt,jobs[pt['task_id']],receipts[pt['task_id']]);fallback=attempt(ft,jobs[ft['task_id']],receipts[ft['task_id']]);previous=next(r for r in oldchoices if r['case_id']==case)
  if prefix['Phi']!=previous['attempts'][0]['Phi'] or prefix['task_id']!=previous['attempts'][0]['task_id']:raise ValueError('Original CIE prefix changed')
  attempts=[dict(prefix,kind='prefix')]
  if gate['trigger']:attempts.append(dict(fallback,kind='fallback'))
  eligible=[x for x in attempts if x['valid'] and x['Phi'] is not None];chosen=min(eligible,key=lambda x:x['Phi']) if eligible else None;directory=a.output/'deployments'/(case+'__CIE_then_GlobalDual');directory.mkdir(parents=True,exist_ok=False)
  row=dict(case_id=case,method='CIE_then_GlobalDual',trigger=gate['trigger'],trigger_source=str(a.original96/'TRIGGER_DECISIONS.json'),selected=chosen,selection='Same originalPhi minimum; exact ties CIE; no GT',attempts=attempts,prefix_full_child_CPU_seconds=prefix['full_child_CPU_seconds'],fallback_full_child_CPU_seconds=fallback['full_child_CPU_seconds'] if gate['trigger'] else 0,precomputed_fallback_trial_full_child_CPU_seconds=fallback['full_child_CPU_seconds'],deployment_full_child_CPU_seconds=add_cost(prefix,fallback,'full_child_CPU_seconds',gate['trigger']),deployment_physical_LU=add_cost(prefix,fallback,'physical_LU',gate['trigger']),deployment_physical_complex_RHS=add_cost(prefix,fallback,'physical_complex_RHS',gate['trigger']),deployment_cost_additive=False,shared_campaign_new_inverse_CPU_seconds=0,truth_access='NONE',registration_context=rule['registration_context'],budget_scope='Triggered view can pay two480CPU allowances; not equal actual CPU or budget to one CIE attempt')
  row.update(deployment_cost_info(prefix,fallback,gate['trigger']))
  if chosen:
   source=resolve(chosen['path']);shutil.copy2(source/'R11_CAP_PREDICTION.npz',directory/'prediction.npz');shutil.copy2(source/'R11_CAP_PREDICTION.json',directory/'prediction.json')
  dump(directory/'SELECTION.json',row);rows.append(row)
 dump(a.output/'DEPLOYMENT_SELECTIONS.json',dict(status='ADDITIONAL48_GLOBAL_PUBLIC_CHOICES_SEALED',rows=rows,truth_access='NONE'))
 dump(a.output/'TOTAL144_SELECTION_REFERENCES.json',dict(original96=str(a.original96/'DEPLOYMENT_SELECTIONS.json'),appended48=str(a.output/'DEPLOYMENT_SELECTIONS.json'),total_choices=144,original96_unchanged=True,truth_access='NONE'))
 dump(a.output/'FROZEN_APPENDIX_USED.json',rule)
 receipt=dict(status='ADDITIONAL48_GLOBAL_PUBLIC_CHOICES_SEALED_BEFORE48_GT',registered_utc=rule['registered_utc'],registration_context=rule['registration_context'],seal_started_utc=started,seal_finished_utc=utc(),additional_choices=48,total_deployment_choices=144,original96_unchanged=True,triggers=sum(r['trigger'] for r in rows),valid_selected=sum(r['selected'] is not None for r in rows),new_inverse_runs=0,new_forward_calls=0,truth_access='NONE',CPU_seconds=time.process_time(),CPU_scope='Absolute whole append-only public selector process_time, imports/copies included through receipt preparation',wall_seconds=time.perf_counter()-wall)
 dump(a.output/'RECEIPT.json',receipt);print(json.dumps({k:receipt[k] for k in ('status','additional_choices','total_deployment_choices','CPU_seconds')}))
if __name__=='__main__':main()
