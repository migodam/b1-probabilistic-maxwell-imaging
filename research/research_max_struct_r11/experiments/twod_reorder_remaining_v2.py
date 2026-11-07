"""Engineeringv2: aggregate prior terminal dispatchers, expand only unstarted tasks."""
import argparse,json,time
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--queue',type=Path,action='append',required=True);ap.add_argument('--manifest',type=Path,required=True);ap.add_argument('--cpus',type=int,nargs='+',required=True);a=ap.parse_args()
 if len(set(a.cpus))!=len(a.cpus) or not set(a.cpus).issubset({0,2,4,6,8,10}) or len(a.cpus)>4:raise ValueError('Root authorized independent Pcores only')
 doc=json.loads(a.manifest.read_text());tasks=doc['tasks'];lookup={t['task_id']:t for t in tasks};done={};qs=[]
 for p in a.queue:
  q=json.loads(p.read_text());qs.append(q)
  if q['status'] not in ('COMPLETED','COMPLETED_WITH_ERRORS','STOPPED_DISPATCH_CHILDREN_FINISHED') or q.get('running'):raise RuntimeError('Prior own queues must finish inflight naturally')
  for j in q['jobs']:
   if j['task_id'] in done:raise ValueError('Do not duplicate paid task')
   done[j['task_id']]=j
 if qs[-1]['status']!='STOPPED_DISPATCH_CHILDREN_FINISHED':raise ValueError('Last own dispatcher must explicitly stop future launch')
 unstarted=set(qs[-1]['not_started_task_ids'])
 if unstarted!=set(lookup)-set(done):raise ValueError('Aggregate paid/unstarted partition mismatch')
 methods=('cie_single','cie_dual','global_dual','fullG19','freezeA','contractedT');ordered=[]
 for case in sorted({t['case_id'] for t in tasks}):
  rotation=int(case.rsplit('_',1)[-1])%6
  for method in methods[rotation:]+methods[:rotation]:
   found=[t for t in tasks if t['case_id']==case and t.get('material_readout',{}).get('label',t['arm'])==method]
   if len(found)!=1:raise ValueError('Balanced frozen case/method required')
   t=found[0]
   if t['task_id'] in unstarted:
    if (REPO/t['output_dir']).exists():raise FileExistsError('Unstarted scientific output exists; preserve and audit')
    ordered.append(t)
 if len(ordered)!=len(unstarted):raise ValueError('Task dropped or duplicated')
 target=ROOT/'configs/twod_independent48_expanded_remaining_tasks.json';receiptpath=ROOT/'delegated/twod/independent48/EXPANDED_DISPATCH_RECEIPT.json'
 for p in (target,receiptpath):
  if p.exists():raise FileExistsError(p)
 target.write_text(json.dumps(dict(status='ENGINEERING_V2_EXPANDED_INDEPENDENT_PCORES_NO_SCIENCE_CHANGE',tasks=ordered),indent=2)+'\n')
 receipt=dict(status='REORDERED_ONLY_UNSTARTED_FOR_EXPANDED_PCORE_QUEUE',engineering_version=2,created_utc=datetime.now(timezone.utc).isoformat(),cpus=a.cpus,threads_per_worker=1,whole_child_CPU_cap_seconds=480,prior_queues=[str(p) for p in a.queue],retained_paid_task_ids=sorted(done),retained_paid=len(done),remaining=len(ordered),original_total=288,new_manifest=str(target),science='Same48 observations, six methods, initialization/prior/stages/selection/caps; no rerun or GT access',dispatch_order='Same within-case case-index-mod6 cyclic order after removing all paid tasks',hardware_scope='Same Windows host and independent Pcores; worker GetProcessTimes includes whole child, phase/core retained',resource_evidence='RESOURCE_UPDATE_4/5/6 and root explicit allocation; no external process interrupted',physical_calls=0,GT_access='NONE',whole_child_CPU_seconds=time.process_time())
 receiptpath.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
