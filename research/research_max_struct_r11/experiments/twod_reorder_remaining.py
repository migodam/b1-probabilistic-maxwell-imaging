"""Reorder only unstarted frozen tasks after explicit no-kill STOP_DISPATCH."""
import argparse,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--old-queue',type=Path,required=True);ap.add_argument('--manifest',type=Path,required=True);a=ap.parse_args();cpu=time.process_time()
 q=json.loads(a.old_queue.read_text());original=json.loads(a.manifest.read_text());tasks=original['tasks'];done={j['task_id'] for j in q['jobs']}
 if q['status']!='STOPPED_DISPATCH_CHILDREN_FINISHED' or q.get('running'):raise RuntimeError('Wait for old own queue inflight to finish naturally')
 if len(done)!=len(q['jobs']):raise ValueError('Duplicate completed task')
 unstarted=set(q['not_started_task_ids']);lookup={t['task_id']:t for t in tasks}
 if unstarted!=set(lookup)-done:raise ValueError('Unstarted inventory mismatch')
 methods=('cie_single','cie_dual','global_dual','fullG19','freezeA','contractedT');cases=sorted({t['case_id'] for t in tasks});ordered=[]
 for case in cases:
  rotation=int(case.rsplit('_',1)[-1])%6;order=methods[rotation:]+methods[:rotation]
  for method in order:
   found=[t for t in tasks if t['case_id']==case and t.get('material_readout',{}).get('label',t['arm'])==method]
   if len(found)!=1:raise ValueError('Expected one frozen case/method')
   if found[0]['task_id'] in unstarted:ordered.append(found[0])
 if {t['task_id'] for t in ordered}!=unstarted or len(ordered)!=len(unstarted):raise ValueError('Reorder duplicated/dropped task')
 for t in ordered:
  p=REPO/t['output_dir']
  if p.exists():raise FileExistsError('Unstarted output already exists; preserve and audit')
 target=ROOT/'configs/twod_independent48_interleaved_remaining_tasks.json'
 receiptpath=ROOT/'delegated/twod/independent48/DISPATCH_REORDER_RECEIPT.json'
 for p in (target,receiptpath):
  if p.exists():raise FileExistsError(p)
 doc=dict(status='FROZEN_EXECUTION_ORDER_CORRECTION_NO_SCIENTIFIC_CHANGE',input_base='.',tasks=ordered)
 target.write_text(json.dumps(doc,indent=2)+'\n')
 receipt=dict(status='REORDERED_UNSTARTED_ONLY',reason='Root requested reduce method/time-temperature-load confounding before confirmation GT access',rule='Within each case rotate six methods by case index modulo6; preserve case order; remove all already-paid tasks',methods=list(methods),previous_dispatch='All baseline144 then candidates144',old_queue=str(a.old_queue),completed_task_ids=sorted(done),retained_completed=len(done),remaining=len(ordered),first_remaining_task_ids=[t['task_id'] for t in ordered[:18]],new_manifest=str(target),scientific_task_contents_unchanged=True,new_materials=0,repeated_inversions=0,GT_access='NONE',physical_calls=0,known_limitation='Already-started front segment retains original baseline-first time exposure; never silently rerun/drop it',whole_child_CPU_seconds=time.process_time(),control_work_CPU_seconds=time.process_time()-cpu)
 receiptpath.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
