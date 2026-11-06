"""Inventory terminal prediction artifacts using JSON/stat only; no GT/arrays."""
from pathlib import Path
from collections import Counter
import json,shutil,time
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]


def main():
    cpu,wall=time.process_time(),time.perf_counter();rows=[];queues=[]
    for cohort,manifest,tag,count in [('baseline','twod_tasks.json','REMOTE',120),
            ('denominator','denominator_confirmation_tasks.json','DENOMINATOR',48)]:
        path=ROOT/('delegated/twod/fetched_'+cohort+'_queue/RECEIPT.json')
        queue=json.loads(path.read_text())
        shutil.copy2(path,ROOT/('delegated/twod/'+tag+'_QUEUE_RECEIPT.json'))
        tasks=json.loads((ROOT/'configs'/manifest).read_text())['tasks']
        if len(tasks)!=count:raise ValueError('Frozen task count mismatch')
        job_ids={v['task_id'] for v in queue['jobs']}
        if job_ids!={v['task_id'] for v in tasks}:raise ValueError('Queue/manifest IDs differ')
        queues.append(dict(cohort=cohort,status=queue['status'],terminal=len(queue['jobs']),
            failed=sum(v['status']!='FINISHED' for v in queue['jobs']),
            whole_child_CPU_seconds=queue['completed_child_CPU_seconds'],
            queue_controller_CPU_seconds=queue['controller_CPU_seconds'],wall_seconds=queue['wall_seconds']))
        for task in tasks:
            p=REPO/task['output_dir'];receipt=json.loads((p/'RECEIPT.json').read_text())
            names=['prediction.npz','initial.npz','selected_before_corrector.npz',
                'HISTORY.json','CORRECTOR.json','BUDGET_STATES.npz','R10_BUDGET_STATES.npz',
                'BASIS_METADATA.json','BASIS_FREQUENCIES.npz','SOURCE_METRICS.npz',
                'AUXILIARY_CANDIDATES.npz','NUMERICAL_SCHEDULE.json','R10_DPCS_RULE.json']
            names += [f'stage_{i:02d}_{kind}.npz' for i in range(4) for kind in ['joint','physical']]
            missing=[name for name in names if not (p/name).is_file() or (p/name).stat().st_size==0]
            prep=p.parent/(p.name+'__adapter_inputs')
            if not (prep/'RECEIVER_DPCS.npz').is_file():missing.append('adapter_inputs/RECEIVER_DPCS.npz')
            corr=receipt.get('corrector') or {}
            rows.append(dict(task_id=task['task_id'],case_id=task['case_id'],cohort=cohort,
                algorithm_label=task.get('material_readout',{}).get('label',task['arm']),
                status=receipt['status'],missing=missing,
                auxiliary_iterations=receipt.get('auxiliary_iterations'),
                corrector_iterations=corr.get('iterations'),corrector_success=corr.get('success'),
                corrector_stop=corr.get('stop'),material_dof=receipt.get('material_dof'),
                physical_cost=receipt.get('physical_model_cost'),truth_access=receipt['truth_access']))
    valid=(len(rows)==168 and all(v['status']=='SEALED_PUBLIC_PREDICTION' and
        not v['missing'] and v['truth_access']=='NONE' for v in rows) and
        all(q['status']=='COMPLETED' and q['failed']==0 for q in queues))
    result=dict(status='ALL_168_LOCAL_PREDICTIONS_SEALED_COMPLETE' if valid else 'INCOMPLETE',
        tasks=len(rows),algorithm_counts=dict(Counter(v['algorithm_label'] for v in rows)),
        queues=queues,rows=rows,auxiliary_current_steps=sum(v['auxiliary_iterations'] or 0 for v in rows),
        corrector_stop_counts=dict(Counter(v['corrector_stop'] for v in rows)),
        whole_child_CPU_seconds=sum(q['whole_child_CPU_seconds'] for q in queues),
        truth_access='NONE',array_payloads_read=0,physical_calls=0,scoring_calls=0,
        GT_authorization='Root decides separately; this inventory is not a GT reader or scoring command',
        verification_scope='Task ID sets, terminal JSON statuses and required nonempty saved artifacts only; tar/gzip extraction already completed; no SHA checks',
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
    (ROOT/'delegated/twod/OUTPUT_SYNC_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','tasks','algorithm_counts','auxiliary_current_steps','whole_child_CPU_seconds','corrector_stop_counts')}))
    if not valid:raise RuntimeError('Incomplete inventory; no scoring gate')


if __name__=='__main__':main()
