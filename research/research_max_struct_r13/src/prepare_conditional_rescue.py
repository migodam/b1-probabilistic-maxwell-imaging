"""Construct a pre-frozen public residual-triggered rescue cohort, no truth."""
import argparse,json,time
from pathlib import Path
import numpy as np


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();r=a.root
    policy=json.loads((r/'configs/conditional_rescue_policy.json').read_text())
    core=json.loads((r/'configs/confirmation_core48.json').read_text())
    q=json.loads((r/'queues/confirmation_core48/RECEIPT.json').read_text())
    if q['pending'] or q['running']:raise RuntimeError('Parent cohort not complete')
    jobs={j['task_id']:j for j in q['jobs']};tasks=[];rows=[]
    dest=r/'data_public/rescue48';dest.mkdir(exist_ok=False)
    for task in core['tasks']:
        if task['method']!=policy['parent_method']:continue
        folder=r/'RESULTS'/policy['parent_campaign']/task['task_id']
        meta=json.loads((folder/'primary.json').read_text());receipt=json.loads((folder/'RECEIPT.json').read_text())
        trigger=bool(meta['residual_over_delta']>1)
        row=dict(case_id=task['case_id'],parent_task_id=task['task_id'],trigger=trigger,
            parent_metadata=meta,parent_receipt=str(folder/'RECEIPT.json'),
            parent_external_CPU=jobs[task['task_id']]['external'].get('child_process_CPU_seconds'),
            parent_status=receipt['status'])
        rows.append(row)
        if not trigger:continue
        with np.load(r/task['input_file'],allow_pickle=False) as z:pub={k:z[k].copy() for k in z.files}
        with np.load(folder/'primary.npz',allow_pickle=False) as z:pub['m_init']=z['m'].copy()
        path=dest/(task['case_id']+'.npz');np.savez_compressed(path,**pub)
        for name,algorithm in policy['methods'].items():
            tasks.append(dict(task_id=task['case_id']+'__'+name,case_id=task['case_id'],method=name,algorithm=algorithm,
                campaign=policy['campaign'],physics_kind='scalar_tm_2d',cohort='fresh48_conditional_rescue',
                input_file=str(path.relative_to(r)),seed=13013,parent_run=task['task_id']))
    result=dict(campaign=policy['campaign'],status='FROZEN_PUBLIC_RULE_INSTANTIATED_NO_GT',
        protocol=policy['protocol'],tasks=tasks,all48_deployment_rows=rows,
        policy_file='configs/conditional_rescue_policy.json',truth_access='NONE',new_forward=0,
        preparation_whole_process_CPU_seconds=time.process_time())
    with (r/'configs/confirmation_rescue48.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(cases=len(rows),triggered=sum(x['trigger'] for x in rows),tasks=len(tasks),truth_access='NONE')))


if __name__=='__main__':main()
