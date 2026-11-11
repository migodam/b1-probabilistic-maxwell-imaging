"""Index completed paid public predictions and disjoint queue costs; no GT/FW."""
from pathlib import Path
import argparse,json,time,math
ROOT=Path(__file__).resolve().parent
METHODS=('cie_single','cie_dual','global_dual','a_plus','fullG19','dr','multistart','nn_direct')

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stage',choices=('dev6','confirm24'),required=True)
    a=ap.parse_args();out=ROOT/('SEALED_'+a.stage.upper()+'_PUBLIC_COHORT.json')
    if out.exists():raise FileExistsError(out)
    manifest=json.loads((ROOT/'data_public'/a.stage/'manifest.json').read_text())
    # Queue terminals, including old startup failures, remain distinct attempts.
    attempts=[];controllers=[]
    for path in (ROOT/'queues').glob('*/RECEIPT.json'):
        q=json.loads(path.read_text());controllers.append(dict(queue=str(path),status=q['status'],controller_CPU_seconds=q.get('controller_CPU_seconds')))
        for job in q['jobs']:attempts.append(dict(queue=str(path),**job))
    byid={}
    for attempt in attempts:
        if attempt['status']=='FINISHED':
            ident=attempt['task_id']
            if ident in byid:raise ValueError('Multiple completed paid predictions for same canonical task; root must explicitly adjudicate provenance: '+ident)
            byid[ident]=attempt
    rows=[]
    import numpy as np
    for case in manifest['records']:
        for method in METHODS:
            ident=case['case_id']+'__'+method+'__240CPU';task=ROOT/'RESULTS'/a.stage/ident
            job=byid.get(ident)
            if job is None:raise ValueError('Cohort incomplete; no finished controller terminal: '+ident)
            receipt=json.loads((task/'RECEIPT.json').read_text())
            baseline=method in METHODS[:5]
            prediction=task/('R11_CAP_PREDICTION.npz' if baseline else 'prediction.npz')
            with np.load(prediction,allow_pickle=False) as z:m=z['m']
            if m.shape!=(576,) or not np.isfinite(m).all():raise ValueError('Nonfinite/invalid paid material; preserved for root review: '+ident)
            phi=receipt.get('cap_physical_objective') if baseline else receipt.get('objective')
            if phi is None or not math.isfinite(float(phi)):raise ValueError('No finite paid public Phi: '+ident)
            h=task/'HISTORY.json'
            aux_iterations=sum('iteration' in r for r in json.loads(h.read_text())) if h.exists() else None
            secondary_path=task/'NOISE_FEASIBLE_MIN_PRIOR.json'
            secondary=None
            if secondary_path.exists():
                secondary_info=json.loads(secondary_path.read_text())
                secondary=dict(status=secondary_info['status'],metadata=str(secondary_path),
                    prediction=str(task/'NOISE_FEASIBLE_MIN_PRIOR.npz') if (task/'NOISE_FEASIBLE_MIN_PRIOR.npz').exists() else None,
                    prefixes=str(task/'NOISE_FEASIBLE_MIN_PRIOR_TIMEPREFIX.npz'),
                    history=str(task/'NOISE_FEASIBLE_MIN_PRIOR_HISTORY.json'),
                    primary_selection_unchanged=True,estimator='CN99_MIN_SAME_PRIOR_AUDIT_ONLY_V1')
            rows.append(dict(task_id=ident,case_id=case['case_id'],method=method,prediction=str(prediction),
                receipt=str(task/'RECEIPT.json'),status=receipt['status'],Phi=phi,
                completed_external_child_CPU_seconds=job.get('child_process_CPU_seconds'),
                controller_terminal_queue=job['queue'],auxiliary_iterations_from_paid_history=aux_iterations,
                secondary_candidate_selection=secondary,
                source_selection='Already frozen paid original-Phi incumbent; no new candidates, optimization, forward or GT',
                cost_prefix_artifacts=[str(p) for p in ((task/'R11_INCUMBENT_TRAJECTORY.npz',task/'R11_CPU_CHECKPOINTS.json') if baseline else (task/'paid_prefix.npz',task/'history.json')) if p.exists()]))
    failures=[dict(task_id=j['task_id'],queue=j['queue'],status=j['status'],returncode=j.get('returncode'),
                   child_CPU_seconds=j.get('child_process_CPU_seconds')) for j in attempts if j['status']!='FINISHED']
    result=dict(status='COMPLETE_PUBLIC_COHORT_SEALED_BEFORE_GT',stage=a.stage,independent_materials=len(manifest['records']),
        methods=list(METHODS),predictions=len(rows),rows=rows,truth_access='NONE',scoring='NOT_RUN',new_FW=0,
        cost=dict(stage_finished_child_CPU_seconds=sum(r['completed_external_child_CPU_seconds'] or 0 for r in rows),
            campaign_all_ended_attempt_CPU_seconds=sum(j.get('child_process_CPU_seconds') or 0 for j in attempts),
            campaign_failed_attempts=failures,controllers=controllers,
            scope='Stage costs and full campaign attempt/controller costs are overlapping views; do not add stages to campaign twice. Controller and transport final tails may be partly unmetered.'),
        seal_process_CPU_seconds_near_write=time.process_time())
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=result['status'],stage=a.stage,predictions=len(rows),path=str(out),GT='NONE')),flush=True)

if __name__=='__main__':main()
