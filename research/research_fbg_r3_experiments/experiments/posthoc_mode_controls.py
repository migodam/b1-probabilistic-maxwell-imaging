"""Post-disclosure component diagnosis; does not revise the fresh protocol."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[k]='1'
from pathlib import Path
import argparse,json,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from e1_current_release import infer,HERE

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--cases',nargs='+',required=True)
    args=ap.parse_args()
    public=HERE/'data/native_confirmation_public'
    cases=[c for c in json.loads((public/'cases.json').read_text()) if c['case_id'] in args.cases]
    tag='e1_posthoc_mode_controls'
    tasks=[dict(**c,method=method,chart='cie',start=start,wall=4.,tag=tag,public_dir=str(public),
               backbone='observable',joint_steps=80,max_lu=120,max_rhs=2400,affine_mode='native_tsom')
           for c in cases for start in (0,1,4) for method in ('n1_random','n1_energy','n1_normal','n1_soft')]
    out=HERE/'results'/tag
    out.mkdir(exist_ok=False)
    (out/'PLAN.json').write_text(json.dumps(dict(status='FROZEN_BEFORE_COMPONENT_PREDICTIONS',
        classification='POST_DISCLOSURE_GAIN_CASE_COMPONENT_DIAGNOSTIC',
        selection_reason='Paid fresh task logs selected random dictionary mode. Test ordinary selectors on same disclosed gain object; not a new blind cohort.',
        no_retuning=True,unchanged_source='audit/NATIVE_CONFIRMATION_SOURCE.zip',truth_access_by_inference=False,tasks=tasks),indent=2))
    rows=[];started=time.perf_counter()
    with ProcessPoolExecutor(max_workers=4) as pool:
        for f in as_completed([pool.submit(infer,t) for t in tasks]):
            row=f.result();rows.append(row)
            with (out/'rows.jsonl').open('a') as h:h.write(json.dumps(row)+'\n')
    summary=dict(status='PREDICTIONS_CLOSED_TRUTH_NOT_READ',classification='POST_DISCLOSURE_GAIN_CASE_COMPONENT_DIAGNOSTIC',
                 rows=rows,count=len(rows),batch_wall=time.perf_counter()-started,workers=4,blas_threads_per_worker=1)
    (out/'SUMMARY.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps({k:v for k,v in summary.items() if k!='rows'}))
if __name__=='__main__':main()

