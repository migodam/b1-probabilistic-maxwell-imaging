"""Replay an existing public-data task plan under a new output tag."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[k]='1'
from pathlib import Path
import argparse,json,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from e1_current_release import infer,HERE

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--plan',required=True)
    ap.add_argument('--public-dir',required=True)
    ap.add_argument('--tag',required=True)
    ap.add_argument('--workers',type=int,default=4)
    args=ap.parse_args()
    plan=json.loads(Path(args.plan).read_text())
    tasks=[{**t,'tag':args.tag,'public_dir':str(Path(args.public_dir).resolve())} for t in plan['tasks']]
    out=HERE/'results'/args.tag
    out.mkdir(exist_ok=False)
    (out/'PLAN.json').write_text(json.dumps({**plan,'tasks':tasks,'replay_of':args.plan,'truth_access':False},indent=2))
    rows=[];started=time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(infer,t) for t in tasks]):
            row=future.result();rows.append(row)
            with (out/'rows.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
    summary=dict(status='PREDICTIONS_CLOSED_TRUTH_NOT_READ',rows=rows,count=len(rows),
                 batch_wall=time.perf_counter()-started,workers=args.workers,blas_threads_per_worker=1,replay_of=args.plan)
    (out/'SUMMARY.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps({k:v for k,v in summary.items() if k!='rows'}))
if __name__=='__main__':main()

