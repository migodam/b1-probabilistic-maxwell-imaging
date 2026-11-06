"""Portable endpoint grading after every reproduction task is terminal.

Uses bundled offline labels only after saving the complete public output list.
No forward/adjoint calls, no candidate selection from truth, no historical paths.
"""
import argparse,csv,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r7.experiments.score_primary import material_metrics,noise_metrics,fixed_prior

def main():
    wall,cpu=time.perf_counter(),time.process_time()
    ap=argparse.ArgumentParser();ap.add_argument('--manifests',nargs='+',required=True,type=Path)
    ap.add_argument('--cohort',choices=['development','confirmation'],required=True);ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args();tasks=[t for p in args.manifests for t in json.loads(p.read_text())['tasks']]
    entries=[]
    for t in tasks:
        folder=Path(t['output_dir']);rec=json.loads((folder/'RECEIPT.json').read_text())
        if rec['status'] not in ('SEALED_PUBLIC_PREDICTION','COMPLETED_PUBLIC_PREDICTION','FAILED'):
            raise ValueError('All requested predictions must be terminal before truth: '+t['task_id'])
        entries.append(dict(task=t,receipt=rec,prediction=str(folder/'prediction.npz')))
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'PUBLIC_SEAL.json').write_text(json.dumps(dict(status='ALL_REQUESTED_OUTPUTS_TERMINAL_BEFORE_GT',entries=entries),indent=2)+'\n')
    truth_dir=ROOT/'sealed_truth'/args.cohort
    meta={r['case_id']:r for r in json.loads((truth_dir/'score_metadata.json').read_text())['records']}
    rows=[]
    for e in entries:
        t,r=e['task'],e['receipt'];row=dict(task_id=t['task_id'],case_id=t['case_id'],status=r['status'],
            CPU_seconds=r.get('r9_total_cpu_seconds',r.get('process_CPU_seconds')))
        if r['status']!='FAILED':
            with np.load(truth_dir/(t['case_id']+'.npz')) as z:truth,support=z['m'],z['mask']
            with np.load(e['prediction']) as z:m,yp=z['m'],z['Y_prediction']
            public={}
            for name in t['public_files']:
                with np.load(name) as z:
                    for k in ('Y','sigma','noise_radius'):
                        if k in z:public[k]=z[k]
            row.update(material_metrics(m,truth,support,meta[t['case_id']]['amplitude']))
            row.update(noise_metrics(yp,public,meta[t['case_id']].get('same_grid',True)))
            row['original_Phi']=.5*row['relative_residual']**2+fixed_prior(m)
        rows.append(row)
    (args.output/'SCORES.json').write_text(json.dumps(rows,indent=2)+'\n')
    keys=sorted(set(k for r in rows for k in r))
    with (args.output/'SCORES.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
    (args.output/'RECEIPT.json').write_text(json.dumps(dict(status='OFFLINE_SCORED',cohort=args.cohort,tasks=len(tasks),
        truth_access='Only after all requested predictions terminal',new_physics_calls=0,CPU_seconds=time.process_time()-cpu,
        wall_seconds=time.perf_counter()-wall),indent=2)+'\n')

if __name__=='__main__':main()
