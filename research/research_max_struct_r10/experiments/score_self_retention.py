"""Offline scoring of two completed, exposed historical DEV controls."""
import argparse,json,time,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r7.experiments.score_primary import material_metrics


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--historical-repo',required=True);args=parser.parse_args()
    r9=Path(args.historical_repo)/'research/research_max_struct_r9'
    started=time.process_time();rows=[]
    for case in ('obs_b20d5a55677d4eca','obs_40ac5b79e8a74ac3'):
        p=ROOT/'RESULTS/denominator_dev'/f'{case}__green_strong__galerkin__retain_self'
        receipt=json.loads((p/'RECEIPT.json').read_text())
        if receipt['status']!='SEALED_PUBLIC_PREDICTION':raise ValueError('Prediction not sealed')
        with np.load(r9/'sealed_truth/development'/f'{case}.npz') as t:
            for state in ('selected_before_corrector','prediction'):
                with np.load(p/(state+'.npz')) as z:m=z['m']
                rows.append(dict(case_id=case,state=state,**material_metrics(m,t['m'],t['mask'],12.)))
    result=dict(status='COMPLETED_POSTSEAL_HISTORICAL_DEV',truth_access='EXPOSED_OFFLINE_DIAGNOSTIC',
        comparison='Retain exact diagonal G while replacing nonlocal remainder by rank19 projection',
        CPU_seconds=time.process_time()-started,new_physics_calls=0,rows=rows)
    (ROOT/'RESULTS/denominator_dev/SELF_RETENTION_SCORE.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
