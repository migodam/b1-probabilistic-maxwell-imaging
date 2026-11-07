"""Postseal scoring of the two exposed, optimizer-only DEV controls."""
import argparse,json,time,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r7.experiments.score_primary import material_metrics


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--historical-repo',required=True);args=parser.parse_args()
    old=Path(args.historical_repo)/'research/research_max_struct_r9'
    cases=('obs_b20d5a55677d4eca','obs_40ac5b79e8a74ac3')
    root=ROOT/'RESULTS/onehop_solver_control';paths={}
    for case in cases:
        p=root/f'{case}__onehop_lbfgs8_corrected_dispatch'
        receipt=json.loads((p/'RECEIPT.json').read_text())
        if receipt['status']!='SEALED_PUBLIC_PREDICTION':raise ValueError('Both controls must seal first')
        schedule=json.loads((p/'NUMERICAL_SCHEDULE.json').read_text())
        if schedule['material_inner_maxiter']!=8:raise ValueError('Not the frozen eight-step control')
        paths[case]=(p,receipt)
    cpu=time.process_time();rows=[]
    for case,(p,receipt) in paths.items():
        with np.load(old/'sealed_truth/development'/f'{case}.npz',allow_pickle=False) as truth:
            for name in ('selected_before_corrector','prediction'):
                with np.load(p/(name+'.npz'),allow_pickle=False) as z:m=z['m']
                rows.append(dict(case_id=case,state=name,**material_metrics(m,truth['m'],truth['mask'],12.),
                    inverse_CPU_seconds=receipt['cpu_seconds'],inverse_cost=receipt['physical_model_cost']))
    result=dict(status='COMPLETED_POSTSEAL_HISTORICAL_DEV',truth_access='EXPOSED_OFFLINE_DIAGNOSTIC',
        question='Does eight-step LBFGS with unchanged onehop metric explain the new denominator rescue?',
        CPU_seconds=time.process_time()-cpu,new_physics_calls=0,rows=rows)
    (root/'SCORE.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
