from pathlib import Path
import sys,json,time,argparse
from datetime import datetime,timezone
import numpy as np
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_max_struct_r3.src.passive_continuation import PassiveModel
from research.research_fbg_r3_experiments.src.full_maxwell_corrector import corrector
from research.research_max_struct_r2.experiments.run_blind import clean

def run(arm):
    dest=ROOT/'results/development'/arm;dest.mkdir(parents=True,exist_ok=True)
    if (dest/'RECEIPT.json').exists():raise FileExistsError('Preserve existing run')
    wall,cpu=time.perf_counter(),time.process_time()
    case='obs_b20d5a55677d4eca';r2=REPO/'research/research_max_struct_r2'
    with np.load(r2/'data_public'/f'{case}.npz') as z:public={k:z[k] for k in z.files}
    with np.load(r2/'predictions'/case/'native.npz') as z:m=z['m'].copy()
    Y=public['Y'];m0=m.copy();stages=[];models=[]
    rec=dict(run_id=arm,case_id=case,classification='PUBLISHED_R2_DEVELOPMENT',truth_access='NONE',
        original_seed='R2 native corrected endpoint',physical_convention='exp(-i omega t), scalar TM disk pulse VIE',
        mesh=[24,24],sources=12,receivers=48,material_class='real [.02,20] final; auxiliary Im contrast=gamma>=0',
        prior='same R2 smooth and TV',online_selection_rule='original full Maxwell Phi only',
        code_hash=None,config_hash=None,hash_reason='User requested no hash checking',
        started_utc=datetime.now(timezone.utc).isoformat())
    try:
        schedule=[(4.,250),(1.,250),(.25,250),(0.,1000)] if arm=='passive_path' else [(0.,1750)]
        for gamma,iterations in schedule:
            model=PassiveModel(public,gamma);models.append(model)
            out=corrector(model,Y,m,dict(method='lbfgs',max_nfev=iterations,box=(.02,20.)))
            m=out['m'];stages.append(dict(gamma=gamma,iterations=iterations,stop=out['stop'],cost=dict(model.cost),
                wall=out['wall'],history=out['history']))
            np.savez_compressed(dest/f'gamma_{gamma}.npz',m=m)
        original=PassiveModel(public,0.);models.append(original)
        candidates=[m0,m];values=[original.physical_score(x,Y) for x in candidates]
        best=candidates[int(np.argmin(values))];state=original.state(best)
        value,g=original.data_value_gradient(best,Y)
        pg=best-np.clip(best-g,.02,20.)
        np.savez_compressed(dest/'prediction.npz',m=best,Y_prediction=state['Y'],j=state['j'],Q=state['Q'])
        rec.update(status='SEALED_DEVELOPMENT_PREDICTION',stages=stages,physical_objective=value,
            physical_residual=float(np.linalg.norm(state['Y']-Y)/np.linalg.norm(Y)),
            physical_whitened_residual=float(np.linalg.norm(state['Y']-Y)/float(public['sigma'])),
            projected_gradient_inf=float(np.max(abs(pg))),candidate_ids=['R2_native','path_endpoint'],
            constraint_violation=float(max(0,.02-best.min(),best.max()-20)),constitutive_defect=state['state_relative_residual'])
    except Exception as ex:
        rec.update(status='FAILED',error=repr(ex));raise
    finally:
        keys=set().union(*(q.cost for q in models))
        rec.update(cost={k:sum(q.cost.get(k,0) for q in models) for k in keys},wall_seconds=time.perf_counter()-wall,
            cpu_seconds=time.process_time()-cpu,finished_utc=datetime.now(timezone.utc).isoformat())
        (dest/'RECEIPT.json').write_text(json.dumps(clean(rec),indent=2,allow_nan=False))
    print(json.dumps({k:rec[k] for k in ['run_id','status','physical_residual','projected_gradient_inf','cost','wall_seconds']}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('arm',choices=['passive_path','extra_original']);args=p.parse_args()
    with threadpool_limits(limits=1):run(args.arm)
