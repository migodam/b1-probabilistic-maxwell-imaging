"""Historical public-state counterfactuals, no truth access or inverse replay."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,json,time,argparse,traceback
from collections import Counter
import numpy as np
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.src.lowrank_denominator import prepare_bases,LowRankDenominator,matrix_material_readout
from research.research_max_struct_r4.src.cie_onehop import OneHopAbsorption
from research.research_max_struct_r7.src.aplus_readout import mobius_material_readout
from research.research_fbg_r3_experiments.src.physics import MaxwellModel


def dump(p,v):
    Path(p).write_text(json.dumps(v,indent=2,allow_nan=False,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist())+'\n')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--historical-repo',required=True);args=parser.parse_args()
    old=Path(args.historical_repo);out=ROOT/'RESULTS/denominator_fixed_state';out.mkdir(parents=True,exist_ok=False)
    wall,cpu=time.perf_counter(),time.process_time();receipt=dict(status='RUNNING',truth_access='NONE',data_status='EXPOSED_HISTORICAL_PUBLIC_STATES',purpose='fixed-current counterfactual material updates; not imaging superiority')
    dump(out/'RECEIPT.json',receipt)
    rows=[];basis_setup=[];checks=[];physical=Counter();aux=Counter()
    try:
        with threadpool_limits(limits=1):
            for case in ('obs_b20d5a55677d4eca','obs_40ac5b79e8a74ac3'):
                task=json.loads((old/'research/research_max_struct_r9/RESULTS/development'/f'{case}__global_dual__fixed_001'/'PUBLIC_TASK.json').read_text())
                public={}
                for public_path in task['public_files']:
                    with np.load(public_path) as z:public.update(dict(z))
                if 'E0' not in public:public['E0']=public['E']
                model=MaxwellModel(public);Y=public['Y'];n=model.n
                bases,info=prepare_bases(model.G,model.S);basis_setup.append(dict(case_id=case,**info))
                np.savez_compressed(out/f'{case}_bases.npz',**bases)
                Os=np.diag(1/np.sum(abs(model.E)**2,axis=0));Ls=np.sqrt(Os)
                lo,hi=np.full(n,.02),np.full(n,20.)
                reg=lambda x:model.reg.evaluate(x)[:2]
                for origin in ('global_dual','a_plus'):
                    folder=old/'research/research_max_struct_r9/RESULTS/development'/f'{case}__{origin}__fixed_001'
                    for stage in (0,1):
                        with np.load(folder/f'stage_{stage:02d}_joint.npz') as z:m,J,Q=[z[k].copy() for k in ('m','J','Q')]
                        beta=1. if stage==0 else .2
                        phi,gphi=model.data_value_gradient(m,Y)
                        base=dict(case_id=case,origin=origin,stage=stage,beta=beta,old_Phi=phi)
                        cfg=dict(material_inner_maxiter=60,material_inner_ftol=1e-12,material_inner_gtol=1e-8,aplus_regularizer=model.reg,
                                 aplus_sweeps=2,aplus_arc_candidates=True,aplus_inherited_gtol=True,matrix_inner_maxiter=8)
                        for kind in ('a_plus','twofold_galerkin','green_equal_rank_galerkin','twofold_right'):
                            count=Counter();blocks=OneHopAbsorption(model.G,model.shape,count)
                            started=time.process_time()
                            if kind=='a_plus':
                                new,detail=mobius_material_readout(blocks,J,Q,m,beta,Os,Ls,1.,(lo,hi),reg,cfg,count)
                            else:
                                basis_name='green_equal_rank' if kind.startswith('green') else 'twofold'
                                projection='right' if kind.endswith('right') else 'galerkin'
                                op=LowRankDenominator(model.G,bases[basis_name],projection,count);cfg['matrix_denominator']=op
                                if len(checks)<2 and kind=='twofold_galerkin':
                                    x=.9*m+.1;direction=np.sin(np.arange(n)+.137);direction/=np.linalg.norm(direction)
                                    f,g=op.value_gradient(J,Q,x,beta,Os,1.,reg)
                                    h=1e-5;fp=op.value_gradient(J,Q,x+h*direction,beta,Os,1.,reg)[0];fm=op.value_gradient(J,Q,x-h*direction,beta,Os,1.,reg)[0]
                                    fd=(fp-fm)/(2*h);exact=float(g@direction)
                                    checks.append(dict(**base,finite_difference=fd,analytic=exact,relative_error=abs(fd-exact)/max(abs(fd),abs(exact),1e-12)))
                                new,detail=matrix_material_readout(blocks,J,Q,m,beta,Os,Ls,1.,(lo,hi),reg,cfg,count)
                            elapsed=time.process_time()-started;step=new-m
                            value=model.physical_score(new,Y);dot=float(step@gphi)
                            row=dict(**base,method=kind,new_Phi=value,delta_Phi=value-phi,linear_prediction=dot,
                                     linear_error=value-phi-dot,step_norm=float(np.linalg.norm(step)),
                                     cosine_negative_physical_gradient=-dot/max(np.linalg.norm(step)*np.linalg.norm(gphi),1e-300),
                                     readout_cpu_seconds=elapsed,material_readout=detail,cost=dict(count))
                            rows.append(row);aux.update(count)
                            np.savez_compressed(out/f'{case}_{origin}_s{stage}_{kind}.npz',m_before=m,m_after=new,step=step)
                            dump(out/'RESULTS.json',rows)
                physical.update(model.cost)
        receipt['status']='COMPLETED'
    except Exception:
        receipt.update(status='FAILED',traceback=traceback.format_exc());raise
    finally:
        receipt.update(wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,physical_cost=dict(physical),auxiliary_cost=dict(aux),basis_setup=basis_setup,derivative_wiring=checks,observations=len(rows))
        dump(out/'RECEIPT.json',receipt)


if __name__=='__main__':main()
