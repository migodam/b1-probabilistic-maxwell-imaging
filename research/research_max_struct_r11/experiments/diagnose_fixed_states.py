"""Ten frozen exposed endpoint counterfactuals; GT is a final offline scoring pass."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,json,time,traceback
from collections import Counter
import numpy as np
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).absolute().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r10.src.lowrank_denominator import LowRankDenominator
from research.research_max_struct_r2.src.cie_global import _baseline_readout
from research.research_max_struct_r11.src.material_readouts import full_material_readout,lagged_material_readout,frozen_value_gradient
from research.research_max_struct_r11.src.diagnostic_core import measure_start,cosine,defect_stats
OUT=ROOT/'delegated/diagnostics'

def dump(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist())+'\n')

def difference(a,b):return {k:a.get(k,0)-b.get(k,0) for k in set(a)|set(b) if a.get(k,0)!=b.get(k,0)}

def main():
    wall,cpu=time.perf_counter(),time.process_time()
    receipt=dict(status='RUNNING',historical_data_status='EXPOSED_PUBLIC_R10',truth_access='NONE',
       scope='10 states; 3 candidate readouts per state; no inverse replay; no threshold or rule selection',
       ordinary='scalar CIE full-real readout at same fixed J/Q/m/beta; mature inherited maxiter60',
       full_lagged_budget='same scalar seed and maxiter8 matrix inner; no rescue budget')
    dump(OUT/'DIAGNOSTIC_RECEIPT.json',receipt)
    inventory=json.loads((OUT/'STATE_INVENTORY.json').read_text())
    starts=[];rows=[];aux_total=Counter();physical_total=Counter()
    try:
        with threadpool_limits(limits=1):
            for obj in inventory['objects']:
                folder=Path(obj['folder']);public={}
                for path in obj['portable_public_files']:
                    with np.load(path) as z:public.update(dict(z))
                model=MaxwellModel(public);Y=public['Y'];reg=lambda x:model.reg.evaluate(x)[:2]
                with np.load(folder/'SOURCE_METRICS.npz') as z:Os=np.diag(z['state_weights']);Ls=np.diag(np.sqrt(z['state_weights']))
                with np.load(folder/'DENOMINATOR_BASIS.npz') as z:U=z['U'].copy();B=z['B'].copy()
                lo,hi=np.full(model.n,.02),np.full(model.n,20.)
                for state in obj['stages']:
                    stage=state['stage'];sid=f"{obj['case_id']}_s{stage}"
                    with np.load(state['path']) as z:m,J,Q=[z[k].copy() for k in ('m','J','Q')]
                    beta=state['metadata']['beta_material'];counts=Counter()
                    op=LowRankDenominator(model.G,U,'galerkin',counts);op.B=B.copy()
                    before=dict(model.cost);old_phi,gphi=model.data_value_gradient(m,Y)
                    start_physics=difference(dict(model.cost),before)
                    meta,arrays=measure_start(op,J,Q,m,beta,Os,1.,reg,gphi,(lo,hi))
                    consistency=float(np.linalg.norm(Q-(model.E+model.G@J))/max(np.linalg.norm(Q),1e-300))
                    smeta=dict(state_id=sid,case_id=obj['case_id'],stage=stage,mf=state['metadata']['basis_info']['mf'],
                        old_Phi=old_phi,physical_cost=start_physics,saved_Q_consistency_relative=consistency,
                        diagnostic_cost=dict(counts),**meta)
                    starts.append(smeta)
                    np.savez_compressed(OUT/f'{sid}_START.npz',m=m,J=J,Q=Q,**arrays)
                    aux_total.update(counts)
                    reference_factor=op.factor(m,beta)
                    cfg=dict(material_inner_maxiter=60,material_inner_ftol=1e-12,material_inner_gtol=1e-8,matrix_inner_maxiter=8,matrix_denominator=op)
                    for method in ('ordinary','full','lagged'):
                        count=Counter();op.cost=count;beg=time.process_time()
                        if method=='ordinary':
                            new,detail=_baseline_readout(J@Ls,Q@Ls,m,beta,np.ones(J.shape[1]),1.,(lo,hi),reg,True,cfg,count)
                        else:
                            fun=full_material_readout if method=='full' else lagged_material_readout
                            new,detail=fun(None,J,Q,m,beta,Os,Ls,1.,(lo,hi),reg,cfg,count)
                        elapsed=time.process_time()-beg;step=new-m
                        # Freeze candidate arrays before any truth file is opened.
                        np.savez_compressed(OUT/f'{sid}_{method}_CANDIDATE.npz',m_before=m,m_after=new,step=step)
                        before=dict(model.cost);new_phi=model.physical_score(new,Y)
                        phys=difference(dict(model.cost),before)
                        e=J-new[:,None]*Q;w=op.solve(op.factor(new,beta),e)
                        frozen_w=op.solve(reference_factor,e)
                        prior=float(reg(new)[0])
                        moving_loss=.5*beta**2*float(np.vdot(w,w@Os).real)+prior
                        frozen_loss=.5*beta**2*float(np.vdot(frozen_w,frozen_w@Os).real)+prior
                        linear=float(gphi@step)
                        row=dict(state_id=sid,case_id=obj['case_id'],stage=stage,method=method,
                            old_Phi=old_phi,new_Phi=new_phi,delta_Phi=new_phi-old_phi,linear_prediction=linear,
                            nonlinear_remainder=new_phi-old_phi-linear,step_norm=float(np.linalg.norm(step)),
                            step_inf=float(np.linalg.norm(step,np.inf)),cos_step_negative_Phi_gradient=cosine(step,-gphi),
                            box_lower_after=int(np.count_nonzero(new<=lo+1e-12)),box_upper_after=int(np.count_nonzero(new>=hi-1e-12)),
                            moving_material_loss=moving_loss,frozen_material_loss=frozen_loss,prior_value=prior,
                            defect_after=defect_stats(e,w,Os),frozen_weighted_defect_after_squared=float(np.vdot(frozen_w,frozen_w@Os).real),
                            material_readout=detail,readout_cpu_seconds=elapsed,auxiliary_cost=dict(count),physical_cost=phys,
                            truth_access_at_candidate_generation='NONE',GT_status='PENDING_OFFLINE_DIAGNOSTIC')
                        rows.append(row);aux_total.update(count)
                        dump(OUT/'START_DIAGNOSTICS.json',starts);dump(OUT/'FIXED_STATE_COUNTERFACTUALS.json',rows)
                physical_total.update(model.cost)
        dump(OUT/'CANDIDATES_SEALED.json',dict(status='GENERATED_BEFORE_TRUTH_ACCESS',candidates=len(rows),selection='NONE; all retained'))
        # This pass must not affect candidates, steps, optimizer settings or rules.
        truth_schema={}
        for obj in inventory['objects']:
            with np.load(obj['portable_truth_path_OFFLINE_ONLY']) as z:
                truth_schema[obj['case_id']]={k:list(z[k].shape) for k in z.files}
                key=next((k for k in ('m','m_true','truth','material','m_gt') if k in z.files),None)
                if key is None:raise KeyError('No explicit material truth key in '+str(z.files))
                gt=np.asarray(z[key],float).reshape(-1)
            for row in rows:
                if row['case_id']!=obj['case_id']:continue
                with np.load(OUT/f"{row['state_id']}_{row['method']}_CANDIDATE.npz") as z:
                    old=np.linalg.norm(z['m_before']-gt)/np.linalg.norm(gt)
                    new=np.linalg.norm(z['m_after']-gt)/np.linalg.norm(gt)
                row.update(GT_status='OFFLINE_DIAGNOSTIC',NRMSE_definition='Relative L2 error ||m-GT||/||GT||; not amplitude NRMSE',old_NRMSE=float(old),new_NRMSE=float(new),delta_NRMSE=float(new-old))
        dump(OUT/'FIXED_STATE_COUNTERFACTUALS.json',rows)
        receipt.update(status='COMPLETE',truth_access='OFFLINE_DIAGNOSTIC_AFTER_ALL_CANDIDATES_SAVED',truth_schema=truth_schema)
    except Exception:
        receipt.update(status='FAILED',traceback=traceback.format_exc());raise
    finally:
        receipt.update(wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,states=len(starts),candidates=len(rows),physical_cost=dict(physical_total),auxiliary_cost=dict(aux_total))
        dump(OUT/'DIAGNOSTIC_RECEIPT.json',receipt)
if __name__=='__main__':main()
