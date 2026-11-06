"""Correct ordinary has_reg wiring without repeating full/lagged or original gradient."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,time
from collections import Counter
import numpy as np
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).absolute().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r11.experiments.diagnose_fixed_states import dump,difference,MaxwellModel,LowRankDenominator,_baseline_readout,cosine,defect_stats,OUT
from research.research_max_struct_r11.src.diagnostic_core import projected_gradient

def main():
    if (OUT/'ORDINARY_REPAIR_RECEIPT.json').exists():
        raise RuntimeError('Repair already completed; refuse to overwrite audit history')
    wall,cpu=time.perf_counter(),time.process_time()
    rows=json.loads((OUT/'FIXED_STATE_COUNTERFACTUALS.json').read_text());starts=json.loads((OUT/'START_DIAGNOSTICS.json').read_text())
    inventory=json.loads((OUT/'STATE_INVENTORY.json').read_text());oldreceipt=json.loads((OUT/'DIAGNOSTIC_RECEIPT.json').read_text())
    dump(OUT/'SEED_ONLY_ATTEMPT_COUNTERFACTUALS.json',rows);dump(OUT/'SEED_ONLY_ATTEMPT_RECEIPT.json',oldreceipt)
    physical=Counter();aux=Counter();newrows=[]
    with threadpool_limits(limits=1):
      for obj in inventory['objects']:
        folder=Path(obj['folder']);public={}
        for p in obj['portable_public_files']:
            with np.load(p) as z:public.update(dict(z))
        model=MaxwellModel(public);Y=public['Y'];reg=lambda x:model.reg.evaluate(x)[:2]
        with np.load(folder/'SOURCE_METRICS.npz') as z:Os=np.diag(z['state_weights']);Ls=np.diag(np.sqrt(z['state_weights']))
        with np.load(folder/'DENOMINATOR_BASIS.npz') as z:U=z['U'].copy();B=z['B'].copy()
        lo,hi=np.full(model.n,.02),np.full(model.n,20.)
        for state in obj['stages']:
          sid=f"{obj['case_id']}_s{state['stage']}";oldrow=next(r for r in rows if r['state_id']==sid and r['method']=='ordinary')
          # Preserve actually executed unregularized scalar seed control.
          source=OUT/f'{sid}_ordinary_CANDIDATE.npz';source.rename(OUT/f'{sid}_ordinary_seed_unregularized_CANDIDATE.npz')
          with np.load(OUT/f'{sid}_START.npz') as z:m,J,Q,gphi=[z[k].copy() for k in ('m','J','Q','physical_gradient')]
          beta=state['metadata']['beta_material'];count=Counter();op=LowRankDenominator(model.G,U,'galerkin',count);op.B=B.copy()
          fac=op.factor(m,beta)
          cfg=dict(material_inner_maxiter=60,material_inner_ftol=1e-12,material_inner_gtol=1e-8)
          beg=time.process_time();new,detail=_baseline_readout(J@Ls,Q@Ls,m,beta,np.ones(J.shape[1]),1.,(lo,hi),reg,True,cfg,count)
          elapsed=time.process_time()-beg;step=new-m
          np.savez_compressed(source,m_before=m,m_after=new,step=step)
          before=dict(model.cost);phi=model.physical_score(new,Y);phys=difference(dict(model.cost),before)
          e=J-new[:,None]*Q;w=op.solve(op.factor(new,beta),e);fw=op.solve(fac,e);prior=float(reg(new)[0]);linear=float(gphi@step)
          newrows.append(dict(state_id=sid,case_id=obj['case_id'],stage=state['stage'],method='ordinary',old_Phi=oldrow['old_Phi'],new_Phi=phi,delta_Phi=phi-oldrow['old_Phi'],linear_prediction=linear,
            nonlinear_remainder=phi-oldrow['old_Phi']-linear,step_norm=float(np.linalg.norm(step)),step_inf=float(np.linalg.norm(step,np.inf)),cos_step_negative_Phi_gradient=cosine(step,-gphi),
            box_lower_after=int(np.count_nonzero(new<=lo+1e-12)),box_upper_after=int(np.count_nonzero(new>=hi-1e-12)),
            moving_material_loss=.5*beta**2*float(np.vdot(w,w@Os).real)+prior,frozen_material_loss=.5*beta**2*float(np.vdot(fw,fw@Os).real)+prior,prior_value=prior,
            defect_after=defect_stats(e,w,Os),frozen_weighted_defect_after_squared=float(np.vdot(fw,fw@Os).real),material_readout=detail,readout_cpu_seconds=elapsed,auxiliary_cost=dict(count),physical_cost=phys,
            truth_access_at_candidate_generation='Previously exposed offline GT exists; no truth loaded or used in this repair generation pass',GT_status='PENDING_OFFLINE_DIAGNOSTIC'))
          # Matched scalar CIE gradient plus same prior, with no additional FW.
          R=beta*(J-m[:,None]*Q)/(1+beta*m[:,None]);dR=-beta*(Q+beta*J)/(1+beta*m[:,None])**2
          go=np.real(np.sum(dR.conj()*(R@Os),axis=1))+reg(m)[1];pg=projected_gradient(m,go,lo,hi)
          smeta=next(r for r in starts if r['state_id']==sid)
          smeta['gradients']['g_ordinary']=dict(norm=float(np.linalg.norm(go)),inf=float(np.linalg.norm(go,np.inf)),descent_cos_negative_physical_gradient=cosine(-go,-gphi),
            box_projected_gradient_norm=float(np.linalg.norm(pg)),box_projected_gradient_inf=float(np.linalg.norm(pg,np.inf)),projected_descent_cos_negative_physical_gradient=cosine(-pg,-gphi))
          aux.update(count)
        physical.update(model.cost)
    dump(OUT/'ORDINARY_REPAIR_CANDIDATES_SEALED.json',dict(candidates=len(newrows),status='All saved before repair scoring pass',settings='has_reg=True; maxiter60; inherited prior and bounds'))
    for obj in inventory['objects']:
      with np.load(obj['portable_truth_path_OFFLINE_ONLY']) as z:gt=z['m']
      for row in newrows:
        if row['case_id']!=obj['case_id']:continue
        with np.load(OUT/f"{row['state_id']}_ordinary_CANDIDATE.npz") as z:
            old=float(np.linalg.norm(z['m_before']-gt)/np.linalg.norm(gt));new=float(np.linalg.norm(z['m_after']-gt)/np.linalg.norm(gt))
        row.update(GT_status='OFFLINE_DIAGNOSTIC',NRMSE_definition='Relative L2 error ||m-GT||/||GT||; not amplitude NRMSE',old_NRMSE=old,new_NRMSE=new,delta_NRMSE=new-old)
    final=[next(x for x in newrows if x['state_id']==r['state_id']) if r['method']=='ordinary' else r for r in rows]
    dump(OUT/'FIXED_STATE_COUNTERFACTUALS.json',final);dump(OUT/'START_DIAGNOSTICS.json',starts)
    repair=dict(status='COMPLETE',reason='Initial ordinary has_reg=False produced scalar seed only; corrected to True without repeating full/lagged or original Phi gradients',
       states=10,candidates=10,physical_cost=dict(physical),auxiliary_cost=dict(aux),wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu)
    dump(OUT/'ORDINARY_REPAIR_RECEIPT.json',repair)
    pc=Counter(oldreceipt['physical_cost']);pc.update(physical);ac=Counter(oldreceipt['auxiliary_cost']);ac.update(aux)
    oldreceipt.update(physical_cost=dict(pc),auxiliary_cost=dict(ac),correction='See ORDINARY_REPAIR_RECEIPT; original 10 scalar-seed FW retained in total',ordinary='scalar CIE full-real inherited prior readout, has_reg=True; maxiter60',total_executed_candidates_including_seed_attempt=40)
    dump(OUT/'DIAGNOSTIC_RECEIPT.json',oldreceipt)
if __name__=='__main__':main()
