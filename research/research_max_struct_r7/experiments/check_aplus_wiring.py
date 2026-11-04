"""Development wiring audit using already stored public current freezes; no GT."""
from pathlib import Path
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import sys,json,time
import numpy as np
from collections import Counter
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from research.research_max_struct_r6.experiments.run_fft_case import load_public
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r4.src.cie_onehop import OneHopAbsorption,onehop_material_readout,onehop_material_value_gradient
from research.research_max_struct_r7.src.aplus_readout import mobius_material_readout,projected_material_readout
from threadpoolctl import threadpool_limits
OUT=ROOT/'research/research_max_struct_r7/delegated/aplus';OUT.mkdir(exist_ok=True,parents=True)

def main():
 rows=[];derivative_rows=[]
 with threadpool_limits(limits=1):
  for case in ('obs_b20d5a55677d4eca','obs_7b8d5df517f1497d'):
   candidates=list((ROOT/'research/research_max_struct_r6/RESULTS/primary').glob(case+'__a_material'))
   if not candidates:
    print('Unavailable public freeze',case);continue
   folder=candidates[0];task=json.loads((folder/'PUBLIC_TASK.json').read_text());model=MaxwellModel(load_public(task['public_files']))
   for stage in range(3):
    freeze=folder/f'stage_{stage:02d}_joint.npz'
    if not freeze.is_file():continue
    data=np.load(freeze);m,J,Q=data['m'],data['J'],data['Q'];beta=1. if stage==0 else .2
    weights=1/np.sum(abs(model.E)**2,axis=0);Os=np.diag(weights);Ls=np.diag(np.sqrt(weights));bounds=(np.full(model.n,.02),np.full(model.n,20.))
    audit_cost=Counter();audit_blocks=OneHopAbsorption(model.G,model.shape,audit_cost);idx=audit_blocks.indices;b=beta-audit_blocks.diag_g;z=m[idx]/(1+b[idx]*m[idx]);Jb=(J@Ls)[idx];H=((Q+b[:,None]*J)@Ls)[idx];V=Jb-z[:,:,None]*H;R=beta*(V+z[:,:,None]*(audit_blocks.goff@V));state_gradient=np.zeros(model.n);arc_errors=[]
    for color in range(4):
     B=-beta*z[:,:,None]*audit_blocks.goff[:,:,color,None]*H[:,color,None,:];B[:,color,:]=beta*(-H[:,color,:]+(audit_blocks.goff@V)[:,color,:]);A=R-z[:,color,None,None]*B;q=np.sum(abs(B)**2,axis=(1,2));c=np.sum(B.conj()*A,axis=(1,2));cells=idx[:,color];bc=b[cells];mc=m[cells];n1=c.real;n2=.5*q+(c*bc).real;d1=2*bc.real;d2=abs(bc)**2;state_gradient[cells]=(n1+2*n2*mc+(n2*d1-n1*d2)*mc**2)/abs(1+bc*mc)**4
     test=m.copy();cell=int(cells[0]);test[cell]=12.;newz=12./(1+b[cell]*12.);pred=R.copy();pred[0]+=(newz-z[0,color])*B[0];actual=(audit_blocks.apply(test,beta,J-test[:,None]*Q)@Ls)[idx];arc_errors.append(float(np.linalg.norm(actual-pred)/max(np.linalg.norm(actual),1e-300)))
    actual_f,actual_g=onehop_material_value_gradient(audit_blocks,J,Q,m,beta,Os,1.,lambda zz:model.reg.evaluate(zz)[:2],audit_cost);new_g=state_gradient+model.reg.evaluate(m)[1];gerror=float(np.max(abs(actual_g-new_g)));assert gerror<1e-10+1e-6*np.max(abs(actual_g)),gerror;assert max(arc_errors)<1e-12,arc_errors
    derivative_rows.append(dict(case_id=case,stage=stage,gradient_max_abs_error=gerror,gradient_reference_max_abs=float(np.max(abs(actual_g))),conditional_arc_relative_errors=arc_errors))
    for name,fn in [('original_a',onehop_material_readout),('mobius_mm',mobius_material_readout),('projected_adjoint',projected_material_readout),('mobius_guarded',mobius_material_readout),('mobius_arc_candidates',mobius_material_readout)]:
     cost=Counter();blocks=OneHopAbsorption(model.G,model.shape,cost);cfg=dict(aplus_regularizer=model.reg,aplus_inherited_gtol=name in ('mobius_guarded','mobius_arc_candidates'),aplus_arc_candidates=name=='mobius_arc_candidates',cache_material_value_gradient=True,material_inner_maxiter=60,material_inner_ftol=1e-12,material_inner_gtol=1e-8)
     cpu=time.process_time();new,info=fn(blocks,J,Q,m,beta,Os,Ls,1.,bounds,lambda z:model.reg.evaluate(z)[:2],cfg,cost)
     row=dict(case_id=case,stage=stage,kind=name,cpu_seconds=time.process_time()-cpu,info=info,cost=dict(cost),finite=bool(np.isfinite(new).all()),in_box=bool(np.all((new>=.02)&(new<=20))),descent=bool(info['material_final_block_loss']<=info['material_old_block_loss']+1e-12))
     rows.append(row)
     assert row['finite'] and row['in_box'] and row['descent'],row
     np.savez_compressed(OUT/f'freeze_{case}_s{stage}_{name}.npz',m_old=m,m_new=new,J=J,Q=Q)
 (OUT/'WIRING_CHECKS.json').write_text(json.dumps(dict(status='PASS',truth_access='NONE',physical_LU=0,physical_RHS=0,gradient_wiring=derivative_rows,rows=rows),indent=2))
 print(json.dumps([dict(case=r['case_id'],stage=r['stage'],kind=r['kind'],cpu=r['cpu_seconds'],loss=r['info']['material_final_block_loss'],old=r['info']['material_old_block_loss']) for r in rows],indent=2))
if __name__=='__main__':main()
