"""Single actual contracted-readout/prior wiring directional derivative check; zero FW."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,time
from collections import Counter
import numpy as np
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).absolute().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r11.src import material_readouts as mod
from research.research_max_struct_r10.src.lowrank_denominator import LowRankDenominator
from research.research_fbg_r3_experiments.src.physics import MaxwellModel

def main():
    out=ROOT/'delegated/diagnostics';inv=json.loads((out/'STATE_INVENTORY.json').read_text())['objects'][0];folder=Path(inv['folder'])
    with threadpool_limits(limits=1):
      public={}
      for p in inv['portable_public_files']:
       with np.load(p) as z:public.update(dict(z))
      model=MaxwellModel(public);reg=lambda x:model.reg.evaluate(x)[:2]
      with np.load(folder/'stage_00_joint.npz') as z:m,J,Q=[z[k] for k in ('m','J','Q')]
      with np.load(folder/'SOURCE_METRICS.npz') as z:Os=np.diag(z['state_weights']);Ls=np.diag(np.sqrt(z['state_weights']))
      with np.load(folder/'DENOMINATOR_BASIS.npz') as z:U=z['U'];B=z['B']
      c=Counter();op=LowRankDenominator(model.G,U,'galerkin',c);op.B=B.copy();lo=np.full(model.n,.02);hi=np.full(model.n,20.)
      saved=mod.minimize;checks=[]
      def observed(fun,start,*args,**kw):
        b=np.asarray(kw['bounds']);x=.85*start+.15*b.mean(axis=1)
        d=np.sin(np.arange(len(x))+.137);d/=np.linalg.norm(d)
        h=1e-6;v,g=fun(x);fd=(fun(x+h*d)[0]-fun(x-h*d)[0])/(2*h);analytic=float(g@d)
        checks.append(dict(coordinate='z',value=float(v),finite_difference=float(fd),analytic=analytic,relative_error=abs(fd-analytic)/max(abs(fd),abs(analytic),1e-12),step=h,
          has_original_material_prior=True,prior_value=float(reg(x/(1-x))[0])))
        return saved(fun,start,*args,**kw)
      mod.minimize=observed
      try:new,detail=mod.contracted_lagged_material_readout(None,J,Q,m,1.,Os,Ls,1.,(lo,hi),reg,dict(matrix_denominator=op,matrix_inner_maxiter=8,material_inner_gtol=1e-8),c)
      finally:mod.minimize=saved
      receipt=dict(status='PASS' if checks[0]['relative_error']<1e-5 else 'FAIL',case_id=inv['case_id'],stage=0,physical_forward=0,physical_adjoint=0,truth_access='NONE',checks=checks,
        readout=detail,auxiliary_cost=dict(c),result_box_lower_excess=float(max(0,np.max(lo-new))),result_box_upper_excess=float(max(0,np.max(new-hi))))
      (out/'CONTRACTED_WIRING_CHECK.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
