"""One actual H adjoint/directional derivative check; no Maxwell FW or GT."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,time
from collections import Counter
import numpy as np
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).absolute().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r11.src.cie_feedback_current import feedback_current_maps,NativeFourierAPCS,DiskPulseToeplitzGreen
from research.research_max_struct_r10.src.lowrank_denominator import LowRankDenominator
OUT=ROOT/'delegated/diagnostics'

def main():
    wall,cpu=time.perf_counter(),time.process_time();obj=json.loads((OUT/'STATE_INVENTORY.json').read_text())['objects'][0];f=Path(obj['folder']);pub={};c=Counter()
    with threadpool_limits(limits=1):
      for p in obj['portable_public_files']:
       with np.load(p) as z:pub.update(dict(z))
      with np.load(f/'stage_00_joint.npz') as z:m=z['m'];A=z['A']
      with np.load(f/'DENOMINATOR_BASIS.npz') as z:Ub=z['U'];B=z['B']
      with np.load(f/'SOURCE_METRICS.npz') as z:Os=np.diag(z['state_weights']);Od=np.diag(z['data_weights'])
      with np.load(f.parent/(f.name+'__adapter_inputs')/'RECEIVER_DPCS.npz') as z:Vd=z['Vd'];JD=z['JD']
      G,E,S,Y=pub['G'],pub['E0'],pub['S'],pub['Y'];shape=(24,24)
      op=LowRankDenominator(G,Ub,'galerkin',c);op.B=B.copy()
      U=NativeFourierAPCS(shape,3,Vd,c);green=DiskPulseToeplitzGreen(pub['coords'],shape,float(pub['radius']),float(pub['k']),c)
      H,HH,T=feedback_current_maps(op,U,green,m,1.,JD,E,green.apply(JD),c)
      rng=np.random.default_rng(2026110708);d=rng.normal(size=A.shape)+1j*rng.normal(size=A.shape);d/=np.linalg.norm(d)
      z=rng.normal(size=JD.shape)+1j*rng.normal(size=JD.shape);z/=np.linalg.norm(z)
      lhs=np.vdot(H(d),z);rhs=np.vdot(d,HH(z));adjerr=float(abs(lhs-rhs)/max(abs(lhs),abs(rhs),1e-300))
      def fun(a):
        D=S@JD-Y+S@U.apply(a);R=T+H(a)
        v=.5*float(np.vdot(D,D@Od).real+np.vdot(R,R@Os).real)
        g=U.adjoint(S.conj().T@(D@Od))+HH(R@Os)
        return v,g
      v,g=fun(A);h=1e-5;fd=(fun(A+h*d)[0]-fun(A-h*d)[0])/(2*h);analytic=float(np.vdot(g,d).real)
      derr=float(abs(fd-analytic)/max(abs(fd),abs(analytic),1e-12))
      rec=dict(status='PASS' if max(adjerr,derr)<1e-5 else 'FAIL',case_id=obj['case_id'],stage=0,physical_forward=0,physical_adjoint=0,truth_access='NONE',
        actual_production_H='src/cie_feedback_current.py::feedback_current_maps',adjoint_relative_error=adjerr,
        directional_derivative=dict(finite_difference=float(fd),analytic=analytic,relative_error=derr,step=h),
        prior='Original material prior is constant in this fixed-m current derivative and therefore contributes zero',
        operation_cost=dict(c),CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
      (OUT/'BOTHBLOCKS_WIRING_CHECK.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec))
if __name__=='__main__':main()
