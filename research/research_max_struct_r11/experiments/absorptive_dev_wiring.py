"""One public-state check of new complex current adjoint and real-material gradient;0FW."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,time
from collections import Counter
import numpy as np
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).absolute().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO));OUT=ROOT/'delegated/diagnostics/absorptive_dev'
from research.research_max_struct_r11.src.absorptive_denominator import AbsorptiveDenominator
from research.research_max_struct_r11.src.cie_feedback_current import feedback_current_maps,NativeFourierAPCS,DiskPulseToeplitzGreen
from research.research_z_alg_r4_diversity.src.physics import Regularizer

def main():
 wall,cpu=time.perf_counter(),time.process_time();inv=json.loads((ROOT/'delegated/diagnostics/STATE_INVENTORY.json').read_text())['objects'][0];folder=Path(inv['folder']);public={}
 for p in inv['portable_public_files']:
  with np.load(p) as z:public.update(dict(z))
 with np.load(folder/'stage_00_joint.npz') as z:oldm=z['m'];J=z['J'];Q=z['Q'];A=z['A']
 with np.load(folder/'DENOMINATOR_BASIS.npz') as z:Ug=z['U'].copy();Bg=z['B'].copy()
 with np.load(folder/'SOURCE_METRICS.npz') as z:Os=np.diag(z['state_weights'])
 with np.load(folder.parent/(folder.name+'__adapter_inputs')/'RECEIVER_DPCS.npz') as z:Vd=z['Vd'];JD=z['JD']
 m=.9*oldm+.01;regobj=Regularizer((24,24));reg=lambda x:regobj.evaluate(x)[:2];checks=[];total=Counter();rng=np.random.default_rng(202611071122)
 d=rng.normal(size=A.shape)+1j*rng.normal(size=A.shape);d/=np.linalg.norm(d);z=rng.normal(size=J.shape)+1j*rng.normal(size=J.shape);z/=np.linalg.norm(z);direction=np.sin(np.arange(len(m))+.137);direction/=np.linalg.norm(direction)
 with threadpool_limits(limits=1):
  for name,green in [('absorptiveGreenBothBlocks',True),('absorptiveScalarBothBlocks',False)]:
   count=Counter();op=AbsorptiveDenominator(public['G'],Ug if green else np.empty((len(m),0),complex),'galerkin',count)
   if green:op.B=Bg.copy()
   U=NativeFourierAPCS((24,24),3,Vd,count);fft=DiskPulseToeplitzGreen(public['coords'],(24,24),float(public['radius']),float(public['k']),count)
   H,HH,T=feedback_current_maps(op,U,fft,m,1.,JD,public['E0'],fft.apply(JD),count);lhs=np.vdot(H(d),z);rhs=np.vdot(d,HH(z));adjerr=float(abs(lhs-rhs)/max(abs(lhs),abs(rhs),1e-300))
   f,g=op.value_gradient(J,Q,m,1.,Os,1.,reg);h=1e-5;fd=(op.value_gradient(J,Q,m+h*direction,1.,Os,1.,reg)[0]-op.value_gradient(J,Q,m-h*direction,1.,Os,1.,reg)[0])/(2*h);exact=float(g@direction);rel=float(abs(fd-exact)/max(abs(fd),abs(exact),1e-12))
   checks.append(dict(method=name,rank=op.rank,current_actual_H_adjoint_relative_error=adjerr,material_finite_difference=float(fd),material_analytic=exact,material_relative_error=rel,material_difference_step=h,includes_original_R=True,status='PASS' if adjerr<1e-10 and rel<1e-5 else 'FAIL',counts=dict(count)));total.update(count)
 receipt=dict(status='PASS' if all(x['status']=='PASS' for x in checks) else 'FAIL',checks=checks,public_state=str(folder/'stage_00_joint.npz'),material_test_point='0.9*saved_m+0.01, interior realbox, no GT',beta_magnitude=1.,shift='c=-1j*beta',outer_prefactor='positive real beta',physical_forward=0,physical_adjoint=0,GT_access='NONE',counts=dict(total),CPU_seconds=time.process_time()-cpu,near_exit_process_CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-wall)
 (OUT/'WIRING_CHECK.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
 if receipt['status']!='PASS':raise RuntimeError('New wiring failed; do not run8 until root notified')
if __name__=='__main__':main()
