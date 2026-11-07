"""Exposed actual reset-current directions: unresolved G-K transport, no FW or GT."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import sys,json,time
from pathlib import Path
from collections import Counter
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r10.src.lowrank_denominator import LowRankDenominator

def load(p):
 with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def energy(x,omega):return float(np.sum(abs(x)**2*omega[None,:]))
def main():
 wall,cpu=time.perf_counter(),time.process_time();counts=Counter();rows=[]
 inventory=json.loads((ROOT/'delegated/diagnostics/STATE_INVENTORY.json').read_text())
 for obj in inventory['objects'][:3]:
  pub={}
  for p in obj['portable_public_files']:pub.update(load(p))
  G=pub['G'];omega=load(Path(obj['folder'])/'SOURCE_METRICS.npz')['state_weights']
  kb=load(Path(obj['folder'])/'DENOMINATOR_BASIS.npz')
  op=LowRankDenominator(G,kb['U'],'galerkin',counts);op.B=kb['B']
  for st in obj['stages']:
   saved=load(ROOT/'RESULTS/current_curvature_exposed'/f"{obj['case_id']}_s{st['stage']}.npz")
   m=saved['m'];beta=float(st['metadata']['beta_current']);fac=op.factor(m,beta)
   for direction in ('scalar','green'):
    dJ=saved[direction+'_current_step'];gd=G@dJ;kd=op.action(dJ)
    counts['dense_G_diagnostic_complex_RHS']+=dJ.shape[1]
    resolvent_part=beta*dJ-beta**2*op.solve(fac,m[:,None]*dJ)
    remainder=-beta*op.solve(fac,m[:,None]*(gd-kd))
    h=resolvent_part+remainder
    en=energy(h,omega);cross=float(2*np.real(np.sum(resolvent_part.conj()*remainder*omega[None,:])))
    rows.append(dict(case_id=obj['case_id'],stage=st['stage'],direction=direction,beta=beta,
     current_step_energy=energy(dJ,omega),current_residual_step_energy=en,
     resolvent_part_energy=energy(resolvent_part,omega),unrepresented_remainder_energy=energy(remainder,omega),
     signed_cross_energy=cross,remainder_over_total_norm=float(np.sqrt(energy(remainder,omega)/en)),
     GminusK_over_G_step_norm=float(np.sqrt(energy(gd-kd,omega)/energy(gd,omega))),
     scope='Same exposed saved orthonormal-state reset direction; not native PR or an online predictor'))
 out=ROOT/'RESULTS/current_remainder_exposed';out.mkdir(parents=True,exist_ok=False)
 (out/'METRICS.json').write_text(json.dumps(rows,indent=2)+'\n')
 (out/'RECEIPT.json').write_text(json.dumps(dict(status='COMPLETE',truth_access='NONE',physical_forward=0,physical_adjoint=0,
  CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,auxiliary_counts=dict(counts),
  question='Does the finite-rank G-K remainder dominate the new current direction metric, beyond uniform full-G saturation?'),indent=2)+'\n')
 print(json.dumps(rows,indent=2))
if __name__=='__main__':main()
