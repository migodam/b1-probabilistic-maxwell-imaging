"""One fixed public-input derivative/adjoint wiring check; no physics or optimizer."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import time
cpu,wall=time.process_time(),time.perf_counter()
import sys,json,csv,traceback
from pathlib import Path
from collections import Counter
R=Path(__file__).resolve().parents[2];W=R.parents[1]
sys.path.insert(0,str(W))
import numpy as np
from research.research_max_struct_r15.src.twofold_denominator import TwofoldDenominator,matrix_material_readout,twofold_material_readout
from research.research_max_struct_r10.src.lowrank_denominator import LowRankDenominator
from research.research_z_alg_r4_diversity.src.physics import Regularizer
out=Path(__file__).resolve().parent
receipt=dict(status='RUNNING',truth_access='NONE',physics_calls=0,physical_LU=0,physical_RHS=0,new_SVD=0,
             online_inverse_runs=0,optimizers_run=0,errors=[],qualification='No all-box contraction guarantee for truncated R correction')
try:
 public=W/'research/research_max_struct_r14/data_public/independent24/r14_new_0000.npz'
 spaces=R/'audit/rank_characterization/current_complex_v1/operator_group_00/SPACES.npz'
 with np.load(public,allow_pickle=False) as z:G=z['G'];E=z['E0'];shape=tuple(z['grid_shape'])
 with np.load(spaces,allow_pickle=False) as z:U=z['U11'];Ug=z['G_right_frame'][:,:19];Us=z['S_right_frame'][:,:19]
 n,ns=E.shape;index=np.arange(n,dtype=float)
 # Fixed auxiliary inputs, not a physical PDE solution or saved GT candidate.
 J=.1*E;Q=E+G@J
 Os=np.diag(1/np.sum(abs(E)**2,axis=0))
 regularizer=Regularizer(shape);reg=lambda m:regularizer.evaluate(m)[:2]
 materials={'uniform':np.full(n,.1),'deterministic_nonuniform':.1+.035*np.sin(.17*index)+.02*np.cos(.31*index)}
 directions=[np.sin(.23*(index+1)),np.cos(.37*(index+1))]
 directions=[h/np.linalg.norm(h) for h in directions]
 beta=1.;rho=1.;epsilon=1e-4
 cost=Counter(auxiliary_public_G_actions=1,auxiliary_public_G_complex_RHS=ns)
 rows=[];equality=[];adjoints=[];operator_receipts=[]
 for name,basis in [('U11',U),('rank0',np.empty((n,0),complex))]:
  op=TwofoldDenominator(G,basis,shape,one_hop=True,cost=cost)
  operator_receipts.append(dict(basis=name,**op.setup_receipt))
  for mname,m in materials.items():
   v,g=op.value_gradient(J,Q,m,beta,Os,rho,reg)
   for di,h in enumerate(directions):
    vp,_=op.value_gradient(J,Q,m+epsilon*h,beta,Os,rho,reg)
    vm,_=op.value_gradient(J,Q,m-epsilon*h,beta,Os,rho,reg)
    analytic=float(g@h);numeric=(vp-vm)/(2*epsilon);err=abs(analytic-numeric)
    rows.append(dict(basis=name,material=mname,direction=di,epsilon=epsilon,value=v,analytic=analytic,central_FD=numeric,absolute_error=err,relative_error=err/max(abs(analytic),abs(numeric),1e-15)))
   # Deterministic complex dual pairing; no random physics or additional basis.
   X=E[:,:2]*(1+.01j*index[:,None]);Y=E[:,-2:]*(.4+.003j*index[:,None])
   TX=op.transform(m,beta,X);aY=op.adjoint(m,beta,Y)
   lhs=np.vdot(Y,TX);rhs=np.vdot(aY,X)
   adjoints.append(dict(basis=name,material=mname,absolute_error=float(abs(lhs-rhs)),relative_error=float(abs(lhs-rhs)/max(abs(lhs),abs(rhs),1e-15))))
 for name,basis in [('U11',U),('G19_control',Ug),('S19_control',Us)]:
  op=TwofoldDenominator(G,basis,shape,one_hop=False,cost=cost)
  original=LowRankDenominator(G,basis,projection='galerkin',cost=cost)
  cost['reference_setup_G_times_U_complex_columns']+=basis.shape[1]
  cost['reference_setup_basis_products']+=2
  for mname,m in materials.items():
   v,g=op.value_gradient(J,Q,m,beta,Os,rho,reg);vo,go=original.value_gradient(J,Q,m,beta,Os,rho,reg)
   equality.append(dict(basis=name,material=mname,value_exact_equal=v==vo,gradient_array_equal=bool(np.array_equal(g,go)),value_difference=v-vo,gradient_Linf_difference=float(np.max(abs(g-go)))))
 # Check shared-counter reassignment used by the unchanged inherited readout.
 shared=Counter();op.cost=shared
 if op.blocks is not None:raise AssertionError('Expected final no-hop control')
 sharing_setup=Counter()
 sharing=TwofoldDenominator(G,U,shape,one_hop=True,cost=sharing_setup);sharing.cost=shared
 if sharing.blocks.cost is not shared:raise AssertionError('Counter reassignment not shared')
 np.savez_compressed(out/'FIXED_WIRING_INPUTS.npz',J=J,Q=Q,Os=Os,U11=U,uniform=materials['uniform'],deterministic_nonuniform=materials['deterministic_nonuniform'],direction0=directions[0],direction1=directions[1])
 with (out/'DERIVATIVE_VALUES.csv').open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
 dump=dict(derivative_values=rows,nohop_original_equality=equality,transform_adjoint_pairings=adjoints)
 (out/'VALUES.json').write_text(json.dumps(dump,indent=2)+'\n')
 passed=all(r['relative_error']<2e-6 or r['absolute_error']<1e-10 for r in rows) and all(r['value_exact_equal'] and r['gradient_array_equal'] for r in equality) and max(r['relative_error'] for r in adjoints)<1e-11
 receipt.update(status='PASS' if passed else 'CHECK_FAILED',public_input=str(public),public_keys_read=['G','E0','grid_shape'],basis_source=str(spaces),basis_keys_read=['U11','G_right_frame','S_right_frame'],U11_rank=U.shape[1],
  auxiliary_inputs='J=.1 E0, Q=E0+G J, fixed throughout material derivatives; not physical PDE solution',
  metric='Os=diag(1/sum |E0|² per source), fixed; original Regularizer defaults smooth=tv=1e-4 epsilon=.01',
  beta=beta,rho=rho,epsilon=epsilon,legal_material_bounds=[.02,20.],material_ranges={k:[float(m.min()),float(m.max())] for k,m in materials.items()},
  directions='unit deterministic sin(.23*(index+1)), cos(.37*(index+1)); same two directions for both materials',
  derivative_max_relative_error=max(r['relative_error'] for r in rows),derivative_max_absolute_error=max(r['absolute_error'] for r in rows),
  nohop_original_bitwise_equal=all(r['value_exact_equal'] and r['gradient_array_equal'] for r in equality),
  adjoint_max_relative_error=max(r['relative_error'] for r in adjoints),operation_cost=dict(cost),operator_setups=operator_receipts,
  inherited_material_readout_identity=twofold_material_readout is matrix_material_readout,shared_counter_reassignment='PASS',
  cost_limits='Counts inherited small LU/triangular solves and block actions once in original counters; separate request/action counters are labels and must not be summed as extra operations. Shared-counter construction diagnostic setup is in separate count below.',
  counter_reassignment_setup_cost=dict(sharing_setup),BLAS_threads=1)
except BaseException as exc:
 receipt.update(status='FAILED',errors=[repr(exc)],traceback=traceback.format_exc());raise
finally:
 receipt.update(CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,whole_process_CPU_seconds_at_receipt=time.process_time())
 (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['status','U11_rank','derivative_max_relative_error','nohop_original_bitwise_equal','adjoint_max_relative_error','CPU_seconds','wall_seconds']},indent=2))
if receipt['status']!='PASS':raise SystemExit(1)
