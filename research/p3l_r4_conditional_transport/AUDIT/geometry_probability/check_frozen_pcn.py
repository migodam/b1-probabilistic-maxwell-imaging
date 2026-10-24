"""Pure frozen-matrix prior-reversibility and material-RMS normalization checks."""
from pathlib import Path
import json,time,resource,sys
import numpy as np
from p3r4.frozen_geometry_pcn import *
wall,cpu=time.monotonic(),time.process_time();rng=np.random.default_rng(20261025112)
root=Path('research/p3l_r4_conditional_transport/AUDIT/geometry_probability/PREPARED_V1_VALIDATED');protocol=json.loads((root/'FROZEN_PROTOCOL.json').read_text());results=[];maxerror=0.
for row in protocol['rows']:
 with np.load(row['teacher_file']['path']) as z:D=real_decoder_metric(z['Dchi']);x0=z['x0']
 for name in KERNELS:
  with np.load(row['kernels'][name]['file']) as z:kernel=FrozenMatrixPCN(z['Q'])
  assert kernel.values.min()>0 and kernel.values.max()<=Q_CAP*(1+1e-10)
  actual=kernel.material_linear_expected_rms(D,x0);assert abs(actual-RMS_TARGET)<1e-10
  covariance_error=float(np.max(np.abs(kernel.rho@kernel.rho.T+kernel.Q-np.eye(66))));assert covariance_error<1e-12
  for t in range(10):
   x=rng.normal(size=66);candidate=kernel.propose(x,rng.normal(size=66))
   ratio=kernel.log_proposal(x,candidate)-kernel.log_proposal(candidate,x)
   expected=.5*(candidate@candidate-x@x);error=abs(ratio-expected);maxerror=max(maxerror,error);assert error<1e-7
   phi0,phi1=rng.uniform(0,100,size=2);explicit=-phi1+phi0-.5*(candidate@candidate-x@x)+ratio
   assert abs(explicit-(phi0-phi1))<1e-7
  results.append(dict(parent=row['parent'],kernel=name,min_eigen=float(kernel.values.min()),max_eigen=float(kernel.values.max()),expected_RMS=actual,covariance_identity_max_error=covariance_error))
# Explicit cap binding must retain the achieved scale and reject matched claim.
Dsmall=np.eye(1152,66)*1e-8;_,capped=normalized_kernel(np.eye(66),Dsmall,np.zeros(66))
assert capped['cap_binds'] and not capped['matched_RMS'] and abs(capped['lambda_max']-Q_CAP)<1e-12
# Degenerate/rank-deficient proposals may never masquerade as full66 kernels.
try:FrozenMatrixPCN(np.diag(np.r_[np.ones(65)*.001,0.]))
except ValueError:pass
else:raise AssertionError('rank deficient Q accepted')
result=dict(status='PASS',prepared_kernels=len(results),reversibility_pairs=320,max_logproposal_prior_cancellation_error=maxerror,
 all66_eigenvalues_positive=True,all_material_RMS_matched_in_actual_preparation=True,cap_binding_disclosure='PASS',rank_deficiency_rejected='PASS',rows=results,
 physics_RHS=0,physics_states=0,training_updates=0,wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,
 peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024))
Path('research/p3l_r4_conditional_transport/AUDIT/geometry_probability/KERNEL_MATH_CHECK.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='rows'}))
