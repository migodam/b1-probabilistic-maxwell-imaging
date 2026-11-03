"""One spectral evaluation of the newly derived hard-lock noise distribution.
No Monte Carlo noise sweep, no cutoff tuning; GT is disclosed R2 diagnostic only.
"""
from pathlib import Path
import sys,json,time
from datetime import datetime,timezone
import numpy as np
from scipy.linalg import svd,lu_solve,cho_factor,cho_solve,eigvalsh
from scipy.stats import chi2
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel

def main():
    tick,cpu=time.perf_counter(),time.process_time();case='obs_b20d5a55677d4eca';r2=REPO/'research/research_max_struct_r2'
    dest=ROOT/'results/diagnostics/NOISE_FLOOR_DISTRIBUTION.json'
    if dest.exists():raise FileExistsError('preserve realization')
    with np.load(r2/'data_public'/f'{case}.npz') as z:public={k:z[k]for k in z.files}
    with np.load(r2/'sealed_truth'/f'{case}.npz') as z:m=z['m'];noise=z['noise']
    model=MaxwellModel(public);st=model.state(m);U,s,Vh=svd(model.S,full_matrices=False)
    r=np.count_nonzero(s>1e-3*s[0]);D=Vh[:r];P=1/(1+m)
    with model._physics.operation('DPCS_NOISE_TRANSFER',rhs=int(r),receiver_transpose_RHS=int(r)):
        A=lu_solve(st['lu'],D.T,trans=1,check_finite=False).T/P
    Gram=A@A.conj().T;fac=cho_factor(Gram,lower=True)
    invsig=np.diag(1/s[:r]);C=float(public['sigma'])**2*invsig@cho_solve(fac,invsig)
    ce=eigvalsh(C);ws=1/np.sum(abs(model.E)**2,axis=0);be=(ws[:,None]*ce[None,:]/2).ravel()
    expected=float(be.sum());t=np.log(100);bound=float(expected+np.sqrt(2*np.sum(be**2)*t)+be.max()*t)
    b=(U[:,:r].conj().T@noise)/s[:r,None]
    observed=float(.5*np.sum(ws*np.real(np.sum(b.conj()*cho_solve(fac,b),axis=0))))
    rec=dict(status='COMPLETED',classification='OFFLINE_DIAGNOSTIC_PUBLIC_GT',case_id=case,truth_access='PUBLISHED_R2_ONLY',
        retained_rank=int(r),complex_noise_sd=float(public['sigma']),covariance_eigenvalues=ce.tolist(),
        expected_hard_lock_state_term=expected,observed_hard_lock_state_term=observed,upper_probability_99=bound,
        retained_whitened_noise_norm=float(np.linalg.norm(U[:,:r].conj().T@noise)/float(public['sigma'])),
        retained_noise_radius99=float(np.sqrt(chi2.ppf(.99,2*r*model.ns)/2)),
        physical_cost=dict(model.cost),auxiliary_cost=dict(receiver_SVD=1,Gram_Cholesky=1,Gram_RHS=int(r+model.ns),Gram_eigenproblem_dimension=int(r)),
        validation='Distribution conditional on fixed GT,H,modes; no claim for data-adapted selected model or continuous PDE',
        wall_seconds=time.perf_counter()-tick,cpu_seconds=time.process_time()-cpu,finished_utc=datetime.now(timezone.utc).isoformat())
    dest.write_text(json.dumps(rec,indent=2));print(json.dumps({k:rec[k]for k in ['expected_hard_lock_state_term','observed_hard_lock_state_term','upper_probability_99','retained_whitened_noise_norm','retained_noise_radius99','physical_cost','wall_seconds']}))
if __name__=='__main__':
    with threadpool_limits(limits=1):main()
