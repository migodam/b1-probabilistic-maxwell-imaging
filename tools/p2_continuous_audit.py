"""Direct fresh-panel CDF integrals, without physics or node-bin broadening.

This retains the frozen model evaluation. Its original resampled-node CDFs
remain diagnostics; this file supplies a genuinely continuous reference.
"""
from pathlib import Path
import sys,time
import numpy as np
from scipy.special import logsumexp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.b1_joint.reference import Reference
from src.b1_joint.response import Response
from src.b1_joint.material import RMAX
from src.b1_joint.likelihood import log_likelihood,log_cdf_difference
from src.b1_joint.accounting import dump

def main():
    out=ROOT/'research/joint_p2/results';dest=out/'FRESH_CONTINUOUS_REFERENCE.npz'
    if dest.exists():raise RuntimeError('completed continuous audit is retained')
    data=np.load(out/'FRESH_FINAL_OBSERVATIONS.npz');y=data['y'][:,:32];sd=data['sigma'][:32]
    ref=Reference(Response(out/'COMMON_RESPONSE_BANK.npz'));n=len(y)
    results=[ref.solve(a,sd,order=13) for a in y]
    fine=[ref.solve(a,sd,order=17) for a in y[:8]]
    bounds=np.array([[0,RMAX[0]],[0,RMAX[1]],[2.2,2.8],[.04,.16]])
    grid=np.array([np.linspace(a,b,65) for a,b in bounds]);cdf=np.empty((n,4,65));conditional=np.empty((n,12,4,65))
    # Shared cells depend on the target and threshold, never observation/truth.
    for axis in range(4):
        for ti,cut in enumerate(grid[axis]):
            logs=np.full((n,12),-np.inf)
            for event in range(12):
                cell=ref.cell(event,13,cut=(axis,float(cut)))
                if cell is None:continue
                lp=np.log(cell['prior'])
                # Fixed noise permits shared Gaussian sufficient statistics;
                # this is the exact same bounded-gain marginal, batched over y.
                f=cell['F1']/sd;ff=np.sum(f*f,axis=1);yy=y/sd;dot=f@yy.T;den=1+.15**2*ff
                q=np.sum(yy*yy,axis=1)[None,:]-2*dot+ff[:,None]-.15**2*(dot-ff[:,None])**2/den[:,None]
                ll=-.5*(q+(2*np.log(sd).sum()+np.log(den))[:,None]+32*np.log(2*np.pi))
                v=1/(ff+1/.15**2);m=v[:,None]*(dot+1/.15**2);ss=np.sqrt(v)[:,None]
                ll+=log_cdf_difference((.7-m)/ss,(1.3-m)/ss)-log_cdf_difference(-2.,2.)
                # The first cell is independently compared with the scalar law.
                if axis==0 and ti==0 and event==0:
                    if abs(ll[:,0]-log_likelihood(cell['F1'],y[0],sd)).max()>1e-9:raise FloatingPointError('batched likelihood differs')
                logs[:,event]=logsumexp(ll+lp[:,None],axis=0)
            den=np.array([r['log_evidence'] for r in results]);ep=np.array([r['log_event_mass'] for r in results])
            cdf[:,axis,ti]=np.exp(logsumexp(logs,axis=1)-den)
            conditional[:,:,axis,ti]=np.exp(logs-den[:,None]-ep)
        print(f'continuous CDF axis {axis} complete',flush=True)
    # Direct inversion, not inversion of a binned/node CDF, on a fixed subset.
    summaries=[];checks=[]
    for i in range(8):
        a=ref.continuous_summary(results[i],grid_size=9)
        b=ref.continuous_summary(fine[i],grid_size=9)
        checks.append(dict(observation=i,CDF_max=float(abs(np.array(a[0]['CDF'])-b[0]['CDF']).max()),quantile_fraction=float((abs(np.array(a[0]['quantiles_05_50_95'])-b[0]['quantiles_05_50_95'])/(bounds[:,1]-bounds[:,0])[:,None]).max()),coherence=bool(a[0]['coherence_pass'] and b[0]['coherence_pass'])))
        summaries.append(b[0])
    np.savez_compressed(dest,grid=grid,CDF=cdf,conditional_CDF=conditional,mean=np.array([r['mean'] for r in results]),SD=np.array([r['SD'] for r in results]),event_mass=np.array([r['event_mass'] for r in results]),parent=data['parent'],target_identity=np.array('joint-p2-radial-v1:C:tensor-Chebyshev-five-per-physical-axis-radial-v1'))
    dump(out/'FRESH_CONTINUOUS_AUDIT.json',dict(observations=n,order=13,grid_points=65,measure='original continuous coefficient prior with geometry pushforward atoms; truncated integrals, no bin widening',old_sample_CDF='preserved as node-resampling diagnostic, not continuous reference',direct_refinement_checks=checks,first_eight_direct_quantile_summaries=summaries,monotone_violation=float(max(0.,-np.diff(cdf,axis=-1).min())),upper_error=float(abs(cdf[:,:,-1]-1).max()),new_physics=0,limitations='grid max CDF is not uniform supremum; 13/17 refinement subset empirical; no continuum/real-world calibration'))
if __name__=='__main__':main()
