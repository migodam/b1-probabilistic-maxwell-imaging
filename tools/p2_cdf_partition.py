"""Positive disjoint-interval integration of the SAME continuous density.

No probability is spread uniformly inside bins. Only CDF grid values are
computed by adding independently integrated adjacent interval masses. The
normalizer of each partition is compared with the original tensor integral.
Directly inverted continuous quantiles remain a separate saved calculation.
"""
from pathlib import Path
import sys,json
import numpy as np
from scipy.special import logsumexp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.b1_joint.reference import Reference
from src.b1_joint.response import Response
from src.b1_joint.material import RMAX,SCALES,event_intervals,BOUNDS
from src.b1_joint.likelihood import log_likelihood,log_cdf_difference
from src.b1_joint.accounting import dump

def batch(F,y,sd):
    f=F/sd;ff=np.sum(f*f,axis=1);yy=y/sd;dot=f@yy.T;den=1+.15**2*ff
    q=np.sum(yy*yy,axis=1)[None,:]-2*dot+ff[:,None]-.15**2*(dot-ff[:,None])**2/den[:,None]
    ll=-.5*(q+(2*np.log(sd).sum()+np.log(den))[:,None]+32*np.log(2*np.pi))
    v=1/(ff+1/.15**2);m=v[:,None]*(dot+1/.15**2);s=np.sqrt(v)[:,None]
    return ll+log_cdf_difference((.7-m)/s,(1.3-m)/s)-log_cdf_difference(-2.,2.)

def integrate(ref,y,sd,grid,order):
    n=len(y);intervals=event_intervals();cdf=np.zeros((n,4,65));cond=np.zeros((n,12,4,65));norm=np.zeros((n,4));event=np.zeros((n,4,12))
    for axis in range(4):
        logs=np.full((65,n,12),-np.inf)
        for event_id in range(12):
            c,s=divmod(event_id,3);box=[intervals[0][c],intervals[1][s],tuple(BOUNDS[2]),tuple(BOUNDS[3])]
            if axis<2 and box[axis] is None:
                cell=ref.cell(event_id,order,box_override=box);logs[0,:,event_id]=logsumexp(batch(cell['F1'],y,sd)+np.log(cell['prior'])[:,None],axis=0);continue
            for ti in range(1,65):
                lo,hi=grid[axis,ti-1:ti+1]
                if axis<2:lo,hi=-(1-(hi/SCALES[axis])**2)**-2,-(1-(lo/SCALES[axis])**2)**-2
                lo=max(lo,box[axis][0]);hi=min(hi,box[axis][1])
                if hi<=lo:continue
                part=list(box);part[axis]=(lo,hi);cell=ref.cell(event_id,order,box_override=part)
                logs[ti,:,event_id]=logsumexp(batch(cell['F1'],y,sd)+np.log(cell['prior'])[:,None],axis=0)
        den=logsumexp(logs,axis=(0,2));ed=logsumexp(logs,axis=0);norm[:,axis]=den;event[:,axis]=np.exp(ed-den[:,None])
        masses=np.exp(logs-den[None,:,None]);cdf[:,axis]=np.cumsum(masses.sum(2),axis=0).T
        cond[:,:,axis]=np.cumsum(np.exp(logs-ed[None,:,:]),axis=0).transpose(1,2,0)
        print(f'positive CDF order{order} axis{axis}',flush=True)
    return cdf,cond,norm,event

def main():
    out=ROOT/'research/joint_p2/results';old=np.load(out/'FRESH_CONTINUOUS_REFERENCE.npz');data=np.load(out/'FRESH_FINAL_OBSERVATIONS.npz');y=data['y'][:,:32];sd=data['sigma'][:32];ref=Reference(Response(out/'COMMON_RESPONSE_BANK.npz'));grid=old['grid'];a,b,za,ea=integrate(ref,y,sd,grid,13);fine,fb,zf,ef=integrate(ref,y[:8],sd,grid,17)
    base=np.array([ref.solve(v,sd,order=13)['log_evidence'] for v in y]);arrays={k:old[k] for k in old.files};arrays.update(CDF=a,conditional_CDF=b,partition_log_evidence=za,partition_event_mass=ea)
    np.savez_compressed(out/'FRESH_CONTINUOUS_REFERENCE_PARTITION.npz',**arrays)
    audit=json.loads((out/'FRESH_CONTINUOUS_AUDIT.json').read_text());audit.update(status='POSITIVE_PARTITION_INTEGRAL_REFERENCE',old_truncated_monotonicity_violation=audit['monotone_violation'],monotone_violation=float(max(0.,-np.diff(a,axis=-1).min())),important_conditional_monotone_violation=float(max(0.,-np.diff(b,axis=-1).min())),upper_error=float(abs(a[:,:,-1]-1).max()),max_partition_log_evidence_disagreement=float(abs(za-base[:,None]).max()),max_partition_event_mass_disagreement=float(abs(ea-old['event_mass'][:,None,:]).max()),fixed_first_eight_partition_CDF_13_17_difference=float(abs(a[:8]-fine).max()),fixed_first_eight_partition_conditional_CDF_13_17_difference=float(abs(b[:8]-fb).max()),max_old_new_CDF_difference=float(abs(a-old['CDF']).max()),new_measure='same continuous coefficient-prior likelihood; positive quadrature on adjacent disjoint intervals; no within-bin uniform density or bin quantile',normalization_scope='whole support per integration partition; different numerical approximations to same exact evidence; discrepancies explicitly shown',first_eight_partition_fine_CDF=fine.tolist(),reference_quantile_scope='first8 direct continuous inversion unchanged; no quantile interpolation from grid')
    dump(out/'FRESH_CONTINUOUS_PARTITION_AUDIT.json',audit)
if __name__=='__main__':main()
