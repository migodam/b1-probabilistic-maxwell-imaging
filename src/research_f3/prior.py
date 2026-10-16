import copy
import numpy as np
from prior_constrained import prior

LIMIT=float(np.log(9/4))
SCALE=.00527555

class QuotientPrior:
    """q=(six original uniform shape coordinates,t); rho=(2t-1)log(9/4).

    Common weight scale is integrated out, not assigned a uniform gauge-fixed
    prior. Gaussian-label interchange remains a discrete decoder symmetry.
    """
    lo=np.r_[np.full(4,-.13),np.full(2,.05)]
    hi=np.r_[np.full(4,.13),np.full(2,.095)]
    def __init__(self,meta):self.meta=meta;self.row=copy.deepcopy(meta['generator']);self.bg=meta['background_id']
    def from_eight(self,p):
        w=.4+.5*np.asarray(p)[6:];rho=np.log(w[0]/w[1]);return np.r_[p[:6],(rho+LIMIT)/(2*LIMIT)]
    def decode(self,q):
        q=np.asarray(q);r=copy.deepcopy(self.row)
        if q.shape!=(7,) or np.any(q<0) or np.any(q>1):raise ValueError('outside original quotient support')
        x=self.lo+(self.hi-self.lo)*q[:6];rho=(2*q[6]-1)*LIMIT
        w=np.clip(np.array([.6*np.exp(rho/2),.6*np.exp(-rho/2)]),.4,.9)
        if q[6]==0:w=np.array([.4,.9])
        if q[6]==1:w=np.array([.9,.4])
        d=r['geometry'];d['gaussian_centers'][:2]=x[:4].reshape(2,2).tolist();d['gaussian_widths'][:2]=x[4:].tolist();d['gaussian_weights'][:2]=w.tolist()
        return prior.decode(r,background_id=self.bg),r
    def logprior(self,q):
        q=np.asarray(q)
        if np.any(q<0) or np.any(q>1):return -np.inf
        gap=2*LIMIT*min(q[6],1-q[6]);r=LIMIT-gap
        if gap<=0:return -np.inf
        return float(np.log(2*LIMIT)+2*np.log(.9)-r+np.log(-np.expm1(-2*gap))-np.log(2*.5**2))
    def derivative(self,q):
        _,row=self.decode(q);d=row['geometry'];p=prior.geometry_points()-np.array(d['center']);co,si=np.cos(d['rotation']),np.sin(d['rotation'])
        xx=(co*p[:,0]+si*p[:,1])/d['scale'];yy=(-si*p[:,0]+co*p[:,1])/(d['scale']*d['axis_ratio'])
        vals=[];D=np.zeros((len(p),7))
        for i in range(2):
            gx,gy=d['gaussian_centers'][i];w=d['gaussian_weights'][i];s=d['gaussian_widths'][i];dx,dy=xx-gx,yy-gy;dist=dx*dx+dy*dy
            z=w*np.exp(-dist/(2*s*s));vals.append(z)
            D[:,2*i]=z*dx/s**2*.26;D[:,2*i+1]=z*dy/s**2*.26;D[:,4+i]=z*dist/s**3*.045
            D[:,6]+=z*LIMIT*(1 if i==0 else -1)
        z=np.sum(vals,axis=0);j=np.argmax(z);return 1.2*prior.V*(D/z[j]-z[:,None]*D[j]/z[j]**2)

def reflected(x):
    """Coordinate reflection. Symmetric for a fixed diagonal Gaussian only."""
    z=np.mod(x,2);return np.where(z<=1,z,2-z)

def finite_number(x):return float(x) if np.isfinite(x) else None
