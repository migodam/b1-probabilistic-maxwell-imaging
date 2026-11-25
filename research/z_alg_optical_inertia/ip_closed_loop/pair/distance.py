"""Distance dual for two conditional complex output quadrics."""
import numpy as np
from scipy.linalg import eigh, cho_factor, cho_solve, block_diag
from scipy.optimize import minimize

def herm(a):
    return (a+a.conj().T)*0.5

def translate(q, input_basis, reference):
    n, p = input_basis.shape[0], reference.shape[0]
    transform = np.block([[input_basis,np.zeros((n,p))],
                          [reference@input_basis,np.eye(p)]])
    return herm(transform.conj().T@q@transform)

def solve_weights(q0,q1,source_count,weights):
    l,h = weights
    n = source_count
    identity = np.eye(len(q0)-n)
    r = np.block([[identity+l*q0[n:,n:],-identity],[-identity,identity+h*q1[n:,n:]]])
    b = np.vstack((l*q0[n:,:n],h*q1[n:,:n]))
    cf=cho_factor(herm(r),check_finite=False)
    return herm(l*q0[:n,:n]+h*q1[:n,:n]-b.conj().T@cho_solve(cf,b,check_finite=False))

def optimize(q0,q1,n):
    def evaluate(t, vectors=False):
        l=np.exp(np.asarray(t))
        try:
            h=solve_weights(q0,q1,n,l)
            if vectors:
                e,v=eigh(h,subset_by_index=(n-1,n-1),check_finite=False)
                return float(e[0]),v[:,-1]
            return float(eigh(h,subset_by_index=(n-1,n-1),eigvals_only=True,check_finite=False)[0])
        except np.linalg.LinAlgError:
            return (-np.inf,None) if vectors else -np.inf
    starts = [(j+k,j-k) for j in (-11.,-6.,-2.,1.,4.) for k in (-1.75,0.,1.75)]
    ranked=sorted(starts,key=evaluate,reverse=True)
    best_t,best=ranked[0],evaluate(ranked[0])
    if not np.isfinite(best):
        return 0.,None,None
    scale=max(abs(best),1.e-14)
    def cost(v):
        score=evaluate(v)
        return -score/scale if np.isfinite(score) else 1.e8
    for t in ranked[:2]:
        ans=minimize(cost,t,method="Nelder-Mead",bounds=((-23.,14.),(-23.,14.)),
                     options={"maxiter":80,"xatol":.006,"fatol":.001})
        v=evaluate(ans.x)
        if v>best:
            best,best_t=v,ans.x
    score,vec=evaluate(best_t,True)
    return score,vec,np.exp(best_t)

def direction_support(q,a,w):
    n=len(a)
    r=q[n:,n:]
    try:
        cf=cho_factor(r,check_finite=False)
        center=-cho_solve(cf,q[n:,:n]@a,check_finite=False)
        gap=-np.vdot(a,q[:n,:n]@a).real-np.vdot(center,q[n:,:n]@a).real
        reach=np.sqrt(max(0.,gap)*np.vdot(w,cho_solve(cf,w,check_finite=False)).real)
        dot=np.vdot(w,center).real
        return dot-reach,dot+reach
    except np.linalg.LinAlgError:
        return -np.inf,np.inf

def observed_lower(q,a,y,radius):
    """Lower bound for min_{||e||<=radius} q(a,y+e), for fixed dual scalar."""
    n=len(a); R=q[n:,n:]; b=q[n:,:n]@a+R@y
    z=np.r_[a,y]; value=float(np.vdot(z,q@z).real)
    ev,V=eigh(R,check_finite=False); coeff=abs(V.conj().T@b)**2
    lo=max(0.,-float(ev[0]))+1.e-12
    def deriv(t):
        return float(np.sum(coeff/(ev+t)**2)-radius**2)
    if deriv(lo)>0:
        from scipy.optimize import brentq
        hi=max(1.,2*lo)
        while deriv(hi)>0:hi*=2
        multiplier=brentq(deriv,lo,hi)
    else:multiplier=lo
    return value-multiplier*radius**2-float(np.sum(coeff/(ev+multiplier))),float(multiplier)
