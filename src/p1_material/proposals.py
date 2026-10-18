"""Material-only controls. Retraction outputs candidates, never posterior samples."""
import numpy as np
from scipy.optimize import brentq,minimize
from .prior import structural_change


def normalized_direction(prior,q,rng):
    v=rng.normal(size=66);v[64:]*=.25
    for i in range(64,66):
        if (q[i]<=0 and v[i]<0) or (q[i]>=1 and v[i]>0):v[i]*=-1
    w=prior.derivative(q)@v;n=np.linalg.norm(w)/24
    if n==0:raise ValueError('no material velocity')
    return v/n,w/n


def frozen_latent_step(prior,q,v,target):
    x0,_=prior.decode(q);cap=np.inf
    for k in [64,65]:
        if v[k]>0:cap=min(cap,(1-q[k])/v[k])
        if v[k]<0:cap=min(cap,-q[k]/v[k])
    hi=min(float(target),cap*(1-1e-12));lo=0
    if hi<=0:return q.copy(),False
    f=lambda t:prior.material_rms(prior.decode(q+t*v)[0],x0)-target
    for _ in range(20):
        if f(hi)>=0:
            step=brentq(f,lo,hi,xtol=1e-12);return q+step*v,True
        lo=hi;hi=min(2*hi,cap*(1-1e-12))
        if hi<=lo:break
    return q+hi*v,False


def material_retraction(prior,q0,v,w,target):
    x0,_=prior.decode(q0);goal=x0+target*w;seed,matched=frozen_latent_step(prior,q0,v,target)
    def objective(q):
        d=prior.decode(q)[0]-goal;dq=q-q0
        return float(np.vdot(d,d).real/576+1e-8*np.dot(dq,dq))
    def jac(q):return 2*np.real(prior.derivative(q).conj().T@(prior.decode(q)[0]-goal))/576+2e-8*(q-q0)
    def sphere(q):d=prior.decode(q)[0]-x0;return float(np.vdot(d,d).real/576-target*target)
    def sphere_jac(q):return 2*np.real(prior.derivative(q).conj().T@(prior.decode(q)[0]-x0))/576
    def progress(q):return float(np.vdot(w,prior.decode(q)[0]-x0).real/576-.5*target)
    def progress_jac(q):return np.real(prior.derivative(q).conj().T@w)/576
    fit=minimize(objective,seed,jac=jac,bounds=[(None,None)]*64+[(0,1)]*2,method='SLSQP',constraints=[dict(type='eq',fun=sphere,jac=sphere_jac),dict(type='ineq',fun=progress,jac=progress_jac)],options=dict(ftol=1e-11,maxiter=150))
    q=fit.x;x,_=prior.decode(q);rms=prior.material_rms(x,x0)
    return q,dict(optimizer_success=bool(fit.success),message=str(fit.message),iterations=int(fit.nit),actual_RMS=rms,target_matched=abs(rms-target)<=.02*target,progress=float(progress(q)+.5*target),structural_change=structural_change(prior,q0,q),candidate_only=True,full_Maxwell_validated=False)
