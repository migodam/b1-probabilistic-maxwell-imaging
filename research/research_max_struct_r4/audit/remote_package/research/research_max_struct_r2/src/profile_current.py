"""Known variable projection applied to full-current CIE; not a new theorem.

Uniform source metrics use sums normalized by total data/incident energy.
Exact current elimination, including derivative of beta/(1+beta*m).
"""
import time
from collections import Counter
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize

def run_profile(model,Y,m0,maxiter=250,beta=1.,gradient_check=False):
    start=time.perf_counter(); cpu=time.process_time(); cost=Counter(); history=[]
    G,E,S=model.G,model.E,model.S
    yd=float(np.vdot(Y,Y).real); ed=float(np.vdot(E,E).real)
    SS=S.conj().T@S/yd; SY=S.conj().T@Y/yd
    cost['S_Gram_builds']=1
    def objective(m,details=False):
        w=beta/(1+beta*m); wp=-w*w
        L=np.eye(model.n)-m[:,None]*G
        H=w[:,None]*L; b=(w*m)[:,None]*E
        A=SS+H.conj().T@H/ed
        fac=cho_factor(A,lower=True,check_finite=False)
        J=cho_solve(fac,SY+H.conj().T@b/ed,check_finite=False)
        cost['current_Cholesky']+=1; cost['current_RHS']+=model.ns
        cost['dense_Gram_builds']+=1; cost['G_columns']+=model.ns
        Q=E+G@J; e=J-m[:,None]*Q; r=w[:,None]*e; d=S@J-Y
        rv,rg,_=model.reg.evaluate(m)
        value=.5*np.vdot(d,d).real/yd+.5*np.vdot(r,r).real/ed+rv
        grad=np.real(np.sum(r.conj()*(wp[:,None]*e-w[:,None]*Q),axis=1))/ed+rg
        history.append(dict(auxiliary=float(value),data=float(.5*np.vdot(d,d).real/yd),
                            state=float(.5*np.vdot(r,r).real/ed),wall=time.perf_counter()-start))
        return (float(value),grad,J) if details else (float(value),grad)
    if gradient_check:
        h=np.cos(np.arange(model.n)+.37); h/=np.linalg.norm(h); step=1e-5
        value,g=objective(m0)
        finite=(objective(m0+step*h)[0]-objective(m0-step*h)[0])/(2*step)
        return dict(analytic=float(g@h),finite=float(finite),absolute_error=float(abs(g@h-finite)),
                    cost=dict(cost),wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu)
    opt=minimize(objective,m0,jac=True,method='L-BFGS-B',bounds=[(.02,20.)]*model.n,
                 options=dict(maxiter=maxiter,maxfun=4*maxiter,ftol=1e-12,gtol=1e-8,maxls=30))
    value,grad,J=objective(opt.x,True)
    cost.update(wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu)
    return dict(m=opt.x,J=J,history=history,cost=dict(cost),stop=str(opt.message),
                nit=int(opt.nit),auxiliary=value,gradient_norm=float(np.linalg.norm(grad)))
