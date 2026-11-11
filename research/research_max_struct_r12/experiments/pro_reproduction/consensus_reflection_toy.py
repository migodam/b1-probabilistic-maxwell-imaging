"""Exploratory scalar-TM finite VIE feasibility study, 2026-11-08.

Purpose: test exact cellwise current/field projection and reflected feasibility
iteration, not benchmark against mature CIE or claim continuum accuracy.
All four cases and algorithm parameters are fixed before running.
Requires numpy, scipy. Run with OPENBLAS_NUM_THREADS=1.
"""
from pathlib import Path
import json, time
import numpy as np
from scipy.special import hankel1, j1
from scipy.linalg import cho_factor, cho_solve, lu_factor, lu_solve, eigh
from scipy.optimize import minimize, brentq
from scipy.stats import chi2

OUT = Path(__file__).parent

def setup(n=10, ns=6, nr=24):
    k = 2*np.pi
    side = 1.0
    dx = side/n
    c = (np.arange(n)+.5)*dx-side/2
    xx, yy = np.meshgrid(c,c,indexing='xy')
    xy = np.stack([xx.ravel(), yy.ravel()],axis=1)
    a = dx/np.sqrt(np.pi)
    r = np.linalg.norm(xy[:,None,:]-xy[None,:,:],axis=-1)
    np.fill_diagonal(r,1.)
    factor = 1j*np.pi*k*a/2*j1(k*a)
    G = factor*hankel1(0,k*r)
    np.fill_diagonal(G,1j*np.pi*k*a/2*hankel1(1,k*a)-1)
    ts = np.arange(ns)*2*np.pi/ns
    dirs = np.stack([np.cos(ts),np.sin(ts)],axis=1)
    E = np.exp(1j*k*xy@dirs.T)
    tr = np.arange(nr)*2*np.pi/nr
    rx = 3*np.stack([np.cos(tr),np.sin(tr)],axis=1)
    S = factor*hankel1(0,k*np.linalg.norm(rx[:,None,:]-xy[None,:,:],axis=-1))
    return xy,G,E,S

class Physics:
    def __init__(self,G,E,S,Y):
        self.G,self.E,self.S,self.Y=G,E,S,Y
        self.norm=float(np.vdot(Y,Y).real)
        self.lu_count=0
        self.rhs_count=0
    def forward(self,m,grad=False):
        fac=lu_factor(np.eye(len(m))-m[:,None]*self.G)
        self.lu_count+=1
        J=lu_solve(fac,m[:,None]*self.E)
        self.rhs_count+=self.E.shape[1]
        Q=self.E+self.G@J
        r=self.S@J-self.Y
        f=.5*float(np.vdot(r,r).real)/self.norm
        if not grad:return f,J,Q
        adj=lu_solve(fac,self.S.conj().T@r,trans=2)
        self.rhs_count+=self.E.shape[1]
        g=np.real(np.sum(np.conj(adj)*Q,axis=1))/self.norm
        return f,g

class Projectors:
    def __init__(self,G,E,S,Y,delta,scale=6.,mmax=12.):
        self.G,self.E,self.S,self.Y=G,E,S,Y
        self.delta,self.scale,self.amax=delta,scale,mmax/scale
        self.B,self.T=scale*G,scale*S
        self.H=cho_factor(np.eye(len(G))+self.B.conj().T@self.B)
        self.V=cho_solve(self.H,self.T.conj().T)
        K=self.T@self.V
        self.val,self.U=eigh((K+K.conj().T)/2)
        self.val=np.maximum(self.val,0)
        self.countA,self.countC=0,0
    def affine_ball(self,X,Q):
        self.countA+=1
        xu=cho_solve(self.H,X+self.B.conj().T@(Q-self.E))
        r=self.T@xu-self.Y
        if np.linalg.norm(r)<=self.delta:return xu,self.E+self.B@xu
        z=self.U.conj().T@r
        power=np.sum(abs(z)**2,axis=1)
        def discrepancy(mu):
            return float(np.sum(power/(1+mu*self.val)**2)-self.delta**2)
        hi=1.
        while discrepancy(hi)>0:
            hi*=4
            if hi>1e30:raise RuntimeError('Data ball projection inconsistent or numerically rank deficient')
        mu=brentq(discrepancy,0,hi,xtol=1e-11,rtol=1e-12)
        rmu=self.U@(z/(1+mu*self.val)[:,None])
        x=xu-mu*self.V@rmu
        return x,self.E+self.B@x
    def material(self,X,Q):
        self.countC+=1
        aa=np.sum(abs(X)**2,axis=1)
        bb=np.sum(abs(Q)**2,axis=1)
        cc=np.real(np.sum(np.conj(X)*Q,axis=1))
        d=bb-aa
        disc=np.sqrt(d*d+4*cc*cc)
        root1=np.divide(-d+disc,2*cc,out=np.full_like(d,np.nan),where=abs(cc)>1e-25)
        root2=np.divide(-d-disc,2*cc,out=np.full_like(d,np.nan),where=abs(cc)>1e-25)
        candidates=np.stack([np.zeros_like(d),np.full_like(d,self.amax),root1,root2],axis=1)
        valid=np.isfinite(candidates)&(candidates>=0)&(candidates<=self.amax)
        clean=np.nan_to_num(candidates,nan=0,posinf=0,neginf=0)
        values=(aa[:,None]-2*clean*cc[:,None]+clean**2*bb[:,None])/(1+clean**2)
        values[~valid]=np.inf
        idx=np.argmin(values,axis=1)
        alpha=clean[np.arange(len(X)),idx]
        q=(Q+alpha[:,None]*X)/(1+alpha**2)[:,None]
        return alpha[:,None]*q,q,alpha*self.scale

def stats(m,truth,f):
    return {'relative_material_l2':float(np.linalg.norm(m-truth)/np.linalg.norm(truth)),
            'relative_data_l2':float(np.sqrt(2*f))}

def run_case(name,truth,G,E,S,seed):
    base=Physics(G,E,S,np.ones((len(S),E.shape[1]),complex))
    _,jt,qt=base.forward(truth)
    yclean=S@jt
    rng=np.random.default_rng(seed)
    sigma=.01*np.linalg.norm(yclean)/np.sqrt(yclean.size)
    Y=yclean+sigma/np.sqrt(2)*(rng.normal(size=yclean.shape)+1j*rng.normal(size=yclean.shape))
    delta=sigma*np.sqrt(chi2.ppf(.99,2*Y.size)/2)
    m0=np.full(len(truth),.1)
    record={'case':name,'sigma_complex_rms':float(sigma),'delta_CN99':float(delta),
            'relative_noise_bound':float(delta/np.linalg.norm(Y)),'methods':{}}
    arrays={'truth':truth,'Y':Y,'yclean':yclean}
    p=Physics(G,E,S,Y)
    t=time.perf_counter()
    sol=minimize(lambda m:p.forward(m,True),m0,jac=True,method='L-BFGS-B',bounds=[(0,12)]*len(truth),
                 options={'maxiter':600,'maxls':30,'ftol':1e-14,'gtol':1e-10})
    f,j,q=p.forward(sol.x)
    record['methods']['direct']={**stats(sol.x,truth,f),'seconds':time.perf_counter()-t,
            'lu_count':p.lu_count,'rhs_count':p.rhs_count,'iterations':int(sol.nit),'status':str(sol.message)}
    arrays['direct']=sol.x
    for mode in ['alternating_projection','douglas_rachford']:
        p=Physics(G,E,S,Y)
        t=time.perf_counter()
        pr=Projectors(G,E,S,Y,delta)
        f,J,Q=p.forward(m0)
        X=J/pr.scale
        best=(f,m0.copy(),0)
        hist=[]
        for it in range(1,3001):
            cx,cq,cm=pr.material(X,Q)
            if mode=='douglas_rachford':
                ax,aq=pr.affine_ball(2*cx-X,2*cq-Q)
                gap=np.hypot(np.linalg.norm(ax-cx),np.linalg.norm(aq-cq))
                X=X+ax-cx;Q=Q+aq-cq
            else:
                ax,aq=pr.affine_ball(cx,cq)
                gap=np.hypot(np.linalg.norm(ax-cx),np.linalg.norm(aq-cq))
                X,Q=ax,aq
            if it%50==0:
                ff,_,_=p.forward(cm)
                if ff<best[0]:best=(ff,cm.copy(),it)
                hist.append([it,ff,float(gap),float(np.linalg.norm(cm-truth)/np.linalg.norm(truth))])
        elapsed=time.perf_counter()-t
        ff,mm,at=best
        record['methods'][mode]={**stats(mm,truth,ff),'seconds':elapsed,
               'lu_count':p.lu_count,'rhs_count':p.rhs_count,'projection_A':pr.countA,'projection_C':pr.countC,
               'selected_iteration_by_data':at,'last_shadow':stats(cm,truth,hist[-1][1]),'last_gap':hist[-1][2]}
        arrays[mode]=mm
        arrays[mode+'_history']=np.asarray(hist)
        # Same 200-step physical polishing control, reported separately and charged.
        before_lu,before_rhs=p.lu_count,p.rhs_count
        tp=time.perf_counter()
        polished=minimize(lambda m:p.forward(m,True),mm,jac=True,method='L-BFGS-B',bounds=[(0,12)]*len(truth),
                options={'maxiter':200,'maxls':30,'ftol':1e-14,'gtol':1e-10})
        fpol,_,_=p.forward(polished.x)
        record['methods'][mode+'_plus_200_direct']={**stats(polished.x,truth,fpol),
                'extra_seconds':time.perf_counter()-tp,'extra_lu':p.lu_count-before_lu,
                'extra_rhs':p.rhs_count-before_rhs,'iterations':int(polished.nit)}
        arrays[mode+'_polished']=polished.x
    np.savez_compressed(OUT/(name+'.npz'),**arrays)
    return record

def main():
    xy,G,E,S=setup()
    x,y=xy[:,0],xy[:,1]
    rad=np.hypot(x,y)
    ring=(rad>.19)&(rad<.36)
    shell=ring&~((x>.0)&(abs(y)<.11))
    bridge=((abs(x)>.19)&(abs(x)<.36)&(abs(y)<.31))|((abs(x)<.31)&(abs(y)<.065))
    disk=rad<.29
    cases=[('ring_high',8.*ring),('split_shell_high',8.*shell),('bridge_high',8.*bridge),('disk_low',2.*disk)]
    all_results={'scope':'Exploratory fixed finite scalar TM model; not CIE comparison, no equal-cost performance claim.',
            'parameters':{'n':10,'ns':6,'nr':24,'side_wavelength':1.,'box':[0,12],'scale':6.,'projection_steps':3000,
                          'physical_audit_every':50,'noise_relative_rms':.01,'initial_material':.1,'prior':'none'},'cases':[]}
    for i,(name,truth) in enumerate(cases):
        result=run_case(name,truth,G,E,S,7301+i)
        all_results['cases'].append(result)
        print(json.dumps(result,ensure_ascii=False),flush=True)
        (OUT/'toy_results.json').write_text(json.dumps(all_results,indent=2,ensure_ascii=False))
    print('DONE',flush=True)

if __name__=='__main__':main()
