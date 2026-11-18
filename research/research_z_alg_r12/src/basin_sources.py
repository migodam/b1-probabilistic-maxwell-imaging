"""Stationary box faces under coherent reweighting of previously measured transmitters."""
import argparse
import json
import os
from pathlib import Path
from time import perf_counter, process_time

for _key in ("OMP_NUM_THREADS","VECLIB_MAXIMUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np
from scipy.linalg import eigh, lu_factor, lu_solve
from scipy.optimize import minimize
from scipy.spatial.distance import cdist
from scipy.special import j0

from high_endpoints import load_endpoint, physics, bounds, noise_radius

OUT = Path(__file__).resolve().parents[1] / "results"
SEED = 121407


def hermitian(a):
    return (a + a.conj().swapaxes(-1,-2)) / 2


class Wave:
    def __init__(self, public):
        self.G, self.B, self.S, self.Y = (public[key] for key in ("G","E0","S","Y"))
        self.n, self.s = self.B.shape
        self.den = float(np.vdot(self.Y,self.Y).real)
        self.gram = j0(float(public['k'])*cdist(public['sources'],public['sources']))
        ev, vec = eigh(self.gram)
        self.white = (vec / np.sqrt(ev)) @ vec.T
        self.reg = physics(public)[1]
        self.lo, self.hi = bounds(public)
        self.lu = self.rhs = 0

    def state(self, m, with_gradient=True, with_kernel=False):
        # E = B+GJ, J=mE. P propagates a local current source to each receiver.
        lu = lu_factor(np.eye(self.n)-m[:,None]*self.G,check_finite=False)
        j = lu_solve(lu,m[:,None]*self.B,check_finite=False)
        self.lu += 1
        self.rhs += self.s
        E = self.B + self.G @ j
        res = self.S @ j-self.Y
        P = K = None
        if with_gradient:
            P = lu_solve(lu,self.S.T,trans=1,check_finite=False).T
            self.rhs += self.S.shape[0]
        if with_kernel:
            K = lu_solve(lu,self.G.T,trans=1,check_finite=False).T
            self.rhs += self.n
        return res,E,P,K

    def objective(self,m,W):
        r,E,P,_ = self.state(m)
        rw = r @ W
        reg, gp, _ = self.reg.evaluate(m)
        f = .5*float(np.vdot(r,rw).real)/self.den + reg
        g = np.sum((P.T @ rw.conj())*E,axis=1).real/self.den + gp
        return f,g

    def full_hessian(self,m):
        r,E,P,K = self.state(m,with_kernel=True)
        reg, gp, _ = self.reg.evaluate(m)
        g = np.sum((P.T @ r.conj())*E,axis=1).real/self.den+gp
        jtj = (P.conj().T @ P) * (E.conj() @ E.T)
        hres = (P.T @ r.conj() @ E.T)*K
        D = self.reg.D
        dm = D @ m
        w = (2*self.reg.smooth + self.reg.tv*self.reg.epsilon**2/(dm*dm+self.reg.epsilon**2)**1.5)/self.reg.edges
        H = (jtj+hres+hres.T).real/self.den+(D.T @ D.multiply(w[:,None])).toarray()
        return r,E,P,g,(H+H.T)/2, .5*np.vdot(r,r).real/self.den+reg

    def score(self,m):
        r,_,_,_= self.state(m,False)
        rel = np.linalg.norm(r)/np.sqrt(self.den)
        return {'phi':float(.5*rel**2+self.reg.evaluate(m)[0]),'relative_residual':float(rel)}

    def fit(self, m, W, niter, free=None):
        a = m if free is None else m[free]
        def fun(x):
            if free is None:
                v = x
            else:
                v = m.copy()
                v[free]=x
            f,g = self.objective(v,W)
            return 1e5*f,1e5*(g if free is None else g[free])
        fit = minimize(fun,a,method='L-BFGS-B',jac=True,
                       bounds=[(self.lo,self.hi)]*len(a) if free is None else None,
                       options={'maxiter':niter,'gtol':2e-7,'ftol':2e-15,'maxcor':20,'maxls':30})
        v=m.copy()
        if free is None:
            v=fit.x
        else:
            v[free]=fit.x
        return v,dict(iterations=fit.nit,evaluations=fit.nfev,status=str(fit.message))


def source_directions(w,m):
    r,E,P,g,H,phi=w.full_hessian(m)
    s=w.s
    free = (m>w.lo+1e-8) & (m<w.hi-1e-8)
    active = np.flatnonzero(~free)
    signs = np.where(m[active]<w.lo+1e-8,1.,-1.)
    fe,fu = eigh(H[np.ix_(free,free)])
    # Lagrange multiplier on a stationary face, and its derivative under a new data metric.
    v = np.eye(w.n)[:,active]
    v[free,:] = -np.linalg.lstsq(H[np.ix_(free,free)],H[np.ix_(free,active)],rcond=1e-11)[0]
    v *= signs
    raw = np.einsum('is,ia,it->ast', (r.conj().T @ P).T,v,E,optimize=True)/w.den
    mat = hermitian(raw)
    M = w.white @ mat @ w.white
    ev,U=np.linalg.eigh(M)
    old = np.trace(mat,axis1=1,axis2=2).real
    mu=signs*g[active]
    slope=s*ev[:,0]-old
    ratio=np.where((mu>0)&(slope<0),mu/np.maximum(-slope,1e-30),np.inf)
    i=int(np.argmin(ratio))
    a=np.sqrt(s)*w.white @ U[i,:,0]
    rng=np.random.default_rng(SEED)
    random=rng.standard_normal(s)+1j*rng.standard_normal(s)
    random=np.sqrt(s)*w.white @ (random/np.linalg.norm(random))
    _,resvec=eigh(w.white @ (r.conj().T @ r) @ w.white)
    residual=np.sqrt(s)*w.white @ resvec[:,-1]
    # Ignore free-material response to weights as a mechanism control.
    raw_fixed=hermitian(np.einsum('is,it->ist',(r.conj().T @ P).T[active],E[active]))*signs[:,None,None]/w.den
    val0,vec0=np.linalg.eigh(w.white @raw_fixed@w.white)
    slope0=s*val0[:,0]-np.trace(raw_fixed,axis1=1,axis2=2).real
    rat0=np.where((mu>0)&(slope0<0),mu/np.maximum(-slope0,1e-30),np.inf)
    j=int(np.argmin(rat0))
    unprofiled=np.sqrt(s)*w.white@vec0[j,:,0]
    weights={'profile':a,'fixed_face_ignored':unprofiled,'largest_residual':residual,'random':random}
    def slopes(beam):
        return np.einsum('s,ast,t->a',beam.conj(),mat,beam).real-old
    stat={'free':int(sum(free)),'active':len(active),'phi':float(phi),'free_hessian_min':float(fe[0]),
        'free_hessian_max':float(fe[-1]),'free_gradient_inf':float(max(abs(g[free]))),
        'min_multiplier':float(mu.min()),'index':int(active[i]),'prior_multiplier':float(mu[i]),
        'profile_rate':float(slope[i]),'predicted_weight':float(ratio[i]),
        'comparison_index_unprofiled':int(active[j]),'predictions':{}}
    for name,beam in weights.items():
        rate=slopes(beam)
        tau=np.where((mu>0)&(rate<0),mu/np.maximum(-rate,1e-30),np.inf)
        ni=int(np.argmin(tau))
        stat['predictions'][name]={'first_index':int(active[ni]),'first_weight':float(tau[ni]),
             'chosen_index_rate':float(rate[i]),'power':float(np.vdot(beam,w.gram@beam).real),
             'feed_norm':float(np.linalg.norm(beam))}
    return stat, weights, (r,E,P,g,H,free,active,mu,mat)


def evaluate_labels(labels,steps):
    for label in labels:
        _,pub,m=load_endpoint(label)
        w=Wave(pub)
        start=perf_counter(); cpu=process_time()
        geom,weights,intern=source_directions(w,m)
        r,E,P,g,H,free,active,mu,_=intern
        I=np.eye(w.s)
        dist=1e-5
        h=np.sin(np.arange(w.n)+1.); h/=np.linalg.norm(h)
        fd=(w.objective(m+dist*h,I)[1]-w.objective(m-dist*h,I)[1])/(2*dist)
        geom['hessian_vector_fd_relative']=float(np.linalg.norm(fd-H@h)/np.linalg.norm(fd))
        geom['geometry_cost']={'wall_s':perf_counter()-start,'cpu_s':process_time()-cpu,'LU':w.lu,'rhs':w.rhs}
        rows={}
        mats={'initial':m}
        for name,a in (weights.items() if steps else []):
            W=np.outer(a,a.conj())
            start=perf_counter();cpu=process_time();f0,r0=w.lu,w.rhs
            v,m1=w.fit(m,.04*I+.96*W,steps)
            q,m2=w.fit(v,I,steps*2)
            rows[name]={'temporary':{**m1,**w.score(v)},'restored':{**m2,**w.score(q)},
                        'moved_rms':float(np.linalg.norm(m-q)/np.sqrt(w.n)),
                        'wall_s':perf_counter()-start,'cpu_s':process_time()-cpu,'LU':w.lu-f0,'rhs':w.rhs-r0}
            mats[name]=q; mats[name+'_temporary']=v
        if steps:
            start=perf_counter();cpu=process_time();f0,r0=w.lu,w.rhs
            q,info=w.fit(m,I,steps*3)
            rows['continue']={'restored':{**info,**w.score(q)},'wall_s':perf_counter()-start,
                     'cpu_s':process_time()-cpu,'LU':w.lu-f0,'rhs':w.rhs-r0}
            mats['continue']=q
        np.savez_compressed(OUT/f'basin_{label}.npz',**mats,**{'beam_'+k:v for k,v in weights.items()})
        record={'label':label,'steps':steps,'geometry':geom,'branches':rows}
        (OUT/f'basin_{label}.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
        print(label,geom['free'],geom['free_hessian_min'],geom['predicted_weight'],
              {k:round(v['restored']['phi'],6) for k,v in rows.items()},flush=True)


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--labels',nargs='+',default=['ring','shell','bridge'])
    cli.add_argument('--steps',type=int,default=55)
    args=cli.parse_args()
    evaluate_labels(args.labels,args.steps)
