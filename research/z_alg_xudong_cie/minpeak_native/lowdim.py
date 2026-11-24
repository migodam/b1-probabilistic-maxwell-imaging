"""Use a projected wavefield to propose low-cost starts for the original native inverse."""

import argparse
from collections import Counter
import json

import numpy as np
from scipy.linalg import lu_factor, lu_solve, solve
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits

from search import HERE, N, Model, full_ring_inputs, groups, elapsed, clocks


def projected_model(spectrum,sensor,incident,simple_y,side):
    index=groups(side)
    dim=side**2
    counts=np.bincount(index)
    basis=np.eye(dim)[index]/np.sqrt(counts[index,None])
    native=Model(spectrum,sensor,incident,simple_y)
    g=np.empty((dim,dim),complex)
    for k in range(0,dim,20):
        g[:,k:k+20]=basis.T @ native.green(basis[:,k:k+20])
    return index,g,sensor@basis,basis.T@incident,native.count


class ProjectedFit:
    def __init__(self,g,s,b,y):
        self.g,self.s,self.b,self.y=g,s,b,y
        self.scale=np.linalg.norm(y)**2
        self.count=Counter()

    def calculate(self,m,weight):
        A=-m[:,None]*self.g
        A.flat[::len(m)+1]+=1
        factor=lu_factor(A,check_finite=False)
        state=lu_solve(factor,m[:,None]*self.b,check_finite=False)
        adjoint=lu_solve(factor,self.s.conj().T,trans=2,check_finite=False)
        error=self.y-self.s@state
        self.count['projected_dense_lu']+=1
        self.count['projected_physical_state_rhs']+=self.b.shape[1]
        self.count['projected_receiver_adjoint_rhs']+=self.s.shape[0]
        physical=float(np.linalg.norm(error))
        if weight is None:
            g=np.sum(((self.b+self.g@state)*(adjoint @ error).conj()).real,axis=1)
            return .5*physical**2/self.scale, -g/self.scale,physical
        Dm=adjoint.conj().T@adjoint
        Dm.flat[::len(Dm)+1]+=weight
        correction=adjoint @ solve(Dm,error,assume_a='pos',check_finite=False)
        state += lu_solve(factor,correction,check_finite=False)
        self.count['projected_relaxed_state_rhs']+=self.b.shape[1]
        loss=np.linalg.norm(self.s@state-self.y)**2 + weight*np.linalg.norm(correction)**2
        gradient=-weight*np.sum(((self.b+self.g@state)*correction.conj()).real,axis=1)
        return .5*loss/self.scale,gradient/self.scale,physical


def one(projected,spectrum,sensor,incident,obs,args):
    start=clocks()
    index,g,s,b,_=projected
    side=args.side
    yy,xx=(np.indices((side,side))+.5)*2/side-1.
    base=np.exp(-.5*((xx/.38)**2+(yy/.38)**2)).ravel()
    fit=ProjectedFit(g,s,b,obs.y)
    native=Model(spectrum,sensor,incident,obs.y)
    best=np.inf
    audits=[]
    best_m=None
    for amp in args.amplitudes:
        m=amp*base
        for weight in [1.e-5,1.e-4,1.e-3,1.e-2,.1,None]:
            stage=clocks()
            def loss(m):
                f,d,_=fit.calculate(m,weight)
                return f,d
            result=minimize(loss,m,jac=True,method='L-BFGS-B',bounds=[(0,4)]*len(m),
                            options={'ftol':2e-10,'gtol':2e-9,'maxfun':args.calls,
                                     'maxiter':args.calls,'maxcor':12})
            m=result.x
            f,g0,projected_residual=fit.calculate(m,weight)
            if weight in (.01,.1,None):
                fr,_,_=native.full(m[index],gradient=False)
                if fr<best:
                    best,best_m=float(fr),m[index].copy()
            else:
                fr=None
            row=dict(amplitude=amp,penalty=weight,physical_projected_residual=projected_residual,
                     native_residual=float(fr) if fr is not None else None,peak=float(m.max()),
                     objective=f,iterations=int(result.nit),evaluations=int(result.nfev),
                     status=str(result.message),**elapsed(stage))
            print(obs.id,row,flush=True)
            audits.append(row)
    np.savez_compressed(args.out/(obs.id+'.npz'),m_seed=best_m)
    (args.out/(obs.id+'.json')).write_text(json.dumps({'observation':obs.id,'native_residual':best,
        'details':audits,'native_counts':dict(native.count),'projected_counts':dict(fit.count),
        **elapsed(start)},indent=2,allow_nan=False)+'\n')
    print('DONE',obs.id,best,flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=type(HERE),default=HERE/'wri')
    parser.add_argument('--object',action='append')
    parser.add_argument('--noise',type=float,action='append')
    parser.add_argument('--amplitudes',type=float,nargs='+',default=[1.5,2.4,3.3])
    parser.add_argument('--side',type=int,default=16)
    parser.add_argument('--calls',type=int,default=110)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    s,p,b,data=full_ring_inputs()
    start=clocks()
    projected=projected_model(s,p,b,data[0].y,args.side)
    print('SHARED PROJECTOR',dict(projected[-1]),elapsed(start),flush=True)
    for obs in data:
        if obs.noise_pct in (args.noise or [1.,2.]) and (not args.object or obs.object_id in args.object):
            one(projected,s,p,b,obs,args)


if __name__=='__main__':
    with threadpool_limits(1):
        main()
