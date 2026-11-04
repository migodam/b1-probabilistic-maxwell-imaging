"""Exact finite-background material actions; no series truncation or learning.

The small port objective contracts source/receiver matrices once. A template
with k changed scalar cells has state rank k even with only two amplitudes.
"""
import time
import numpy as np
from scipy.linalg import lu_solve
from scipy.optimize import minimize

class RationalAnchor:
    def __init__(self,model,Y,m,pool):
        self.model,self.Y,self.m=model,Y,m.copy()
        self.state=model.state(m);self.scale=float(np.vdot(Y,Y).real)
        self.error=self.state['Y']-Y;self.E=self.state['Q'];self.pool=np.asarray(pool,int)
        self.base_data=.5*float(np.vdot(self.error,self.error).real)/self.scale
        self.base_prior=model.reg.evaluate(m)[0]
        start,cpu=time.perf_counter(),time.process_time()
        with model._physics.operation('RATIONAL_RECEIVER_TRANSPOSE',rhs=model.nr,receiver_transpose_RHS=model.nr):
            self.O=lu_solve(self.state['lu'],model.S.T,trans=1,check_finite=False).T
        P=np.eye(model.n,dtype=complex)[:,self.pool]
        with model._physics.operation('RATIONAL_SELECTED_COLUMNS',rhs=len(self.pool),selected_resolvent_RHS=len(self.pool)):
            self.C=model.G@lu_solve(self.state['lu'],P,check_finite=False)
        self.setup_wall=time.perf_counter()-start;self.setup_cpu=time.process_time()-cpu
        self.lookup={int(i):j for j,i in enumerate(self.pool)}

    def action(self,groups,weights=None):return FiniteAction(self,groups,weights)

class FiniteAction:
    def __init__(self,anchor,groups,weights=None):
        self.anchor=anchor;self.groups=[np.array(g,int) for g in groups]
        self.indices=np.concatenate(self.groups);self.rank=len(self.indices)
        self.assignment=np.zeros((self.rank,len(groups)));start=0
        for j,g in enumerate(groups):
            self.assignment[start:start+len(g),j]=1. if weights is None else weights[j]
            start+=len(g)
        col=[anchor.lookup[int(i)] for i in self.indices]
        self.K=anchor.C[np.ix_(self.indices,col)]
        self.O=anchor.O[:,self.indices];self.E=anchor.E[self.indices]
        self.A=self.O.conj().T@self.O/anchor.scale
        self.B=self.E@self.E.conj().T
        self.Z=self.O.conj().T@anchor.error@self.E.conj().T/anchor.scale
        self.bounds=[]
        for j,g in enumerate(self.groups):
            w=self.assignment[self.assignment[:,j]!=0,j]
            lo=(.02-anchor.m[g])/w;hi=(20.-anchor.m[g])/w
            self.bounds.append((float(np.maximum(np.minimum(lo,hi),-1e20).max()),float(np.minimum(np.maximum(lo,hi),1e20).min())))
        D=anchor.model.reg.D[:,self.indices].tocsr();used=np.flatnonzero(np.asarray(abs(D).sum(axis=1)).ravel())
        self.edge_base=(anchor.model.reg.D@anchor.m)[used]
        self.edge_action=(D[used]@self.assignment)
        self.queries=0
        if self.rank==2 and weights is not None:raise ValueError('Weighted rank2 Gram not enabled; use full small resolvent')
        if self.rank==2:
            k=self.K;O=self.O;E=self.E;e=anchor.error/np.sqrt(anchor.scale)
            U=O[:,0,None]*E[0,None,:]/np.sqrt(anchor.scale)
            V=O[:,1,None]*E[1,None,:]/np.sqrt(anchor.scale)
            H=-k[1,1]*U-k[0,0]*V+(k[0,1]*O[:,0,None]*E[1,None,:]+k[1,0]*O[:,1,None]*E[0,None,:])/np.sqrt(anchor.scale)
            c=k[0,0]*k[1,1]-k[0,1]*k[1,0]
            T=np.stack((e,U-k[0,0]*e,V-k[1,1]*e,H+c*e)).reshape(4,-1)
            self.gram=(T.conj()@T.T).real;self.det_coeff=np.array([1,-k[0,0],-k[1,1],c])

    def prior(self,x):
        reg=self.anchor.model.reg
        d=self.edge_base+self.edge_action@x
        def f(a):return reg.smooth*np.sum(a*a)+reg.tv*np.sum(np.sqrt(a*a+reg.epsilon**2)-reg.epsilon)
        return self.anchor.base_prior+(f(d)-f(self.edge_base))/reg.edges

    def evaluate(self,x):
        self.queries+=1;x=np.asarray(x,float)
        if self.rank==2:
            p=np.array([1,x[0],x[1],x[0]*x[1]])
            determinant=p@self.det_coeff
            data=.5*(p@self.gram@p)/abs(determinant)**2
        else:
            delta=self.assignment@x
            T=np.linalg.solve(np.eye(self.rank)-delta[:,None]*self.K,np.diag(delta))
            data=self.anchor.base_data+np.vdot(self.Z,T).real+.5*np.trace(T.conj().T@self.A@T@self.B).real
        return float(data+self.prior(x))

    def material(self,x):
        m=self.anchor.m.copy();m[self.indices]+=self.assignment@x
        return m

    def port_delta(self,x):
        delta=self.assignment@x
        return self.O@np.linalg.solve(np.eye(self.rank)-delta[:,None]*self.K,delta[:,None]*self.E)

    def search(self,points=9):
        # Candidate domain search, fixed before endpoint/GT inspection.
        axes=[np.unique(np.r_[0,np.linspace(a,b,points)]) for a,b in self.bounds]
        xx=np.array(np.meshgrid(*axes,indexing='ij')).reshape(2,-1).T
        values=np.array([self.evaluate(x) for x in xx])
        local=minimize(self.evaluate,np.zeros(2),method='L-BFGS-B',bounds=self.bounds,options=dict(maxiter=100,ftol=1e-13,gtol=1e-9))
        best=(float(local.fun),local.x.copy(),'local_exact_from_origin')
        for idx in np.argsort(values)[:3]:
            z=minimize(self.evaluate,xx[idx],method='L-BFGS-B',bounds=self.bounds,options=dict(maxiter=100,ftol=1e-13,gtol=1e-9))
            if z.fun<best[0]:best=(float(z.fun),z.x.copy(),'finite_grid_then_exact_local')
        return dict(best_Phi=best[0],best_x=best[1],selected_kind=best[2],origin_Phi=self.evaluate([0,0]),local_origin_Phi=float(local.fun),local_origin_x=local.x,grid_points=xx,grid_Phi=values,queries=self.queries,state_rank=self.rank,amplitude_dimension=2,groups=[g.tolist() for g in self.groups])
