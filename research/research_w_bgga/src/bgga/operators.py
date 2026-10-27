"""Ordinary exact fixed-space absorption; not a BGGA novelty claim.

Scattered-field convention: (L_W - V_X) u = b_W + b_X.
The RHS b_X is assembled with analytic incident fields, never V_X times a
nodal interpolant. Both source and receiver terms are retained.
"""
import time
import numpy as np
from scipy.linalg import lu_factor,lu_solve
from scipy.sparse.linalg import splu,gmres,LinearOperator

def relative_backward(L,u,b):
    residual=np.linalg.norm(L@u-b)
    denom=np.linalg.norm(L.data if hasattr(L,'tocsc') else L)*np.linalg.norm(u)+np.linalg.norm(b)
    return dict(absolute=float(residual),backward=float(residual/denom) if denom else 0.,
        rhs_relative=float(residual/np.linalg.norm(b)) if np.linalg.norm(b) else None)

class FixedSpaceCache:
    def __init__(self,Lw,bw,R,direct,U,book,phase='G0'):
        self.L=Lw.astype(np.complex128).tocsc();self.b=np.asarray(bw,complex);self.R=R;self.direct=np.asarray(direct)
        self.U=np.asarray(U,int);self.book=book;self.phase=phase
        self.book.reserve(phase,'ordinary background sparse LU',factorizations=1)
        self.factor=splu(self.L,permc_spec='COLAMD')
        self.book.reserve(phase,'ordinary background primal',rhs=self.b.shape[1])
        self.uw=self.factor.solve(self.b);self.wall_data=self.direct+(R@self.uw).T
        self.images=None;self.C=None;self.A=None
    def build_ports(self,max_bytes=1500000000,batch=24):
        n,m=self.L.shape[0],len(self.U)
        required=16*(n*m+3*m*m)
        if required>max_bytes:raise MemoryError(f'public ports need at least {required} bytes; refusing allocation')
        self.book.reserve(self.phase,'public ROI Green ports',rhs=m,scope='all target ROI, not current target mask')
        G=np.empty((n,m),complex)
        for j in range(0,m,batch):
            z=np.zeros((n,min(batch,m-j)),complex,order='F');z[self.U[j:j+batch],np.arange(z.shape[1])]=1
            G[:,j:j+batch]=self.factor.solve(z)
        self.images=G;self.C=G[self.U];self.A=np.asarray(self.R@G)
        return dict(full_dofs=n,public_ports=m,estimated_bytes=required,setup_rhs=m)
    def full(self,V,bx):
        L=(self.L-V).tocsc();b=self.b+bx
        self.book.reserve(self.phase,'ordinary target-present sparse LU',factorizations=1)
        f=splu(L,permc_spec='COLAMD')
        self.book.reserve(self.phase,'ordinary full primal',rhs=b.shape[1])
        u=f.solve(b);return u,self.direct+(self.R@u).T,relative_backward(L,u,b)
    def woodbury(self,V,bx):
        if self.images is None:raise RuntimeError('explicit public port setup required')
        outside=V.copy();inside=np.zeros(V.shape[0],bool);inside[self.U]=True
        if np.any(np.asarray(np.abs(outside[~inside]).sum(axis=1)).ravel()>0) or np.any(bx[~inside]!=0):
            raise ValueError('target update lies outside frozen public ROI')
        v=V[self.U][:,self.U].toarray()
        mu=bx+V@self.uw
        # delta = G_W(mu + V delta), exact finite material feedback.
        core=np.eye(len(self.U))-self.C@v
        self.book.dense(self.phase,'ordinary exact resolvent core',len(self.U))
        q=lu_solve(lu_factor(core),self.C@mu[self.U])
        injected=mu[self.U]+v@q;delta=self.images@injected
        u=self.uw+delta;data=self.wall_data+(self.A@injected).T
        return u,data,dict(**relative_backward(self.L-V,u,self.b+bx),
            core_backward=relative_backward(core,q,self.C@mu[self.U]),
            target_structural_ports=len(self.U))
    def schur(self,V,bx,max_bytes=1500000000,batch=24):
        """Independent ordinary exterior elimination; f_X is explicit."""
        U=self.U;C=np.setdiff1d(np.arange(self.L.shape[0]),U)
        if 16*(len(C)*len(U)+3*len(U)**2)>max_bytes:raise MemoryError('Schur allocation refused')
        lc=self.L[C][:,C].tocsc();luc=self.L[U][:,C];lcu=self.L[C][:,U].tocsc()
        self.book.reserve(self.phase,'ordinary exterior sparse LU',factorizations=1)
        fc=splu(lc,permc_spec='COLAMD')
        self.book.reserve(self.phase,'ordinary exterior port RHS',rhs=len(U)+self.b.shape[1])
        Z=np.empty((len(C),len(U)),complex)
        for j in range(0,len(U),batch):Z[:,j:j+batch]=fc.solve(lcu[:,j:j+batch].toarray())
        z=fc.solve(self.b[C]);K=self.L[U][:,U].toarray()-luc@Z
        f=self.b[U]-luc@z;rw=np.asarray(self.R[:,U]-self.R[:,C]@Z)
        d=self.direct+(self.R[:,C]@z).T
        vx=V[U][:,U].toarray();core=K-vx
        self.book.dense(self.phase,'ordinary exact Schur core',len(U))
        a=lu_solve(lu_factor(core),f+bx[U]);u=np.empty_like(self.b)
        u[U]=a;u[C]=z-Z@a
        data=d+(rw@a).T
        return u,data,relative_backward(self.L-V,u,self.b+bx)
    def schur_from_paid_ports(self,V,bx):
        """Ordinary Schur from the exact same paid port cache, no extra RHS.

        C is the restricted wall inverse. K=C^-1 is the full exterior-eliminated
        operator when both inverses exist. No SVD/pinv or stability assertion.
        """
        if self.images is None:raise RuntimeError('paid public ports required')
        m=len(self.U)
        self.book.dense(self.phase,'ordinary Schur restricted inverse LU',m)
        cf=lu_factor(self.C);K=lu_solve(cf,np.eye(m))
        self.book.dense(self.phase,'ordinary Schur target core LU',m)
        q=lu_solve(lu_factor(K-V[self.U][:,self.U].toarray()),K@self.uw[self.U]+bx[self.U])
        inj=bx[self.U]+V[self.U][:,self.U]@q
        u=self.uw+self.images@inj
        return u,self.wall_data+(self.A@inj).T,relative_backward(self.L-V,u,self.b+bx)
    def krylov(self,V,bx):
        L=(self.L-V).tocsc();b=self.b+bx;iterations=[];u=[]
        for col in range(b.shape[1]):
            count=[0]
            def apply(x):
                self.book.reserve(self.phase,'ordinary background-preconditioner triangular RHS',rhs=1)
                return self.factor.solve(x)
            def cb(x):count[0]+=1
            x,info=gmres(L,b[:,col],M=LinearOperator(L.shape,matvec=apply,dtype=complex),
                rtol=1e-11,atol=0,maxiter=100,callback=cb,callback_type='pr_norm')
            if info:raise RuntimeError(f'ordinary GMRES NONCONVERGED {info}')
            u.append(x);iterations.append(count[0])
        u=np.column_stack(u)
        return u,self.direct+(self.R@u).T,dict(**relative_backward(L,u,b),iterations=iterations)

def goal_galerkin(L,b,R,direct,Q,book,phase='G0'):
    """Same-space ordinary goal ROM; Q provenance/setup costs are caller-owned."""
    book.dense(phase,'ordinary goal Galerkin core',Q.shape[1])
    core=Q.conj().T@(L@Q);a=lu_solve(lu_factor(core),Q.conj().T@b);u=Q@a
    return u,direct+(R@u).T,relative_backward(L,u,b)
