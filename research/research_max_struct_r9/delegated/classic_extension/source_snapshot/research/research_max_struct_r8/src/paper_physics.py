"""Paper acquisition and equal-area disk MoM; FFT forward, dense small inverse.

The disk is a pulse integration approximation to a square of the same area,
not a non-overlapping physical disk. No .45h empirical filling fraction.
"""
from collections import Counter
import time
import numpy as np
from scipy.fft import fftn, next_fast_len
from scipy.special import hankel1, jv
from scipy.sparse.linalg import LinearOperator, gmres, bicgstab
from scipy.linalg import lu_solve, lu_factor
from research.research_max_struct_r6.src.fft_tsom_ops import DiskPulseToeplitzGreen
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_z_alg_r4_diversity.src.physics import Regularizer
from research.research_max_struct_r4.src.cie_onehop import OneHopAbsorption


class EqualAreaGreen(DiskPulseToeplitzGreen):
    def __init__(self, centers, shape, radius, k, cost=None):
        self.shape=tuple(shape);self.n=int(np.prod(shape));self.cost=Counter() if cost is None else cost
        c=np.asarray(centers).reshape(self.shape+(2,));v0=c[1,0]-c[0,0];v1=c[0,1]-c[0,0]
        self.pad=tuple(next_fast_len(2*s-1) for s in self.shape)
        a,b=np.meshgrid(np.arange(1-self.shape[0],self.shape[0]),np.arange(1-self.shape[1],self.shape[1]),indexing='ij')
        d=np.linalg.norm(a[...,None]*v0+b[...,None]*v1,axis=-1);zero=(a==0)&(b==0)
        kernel=.5j*np.pi*k*radius*jv(1,k*radius)*hankel1(0,k*np.where(zero,1.,d))
        self.diagonal=.5j*np.pi*k*radius*hankel1(1,k*radius)-1
        kernel[zero]=self.diagonal;embed=np.zeros(self.pad,complex);embed[a%self.pad[0],b%self.pad[1]]=kernel
        self.spectrum=fftn(embed);self.cost['Green_kernel_FFT_setup']+=1


def acquisition(side, frequency, dense=True):
    h=2/side;axis=(np.arange(side)+.5)*h-1
    yy,xx=np.meshgrid(axis,axis,indexing='ij');coords=np.c_[xx.ravel(),yy.ravel()]
    k=2*np.pi*frequency/3e8;radius=h/np.sqrt(np.pi)
    a=2*np.pi*np.arange(20)/20;directions=np.c_[np.cos(a),np.sin(a)]
    E=np.exp(1j*k*(coords@directions.T));b=2*np.pi*np.arange(40)/40
    receivers=3.75*np.c_[np.cos(b),np.sin(b)]
    def external(d):return .5j*np.pi*k*radius*jv(1,k*radius)*hankel1(0,k*d)
    S=external(np.linalg.norm(receivers[:,None]-coords[None,:],axis=2))
    public=dict(E0=E,E=E,S=S,centers=coords,radius=radius,k=k,grid_shape=(side,side))
    if dense:
        d=np.linalg.norm(coords[:,None]-coords[None,:],axis=2);np.fill_diagonal(d,1.)
        G=external(d);np.fill_diagonal(G,.5j*np.pi*k*radius*hankel1(1,k*radius)-1);public['G']=G
    return public


def austria(coords, epsilon_r):
    x,y=np.asarray(coords).T
    disks=((x-.3)**2+(y-.6)**2<=.2**2)|((x+.3)**2+(y-.6)**2<=.2**2)
    r2=x*x+(y+.2)**2;return ((disks|((r2>=.3**2)&(r2<=.6**2)))*(epsilon_r-1)).astype(float)


def fft_forward(public,m,rtol=1e-9,maxiter=100,restart=100):
    t,c=time.perf_counter(),time.process_time();counter=Counter();shape=tuple(public['grid_shape'])
    green=EqualAreaGreen(public['centers'],shape,public['radius'],public['k'],counter)
    E,S=public['E0'],public['S'];m=np.asarray(m);n=len(m)
    def mv(z):return z-m*green.apply(z[:,None])[:,0]
    A=LinearOperator((n,n),matvec=mv,dtype=complex)
    diagonal=1-m*green.diagonal;P=LinearOperator((n,n),matvec=lambda z:z/diagonal,dtype=complex)
    J=np.empty_like(E);solves=[]
    for s in range(E.shape[1]):
        trace=[];b=m*E[:,s]
        J[:,s],info=gmres(A,b,M=P,rtol=rtol,atol=0.,restart=restart,maxiter=maxiter,callback=lambda v:trace.append(float(v)),callback_type='pr_norm')
        residual=float(np.linalg.norm(mv(J[:,s])-b)/max(np.linalg.norm(b),1e-300))
        solves.append(dict(source=s,info=int(info),iterations=len(trace),state_relative_residual=residual))
        if info:raise RuntimeError('Forward GMRES did not converge: '+str(solves[-1]))
    return S@J,dict(wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c,LU=0,physical_complex_RHS=E.shape[1],counts=dict(counter),solves=solves)


def make_model(public,prior=True):
    model=MaxwellModel(public);model.reg=Regularizer(model.shape,smooth=1e-4 if prior else 0.,tv=1e-4 if prior else 0.)
    return model


class ShapeOnlyGreen:
    def __init__(self,n,public=None):self.shape=(n,n);self.public=public


class PaperOneHopAbsorption(OneHopAbsorption):
    def __init__(self,G,shape,cost):
        if not isinstance(G,ShapeOnlyGreen):
            super().__init__(G,shape,cost);return
        self.cost=cost;grid=np.arange(np.prod(shape)).reshape(shape)
        self.indices=np.array([grid[i:i+2,j:j+2].ravel() for i in range(0,shape[0],2) for j in range(0,shape[1],2)])
        p=G.public;c=p['centers'][self.indices];d=np.linalg.norm(c[:,:,None]-c[:,None,:],axis=-1)
        k,r=p['k'],p['radius'];np.put_along_axis(d,np.broadcast_to(np.arange(4)[None,:,None],(len(d),4,1)),1.,axis=2)
        self.goff=.5j*np.pi*k*r*jv(1,k*r)*hankel1(0,k*d);self.goff[:,np.arange(4),np.arange(4)]=0
        self.n=int(np.prod(shape));self.nblocks=len(self.indices);self.diag_g=np.full(self.n,.5j*np.pi*k*r*hankel1(1,k*r)-1)
        cost['block_Green_submatrices_generated_without_dense_G']+=self.nblocks
        cost['block_offdiagonal_nonzero_slots']+=12*self.nblocks;cost['self_Green_diagonal_entries_extracted']+=self.n


class MatrixFreeModel:
    """Same disk-pulse VIE, iterative full-physics state and analytic adjoint.

    Supports scalar CIE/NIE/LSIE arm only; no fake dense G for block methods.
    The observed residual is checked for every source/adjoint solve.
    """
    def __init__(self,public,prior=True):
        self.E=public['E0'];self.S=public['S'];self.shape=tuple(public['grid_shape']);self.n,self.ns=self.E.shape
        self.cost=Counter();self.G=ShapeOnlyGreen(self.n,public);self._current=None;self._solve_initial={}
        self.paper_geometry=dict(coords=public['centers'],radius=public['radius'],k=public['k'])
        self.green=EqualAreaGreen(public['centers'],self.shape,public['radius'],public['k'],self.cost)
        self.reg=Regularizer(self.shape,smooth=1e-4 if prior else 0.,tv=1e-4 if prior else 0.)
        self._dense_m=None;self._dense_lu=None
    def dense_fallback(self,m,B,adjoint=False):
        """Exactsame discreteoperator; retained libraryLU for all source RHS.

        Row-block assembly avoids an n*n*2 displacement tensor. Count failed
        iterative work separately and charge complete LU/setup/all RHS.
        """
        if self._dense_m is None or not np.array_equal(m,self._dense_m):
            t,c=time.perf_counter(),time.process_time();p=self.G.public;k,r=p['k'],p['radius'];coords=p['centers']
            A=np.empty((self.n,self.n),dtype=complex,order='F')
            for i in range(0,self.n,128):
                end=min(i+128,self.n);d=np.linalg.norm(coords[i:end,None]-coords[None,:],axis=2)
                d[np.arange(end-i),np.arange(i,end)]=1
                A[i:end]=(-m[i:end,None])*(.5j*np.pi*k*r*jv(1,k*r)*hankel1(0,k*d))
            np.fill_diagonal(A,1-m*self.green.diagonal)
            self._dense_lu=lu_factor(A,overwrite_a=True,check_finite=False);self._dense_m=m.copy()
            self.cost['LU']+=1;self.cost['dense_fallback_operator_assemblies']+=1;self.cost['dense_fallback_setup_wall_seconds']+=time.perf_counter()-t;self.cost['dense_fallback_setup_CPU_seconds']+=time.process_time()-c
        self.cost['dense_fallback_calls']+=1;self.cost['actual_complex_RHS']+=B.shape[1];self.cost['dense_fallback_complex_RHS']+=B.shape[1]
        Z=lu_solve(self._dense_lu,B,trans=2 if adjoint else 0,check_finite=False)
        residual=Z-self.green.adjoint(m[:,None]*Z)-B if adjoint else Z-m[:,None]*self.green.apply(Z)-B
        rel=np.linalg.norm(residual,axis=0)/np.maximum(np.linalg.norm(B,axis=0),1e-300)
        if not np.all(np.isfinite(Z)) or np.max(rel)>2e-9:raise RuntimeError('Matrix-free physical solve failed densefallback residual='+str(float(np.max(rel))))
        self.cost['maximum_dense_fallback_relative_residual']=max(self.cost.get('maximum_dense_fallback_relative_residual',0),float(np.max(rel)))
        self._solve_initial[adjoint]=Z.copy();return Z
    def solve(self,m,B,adjoint=False):
        if self._dense_m is not None and np.array_equal(m,self._dense_m):
            return self.dense_fallback(m,B,adjoint)
        def mv(v):
            return v-self.green.adjoint((m*v)[:,None])[:,0] if adjoint else v-m*self.green.apply(v[:,None])[:,0]
        A=LinearOperator((self.n,self.n),matvec=mv,dtype=complex)
        diagonal=1-m*(np.conj(self.green.diagonal) if adjoint else self.green.diagonal)
        P=LinearOperator((self.n,self.n),matvec=lambda v:v/diagonal,dtype=complex);Z=np.empty_like(B)
        for j in range(B.shape[1]):
            old=self._solve_initial.get(adjoint);x0=None if old is None or old.shape!=B.shape else old[:,j]
            it=[];Z[:,j],info=gmres(A,B[:,j],M=P,x0=x0,rtol=1e-9,atol=0,restart=100,maxiter=10,callback=lambda v:it.append(v),callback_type='pr_norm')
            self.cost['GMRES_iterations']+=len(it);self.cost['actual_complex_RHS']+=1
            if info:
                self.cost['failed_GMRES_attempts']+=1
                return self.dense_fallback(m,B,adjoint)
            res=float(np.linalg.norm(mv(Z[:,j])-B[:,j])/max(np.linalg.norm(B[:,j]),1e-300))
            self.cost['maximum_GMRES_relative_residual']=max(self.cost.get('maximum_GMRES_relative_residual',0),res)
            if info or res>2e-9:return self.dense_fallback(m,B,adjoint)
        self._solve_initial[adjoint]=Z.copy();return Z
    def state(self,m):
        if self._current is not None and np.array_equal(self._current['m'],m):return self._current
        self.cost['FULL_FORWARD']+=1;j=self.solve(m,m[:,None]*self.E);Q=self.E+self.green.apply(j)
        self._current=dict(m=m.copy(),j=j,Q=Q,Y=self.S@j,state_relative_residual=float(np.linalg.norm(j-m[:,None]*Q)/max(np.linalg.norm(j),1e-300)))
        return self._current
    def physical_score(self,m,Y):
        r=self.state(m)['Y']-Y;return .5*float(np.vdot(r,r).real/np.vdot(Y,Y).real)+self.reg.evaluate(m)[0]
    def data_value_gradient(self,m,Y):
        st=self.state(m);r=st['Y']-Y;scale=np.vdot(Y,Y).real;a=self.solve(m,self.S.conj().T@r,adjoint=True)
        prior,gprior,_=self.reg.evaluate(m);self.cost['DATA_ADJOINT']+=1
        return .5*float(np.vdot(r,r).real/scale)+prior,np.real(np.sum(a.conj()*st['Q'],axis=1))/scale+gprior
