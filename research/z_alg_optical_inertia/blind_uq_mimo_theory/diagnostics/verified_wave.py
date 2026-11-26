"""A posteriori wave bounds for exact stored Fourier and multiport coefficients.

The validation transform uses explicit unitary DFT matrices with Arb-enclosed
twiddles. FFTs and GMRES are used freely to propose numerical vectors.
"""
from pathlib import Path
from time import perf_counter, process_time
import os, time
for name in ('OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OMP_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.fft import fft2,ifft2
from scipy.linalg import lu_factor,lu_solve
from flint import arb,acb,ctx
ctx.prec=128
P=Path(__file__).resolve().parents[1]
EPS=2.**-53
def up(x):
    return float(np.nextafter(float(x)*(1+128*EPS)+1e-290,np.inf))
def fro(x):
    z=np.asarray(x); return up(np.linalg.norm(z.ravel())*(1+EPS*16*z.size))
def opnorm(x):
    return up(min(fro(x),np.sqrt(np.linalg.norm(x,1)*np.linalg.norm(x,np.inf))*(1+32*EPS*sum(x.shape))))
def rd(n):
    return up(32*(n+3)*EPS/(1-32*(n+3)*EPS))
def arb_upper(x):return float(np.nextafter(float(x.upper()),np.inf))
class VerifiedWave:
    def __init__(self,op):
        self.sym=op['spectrum']; self.side=int(op['grid_shape'][0]);self.n=self.side**2
        self.p=self.sym.shape[0];assert self.sym.shape==(self.p,self.p)
        self.S=op['S']; self.B=op['E0']; self.L=op['low_rank_left'];self.R=op['low_rank_right']
        self.snorm=fro(self.S);self.g0=up(np.max(abs(self.sym))*(1+8*EPS))
        self.lnorm=opnorm(self.L);self.rnorm=opnorm(self.R);self.gnorm=up(self.g0+self.lnorm*self.rnorm)
        d=self.p; w=acb(0,2*arb.pi()/d).exp();sq=arb(d).sqrt()
        integers=np.arange(d); idx=(integers[:,None]*integers)%d
        exact=[w**i/sq for i in range(d)]
        num=np.array([complex(a) for a in exact])
        self.tw_error=max(arb_upper(abs(a-acb(float(z.real),float(z.imag)))) for a,z in zip(exact,num))
        self.U=num[idx] # unitary inverse DFT
        self.one=up(d*self.tw_error+rd(d)*(np.sqrt(d)+d*self.tw_error))
        self.two=up(2*self.one+self.one**2)
        self.full_green_err=up((2*self.two+self.two**2+rd(2)*(1+self.two)**2)*self.g0+
                 (rd(self.n)+rd(self.L.shape[1])+rd(self.n)*rd(self.L.shape[1]))*fro(self.L)*fro(self.R)
                 +rd(2)*self.gnorm)
        # Exact lag samples of inverse DFT, evaluated at arbitrary precision.
        v=[[acb(float(z.real),float(z.imag)) for z in row] for row in self.sym]
        cols=zip(*(acb.dft(row,inverse=True) for row in v))
        rows=list(zip(*(acb.dft(list(col),inverse=True) for col in cols)))
        self.lag=np.array([[complex(z) for z in row] for row in rows])
        self.lagerr=max(arb_upper(abs(z-acb(float(c.real),float(c.imag)))) for row,rr in zip(rows,self.lag) for z,c in zip(row,rr))
    def transform(self,arr,reverse=False):
        # Channels first for independent BLAS matrix products.
        D=self.U if reverse else self.U.conj()
        return (D@np.asarray(arr,order='C'))@D.T
    def apply(self,q,strict=False):
        flat=q.ndim==1
        if flat:q=q[:,None]
        if strict:
            a=np.zeros((q.shape[1],self.p,self.p),complex)
            a[:,:self.side,:self.side]=q.T.reshape((-1,self.side,self.side))
            t=self.transform(a)*self.sym
            b=self.transform(t,True)[:,:self.side,:self.side].reshape(q.shape[1],-1).T
        else:
            t=fft2(q.reshape(self.side,self.side,-1),s=self.sym.shape,axes=(0,1))*self.sym[...,None]
            b=ifft2(t,axes=(0,1))[:self.side,:self.side].reshape(q.shape)
        ret=b+self.L@(self.R@q)
        return ret[:,0] if flat else ret
    def matrix(self,idx):
        y,x=divmod(idx,self.side)
        g=self.lag[(y[:,None]-y[None,:])%self.p,(x[:,None]-x[None,:])%self.p]
        return g+self.L[idx]@self.R[:,idx]
    def active_inv(self,x,idx):
        h=len(idx);u=float(np.max(abs(x)));g=self.matrix(idx)
        gerr=up(h*self.lagerr+rd(self.R.shape[0])*fro(self.L[idx])*fro(self.R[:,idx])+rd(1)*fro(g))
        A=np.eye(h)-x[idx,None]*g; Aerr=up(u*gerr+rd(2)*(np.sqrt(h)+u*fro(g)))
        F=lu_factor(A,check_finite=False);Q=lu_solve(F,np.eye(h),check_finite=False)
        nq=opnorm(Q)
        rho=up(fro(np.eye(h)-A@Q)+rd(h)*(np.sqrt(h)+fro(A)*fro(Q))+Aerr*nq)
        if rho>=1:raise ArithmeticError('Direct inverse defect >= 1')
        active=up(nq/(1-rho))
        # Block triangular I-XG has identity rows on exact zero cells.
        c=active if h==self.n else up(1+active*(1+u*self.gnorm))
        return c,dict(method='active full inverse residual',active=h,defect=rho,operator_error=gerr)
    def fourier_inv(self,x,q0=.78):
        # A=H-U V*, ||I-H|| <= q. Find Z with A Z approximately U.
        # Then A(I+Z V*)=H+(AZ-U)V*, and obtain a full inverse bound.
        u=float(x.max()); select=np.flatnonzero(abs(self.sym.ravel())>q0/u)
        kt=len(select);rank=kt+self.R.shape[0];modes=np.unravel_index(select,self.sym.shape)
        t=self.sym.copy();t.flat[select]=0
        q=up(u*np.max(abs(t))*(1+8*EPS))
        if q>=1:raise ArithmeticError('Fourier tail is too large')
        y,gx=np.indices((self.side,self.side));y=y.ravel();gx=gx.ravel()
        yy,xx=modes
        # Columns of compression of the padded exact two-dimensional inverse DFT.
        f=self.U[y[:,None],yy[None,:]]*self.U[gx[:,None],xx[None,:]]
        left=x[:,None]*np.column_stack((f*self.sym.ravel()[select],self.L))
        nuerr=up(u*(np.sqrt(self.n*kt)*self.g0*(2*self.tw_error/np.sqrt(self.p)+self.tw_error**2+rd(3)/self.p)+rd(3)*fro(self.L)))
        vn=up(np.sqrt(1+self.rnorm**2))
        def vstar(mat,strict):
            if strict:
                a=np.zeros((mat.shape[1],self.p,self.p),complex)
                a[:,:self.side,:self.side]=mat.T.reshape((-1,self.side,self.side))
                values=self.transform(a)[:,yy,xx].T
            else:
                values=fft2(mat.reshape(self.side,self.side,-1),s=self.sym.shape,axes=(0,1))[yy,xx]/self.p
            return np.vstack((values,self.R@mat))
        def tail(mat):
            values=fft2(mat.reshape(self.side,self.side,-1),s=self.sym.shape,axes=(0,1))*t[...,None]
            return x[:,None]*ifft2(values,axes=(0,1))[:self.side,:self.side].reshape(mat.shape)
        D=np.empty_like(left); reduced=np.empty((rank,rank),complex)
        chunk=12
        for j in range(0,rank,chunk):
            sl=slice(j,j+chunk); z=left[:,sl].copy();v=z.copy()
            for _ in range(32):
                v=tail(v);z+=v
                if fro(v)<1e-9:break
            D[:,sl]=z;reduced[:,sl]=vstar(z,False)
        inverse=lu_solve(lu_factor(np.eye(rank)-reduced,check_finite=False),np.eye(rank),check_finite=False)
        Z=D@inverse
        sq=0.
        for j in range(0,rank,chunk):
            sl=slice(j,j+chunk);z=Z[:,sl];g=self.apply(z,True);r=z-x[:,None]*g-left[:,sl]
            err=up(u*self.full_green_err*fro(z)+rd(4)*(fro(z)+u*fro(g)+fro(left[:,sl])))
            sq=up(sq+up(fro(r)+err)**2)
        res=up(np.sqrt(sq)+nuerr); gap=1-q-res*vn
        if gap<=0:raise ArithmeticError(f'Fourier residual too large: q={q} res={res} vn={vn}')
        c=up((1+fro(Z)*vn)/(gap*(1-128*EPS)))
        return c,dict(method='native Fourier tail with independently evaluated DFT residual',rank=rank,q=q,
                      approximate_block_defect=res,verified_gap_lower=gap*(1-128*EPS),block_size=chunk)
    def inverse(self,x):
        i=np.flatnonzero(x);t=time.perf_counter();cpu=time.process_time()
        if len(i)==0:c,log=1.,dict(method='air')
        elif self.side==32 or len(i)<=1600:c,log=self.active_inv(x,i)
        else:c,log=self.fourier_inv(x)
        return c,dict(**log,inverse_norm_upper=c,wall_s=time.perf_counter()-t,cpu_s=time.process_time()-cpu)
    def solve_error(self,x,y,b,binerr,inv):
        g=self.apply(y,True); res=y-x[:,None]*g-b
        ue=up(fro(res)+max(x)*self.full_green_err*fro(y)+rd(4)*(fro(y)+max(x)*fro(g)+fro(b))+binerr)
        return up(inv*ue)
