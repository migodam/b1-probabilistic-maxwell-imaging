"""Free-space volume-weighted dyadic DDA. J=chi E, cell-major Cartesian layout."""
import os
os.environ.setdefault("OMP_NUM_THREADS","1")
from collections import Counter
import numpy as np
import torch
from pathlib import Path
import time,json
ROOT=Path(__file__).resolve().parent
K=2*np.pi
WIDTH=.6
DEV=torch.device(os.environ.get("SRAP_DEVICE","cuda"))
CD=torch.complex64
torch.set_num_threads(1)
def t(a):return torch.as_tensor(np.asarray(a),device=DEV,dtype=CD)
def npy(a):return a.detach().cpu().numpy()
def write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,allow_nan=False))
def green(dv):
    d=np.linalg.norm(dv,axis=-1);u=dv/d[...,None];nn=u[...,None,:]*u[...,:,None];eye=np.eye(3)
    return (np.exp(1j*K*d)/(4*np.pi*d))[...,None,None]*(K**2*(eye-nn)+(1j*K/d-1/d**2)[...,None,None]*(eye-3*nn))
def ports(n,geom,phase):
    if geom=='plane_grid':
        # Shifted 2x4 transmitter and 4x4 receiver rectangular lattices.
        if n==8:
            xs=np.array([-3,-1,1,3])*.09
            ys=np.array([-1,1])*np.sqrt(2*.08**2-.09**2)
        elif n==16:
            xs=ys=np.array([-3,-1,1,3])*.08
        else:raise ValueError(n)
        x,y=np.meshgrid(xs,ys,indexing='xy')
        z=np.sqrt(.9**2-2*.08**2)
        pos=np.stack([x.ravel(),y.ravel(),np.full(n,z)],axis=1)
        pol=np.tile(np.eye(3)[:2],(n,1,1))
        return pos,pol
    if geom!='hemi':raise ValueError(geom)
    z=(np.arange(n)+.5)/n
    az=2*np.pi*((np.arange(n)*.61803398875+phase/n)%1)
    c=np.sqrt(1-z*z)
    dirs=np.stack([c*np.cos(az),c*np.sin(az),z],axis=1)
    p1=np.stack([-np.sin(az),np.cos(az),0*az],axis=1);p2=np.cross(dirs,p1)
    return .9*dirs,np.stack([p1,p2],axis=1)
class Model:
    def __init__(self,n,geom):
        self.n=n;self.h=WIDTH/n;self.N=n**3;self.M=3*self.N;self.count=Counter()
        idx=np.indices((n,n,n)).reshape(3,-1).T;self.pts=(idx-(n-1)/2)*self.h
        p=2*n;v=np.arange(p);v=np.where(v<=n,v,v-p)*self.h
        dv=np.stack(np.meshgrid(v,v,v,indexing='ij'),axis=-1);dv[0,0,0]=[1,0,0]
        ker=self.h**3*green(dv);self.selfterm=-1/3+1j*K**3*self.h**3/(6*np.pi);ker[0,0,0]=self.selfterm*np.eye(3)
        self.kernel=torch.fft.fftn(t(ker.transpose(3,4,0,1,2).copy()),dim=(-3,-2,-1))
        src,sp=ports(8,geom,.18);rx,rp=ports(16,geom,.63)
        self.E=t(np.einsum('nsij,soj->niso',green(self.pts[:,None,:]-src[None,:,:]),sp).reshape(self.M,16))
        self.S=t((self.h**3*np.einsum('roi,rnij->ronj',rp,green(rx[:,None,:]-self.pts[None,:,:]))).reshape(32,self.M).copy())
    def G(self,x,adj=False):
        self.count['Gadj_rhs' if adj else 'G_rhs']+=x.shape[1]
        n=self.n;rhs=x.shape[1];a=x.T.reshape(rhs,n,n,n,3).permute(0,4,1,2,3).contiguous()
        z=torch.fft.fftn(a,s=(2*n,)*3,dim=(-3,-2,-1));ker=self.kernel.conj().transpose(0,1) if adj else self.kernel
        # Explicit three-component contraction keeps FFT batch layout contiguous.
        out=torch.stack([sum(ker[i,j]*z[:,j] for j in range(3)) for i in range(3)],dim=1)
        y=torch.fft.ifftn(out,dim=(-3,-2,-1))[:,:,:n,:n,:n]
        return y.permute(2,3,4,1,0).reshape(self.M,rhs)
    def solve(self,m,rhs,adj=False,tol=2e-5):
        """Per-source preconditioned BiCGSTAB; original-equation certification."""
        d=m.repeat_interleave(3)[:,None];diag=1-d*self.selfterm
        if adj:diag=diag.conj()
        def A(v):return v-self.G(d*v,True) if adj else v-d*self.G(v)
        self.count['adjoint_solve_rhs' if adj else 'forward_solve_rhs']+=rhs.shape[1]
        x=torch.zeros_like(rhs);r=rhs.clone();rh=r.clone();v=torch.zeros_like(r);p=v.clone()
        sz=torch.linalg.vector_norm(rhs,dim=0).clamp_min(1e-18);one=torch.ones_like(sz,dtype=CD)
        rho=one;alpha=one;omega=one
        def divide(a,b):return torch.where(b.abs()>1e-30,a/b,torch.zeros_like(a))
        for it in range(160):
            active=(torch.linalg.vector_norm(r,dim=0)>tol*sz)
            if not active.any().item():break
            rn=(rh.conj()*r).sum(0);beta=divide(rn,rho)*divide(alpha,omega)
            p=r+beta*(p-omega*v);ph=p/diag;v=A(ph);alpha=divide(rn,(rh.conj()*v).sum(0))
            s=r-alpha*v;sh=s/diag;u=A(sh);omega=divide((u.conj()*s).sum(0),(u.conj()*u).sum(0))
            x=x+torch.where(active,alpha*ph+omega*sh,0);r=torch.where(active,s-omega*u,0);rho=rn
        err=(torch.linalg.vector_norm(A(x)-rhs,dim=0)/sz).max().item()
        self.count['bicg_iterations']+=it+1
        if not np.isfinite(err) or err>5*tol:raise ArithmeticError(('BiCGSTAB',adj,err,it))
        return x,err
    def physical(self,m,Y=None,grad=False):
        d=m.repeat_interleave(3)[:,None];J,err=self.solve(m,d*self.E);pred=self.S@J
        if Y is None:return pred,err
        diff=pred-Y;weights=1/(abs(Y)**2).sum(0);phi=((abs(diff)**2).sum(0)*weights).mean()
        rel=(torch.linalg.norm(diff)/torch.linalg.norm(Y)).item()
        if not grad:return float(phi),rel
        lam,_=self.solve(m,self.S.mH@(diff*weights/Y.shape[1]),True)
        g=2*(lam.conj()*(self.E+self.G(J))).real.reshape(self.N,3,-1).sum((1,2))
        return float(phi),npy(g).astype(float)
