"""CPU 3D vector CM+RR DDA adapter, exp(-i omega t), lambda0=Z0=1.

Fields are columns (3N,P), cell-major/component-minor. Physical epsilon_r
is a symmetric tensor; mu_r=1. This is a point-dipole discretization, not a
compatible discrete de Rham complex or a continuum error certificate.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
from pathlib import Path
import time
import warnings
import numpy as np
from scipy import linalg as la
try:
    from .a9_green import build_goff, dyadic_green
except ImportError:
    from a9_green import build_goff, dyadic_green


def block_diagonal(tensors):
    t = np.asarray(tensors, complex)
    out = np.zeros((3*len(t),3*len(t)),complex)
    for i,b in enumerate(t): out[3*i:3*i+3,3*i:3*i+3] = b
    return out


def plane_waves(points,k,dirs,pols):
    dirs=np.asarray(dirs,float).reshape(-1,3)
    pols=np.asarray(pols,complex).reshape(-1,3)
    if len(dirs)!=len(pols) or np.any(abs(la.norm(dirs,axis=1)-1)>1e-12):
        raise ValueError('directions must be unit and match polarizations')
    if np.any(abs(np.sum(dirs*pols,axis=1))>1e-12):
        raise ValueError('polarization must be transverse')
    phase=np.exp(1j*k*(np.asarray(points)@dirs.T))
    return (phase[:,None,:]*pols.T[None,:,:]).reshape(3*len(points),len(dirs))


def magnetic_green(dv,k):
    """H from an electric dipole p: curl(G_E p)/(ik), Z0=1.

    The odd kernel has zero symmetric-cell principal self contribution.
    """
    d=la.norm(dv,axis=-1)
    if np.any(d<=0): raise ValueError('external magnetic Green requires separation')
    n=dv/d[...,None]
    cross=np.zeros(dv.shape[:-1]+(3,3),complex)
    cross[...,0,1]=-n[...,2];cross[...,0,2]=n[...,1]
    cross[...,1,0]=n[...,2];cross[...,1,2]=-n[...,0]
    cross[...,2,0]=-n[...,1];cross[...,2,1]=n[...,0]
    return (k*k*np.exp(1j*k*d)/(4*np.pi*d)*(1+1j/(k*d)))[...,None,None]*cross


@dataclass
class VectorState:
    model: object
    epsilon: np.ndarray
    directions: np.ndarray
    polarizations: np.ndarray
    E: np.ndarray
    H: np.ndarray
    j: np.ndarray
    J: np.ndarray
    p: np.ndarray
    E_exc: np.ndarray
    incident: np.ndarray
    matrix: np.ndarray
    factor: tuple
    receipt: dict

    def near(self,receivers):
        """Full external near fields; receivers must lie outside all cells."""
        r=np.asarray(receivers,float).reshape(-1,3)
        # The adapter contract excludes cell interiors, not only exact centers.
        if np.any(np.max(abs(r[:,None,:]-self.model.points[None,:,:]),axis=-1)
                  <= self.model.volume**(1/3)/2):
            raise ValueError('near receivers must be outside occupied voxel cells')
        E_s=self.model.external_electric_operator(r)@self.j
        H_s=self.model.external_magnetic_operator(r)@self.j
        E_i=plane_waves(r,self.model.k,self.directions,self.polarizations)
        hpol=np.cross(self.directions,self.polarizations)
        H_i=plane_waves(r,self.model.k,self.directions,hpol)
        return dict(E=E_i+E_s,H=H_i+H_s,E_scattered=E_s,H_scattered=H_s)

    def far(self,directions):
        """Amplitudes F: E_sc(r n)=exp(ikr)/r F_E(n)+O(r^-2)."""
        d=np.asarray(directions,float).reshape(-1,3)
        if np.any(abs(la.norm(d,axis=1)-1)>1e-12): raise ValueError('unit far directions required')
        phase=np.exp(-1j*self.model.k*(d@self.model.points.T))
        dip=np.einsum('rn,ncp->rcp',phase,self.p.reshape(-1,3,self.E.shape[1]))
        E=self.model.k**2/(4*np.pi)*(dip-d[:,:,None]*np.einsum('rc,rcp->rp',d,dip)[:,None,:])
        H=np.cross(d[:,:,None],E,axisa=1,axisb=1,axisc=1)
        return dict(E=E.reshape(3*len(d),-1),H=H.reshape(3*len(d),-1))


class TensorDDA:
    def __init__(self,points,volume,k=2*np.pi):
        start=time.perf_counter(); cpu=time.process_time()
        self.points=np.asarray(points,float).reshape(-1,3)
        self.volume=float(volume); self.k=float(k); self.n=3*len(self.points)
        if self.n==0 or self.volume<=0 or self.k<=0: raise ValueError('invalid grid')
        self.Goff=build_goff(self.points,self.k)
        self.self_density=-1/3+1j*self.k**3*self.volume/(6*np.pi)
        self.G=self.volume*self.Goff+self.self_density*np.eye(self.n)
        self.setup_receipt=dict(wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu,
                                N=len(self.points),LU=0,RHS=0)

    def external_electric_operator(self,receivers):
        r=np.asarray(receivers,float).reshape(-1,3)
        g=dyadic_green(r[:,None,:]-self.points[None,:,:],self.k)
        return self.volume*g.transpose(0,2,1,3).reshape(3*len(r),self.n)

    def external_magnetic_operator(self,receivers):
        r=np.asarray(receivers,float).reshape(-1,3)
        g=magnetic_green(r[:,None,:]-self.points[None,:,:],self.k)
        return self.volume*g.transpose(0,2,1,3).reshape(3*len(r),self.n)

    def _internal_magnetic_operator(self):
        dv=self.points[:,None,:]-self.points[None,:,:]
        mask=np.eye(len(self.points),dtype=bool)
        dv[mask]=[1.,0,0]
        g=magnetic_green(dv,self.k); g[mask]=0
        return self.volume*g.transpose(0,2,1,3).reshape(self.n,self.n)

    def solve(self,epsilon,directions,polarizations,receipt_path=None):
        start=time.perf_counter();cpu=time.process_time()
        receipt=dict(status='RUNNING',LU=0,RHS=0,N=len(self.points),P=len(directions),
                     precision='complex128',device='CPU',setup=self.setup_receipt)
        try:
            eps=np.asarray(epsilon,complex)
            if eps.shape!=(len(self.points),3,3): raise ValueError('epsilon must be (N,3,3)')
            if not np.all(np.isfinite(eps)): raise ValueError('nonfinite epsilon')
            if not np.allclose(eps,eps.transpose(0,2,1),rtol=0,atol=1e-12):
                raise ValueError('reciprocal epsilon must be complex symmetric')
            chi=block_diagonal(eps-np.eye(3))
            inc=plane_waves(self.points,self.k,directions,polarizations)
            L=np.eye(self.n)-self.G@chi
            receipt['LU']=1
            with warnings.catch_warnings():
                warnings.simplefilter('error',la.LinAlgWarning)
                fac=la.lu_factor(L,check_finite=True)
            receipt['RHS']=inc.shape[1]
            E=la.lu_solve(fac,inc,check_finite=True)
            j=chi@E;p=self.volume*j
            H=plane_waves(self.points,self.k,directions,np.cross(directions,polarizations))
            H+=self._internal_magnetic_operator()@j
            E_exc=E-(self.self_density/self.volume)*p
            receipt.update(status='COMPLETED',linear_relative_residual=float(la.norm(L@E-inc)/la.norm(inc)),
                           constitutive_relative_residual=float(la.norm(j-chi@E)/max(la.norm(j),1e-300)))
            return VectorState(self,eps,np.asarray(directions),np.asarray(polarizations),E,H,j,
                               -1j*self.k*j,p,E_exc,inc,L,fac,receipt)
        except Exception as exc:
            receipt.update(status='FAILED',error=repr(exc));raise
        finally:
            receipt.update(wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu)
            if receipt_path:
                path=Path(receipt_path);path.parent.mkdir(parents=True,exist_ok=True)
                path.write_text(json.dumps(receipt,indent=2,allow_nan=False))


def sphere_grid(radius,h):
    """Cell-center staircase sphere on cubic cells of side h."""
    a=np.arange(-int(np.ceil(radius/h)),int(np.ceil(radius/h))+1)*h
    p=np.stack(np.meshgrid(a,a,a,indexing='ij'),axis=-1).reshape(-1,3)
    return p[np.sum(p*p,axis=1)<=radius*radius*(1+1e-14)],h**3
