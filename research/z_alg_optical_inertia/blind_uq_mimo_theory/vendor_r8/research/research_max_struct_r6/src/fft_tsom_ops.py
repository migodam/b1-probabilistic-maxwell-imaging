"""Matrix-free native Fourier-TSOM wiring; no historical solver is modified.

The Fourier coefficients keep their native Euclidean metric. Receiver
projection remains low rank, reception remains dense, and the finite-domain
scalar-TM Green action uses an exact padded disk-pulse Toeplitz embedding.
"""
from collections import Counter
import numpy as np
from scipy.fft import fftn, ifftn, next_fast_len
from scipy.special import hankel1, jv


class NativeFourierAPCS:
    def __init__(self, shape, mf, Vd, cost=None):
        self.shape = tuple(map(int, shape))
        mf=int(mf)
        if len(self.shape) != 2 or not 1 <= int(mf) <= min(self.shape)//2:
            raise ValueError('Two-dimensional grid and legal four-corner mask required')
        axes = [np.unique(np.r_[np.arange(mf), np.arange(s-mf, s)]) for s in self.shape]
        self.frequencies = np.array(np.meshgrid(*axes,indexing='ij')).reshape(2,-1).T
        self.Vd = np.asarray(Vd,complex)
        self.n = int(np.prod(self.shape)); self.k = len(self.frequencies)
        if self.Vd.ndim != 2 or self.Vd.shape[0] != self.n: raise ValueError('DPCS grid mismatch')
        self.cost = Counter() if cost is None else cost

    def project(self, x):
        self.cost['DPCS_projection_RHS'] += x.shape[1]
        return x-self.Vd@(self.Vd.conj().T@x)

    def apply(self, coefficients):
        a=np.asarray(coefficients,complex)
        if a.ndim != 2 or a.shape[0] != self.k: raise ValueError('Expected mask coefficients by source')
        packed=np.zeros(self.shape+(a.shape[1],),complex)
        packed[tuple(self.frequencies.T)] = a
        self.cost['coefficient_IFFT_calls'] += 1
        self.cost['coefficient_IFFT_2D_transforms'] += a.shape[1]
        field=ifftn(packed,axes=(0,1),norm='ortho').reshape(self.n,-1)
        return self.project(field)

    def adjoint(self, current):
        z=np.asarray(current,complex)
        if z.ndim != 2 or z.shape[0] != self.n: raise ValueError('Expected cell by source current')
        packed=self.project(z).reshape(self.shape+(z.shape[1],))
        self.cost['coefficient_FFT_calls'] += 1
        self.cost['coefficient_FFT_2D_transforms'] += z.shape[1]
        spectrum=fftn(packed,axes=(0,1),norm='ortho')
        return spectrum[tuple(self.frequencies.T)]

    def transport(self, old, old_frequencies):
        """Exact nested-mask zero padding without changing coefficient metric."""
        lookup={tuple(f): i for i,f in enumerate(old_frequencies)}
        out=np.zeros((self.k,old.shape[1]),complex)
        for i,f in enumerate(self.frequencies):
            if tuple(f) in lookup: out[i]=old[lookup[tuple(f)]]
        return out


class DiskPulseToeplitzGreen:
    def __init__(self, centers, shape, radius, k, cost=None):
        self.shape=tuple(map(int,shape)); self.n=int(np.prod(self.shape))
        if len(self.shape)!=2 or min(self.shape)<2: raise ValueError('Uniform two-dimensional grid required')
        coords=np.asarray(centers,float).reshape(self.shape+(2,))
        v0=coords[1,0]-coords[0,0]; v1=coords[0,1]-coords[0,0]
        ii,jj=np.indices(self.shape)
        expected=coords[0,0]+ii[...,None]*v0+jj[...,None]*v1
        if not np.allclose(coords,expected,rtol=1e-12,atol=1e-13):
            raise ValueError('Centers must be a full affine uniform grid in native flattening order')
        if abs(np.dot(v0,v1))>1e-12*np.linalg.norm(v0)*np.linalg.norm(v1):
            raise ValueError('Cartesian grid required')
        if float(radius)<=0 or float(k)<=0 or min(np.linalg.norm(v0),np.linalg.norm(v1))<=2*radius:
            raise ValueError('Positive wavenumber/nonoverlapping pulse disks required')
        self.cost=Counter() if cost is None else cost
        self.pad=tuple(next_fast_len(2*s-1) for s in self.shape)
        d0=np.arange(-(self.shape[0]-1),self.shape[0]); d1=np.arange(-(self.shape[1]-1),self.shape[1])
        a,b=np.meshgrid(d0,d1,indexing='ij')
        dist=np.linalg.norm(a[...,None]*v0+b[...,None]*v1,axis=-1)
        zero=(a==0)&(b==0); safe=np.where(zero,1.,dist)
        kernel=.5j*np.pi*k*radius*jv(1,k*radius)*hankel1(0,k*safe)
        kernel[zero]=.5j*np.pi*k*radius*hankel1(1,k*radius)-1
        embedded=np.zeros(self.pad,complex)
        embedded[a%self.pad[0],b%self.pad[1]]=kernel
        self.spectrum=fftn(embedded)
        self.cost['Green_kernel_FFT_setup']+=1
        self.cost['Green_nonperiodic_lags_generated']+=kernel.size

    def apply(self, current, adjoint=False):
        x=np.asarray(current,complex)
        if x.ndim!=2 or x.shape[0]!=self.n: raise ValueError('Expected cell by source current')
        work=np.zeros(self.pad+(x.shape[1],),complex)
        work[:self.shape[0],:self.shape[1]]=x.reshape(self.shape+(x.shape[1],))
        spectrum=self.spectrum.conj() if adjoint else self.spectrum
        out=ifftn(fftn(work,axes=(0,1))*spectrum[...,None],axes=(0,1))
        self.cost['Green_adjoint_calls' if adjoint else 'Green_forward_calls']+=1
        self.cost['Green_FFT_2D_transforms']+=x.shape[1]
        self.cost['Green_IFFT_2D_transforms']+=x.shape[1]
        return out[:self.shape[0],:self.shape[1]].reshape(self.n,-1)

    def adjoint(self,current): return self.apply(current,adjoint=True)


class FixedCurrentBlock:
    """Composition interface for a fixed-material residual/readout plugin.

    transform and transform_adjoint act only in cell/current space. A can
    supply its exact fixed-m P and P*, independently of the FFT propagation.
    readout receives J and Q from the updated current, then owns material
    bounds/prior/derivative handling; this class never substitutes its prior.
    """
    def __init__(self,basis,green,S,E,Y,JD,m,Od,Os,rho=1.,transform=None,transform_adjoint=None):
        self.U,self.G=basis,green
        self.S,self.E,self.Y,self.JD=[np.asarray(x,complex) for x in (S,E,Y,JD)]
        self.m=np.asarray(m,float); self.Od=np.asarray(Od,complex); self.Os=np.asarray(Os,complex)
        self.rho=float(rho)
        if (transform is None)!=(transform_adjoint is None): raise ValueError('Both residual actions required')
        self.P=transform or (lambda x:x); self.PH=transform_adjoint or (lambda x:x)

    def evaluate(self,A):
        J=self.JD+self.U.apply(A); Q=self.E+self.G.apply(J)
        D=self.S@J-self.Y; R=self.P(J-self.m[:,None]*Q)
        weighted=self.PH(R@self.Os)
        back=self.S.conj().T@(D@self.Od)+self.rho*(weighted-self.G.adjoint(self.m[:,None]*weighted))
        grad=self.U.adjoint(back)
        value=.5*np.vdot(D,D@self.Od).real+.5*self.rho*np.vdot(R,R@self.Os).real
        return float(value),grad,J,Q

    def step(self,A,old_gradient=None,old_direction=None,readout=None):
        """One global PR+ step with exact fixed-material quadratic line search."""
        value,grad,_,_=self.evaluate(A)
        pr=0. if old_gradient is None else max(0.,float(np.vdot(grad,grad-old_gradient).real/max(np.vdot(old_gradient,old_gradient).real,1e-300)))
        direction=-grad if old_gradient is None else -grad+pr*old_direction
        if np.vdot(grad,direction).real>=0: direction=-grad; pr=0.
        dJ=self.U.apply(direction); dQ=self.G.apply(dJ)
        dD=self.S@dJ; dR=self.P(dJ-self.m[:,None]*dQ)
        denom=float(np.vdot(dD,dD@self.Od).real+self.rho*np.vdot(dR,dR@self.Os).real)
        alpha=-float(np.vdot(grad,direction).real)/denom if denom>1e-300 else 0.
        newA=A+alpha*direction
        newvalue,_,J,Q=self.evaluate(newA)
        material=None if readout is None else readout(J,Q)
        return dict(coefficients=newA,current=J,total_field=Q,gradient=grad,direction=direction,
                    alpha=alpha,pr=pr,before=value,after=newvalue,material=material)
