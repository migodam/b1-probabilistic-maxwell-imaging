"""Uniform pulse/collocation electric VIE cube integrals, additive R9 backend.

Standard fixed quadrature, not a new theorem or passivity certificate. The
dyadic_green integrand already includes k². No extra volume factor is applied
to the returned G or external operator. Point incident fields are collocated.
"""
import time
import numpy as np
from numpy.polynomial.legendre import leggauss
from research.research_max_struct_r1.src.a9_green import dyadic_green


def cube_rule(h,order):
    nodes,weights=leggauss(int(order));a=float(h)/2
    q=np.stack(np.meshgrid(nodes,nodes,nodes,indexing='ij'),-1).reshape(-1,3)*a
    w=np.prod(np.stack(np.meshgrid(weights,weights,weights,indexing='ij'),-1),axis=-1).ravel()*a**3
    return q,w


def self_scalar(h,k,order=16):
    """PV/delta cubic symmetry: -1/3 + k²/(6π) ∫cube exp(ikr)/r.

    Eight octants and three max-coordinate parts; no singular quadrature node.
    Includes the full complex integral, including the dynamic real self term.
    """
    q,w=leggauss(int(order));q=(q+1)/2;w=w/2;a=float(h)/2
    t,u,v=np.meshgrid(q,q,q,indexing='ij');weights=np.prod(np.stack(np.meshgrid(w,w,w,indexing='ij'),-1),axis=-1)
    root=np.sqrt(1+u*u+v*v)
    integral=24*a*a*np.sum(weights*t*np.exp(1j*k*a*t*root)/root)
    return complex(-1/3+k*k/(6*np.pi)*integral)


def integrated_blocks(displacements,h,k,order,chunk_lags=256):
    """Offself ∫cube dyadic_green(d-s,k) ds at fixed arbitrary d.

    Caller must exclude self/interior receivers. Chunk size changes memory only.
    """
    d=np.asarray(displacements,float).reshape(-1,3);q,w=cube_rule(h,order)
    result=np.empty((len(d),3,3),complex)
    chunk=max(1,min(int(chunk_lags),int(150000/max(len(q),1))))
    for start in range(0,len(d),chunk):
        stop=min(start+chunk,len(d));blocks=dyadic_green(d[start:stop,None,:]-q[None,:,:],k)
        result[start:stop]=np.einsum('lqij,q->lij',blocks,w)
    return result


class CellIntegratedGreen:
    """Cached integer-lag LUT and dense/receiver interfaces for uniform grids.

    Fixed primary orders8/3/16 and near Chebyshev radius2 are protocol choices.
    No material, observations, losses or GT enter operator construction.
    """
    def __init__(self,h,k=2*np.pi,near_order=8,far_order=3,self_order=16,near_radius=2):
        self.h=float(h);self.k=float(k);self.near_order=int(near_order);self.far_order=int(far_order)
        self.self_order=int(self_order);self.near_radius=int(near_radius)
        if self.h<=0 or self.k<=0:raise ValueError('Positive h/k required')
        self.self_density=self_scalar(self.h,self.k,self.self_order)
        self.lut=None;self.limits=None;self.receipt={}

    def build_lut(self,limits):
        wall,cpu=time.perf_counter(),time.process_time();limits=np.asarray(limits,int)
        if limits.shape!=(3,) or np.any(limits<0):raise ValueError('Three nonnegative integer lag limits required')
        shape=2*limits+1;offsets=np.indices(tuple(shape)).reshape(3,-1).T-limits
        norm=np.max(abs(offsets),axis=1);zero=norm==0;near=(norm>0)&(norm<=self.near_radius);far=norm>self.near_radius
        blocks=np.empty((len(offsets),3,3),complex);blocks[zero]=self.self_density*np.eye(3)
        blocks[near]=integrated_blocks(offsets[near]*self.h,self.h,self.k,self.near_order)
        blocks[far]=integrated_blocks(offsets[far]*self.h,self.h,self.k,self.far_order)
        self.lut=blocks.reshape(tuple(shape)+(3,3));self.limits=limits
        self.receipt=dict(kind='pulse_collocation_fixed_cube_integral',h=self.h,k=self.k,near_radius=self.near_radius,
            near_order=self.near_order,far_order=self.far_order,self_order=self.self_order,
            unique_integer_lags=len(offsets),near_nonself_lags=int(near.sum()),far_lags=int(far.sum()),
            dyadic_quadrature_evaluations=int(near.sum())*self.near_order**3+int(far.sum())*self.far_order**3,
            self_density=[self.self_density.real,self.self_density.imag],LUT_bytes=self.lut.nbytes,
            LUT_CPU_seconds=time.process_time()-cpu,LUT_wall_seconds=time.perf_counter()-wall,
            no_material_or_data_dependency=True,continuum_passivity_certificate='NOT_CLAIMED')
        return self.lut

    def assemble_dense(self,integer_coordinates,row_chunk=64):
        wall,cpu=time.perf_counter(),time.process_time();raw=np.asarray(integer_coordinates)
        if raw.ndim!=2 or raw.shape[1]!=3 or not np.array_equal(raw,np.rint(raw)):raise ValueError('Exact uniform-grid integer coordinates required')
        coordinates=np.asarray(raw,int)
        if len(np.unique(coordinates,axis=0))!=len(coordinates):raise ValueError('Duplicate cells')
        limits=coordinates.max(0)-coordinates.min(0)
        if self.lut is None or not np.array_equal(self.limits,limits):self.build_lut(limits)
        n=len(coordinates);G=np.empty((3*n,3*n),complex)
        for start in range(0,n,row_chunk):
            stop=min(start+row_chunk,n);delta=coordinates[start:stop,None,:]-coordinates[None,:,:]+self.limits
            b=self.lut[delta[...,0],delta[...,1],delta[...,2]]
            G[3*start:3*stop]=b.transpose(0,2,1,3).reshape(3*(stop-start),3*n)
        self.receipt.update(active_cells=n,vector_dimension=3*n,dense_G_bytes=G.nbytes,
            assemble_total_CPU_seconds=time.process_time()-cpu,assemble_total_wall_seconds=time.perf_counter()-wall)
        return G

    def external_electric_operator(self,points,receivers,row_chunk=16):
        """Uniform-source-cell integral, fixed Gauss3; adjoint is S.conj().T.

        Background receivers must lie strictly outside all occupied source cubes.
        """
        p=np.asarray(points,float).reshape(-1,3);r=np.asarray(receivers,float).reshape(-1,3)
        if np.any(np.max(abs(r[:,None,:]-p[None,:,:]),axis=-1)<=self.h/2):raise ValueError('Receiver inside a source cube')
        S=np.empty((3*len(r),3*len(p)),complex)
        for start in range(0,len(r),row_chunk):
            stop=min(start+row_chunk,len(r));d=r[start:stop,None,:]-p[None,:,:]
            b=integrated_blocks(d.reshape(-1,3),self.h,self.k,self.far_order).reshape(stop-start,len(p),3,3)
            S[3*start:3*stop]=b.transpose(0,2,1,3).reshape(3*(stop-start),3*len(p))
        return S

    @staticmethod
    def receiver_adjoint(S,residual):return np.asarray(S).conj().T@residual

    def far_field(self,points,current,directions):
        """Uniform-cell current far amplitude, including exact cube form factor."""
        p=np.asarray(points,float).reshape(-1,3);d=np.asarray(directions,float).reshape(-1,3);j=np.asarray(current,complex)
        if j.ndim==1:j=j[:,None]
        if j.shape[0]!=3*len(p) or np.any(abs(np.linalg.norm(d,axis=1)-1)>1e-12):raise ValueError('Invalid current/directions')
        phase=np.exp(-1j*self.k*d@p.T)
        dip=np.einsum('rn,ncp->rcp',phase,self.h**3*j.reshape(len(p),3,-1))
        factor=np.prod(np.sinc(self.k*self.h*d/(2*np.pi)),axis=1)
        transverse=dip-d[:,:,None]*np.einsum('rc,rcp->rp',d,dip)[:,None,:]
        return self.k**2/(4*np.pi)*factor[:,None,None]*transverse
