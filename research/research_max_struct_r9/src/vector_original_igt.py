"""Separate original3lambda pulse/collocation IGT adapter; frozen CM untouched.

Uses CellIntegratedGreen for self, every usable integer-lag block, external S,
and the exposed local A block. Incident dipoles remain center-collocated.
"""
import time
import numpy as np
from scipy import fft
from .cell_integrated_green import CellIntegratedGreen
from .vector_original_extension import OriginalVectorModel


class CellIntegratedOriginalModel(OriginalVectorModel):
    def __init__(self,side,*,source_indices=(59,),active_indices=None,generation_only=False,crop_generation_fft=False):
        wall,cput=time.perf_counter(),time.process_time()
        self.integrator=CellIntegratedGreen(3./side,2*np.pi,near_order=8,far_order=3,self_order=16,near_radius=2)
        super().__init__(side,active_indices=active_indices,generation_only=generation_only,crop_generation_fft=crop_generation_fft,
            gmres_rtol=2e-6,gmres_restart=60,gmres_cycles=20,block_rhs=4)
        self.source_indices=np.asarray(source_indices,int)
        if self.source_indices.ndim!=1 or np.any(self.source_indices<0) or np.any(self.source_indices>=60):raise ValueError('Original source indices0..59')
        self.inc=self.inc[:,self.source_indices].copy();self.P=len(self.source_indices)
        external_wall,external_cpu=time.perf_counter(),time.process_time()
        # Already includes source-cell volume. Do not multiply by h^3 again.
        self.S=self.integrator.external_electric_operator(self.points,self.receivers)
        self.counts.update(IGT_external_operator_builds=1,discarded_CM_point_external_operator_builds=1)
        external=dict(CPU_seconds=time.process_time()-external_cpu,wall_seconds=time.perf_counter()-external_wall,
            source_cells=self.N,receiver_points=60,receiver_vector_rows=180,quadrature_order=3,
            dyadic_quadrature_evaluations=self.N*60*3**3,bytes=self.S.nbytes)
        self.setup.update(CPU_seconds=time.process_time()-cput,wall_seconds=time.perf_counter()-wall,
            physical_model='Independent uniform-pulse/collocation source-cell integrated Green; not author software or continuum certificate',
            self_real=float(self.s.real),self_imag=float(self.s.imag),incident_bytes=self.inc.nbytes,receiver_bytes=self.S.nbytes,
            source_indices=self.source_indices.tolist(),stored_incident_source_fields=self.P,
            incident_convention='Point electric dipole field at cell center; no source averaging',
            IGT_kernel=self.integrator.receipt,IGT_external=external,
            local_A_block='Exact off-diagonal2x2x2 blocks from this same integrated integer-lag LUT',
            original_inverse_DOFS_unchanged=81000 if side==30 else None)

    def configure_embedding(self,crop=False):
        if crop and not self.generation_only:raise ValueError('Known-support crop generation-only')
        coordinates=np.array(np.unravel_index(self.active,self.shape)).T
        self.embedding_origin=coordinates.min(0) if crop else np.zeros(3,int)
        self.convolution_shape=tuple((coordinates.max(0)-self.embedding_origin+1).tolist()) if crop else self.shape
        self.convolution_N=int(np.prod(self.convolution_shape))
        embedding=np.ravel_multi_index((coordinates-self.embedding_origin).T,self.convolution_shape)
        self.green_rows=(3*embedding[:,None]+np.arange(3)).ravel()
        self.pad=tuple(fft.next_fast_len(2*d) for d in self.convolution_shape) if crop else tuple(2*d for d in self.shape)
        limits=np.asarray(self.convolution_shape)-1;self.integrator.build_lut(limits);self.s=self.integrator.self_density
        lagaxes=[np.where(np.arange(p)<=p//2,np.arange(p),np.arange(p)-p) for p in self.pad]
        lag=np.stack(np.meshgrid(*lagaxes,indexing='ij'),-1)
        usable=np.all(abs(lag)<=limits,axis=-1)
        # +/-shape Nyquist entries never represent an occupied-cell difference.
        # They remain exactly zero; no impossible displacement is integrated.
        kernel=np.zeros(self.pad+(3,3),complex);indices=lag[usable]+limits
        kernel[usable]=self.integrator.lut[indices[:,0],indices[:,1],indices[:,2]]
        self.kernel=fft.fftn(kernel,axes=(0,1,2));self.counts.update(Green_kernel_builds=1,IGT_integer_lag_LUT_builds=1)

    def local_offdiagonal_block(self):
        c=np.indices((2,2,2)).reshape(3,-1).T;d=c[:,None]-c[None,:]+self.integrator.limits
        block=self.integrator.lut[d[...,0],d[...,1],d[...,2]].copy()
        block[np.arange(8),np.arange(8)]=0
        return block.transpose(0,2,1,3).reshape(24,24)
