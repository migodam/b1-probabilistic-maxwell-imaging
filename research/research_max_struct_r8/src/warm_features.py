"""Observation-only Born and DPCS features, independent of material truth."""
import time
import numpy as np
from scipy.linalg import svd

FEATURE_NAMES=['born_real','born_imag','dpcs_readout_real','dpcs_readout_imag',
               'log_current_energy','log_total_field_ratio','source_aperture','sources_over_12',
               'x_over_width','y_over_width']

class FeatureOperator:
    def __init__(self,public):
        t,c=time.perf_counter(),time.process_time()
        self.G=np.asarray(public['G']);self.S=np.asarray(public['S']);self.E=np.asarray(public['E0'])
        self.shape=tuple(public['grid_shape']);self.n=self.E.shape[0]
        U,s,Vh=svd(self.S,full_matrices=False,check_finite=False)
        rank=int(np.count_nonzero(s>1e-3*s[0]));self.pinv=(Vh[:rank].conj().T/s[:rank])@U[:,:rank].conj().T
        self.diag=np.sum(abs(self.S)**2,axis=0)*np.sum(abs(self.E)**2,axis=1)
        self.energy=np.mean(abs(self.E)**2,axis=1)
        self.coords=np.asarray(public.get('centers',public.get('coords')))
        self.aperture=float(public.get('source_aperture_fraction',1.))
        self.setup=dict(wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c,
                        receiver_SVD=1,dpcs_rank=rank,LU=0,complex_RHS=0)

    def __call__(self,Y):
        t,c=time.perf_counter(),time.process_time()
        born=np.sum(self.E.conj()*(self.S.conj().T@Y),axis=1)/self.diag
        J=self.pinv@Y;Q=self.E+self.G@J
        read=np.sum(Q.conj()*J,axis=1)/np.maximum(np.sum(abs(Q)**2,axis=1),1e-20)
        je=np.mean(abs(J)**2,axis=1)/np.maximum(self.energy,1e-20)
        qe=np.mean(abs(Q)**2,axis=1)/np.maximum(self.energy,1e-20)
        # Bounded robust encoding is a fixed public transform, never fitted on test labels.
        x=[np.arcsinh(born.real),np.arcsinh(born.imag),np.arcsinh(read.real),np.arcsinh(read.imag),
           np.log1p(je),np.log(np.maximum(qe,1e-12)),np.full(self.n,self.aperture),
           np.full(self.n,self.E.shape[1]/12),self.coords[:,0]/.48,self.coords[:,1]/.48]
        features=np.stack(x).reshape(len(x),*self.shape).astype(np.float32)
        return features,dict(born=np.clip(born.real,.02,20.),dpcs=np.clip(read.real,.02,20.)),dict(
            wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c,
            Green_actions_columns=self.E.shape[1],receiver_adjoint_columns=self.E.shape[1],
            DPCS_pseudoinverse_columns=self.E.shape[1],LU=0,complex_RHS=0,truth_access='NONE')
