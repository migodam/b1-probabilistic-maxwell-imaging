"""R10 mature physical LU, only retaining CUDA allocator cache between m's."""
import time
import numpy as np
import torch
from research.research_max_struct_r10.src.vector_cuda_dense import DenseResidentVectorModel

class CachedDenseResidentVectorModel(DenseResidentVectorModel):
    def __init__(self,cpu,gate):
        super().__init__(cpu,gate)
        self.setup['R11_allocator_policy']='Release old dense_state; retain CUDA allocator cache; no per-material empty_cache'

    def _factor(self,chi):
        key=chi.detach().cpu().numpy()
        if self.dense_state is not None and np.array_equal(key,self.dense_state['m']):
            self.counts['retained_dense_factor_hits']+=1;return self.dense_state
        self.dense_state=None
        wall,cput=time.perf_counter(),time.process_time()
        d=torch.repeat_interleave(chi,3)[:,None];K=-d*self.G;K.diagonal().add_(1)
        factor=torch.linalg.lu_factor_ex(K,check_errors=False)
        info=int(factor.info.item());del K
        self.counts['physical_LU']+=1;self.counts['physical_matrix_assemblies']+=1
        self.counts['allocator_cache_retained_factor_changes']+=1
        if info:raise RuntimeError('Qualified mature CUDA LU failed info='+str(info))
        self.dense_state=dict(m=key.copy(),LU=factor.LU,pivots=factor.pivots,J=None)
        self.batch_records.append(dict(operation='retained_CUDA_LU',LU=1,primal_RHS=0,adjoint_RHS=0,
            CPU_seconds=time.process_time()-cput,wall_seconds=time.perf_counter()-wall,additive=False,allocator_cache_retained=True))
        return self.dense_state
