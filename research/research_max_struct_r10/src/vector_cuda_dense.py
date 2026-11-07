"""Thin retained-dense CUDA physical adapter; inherited vector FFT auxiliary.

Same R1/R9 point-DDA G, mature Torch LU/lu_solve, all60 RHS. No solver kernel
or new physical model. Primal API intentionally means only chi*Einc; adjoint
uses the cotangent. Material-dependent cache is reset between algorithm arms.
"""
import time,json
import numpy as np
import torch
from research.research_max_struct_r1.src.a9_green import build_goff
from .vector_cuda_model import ResidentVectorModel


class DenseResidentVectorModel(ResidentVectorModel):
    def __init__(self,cpu,gate):
        super().__init__(cpu,gate)
        wall,cput=time.perf_counter(),time.process_time()
        matrix=build_goff(cpu.points,cpu.k);matrix*=cpu.h**3
        matrix[np.diag_indices(self.M)]+=cpu.s
        self.G=torch.as_tensor(matrix,device=self.device,dtype=torch.complex128)
        torch.cuda.synchronize();del matrix
        self.dense_state=None;self.batch_records=[]
        self.counts.update(dense_G_assemblies=1,dense_G_uploads=1)
        self.setup.update(backend='Same CM+RR dyadic matrix; retained Torch CUDA default LU/all60 source+sameLU adjoint; inherited packed16 FFT auxiliary',
            dense_G_bytes=self.G.numel()*self.G.element_size(),dense_matrix_CPU_seconds=time.process_time()-cput,
            dense_matrix_wall_seconds=time.perf_counter()-wall,linalg_library=str(torch.backends.cuda.preferred_linalg_library()))

    def _factor(self,chi):
        key=chi.detach().cpu().numpy()
        if self.dense_state is not None and np.array_equal(key,self.dense_state['m']):
            self.counts['retained_dense_factor_hits']+=1;return self.dense_state
        # Release the obsolete material's factor before creating the next one.
        # No live optimizer/current history is reset by this storage reuse.
        self.dense_state=None;torch.cuda.empty_cache()
        wall,cput=time.perf_counter(),time.process_time()
        d=torch.repeat_interleave(chi,3)[:,None];K=-d*self.G;K.diagonal().add_(1)
        factor=torch.linalg.lu_factor_ex(K,check_errors=False)
        info=int(factor.info.item());del K
        self.counts['physical_LU']+=1;self.counts['physical_matrix_assemblies']+=1
        if info:raise RuntimeError('Qualified mature CUDA LU failed info='+str(info))
        self.dense_state=dict(m=key.copy(),LU=factor.LU,pivots=factor.pivots,J=None)
        self.batch_records.append(dict(operation='retained_CUDA_LU',LU=1,primal_RHS=0,adjoint_RHS=0,
            CPU_seconds=time.process_time()-cput,wall_seconds=time.perf_counter()-wall,additive=False))
        return self.dense_state

    def solve(self,chi,rhs,*,adjoint=False,label='',progress=None,checkpoint_dir=None):
        if checkpoint_dir is not None:raise ValueError('Physical dense adapter has no generation checkpoint replay')
        self.gate.wait('before_dense_'+label,torch.cuda.memory_reserved())
        chi=torch.as_tensor(chi,device=self.device,dtype=torch.float64)
        d=torch.repeat_interleave(chi,3)[:,None]
        if not adjoint and (rhs.shape!=self.inc.shape or not torch.equal(rhs,d*self.inc)):
            raise ValueError('Dense primal scope is specifically chi Einc, not arbitrary linear RHS')
        state=self._factor(chi);wall,cput=time.perf_counter(),time.process_time()
        if not adjoint and state['J'] is not None:
            self.counts['retained_dense_state_hits']+=1;return state['J']
        out=torch.linalg.lu_solve(state['LU'],state['pivots'],rhs,adjoint=adjoint)
        if adjoint:
            residual=out-self.green(d*out,True)-rhs
            self.counts['adjoint_RHS']+=rhs.shape[1];self.counts['physical_adjoint_RHS']+=rhs.shape[1]
        else:
            residual=out-d*(self.inc+self.green(out))
            self.counts['primal_RHS']+=rhs.shape[1];self.counts['physical_source_RHS']+=rhs.shape[1]
            state['J']=out
        relative=(torch.linalg.vector_norm(residual,dim=0)/torch.clamp(torch.linalg.vector_norm(rhs,dim=0),min=1e-300)).cpu().numpy()
        worst=float(relative.max());self.max_residual=max(self.max_residual,worst)
        rows=[dict(label=label,source=i,adjoint=adjoint,relative_residual=float(relative[i]),iterations=0,
            solver='same_retained_CUDA_dense_LU',CPU_seconds=None,wall_seconds=None,timing_scope='All60 batch, not per-source additive') for i in range(rhs.shape[1])]
        self.solve_records.extend(rows)
        if getattr(self,'progress_file',None):
            with self.progress_file.open('a') as handle:handle.writelines(json.dumps(r)+'\n' for r in rows)
        if progress:
            for row in rows:progress(row)
        self.batch_records.append(dict(operation='sameLU_adjoint_allRHS' if adjoint else 'allRHS_forward',
            LU=0,primal_RHS=0 if adjoint else rhs.shape[1],adjoint_RHS=rhs.shape[1] if adjoint else 0,
            CPU_seconds=time.process_time()-cput,wall_seconds=time.perf_counter()-wall,maximum_original_residual=worst,additive=False))
        if not np.isfinite(relative).all() or worst>self.rtol:raise RuntimeError('Dense original residual exceeds same frozen2e-6 tolerance')
        return out
