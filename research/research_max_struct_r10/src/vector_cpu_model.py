"""Prepared CPU fallback for the same vector model; NOT RUN by GPU queue.

Uses the existing R9 SciPy qualified GMRES/FFT on CPU. Host Torch arrays wrap
NumPy values without a CUDA context. Current/readout/budget code is shared.
CPU MGS and CUDA CGS2 paths need not be bitwise identical; residual is qualified.
"""
import numpy as np
import torch
import time
from scipy import linalg as la


class NoCudaGate:
    helper_CPU_seconds=0.;query_count=0;unmeasured_helper_invocations=0
    def wait(self,boundary,memory_bytes=None):
        return dict(action='CPU_ONLY_NO_CUDA_WORK',boundary=boundary,own_memory_bytes=0)


class CpuVectorModel:
    def __init__(self,cpu):
        self.cpu=cpu;self.device=torch.device('cpu');self.N=cpu.N;self.M=cpu.M
        self.P=cpu.P;self.side=cpu.side;self.s=cpu.s;self.block_rhs=16
        self.fft_layout='packed';self.inc=self.array(cpu.inc);self.S=self.array(cpu.S)
        self.counts=cpu.counts;self.solve_records=cpu.solve_records;self.gate=NoCudaGate()
        self.setup=dict(cpu.setup,backend='Inherited CPU scipy GMRES/MGS + scipy padded dyadic FFT; CPU Torch material/current arrays')

    @staticmethod
    def array(value):return torch.as_tensor(np.asarray(value),dtype=torch.complex128)

    def green(self,value,adjoint=False):
        return torch.from_numpy(self.cpu.green(value.detach().numpy(),adjoint))

    def solve(self,chi,rhs,*,adjoint=False,label='',progress=None,checkpoint_dir=None):
        begin=len(self.solve_records)
        out=self.cpu.solve(chi.detach().numpy(),rhs.detach().numpy(),adjoint=adjoint,label=label,progress=progress)
        if any(row.get('relative_residual',float('inf'))>self.cpu.rtol for row in self.solve_records[begin:]):
            raise RuntimeError('CPU original residual exceeded the same frozen tolerance; no silent tolerance change')
        return torch.from_numpy(out)


class CpuDenseVectorModel(CpuVectorModel):
    """R9 paid_density_state retained LU/all60 RHS, same CM+RR grid physics.

    Auxiliary current propagation stays on the inherited padded FFT. The dense
    matrix is the existing R1 A9 dyadic primitive with the same volume/self.
    No fresh factorization occurs for a gradient at the same paid material.
    """
    def __init__(self,cpu):
        super().__init__(cpu)
        from research.research_max_struct_r1.src.a9_green import build_goff
        wall,cput=time.perf_counter(),time.process_time()
        self.G=build_goff(cpu.points,cpu.k);self.G*=cpu.h**3
        self.G[np.diag_indices(self.M)]+=cpu.s
        self.counts['dense_G_assemblies']+=1;self.dense_state=None
        self.setup.update(backend='R9 paid_density_state dense SciPy LU/all60RHS with retained adjoint; inherited FFT auxiliary propagation',
            dense_G_bytes=self.G.nbytes,dense_assembly_CPU_seconds=time.process_time()-cput,
            dense_assembly_wall_seconds=time.perf_counter()-wall)

    def state(self,m):
        from research.research_max_struct_r9.src.vector_extension import paid_density_state
        m=np.asarray(m,float)
        if self.dense_state is not None and np.array_equal(m,self.dense_state['m']):
            self.counts['retained_dense_state_hits']+=1;return self.dense_state
        before=self.counts['physical_LU'];wall,cput=time.perf_counter(),time.process_time()
        self.dense_state=paid_density_state(self.G,self.cpu.inc,m,self.cpu.S,self.counts)
        self.counts['primal_RHS']+=self.P
        st=self.dense_state;residual=st['J']-np.repeat(m,3)[:,None]*st['E']
        denom=la.norm(np.repeat(m,3)[:,None]*self.cpu.inc,axis=0)
        relative=la.norm(residual,axis=0)/np.maximum(denom,1e-300)
        self.last_dense_batch=dict(operation='shared_LU_and_all_source_forward',CPU_seconds=time.process_time()-cput,
            wall_seconds=time.perf_counter()-wall,LU=int(self.counts['physical_LU']-before),primal_RHS=self.P,
            maximum_relative_residual=float(relative.max()),timing_scope='One batch, not per-source additive')
        for source in range(self.P):self.solve_records.append(dict(label='dense_all60_forward',source=source,adjoint=False,
            relative_residual=float(relative[source]),solver='retained_dense_LU',iterations=0,
            CPU_seconds=None,wall_seconds=None,timing_scope='Included in shared all60 batch'))
        return st

    def solve(self,chi,rhs,*,adjoint=False,label='',progress=None,checkpoint_dir=None):
        """Primal scope is specifically L^-1 (chi Einc), not arbitrary RHS.

        All actual calls here are the physical state/generator at that fixed
        RHS. The adjoint accepts the loss cotangent and reuses the same LU.
        """
        m=chi.detach().numpy();st=self.state(m)
        if not adjoint:
            expected=torch.repeat_interleave(chi,3)[:,None]*self.inc
            if rhs.shape!=expected.shape or not torch.equal(rhs,expected):
                raise ValueError('Dense primal adapter supports only frozen chi Einc physical RHS')
            return torch.from_numpy(st['J'])
        wall,cput=time.perf_counter(),time.process_time()
        out=la.lu_solve(st['factor'],rhs.detach().numpy(),trans=2,check_finite=False)
        self.counts['physical_adjoint_RHS']+=rhs.shape[1];self.counts['adjoint_RHS']+=rhs.shape[1]
        d=np.repeat(m,3)[:,None]
        # One qualified residual of the actual adjoint batch; G*D placement retained.
        residual=out-self.G.conj().T@(d*out)-rhs.detach().numpy()
        relative=la.norm(residual,axis=0)/np.maximum(la.norm(rhs.detach().numpy(),axis=0),1e-300)
        self.last_dense_adjoint_batch=dict(operation='all_source_adjoint_same_LU',CPU_seconds=time.process_time()-cput,
            wall_seconds=time.perf_counter()-wall,LU=0,adjoint_RHS=rhs.shape[1],maximum_relative_residual=float(relative.max()))
        for source in range(rhs.shape[1]):self.solve_records.append(dict(label=label,source=source,adjoint=True,
            relative_residual=float(relative[source]),solver='same_retained_dense_LU',iterations=0,
            CPU_seconds=None,wall_seconds=None,timing_scope='Included in shared all60 adjoint batch'))
        if relative.max()>2e-6:raise RuntimeError('Dense adjoint residual violates frozen tolerance')
        return torch.from_numpy(out)
