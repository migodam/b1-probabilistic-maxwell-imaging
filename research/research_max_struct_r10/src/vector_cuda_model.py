"""Resident R9 CUDA physics with cooperative boundary checks; no new solver."""
import time
import numpy as np
from scipy import linalg as la
from scipy.stats import chi2
import torch
from research.research_max_struct_r9.src.vector_original_cuda import CudaOriginalModel


class ResidentVectorModel(CudaOriginalModel):
    def __init__(self, cpu, gate):
        gate.wait('before_CUDA_context')
        super().__init__(cpu, fft_layout='packed', block_rhs=16)
        self.gate = gate
        gate.wait('after_shared_upload', torch.cuda.memory_reserved())

    def solve(self, chi, rhs, *, adjoint=False, label='', progress=None, checkpoint_dir=None):
        self.gate.wait('before_'+label, torch.cuda.memory_reserved())
        def checkpoint(row):
            if getattr(self,'progress_file',None):
                import json
                with self.progress_file.open('a') as handle:handle.write(json.dumps(row)+'\n')
            if progress: progress(row)
            if (row.get('source', -1)+1) % 20 == 0:
                self.gate.wait(label+'_source'+str(row['source']), torch.cuda.memory_reserved())
        return super().solve(chi, rhs, adjoint=adjoint, label=label,
            progress=checkpoint, checkpoint_dir=checkpoint_dir)


class ExactReceiverGeometry:
    """One actual economy SVD of S, not sqrt of a noisy Gram spectrum."""
    def __init__(self, model):
        wall, cpu = time.perf_counter(), time.process_time()
        self.U, self.s, self.Vh = la.svd(model.cpu.S, full_matrices=False, check_finite=False)
        threshold = np.finfo(float).eps * max(model.cpu.S.shape) * self.s[0]
        self.supported = int(np.sum(self.s > threshold))
        self.info = dict(algorithm='SciPy economy SVD of actual S',
            supported_rank=self.supported, numerical_threshold=float(threshold),
            singular_values=self.s.tolist(), vector_domain_dimension=model.M,
            data_range_dimension=180, CPU_seconds=time.process_time()-cpu,
            wall_seconds=time.perf_counter()-wall, observation_independent=True)
        model.counts['receiver_exact_SVD'] += 1


class NoiseReceiverSubspace:
    def __init__(self, model, geometry, Y, noise):
        wall, cpu = time.perf_counter(), time.process_time()
        if noise.get('model') != 'iid_circular_complex' or noise.get('whitening') not in (None, False):
            raise ValueError('Only raw iid CN scalar-sigma data supported; no silent normalization')
        sigma = np.asarray(noise['sigma'])
        if sigma.size != 1 or np.iscomplexobj(sigma) or not np.isfinite(sigma).all() or float(sigma)<0:
            raise ValueError('Scalar public complex noise sigma required')
        U, singular, Vh = geometry.U, geometry.s, geometry.Vh
        coordinate = U.conj().T @ Y
        outside = float(la.norm(Y-U@coordinate)**2)
        tail = outside + np.r_[np.cumsum(np.sum(abs(coordinate)**2, axis=1)[::-1])[::-1], 0.]
        radius2 = float(sigma)**2 * chi2.ppf(.99, 2*Y.size)/2
        legal = np.flatnonzero(tail[:geometry.supported+1] <= radius2)
        if not len(legal): raise ValueError('No numerically supported prefix meets public full noise ball')
        self.rank = int(legal[0]); r = self.rank
        self.V = model.array(Vh[:r].conj().T)
        self.Jd = model.array(Vh[:r].conj().T @ (coordinate[:r]/singular[:r,None]))
        self.Qd = model.inc + model.green(self.Jd)
        self.Dd = model.S @ self.Jd - model.array(Y)
        self.wd = 1/torch.sum(abs(model.array(Y))**2, dim=0)
        self.ws = 1/torch.sum(abs(model.inc)**2, dim=0)
        self.info = dict(rank=r, rule='public_iid_CN_full_noise_ball99', sigma=float(sigma),
            radius_squared=radius2, real_noise_dof=2*Y.size,
            fixed_010_diagnostic_rank=int(np.sum(singular >= .1*singular[0])),
            discarded_energy_squared_by_rank=tail.tolist(),
            supported_rank=geometry.supported, coordinate_system='complex3N vector current right domain',
            noise_scope='One full-noise event for all projectors; no per-r significance; DOES NOT bound crossgrid model error',
            CPU_seconds=time.process_time()-cpu, wall_seconds=time.perf_counter()-wall)
        model.counts['observation_noise_DPCS_fit'] += 1

    def project(self, X):
        return X-self.V@(self.V.conj().T@X)
