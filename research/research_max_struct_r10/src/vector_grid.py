"""R10 geometry adapter of the frozen R9 full-vector CM+RR FFT operator.

No material generation, solve, CUDA allocation, or truth access on import.
Only grid width and explicitly declared external acquisition are generalized.
"""
from collections import Counter
import time
import numpy as np
from research.research_max_struct_r1.src.a9_green import dyadic_green
from research.research_max_struct_r9.src.vector_original_extension import OriginalVectorModel


def grid_points(side, width):
    ijk = np.indices((side,) * 3).reshape(3, -1).T
    return (ijk - (side - 1) / 2) * (width / side)


def external_geometry(contract):
    """Three circles, independent of inverse/generation resolution."""
    angle = 2 * np.pi * np.arange(20) / 20
    circles, polarizations = [], []
    for axis, center in zip((2, 0, 1), contract['circle_centers']):
        q = np.zeros((20, 3)); transverse = [a for a in range(3) if a != axis]
        q[:, transverse[0]] = contract['circle_radius'] * np.cos(angle)
        q[:, transverse[1]] = contract['circle_radius'] * np.sin(angle)
        q += np.asarray(center)
        pol = np.zeros_like(q); pol[:, axis] = 1
        circles.append(q); polarizations.append(pol)
    return np.concatenate(circles), np.concatenate(polarizations)


class VectorGrid(OriginalVectorModel):
    """Inherits R9 embedding, Green/adjoint and qualified CPU GMRES unchanged.

    Compatible with the existing CudaOriginalModel resident packed16 wrapper.
    All material cells remain unknown; no known-support restriction is exposed.
    """
    def __init__(self, side, width, acquisition, *, gmres_rtol=2e-6,
                 gmres_restart=60, gmres_cycles=20):
        wall, cpu = time.perf_counter(), time.process_time()
        if side % 2 or side < 4 or width <= 0:
            raise ValueError('Positive width and even grid required')
        self.side = int(side); self.width = float(width); self.h = width / side
        self.k = 2 * np.pi; self.shape = (side,) * 3; self.pad = (2 * side,) * 3
        self.full_N = side ** 3; self.active = np.arange(self.full_N)
        self.generation_only = False; self.N = self.full_N; self.M = 3 * self.N
        self.row_indices = np.arange(self.M); self.points = grid_points(side, width)
        self.P = 60; self.nr = 180; self.block_rhs = 4
        self.rtol = gmres_rtol; self.restart = gmres_restart; self.cycles = gmres_cycles
        self.s = -1 / 3 + 1j * self.k ** 3 * self.h ** 3 / (6 * np.pi)
        self.counts = Counter(); self.solve_records = []; self.max_residual = 0.
        self._physical_cache = None
        self.sources, self.polarizations = external_geometry(acquisition)
        self.receivers = self.sources.copy()
        if np.any(np.max(abs(self.sources), axis=1) <= width / 2):
            raise ValueError('A source/receiver entered unknown-material box')
        self.configure_embedding(False)
        self.inc = np.empty((self.M, 60), complex)
        self.S = np.empty((180, self.M), complex)
        for i in range(0, 60, self.block_rhs):
            j = min(i + self.block_rhs, 60)
            g = dyadic_green(self.points[:, None] - self.sources[None, i:j], self.k)
            self.inc[:, i:j] = np.einsum('nspq,sq->nps', g, self.polarizations[i:j]).reshape(self.M, j-i)
            g = dyadic_green(self.receivers[i:j, None] - self.points[None], self.k)
            self.S[3*i:3*j] = self.h ** 3 * g.transpose(0, 2, 1, 3).reshape(3*(j-i), self.M)
        self.counts.update(incident_dipole_fields=60, receiver_vector_rows=180)
        self.setup = dict(CPU_seconds=time.process_time()-cpu,
            wall_seconds=time.perf_counter()-wall, full_cells=self.full_N,
            stored_current_cells=self.N, complex_current_components=self.M,
            known_support_restriction=False, block_rhs=4,
            kernel_bytes=self.kernel.nbytes, incident_bytes=self.inc.nbytes,
            receiver_bytes=self.S.nbytes, geometry_width=width,
            source_radius=acquisition['circle_radius'], source_centers=acquisition['circle_centers'],
            source_geometry_resolution_independent=True,
            self_real=float(self.s.real), self_imag=float(self.s.imag),
            convolution_shape=list(self.shape), convolution_padding=list(self.pad),
            generation_fft_crop=False, embedding_origin_indices=[0, 0, 0],
            physical_model='Inherited R1/R8/R9 dyadic CM+RR point DDA; not new physics or author reproduction',
            implementation_provenance='R10 geometry-only constructor; inherited OriginalVectorModel.configure_embedding/green/solve')
