"""DCTMC-style development control; known lifting, not a new algorithm.

Levinson--Markel, PRE 94,043317 (2016), arXiv:1401.3319v3, Eqs.10,11,23,24,26.
This variation starts at a public host, extracts real diagonal V into a box,
and truncates product singular coordinates at a declared relative noise-level
cutoff. It does not reproduce the paper's numerical experiments or weighted
off-diagonal force-diagonalization. No symmetry is imposed on intermediate T/V.
"""
import warnings
import numpy as np
from scipy import linalg as la


def relative_norm(a, reference):
    return float(la.norm(a) / max(float(la.norm(reference)), np.finfo(float).tiny))


class TMatrixCompletionControl:
    def __init__(self, model, Y, product_rtol=1e-3, box=(.02, 20.)):
        self.model, self.Y = model, np.asarray(Y, complex)
        self.box = tuple(map(float, box))
        self.product_rtol = float(product_rtol)
        self.eye = np.eye(model.n, dtype=complex)
        self.operation = model._physics.operation
        with self.operation('TMATRIX_DATA_SVD', receiver_SVD=1, incident_SVD=1):
            self.Ua, self.sa, vha = la.svd(model.S, full_matrices=False, check_finite=True)
            self.Ub, self.sb, vhb = la.svd(model.E, full_matrices=False, check_finite=True)
            self.Va, self.Vb = vha.conj().T, vhb.conj().T
            self.Ysv = self.Ua.conj().T @ self.Y @ self.Vb
            self.products = self.sa[:, None] * self.sb[None, :]
            self.cutoff = self.product_rtol * float(self.products.max())
            self.mask = self.products > self.cutoff
            if not np.any(self.mask):
                raise ValueError('No product singular coordinate retained')
            self.known = np.zeros_like(self.Ysv)
            self.known[self.mask] = self.Ysv[self.mask] / self.products[self.mask]
        self.data_info = dict(
            receiver_singular_values=self.sa.tolist(), incident_singular_values=self.sb.tolist(),
            product_rtol=self.product_rtol, product_absolute_cutoff=self.cutoff,
            retained_complex_coordinates=int(self.mask.sum()), total_measured_complex_coordinates=int(self.mask.size),
            unconstrained_T_complex_coordinates=int(model.n**2-self.mask.sum()),
            discarded_observed_data_relative_norm=float(la.norm(self.Ysv[~self.mask])/la.norm(self.Y)),
            cutoff_semantics='DEV noise-level truncation, not certified numerical rank or noise-optimal threshold',
            retained_mask=self.mask.tolist(),
            reciprocal_green_relative_transpose_defect=relative_norm(model.G-model.G.T, model.G),
            symmetry_policy='No T/V symmetrization; physical local real diagonal material only')

    def physical(self, material):
        """One original full Maxwell candidate and all n columns of its T."""
        m = np.asarray(material, float)
        state = self.model.state(m)
        with self.operation('TMATRIX_PHYSICAL_DENSE_T', rhs=self.model.n,
                            T_dense_RHS=self.model.n, T_dense_reconstructions=1,
                            retained_LU_dense_solve_dimension=self.model.n):
            T = la.lu_solve(state['lu'], np.diag(m).astype(complex), check_finite=False)
        if not np.all(np.isfinite(T)):
            raise FloatingPointError('Nonfinite physical T')
        with self.operation('TMATRIX_PHYSICAL_SCORE', physical_candidate_scores=1):
            residual = state['Y']-self.Y
            prior = self.model.reg.evaluate(m)[0]
            score = float(.5*la.norm(residual)**2/la.norm(self.Y)**2+prior)
        row = dict(physical_Phi=score, physical_relative_residual=relative_norm(residual,self.Y),
                   prior=float(prior), material_min=float(m.min()), material_max=float(m.max()),
                   material_lower_count=int(np.count_nonzero(m == self.box[0])),
                   material_upper_count=int(np.count_nonzero(m == self.box[1])),
                   physical_T_relative_transpose_defect=relative_norm(T-T.T,T))
        return T, state, row

    def overwrite(self, T):
        """Exact retained-coordinate overwrite, leaving other coordinates alone."""
        with self.operation('TMATRIX_DATA_OVERWRITE', data_overwrites=1,
                            T_to_singular_products=2, singular_to_T_products=2,
                            completion_data_products=2):
            coordinates = self.Va.conj().T @ T @ self.Ub
            delta = np.zeros_like(coordinates)
            delta[self.mask] = self.known[self.mask]-coordinates[self.mask]
            completed = T + self.Va @ delta @ self.Ub.conj().T
            if not np.all(np.isfinite(completed)):
                raise FloatingPointError('Nonfinite overwritten T')
            prediction = self.model.S @ (completed @ self.model.E)
            sv_residual = self.Ua.conj().T @ (prediction-self.Y) @ self.Vb
            row = dict(completed_full_data_relative_residual=relative_norm(prediction-self.Y,self.Y),
                       completed_retained_data_relative_residual=float(la.norm(sv_residual[self.mask])/la.norm(self.Y)),
                       completed_discarded_data_relative_residual=float(la.norm(sv_residual[~self.mask])/la.norm(self.Y)),
                       overwrite_relative_change=relative_norm(completed-T,T),
                       completed_T_relative_transpose_defect=relative_norm(completed-completed.T,completed))
        return completed, row

    def localize(self, T):
        """Full n-RHS nonlocal V solve, followed by the frozen legal projection."""
        with self.operation('TMATRIX_NONLOCAL_LOCALITY_SOLVE', rhs=self.model.n,
                            LU=1, locality_LU=1, locality_dense_RHS=self.model.n,
                            locality_matrix_products=1, locality_factorization_dimension=self.model.n):
            A = self.eye + T @ self.model.G
            with warnings.catch_warnings():
                warnings.simplefilter('error', la.LinAlgWarning)
                factor = la.lu_factor(A, check_finite=False)
                V = la.lu_solve(factor, T, check_finite=False)
            if not np.all(np.isfinite(V)):
                raise FloatingPointError('Nonfinite nonlocal V')
        with self.operation('TMATRIX_LEGAL_LOCAL_PROJECTION', legal_local_projections=1):
            diagonal = np.diag(V).copy()
            offdiag = V.copy()
            np.fill_diagonal(offdiag, 0)
            m = np.clip(diagonal.real, *self.box)
            row = dict(nonlocal_V_fro=float(la.norm(V)), diagonal_V_fro=float(la.norm(diagonal)),
                       offdiagonal_V_fro=float(la.norm(offdiag)),
                       offdiagonal_relative_to_full=relative_norm(offdiag,V),
                       offdiagonal_relative_to_diagonal=relative_norm(offdiag,diagonal),
                       imaginary_diagonal_fro=float(la.norm(diagonal.imag)),
                       imaginary_diagonal_relative=relative_norm(diagonal.imag,diagonal),
                       legal_local_projection_relative_change=relative_norm(V-np.diag(m),V),
                       nonlocal_V_relative_transpose_defect=relative_norm(V-V.T,V),
                       clipped_lower=int(np.count_nonzero(diagonal.real<self.box[0])),
                       clipped_upper=int(np.count_nonzero(diagonal.real>self.box[1])))
        return m, V, row
