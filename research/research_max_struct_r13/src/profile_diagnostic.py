"""Public-only profiles of ONE fixed scalar-CIE current block.

This is a diagnostic of nested affine current families, not an inversion
algorithm or common-objective interpretation of dual-beta G19 alternation.
Stable thin QR and a small receiver SVD solve every source's weighted least
squares without forming a normal matrix or changing source normalization.
"""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_key] = '1'
from collections import Counter
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np
from scipy.linalg import qr, svd, solve_triangular, lstsq
from research.research_max_struct_r6.src.fft_tsom_ops import NativeFourierAPCS
from research.research_max_struct_r10.src.twod_adapter import noise_ball_rank


class ScalarCurrentProfiles:
    """Profiles of .5||SJ-Y||_Od²+.5 rho||p(J-m(E+GJ))||_Os².

    beta, source weights, DPCS and current bases stay fixed across material
    probes. The common material prior is reported separately and cancels from
    representation penalties. These profiles do not use G19's material metric.
    """
    def __init__(self, model, Y, Vd, JD, *, mf=3, beta=1., rho=1.):
        self.model, self.Y = model, np.asarray(Y, complex)
        self.Vd, self.JD = np.asarray(Vd, complex), np.asarray(JD, complex)
        self.n, self.ns = model.E.shape
        self.beta, self.rho = float(beta), float(rho)
        if self.beta <= 0 or self.rho <= 0:
            raise ValueError('Positive fixed beta/rho required')
        self.cost = Counter()
        self.wd = 1/np.maximum(np.sum(abs(Y)**2, axis=0), 1e-300)
        self.ws = 1/np.maximum(np.sum(abs(model.E)**2, axis=0), 1e-300)
        native = NativeFourierAPCS(model.shape, mf, Vd, self.cost)
        F = native.apply(np.eye(native.k, dtype=complex))
        u, singular, _ = svd(F, full_matrices=False, check_finite=False)
        self.cost['restricted_basis_setup_SVD'] += 1
        keep = singular > 1e-12*max(singular[0], 1e-300)
        complete, _ = qr(Vd, mode='full', check_finite=False)
        self.cost['full_APCS_complement_QR'] += 1
        self.bases = dict(restricted=u[:, keep], full_apcs=complete[:, Vd.shape[1]:],
                          full_current=np.eye(self.n, dtype=complex))
        self.info = dict(beta_current=self.beta, rho=self.rho, restricted_mf=int(mf),
            DPCS_rank=Vd.shape[1], independent_complex_spatial_ranks={k:U.shape[1] for k,U in self.bases.items()},
            actual_real_current_DOF={k:2*U.shape[1]*self.ns for k,U in self.bases.items()},
            common_offset='original receiver DPCS JD, also retained in full-current coordinates',
            source_data_weights=self.wd.tolist(), source_state_weights=self.ws.tolist(),
            material_prior='separate, common across families; excluded from penalty differences')

    def evaluate(self, m, retain_current=False):
        m = np.asarray(m, float)
        p = self.beta/(1+self.beta*m)
        C = p[:, None]*(np.eye(self.n)-m[:, None]*self.model.G)
        target = p[:, None]*m[:, None]*self.model.E
        s0 = self.model.S@self.JD-self.Y
        c0 = C@self.JD-target
        self.cost['profile_material_evaluations'] += 1
        self.cost['offset_dense_C_RHS'] += self.ns
        self.cost['offset_receiver_RHS'] += self.ns
        rows = {}
        for family, U in self.bases.items():
            started = time.process_time()
            CU, SU = C@U, self.model.S@U
            self.cost['profile_dense_C_basis_RHS'] += U.shape[1]
            self.cost['profile_receiver_basis_RHS'] += U.shape[1]
            Q, R, piv = qr(CU, mode='economic', pivoting=True, check_finite=False)
            self.cost['profile_QR'] += 1
            diagonal = abs(np.diag(R))
            qrank = int(np.count_nonzero(diagonal > 1e-12*max(diagonal.max(initial=0.), 1e-300)))
            if qrank != U.shape[1]:
                # Rare near-pole/rank failure: retain stable direct QR LS with
                # exact source weights rather than normal equations or ridge.
                coeff = np.empty((U.shape[1], self.ns), complex)
                for s in range(self.ns):
                    A = np.vstack([np.sqrt(self.wd[s])*SU,
                        np.sqrt(self.rho*self.ws[s])*CU])
                    b = -np.r_[np.sqrt(self.wd[s])*s0[:, s],
                        np.sqrt(self.rho*self.ws[s])*c0[:, s]]
                    coeff[:, s], _, _, _ = lstsq(A, b, cond=1e-12,
                        lapack_driver='gelsy', check_finite=False)
                    self.cost['fallback_source_direct_QR_LS'] += 1
                strategy = 'rank-deficient fallback: original stacked weighted QR least squares per source'
            else:
                # CU[:,piv]=QR, x=R coeff[piv]. Data acts through B=SU[:,piv]/R.
                B = solve_triangular(R.T, SU[:, piv].T, lower=True,
                                     check_finite=False).T
                x0 = -Q.conj().T@c0
                L, sig, Vh = svd(B, full_matrices=False, check_finite=False)
                self.cost['profile_small_receiver_SVD'] += 1
                ratio = self.wd/(self.rho*self.ws)
                factor = sig[:, None]*ratio[None, :]/(1+sig[:, None]**2*ratio[None, :])
                x = x0-Vh.conj().T@(factor*(L.conj().T@(s0+B@x0)))
                coeff = np.empty_like(x)
                coeff[piv] = solve_triangular(R, x, check_finite=False)
                self.cost['profile_triangular_RHS'] += B.shape[0]+self.ns
                strategy = 'pivoted thin QR of C U plus receiver SVD; exact sourcewise weights'
            J = self.JD+U@coeff
            Qtotal = self.model.E+self.model.G@J
            D = self.model.S@J-self.Y
            residual = p[:, None]*(J-m[:, None]*Qtotal)
            data = .5*float(np.sum(abs(D)**2*self.wd[None, :]))
            state = .5*self.rho*float(np.sum(abs(residual)**2*self.ws[None, :]))
            # Envelope derivative: current is optimized in a fixed affine
            # family. d[p(J-mQ)]/dm=-beta*(Q+beta J)/(1+beta m)^2.
            derivative = -self.beta*(Qtotal+self.beta*J)/(1+self.beta*m[:, None])**2
            gradient = self.rho*np.real(np.sum(residual.conj()*derivative*self.ws[None, :], axis=1))
            back = self.model.S.conj().T@(D*self.wd[None, :])
            back += self.rho*C.conj().T@(residual*self.ws[None, :])
            tangent_stationarity = U.conj().T@back
            self.cost['profile_final_dense_G_RHS'] += self.ns
            self.cost['profile_final_receiver_RHS'] += self.ns
            self.cost['profile_stationarity_dense_CH_RHS'] += self.ns
            self.cost['profile_stationarity_receiver_adjoint_RHS'] += self.ns
            self.cost['profile_exact_LS_solutions'] += 1
            row = dict(value=data+state, data_term=data, state_term=state, gradient=gradient,
                independent_complex_spatial_rank=U.shape[1], QR_rank=qrank,
                QR_min_abs_diagonal=float(diagonal.min(initial=np.inf)),
                QR_max_abs_diagonal=float(diagonal.max(initial=0.)),
                tangent_stationarity_norm=float(np.linalg.norm(tangent_stationarity)),
                strategy=strategy, CPU_seconds=time.process_time()-started)
            if retain_current:
                row.update(J=J, Q=Qtotal)
            rows[family] = row
        rows['penalties'] = dict(Fourier=rows['restricted']['value']-rows['full_apcs']['value'],
            fixed_DPCS=rows['full_apcs']['value']-rows['full_current']['value'],
            total=rows['restricted']['value']-rows['full_current']['value'])
        return rows


def material_directions(model, Y, m, profiles):
    """Five public directions; no material truth, target IoU or selection."""
    before = dict(model.cost)
    value, gradient = model.data_value_gradient(m, Y)
    state = model.state(m)
    residual = state['Y']-Y
    port = model.S.conj().T@residual
    port_gradient = np.real(np.sum(state['Q'].conj()*port, axis=1))/float(np.vdot(Y, Y).real)
    raw = dict(physical_Phi_descent=-gradient,
        full_profile_descent=-profiles['full_current']['gradient'],
        Fourier_penalty_descent=-(profiles['restricted']['gradient']-profiles['full_apcs']['gradient']),
        DPCS_penalty_descent=-(profiles['full_apcs']['gradient']-profiles['full_current']['gradient']),
        linear_port_backprojection=-port_gradient)
    directions, norms = {}, {}
    for name, d in raw.items():
        norm = float(np.linalg.norm(d, np.inf))
        directions[name] = d/max(norm, 1e-300)
        norms[name] = dict(raw_inf_norm=norm, zero_direction=norm == 0,
            normalization='unit material infinity norm', uses_truth=False)
    return directions, dict(physical_Phi=float(value), direction_info=norms,
        physical_model_delta={k:model.cost.get(k, 0)-before.get(k, 0) for k in set(model.cost)|set(before)},
        extra_linear_port_receiver_adjoint_complex_RHS=model.ns,
        linear_port_qualification='unpropagated S* residual material proxy; not Maxwell adjoint or GN direction')


def _compact(rows):
    return {family:{k:v for k,v in row.items() if k not in ('gradient', 'J', 'Q')}
            for family,row in rows.items()}


def run_public_input(source, out, *, step=.02, mf=3):
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel
    wall, cpu = time.perf_counter(), time.process_time()
    out.mkdir(parents=True, exist_ok=False)
    with np.load(source, allow_pickle=False) as z:
        forbidden = [k for k in z.files if any(word in k.lower() for word in ('truth', 'm_true', 'ground', 'label'))]
        if forbidden:
            raise ValueError('Public-only profile input contains forbidden truth-like fields')
        public = {k:z[k].copy() for k in z.files}
    m = public.pop('m_init')
    Y = public['Y']
    model = MaxwellModel(public)
    decision = noise_ball_rank(model.S, Y, float(public['sigma']), out/'RECEIVER_DPCS.npz')
    if decision['rank'] is None:
        raise ValueError('No supported public DPCS prefix')
    with np.load(out/'RECEIVER_DPCS.npz', allow_pickle=False) as z:
        operator = ScalarCurrentProfiles(model, Y, z['Vd'], z['JD'], mf=mf)
    base = operator.evaluate(m, retain_current=True)
    directions, direction_info = material_directions(model, Y, m, base)
    arrays = dict(m=m, direction_names=np.array(list(directions)), directions=np.array(list(directions.values())))
    for family in operator.bases:
        arrays['J_'+family] = base[family]['J']
        arrays['gradient_'+family] = base[family]['gradient']
    lower, upper = .02, 20.
    probes = []
    for name, d in directions.items():
        for sign in (-1, 1):
            candidate = np.clip(m+sign*step*d, lower, upper)
            delta = candidate-m
            rows = operator.evaluate(candidate)
            derivative_prediction = {family:float(np.dot(base[family]['gradient'], delta)) for family in operator.bases}
            probes.append(dict(direction=name, sign=sign, step=step, material_inf_change=float(np.linalg.norm(delta, np.inf)),
                material_l2_change=float(np.linalg.norm(delta)), clipped_cells=int(np.count_nonzero(candidate != m+sign*step*d)),
                profiles=_compact(rows), predicted_first_order_profile_change=derivative_prediction,
                actual_profile_change={family:rows[family]['value']-base[family]['value'] for family in operator.bases},
                penalty_change={k:rows['penalties'][k]-base['penalties'][k] for k in rows['penalties']}))
            arrays[name+('_minus' if sign < 0 else '_plus')] = candidate
    # One small one-sided envelope derivative wiring check, with the same box
    # retraction. It validates new diagnostic code, not a frozen physics identity.
    eps = 1e-6
    tiny = np.clip(m+eps*directions['full_profile_descent'], lower, upper)-m
    tiny_rows = operator.evaluate(m+tiny)
    derivative_check = {family:dict(actual=tiny_rows[family]['value']-base[family]['value'],
        predicted=float(np.dot(base[family]['gradient'], tiny)),
        remainder=float(tiny_rows[family]['value']-base[family]['value']-np.dot(base[family]['gradient'], tiny)))
        for family in operator.bases}
    report = dict(status='PUBLIC_FIXED_CURRENT_BLOCK_PROFILE_DIAGNOSTIC', source=str(source), truth_access='NONE',
        qualification='Exact diagnostic of ONE scalar-CIE current block at beta=1; not common-objective dual-beta/G19 theory, not a new algorithm',
        material_source='public exposed initializer m_init; no recovered/truth material selection',
        current_block=operator.info, public_DPCS_rule=decision, base=_compact(base),
        common_material_prior=float(model.reg.evaluate(m)[0]), direction_construction=direction_info,
        probes=probes, envelope_derivative_wiring=derivative_check,
        profile_cost=dict(operator.cost), physical_model_cost=dict(model.cost),
        CPU_seconds=time.process_time()-cpu, wall_seconds=time.perf_counter()-wall,
        execution='local existing NumPy/SciPy environment, BLAS1; input/operator/model setup included; final output serialization excluded')
    np.savez_compressed(out/'PROFILE_MATERIALS_AND_BASE_CURRENTS.npz', **arrays)
    (out/'PROFILE.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', type=Path, action='append', required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--step', type=float, default=.02)
    p.add_argument('--mf', type=int, default=3)
    a = p.parse_args()
    if len(a.input) > 3:
        raise ValueError('Bounded public diagnostic accepts at most three inputs')
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        for source in a.input:
            report = run_public_input(source.resolve(), a.output/source.stem, step=a.step, mf=a.mf)
            print(json.dumps(dict(case=source.stem, CPU_seconds=report['CPU_seconds'],
                base_penalties=report['base']['penalties'], physical_model_cost=report['physical_model_cost'])), flush=True)


if __name__ == '__main__':
    main()
