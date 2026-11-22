#!/usr/bin/env python3
"""Independently certify selected MAX-STRUCT profiles on the archived finite model."""
import argparse
from collections import Counter, defaultdict
from pathlib import Path

from scan import HERE, R10_REFERENCE, adjoint, clock, hermitian, load_public, read, timing, write

import numpy as np
from scipy.linalg import cholesky, eigh
from flint import acb, acb_mat, arb, ctx
from threadpoolctl import threadpool_limits

from core.certificate import ball, gershgorin, safe_float, star

ARMS = ('pure_joint', 'optical_joint', 'definite_joint', 'indefinite_joint', 'single_original', 'single_coherent')


def hi(x):
    return safe_float(x, upper=True)


def lo(x):
    return safe_float(x, upper=False)


def frob_upper(matrix):
    return sum((abs(entry).upper()**2 for entry in matrix.entries()), arb(0)).upper().sqrt().upper()


def eye(n):
    result = acb_mat(n, n)
    for i in range(n):
        result[i,i] = 1
    return result


def positive_norm_bounds(scores):
    lower = sum((max(arb(0), x.lower()) for x in scores), arb(0)).lower()
    upper = sum((max(arb(0), x.upper()) for x in scores), arb(0)).upper()
    return [max(0.,lo(max(arb(0), lower).sqrt())), max(0.,hi(max(arb(0),upper).sqrt()))]


def numeric_score(row, obs_i, arm):
    nu = row['nu']
    if nu is None or nu >= 12:
        return 0.
    if ((arm == 'pure_joint' and row['lambda_over_mu'] != 0)
        or (arm == 'optical_joint' and row['lambda_over_mu'] == 0)
        or (arm in ('definite_joint', 'single_original', 'single_coherent') and nu != 0)
        or (arm == 'indefinite_joint' and nu == 0)):
        return 0.
    terms = row['case_matrices'][obs_i]
    if arm == 'single_original':
        return float(np.sqrt(max(0.,*terms['diagonal'])))
    if arm == 'single_coherent':
        return float(np.sqrt(max(0.,terms['eig'][0])))
    return float(np.sqrt(np.maximum(terms['eig'][nu:],0.).sum()))


def select(grid, caps):
    index = {row['profile']:row for row in grid['values']}
    jobs = defaultdict(set)
    reasons = defaultdict(set)
    ranks = {}
    for cap in caps:
        cap_groups = [r for r in grid['values'] if r['cap'] == cap]
        for row in cap_groups:
            if (cap, row['mu'], row['lambda_over_mu']) == R10_REFERENCE:
                jobs[row['profile']].update(range(len(grid['cases'])))
                reasons[row['profile']].add('frozen_R10_weight')
        best_min = min((r for r in cap_groups if r['nu'] is not None),
                       key=lambda r: (r['nu'],r['profile']))
        if best_min['nu'] >= 12:
            jobs[best_min['profile']]
            reasons[best_min['profile']].add('numeric_minimum_inertia_when_all_profiles_exhaust_sources')
        for obs_i, case in enumerate(grid['cases']):
            ranks[case, cap] = {}
            for arm in ARMS:
                ranked = [(numeric_score(r, obs_i, arm),r['profile']) for r in cap_groups]
                committed = sorted((v for v in ranked if v[0] > 0.),
                                   key=lambda v: (-v[0],v[1]))[:2]
                ranks[case,cap][arm] = [dict(profile=pid, numeric=score) for score,pid in committed]
                for score,pid in committed:
                    jobs[pid].add(obs_i)
                    reasons[pid].add(arm)
    return index, jobs, reasons, ranks


class SharedGeometry:
    def __init__(self, G, S, B):
        self.S = S
        self.B = B
        self.S2 = hermitian(adjoint(S) @ S)
        self.H = hermitian(G)
        self.O = (G-adjoint(G))/(2j)
        self.n = len(G)
        self.I = np.eye(self.n)
        self.ia = eye(self.n)
        Ga, self.sa, self.ba = map(ball,(G,S,B))
        self.s2a = star(self.sa)*self.sa
        self.ha = (Ga+star(Ga))/2
        self.oa = (Ga-star(Ga))/acb(0,2)
        self.sta = star(self.sa)

    def certify(self, profile):
        u, mu, t = [profile[k] for k in ('cap','mu','lambda_over_mu')]
        started = clock()
        K = self.S2+mu*(self.I-u*self.H+t*self.O)
        d, basis = eigh(K,check_finite=False,driver='evd')
        shift = float(np.min(abs(d)) / 2.)
        if not shift > 0:
            raise ArithmeticError('No nonzero gap proposal.')
        nu = int(np.count_nonzero(d<0))
        negative = basis[:,:nu]
        W = negative * np.sqrt(-2.*d[:nu])
        F = cholesky(K+W@adjoint(W)-shift*self.I,lower=True,check_finite=False)
        double_cost = timing(started)
        check = clock()
        ae_u, ae_mu, ae_t = [acb(x) for x in (u,mu,t)]
        Ka = self.s2a + ae_mu*(self.ia-ae_u*self.ha+ae_t*self.oa)
        Fa = ball(F)
        defect = Ka - acb(shift)*self.ia - Fa*star(Fa)
        if nu:
            Wa = ball(W)
            defect += Wa*star(Wa)
        upper_defect = frob_upper(defect)
        gamma = (arb(shift)-upper_defect).lower()
        if not gamma > 0:
            raise ArithmeticError('Positive shifted Cholesky defect cannot establish the gap.')
        negative_ceiling = None
        if nu:
            Na = ball(negative)
            neg_edge = gershgorin(star(Na)*(Ka*Na+acb(shift)*Na),upper=True)
            negative_ceiling = hi(neg_edge)
            if not neg_edge < 0:
                raise ArithmeticError('Negative shifted subspace is not verified.')
        z = (ae_u+acb(0,1)*ae_t)*ae_mu/2
        zb_numeric = (.5*mu*(u+1j*t))*self.B
        report = dict(proof_nu=nu, proof_gap_lower=lo(gamma),
                      negative_projected_max_gershgorin_upper=negative_ceiling,
                      positive_shifted_cholesky_defect_upper=hi(upper_defect),
                      geometry_double=double_cost, geometry_Arb=timing(check), **timing(started))
        return (d,basis,Ka,gamma, z,zb_numeric),report

    def sources(self, data, ys, y_arb, y2_arb, sy_arb, radii):
        d,basis,Ka,gamma,z,zb_numeric = data
        nu = int(np.count_nonzero(d<0))
        if nu >= 12:
            return dict(joint_interval=[0.,0.], single_original_interval=[0.,0.],
                        single_coherent_interval=[0.,0.], exhausted=True, joint_pass=False,
                        single_original_pass=False, single_coherent_pass=False)
        started = clock()
        radii = float(radii)
        L = self.S.conj().T @ ys + zb_numeric
        X = basis @ ((basis.conj().T @ L)/d[:,None])
        Xa = ball(X)
        La = sy_arb + z*self.ba
        KX = Ka*Xa
        D = KX-La
        correction = (frob_upper(D)**2/gamma).upper()
        H0 = y2_arb-star(La)*Xa-star(Xa)*La+star(Xa)*KX
        center = np.array([[complex(item) for item in row] for row in H0.tolist()])
        vals,V = eigh(hermitian(center),check_finite=False,driver='evd')
        V = ball(V)
        Vstar = star(V)
        gram_error = frob_upper(Vstar*V-eye(len(vals)))
        diag = acb_mat(len(vals),len(vals))
        for j,x in enumerate(vals):
            diag[j,j] = acb(float(x))
        spectral_error = frob_upper(Vstar*H0*V-diag)
        if not gram_error < 1:
            raise ArithmeticError('Source eigenvectors not numerically independent.')
        envelope = ((spectral_error+gram_error*arb(float(np.max(abs(vals)))))/(1-gram_error)
                    +correction).upper()
        eigenbounds = [arb(float(x),envelope) for x in vals[::-1]]
        joint = positive_norm_bounds(eigenbounds[nu:])
        if nu:
            original = coherent = [0.,0.]
            per_source = None
        else:
            per_source = [positive_norm_bounds([item.real+arb(0.,correction)])
                          for item in (H0[j,j] for j in range(len(vals)))]
            original = [max(a[i] for a in per_source) for i in (0,1)]
            coherent = positive_norm_bounds(eigenbounds[:1])
        return dict(joint_interval=joint, single_original_interval=original,
                    single_coherent_interval=coherent, original_source_intervals=per_source,
                    H_position_error_upper=hi(envelope), inverse_correction_upper=hi(correction),
                    joint_pass=joint[0]>radii, single_original_pass=original[0]>radii,
                    single_coherent_pass=coherent[0]>radii, **timing(started))


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--public',type=Path,default=HERE/'public_inputs.npz')
    args = cli.parse_args()
    started = clock()
    frozen = read(HERE/'frozen.json')
    grid = read(HERE/'numerical_grid.json')
    index,jobs,reasons,ranks = select(grid,frozen['caps'])
    G,S,B,Y,radii = load_public(args.public, frozen)
    ctx.prec,ctx.threads = 128,1
    geo = SharedGeometry(G,S,B)
    y_arb = list(map(ball,Y))
    y2_arb = [star(y)*y for y in y_arb]
    sy_arb = [geo.sta*y for y in y_arb]
    checked = {}
    cost = Counter()
    setup = timing(started)
    for pid in sorted(jobs, key=lambda i: (not index[i]['historical_r10_profile'],index[i]['cap'],i)):
        p = index[pid]
        row = {k:v for k,v in p.items() if k!='case_matrices'}
        row['selected_for']=sorted(reasons[pid])
        cost['geometry_profiles_attempted'] += 1
        try:
            data, proof = geo.certify(p)
            row.update(status='geometry_verified', **proof, per_case={})
            cost['geometry_profiles_verified'] += 1
            for obs_i in sorted(jobs[pid]):
                case = frozen['cases'][obs_i]
                cost['source_proofs_attempted'] += 1
                try:
                    row['per_case'][case] = geo.sources(data,Y[obs_i],y_arb[obs_i],
                                                        y2_arb[obs_i],sy_arb[obs_i],radii[obs_i])
                    cost['source_proofs_verified'] += 1
                    if proof['proof_nu'] < 12:
                        cost['certification_source_rhs'] += 12
                except (ArithmeticError, ValueError, np.linalg.LinAlgError) as exc:
                    row['per_case'][case] = dict(status='source_interval_unresolved',reason=str(exc))
        except (ArithmeticError,ValueError,np.linalg.LinAlgError) as exc:
            row.update(status='geometry_unresolved',reason=str(exc))
        checked[pid]=row
        if len(checked)%10==0 or len(checked)==len(jobs) or p['historical_r10_profile']:
            report = dict(receipts=checked, jobs_planned=len(jobs), setup=setup,
                          precision=ctx.prec, Arb_threads=ctx.threads, operations=dict(cost),
                          numerical_proposals=[dict(case=k[0],cap=k[1],ranks=v) for k,v in ranks.items()],
                          **timing(started))
            write(HERE/'verified.json',report)
            print('verified',len(checked),'/',len(jobs),'last cap',p['cap'], 'time',
                  round(report['wall_seconds'],2),'s',flush=True)


if __name__=='__main__':
    with threadpool_limits(1):
        main()
