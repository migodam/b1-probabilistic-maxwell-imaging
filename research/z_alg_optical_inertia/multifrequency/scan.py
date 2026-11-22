#!/usr/bin/env python3
"""Compare independent single-frequency and coupled shared-material residual bounds."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'

from collections import Counter
from pathlib import Path
import sys
import json

import numpy as np
from scipy.linalg import block_diag, ldl, solve_triangular

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = ROOT / 'results/multifrequency'
sys.path.insert(0, str(ROOT))
from core.inertia import adjoint, clock, hermitian, result_row, spectral_tolerance, timing

CAP = 0.5
MUS = (1.e-4, 1.e-3, 1.e-2, 1.e-1)
SINGLE_RATIOS = (0., .3, 1., 3.)
JOINT_RATIOS = (.3, 1.)
COUPLINGS = {'block_diagonal': (0., 0.), 'same_sign': (.45, .45),
             'opposite_sign': (.45, -.45), 'optical_off_diagonal_only': (0., .45)}


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def save(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / name
    pending = path.with_suffix('.partial')
    pending.write_text(json.dumps(value, indent=1, ensure_ascii=False) + '\n', encoding='utf-8')
    pending.replace(path)


class Frequency:
    def __init__(self, G, S, B, Y, obs):
        self.G, self.S, self.B, self.Y, self.obs = G, S, B, Y, obs
        self.n, self.s = B.shape
        self.S2 = hermitian(adjoint(S)@S)
        self.SY = [adjoint(S)@y for y in Y]
        self.Y2 = [hermitian(adjoint(y)@y) for y in Y]

    def K(self, mu, ratio):
        mat = self.S2 - mu*hermitian((CAP+1j*ratio)*self.G)
        mat.flat[::self.n+1] += mu
        return mat

    def L(self, mu, ratio):
        return [sy + .5*mu*(CAP+1j*ratio)*self.B for sy in self.SY]


def load():
    p400 = ROOT / 'synthetic/public'
    p500 = HERE / 'data/public'
    m400, m500 = read(p400 / 'manifest.json'), read(p500 / 'manifest.json')
    ids = [row['id'] for row in m500['observations']]
    notes = {row['id']:row for row in m400['observations']}
    obs400 = [notes[k] for k in ids]
    obs500 = m500['observations']
    half = m400['geometries']['PW180_10x20']
    with np.load(p400/half['green_file'], allow_pickle=False) as vg, np.load(p400/half['file'], allow_pickle=False) as vv:
        a = {'G':vg['G'], 'S':vv['S'], 'B':vv['B']}
    with np.load(p500/m500['geometry'], allow_pickle=False) as v:
        b = {k:v[k] for k in ('G','S','B')}
    return ids, [Frequency(**ops, Y=read_observations(p,obs), obs=obs)
                 for ops, p, obs in ((a,p400,obs400),(b,p500,obs500))]


def read_observations(folder, rows):
    values = []
    for row in rows:
        with np.load(folder/row['file'], allow_pickle=False) as z:
            values.append(z['Y'])
    return values


def evaluate(K, Ls, YY, sources):
    started = clock()
    lu, d, perm = ldl(K, check_finite=False)
    pivots = np.flatnonzero(np.diag(d,1))
    plain = np.setdiff1d(np.arange(len(d)), np.r_[pivots,pivots+1])
    diagonal = d.diagonal()[plain].real
    pivots_values = [np.linalg.eigvalsh(d[p:p+2,p:p+2]) for p in pivots]
    values = np.r_[diagonal, np.asarray(pivots_values).ravel()]
    gap = float(np.min(np.abs(values)))
    floor = spectral_tolerance(values)
    if gap <= floor:
        return dict(status='numeric_LDL_unresolved', pivot_min_abs=gap, results=None, **timing(started))
    nu = int(np.count_nonzero(values < 0))
    if nu >= sources:
        results = [result_row(nu,sources) for _ in YY]
        return dict(status='inertia_exhausted', nu=nu, pivot_min_abs=gap, results=results, **timing(started))
    T = lu[perm]
    L = np.concatenate(Ls, axis=1)
    work = solve_triangular(T,L[perm],lower=True,unit_diagonal=True,check_finite=False)
    work[plain] /= diagonal[:,None]
    for p in pivots:
        work[p:p+2] = np.linalg.solve(d[p:p+2,p:p+2],work[p:p+2])
    X = np.empty_like(L)
    X[perm] = solve_triangular(T,work,trans='C',lower=True,unit_diagonal=True,check_finite=False)
    rel = float(np.linalg.norm(K@X-L)/np.linalg.norm(L))
    if rel > 1.e-7:
        return dict(status='numeric_solve_unresolved', nu=nu, pivot_min_abs=gap,
                    relative_solve_defect=rel, results=None, **timing(started))
    values = [result_row(nu,sources, gram-adjoint(l)@x)
              for gram,l,x in zip(YY,Ls,np.split(X,len(Ls),axis=1),strict=True)]
    return dict(status='numeric_only', nu=nu, pivot_min_abs=gap, relative_solve_defect=rel,
                results=values, rhs=L.shape[1], **timing(started))


def joint(f, g, mu, t, r, q):
    z = .5*mu*(CAP*r+1j*t*q)
    K12 = -(z*g.G+z.conjugate()*adjoint(f.G))
    K12.flat[::len(K12)+1] += mu*r
    K = np.block([[f.K(mu,t),K12],[adjoint(K12),g.K(mu,t)]])
    L = [np.block([[l1,z*g.B],[z*f.B,l2]]) for l1,l2 in zip(f.L(mu,t),g.L(mu,t),strict=True)]
    grams = [block_diag(y1,y2) for y1,y2 in zip(f.Y2,g.Y2,strict=True)]
    return evaluate(K,L,grams,f.s+g.s)


def best(rows, index):
    candidates = [(i,row) for i,row in enumerate(rows) if row['results'] is not None]
    i, row = max(candidates, key=lambda x:(x[1]['results'][index]['lower_bound_numeric'], -x[1]['nu']))
    return dict(lower=row['results'][index]['lower_bound_numeric'], nu=row['nu'], row_index=i,
                mu=row['mu'], ratio=row['ratio'])


def summary_by_observation(ids, frequencies, single, coupled):
    results = {}
    for k, obs_id in enumerate(ids):
        terms = [best(rows,k) for rows in single]
        lower_sum = np.hypot(*(term['lower'] for term in terms))
        sample_radius = np.hypot(*(freq.obs[k]['radius_noise_fro'] for freq in frequencies))
        joint_winners = {name:best(values,k) for name,values in coupled.items()}
        for row in joint_winners.values():
            row['excludes_sample_noise_fro_numeric'] = bool(row['lower'] > sample_radius)
            row['gain_over_independent'] = float(row['lower']-lower_sum)
        results[obs_id] = dict(sample_radius_fro=float(sample_radius), independent_terms=terms,
                               independent_sum_norm=float(lower_sum),
                               independent_excludes_sample_noise_fro_numeric=bool(lower_sum>sample_radius),
                               joint=joint_winners)
    return results


def score_proxy(result):
    p400, p500 = ROOT/'synthetic/evaluation', HERE/'data/evaluation'
    proxy_notes = {}
    for obs, v in result['observation_scores'].items():
        terms = {}
        for f,p in ((400,p400),(500,p500)):
            with np.load(p/(obs+'.npz'),allow_pickle=False) as archive:
                terms[f] = {k:float(archive[k]) for k in ('discrepancy_fro','proxy_error_fro')}
        discrepancy = float(np.hypot(terms[400]['discrepancy_fro'],terms[500]['discrepancy_fro']))
        effective_radius = discrepancy+v['sample_radius_fro']
        proxy_notes[obs] = {'terms':terms,'conditional_fine_to_coarse_proxy_fro':discrepancy,
                            'conditional_noise_plus_proxy_radius':effective_radius,
                            'direct_joint_proxy_error':float(np.hypot(terms[400]['proxy_error_fro'],terms[500]['proxy_error_fro'])),
                            'independent_excludes_conditional_proxy_numeric':v['independent_sum_norm']>effective_radius,
                            'joint_excludes_conditional_proxy_numeric':{name:row['lower']>effective_radius
                                                                         for name,row in v['joint'].items()}}
    save('proxy_only.json', proxy_notes)


def main():
    total = clock()
    ids, frequencies = load()
    single, coupled, sections = [], {}, Counter()
    for freq in frequencies:
        rows = []
        for t in SINGLE_RATIOS:
            for mu in MUS:
                rows.append(dict(mu=mu, ratio=t, **evaluate(freq.K(mu,t),freq.L(mu,t),freq.Y2,freq.s)))
        sections['single_1024_LDL'] += len(rows)
        single.append(rows)
    for name,(r,q) in COUPLINGS.items():
        rows = []
        for t in JOINT_RATIOS:
            for mu in MUS:
                rows.append(dict(mu=mu, ratio=t, rho_mu=r, rho_lambda=q, **joint(*frequencies,mu,t,r,q)))
        sections['joint_2048_LDL'] += len(rows)
        coupled[name] = rows
        print(name, 'ν', [row.get('nu') for row in rows], flush=True)
    gauge_check = joint(*frequencies,MUS[1],JOINT_RATIOS[0],-.45,-.45)
    gauge_reference = next(row for row in coupled['same_sign'] if row['mu']==MUS[1] and row['ratio']==JOINT_RATIOS[0])
    maximum_difference = max(abs(a['lower_bound_numeric']-b['lower_bound_numeric'])
        for a,b in zip(gauge_check['results'],gauge_reference['results'],strict=True))
    assert gauge_check['nu']==gauge_reference['nu'] and maximum_difference < 1e-8
    sections['joint_2048_LDL_gauge_control'] += 1
    summarize = summary_by_observation(ids,frequencies,single,coupled)
    result = dict(cap=CAP, frequencies_MHz=[400,500], object_count=6, observation_order=ids,
        parameter_grid=dict(mu=MUS,single_lambda_over_mu=SINGLE_RATIOS,joint_lambda_over_mu=JOINT_RATIOS,
                            joint_off_diagonal_relative=COUPLINGS),
        gauge_check=dict(nu=gauge_check['nu'],max_abs_lower_difference=maximum_difference),
        observations_by_frequency=single, joint_parameter_results=coupled,
        observation_scores=summarize, operations=dict(sections), **timing(total))
    save('numerical.json', result)
    score_proxy(result)


if __name__ == '__main__':
    main()
