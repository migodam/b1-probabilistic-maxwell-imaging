"""Check full-menu overlap under a common gain, and under a known unit gain."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from time import process_time

import numpy as np
from scipy.optimize import brentq, minimize

from data14 import Fit, ROOT, load_cases, ri
from pivot_menu_witnesses import MenuProof, OFFSETS

SELECTED = ('spectral-12', 'fresh6-02', 'fresh6-06', 'fresh6-10', 'fresh6-16', 'confirm6-06')
MARGIN = .997
MODES = ('independent', 'common', 'known_unit')


def gain_ri(z):
    return [float(z.real), float(z.imag)]


def common_gain(case, f):
    d = [case.old-x for x in f]
    coef = np.array([np.vdot(x,x).real + (case.additive_radius / case.gain_radius)**2 for x in f])
    centers = [np.vdot(x,r)/c for x,r,c in zip(f,d,coef)]
    norm_base = np.array([np.vdot(r,r).real for r in d]) - coef*np.abs(centers)**2

    def score(t):
        z = (1-t)*centers[0] + t*centers[1]
        return (norm_base + coef*np.abs(z-centers)**2) / case.additive_radius**2

    fn = lambda t: score(t)[0] - score(t)[1]
    if fn(0.) >= 0.:
        t = 0.
    elif fn(1.) <= 0.:
        t = 1.
    else:
        t = brentq(fn, 0., 1.)
    z = (1-t)*centers[0] + t*centers[1]
    hard = case.hard_gain_radius or case.gain_radius
    if abs(z) > hard:
        z *= .999*hard/abs(z)
    return gain_ri(z)


class ControlledSearch:
    def __init__(self, proof):
        self.proof = proof
        self.case, self.model = proof.case, proof.model
        self.E, self.G = self.case.additive_radius, self.case.gain_radius
        self.r = proof.radius
        self.box = self.case.box
        self.hard = self.case.hard_gain_radius or self.G
        self.material_ids = np.r_[np.arange(6), 8+np.arange(6)]

    def responses(self, point):
        x, u = point[:6], point[6:]
        g = 1+self.G*complex(*u)
        f, j = zip(*(self.model.response(x,float(d)) for d in OFFSETS))
        f, j = np.asarray(f), np.asarray(j)
        v = g*f
        arr = np.concatenate((g*j, self.G*f[...,None], 1j*self.G*f[...,None]), axis=-1)
        J = np.stack((arr.real, arr.imag), axis=-2).reshape(len(OFFSETS), -1, 8)
        return v, J

    def solve(self, materials, gains, mode, maxiter, gap_cap=1.):
        start = process_time()
        n0 = self.model.factorizations
        indices = self.material_ids
        Z = np.zeros((16, 12 if mode=='known_unit' else 14 if mode=='common' else 16))
        if mode=='independent':
            Z[:] = np.eye(16)
        else:
            Z[indices,:12] = np.eye(12)
            if mode=='common':
                Z[[6,7,14,15],12:] = np.r_[np.eye(2),np.eye(2)]
        x0 = np.r_[materials[0], np.asarray(gains[0])/self.G,
                   materials[1], np.asarray(gains[1])/self.G]
        z0 = OFFSETS.index(0.)
        _, aJ = self.responses(x0[:8])
        _, bJ = self.responses(x0[8:])
        W = [np.vstack((J[z0]/self.E, J.reshape(-1,8)/(self.r*np.sqrt(len(OFFSETS))), np.eye(8)[6:]))
             for J in (aJ,bJ)]
        jac = np.block([[W[0],np.zeros_like(W[0])],[np.zeros_like(W[1]),W[1]]])
        _, s, vh = np.linalg.svd(jac @ Z,full_matrices=False)
        T = Z @ (vh.T / np.maximum(s,.002))
        bounds_M = T[indices]
        checked = []
        cache = None
        cap = MARGIN * (self.hard/self.G)**2
        scale = max(.02, np.mean((x0[:6]-x0[8:14])**2))

        def objective(y):
            q = x0 + T@y
            diff = q[:6] - q[8:14]
            g = np.r_[diff,0.,0.,-diff,0.,0.]
            return -np.dot(diff,diff)/(6*scale), -(T.T@g)/(3*scale)

        def con(y):
            nonlocal cache
            if cache is not None and np.array_equal(cache[0], y):
                return cache[1:]
            q = x0 + T@y
            A, AJ = self.responses(q[:8])
            B, BJ = self.responses(q[8:])
            oa, ob = [(x[z0]-self.case.old)/self.E for x in (A,B)]
            da, db = ri(oa), ri(ob)
            h1, h2 = q[6:8], q[14:16]
            delta = (A-B)/(2*self.r)
            d = np.stack((delta.real,delta.imag),axis=-1).reshape(len(OFFSETS),-1)
            Js = np.hstack(((-np.einsum('fs,fsj->fj', d, AJ)/self.r),
                             (np.einsum('fs,fsj->fj', d, BJ)/self.r)))
            old = np.zeros((4,16))
            old[0,:8] = -2*da@AJ[z0]/self.E
            old[1,8:] = -2*db@BJ[z0]/self.E
            old[0,6:8] -= 2*h1
            old[1,14:16] -= 2*h2
            old[2,6:8], old[3,14:16] = -2*h1,-2*h2
            vals = np.r_[MARGIN-da@da-h1@h1,MARGIN-db@db-h2@h2,
                         cap-h1@h1,cap-h2@h2,MARGIN*gap_cap**2-np.sum(d*d,axis=1),
                         q[indices]-self.box[0], self.box[1]-q[indices]]
            deriv = np.vstack((np.vstack((old, Js))@T,bounds_M,-bounds_M))
            if np.min(vals)>=-1e-8:
                checked.append((float(np.sqrt(np.mean((q[:6]-q[8:14])**2))),
                                q[indices].reshape(2,6).copy(), (q[[6,7,14,15]]*self.G).reshape(2,2).tolist()))
            cache = y.copy(), vals, deriv
            return vals, deriv

        res = minimize(objective,np.zeros(T.shape[1]),jac=True,method='SLSQP',
                       constraints=[dict(type='ineq',fun=lambda y:con(y)[0],jac=lambda y:con(y)[1])],
                       options=dict(maxiter=maxiter,ftol=1e-10))
        final_vals, _ = con(res.x)
        verified = None
        for _, x, zz in sorted(checked,key=lambda el:-el[0]):
            r = self.proof.verify(x, zz)
            if r['status']=='all_offsets_overlap' and r['peak_over_2r_upper'] < gap_cap:
                verified = r
                break
        return verified, dict(status=int(res.status), iter=int(res.nit), final_min_constraint=float(np.min(final_vals)),
                              feasible_numeric_iterates=len(checked), cpu_s=process_time()-start,
                              float_forward_derivative_calls=self.model.factorizations-n0)


def test(job):
    case, menu, maxiter = job
    start = process_time()
    proof = MenuProof(case)
    pairs = {}
    for v in menu['trials'].values():
        if 'materials' in v:
            token = tuple(tuple(x) for x in v['materials'])
            pairs.setdefault(token,v)
    outcomes = {k: [] for k in MODES}
    targeted = case.case_id in SELECTED
    for value in pairs.values():
        pts = np.array(value['materials'])
        responses = [proof.model.response(x,0.)[0] for x in pts]
        z = common_gain(case,responses)
        outcomes['common'].append((pts, [z,z], proof.verify(pts,[z,z])))
        zero = [[0.,0.],[0.,0.]]
        outcomes['known_unit'].append((pts,zero,proof.verify(pts,zero)))
        outcomes['independent'].append((pts,value['gains_zeta_ri'],value))
    out = {}
    search = ControlledSearch(proof)
    for mode, trials in outcomes.items():
        valid = [c for _,_,c in trials if c['status']=='all_offsets_overlap']
        prebest = max(valid,key=lambda c:c['material_rms_lower'], default=None)
        attempts=[]
        if targeted:
            latest = sorted(trials,key=lambda q:-float(np.linalg.norm(q[0][0]-q[0][1])))[:1]
            viable = [(p,z,c) for p,z,c in trials if c['status']=='all_offsets_overlap']
            if viable:
                front = max(viable,key=lambda v:v[2]['material_rms_lower'])
                if not np.array_equal(front[0],latest[0][0]):
                    latest.append(front)
            for pts, gains, _ in latest:
                pts=pts.copy()
                if mode=='known_unit':
                    fit=Fit(proof.model,[0.],[case.old])
                    pts=np.array([fit.solve(x,110)[0] if np.linalg.norm(fit.evaluate(x)[0])>case.additive_radius else x
                                  for x in pts])
                witness, result=search.solve(pts,gains,mode,maxiter)
                attempts.append(result)
                if witness:
                    valid.append(witness)
        out[mode]=dict(best=max(valid,key=lambda c:c['material_rms_lower'],default=None),
                       before_optimization=prebest, trials=[t[2] for t in trials], optimizations=attempts)
    return case.case_id,dict(cohort=case.cohort,modes=out,cpu_s=process_time()-start,
                             arb_factors=proof.solves,arb_rhs=3*proof.solves,
                             forward_derivative_calls=proof.model.factorizations)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workers',type=int,default=3)
    p.add_argument('--maxiter',type=int,default=65)
    p.add_argument('--slack-control',type=float)
    a=p.parse_args()
    root=ROOT/'results'
    if a.slack_control is not None:
        forced_slack(a.slack_control,a.maxiter,root)
        return
    archived=json.loads((root/'pivot_menu_witnesses.json').read_text())
    jobs=[(case,archived['cases'][case.case_id],a.maxiter) for case in load_cases()]
    rows={}
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for case,data in pool.map(test,jobs):
            rows[case]=data
            print(case,{m:(round(v['best']['material_rms_lower'],7),round(v['best']['peak_over_2r_upper'],7))
                        if v['best'] else None for m,v in data['modes'].items()},flush=True)
    (root/'pivot_menu_gain_controls.json').write_text(json.dumps(dict(offsets=OFFSETS,targeted_optimization_cases=SELECTED,
                    cases=rows),indent=2,allow_nan=False)+'\n')


def forced_slack(target,maxiter,root):
    chosen = ('fresh6-02','fresh6-16','confirm6-06')
    controls = json.loads((root/'pivot_menu_gain_controls.json').read_text())['cases']
    out={}
    for case in load_cases():
        if case.case_id not in chosen:
            continue
        proof=MenuProof(case)
        search=ControlledSearch(proof)
        rows={}
        for mode in ('common','known_unit'):
            anchor=controls[case.case_id]['modes'][mode]['best']
            if anchor is None or anchor['material_rms_lower']<.01:
                continue
            m,z=anchor['materials'],anchor['gains_zeta_ri']
            if anchor['peak_over_2r_upper']<target:
                witness, attempts = anchor, None
            else:
                witness, attempts=search.solve(m,z,mode,maxiter,gap_cap=target)
            rows[mode]=dict(best=witness,optimization=attempts)
        out[case.case_id]=dict(modes=rows,arb_factors=proof.solves,arb_rhs=3*proof.solves)
        print(case.case_id,{m:v['best']['material_rms_lower'] if v['best'] else None for m,v in rows.items()},flush=True)
    (root/'pivot_menu_gain_slack.json').write_text(json.dumps(dict(offsets=OFFSETS,target=target,
                    cases=out),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
