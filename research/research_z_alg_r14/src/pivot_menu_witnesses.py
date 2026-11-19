"""Find fixed material pairs whose mean gap is small across the entire frequency menu."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from time import perf_counter, process_time

import numpy as np
from flint import acb_mat, arb

from data14 import Action, FrequencyTM, MENU, ROOT, load_cases
from pivot_bounds import ball, lower_float, mat, row, sqnorm, upper_float

OFFSETS = tuple(sorted({0., *(delta for action in MENU for delta in action.deltas)}))
SINGLE = Action('unit_illumination_block', (0.,), (1.,))


class MenuProof:
    def __init__(self, case):
        self.case = case
        self.model = FrequencyTM(case.geometry)
        self.radius = case.new_noise_radius(SINGLE)
        self.eps, self.gamma, self.r = (arb(float(x)) for x in (case.additive_radius, case.gain_radius, self.radius))
        self.hard = arb(float(case.hard_gain_radius or case.gain_radius))
        self.archive = {}
        self.solves = 0
        self.operators = {d: self.freeze(d) for d in OFFSETS}

    def freeze(self, offset):
        H, B, R, *_ = self.model.kernels(offset)
        return mat(H), mat(B), row(np.ravel(R))

    def forward(self, point, offset):
        key = (tuple(point), offset)
        if key not in self.archive:
            H, B, R = self.operators[offset]
            T = acb_mat([[(1 / ball(m) if i == j else 0) - H[i,j]
                          for j in range(len(point))] for i, m in enumerate(point)])
            self.archive[key] = (R * T.solve(B)).entries()
            self.solves += 1
        return self.archive[key]

    def verify(self, materials, gains):
        if np.min(materials) < self.case.box[0] or np.max(materials) > self.case.box[1]:
            return dict(status='material_outside_box')
        zs = [ball(complex(*z)) for z in gains]
        if any(abs(z).upper() > self.hard for z in zs):
            return dict(status='hard_gain_outside')

        def signal(p, z, freq):
            return [(1 + z) * f for f in self.forward(p, freq)]

        old = [[ball(y) - pred for y, pred in zip(self.case.old, signal(p, z, 0.))]
               for p, z in zip(materials, zs)]
        old_costs = [sqnorm(u) / self.eps**2 + sqnorm([z]) / self.gamma**2 for u, z in zip(old, zs)]
        if not all(s < 1 for s in old_costs):
            return dict(status='old_noise_outside', old_score_upper=[upper_float(q) for q in old_costs])
        bounds, all_pass = [], True
        for delta in OFFSETS:
            a, b = [signal(p, z, delta) for p, z in zip(materials, zs)]
            gap2 = sqnorm([x-y for x,y in zip(a,b)])
            all_pass &= gap2 < (2*self.r)**2
            bounds.append(upper_float(gap2.sqrt()))
        dist = (sum(((arb(float(p)) - arb(float(q)))**2 for p, q in zip(*materials)), arb(0))/len(materials[0])).sqrt()
        peak = max(bounds)
        worst = OFFSETS[bounds.index(peak)]
        ratio_upper = upper_float(arb(peak) / (2*self.r))
        return dict(status='all_offsets_overlap' if all_pass else 'full_menu_overlap_unproved',
                    materials=np.asarray(materials).tolist(), gains_zeta_ri=gains,
                    old_score_upper=[upper_float(q) for q in old_costs], material_rms_lower=lower_float(dist),
                    gap_upper_by_offset=bounds, peak_gap_upper=peak, peak_offset=worst,
                    peak_over_2r_upper=ratio_upper,
                    slack_fraction_lower=lower_float(1-arb(ratio_upper)))


def screen(job):
    case, record = job
    start, wall = process_time(), perf_counter()
    proofer = MenuProof(case)
    trials = {}
    for action, run in record['actions'].items():
        a = run['witness']
        if a is not None:
            trials[action] = proofer.verify(a['materials'], a['gains_zeta_ri'])
    eligible = [(name, candidate) for name, candidate in trials.items() if candidate['status'] == 'all_offsets_overlap']
    best = max(eligible, key=lambda x:x[1]['material_rms_lower'], default=None)
    record = dict(case_id=case.case_id, cohort=case.cohort, single_block_noise_radius=proofer.radius,
                  best_source=best[0] if best else None, best=best[1] if best else None, trials=trials,
                  arb_forward_factors=proofer.solves, arb_forward_rhs=3*proofer.solves,
                  cpu_s=process_time()-start, wall_s=perf_counter()-wall)
    return case.case_id, record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workers', type=int, default=3)
    p.add_argument('--case', action='append')
    a = p.parse_args()
    saved = json.loads((ROOT / 'results/pivot_diameter_witnesses.json').read_text())['cases']
    jobs = [(c,saved[c.case_id]) for c in load_cases() if c.case_id in saved and (not a.case or c.case_id in a.case)]
    rows = {}
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for key, val in pool.map(screen, jobs):
            rows[key] = val
            best = val['best']
            if best:
                print(key, val['best_source'], best['material_rms_lower'], best['peak_over_2r_upper'], flush=True)
    out = dict(offsets=OFFSETS, target='per-offset three complex channels, fixed known source/receiver geometry, any nonnegative frequency energy allocation with sum at most one',
               status='strict for the recorded discrete kernels and single-block noise radius', cases=rows)
    (ROOT/'results/pivot_menu_witnesses.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
