"""Old-data material outer regions under an equal CPU allowance per method."""
import argparse
from collections import Counter, deque
from concurrent.futures import ProcessPoolExecutor
import json
from math import fsum, prod
from time import perf_counter, process_time

import numpy as np

from data14 import FrequencyTM, ROOT, load_cases
from pivot_bounds import OldCase, midpoint
from pivot_controls import Cofactors, direct_bernstein, quadratic_positive, split

METHODS = ('r13', 'r10_scalar', 'r10_shared', 'direct')


def make_case(public):
    H, B, R, *_ = FrequencyTM(public.geometry).kernels(0.)
    lo, hi = [np.full(H.shape[0], x, dtype=float) for x in public.box]
    return OldCase(H, R, B, np.zeros(B.shape[1], dtype=complex), public.old,
                   public.additive_radius, public.gain_radius, lo, hi)


def metrics(queue, finished, lo0, hi0):
    rows = [(lo, hi) for lo, hi, *_ in queue] + finished
    if not rows:
        return dict(retained_boxes=0, unresolved_relative_volume=0., outer_hull=None,
                    normalized_hull_diagonal=0., mean_relative_hull_width=0.,
                    volume_weighted_local_diagonal=0.), np.empty((0, 2, len(lo0)))
    boxes = np.asarray(rows)
    low, high = boxes[:, 0], boxes[:, 1]
    initial = hi0 - lo0
    widths = (high - low) / initial
    volume = np.prod(widths, axis=1)
    hull = np.asarray([np.min(low, axis=0), np.max(high, axis=0)])
    hullwidth = (hull[1] - hull[0]) / initial
    norms = np.linalg.norm(widths, axis=1) / np.sqrt(len(lo0))
    return dict(retained_boxes=len(boxes), unresolved_relative_volume=float(volume.sum()),
                outer_hull=hull.tolist(), normalized_hull_diagonal=float(np.linalg.norm(hullwidth) / np.sqrt(len(lo0))),
                mean_relative_hull_width=float(np.mean(hullwidth)),
                volume_weighted_local_diagonal=float(np.average(norms, weights=volume)) if np.any(volume) else 0.), boxes


def run_arm(public, method, budget, tail_depth, trace=False):
    cpu0, wall0 = process_time(), perf_counter()
    case = make_case(public)
    lo, hi = case.lower, case.upper
    cofactor = None
    controls = None
    if method in ('r10_scalar', 'r10_shared'):
        cofactor = Cofactors(case, tail_depth)
        controls = cofactor.controls
    elif method == 'direct':
        controls = direct_bernstein(case)
    setup = process_time() - cpu0
    axes = tuple(range(case.k)) if method == 'direct' else case.pivots
    pending = deque([(lo, hi, 0, controls, 'untested')])
    finished = []
    reasons = Counter()
    evaluated_dark_axes = 0
    max_depth = 0
    ntest = 0
    deferred_status = Counter()
    targets = iter((.5, .1, .01, .001))
    target = next(targets)
    attained = {}
    removed = []
    total_removed = 0.
    original_volume = prod(hi - lo)
    while pending and process_time() - cpu0 < budget:
        a, b, depth, C, previous = pending.popleft()
        ntest += 1
        max_depth = max(depth, max_depth)
        if method == 'r13':
            bound, status, dark, q = case.cover(a, b)
            evaluated_dark_axes += dark
        elif cofactor is not None:
            bound = cofactor.bound(a, b, C)
            status = 'cofactor_exclusion' if bound is None else 'cofactor_open'
            if bound is not None and method == 'r10_shared' and cofactor.joint_bound(*bound):
                bound = None
                status = 'shared_error_exclusion'
        else:
            reject = np.all(C > 0) or np.all(quadratic_positive(C))
            bound = None if reject else (a, b)
            status = 'direct_exclusion' if reject else 'direct_open'
        reasons[status] += 1
        if trace and target:
            clipped = 0. if bound is None else prod(bound[1] - bound[0])
            parent_volume = prod(b - a)
            if clipped < parent_volume:
                removed.append((parent_volume - clipped) / original_volume)
            if ntest % 16 == 0 and removed:
                total_removed = fsum((total_removed, *removed))
                removed.clear()
                unresolved = 1. - total_removed
                while target and unresolved + 1e-12 <= target:
                    attained[str(target)] = dict(cpu_s=process_time()-cpu0, evaluated_boxes=ntest,
                                                 unresolved_relative_volume=unresolved)
                    target = next(targets, None)
        if bound is None:
            continue
        aa, bb = bound
        i = max(range(len(axes)), key=lambda ix: (bb - aa)[axes[ix]])
        j = axes[i]
        m = midpoint(aa[j], bb[j])
        if depth == 48 or not aa[j] < m < bb[j]:
            finished.append((aa, bb))
            deferred_status[status] += 1
            continue
        bleft, aright = bb.copy(), aa.copy()
        bleft[j] = aright[j] = m
        if C is None:
            L = R = None
        else:
            L, R = split(C, i + (2 if cofactor is not None else 0))
        pending.append((aa, bleft, depth + 1, L, status))
        pending.append((aright, bb, depth + 1, R, status))
    deferred_status.update(status for *_, status in pending)
    cpu, wall = process_time() - cpu0, perf_counter() - wall0
    summary, boxes = metrics(pending, finished, lo, hi)
    # The unresolved set volume is an upper bound; the complement is proved excluded.
    summary.update(cpu_s=cpu, wall_s=wall, setup_cpu_s=setup, evaluated_boxes=ntest,
                   evaluated_dark_axes=evaluated_dark_axes, max_depth=max_depth,
                   tail_clip_bins=cofactor.tail_bins if cofactor is not None else None,
                   last_parent_status=dict(deferred_status), decisions=dict(reasons),
                   pivots=list(case.pivots) if method != 'direct' else [],
                   midpoint_chart_spectral_q_proxy=case.chart_q if method != 'direct' else None)
    if trace:
        summary['cpu_to_unresolved_volume_upper_bound'] = attained
    return summary, boxes


def perform(args):
    public, budget, tail_depth, trace = args
    rows, boxes = {}, {}
    for method in METHODS:
        rows[method], boxes[method] = run_arm(public, method, budget, tail_depth, trace)
    return public.case_id, public.cohort, rows, boxes


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--budget', type=float, default=1.0)
    p.add_argument('--workers', type=int, default=3)
    p.add_argument('--tail-depth', type=int, default=4)
    p.add_argument('--group', choices=('spectral', 'natural', 'aperture'), action='append')
    p.add_argument('--case', action='append')
    p.add_argument('--trace', action='store_true')
    p.add_argument('--metrics-only', action='store_true')
    p.add_argument('--limit', type=int)
    p.add_argument('--tag', default='pivot_old')
    a = p.parse_args()
    public = [ca for ca in load_cases() if (not a.group or ca.cohort in a.group)
              and (not a.case or ca.case_id in a.case)]
    if a.limit is not None:
        public = public[:a.limit]
    jobs = [(ca, a.budget, a.tail_depth, a.trace) for ca in public]
    results, material = {}, {}
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for name, cohort, rows, boxes in pool.map(perform, jobs):
            results[name] = dict(cohort=cohort, methods=rows)
            if not a.metrics_only:
                material.update({f'{name}/{k}': v for k, v in boxes.items()})
            print(name, {m: (v['evaluated_boxes'], round(v['unresolved_relative_volume'], 5),
                            round(v['mean_relative_hull_width'], 4), v['decisions']) for m, v in rows.items()}, flush=True)
    target = ROOT / 'results' / a.tag
    if not a.metrics_only:
        np.savez_compressed(target.with_name(target.name + '_regions.npz'), **material)
    target.with_suffix('.json').write_text(json.dumps(dict(cpu_cap_per_case_method_s=a.budget,
        input_scope='observed old reads, public geometry and public noise radii only',
        arb_prec_bits=96, comparison_arms=list(METHODS), cases=results), indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
