"""Equal-border 100/1000 box cost benchmark with direct Arb cross-checks."""
import argparse

from run import (
    PUBLIC, RESULTS, Contract, decode_box, encode_box, gate, get_freeze,
    load_geometry, members, write,
)
from certify import (
    FixedKernel, ball, clock, elapsed, full_ball, full_border, full_interval_border,
    herm, interval_inertia, residual_proposal, star, threshold,
)

import numpy as np
from flint import arb, ctx
from scipy.linalg import eigvalsh
from threadpoolctl import threadpool_limits


def ranges(geo, b, count, seed):
    """Positive distinct threshold vectors on 64 placements with varying endpoints."""
    rng = np.random.default_rng(seed)
    centers = members((2, 9, 2, 9))
    seen = set()
    boxes = []
    while len(boxes) < count:
        c = centers[int(rng.integers(len(centers)))]
        lo = int(rng.integers(18, 53))/32
        hi = lo + int(rng.integers(1, 21))/32
        q = Contract(geo.side, 'single', (2, 9, 2, 9), b, lo, hi, mode='amplitude')
        e, u = q.leaf_box(c)
        rec = encode_box(e, u, b)
        if rec['box_id'] not in seen:
            seen.add(rec['box_id'])
            boxes.append(dict(center=c, amplitude=[lo, hi], box=rec))
    return boxes


class DirectFull:
    """Full n-state decomposition with cached base, trace rows and input Grams."""

    def __init__(self, g, y, b, w):
        n = g.n
        self.K = herm(star(g.sensor)@g.sensor + w*(np.eye(n)-.5*b*(g.green+star(g.green))))
        self.L = star(g.sensor)@y + .5*b*w*g.incident
        self.C = herm(star(y)@y)
        self.T = np.vstack((np.eye(n), g.green))
        self.E = np.vstack((np.zeros_like(g.incident), g.incident))
        self.n, self.b, self.w = n, b, w

    def update(self, ell, cap):
        ids = np.flatnonzero((ell != 0) | (cap != self.b))
        rows = np.r_[ids, self.n+ids]
        k = len(ids)
        T, E = self.T[rows], self.E[rows]
        q = np.arange(k)
        change = np.zeros((2*k, 2*k))
        change[q, k+q] = change[k+q, q] = -self.w*(ell[ids]+cap[ids]-self.b)/2
        change[k+q, k+q] = self.w*ell[ids]*cap[ids]
        return (herm(self.K+star(T)@(change@T)),
                self.L-star(T)@(change@E), herm(self.C+star(E)@(change@E)))


def signed(a, shift=0):
    eigen = eigvalsh(herm(a), check_finite=False)
    near = 2.e-13 * max(abs(eigen))
    if np.any(abs(eigen) <= near):
        return None
    neg = np.count_nonzero(eigen < 0)
    return [len(eigen)-int(neg)+shift, int(neg)]


def benchmark(counts):
    freeze = get_freeze()
    setup = clock()
    geo = load_geometry(PUBLIC/'small_geometry.npz')
    data = freeze['cases'][0]
    with np.load(PUBLIC/f"{data['id']}.npz", allow_pickle=False) as z:
        Y, radius = z['Y'], float(z['radius'])
    acquisition = elapsed(setup)
    acquisition_log = __import__('json').loads((PUBLIC/'acquisition.json').read_text())
    b = freeze['benchmark_default_b']
    w = float(np.max(np.sum(abs(geo.sensor)**2, axis=0))*2.**freeze['weight_exponents'][0])
    p0 = clock()
    candidates = ranges(geo, b, max(counts), freeze['benchmark_random_seed'])
    collection = elapsed(p0)
    for count in counts:
        run_count(count, candidates[:count], geo, Y, radius, w, b, acquisition,
                  acquisition_log['small']['geometry_setup_with_persistence'], collection,
                  freeze['benchmark_full_arb_counts'][list(freeze['benchmark_box_counts']).index(count)],
                  freeze['benchmark_random_seed']+count)


def run_count(count, boxes, geo, Y, radius, w, b, io, construction, generation, cross_check, seed):
    endpoints = [decode_box(row['box'], geo.n) for row in boxes]
    start = clock()
    full = DirectFull(geo, Y, b, w)
    full_setup = elapsed(start)
    start = clock()
    fix = FixedKernel(geo, Y, b, w)
    fixed_setup = elapsed(start)

    fspec, fsign = [], []
    start = clock()
    for (ell, cap) in endpoints:
        K, L, C = full.update(ell, cap)
        fspec.append(residual_proposal(K, -L, C))
        fsign.append(signed(full_border(K, L, C, radius**2)))
    full_time = elapsed(start)

    rspec, rsign = [], []
    start = clock()
    for (ell, cap) in endpoints:
        R, N, H, k = fix.float_reduced(ell, cap)
        rspec.append(residual_proposal(R, N, H, k))
        border = np.block([[R, N], [star(N), H-(radius**2)*np.eye(H.shape[0])]])
        inertia = signed(border, geo.n-k)
        rsign.append([inertia[0], inertia[1]-k] if inertia is not None else None)
    fixed_time = elapsed(start)

    comparable, differences, margin_ambiguous = 0, [], 0
    for row, a, c, q, p in zip(boxes, fspec, rspec, fsign, rsign, strict=True):
        if a is None or c is None or q is None or p is None:
            margin_ambiguous += 1
            continue
        if q != p or a[0] != c[0]:
            raise AssertionError('equivalent borders disagree')
        comparable += 1
        differences.append(float(np.max(abs(a[1]-c[1]))))
        row['numeric_full_inertia_at_r_sq'] = q
        row['numeric_fixed_inertia_at_r_sq'] = p
        row['numeric_source_spectrum_max_abs_error'] = differences[-1]

    start = clock()
    fix.prepare_arb()
    arb_setup = elapsed(start)
    start = clock()
    for row, (ell, cap), hint in zip(boxes, endpoints, rspec, strict=True):
        j, t, _, _ = threshold(hint, radius)
        tested_t = t if t is not None else radius**2
        matrix, k = fix.interval_reduced(ell, cap, tested_t)
        proof = interval_inertia(matrix, full_count=True)
        strict = bool(j is not None and proof is not None and proof['verified_positive'] >= k+j
                      and arb(j)*arb(t) >= arb(radius)**2)
        row['fixed_arb'] = dict(j=j, t=tested_t, counts=proof, strictly_excluded=strict)
    arb_all = elapsed(start)

    start = clock()
    verified_cross_check = 0
    selected = np.random.default_rng(seed).choice(count, size=cross_check, replace=False)
    for idx in sorted(selected):
        row = boxes[idx]
        ell, cap = endpoints[idx]
        af = row['fixed_arb']
        g, s, bmat = geo.restrict_arb(np.arange(geo.n))
        exact = full_interval_border(*full_ball(g, s, bmat, ball(Y), ell, cap, arb(w)), af['t'])
        full_proof = interval_inertia(exact, full_count=True)
        rec = af['counts']
        if rec is not None and full_proof is not None:
            k = int(np.count_nonzero((ell != 0) | (cap != b)))
            if (full_proof['verified_positive'] != geo.n + rec['verified_positive']-k
                    or full_proof['verified_negative'] != rec['verified_negative']-k):
                raise AssertionError('independent exact full replay disagrees')
            verified_cross_check += 1
        row['independent_full_arb'] = full_proof
    cross_check_time = elapsed(start)

    strictly_excluded = sum(row['fixed_arb']['strictly_excluded'] for row in boxes)
    record = dict(count=count, genuinely_distinct_boxes=len({row['box']['box_id'] for row in boxes}),
                  shape_contract='single disk; spatially independent profile-scaled intervals, weak background',
                  geometry_states=geo.n, source_states=Y.shape[1], exact_trace_matrix_dimension=18+Y.shape[1],
                  b=b, common_sector_weight=w, optical_multiplier=0, precision_bits=ctx.prec,
                  all_box_timing_scope='both numeric arms use cached common full terms and verify all threshold signs; Arb fixed on all boxes; separately sampled independent direct Arb from original matrices',
                  direct_full_arb_scope='not extrapolated to all boxes',
                  common_geometry_setup_with_persistence=construction, common_reading_and_model_load=io,
                  candidate_generation_for_largest_batch=generation,
                  full_float_cold_setup=full_setup, fixed_float_cold_setup=fixed_setup,
                  full_float_all_assembly_factorizations=full_time,
                  fixed_float_all_reduced_examinations=fixed_time,
                  fixed_arb_cold_setup=arb_setup, fixed_arb_all_nodes=arb_all,
                  independent_direct_arb_replay=cross_check_time, independent_replayed_boxes=cross_check,
                  independent_matching_full_and_reduced_arb_inertias=verified_cross_check,
                  full_float_matching_boxes=comparable, numerically_ambiguous_boxes=margin_ambiguous,
                  max_numeric_source_spectrum_absolute_difference=max(differences, default=None),
                  fixed_arb_full_sign_counts=sum(row['fixed_arb']['counts'] is not None for row in boxes),
                  strictly_excluded=strictly_excluded,
                  remaining_without_witness_search=count-strictly_excluded,
                  boxes=boxes)
    write(RESULTS/f'benchmark_{count}.json', record)
    print('B benchmark', count, dict(full_numeric=full_time, fixed_numeric=fixed_time,
                                    fixed_arb=arb_all, full_arb_sample=cross_check_time,
                                    matching_arb=verified_cross_check, strict=strictly_excluded), flush=True)


if __name__ == '__main__':
    arg = argparse.ArgumentParser(description=__doc__)
    arg.add_argument('--counts', nargs='+', type=int, default=[100, 1000])
    parsed = arg.parse_args()
    gate()
    ctx.prec, ctx.threads = 128, 1
    with threadpool_limits(limits=1):
        benchmark(parsed.counts)
