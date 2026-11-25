"""Scan translated material boxes on frozen Maxwell state and source spaces."""
import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import eigh, eigvalsh, solve
from scipy.sparse.linalg import svds
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
LEGACY = STUDY.parent / 'z_alg_secant_tasks/highdim'
NATIVE = STUDY / 'physical_revision/array/observations'
FINITE = STUDY / 'synthetic/public'
CONTROL = STUDY / 'ip_unified/numerics/feasibility/strict_endpoints'
CONTROL_ID = 'control_disk075__pw360_32x32__pair2pct'
SMALL_ID = 'shiftbox32_01__pw360_10x20__n02'
LEGACY_ID = 'r10_twod_0001'
NATIVE_IDS = ('optfresh_01__pw360_32x32__pair2pct', 'optfresh_02__pw360_32x32__pair2pct')
WIDTHS = (.25, .5, .75, 1., 1.5, 2.)
PEAKS = (2.5, 4.)
FACTORS = tuple(2.**k for k in (-16, -12, -8, -4, 0))


def clocks():
    return perf_counter(), process_time()


def timing(start):
    return dict(zip(('wall_seconds', 'cpu_seconds'),
                    (a - b for a, b in zip(clocks(), start, strict=True)), strict=True))


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    temp = path.with_suffix('.part')
    temp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    temp.replace(path)


def star(x):
    return x.conj().T


def herm(x):
    return (x + star(x)) * .5


def coordinates(side):
    yy, xx = (np.indices((side, side)) + .5) * (2/side) - 1.
    return np.column_stack((xx.ravel(), yy.ravel()))


@dataclass
class Geometry:
    S: np.ndarray
    B: np.ndarray
    side: int
    G: np.ndarray | None = None
    spectrum: np.ndarray | None = None

    def index(self, active):
        y, x = np.divmod(active, self.side)
        pad = self.spectrum.shape[0]
        return ((y[:, None] - y[None, :]) % pad, (x[:, None] - x[None, :]) % pad)

    def restricted(self, active):
        if self.G is None:
            G = np.fft.ifft2(self.spectrum)[self.index(active)]
        else:
            G = self.G[np.ix_(active, active)]
        return G, self.S[:, active], self.B[active]


@dataclass
class Observation:
    id: str
    geometry: Geometry
    Y: np.ndarray
    radius: float
    source: str


def native_geometry():
    with np.load(FINITE / 'fine_geometry.npz', allow_pickle=False) as data:
        spectrum = data['spectrum']
    with np.load(NATIVE / 'arrays.npz', allow_pickle=False) as data:
        return Geometry(data['S32'], data['B32'], 96, spectrum=spectrum)


def small_geometry():
    with np.load(FINITE / 'coarse_green.npz', allow_pickle=False) as data:
        G = data['G']
    with np.load(FINITE / 'pw360_10x20.npz', allow_pickle=False) as data:
        return Geometry(data['S'], data['B'], 32, G=G)


def make_small_observation(out):
    """The generating material is isolated in a separate offline file."""
    start = clocks()
    geo = small_geometry()
    x, y = coordinates(geo.side).T
    first = ((x + .18)/.30)**2 + ((y - .09)/.20)**2 <= 1
    cut = (x > -.26) & (x < -.15) & (y > .04)
    second = ((x - .29)/.17)**2 + ((y + .13)/.24)**2 <= 1
    m = np.maximum(3.25 * (first & ~cut) * (1. - .2 * (y > .1)), 2.25 * second)
    ids = np.flatnonzero(m)
    G, S, B = geo.restricted(ids)
    A = np.eye(len(ids)) - m[ids, None] * G
    J = solve(A, m[ids, None] * B, check_finite=False)
    clean = S @ J
    rng = np.random.default_rng(202611251403)
    noise = rng.standard_normal(clean.shape) + 1j * rng.standard_normal(clean.shape)
    Y = clean + .02 * np.linalg.norm(clean) * noise / np.linalg.norm(noise)
    r = np.nextafter(np.linalg.norm(Y - clean) + 1.e-6 * np.linalg.norm(Y), np.inf)
    np.savez_compressed(out / (SMALL_ID + '.npz'), Y=Y)
    np.savez_compressed(out / (SMALL_ID + '__offline.npz'), material=m, Y_clean=clean)
    row = dict(id=SMALL_ID, radius=float(r), noise_seed=202611251403, solve_active_cells=len(ids),
               source_rhs=B.shape[1], noise_denominator='||binary64 clean field||_F',
               noise_fraction=.02, forward_guard=1.e-6, **timing(start))
    write(out / (SMALL_ID + '.json'), row)


def load(out, which=None):
    wanted = {LEGACY_ID, SMALL_ID, *NATIVE_IDS, CONTROL_ID} if which is None else set(which)
    result = {}
    if LEGACY_ID in wanted:
        with np.load(LEGACY / 'public_vie_inputs.npz', allow_pickle=False) as data:
            green_key = str(data[f'{LEGACY_ID}__G_key'].item())
            geo = Geometry(data[f'{LEGACY_ID}__S'], data[f'{LEGACY_ID}__E0'], 24, G=data[green_key])
        with np.load(LEGACY / 'results' / (LEGACY_ID + '.npz'), allow_pickle=False) as data:
            Y = data['Y_archive']
        r = read(LEGACY / 'observation_radii.json')[LEGACY_ID]['frozen_observation_ball']
        result[LEGACY_ID] = Observation(LEGACY_ID, geo, Y, r, 'z_alg_secant_tasks/highdim')
    if SMALL_ID in wanted:
        if not (out / (SMALL_ID + '.npz')).exists():
            make_small_observation(out)
        with np.load(out / (SMALL_ID + '.npz'), allow_pickle=False) as data:
            Y = data['Y']
        r = read(out / (SMALL_ID + '.json'))['radius']
        result[SMALL_ID] = Observation(SMALL_ID, small_geometry(), Y, r, 'boxes same-grid forward')
    if wanted.intersection((*NATIVE_IDS, CONTROL_ID)):
        geo = native_geometry()
        rows = {v['id']: v for v in read(NATIVE / 'manifest.json')['observations']}
        for key in NATIVE_IDS:
            if key not in wanted:
                continue
            with np.load(NATIVE / rows[key]['file'], allow_pickle=False) as data:
                Y = data['Y']
            result[key] = Observation(key, geo, Y, rows[key]['padded_radius'], 'v3 original native array')
        if CONTROL_ID in wanted:
            with np.load(CONTROL / (CONTROL_ID + '.npz'), allow_pickle=False) as data:
                Y = data['Y']
            r = read(CONTROL / (CONTROL_ID + '.json'))['radius']
            result[CONTROL_ID] = Observation(CONTROL_ID, geo, Y, r, 'unified physical disk control')
    return result


def placements(obs):
    geo = obs.geometry
    xy = coordinates(geo.side)
    score = np.real(np.einsum('ij,ji->i', star(geo.S) @ obs.Y, star(geo.B)))
    window = np.flatnonzero(np.max(np.abs(xy), axis=1) <= .40)
    match = xy[window[np.argmax(score[window])]]
    return {'origin': np.zeros(2), 'born_peak': match}, score


def box(geo, center, radius, peak, width):
    rad2 = np.sum((coordinates(geo.side) - center)**2, axis=1)
    active = np.flatnonzero(rad2 <= radius**2 + 1.e-15)
    center_m = peak * (1. - .25 * rad2[active] / radius**2)
    lower = np.maximum(0., center_m - width)
    upper = center_m + width
    return active, lower, upper


def ball_parameters(lower, upper):
    # Binary64 proposals.  Strict proof recomputes both from exact endpoints.
    return (lower + upper) * .5, (upper - lower) * .5


def arrays(G, S, B, Y, m, r, W=None):
    if W is None:
        W = np.ones(len(G))
    A = np.eye(len(G)) - m[:, None] * G
    RG, RB = r[:, None] * G, r[:, None] * B
    F = m[:, None] * B
    C = herm(star(A) @ (W[:, None] * A) - star(RG) @ (W[:, None] * RG))
    P = star(A) @ (W[:, None] * F) + star(RG) @ (W[:, None] * RB)
    D = herm(star(F) @ (W[:, None] * F) - star(RB) @ (W[:, None] * RB))
    return A, RG, C, P, D, herm(star(S) @ S), star(S) @ Y, herm(star(Y) @ Y)


def scalar_dimensions(S2, A):
    return float(eigvalsh(S2, subset_by_index=(len(S2) - 1, len(S2) - 1),
                         check_finite=False)[0] * len(A) / np.linalg.norm(A)**2)


def scan_box(obs, active, lower, upper, count, *, factors=FACTORS, W=None):
    start = clocks()
    G, S, B = obs.geometry.restricted(active)
    m, r = ball_parameters(lower, upper)
    A, RG, C, P, D, S2, SY, YY = arrays(G, S, B, obs.Y, m, r, W)
    T = solve(A.T, RG.T, check_finite=False).T
    count['double_dense_A_solves'] += 1
    count['double_dense_A_RHS'] += len(A)
    if len(A) < 160:
        sigma = svds(T, k=1, which='LM', return_singular_vectors=False, rng=117)[0]
    else:
        # The largest singular value is obtained without forming the source
        # Gram. A separate Arb vector test is required for a rigorous crossing.
        sigma = svds(T, k=1, which='LM', return_singular_vectors=False, rng=117, tol=1.e-9)[0]
    count['double_T_largest_singular_value'] += 1
    scale = scalar_dimensions(S2, A)
    by_mu = []
    for factor in factors:
        mu = float(factor * scale)
        K, L = herm(S2 + mu * C), SY + mu * P
        d, V = eigh(K, check_finite=False, driver='evd')
        count['double_K_eigensolves'] += 1
        gap = float(np.min(np.abs(d)))
        tolerance = 128 * len(d) * np.finfo(float).eps * float(max(abs(d[0]), abs(d[-1])))
        row = dict(factor=factor, mu=mu, nu=int(np.count_nonzero(d < 0)), K_gap_numeric=gap,
                   H_eigenvalues_descending=None, distance_over_radius_numeric=None)
        if gap > tolerance:
            if row['nu'] < obs.Y.shape[1]:
                X = V @ ((star(V) @ L) / d[:, None])
                H = herm(YY + mu * D - star(L) @ X)
                h = eigvalsh(H, check_finite=False)[::-1]
                count['double_H_eigensolves'] += 1
                count['stationary_source_rhs'] += obs.Y.shape[1]
                row['H_eigenvalues_descending'] = h.tolist()
                row['distance_over_radius_numeric'] = float(
                    np.sqrt(np.maximum(h[row['nu']:], 0.).sum()) / obs.radius)
            else:
                row['distance_over_radius_numeric'] = 0.
        by_mu.append(row)
    return dict(cells=len(active), peak_lower=float(lower.max()), peak_upper=float(upper.max()),
                norm_T_numeric=float(sigma), by_mu=by_mu, **timing(start))


def table(path, cases):
    fields = ('observation', 'placement', 'radius_fraction', 'peak', 'width', 'cells',
              'peak_lower', 'peak_upper', 'norm_T_numeric', 'max_indefinite_numeric',
              'max_positive_numeric', 'box_cpu_seconds')
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for name, case in cases.items():
            for b in case['boxes']:
                row = {k: b[k] for k in fields if k in b}
                row.update(observation=name, box_cpu_seconds=b['cpu_seconds'])
                for i in (0, 1):
                    values = [v['distance_over_radius_numeric'] for v in b['by_mu']
                              if v['distance_over_radius_numeric'] is not None and (v['nu'] > 0) == bool(i)]
                    row['max_indefinite_numeric' if i else 'max_positive_numeric'] = max(values, default=None)
                w.writerow(row)


def main():
    from collections import Counter
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, default=HERE / 'results')
    ap.add_argument('--case', action='append', help='exact observation ID')
    ap.add_argument('--scan-file', default='scan.json')
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    start = clocks()
    cases, maps = {}, {}
    for name, obs in load(args.out, args.case).items():
        started = clocks()
        loc, score = placements(obs)
        maps[name + '__adjoint_score'] = score.reshape(obs.geometry.side, obs.geometry.side)
        maps[name + '__Y'] = obs.Y
        radii = (.45, .62) if obs.geometry.side == 24 else (.245, .32)
        tests = [(label, rad, peak, width) for label in loc for rad in radii
                 for peak in PEAKS for width in WIDTHS]
        if name == CONTROL_ID:
            tests += [('born_peak', .40, None, .375)]
        count, boxes = Counter(), []
        for label, rad, peak, width in tests:
            if peak is None:
                ids, lower, upper = box(obs.geometry, loc[label], rad, .375, .375)
                lower[:], upper[:] = 0., .75
            else:
                ids, lower, upper = box(obs.geometry, loc[label], rad, peak, width)
            block = scan_box(obs, ids, lower, upper, count)
            short = [r for r in block['by_mu'] if r['distance_over_radius_numeric'] is not None]
            best = max((r['distance_over_radius_numeric'] for r in short), default=0.)
            biggest_indef = max((r['distance_over_radius_numeric'] for r in short if 0 < r['nu']), default=0.)
            print(name, label, rad, peak, width, 'norm', f"{block['norm_T_numeric']:.4g}",
                  'best', f'{best:.4g}', 'indef', f'{biggest_indef:.4g}', flush=True)
            boxes.append(dict(placement=label, center=loc[label].tolist(), radius_fraction=rad,
                              peak=peak, width=width, **block))
        cases[name] = dict(side=obs.geometry.side, radius=obs.radius, source=obs.source,
                           port_shape=list(obs.Y.shape), boxes=boxes, operations=dict(count), **timing(started))
        write(args.out / args.scan_file, dict(cases=cases, **timing(start)))
    table(args.out / Path(args.scan_file).with_suffix('.csv'), cases)
    np.savez_compressed(args.out / (Path(args.scan_file).stem + '_maps.npz'), **maps)


if __name__ == '__main__':
    with threadpool_limits(1):
        main()
