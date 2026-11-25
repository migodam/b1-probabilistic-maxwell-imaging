"""Build and verify a full center cover from known geometry and measured ports."""
import os

for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[variable] = '1'

import argparse
from collections import Counter
from dataclasses import dataclass
from hashlib import blake2b
import json
from pathlib import Path
import resource
import sys

import numpy as np
from flint import arb, ctx
from scipy.linalg import eigvalsh
from scipy.special import hankel1, jv
from threadpoolctl import threadpool_info, threadpool_limits

from certify import (
    FixedKernel, Geometry, ball, certify_direct, clock, elapsed, float_forward,
    full_ball, full_border, full_float, full_interval_border, interval_inertia,
    residual_proposal, strict_witness, threshold, witness_search,
)


HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
STUDY = PAPER.parent
FREEZE = HERE / 'freeze.json'
PUBLIC = HERE / 'public'
RESULTS = HERE / 'results'
STATUS = ('EXCLUDED', 'WITNESS_FOUND', 'UNRESOLVED')


def write(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_suffix('.part')
    part.write_text(json.dumps(obj, indent=2, allow_nan=False, ensure_ascii=False)+'\n')
    part.replace(path)


def digest(*arrays):
    h = blake2b(digest_size=20)
    for a in arrays:
        a = np.ascontiguousarray(a)
        h.update(str(a.dtype).encode())
        h.update(str(a.shape).encode())
        h.update(a.view(np.uint8))
    return h.hexdigest()


def get_freeze():
    return json.loads(FREEZE.read_text())


def gate():
    path = PAPER / 'results/p0.json'
    if not path.is_file():
        raise SystemExit("Experiment B waits for the recorded P0 pass.")
    record = json.loads(path.read_text())
    if not (record.get('pass') is True or record.get('passed') is True
            or record.get('status') in ('PASS', 'PASSED')):
        raise SystemExit(f"Experiment B requires a pass in {path}. Keys: {list(record)}")


def coarse_geometry(config):
    n = config['side']
    dx = 2/n
    a = (np.arange(n)+.5)*dx-1
    yy, xx = np.meshgrid(a, a, indexing='ij')
    xy = np.c_[xx.ravel(), yy.ravel()]
    k = 2*np.pi*config['frequency_hz']/3e8
    radius = dx/np.sqrt(np.pi)
    theta = np.deg2rad(config['transmitter_degrees'])
    phi = np.deg2rad(config['receiver_degrees'])
    tx = np.c_[np.cos(theta), np.sin(theta)]
    rx = 3.75*np.c_[np.cos(phi), np.sin(phi)]
    B = np.exp(1j*k*(xy@tx.T))
    constant = .5j*np.pi*k*radius*jv(1, k*radius)
    S = constant*hankel1(0, k*np.linalg.norm(rx[:, None]-xy[None, :], axis=2))
    d = np.linalg.norm(xy[:, None]-xy[None, :], axis=2)
    np.fill_diagonal(d, 1.)
    G = constant*hankel1(0, k*d)
    np.fill_diagonal(G, .5j*np.pi*k*radius*hankel1(1, k*radius)-1)
    return Geometry(n, S, B, green=G)


def native_geometry(config):
    with np.load(STUDY/'synthetic/public/fine_geometry.npz', allow_pickle=False) as data:
        symbol = data['spectrum']
    archive = STUDY/'physical_revision/array/observations/arrays.npz'
    with np.load(archive, allow_pickle=False) as data:
        S = data['S32'][config['receiver_indices']]
        B = data['B32'][:, config['transmitter_indices']]
    return Geometry(96, S, B, spectrum=symbol)


def load_geometry(path):
    with np.load(path, allow_pickle=False) as z:
        return Geometry(int(z['side']), z['S'], z['B'],
                        green=z['G'] if 'G' in z.files else None,
                        spectrum=z['spectrum'] if 'spectrum' in z.files else None)


def grid_masks(side, shape, center):
    """Raster disks and annuli about an actual pixel, in grid coordinates."""
    yy, xx = np.indices((side, side))
    x, y = center
    if shape == 'double':
        squared = [(xx-(x+offset))**2+(yy-y)**2 for offset in (-2, 2)]
        core = (squared[0] <= 1) | (squared[1] <= 1)
        outer = (squared[0] <= 3) | (squared[1] <= 3)
    else:
        squared = (xx-x)**2+(yy-y)**2
        if shape == 'annulus':
            core = (squared >= 2) & (squared <= 4)
            outer = (squared >= 1) & (squared <= 6)
        elif shape == 'native_single':
            core, outer = squared <= 4, squared <= 9
        elif shape == 'single':
            core, outer = squared <= 1, squared <= 3
        else:
            raise ValueError(shape)
    return core.ravel(), outer.ravel()


def members(win):
    x0, x1, y0, y1 = win
    return [(x, y) for y in range(y0, y1+1) for x in range(x0, x1+1)]


def quadrants(win):
    x0, x1, y0, y1 = win
    mx, my = (x0+x1)//2, (y0+y1)//2
    xs = [(x0, mx), (mx+1, x1)] if x1 > x0 else [(x0, x1)]
    ys = [(y0, my), (my+1, y1)] if y1 > y0 else [(y0, y1)]
    return [(l, r, b, t) for b, t in ys for l, r in xs]


@dataclass(frozen=True)
class Contract:
    side: int
    shape: str
    window: tuple
    b: float
    low: float
    cap: float
    mode: str = 'core_outer'

    def leaf_box(self, center):
        core, outer = grid_masks(self.side, self.shape, center)
        ell = np.zeros(self.side**2)
        up = np.full(self.side**2, self.b)
        if self.mode == 'core_outer':
            ell[core], up[outer] = self.low, self.cap
        elif self.mode == 'amplitude':
            profile = np.where(core, 1., np.where(outer, .5, 0.))
            ell[outer] = self.low*profile[outer]
            up[outer] = self.cap*profile[outer]
        elif self.mode == 'arbitrary_high':
            return self.peak_box(center)
        else:
            raise ValueError(self.mode)
        return ell, up

    def peak_box(self, center):
        x, y = center
        ell = np.zeros(self.side**2)
        ell[y*self.side+x] = self.low
        return ell, np.full(self.side**2, self.cap)

    def cover_box(self, win):
        covers = [self.leaf_box(c) for c in members(win)]
        return np.min([ell for ell, _ in covers], axis=0), np.max([up for _, up in covers], axis=0)


def encode_box(ell, cap, b):
    indices = np.flatnonzero((ell != 0.) | (cap != b))
    return dict(default=[0., b], indices=indices.tolist(), ell=ell[indices].tolist(),
                upper=cap[indices].tolist(), box_id=digest(ell, cap))


def decode_box(record, n):
    a, b = record['default']
    ell, cap = np.full(n, a), np.full(n, b)
    ids = record['indices']
    ell[ids], cap[ids] = record['ell'], record['upper']
    if digest(ell, cap) != record['box_id']:
        raise ValueError('box endpoints changed')
    return ell, cap


def setup_design():
    freeze = get_freeze()
    PUBLIC.mkdir(exist_ok=True)
    OFFLINE = HERE/'offline'
    generators = json.loads((OFFLINE/'generator_parameters.json').read_text())
    summary = {}
    for name, config in freeze['grids'].items():
        since = clock()
        geom = native_geometry(config) if name == 'native' else coarse_geometry(config)
        source = geometry_receipt(geom, name, config)
        file = PUBLIC/f'{name}_geometry.npz'
        arrays = dict(side=geom.side, S=geom.sensor, B=geom.incident)
        arrays.update(G=geom.green) if geom.green is not None else arrays.update(spectrum=geom.spectrum)
        np.savez_compressed(file, **arrays)
        summary[name] = dict(**source, geometry_setup_with_persistence=elapsed(since))
    for case in freeze['cases']:
        started = clock()
        source = PUBLIC/f"{case['grid']}_geometry.npz"
        g = load_geometry(source)
        generator = generators[case['id']]
        rng = np.random.default_rng(generator['seed'])
        # Acquisition isolation: online() only opens the public arrays below.
        win = tuple(case['center_window'])
        eligible = [c for c in members(win) if win[0] < c[0] < win[1] and win[2] < c[1] < win[3]]
        center = eligible[int(rng.integers(len(eligible)))]
        core, outer = grid_masks(g.side, case['shape'], center)
        m = np.zeros(g.n)
        amplitude = generator['generator_amplitude_center'] + .0625*(rng.random(int(np.count_nonzero(outer)))-.5)
        m[outer] = amplitude
        m[outer & ~core] *= .5
        y = float_forward(g, m)
        raw = rng.normal(size=y.shape)+1j*rng.normal(size=y.shape)
        noise = raw * (case['fractional_noise']*np.linalg.norm(y)/np.linalg.norm(raw))
        observed = y+noise
        radius = float(np.nextafter(np.linalg.norm(noise)+1e-6*np.linalg.norm(y), np.inf))
        np.savez_compressed(PUBLIC/f"{case['id']}.npz", Y=observed, radius=radius)
        np.savez_compressed(OFFLINE/f"{case['id']}.npz", material=m, center=center)
        summary[case['id']] = dict(observation_id=case['id'], Y_hash=digest(observed, radius),
                                   radius=radius, clean_field_norm=float(np.linalg.norm(y)),
                                   acquisition_io_and_forward=elapsed(started))
    write(PUBLIC/'acquisition.json', summary)


def geometry_receipt(geom, name, config):
    model = geom.spectrum if geom.spectrum is not None else geom.green
    return dict(grid=[geom.side]*2, cells=geom.n, aperture=config, sensor_rows=geom.sensor.shape[0],
                sources_cols=geom.incident.shape[1], model_hash=digest(model, geom.sensor, geom.incident),
                units='m, Hz, dimensionless real contrast; exp(-i omega t); equal-area disk pulse scalar TM',
                source='native v5/v4 padded spectrum and S32/B32 arrays' if name == 'native'
                       else 'equal-area disk pulse kernel, same convention as archived paper_physics.acquisition',
                contract='exact binary64 stored Green, source and sensor entries' if name == 'small'
                         else 'exact DFT of stored binary64 padded spectrum, stored S/B entries exact')


class Engine:
    def __init__(self, geom, Y, radius, freeze):
        self.geom, self.Y, self.radius, self.freeze = geom, Y, radius, freeze
        self.weights = [float(np.max(np.sum(abs(geom.sensor)**2, axis=0)) * 2.**k)
                        for k in freeze['weight_exponents']]
        self.kernels = {}
        self.setup = {}
        self.witnesses = {}
        self.position_checks = []

    def kernel(self, b, w):
        key = (b, w)
        if key not in self.kernels:
            since = clock()
            self.kernels[key] = FixedKernel(self.geom, self.Y, b, w)
            self.setup[str(key)] = dict(float_setup=elapsed(since), arb_setup=None)
        return self.kernels[key]

    def assess(self, ell, up, b, *, vac_budget=None, max_weights=None, full=False):
        started = clock()
        chosen = self.weights[:max_weights] if max_weights else self.weights
        if self.geom.spectrum is not None or full:
            out = certify_direct(self.geom, self.Y, ell, up, self.radius, chosen,
                                 self.geom.n if vac_budget is None else vac_budget)
        else:
            attempts = []
            for w in chosen:
                kernel = self.kernel(b, w)
                # Record the genuinely cold Arb preparation only when a proposal can exclude.
                was_ready = kernel.strict is not None
                since = clock()
                row = kernel.certify(ell, up, self.radius)
                if not was_ready and kernel.strict is not None:
                    self.setup[str((b, w))]["arb_setup_and_first_node"] = elapsed(since)
                attempts.append(dict(weight=w, **row))
                if row["certificate"] is not None:
                    break
            out = dict(**attempts[-1], weights_tried=len(attempts))
            out['best_numeric_tail_over_radius'] = max(
                (a['numeric_tail_over_radius'] for a in attempts if a['numeric_tail_over_radius'] is not None),
                default=None)
            out['proposal_scan'] = [dict(weight=a['weight'], numeric_tail_over_radius=a['numeric_tail_over_radius'],
                                         numeric_state_negative=a['numeric_state_negative']) for a in attempts]
        return dict(**out, node_seconds=elapsed(started))

    def attempt_witness(self, name, center, ell, up):
        started = clock()
        try:
            material, evidence, cost = witness_search(self.geom, self.Y, ell, up, self.radius,
                                                      self.freeze['witness_max_evaluations'])
        except (ArithmeticError, np.linalg.LinAlgError, ValueError) as exc:
            return None, dict(reason=str(exc), **elapsed(started))
        if evidence:
            self.witnesses[name] = material
            evidence = dict(**evidence, material_array=name, material_hash=digest(material),
                            center=list(center))
        return evidence, dict(**cost, **elapsed(started))

    def certify_existing_witness(self, name, source, center, ell, up):
        material = self.witnesses[source]
        started = clock()
        evidence = strict_witness(self.geom, self.Y, ell, up, self.radius, material)
        if evidence is None:
            return None
        self.witnesses[name] = material
        return dict(**evidence, material_array=name, material_hash=digest(material),
                    center=list(center), verified_source_candidate=source,
                    independent_forward_replay=elapsed(started))

    def tree(self, contract, name, witness=True, vac_budget=None, inherited=None):
        nodes, final = [], {}
        target_centers = members(contract.window)
        witness_checks = 0

        def visit(win, ancestry, parent=None):
            nonlocal witness_checks
            lid = f'{name}__{len(nodes):04d}'
            cset = members(win)
            lo, up = contract.cover_box(win)
            if parent is not None:
                if not np.all(parent[0] <= lo) or not np.all(up <= parent[1]):
                    raise AssertionError('child enclosure exceeds parent')
            proof = self.assess(lo, up, contract.b, vac_budget=vac_budget)
            row = dict(id=lid, parent_id=ancestry, contract=name, center_window=list(win),
                       covered_centers=len(cset), box=encode_box(lo, up, contract.b),
                       status='UNRESOLVED', **proof)
            nodes.append(row)
            if proof['certificate'] is not None:
                row['status'] = 'EXCLUDED'
                final.update({c: dict(status='EXCLUDED', backing_node=lid) for c in cset})
                return
            if len(cset) == 1:
                center = cset[0]
                budget = self.freeze['witness_centers_per_contract']
                if inherited and center in inherited:
                    witness_result = self.certify_existing_witness(lid, inherited[center], center, lo, up)
                    if witness_result is None:
                        raise AssertionError('nested weak box lost a valid same-center witness')
                    row['witness'] = witness_result
                    row['status'] = 'WITNESS_FOUND'
                elif witness and witness_checks < budget:
                    witness_checks += 1
                    witness_result, count = self.attempt_witness(lid, center, lo, up)
                    row['witness_search'] = count
                    if witness_result:
                        row['witness'] = witness_result
                        row['status'] = 'WITNESS_FOUND'
                final[center] = dict(status=row['status'], backing_node=lid)
                return
            children = quadrants(win)
            actual = [c for child in children for c in members(child)]
            if Counter(actual) != Counter(cset):
                raise AssertionError('children must partition center labels')
            for child in children:
                visit(child, lid, (lo, up))

        tick = clock()
        visit(contract.window, None)
        if set(final) != set(target_centers):
            raise AssertionError('center cover incomplete')
        count = Counter(row['status'] for row in final.values())
        return dict(id=name, family=contract.mode, shape=contract.shape, b=contract.b, low=contract.low,
                    cap=contract.cap, status_counts={k:count[k] for k in STATUS},
                    nodes=nodes, final_centers=[dict(center=c, **final[c]) for c in target_centers],
                    distinct_boxes=len(set(n['box']['box_id'] for n in nodes)),
                    total_centers=len(target_centers), witnesses_attempted=witness_checks,
                    elapsed=elapsed(tick))

    def arbitrary_high(self, case):
        freeze = self.freeze
        b, u = freeze['peak_test_threshold'], freeze['peak_global_cap']
        n = self.geom.side
        contract = Contract(n, case['shape'], (0, n-1, 0, n-1), u, b, u,
                            mode='arbitrary_high')
        # Union of full native boxes. Point witnesses depend on this singleton's
        # actual cap event even if a numerical candidate fits several other cells.
        rows = []
        tick = clock()
        peak_budget = freeze['peak_witness_centers']
        checks = 0
        previously_certified = list(self.witnesses)
        for ix, center in enumerate(members(contract.window)):
            ell, up = contract.peak_box(center)
            proof = self.assess(ell, up, u, max_weights=1, full=True)
            wid = f"peak_{ix:03d}"
            row = dict(id=wid, center=center, tested_cell=center[1]*n+center[0],
                       box=encode_box(ell, up, u), status='UNRESOLVED', **proof)
            if proof['certificate']:
                row['status'] = 'EXCLUDED'
            else:
                for source in previously_certified:
                    m = self.witnesses[source]
                    if m[row['tested_cell']] >= b and np.all(m <= u):
                        rec = self.certify_existing_witness(wid, source, center, ell, up)
                        if rec:
                            row['witness'] = rec
                            row['status'] = 'WITNESS_FOUND'
                            break
                # Deterministic additional proposals on a stratified budget.
                if row['status'] == 'UNRESOLVED' and checks < peak_budget:
                    if ix % max(1, n*n//peak_budget) == 0:
                        checks += 1
                        rec, cost = self.attempt_witness(wid, center, ell, up)
                        row['witness_search'] = cost
                        if rec:
                            row['witness'] = rec
                            row['status'] = 'WITNESS_FOUND'
            rows.append(row)
        count = Counter(r['status'] for r in rows)
        # Different quantifier: max over entire grid <= q describes one box.
        peak_cap = freeze['local_peak_upper_box_cap']
        ell = np.zeros(self.geom.n)
        cap = np.full(self.geom.n, peak_cap)
        global_low_box = self.assess(ell, cap, peak_cap, full=True)
        return dict(cell_event=f'exists a specified cell with chi_i >= {b}; all other cells in [0,{u}]',
                    union_complete=len(rows) == self.geom.n,
                    status_counts={k:count[k] for k in STATUS}, threshold=b, global_cap=u,
                    witness_checks=checks, rows=rows,
                    global_low_peak_box=dict(cap=peak_cap, box=encode_box(ell, cap, peak_cap),
                                             **global_low_box),
                    elapsed=elapsed(tick))


def online(selected):
    freeze = get_freeze()
    public_log = json.loads((PUBLIC/'acquisition.json').read_text())
    for case in freeze['cases']:
        ident = case['id']
        if selected and ident not in selected:
            continue
        source = PUBLIC/f"{case['grid']}_geometry.npz"
        started = clock()
        geom = load_geometry(source)
        with np.load(PUBLIC/f'{ident}.npz', allow_pickle=False) as data:
            Y, radius = data['Y'], float(data['radius'])
        e = Engine(geom, Y, radius, freeze)
        prepared = elapsed(started)
        native = case['grid'] == 'native'
        started_exact = clock()
        if native:
            geom.prepare_arb()
        native_exact_dft_time = elapsed(started_exact) if native else None
        shape = case['shape']
        win = tuple(case['center_window'])
        b0, b1 = freeze['center_b']
        u, U = ((case['native_core_lower'], case['native_local_cap']) if native
                else (freeze['core_lower'], freeze['local_cap']))
        trees = []
        bands = []
        if native:
            core = Contract(geom.side, shape, win, b0, u, U)
            trees.append(e.tree(core, 'native_vacuum', vac_budget=freeze['native_max_retained_rows']))
        else:
            inherited = {}
            for b in (b0, b1):
                core = Contract(geom.side, shape, win, b, u, U)
                trees.append(e.tree(core, f'core_outer_b{b:.5g}', inherited=inherited))
                inherited = {tuple(row['witness']['center']): row['id'] for row in trees[-1]['nodes']
                             if row['status'] == 'WITNESS_FOUND'}
            for i, (lo, hi) in enumerate(freeze['amplitude_bins']):
                amplitude = Contract(geom.side, shape, win, b0, lo, hi, mode='amplitude')
                bands.append(e.tree(amplitude, f'amplitude_{i}'))
        peak = e.arbitrary_high(case) if not native else dict(
            union_complete=False, reason='native single-cell union has 9216 full-domain weak/strong boxes; outside this budget',
            points_attempted=0, unresolved_hypotheses=geom.n)
        result = dict(id=ident, object_id=case['object_id'], grid=case['grid'],
                      conditions_aggregated=0, center_trees=trees, amplitude_trees=bands,
                      arbitrary_high=peak, operator=public_log[case['grid']],
                      observation=public_log[ident], online_io_and_numeric_initialization=prepared,
                      common_kernel_setups=e.setup, native_arb_idft_setup=native_exact_dft_time,
                      total_online=elapsed(started), weights=e.weights,
                      peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                      pool_threads=threadpool_info())
        write(RESULTS/f'{ident}.json', result)
        np.savez_compressed(RESULTS/f'{ident}_witnesses.npz', **e.witnesses)
        print(ident, [r['status_counts'] for r in trees], [r['status_counts'] for r in bands],
              peak.get('status_counts'), result['total_online'], flush=True)


def score(selected):
    freeze = get_freeze()
    output = {}
    for case in freeze['cases']:
        ident = case['id']
        if selected and ident not in selected:
            continue
        result_file = RESULTS/f'{ident}.json'
        result = json.loads(result_file.read_text())
        with np.load(HERE/'offline'/f'{ident}.npz', allow_pickle=False) as data:
            mat, center = data['material'], data['center'].tolist()
        geom = load_geometry(PUBLIC/f"{case['grid']}_geometry.npz")
        with np.load(PUBLIC/f'{ident}.npz', allow_pickle=False) as data:
            y, r = data['Y'], float(data['radius'])
        q = strict_witness(geom, y, np.zeros(geom.n), np.full(geom.n, float(mat.max())), r, mat)
        if q is None:
            raise AssertionError('generator fails independent strict forward replay')
        core_results = {}
        shape_membership = {}
        for tree in result['center_trees']+result['amplitude_trees']:
            tree_status = next(n for n in tree['final_centers'] if n['center'] == center)
            mode = 'amplitude' if tree['id'].startswith('amplitude') else 'core_outer'
            contract = Contract(geom.side, case['shape'], tuple(case['center_window']),
                                tree['b'], tree['low'], tree['cap'], mode=mode)
            ell, cap = contract.leaf_box(center)
            feasible_material_in_box = bool(np.all(ell <= mat) and np.all(mat <= cap))
            if feasible_material_in_box and tree_status['status'] == 'EXCLUDED':
                raise AssertionError('strict exclusion at a legal generating pair')
            core_results[tree['id']] = dict(center_status=tree_status['status'],
                                            generator_in_box=feasible_material_in_box)
            if mode == 'amplitude':
                shape_membership[tree['id']] = feasible_material_in_box
        for row in result['arbitrary_high'].get('rows', []):
            ell, up = decode_box(row['box'], geom.n)
            if mat[row['tested_cell']] >= result['arbitrary_high']['threshold'] and np.all(mat <= up):
                if row['status'] == 'EXCLUDED':
                    raise AssertionError('strict exclusion at a legal high cell')
        output[ident] = dict(center=center, material_hash=digest(mat), generator_max=float(mat.max()),
                             strict_generator=q, status_at_generator_center=core_results,
                             amplitude_membership=shape_membership)
    write(RESULTS/'offline_score.json', output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('prepare', 'online', 'score'))
    parser.add_argument('select', nargs='*', help='optional complete observation IDs')
    args = parser.parse_args()
    gate()
    ctx.prec, ctx.threads = 128, 1
    with threadpool_limits(limits=1):
        if args.mode == 'prepare':
            setup_design()
        elif args.mode == 'online':
            online(args.select)
        else:
            score(args.select)


if __name__ == '__main__':
    main()
