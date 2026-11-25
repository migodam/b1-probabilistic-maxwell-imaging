"""Re-read fixed inputs and replay saved spatial certificates at Arb128."""
from run import (
    HERE, PUBLIC, RESULTS, Contract, decode_box, digest, gate, get_freeze,
    load_geometry, members, quadrants, write,
)
from certify import (
    FixedKernel, ball, clock, dagger, elapsed, full_ball, full_interval_border,
    identity, interval_inertia, strict_witness,
)

from collections import Counter
import json

import numpy as np
from flint import arb, ctx
from threadpoolctl import threadpool_limits


def verify_cover(tree, cfg, n):
    family = Contract(n, cfg['shape'], tuple(cfg['center_window']), tree['b'], tree['low'],
                      tree['cap'], tree['family'])
    by_id = {row['id']: row for row in tree['nodes']}
    children = {}
    for row in tree['nodes']:
        children.setdefault(row['parent_id'], []).append(row)
        lo, hi = decode_box(row['box'], n*n)
        expected = family.cover_box(tuple(row['center_window']))
        if not np.array_equal(lo, expected[0]) or not np.array_equal(hi, expected[1]):
            raise AssertionError('node box is not actual memberwise enclosure')
    roots = children[None]
    if len(roots) != 1 or tuple(roots[0]['center_window']) != family.window:
        raise AssertionError('incomplete root window')
    final = [(tuple(row['center']), row['status']) for row in tree['final_centers']]
    if Counter(c for c, _ in final) != Counter(members(family.window)):
        raise AssertionError('incomplete final centers')
    for row in tree['nodes']:
        window = tuple(row['center_window'])
        if row['status'] == 'EXCLUDED' or len(members(window)) == 1:
            if row['id'] in children:
                raise AssertionError('children underneath a terminal node')
        elif sorted(tuple(c['center_window']) for c in children[row['id']]) != sorted(quadrants(window)):
            raise AssertionError('missing quadtree branches')
    for row in tree['final_centers']:
        node = by_id[row['backing_node']]
        window = tuple(node['center_window'])
        if row['status'] != node['status'] or tuple(row['center']) not in members(window):
            raise AssertionError('map label has no terminal proof')
        if row['status'] == 'WITNESS_FOUND' and row['center'] != node['witness']['center']:
            raise AssertionError('center label borrows witness from another center')


def direct_replay(geom, ell, cap, Y, proof, full_state, frames):
    ids = np.arange(geom.n) if full_state else np.flatnonzero(cap)
    g, s, b = geom.restrict_arb(ids)
    K, L, C = full_ball(g, s, b, ball(Y), ell[ids], cap[ids], arb(proof['weight']))
    try:
        state_frames = {}
        state = interval_inertia(K, full_count=True, frames=state_frames)
        if state:
            j = state['verified_negative']+proof['j']
            H = C-dagger(L)*K.solve(L)-arb(proof['t'])*identity(C.nrows())
            source_frames = {}
            trace = interval_inertia(H, j, frames=source_frames)
            if trace:
                frames.update({f'full_state_{k}':v for k,v in state_frames.items()})
                frames.update({f'full_source_{k}':v for k,v in source_frames.items()})
                return dict(method='independent_full_state_and_source_schur',
                            full_state_dimension=len(ids), verified_full_positive_at_least=len(ids)+proof['j'],
                            full_state_inertia=state, source_positive=trace)
    except (ValueError, ZeroDivisionError):
        pass
    border = full_interval_border(K, L, C, proof['t'])
    border_frames = {}
    fallback = interval_inertia(border, len(ids)+proof['j'], frames=border_frames)
    frames.update({f'full_border_{k}':v for k,v in border_frames.items()})
    return None if fallback is None else dict(method='independent_full_border',
                                              full_state_dimension=len(ids), **fallback)


def replay():
    cfg = get_freeze()
    out = {}
    frames = {}
    for case_index, case in enumerate(cfg['cases']):
        started = clock()
        ident = case['id']
        archive = json.loads((RESULTS/f'{ident}.json').read_text())
        geom = load_geometry(PUBLIC/f"{case['grid']}_geometry.npz")
        with np.load(PUBLIC/f'{ident}.npz', allow_pickle=False) as file:
            Y, radius = file['Y'], float(file['radius'])
        if digest(Y, radius) != archive['observation']['Y_hash']:
            raise AssertionError('observations differ from certificate')
        model = geom.spectrum if geom.spectrum is not None else geom.green
        if digest(model, geom.sensor, geom.incident) != archive['operator']['model_hash']:
            raise AssertionError('operator differs from certificate')
        tick = clock()
        geom.prepare_arb()
        setup_time = elapsed(tick)
        native = case['grid'] == 'native'
        tree_rows = archive['center_trees']+archive['amplitude_trees']
        for tree in tree_rows:
            verify_cover(tree, case, geom.side)
        rows = [row for tree in tree_rows for row in tree['nodes']]
        peak = archive['arbitrary_high']
        if peak.get('union_complete'):
            points = [row['tested_cell'] for row in peak['rows']]
            if sorted(points) != list(range(geom.n)):
                raise AssertionError('partial arbitrary-high enumeration')
            for row in peak['rows']:
                ell, cap = decode_box(row['box'], geom.n)
                expected = np.zeros(geom.n)
                expected[row['tested_cell']] = peak['threshold']
                if not np.array_equal(ell, expected) or not np.all(cap == peak['global_cap']):
                    raise AssertionError('arbitrary-high event has wrong voxel constraints')
            rows += peak['rows']
            rows += [dict(id='global_low_peak_box', **peak['global_low_peak_box'])]
        excluded = [row for row in rows if row['certificate'] is not None]
        sample = set(np.random.default_rng(cfg['benchmark_random_seed']+1_000+case_index).choice(
            len(excluded), size=min(12, len(excluded)), replace=False))
        notes, direct, witnesses = [], [], []
        kernels = {}
        with np.load(RESULTS/f'{ident}_witnesses.npz', allow_pickle=False) as witness_set:
            for i, row in enumerate(excluded):
                weights = row['certificate']
                if not arb(weights['j'])*arb(weights['t']) >= arb(radius)**2:
                    raise AssertionError('threshold fails exact noise radius')
                ell, cap = decode_box(row['box'], geom.n)
                tick = clock()
                node_frames = {}
                if weights['algorithm'] == 'exact_fixed_trace_augmented_border':
                    key = (row['box']['default'][1], weights['weight'])
                    if key not in kernels:
                        kernels[key] = FixedKernel(geom, Y, *key)
                    A, k = kernels[key].interval_reduced(ell, cap, weights['t'])
                    fixed_frames = {}
                    proof = interval_inertia(A, k+weights['j'], frames=fixed_frames)
                    node_frames.update({f'fixed_{k}':v for k,v in fixed_frames.items()})
                    if proof is None:
                        raise AssertionError(f"reduced point replay failed {row['id']}")
                    notes.append(dict(id=row['id'], **proof, **elapsed(tick)))
                    if i not in sample:
                        frames.update({f"{ident}__{row['id']}__{k}":v for k,v in node_frames.items()})
                        continue
                proof = direct_replay(geom, ell, cap, Y, weights, full_state=not native, frames=node_frames)
                frames.update({f"{ident}__{row['id']}__{k}":v for k,v in node_frames.items()})
                if proof is None:
                    raise AssertionError(f"independent full replay failed {row['id']}")
                direct.append(dict(id=row['id'], **proof, **elapsed(tick)))
            for row in rows:
                q = row.get('witness')
                if not q:
                    continue
                tick = clock()
                ell, cap = decode_box(row['box'], geom.n)
                material = witness_set[q['material_array']]
                checked = strict_witness(geom, Y, ell, cap, radius, material)
                if checked is None or digest(material) != q['material_hash']:
                    raise AssertionError('saved witness failed independent re-evaluation')
                witnesses.append(dict(id=row['id'], **checked, **elapsed(tick)))
        out[ident] = dict(precision_bits=ctx.prec, independent_full_replays=len(direct),
                          fixed_trace_replays=len(notes), excluded_nodes_replayed=len(excluded),
                          witnesses_replayed=len(witnesses), exact_dft_setup=setup_time if native else None,
                          fixed_rows=notes, full_rows=direct, witness_rows=witnesses,
                          all_contract_cover_checks_passed=True, timing=elapsed(started))
        print(ident, 'replayed', len(excluded), 'excluded nodes; direct', len(direct),
              '; feasible materials', len(witnesses), flush=True)
    write(RESULTS/'independent_replay.json', out)
    np.savez_compressed(RESULTS/'independent_replay_frames.npz', **frames)


if __name__ == '__main__':
    gate()
    ctx.prec, ctx.threads = 128, 1
    with threadpool_limits(limits=1):
        replay()
