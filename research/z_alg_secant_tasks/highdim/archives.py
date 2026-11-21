"""Read the original Maxwell acquisitions and legally available solver outputs."""

from dataclasses import dataclass
import json
import os
from pathlib import Path
import sys

import numpy as np
from scipy.special import hankel1, jv
from scipy.spatial.distance import cdist

from replay import OPERATORS, SOURCE, ScalarTM, load_operators

LOCATION = Path(__file__).resolve().parent
WORKSPACE = Path(os.environ.get('HIGHDIM_WORKSPACE', LOCATION.parents[2] if len(LOCATION.parents) > 2 else Path.cwd()))
MAX = WORKSPACE / 'research/delegated/max_struct_r1'
SRAP = WORKSPACE / 'research/srap_cie_array_study_20261121'
R8_TREE = MAX / 'r8_worktree'
R8 = R8_TREE / 'research/research_max_struct_r8'
R10 = MAX / 'r10_worktree/research/research_max_struct_r10'
R12 = MAX / 'r12_worktree/research/research_max_struct_r12'
AUS = 'austria_eps4_400MHz'
RING = 'r10_twod_0001'
RING_PARENT = f'{RING}__cie_single__noise_ball_99'
AUS_RESULTS = R8 / 'RESULTS/paper_benchmarks' / AUS
SRAP_CASES = ('s2d_007__pw180_10x20', 's2d_007__pw360_10x20', 's2d_010__pw180_10x20')

USE_FROZEN_INPUTS = (os.environ.get('HIGHDIM_USE_FROZEN') == '1' or
                     not (R8 / 'src/paper_physics.py').is_file())
if USE_FROZEN_INPUTS:
    Model = ScalarTM
else:
    sys.path.insert(0, str(R8_TREE))
    from research.research_max_struct_r8.src.paper_physics import acquisition
    from research.research_z_alg_r4_diversity.src.physics import Model


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def rel(path):
    return str(path.relative_to(WORKSPACE))


@dataclass
class Endpoint:
    m: np.ndarray
    provenance: dict
    stored_prediction: np.ndarray | None = None


@dataclass
class Observation:
    id: str
    side: int
    width: float
    positions: np.ndarray
    bounds: tuple[float, float]
    Y: np.ndarray
    radius: float
    model: Model
    endpoints: dict[str, Endpoint]
    centers: tuple[str, ...]
    provenance: dict


def material_record(path, field, extra, prediction=None):
    with np.load(path, allow_pickle=False) as archive:
        return Endpoint(
            np.array(archive[field], dtype=float).ravel(),
            {'file': rel(path), 'field': field, **extra},
            None if prediction is None else np.array(archive[prediction], dtype=complex),
        )


def pick(d, names):
    return {key: d[key] for key in names if key in d}


def frozen_observation(case):
    document = read_json(SOURCE / f'{case}.json')
    operators = load_operators(case)
    with np.load(SOURCE / f'{case}.npz', allow_pickle=False) as source:
        y, positions = source['Y_archive'], source['positions']
        endpoints = {arm: Endpoint(source[f'history__{arm}__m'], row['source'],
                                   source[f'history__{arm}__Y'])
                     for arm, row in document['historical_baselines'].items()}
    return Observation(case, document['side'], document['domain_width'], positions,
                       tuple(document['material_box']), y, document['archive_radius'], Model(operators),
                       endpoints, tuple(document['centers']), document['acquisition'])


def r10_observation():
    if USE_FROZEN_INPUTS:
        return frozen_observation(RING)
    folder = R10 / 'data_public/twod'
    path = folder / f'{RING}.npz'
    operator_path = folder / 'operators.npz'
    with np.load(path, allow_pickle=False) as obs, np.load(operator_path, allow_pickle=False) as ops:
        y = obs['Y']
        radius = float(obs['sigma'] * obs['noise_radius'])
        bounds = tuple(map(float, obs['material_bounds']))
        operators = {key: ops[key] for key in ('G', 'E0', 'S')}
        positions = ops['centers']
        k = float(ops['k'])
        pulse_radius = float(ops['radius'])
        sources, receivers = ops['sources'], ops['receivers']

    parent = R10 / 'RESULTS/twod' / RING_PARENT
    parent_receipt = read_json(parent / 'RECEIPT.json')
    endpoints = {
        'r10_cie_single': material_record(parent / 'prediction.npz', 'm', {
            'contract': 'historical CN99, uniform .1 initialization, original sourcewise single-beta CIE',
            'receipt': rel(parent / 'RECEIPT.json'),
            'historical_receipt_fields': pick(parent_receipt, (
                'status', 'elapsed', 'elapsed_wall_seconds', 'elapsed_cpu_seconds', 'wall_seconds',
                'cpu_seconds', 'CPU_seconds', 'started_utc', 'ended_utc', 'termination',
            )),
        }),
    }
    for arm in ('direct', 'g19', 'ap', 'dr'):
        folder = R12 / 'supplement_paired/RESULTS' / f'A__{RING}__{arm}'
        receipt, primary = (read_json(folder / name) for name in ('RECEIPT.json', 'primary.json'))
        endpoints[f'r12_{arm}'] = material_record(folder / 'primary.npz', 'm', {
            'contract': 'historical exposed same-start restart from the saved single-beta CIE endpoint',
            'receipt': rel(folder / 'RECEIPT.json'),
            'historical_receipt_fields': pick(receipt, (
                'whole_child_CPU_seconds', 'wall_seconds', 'CPU_cap', 'phase_limits', 'terminations',
            )),
            'historical_endpoint_fields': primary,
        }, 'Y_prediction')

    side = 24
    width = float(np.ptp(positions[:, 0]) * side / (side - 1))
    return Observation(RING, side, width, positions, bounds, y, radius, Model(operators),
                       endpoints, ('r10_cie_single', 'r12_g19'), {
        'handbook_id': ['r10_twod', RING, RING, 'r10_twod'],
        'cross_reference': ['r12_paired_restart', RING, f'A__{RING}', 'restart-A'],
        'exposure': 'MAX-STRUCT R10 material, scored before this study, R12 selected exposed difficult case',
        'generator': 'MAX-STRUCT R10, original project-generated matched-grid scalar TM simulation',
        'measurement': rel(path),
        'operators': rel(operator_path),
        'fine_side': side,
        'k': k,
        'disk_radius': pulse_radius,
        'sources': sources.tolist(),
        'receivers': receivers.tolist(),
        'radius_rule': 'saved CN99 radius times saved complex standard deviation',
    })


def austria_observation():
    if USE_FROZEN_INPUTS:
        return frozen_observation(AUS)
    path = AUS_RESULTS / 'PUBLIC_DATA.npz'
    with np.load(path, allow_pickle=False) as obs:
        y = obs['Y']
        side, fine, frequency = (int(obs[key]) for key in ('inverse_side', 'forward_side', 'frequency'))
    operators = acquisition(side, frequency)
    endpoints = {}
    for arm in ('cie_dual', 'cie_single', 'global_dual', 'physical_lbfgs', 'paper_lsie', 'paper_nie'):
        folder = AUS_RESULTS / arm
        receipt = read_json(folder / 'RECEIPT.json')
        endpoints[arm] = material_record(folder / 'PREDICTION.npz', 'm', {
            'contract': 'historical R8 paper adaptation; LSIE and NIE use zero material prior',
            'receipt': rel(folder / 'RECEIPT.json'),
            'historical_receipt_fields': pick(receipt, (
                'wall_seconds', 'cpu_seconds', 'relative_data_residual', 'sourcewise',
                'status', 'candidate_selection', 'selected_stage', 'corrector',
            )),
        }, 'Y_prediction')

    return Observation(AUS, side, 2., operators['centers'], (0., 20.), y, float(.12 * np.linalg.norm(y)),
                       Model(operators), endpoints, ('cie_dual',), {
        'handbook_id': ['r8_austria', 'Austria/eps4', AUS, 'fine150-originalR8'],
        'exposure': 'Literature geometry, R8 regenerated and scored, reused in R12 and in this study',
        'generator': 'MAX-STRUCT R8 independent implementation of 2016 NIE geometry',
        'measurement': rel(path),
        'receipt': rel(AUS_RESULTS / 'DATA_RECEIPT.json'),
        'operators': rel(R8 / 'src/paper_physics.py') + '::acquisition',
        'fine_side': fine,
        'frequency_Hz': frequency,
        'tx': '20 full-ring unit plane waves',
        'rx': '40 full-ring sites at 3.75 m radius',
        'radius_rule': 'historical R12 .12 times observed norm; exact R8 noise tested after endpoint seal',
        'radius_reference': rel(R12 / 'delegated/external/CONTRACT.json'),
    })


def srap_observation(case):
    if USE_FROZEN_INPUTS:
        return frozen_observation(case)
    path = SRAP / 'data/public' / f'{case}.npz'
    with np.load(path, allow_pickle=False) as obs:
        positions, E, S, y = (obs[key] for key in ('centers', 'E', 'S', 'Y'))
        k, pulse, radius = (float(obs[key]) for key in ('k', 'radius', 'delta'))
        meta = json.loads(obs['case_json'].item())

    # SRAP uses the R8 equal-area disk operator. E and S are the saved arrays.
    distances = cdist(positions, positions)
    np.fill_diagonal(distances, 1.)
    G = .5j * np.pi * k * pulse * jv(1, k * pulse) * hankel1(0, k * distances)
    np.fill_diagonal(G, .5j * np.pi * k * pulse * hankel1(1, k * pulse) - 1.)

    endpoints = {}
    for arm in ('dual', 'srap'):
        folder = SRAP / 'results/predictions' / f'{case}__{arm}.npz'
        with np.load(folder, allow_pickle=False) as data:
            rows = json.loads(data['rows_json'].item())
        for view in ('raw', 'light'):
            receipt = next(row for row in rows if row['budget'] == 30 and row['view'] == view)
            endpoints[f'{arm}_{view}'] = material_record(folder, f'm_30_{view}', {
                'contract': 'historical same-30-second full-physics audited auxiliary or light endpoint',
                'historical_receipt_fields': receipt,
            })

    return Observation(case, meta['side'], meta['width'], positions, (0., 20.), y, radius,
                       Model({'G': G, 'E0': E, 'S': S}), endpoints, ('dual_light',), {
        'handbook_id': ['srap_main24', meta['object_id'], case, meta['array']],
        'exposure': 'Independent material in original SRAP work; scored before this study',
        'generator': 'SRAP confirmation simulation, independently generated on finer grid',
        'measurement': rel(path),
        'operators': rel(SRAP / 'data2d.py') + '::SyntheticOperators',
        'fine_side': meta['fine'],
        'frequency_Hz': meta['frequency'],
        'acquisition': meta,
        'radius_rule': 'saved norm of the actual fixed-norm 10% complex Gaussian sample',
    })


def observations():
    yield r10_observation()
    yield austria_observation()
    for case in SRAP_CASES:
        yield srap_observation(case)


def truth_path(case):
    if case == RING:
        return R10 / 'sealed_truth/twod' / f'{case}.npz'
    if case == AUS:
        return R8 / 'data/paper_benchmarks_private' / f'{case}.npz'
    return SRAP / 'data/truth' / f'{case}.npz'


def freeze_replay_inputs():
    if USE_FROZEN_INPUTS:
        raise RuntimeError('Freezing operators requires the original historical worktrees')
    matrices, seen_greens = {}, {}
    for obs in observations():
        model = obs.model
        key = next((name for name, green in seen_greens.items()
                    if model.G.shape == green.shape and np.array_equal(model.G, green)), None)
        if key is None:
            key = f'G_{len(seen_greens)}'
            matrices[key] = seen_greens[key] = model.G
        matrices[f'{obs.id}__G_key'] = np.array(key)
        matrices[f'{obs.id}__E0'] = model.E0
        matrices[f'{obs.id}__S'] = model.S
    np.savez_compressed(OPERATORS, **matrices)
    for obs in observations():
        stored = load_operators(obs.id)
        if not all(np.array_equal(stored[name], getattr(obs.model, name)) for name in stored):
            raise ArithmeticError(f'Operator packaging mismatch for {obs.id}')
    print(f'Saved {len(seen_greens)} unique exact Green arrays and five sets of incident and receiving matrices.')


if __name__ == '__main__':
    freeze_replay_inputs()
