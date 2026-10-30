"""Offline generation only: private shapes never enter online payloads.

All arrays use ascending physical y. Public R3 glyph rows are top-to-bottom
and are flipped once on construction. Importing this module does no physics.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

TX = np.linspace(-.7, .7, 8)
TY = np.linspace(1.4, 2.8, 8)
PATTERN_ROWS = {
    1: ('00100','01100','00100','00100','00100','00100','01110'),
    2: ('01110','10001','00001','00010','00100','01000','11111'),
    3: ('11110','00001','00001','01110','00001','00001','11110'),
    4: ('00010','00110','01010','10010','11111','00010','00010'),
    5: ('11111','10000','10000','11110','00001','00001','11110'),
    6: ('01110','10000','10000','11110','10001','10001','01110'),
    7: ('11111','00001','00010','00100','01000','01000','01000'),
    8: ('11111','10001','10001','11111','10001','10001','11111'),
    9: ('01110','10001','10001','01111','00001','00001','01110'),
}
# This ordering is private. The online IDs contain no geometry or digit labels.
CASE_IDS = tuple(f'c{i:03d}' for i in range(13))
DIGIT_BY_CASE = dict(zip(CASE_IDS[1:10], (4, 9, 2, 7, 1, 8, 5, 3, 6)))


def dump(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    def convert(x):
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.generic): return x.item()
        if isinstance(x, Path): return str(x)
        raise TypeError(type(x).__name__)
    path.write_text(json.dumps(value, indent=2, default=convert, allow_nan=False) + '\n')


def digit_mask(digit):
    m = np.zeros((7, 7), dtype=bool)
    m[:, 1:6] = np.asarray([[c == '1' for c in row] for row in PATTERN_ROWS[digit]])
    return np.flipud(m).copy()


def case_shape(case_id):
    """Return private physical-grid support and evaluator-only identity."""
    if case_id in DIGIT_BY_CASE:
        d = DIGIT_BY_CASE[case_id]
        return digit_mask(d), dict(cohort='digit', digit=d, shape='public_R3_strokes')
    yy, xx = np.mgrid[-3:4, -3:4]
    if case_id == 'c000':
        # Independent asymmetric stair-and-island development object; no digit.
        m = np.zeros((7, 7), bool)
        m[1:5, 1] = True; m[1, 1:4] = True
        m[3:5, 4:6] = True; m[5, 5] = True
        return m, dict(cohort='development', digit=None, shape='asymmetric_stair_island')
    if case_id == 'c010':
        u = .866025403784 * xx + .5 * yy
        v = -.5 * xx + .866025403784 * yy
        m = (u / 2.6) ** 2 + (v / 1.8) ** 2 <= 1
        return m, dict(cohort='mechanism_control', digit=None, shape='ellipse')
    if case_id == 'c011':
        m = np.zeros((7, 7), bool)
        m[1:6, 1:3] = True; m[1:3, 3:6] = True; m[4:6, 3:5] = True
        return m, dict(cohort='mechanism_control', digit=None, shape='nonconvex')
    if case_id == 'c012':
        m = (xx / 2.8) ** 2 + (yy / 2.6) ** 2 <= 1
        m[(xx / 1.2) ** 2 + (yy / 1.1) ** 2 <= 1] = False
        return m, dict(cohort='mechanism_control', digit=None, shape='holed')
    raise ValueError(f'Unknown private case {case_id}')


def make_truth(case_id):
    """Deterministic independent worlds; all 110 online coordinates remain free."""
    support, identity = case_shape(case_id)
    i = CASE_IDS.index(case_id)
    rng = np.random.default_rng(43107000 + i)
    eps = rng.uniform(2, 4.5, 6)
    wall_loss = rng.uniform(.02, .12, 6)
    level = 7. if case_id == 'c012' else 3.5
    contrast = np.where(support, level, 0.)
    loss = np.where(support, .06, 0.)
    if case_id == 'c011':
        eps[1] = 1.; wall_loss[1] = 0.
    z = np.r_[eps, wall_loss, contrast.ravel(), loss.ravel()]
    return dict(z=z, mask=support, target_contrast=contrast, target_loss=loss,
                target_x_edges=TX, target_y_edges=TY), identity


def public_payload(config, case_id):
    """Strict generic metadata allowlist; no seed, support, case role or label."""
    physics_keys = ('frequencies', 'heldout_frequencies', 'sigma', 'sigma_real_imag',
                    'wall_x', 'wall_y', 'wall_x_edges', 'wall_y_edges', 'target_x_edges', 'target_y_edges',
                    'h', 'targetstep', 'wall_refine', 'order', 'max_factor_bytes', 'tx', 'rx',
                    'bounds', 'material_box', 'parameter_scale', 'model', 'TX_count', 'RX_count',
                    'source_indices', 'train_rx_indices', 'heldout_rx_indices', 'noise_order',
                    'primary_contract', 'target_degree', 'shape', 'version')
    payload = {k: config[k] for k in physics_keys if k in config}
    payload.update(case_id=case_id, schema='R4_PUBLIC_INPUT_V1', target_shape=[7, 7],
                   parameter_order=['wall_epsilon6', 'wall_loss6', 'target_contrast49', 'target_loss49'],
                   target_x_edges=TX.tolist(), target_y_edges=TY.tolist(),
                   support_threshold=1.0, target_all_cells_free=True)
    return payload


def prepare_private(root, case_ids=CASE_IDS):
    """Prepare truth only. Explicit physics phase separately creates online data."""
    root = Path(root)
    records = []
    for case_id in case_ids:
        out = root / 'inputs' / 'private' / case_id
        if (out / 'truth.npz').exists():
            raise FileExistsError(f'Refusing truth regeneration: {out}')
        truth, identity = make_truth(case_id)
        out.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(out / 'truth.npz', **truth)
        dump(out / 'identity.json', dict(case_id=case_id, **identity))
        records.append(dict(case_id=case_id, **identity, foreground_cells=int(truth['mask'].sum())))
    dump(root / 'inputs' / 'private' / 'CASE_MAP.json', dict(cases=records, truth_private=True))
    return records


def generate_measurements(root, config, predictor, case_ids=CASE_IDS):
    """Run only when orchestrator supplies its accounted predictor(z, config, phase).

    Predictor returns Prediction(data, diagnostics) or (data, diagnostics). The
    parent owns the physics queue, clocks, budget reservation and implementation.
    """
    root = Path(root)
    sigma = float(config.get('sigma_real_imag', config.get('sigma', .01)))
    if sigma <= 0: raise ValueError('sigma must be positive')
    records = []
    for case_id in case_ids:
        private = root / 'inputs' / 'private' / case_id
        online = root / 'inputs' / 'online' / case_id
        if (online / 'observation.npz').exists():
            raise FileExistsError(f'Refusing observation regeneration: {online}')
        with np.load(private / 'truth.npz', allow_pickle=False) as a:
            truth = {k: a[k] for k in a.files}
        pred = predictor(truth['z'], config, f'GENERATE_{case_id}')
        clean, diagnostics = (pred.data, pred.diagnostics) if hasattr(pred, 'data') else pred
        clean = np.asarray(clean)
        if not np.iscomplexobj(clean) or not np.all(np.isfinite(clean)):
            raise ValueError('Generator requires finite actual complex predictions')
        seed = 43108000 + CASE_IDS.index(case_id)
        rng = np.random.default_rng(seed)
        y = clean + sigma * (rng.normal(size=clean.shape) + 1j * rng.normal(size=clean.shape))
        truth['clean'] = clean
        heldout_diag = None
        if config.get('heldout_frequencies'):
            hold_config = dict(config, frequencies=config['heldout_frequencies'])
            hp = predictor(truth['z'], hold_config, f'GENERATE_HELDOUT_{case_id}')
            heldout, heldout_diag = (hp.data, hp.diagnostics) if hasattr(hp, 'data') else hp
            heldout = np.asarray(heldout)
            if not np.iscomplexobj(heldout) or not np.all(np.isfinite(heldout)):
                raise ValueError('Heldout truth requires finite actual complex predictions')
            truth['heldout_clean'] = heldout
        # Private truth is extended before publishing the measurement payload.
        np.savez_compressed(private / 'truth.npz', **truth)
        online.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(online / 'observation.npz', y=y)
        dump(online / 'config.json', public_payload(config, case_id))
        dump(private / 'generation.json', dict(status='ACTUAL_2D_MEASUREMENTS_GENERATED',
             case_id=case_id, diagnostics=diagnostics, heldout_diagnostics=heldout_diag,
             observation_shape=list(y.shape), online_npz_keys=['y'], same_discrete=True))
        records.append(dict(case_id=case_id, status='ACTUAL_2D_MEASUREMENTS_GENERATED'))
        print(json.dumps(records[-1]), flush=True)
    return records


def make_vectors(case_ids=CASE_IDS):
    """Pure offline vectors; returns private truth/identity without I/O or PDE."""
    return {case_id: make_truth(case_id) for case_id in case_ids}


def accounted_predictor_factory(root, config):
    """Parent-accounted FEM adapter, imported only in explicit physics mode."""
    from gsda.accounting import BudgetBook
    from gsda.fem_model import FEMModel
    book = BudgetBook(root=root, job='generate_inputs')
    cache = {}
    def predictor(z, public_config, phase):
        key = tuple(public_config['frequencies'])
        if key not in cache:
            cache[key] = FEMModel(public_config, book, phase=phase)
        return cache[key].predict(z)
    predictor.book = book
    return predictor
