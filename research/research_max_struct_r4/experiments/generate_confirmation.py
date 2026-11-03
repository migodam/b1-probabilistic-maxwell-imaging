"""Write-once six-object confirmation generator. Never import into inversion.

Requires a parent-finalized configs/confirmation_frozen.json and explicit flag.
This module reuses R3 shape_alpha without changing its geometry distribution.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import json
import secrets
import shutil
import sys
import time
import traceback
import uuid
import numpy as np
from scipy.stats import chi2
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
sys.path.insert(0, str(REPO))
from research.research_max_struct_r3.experiments.generate_blind import shape_alpha
from research.research_z_alg_r4_diversity.src.physics import Model
from research.research_max_struct_r4.experiments.package_confirmation import SOURCE_FILES, REMOTE_SOURCE


def utc():
    return datetime.now(timezone.utc).isoformat()


def dump_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)


def generate(config_path):
    wall, cpu = time.perf_counter(), time.process_time()
    cfg = json.loads(config_path.read_text())
    if not (cfg['status'] == 'FROZEN_BEFORE_GENERATION' and
            cfg['algorithm_finalized'] is True and cfg['generation_authorized'] is True):
        raise RuntimeError('Parent-finalized algorithm and generation freeze required')
    families = list(cfg['private_shape_families'])
    if cfg['case_count'] != 6 or sorted(families) != sorted(['ring', 'curved_bar', 'two_islands']*2):
        raise ValueError('Exactly two independent objects per frozen family required')
    public = ROOT/'data_public/confirmation'
    private = ROOT/'sealed_truth/confirmation'
    audit = ROOT/'audit/confirmation'
    if any(path.exists() for path in (public, private, audit)):
        raise FileExistsError('Write-once realization/attempt must be preserved')
    public.mkdir(parents=True)
    private.mkdir(parents=True)
    audit.mkdir(parents=True)
    receipt = dict(status='RUNNING', started_utc=utc(), observations=0,
                   truth_access='GENERATOR_ONLY', source='experiments/generate_confirmation.py',
                   frozen_config=str(config_path.relative_to(ROOT)), operator_source=cfg['operator_source'],
                   operator_assembly='Frozen existing G/E0/S loaded; no reassembly',
                   generator_import_time='Outside this generate timer; process launcher accounts startup',
                   failed_attempts=[], records=[])
    model = None
    try:
        dump_new(audit/'FROZEN_CONFIG.json', cfg)
        # Capture frozen implementation BEFORE any new realization/forward call.
        snapshot = audit/'executed_source_snapshot'
        for name in SOURCE_FILES+[
                'research/research_max_struct_r4/experiments/generate_confirmation.py',
                'research/research_max_struct_r3/experiments/generate_blind.py']:
            target = snapshot/name
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(REPO/name,target)
        for name in ('remote_runner.py','vendor/threadpoolctl.py'):
            target = snapshot/name
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(REMOTE_SOURCE/name,target)
        with np.load(REPO/cfg['operator_source'], allow_pickle=False) as archive:
            ops = {k: archive[k] for k in ('G', 'E0', 'S', 'centers', 'sources', 'receivers', 'k', 'radius')}
        if (ops['G'].shape, ops['E0'].shape, ops['S'].shape) != ((576, 576), (576, 12), (48, 576)):
            raise ValueError('Frozen acquisition mismatch')
        model = Model(ops)
        operator_payload = {k: ops[k] for k in ('G', 'E0', 'S', 'sources', 'receivers', 'k', 'radius')}
        operator_payload.update(coords=ops['centers'], shape=np.array([24, 24]), grid_shape=np.array([24, 24]))
        np.savez_compressed(public/'operators.npz', **operator_payload)
        # Randomize family order independently of each per-object geometry/noise seed.
        secrets.SystemRandom().shuffle(families)
        pending, mapping, seen = [], [], set()
        with threadpool_limits(limits=1):
            for family in families:
                seed = secrets.randbits(128)
                rng = np.random.default_rng(seed)
                ident = 'obs_'+uuid.uuid4().hex[:16]
                if ident in seen:
                    raise RuntimeError('Opaque ID collision; preserve attempt without regeneration')
                seen.add(ident)
                alpha = shape_alpha(family, rng)
                m = cfg['material']['host']+cfg['material']['amplitude']*alpha.ravel()
                tick, ct = time.perf_counter(), time.process_time()
                state = model.forward(m, label='fresh_confirmation_generation')
                clean = state['y']
                sigma = cfg['noise']['relative_complex_amplitude']*np.linalg.norm(clean)/np.sqrt(clean.size)
                noise = sigma/np.sqrt(2)*(rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape))
                radius = float(np.sqrt(chi2.ppf(cfg['noise']['radius_confidence'], 2*clean.size)/2))
                np.savez_compressed(private/(ident+'.npz'), m=m, alpha=alpha, mask=alpha >= .5,
                                    Y_clean=clean, noise=noise)
                payload = dict(Y=clean+noise, sigma=np.array(sigma), noise_radius=np.array(radius),
                               host=np.array(cfg['material']['host']), amplitude=np.array(cfg['material']['amplitude']),
                               material_bounds=np.array(cfg['material']['box']))
                payload.update({k: np.array(v) for k, v in cfg['prior'].items()})
                pending.append((ident, payload))
                mapping.append(dict(observation_id=ident, family=family, seed_PRIVATE=str(seed),
                                    truth_file=ident+'.npz'))
                receipt['records'].append(dict(observation_id=ident, physical_LU=1, source_complex_RHS=12,
                    linear_relative_residual=float(state['state_relative_residual']),
                    wall_seconds=time.perf_counter()-tick, cpu_seconds=time.process_time()-ct))
        for ident, payload in sorted(pending):
            np.savez_compressed(public/(ident+'.npz'), **payload)
        dump_new(private/'mapping_PRIVATE.json', dict(records=mapping,
            policy='Hold until all18 frozen prediction dispositions sealed and parent authorizes scoring'))
        dump_new(public/'manifest.json', dict(status='FRESH_UNSCORED', created_utc=utc(), count=6,
            operators_path='operators.npz', records=[dict(observation_id=ident, path=ident+'.npz') for ident, _ in sorted(pending)],
            operator_keys=sorted(operator_payload), observation_keys=sorted(pending[0][1]),
            truth_access='NONE_TO_INVERSION', family_mapping='PRIVATE', frozen_config='configs/confirmation_frozen.json'))
        receipt.update(status='GENERATED_PRIVATE_TRUTH_HELD', observations=6)
    except Exception as exc:
        receipt.update(status='FAILED', error=repr(exc), traceback=traceback.format_exc())
        receipt['failed_attempts'].append(dict(error=repr(exc), charged=True))
        raise
    finally:
        receipt.update(finished_utc=utc(), physical_counts=dict(model.counts) if model else {},
                       wall_seconds=time.perf_counter()-wall, cpu_seconds=time.process_time()-cpu,
                       timer_scope='operator load/model setup, six forward states, noise, private/public serialization; receipt own write excluded')
        dump_new(audit/'GENERATION_RECEIPT.json', receipt)
    print(json.dumps(dict(status=receipt['status'], observations=6, public_manifest=str(public/'manifest.json'),
                          truth_exposed_to_inversion=False)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--execute-after-freeze', action='store_true')
    args = parser.parse_args()
    if not args.execute_after_freeze:
        parser.error('Explicit parent generation authorization required')
    generate(ROOT/'configs/confirmation_frozen.json')
