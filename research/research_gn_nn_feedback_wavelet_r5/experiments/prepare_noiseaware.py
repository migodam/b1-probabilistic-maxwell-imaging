"""Phase II public-only common all-scales discrepancy-aware GN preparation.

No changed forward/optimizer, truth loading, reference recomputation, or training.
Finite two-worker queue stops submitting at 19 minutes; legitimate in-flight
45s solves finish. All remaining planned inputs are explicitly NOT_RUN.
"""
from __future__ import annotations
import os
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_name] = '1'
import argparse
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import json
import multiprocessing
from pathlib import Path
import shutil
import sys
import time
import numpy as np

STUDY = Path(__file__).resolve().parents[1]
REPO = STUDY.parents[1]
R4 = STUDY.parent/'research_gn_nn_preservation_r4'
sys.path.insert(0, str(REPO))
from research.research_gn_nn_preservation_r4.src.anchors import DEFAULTS, generate_anchor
ARMS = ['pixel_closed_coupled', 'haar_closed_coupled', 'db2_closed_coupled',
        'haar_closed_no_fullfield', 'haar_free_coupled']


def read(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def task(record, output, config):
    """Public-only worker; record contains no truth/support/family information."""
    obs_id = record['obs_id']
    path = Path(output)/f'{obs_id}.npz'
    started = time.perf_counter()
    if path.exists() or path.with_suffix('.json').exists():
        raise FileExistsError('Refusing to overwrite prior new TRAIN work: '+str(path))
    try:
        result = generate_anchor(R4/'data/operators.npz', R4/'data'/record['public_path'],
                                 path, R4/'data/reference_jacobian.npz', weakreg=0., config=config)
        return dict(obs_id=obs_id, parent_id=record['parent_id'], split='train',
                    disposition='GENERATED', status=result['status'],
                    stop_reason=result['stop_reason'], actual_task_wall_seconds=time.perf_counter()-started,
                    solve_wall_seconds=result['solve_wall_seconds'],
                    deployment_anchor_wall_seconds=result['anchor_plus_feature_wall_seconds'],
                    actual_new_physics_counts=result['counts'],
                    original_deployment_counts=result['counts'], output=str(path),
                    truth_access='NONE_PUBLIC_Y_AND_NOISE_ONLY')
    except Exception as exc:
        failure = dict(obs_id=obs_id, parent_id=record['parent_id'], split='train',
                       disposition='FAILED', status='FAILED',
                       error=f'{type(exc).__name__}: {exc}',
                       actual_task_wall_seconds=time.perf_counter()-started,
                       actual_new_physics_counts=None,
                       cost_uncertainty='Uncaught generator error does not expose partial internal counters; elapsed time preserved, counts unknown rather than zero',
                       truth_access='NONE_PUBLIC_Y_AND_NOISE_ONLY')
        # If a late serialization failure left an endpoint, retain it.
        save(path.with_suffix('.failure.json'), failure)
        return failure


def prepare(study):
    start = time.perf_counter()
    output = study/'inputs/noiseaware_anchors'
    if output.exists() and any(output.iterdir()):
        raise FileExistsError('One-run preparation refuses a nonempty output directory')
    output.mkdir(parents=True, exist_ok=True)
    cohort = read(study/'audit/COHORT.json')
    train_ids, val_ids = set(cohort['train_parent_ids']), set(cohort['val_parent_ids'])
    if len(train_ids) != 128 or len(val_ids) != 78:
        raise ValueError('Frozen Phase I cohort must contain128TRAIN/78VAL mothers')
    metadata = read(R4/'data/manifests/evaluation.json')['records']
    train = [r for r in metadata if r['split'] == 'train' and r['parent_id'] in train_ids]
    val = [r for r in metadata if r['split'] == 'val' and r['parent_id'] in val_ids]
    if len(train) != 256 or len(val) != 156:
        raise ValueError('Need original256TRAIN and156VAL observations')
    reference = R4/'data/reference_jacobian.npz'
    with np.load(reference, allow_pickle=False) as bundle:
        reference_shapes = {key: list(bundle[key].shape) for key in bundle.files}
    if reference_shapes.get('J_ref') != [1152, 576]:
        raise ValueError('Existing public reference J shape is incompatible')
    config = dict(DEFAULTS, discrepancy_each_scale=True)
    if config['discrepancy_stop'] or config['stages'] != [8, 12, 24]:
        raise ValueError('No final-only stopping flag or changed scale schedule allowed')
    source_configs = []
    for r in val:
        for extension in ('npz', 'json'):
            if not (R4/'anchors'/f"{r['obs_id']}_discrepancy_allscales.{extension}").exists():
                raise FileNotFoundError('Missing frozen VAL allscales anchor '+r['obs_id'])
        source = read(R4/'anchors'/f"{r['obs_id']}_discrepancy_allscales.json")
        if not source['config'].get('discrepancy_each_scale'):
            raise ValueError('Cached VAL source is not allscales discrepancy-aware')
        source_configs.append(source['config'])
    freeze = dict(schema='r5-noiseaware-phase-II-preparation-v1', created_unix=time.time(),
        status='FROZEN_BEFORE_GENERATION_OR_COPY', sole_scientific_change='Common GN anchor stopping becomes noise-aware at every8/12/24scale for all five arms',
        scope='Known DEVELOPMENT TRAIN/VAL only; no TEST, no blind claim',
        all_five_arms=ARMS, training='All five from scratch; unchanged network/loss/depth/training budget; parent owns trainer',
        train_mothers=128, train_observations=256, val_mothers=78, val_observations=156,
        train_parent_ids=sorted(train_ids), val_parent_ids=sorted(val_ids),
        train_observation_ids=[r['obs_id'] for r in train], val_observation_ids=[r['obs_id'] for r in val],
        config=config, generator='Unchanged R4 src.anchors.generate_anchor/Model; weakreg=0',
        public_inputs='Original observations contain only Y,noise; immutable original G/E/S operators',
        reference=dict(path=str(reference), keys_and_shapes=reference_shapes,
                       setup='Existing original reference reused; no new SVD/Jacobian/reference solve',
                       original_receipt=str(reference.with_suffix('.json'))),
        cached_VAL_source_config=source_configs[0],
        cached_VAL_semantics='Original cache has each_scale=True and discrepancy_stop=True; each_scale enables every scale. New TRAIN each_scale=True and inherited final-onlyflagFalse. Same effective allscales stopping, no frozen source edits',
        cached_VAL='Copy original NPZ and JSON to new plain obs filenames; original deployment costs and newly executed copy costs kept separately',
        workers=2, BLAS_threads=1, device='local Apple Silicon CPU',
        finite_queue='At most two active tasks; stop new submissions after19minutes, finish legitimate in-flight45s solve, NOT_RUN for remaining',
        truth_arrays_read=False, no_hyperparameter_sweep=True)
    freeze_path = study/'audit/NOISEAWARE_PREPARATION_FREEZE.json'
    if freeze_path.exists():
        raise FileExistsError('Phase II preparation already frozen; no unrecorded rerun')
    save(freeze_path, freeze)
    ledger = []
    for record in val:
        tick = time.perf_counter()
        source_base = R4/'anchors'/f"{record['obs_id']}_discrepancy_allscales"
        for extension in ('npz', 'json'):
            shutil.copy2(source_base.with_suffix('.'+extension), output/f"{record['obs_id']}.{extension}")
        seconds = time.perf_counter()-tick
        source = read(source_base.with_suffix('.json'))
        ledger.append(dict(obs_id=record['obs_id'], parent_id=record['parent_id'], split='val',
            disposition='COPIED_FROZEN_CACHE', status=source['status'],
            cached_source_npz=str(source_base.with_suffix('.npz')), cached_source_json=str(source_base.with_suffix('.json')),
            actual_new_copy_wall_seconds=seconds, actual_new_physics_counts={'LU': 0, 'actual_complex_RHS': 0},
            original_deployment_anchor_wall_seconds=source['anchor_plus_feature_wall_seconds'],
            original_deployment_counts=source['counts'], original_feature_counts=source.get('feature_counts'),
            receipt_semantics='Copied JSON preserves original source paths/costs; this ledger alone distinguishes new copy operation',
            output=str(output/f"{record['obs_id']}.npz"), truth_access='NONE'))
    def checkpoint(status):
        done = [r for r in ledger if r['split'] == 'train']
        save(study/'audit/NOISEAWARE_PREPARATION_PROGRESS.json', dict(status=status,
             elapsed_wall_seconds=time.perf_counter()-start, train_finished=len(done), train_planned=256,
             val_copied=156, rows=ledger, truth_arrays_read=False))
    checkpoint('VAL_COPIED_TRAIN_STARTING')
    safe_records = [{k: r[k] for k in ('obs_id', 'parent_id', 'public_path')} for r in train]
    iterator = iter(safe_records)
    unsubmitted = []
    with ProcessPoolExecutor(max_workers=2, mp_context=multiprocessing.get_context('spawn')) as pool:
        active = {}
        for _ in range(2):
            record = next(iterator, None)
            if record:
                active[pool.submit(task, record, str(output), config)] = record
        while active:
            completed, _ = wait(active, timeout=30., return_when=FIRST_COMPLETED)
            for future in completed:
                record = active.pop(future)
                try:
                    row = future.result()
                except Exception as exc:
                    row = dict(obs_id=record['obs_id'], parent_id=record['parent_id'], split='train',
                               disposition='FAILED', status='FAILED', error=f'{type(exc).__name__}: {exc}',
                               actual_new_physics_counts=None, cost_uncertainty='Worker exception; counters unavailable')
                ledger.append(row)
                if time.perf_counter()-start < 19*60:
                    following = next(iterator, None)
                    if following:
                        active[pool.submit(task, following, str(output), config)] = following
            checkpoint('TRAIN_RUNNING')
    unsubmitted.extend(iterator)
    for record in unsubmitted:
        ledger.append(dict(obs_id=record['obs_id'], parent_id=record['parent_id'], split='train',
                           disposition='NOT_RUN', status='NOT_RUN',
                           reason='19minute submission deadline; legitimate per-anchor budgets unchanged',
                           actual_new_physics_counts={'LU': 0, 'actual_complex_RHS': 0}))
    from collections import Counter
    generated = [r for r in ledger if r['disposition'] == 'GENERATED']
    failures = [r for r in ledger if r['disposition'] == 'FAILED']
    notrun = [r for r in ledger if r['disposition'] == 'NOT_RUN']
    counts = {key: sum(r['actual_new_physics_counts'].get(key, 0) for r in generated) for key in ('LU', 'actual_complex_RHS')}
    result = dict(status='COMPLETE' if len(generated) == 256 else 'PARTIAL_RETAINED',
        actual_elapsed_wall_seconds=time.perf_counter()-start, train_generated=len(generated),
        train_failed=len(failures), train_not_run=len(notrun), val_cached_copied=156,
        actual_new_physics_counts=counts, actual_train_task_wall_seconds_sum=sum(r.get('actual_task_wall_seconds', 0.) for r in generated+failures),
        actual_copy_wall_seconds_sum=sum(r.get('actual_new_copy_wall_seconds', 0.) for r in ledger),
        status_counts=dict(Counter(r['status'] for r in ledger)),
        failure_physics_counts_unknown=bool(failures), rows=ledger,
        cached_VAL_original_cost_not_zeroed=True, no_new_reference_setup=True,
        scope='Public anchor preparation only; original cohort and frozen science budgets unchanged; no GT/TEST/training access')
    save(study/'audit/NOISEAWARE_PREPARATION_RECEIPT.json', result)
    checkpoint(result['status'])
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, default=STUDY)
    args = parser.parse_args()
    prepare(args.study.resolve())


if __name__ == '__main__':
    main()
