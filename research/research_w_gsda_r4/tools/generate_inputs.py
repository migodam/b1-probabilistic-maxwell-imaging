#!/usr/bin/env python3
"""Offline input generation. Preparation is default; physics is explicit."""
from pathlib import Path
import argparse
import importlib
import json
import sys
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from gsda.generate import CASE_IDS, prepare_private, generate_measurements, dump


def load_factory(reference):
    module, sep, name = reference.partition(':')
    if not sep: raise ValueError('Use module:callable for accounted predictor factory')
    return getattr(importlib.import_module(module), name)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=ROOT)
    p.add_argument('--config', type=Path, help='Parent-frozen public physics config')
    p.add_argument('--case', action='append', choices=CASE_IDS, dest='cases')
    p.add_argument('--run-physics', action='store_true', help='Explicit actual PDE generation; parent queue permission required')
    p.add_argument('--predictor-factory', default='gsda.generate:accounted_predictor_factory', help='module:callable accepting (root,config), returning accounted predictor(z,config,phase)')
    a = p.parse_args()
    cases = tuple(a.cases or CASE_IDS)
    if not a.run_physics:
        rows = prepare_private(a.root, cases)
        print(json.dumps(dict(status='PRIVATE_TRUTH_PREPARED_PHYSICS_NOT_RUN', cases=len(rows),
                              full_fe=0, rhs=0, online_measurements_generated=False)))
        return
    if not a.config:
        p.error('--run-physics requires --config; caller owns accounted queue')
    config = json.loads(a.config.read_text())
    predictor = load_factory(a.predictor_factory)(a.root, config)
    try:
        rows = generate_measurements(a.root, config, predictor, cases)
    except BaseException as exc:
        if hasattr(predictor, 'book'):
            predictor.book.failure('GENERATION', 'generate_inputs', exc)
            predictor.book.receipt(status='FAILED', error=repr(exc))
        raise
    if hasattr(predictor, 'book'): predictor.book.receipt(status='COMPLETED', cases=list(cases))
    dump(a.root / 'inputs' / 'GENERATION_RECEIPT.json', dict(status='ACTUAL_2D_MEASUREMENTS_GENERATED',
         cases=rows, truth_private=True, online_npz_keys=['y'], independent_finer_validation=False))
    print(json.dumps(dict(status='ACTUAL_2D_MEASUREMENTS_GENERATED', cases=len(rows))))

if __name__ == '__main__':
    main()
