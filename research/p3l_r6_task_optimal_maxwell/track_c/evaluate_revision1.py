"""Explicit stage-allocation arms against previously sealed full96/64 baseline.

Reuses corrected CDEV statistical machinery in an isolated module namespace.
Never overwrites canonical original evaluation_CDEV artifacts, regenerates its
baseline predictions, constructs a Book/Target, or reads locked-D truths.
"""
from __future__ import annotations
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import argparse
import importlib.util
import json
from pathlib import Path
import sys

TRACK=Path(__file__).resolve().parent
OUT=TRACK.parent
sys.path.insert(0,str(OUT))
from track_c.anchors import REVISION1_CONFIGS


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--configs',nargs='+',required=True,choices=list(REVISION1_CONFIGS))
    parser.add_argument('--output',type=Path,default=TRACK/'evaluation_CDEV_revision1')
    parser.add_argument('--predict-only',action='store_true')
    parser.add_argument('--score-only',action='store_true')
    args=parser.parse_args()
    if len(set(args.configs))!=len(args.configs) or args.predict_only and args.score_only:
        raise ValueError('unique explicit configs and at most one phase flag required')
    if not (TRACK/'REVISION1.json').exists():raise RuntimeError('freeze Revision1 before execution')
    destination=args.output.resolve()
    if not destination.is_relative_to(TRACK.resolve()) or destination==(TRACK/'evaluation_CDEV').resolve():
        raise ValueError('separate Track C revision output required')
    specification=importlib.util.spec_from_file_location('track_c_dev_revision_private',TRACK/'evaluate_dev.py')
    evaluator=importlib.util.module_from_spec(specification);specification.loader.exec_module(evaluator)
    parents=evaluator.expected_parents()
    original_destination=TRACK/'evaluation_CDEV'
    old_seal=json.loads((original_destination/'ALL_PREDICTIONS_SEALED.json').read_text())
    # Validate all original sealed bytes before changing the private module's
    # explicit config family. This reads predictions, never old/new truths.
    old_snapshots=evaluator.validate_seal(original_destination,old_seal,parents)
    reports={}
    for config in (evaluator.BASELINE,*args.configs):
        reports[config]=[]
        for parent in parents:
            folder=TRACK/'contexts/CDEV'/parent/config
            report=json.loads((folder/'REPORT.json').read_text())
            if (report.get('parent')!=parent or report.get('GT_read') is not False or
                    report.get('path_complete') is not True or not (folder/'CONTEXT.npz').exists()):
                raise RuntimeError('NOT_READY: complete64 paired revision contexts required')
            reports[config].append(report)
    evaluator.CONFIG_NAMES=tuple(args.configs)
    if args.score_only:
        revision_seal=json.loads((destination/'ALL_PREDICTIONS_SEALED.json').read_text())
    else:
        revision_seal=evaluator.seal_predictions(destination,parents)
    snapshots=evaluator.validate_seal(destination,revision_seal,parents)
    snapshots.update({key:value for key,value in old_snapshots.items() if key[0]==evaluator.BASELINE})
    global_seal=dict(status='REVISION1_ALL_PREDICTIONS_AND_REUSED_BASELINE_SEALED',truth_read=False,
        CDEV_truth_previously_exposed_to_development=True,new_physics_RHS=0,
        new_baseline_prediction_calls=0,parents=parents,revision_configs=args.configs,
        baseline_original_seal=str(original_destination/'ALL_PREDICTIONS_SEALED.json'),
        baseline_original_artifacts=[row for row in old_seal['artifacts'] if row['config']==evaluator.BASELINE],
        revision_prediction_seal=str(destination/'ALL_PREDICTIONS_SEALED.json'),
        revision_artifacts=revision_seal['artifacts'],model_loads=revision_seal['model_loads'])
    evaluator.write(destination/'REVISION1_GLOBAL_SEAL.json',global_seal)
    if args.predict_only:
        print(json.dumps(dict(status=global_seal['status'],truth_read=False)),flush=True);return
    evaluator.CONFIG_NAMES=(evaluator.BASELINE,*args.configs)
    result=evaluator.score(destination,parents,reports,snapshots,global_seal)
    result.update(revision=1,explicit_configuration_subset=args.configs,
        canonical_original_results_unchanged=True,baseline_predictions_reused_without_recompute=True,
        seal=str(destination/'REVISION1_GLOBAL_SEAL.json'),
        scientific_scope='stage allocation development diagnostic after previous CDEV exposure; no novelty or independent confirmation claim')
    evaluator.write(destination/'CDEV_RESULTS.json',result)
    print(json.dumps(dict(status=result['status'],selection=result['selection'],output=str(destination/'CDEV_RESULTS.json'))),flush=True)


if __name__=='__main__':main()
