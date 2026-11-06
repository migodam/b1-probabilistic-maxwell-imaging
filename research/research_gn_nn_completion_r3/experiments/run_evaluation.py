#!/usr/bin/env python3
"""Public finite-candidate evaluation, closure, offline scoring, and plots."""
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np

REPO = Path(__file__).resolve().parents[3]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
STUDY = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads(Path(path).read_text())


def read_npz(path):
    with np.load(path, allow_pickle=False) as z:
        return {k:z[k] for k in z.files}


def evaluate(args):
    batch_started=time.perf_counter();batch_cpu_started=time.process_time();batch_start_unix=time.time()
    from research.research_gn_nn_completion_r3.src.evaluate import evaluate_arm, evaluate_control, evaluate_full45_baseline, save_result
    from research.research_gn_nn_completion_r3.src.score import write_json
    study = Path(args.study).resolve()
    data = Path(args.data).resolve() if args.data else study/'data'
    manifest_path = Path(args.manifest).resolve() if args.manifest else data/'manifests/online.json'
    manifest = read(manifest_path)
    records = [r for r in manifest['records'] if r['split'] in args.splits]
    if args.max_observations is not None:
        records = records[:args.max_observations]
    if args.obs_id:
        requested=set(args.obs_id)
        records=[r for r in records if r['obs_id']in requested]
    if args.partition_count<1 or args.partition_count>2 or not 0<=args.partition_index<args.partition_count:
        raise ValueError('Require one or two deterministic partitions and a valid partition index')
    if args.partition_block_size<1:
        raise ValueError('Positive deterministic partition block size required')
    records=[r for i,r in enumerate(records)if (i//args.partition_block_size)%args.partition_count==args.partition_index]
    operators = read_npz(data/manifest.get('operators_path', 'operators.npz'))
    if args.secondary_mode=='fixed_background_softfilter_box_qp':
        reference=study/'cache/reference_jacobian.npz'
        if reference.exists() and reference.with_suffix('.json').exists():
            cached=read_npz(reference)
            if not np.all(np.asarray(cached['m_ref'])==.1):
                raise ValueError('Fixed-background reference must be public material .1')
            operators['J_ref']=cached['J_ref']
            sidecar=reference.with_suffix('.json')
            operators['reference_setup_receipt']=dict(path=str(reference),receipt=read(sidecar)if sidecar.exists()else None)
    entries = []
    for run in args.runs + args.controls:
        for record in records:
            obs_id = record['obs_id']
            dest = study/'evaluation'/run
            out = dest/f'{obs_id}.json'
            if out.exists() and (dest/f'{obs_id}.npz').exists():
                entries.append(dict(run=run, obs_id=obs_id, status='REUSED_EXISTING'))
                continue
            prediction_path = study/'predictions'/run/f'{obs_id}.npz'
            if run not in args.controls and not prediction_path.exists():
                entries.append(dict(run=run, obs_id=obs_id, status='NOT_RUN_MISSING_PREDICTION'))
                continue
            anchor_path = study/'anchors'/f'{obs_id}.npz'
            if not anchor_path.exists():
                entries.append(dict(run=run, obs_id=obs_id, status='NOT_RUN_MISSING_ANCHOR'))
                continue
            anchor = read_npz(anchor_path)
            public = {**operators, **read_npz(data/record['public_path'])}
            noise = 'clean' if float(record['noise']) == 0 else 'noisy'
            tick = time.perf_counter()
            try:
                if run == 'weakreg45_anchor':
                    weak_path=study/'anchors'/f'{obs_id}_weakreg1e-5.npz'
                    if not weak_path.exists() or not weak_path.with_suffix('.json').exists():
                        entries.append(dict(run=run,obs_id=obs_id,status='NOT_RUN_MISSING_WEAKREG45_ANCHOR'))
                        continue
                    receipt,arrays=evaluate_full45_baseline(public,anchor,read_npz(weak_path),read(weak_path.with_suffix('.json')),noise=noise)
                elif run in args.controls:
                    receipt, arrays = evaluate_control(public, anchor, control=run, noise=noise)
                else:
                    prediction = read_npz(prediction_path)
                    sidecar = prediction_path.with_suffix('.json')
                    info = read(sidecar) if sidecar.exists() else {}
                    receipt, arrays = evaluate_arm(public, anchor, prediction, arm=info.get('arm',run), noise=noise,
                        preparation_cost=info.get('preparation_cost'), secondary_mode=args.secondary_mode)
            except Exception as exc:
                # Preserve malformed/failed inputs as explicit failed attempts,
                # without fabricating a material endpoint or exact missing costs.
                receipt=dict(schema='gn-nn-completion-r3-evaluation-v1',arm=run,noise=noise,
                             status='FAILED_INPUT_OR_SETUP',stop_reason=f'{type(exc).__name__}: {exc}',
                             selection_status='FAIL_NO_FEASIBLE_CANDIDATE',accounting_status='HOLD_INCOMPLETE_FAILURE_COSTS',
                             wall_seconds=time.perf_counter()-tick,costcounts={},candidates=[],trace=[],attempts=[],
                             top1=None,old_score_top1=None,truth_access='NONE_PUBLIC_INPUT_ONLY',
                             failure_cost_scope='elapsed public call is known; missing preparation/setup operation counts remain unknown')
                arrays=dict(anchor=np.asarray(anchor['m0']).reshape(576),candidate_materials=np.empty((0,576)))
            receipt.update(obs_id=obs_id, parent_id=record['parent_id'], split=record['split'],
                           run=run, prediction_path=str(prediction_path)if run not in args.controls else None,
                           anchor_path=str(study/'anchors'/f'{obs_id}.npz'),
                           outer_call_seconds=time.perf_counter()-tick,
                           execution_partition=dict(index=args.partition_index,count=args.partition_count,block_size=args.partition_block_size,rule='(selected_manifest_index // block_size) % partition_count',BLAS_threads=1))
            save_result(dest, obs_id, receipt, arrays)
            entries.append(dict(run=run,obs_id=obs_id,status=receipt['status'],selection_status=receipt['selection_status'],wall_seconds=receipt['wall_seconds']))
            print(json.dumps(entries[-1]),flush=True)
    partition_suffix=f'_P{args.partition_index}OF{args.partition_count}'if args.partition_count>1 else''
    path=study/'evaluation'/f'RUN_MANIFEST_{"_".join(args.splits)}{partition_suffix}_T{time.time_ns()}.json'
    write_json(path,dict(manifest_path=str(manifest_path),runs=args.runs,controls=args.controls,splits=args.splits,
                         max_observations=args.max_observations,partition_index=args.partition_index,partition_count=args.partition_count,partition_block_size=args.partition_block_size,partition_rule='(selected_manifest_index // block_size) % partition_count',BLAS_threads=1,batch_start_unix=batch_start_unix,batch_end_unix=time.time(),batch_wall_seconds=time.perf_counter()-batch_started,batch_cpu_seconds=time.process_time()-batch_cpu_started,batch_cost_scope='full partition process including public input/output IO; excludes already invoiced NN inference and shared GN construction',entries=entries,truth_access='NONE'))
    return path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['evaluate','seal','score','plots'])
    parser.add_argument('--study',default=str(STUDY))
    parser.add_argument('--data')
    parser.add_argument('--manifest')
    parser.add_argument('--runs',nargs='*',default=[])
    parser.add_argument('--controls',nargs='*',choices=['continue_trf','weakreg_trf','weakreg45_anchor'],default=[])
    parser.add_argument('--splits',nargs='+',default=['val','test'])
    parser.add_argument('--max-observations',type=int)
    parser.add_argument('--obs-id',nargs='+')
    parser.add_argument('--partition-count',type=int,default=1,choices=[1,2])
    parser.add_argument('--partition-index',type=int,default=0,choices=[0,1])
    parser.add_argument('--partition-block-size',type=int,default=8)
    parser.add_argument('--secondary-mode',choices=['full_softfilter_box_qp','fixed_background_softfilter_box_qp'])
    parser.add_argument('--seal')
    parser.add_argument('--output')
    parser.add_argument('--oracle',action='store_true')
    parser.add_argument('--replay-frozen',action='store_true',help='Already-disclosed historical replay after relocation; do not assert original blind-seal mtimes')
    args=parser.parse_args()
    study=Path(args.study).resolve()
    if args.command=='evaluate':
        path=evaluate(args)
    elif args.command=='seal':
        from research.research_gn_nn_completion_r3.src.score import seal_predictions
        data=Path(args.data).resolve()if args.data else study/'data'
        manifest=Path(args.manifest).resolve()if args.manifest else data/'manifests/evaluation.json'
        path=seal_predictions(study,manifest,args.runs+args.controls,splits=args.splits,output=args.output,max_observations=args.max_observations)
    elif args.command=='score':
        from research.research_gn_nn_completion_r3.src.score import score_predictions
        path=score_predictions(args.seal or study/'evaluation/PREDICTIONS_SEALED.json',output=args.output,include_oracle=args.oracle,study_override=study,replay_frozen=args.replay_frozen)
    else:
        from research.research_gn_nn_completion_r3.src.plots import make_plots
        path=make_plots(args.output or study/'evaluation/SCORES.json',study=study,replay_frozen=args.replay_frozen)
    print(path)


if __name__=='__main__':
    main()
