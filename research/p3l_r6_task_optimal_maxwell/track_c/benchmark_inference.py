"""Clean CPU1,30-repeat matched online-head timing; never a physical benchmark.

Starts only when central physics and matched training have no active jobs, and
checks between repeats. Warm loaded models predict every CDEV parent singly in
the same shuffled parent order; arm order is shuffled to balance timing drift.
No truth is opened and no prediction is used for scientific risk selection.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'
import argparse,csv,json,sys,time
from pathlib import Path
import numpy as np
TRACK=Path(__file__).resolve().parent;OUT=TRACK.parent;sys.path.insert(0,str(OUT))
from track_c.frozen_heads import load_action,predict_contexts,SEEDS
from track_c.evaluate_dev import expected_parents,load_contexts,write

ARMS={
    'original_full_deterministic':('native96_probit64','deterministic_matched','actual'),
    'ordinary_low48P64_deterministic':('native48_probit64','deterministic_matched','actual'),
    'ordinary_low24P64_deterministic':('native24_probit64','deterministic_matched','actual'),
    'qKKT_low24P64_deterministic':('native24_probit64','deterministic_matched','q'),
    'original_full_distilled':('native96_probit64','distilled_w128_n32768','actual'),
}


def require_idle():
    path=OUT/'RESOURCE_LEDGER_SNAPSHOT.json'
    if path.exists():
        ledger=json.loads(path.read_text())
        if ledger.get('active_jobs') or ledger.get('active_actions'):raise RuntimeError('HOLD: central physics/job active during clean benchmark')
    training=TRACK/'MATCHED_TRAINING_LEDGER.jsonl'
    rows=[json.loads(line) for line in training.read_text().splitlines()] if training.exists() else []
    live={row['job'] for row in rows if row['event']=='BEGIN'}-{row['job'] for row in rows if row['event']=='END'}
    if live:raise RuntimeError('HOLD: matched MPS/CPU training active during clean benchmark')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=TRACK/'clean_inference_benchmark30')
    parser.add_argument('--repetitions',type=int,default=30)
    parser.add_argument('--arms',nargs='+',choices=list(ARMS),default=list(ARMS))
    args=parser.parse_args()
    if args.repetitions!=30:raise ValueError('requested clean benchmark has30 repetitions')
    if not args.output.resolve().is_relative_to(TRACK.resolve()) or args.output.exists():raise FileExistsError('new Track C benchmark destination required')
    require_idle();args.output.mkdir(parents=True)
    parents=expected_parents();contexts={c:load_contexts(parents,c) for c,_,_ in [ARMS[arm] for arm in args.arms]}
    actions={};loads=[]
    for name in args.arms:
        config,method,closure=ARMS[name]
        for seed in SEEDS:
            tick,cpu=time.perf_counter(),time.process_time()
            if closure=='q':
                from track_c.kkt_closure.closure_action import Action
                action=Action(method=method,seed=seed,closure_mode='q',noise_ratio=1.)
            else:action=load_action(method,seed)
            actions[name,seed]=action
            loads.append(dict(arm=name,seed=seed,load_wall_seconds=time.perf_counter()-tick,load_CPU_seconds=time.process_time()-cpu))
            predict_contexts([contexts[config][0]],action,[parents[0]],1)
    rng=np.random.default_rng(202610270630);rows=[];started,cpu0=time.perf_counter(),time.process_time()
    try:
        with (args.output/'TIMINGS.csv').open('x') as stream:
            writer=csv.DictWriter(stream,fieldnames=['repeat','arm','seed','parent','wall_seconds','CPU_seconds']);writer.writeheader()
            for repeat in range(30):
                require_idle();order=rng.permutation(len(parents));armorder=rng.permutation(len(args.arms))
                for k in armorder:
                    name=args.arms[k];config,_,_=ARMS[name]
                    for seed in SEEDS:
                        for i in order:
                            tick,cpu=time.perf_counter(),time.process_time()
                            prediction=predict_contexts([contexts[config][i]],actions[name,seed],[parents[i]],1)[0]
                            wall,elapsed_cpu=time.perf_counter()-tick,time.process_time()-cpu
                            if prediction.shape!=(576,) or not np.isfinite(prediction).all():raise ValueError('finite benchmark prediction required')
                            row=dict(repeat=repeat,arm=name,seed=seed,parent=parents[i],wall_seconds=wall,CPU_seconds=elapsed_cpu)
                            writer.writerow(row);rows.append(row)
                stream.flush()
                if (repeat+1)%5==0:print(json.dumps(dict(event='CLEAN_INFERENCE_BENCHMARK',repeat=repeat+1)),flush=True)
        summary={}
        for name in args.arms:
            selected=[row for row in rows if row['arm']==name]
            wall=np.array([row['wall_seconds'] for row in selected]);cpu=np.array([row['CPU_seconds'] for row in selected])
            repeatmeans=[float(np.mean([row['wall_seconds'] for row in selected if row['repeat']==r])) for r in range(30)]
            summary[name]=dict(parent_calls=len(selected),mean_wall_seconds=float(wall.mean()),median_wall_seconds=float(np.median(wall)),
                p95_wall_seconds=float(np.quantile(wall,.95)),mean_CPU_seconds=float(cpu.mean()),repeat_mean_wall_seconds=repeatmeans)
        write(args.output/'BENCHMARK.json',dict(status='COMPLETE',repetitions=30,parents=64,seeds=list(SEEDS),arms=summary,
            inference_backend='CPU1; single-parent calls; loaded warm checkpoints; fixed tensors/weights',loads_separate=loads,
            new_RHS=0,GT_read=False,job_wall_seconds=time.perf_counter()-started,job_CPU_seconds=time.process_time()-cpu0,
            physical_anchor_time_included=False,scope='head inference only; standalone physical full/low D timing remains separate',
            uncontrolled_other_OS_load='not measured; controlled research physics/training idle at each repeat boundary'))
    except BaseException as exc:
        write(args.output/'PARTIAL.json',dict(status='INCOMPLETE',error=repr(exc),completed_parent_calls=len(rows),new_RHS=0,GT_read=False,
            job_wall_seconds=time.perf_counter()-started,job_CPU_seconds=time.process_time()-cpu0));raise


if __name__=='__main__':main()
