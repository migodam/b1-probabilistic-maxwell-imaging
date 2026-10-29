"""Fresh observation-only four-arm holdout; freeze predictions before scoring."""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import sys,json,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
P=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P/'src'));sys.path.insert(0,str(P/'experiments'))
from physics import Model,ReducedKernel
from validate import CountingBook
from run_quotient_benchmark import algebraic,fullj,qualified,compare
from run_branch_study import write

METHODS=['quotient','six','ReducedFullJ8','ReducedFullJ64']
SOURCES=['src/quotient_solver.py','src/branch_enumerator.py','src/physics.py','src/algebraic_solver.py',
         'experiments/run_quotient_benchmark.py','experiments/run_quotient_holdout.py','src/validate.py']

def snapshot(paths):
    return [dict(path=path,bytes=(P/path).stat().st_size,mtime_ns=(P/path).stat().st_mtime_ns) for path in paths]

def freeze():
    path=P/'experiments/QUOTIENT_HOLDOUT_RUN_CONTRACT.json'
    if path.exists():raise RuntimeError('Holdout run already frozen')
    original=json.loads((P/'experiments/QUOTIENT_HOLDOUT_CONTRACT.json').read_text())
    now={r['path']:r for r in snapshot([r['path'] for r in original['source_snapshot']])}
    if any(now[row['path']]!=row for row in original['source_snapshot']):
        raise RuntimeError('Pre-generation solver source snapshot changed; stop rather than tune')
    high=json.loads((P/'experiments/HIGH4_CONTRACT.json').read_text())
    control=json.loads((P/'experiments/QUOTIENT_CONTROL64_CONTRACT.json').read_text())
    s8=high['initialization']['FullJ8_starts'];s64=control['starts']['high4_E2_HIGH4_01']
    if not np.array_equal(s8,np.asarray(s64)[:8]):raise RuntimeError('Frozen8startprefix mismatch')
    write(path,dict(schema='Z-ALG-R2-FRESH-QUOTIENT-HOLDOUT-RUN-v1',frozen_utc=datetime.now(timezone.utc).isoformat(),
        original_contract='experiments/QUOTIENT_HOLDOUT_CONTRACT.json',methods=METHODS,
        reference_material=[.5]*4,shared_setup='onegeometryandrank4kernelsetup perarm acrossall4targets,6RHSperarm',
        starts8=s8,starts64=s64,source_snapshot=snapshot(SOURCES),
        source_policy='Original4solver sources unchanged frompre-generation snapshot',
        parameters='Frozenquotient,six,FullJTRF andvalidation settings reused withouttuning',
        repetitions=1,truth='No private folder,truth values,seed,orscores read',
        generation_accounting=dict(declared_source_RHS=8,derivation='4targets×2sourceRHS;parentgeneratorreceipt separate',
                                   wall_seconds=None,status='GENERATOR_WALL_NOT_PROVIDED_TO_SOLVER'),
        prediction_freeze='All16predictionfiles andsourcesnapshotreceipt beforeparentscore notification'))
    print('FRESH_HOLDOUT_RUN_FROZEN',flush=True)

def main():
    if '--freeze' in sys.argv:freeze();return
    contract=json.loads((P/'experiments/QUOTIENT_HOLDOUT_RUN_CONTRACT.json').read_text())
    original=json.loads((P/'experiments/QUOTIENT_HOLDOUT_CONTRACT.json').read_text())
    repair_path=P/'experiments/QUOTIENT_HOLDOUT_DRIVER_REPAIR1.json'
    if '--freeze-repair' in sys.argv:
        if repair_path.exists():raise RuntimeError('Driver repair already frozen')
        frozen_core={row['path']:row for row in original['source_snapshot']}
        if any(row!=frozen_core[row['path']] for row in snapshot(list(frozen_core))):
            raise RuntimeError('Frozen solver core changed')
        write(repair_path,dict(status='DRIVER_ONLY_PRE_PHYSICS_REPAIR_FROZEN',
            frozen_utc=datetime.now(timezone.utc).isoformat(),
            reason='No-repeat guard matched parent public generation receipt;now checks actual QH4 predictions only',
            prior_attempt_actual_RHS=0,prior_attempt_predictions=0,prior_attempt_observations_loaded=0,
            scientific_parameters_changed=False,core_solver_sources_changed=False,source_snapshot=snapshot(SOURCES)))
        print('HOLDOUT_DRIVER_REPAIR_FROZEN',flush=True);return
    expected_sources=json.loads(repair_path.read_text())['source_snapshot'] if repair_path.exists() else contract['source_snapshot']
    if snapshot(SOURCES)!=expected_sources:raise RuntimeError('Sources changed afterrunfreeze')
    if list((P/'results').glob('quotient_holdout_QH4_*.json')):raise RuntimeError('Holdout predictions alreadyexist;no repeat')
    manifest=json.loads((P/'experiments/quotient_holdout_public/manifest.json').read_text())
    observations={row['case']:np.load(P/'experiments/quotient_holdout_public'/row['observations'])['y'] for row in manifest}
    all_outputs=[];receipts=[];prediction_paths=[];whole_begin=time.perf_counter()
    for method in METHODS:
        begin=time.perf_counter();cpu=time.process_time();book=CountingBook();model=Model(book,**original['geometry'])
        kernel=ReducedKernel(model,np.full(4,.5),active=np.arange(4))
        for entry in manifest:
            case=entry['case'];before=model.counts.copy();rhs=book.rhs;ev=len(book.events);tick=time.perf_counter()
            row=dict(id=case,case=case,y=observations[case],box=original['box'],fixed=np.full(4,15.25),
                     starts=np.asarray(contract['starts64'] if method=='ReducedFullJ64' else contract['starts8']))
            output=algebraic(model,kernel,row,method) if method in ('quotient','six') else fullj(model,row,kernel)
            output.update(case=case,method=method,box=original['box'],reference_material=[.5]*4,
                exposure='FRESH_OBSERVATION_HOLDOUT_TRUTH_WITHHELD_FROM_SOLVER',
                incremental_RHS=book.rhs-rhs,model_counts=dict(model.counts-before),events=book.events[ev:],
                task_wall_before_logging=time.perf_counter()-tick)
            path=P/'results'/f'quotient_holdout_{case}_{method}.json';write(path,output)
            all_outputs.append(output);prediction_paths.append(str(path.relative_to(P)))
            print('HOLDOUT',case,method,output['status'],'qualified',len(qualified(output)),flush=True)
        receipt=dict(method=method,target_count=4,setup_RHS=6,accounting=book.summary(model.ns,model.nr),
            total_wall_including_target_logging=time.perf_counter()-begin,total_CPU_seconds=time.process_time()-cpu,
            model_counts=dict(model.counts),qualified_roots=sum(len(qualified(o)) for o in all_outputs if o['method']==method),
            failures=sum(o['status']!='COMPLETE' for o in all_outputs if o['method']==method),
            fallbacks=sum(o.get('fallback_used',False) for o in all_outputs if o['method']==method))
        write(P/'results'/f'quotient_holdout_{method}_accounting.json',receipt);receipts.append(receipt)
    comparison=[]
    for entry in manifest:
        selected={o['method']:o for o in all_outputs if o['case']==entry['case']}
        comparison.append(dict(case=entry['case'],qualified_counts={m:len(qualified(o)) for m,o in selected.items()},
            quotient_vs_six=compare(selected['quotient'],selected['six']),
            quotient_vs_ReducedFullJ8=compare(selected['quotient'],selected['ReducedFullJ8']),
            quotient_vs_ReducedFullJ64=compare(selected['quotient'],selected['ReducedFullJ64'])))
    summary=dict(schema='Z-ALG-R2-FRESH-QUOTIENT-HOLDOUT-SUMMARY-v1',status='PREDICTIONS_COMPLETE_UNSCORED',
        methods=receipts,cases=comparison,new_solver_actual_RHS=sum(r['accounting']['total_complex_rhs'] for r in receipts),
        generation=contract['generation_accounting'],solver_whole_wall_including_logging=time.perf_counter()-whole_begin,
        source_policy='Frozen solver sources beforegeneration andunchanged duringholdout',truth_access='NONE',
        certificate_flag='NOT_RIGOROUSLY_CERTIFIED')
    summary['declared_solver_plus_generation_RHS']=summary['new_solver_actual_RHS']+contract['generation_accounting']['declared_source_RHS']
    write(P/'results/quotient_holdout_summary.json',summary)
    source_final=snapshot(SOURCES)
    if source_final!=expected_sources:raise RuntimeError('Sources changed duringholdout;predictionfreeze rejected')
    freeze_receipt=dict(schema='Z-ALG-R2-FRESH-HOLDOUT-PREDICTION-FREEZE-v1',frozen_utc=datetime.now(timezone.utc).isoformat(),
        status='ALL_PREDICTIONS_FROZEN_BEFORE_SCORING',prediction_count=len(prediction_paths),
        predictions=snapshot(prediction_paths),source_snapshot=source_final,truth_access='NONE',
        summary='results/quotient_holdout_summary.json',run_contract='experiments/QUOTIENT_HOLDOUT_RUN_CONTRACT.json',
        driver_repair='experiments/QUOTIENT_HOLDOUT_DRIVER_REPAIR1.json' if repair_path.exists() else None)
    write(P/'results/quotient_holdout_FREEZE_RECEIPT.json',freeze_receipt)
    print('FRESH_HOLDOUT_PREDICTIONS_FROZEN',json.dumps({r['method']:{'roots':r['qualified_roots'],'RHS':r['accounting']['total_complex_rhs'],'wall':r['total_wall_including_target_logging']} for r in receipts}),flush=True)

if __name__=='__main__':main()
