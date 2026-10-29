#!/usr/bin/env python3
"""Replay the frozen, now-exposed general6 observations into a new directory.

Usage: python replay_general6.py --output-dir /absolute/new/output_directory
This is an exposed reproduction, not a fresh blind experiment.  No truth or
scoring file is loaded.  Reproduction costs belong to this new output directory
and are not automatically added to the original research ledger.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import argparse,json,sys,time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'experiments'))
from physics import Model,ReducedKernel
from validate import CountingBook
from run_six_material_holdout import quotient,reduced64,roots,write,snapshot,SOURCES


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',required=True,type=Path,
                        help='New or empty directory; original frozen results are never overwritten')
    args=parser.parse_args();destination=args.output_dir.expanduser().resolve()
    if destination.exists() and any(destination.iterdir()):
        parser.error('--output-dir must be new or empty')
    if destination==ROOT/'results':parser.error('Choose an output directory separate from frozen results')
    destination.mkdir(parents=True,exist_ok=True)
    contract=json.loads((ROOT/'experiments/GENERAL6_CONTRACT.json').read_text())
    frozen=json.loads((ROOT/'results/general6_SOLVER_FREEZE_RECEIPT.json').read_text())
    manifest=json.loads((ROOT/'experiments/general6_public/manifest.json').read_text())
    inputs=[dict(case=r['case'],file=ROOT/'experiments/general6_public'/r['observation_file']) for r in manifest['cases']]
    outputs=[];receipts=[];stopped=False;begin=time.perf_counter()
    write(destination/'REPLAY_CONTRACT.json',dict(exposure='EXPOSED_REPRODUCTION_NOT_FRESH_BLIND',
        frozen_parameters=frozen['parameters'],optimizer=frozen['optimizer'],
        observation_inputs=[dict(case=r['case'],file=str(r['file'])) for r in inputs],
        starts=contract['initialization']['ReducedFullJ64_starts'],source_snapshot_current=snapshot(SOURCES),
        source_note='Portable checkout mtimes are recorded,not required to equal the original machine freeze.',
        truth_access='NONE',cost_note='New replay costs,not part of original research ledger'))
    for method in ['general_quotient20','ReducedFullJ64']:
        if stopped:break
        tick=time.perf_counter();book=CountingBook();model=Model(book,**contract['geometry'])
        kernel=ReducedKernel(model,np.asarray(frozen['parameters']['reference_material']),active=np.arange(6))
        for row in inputs:
            case=row['case'];y=np.load(row['file'])['y'];rhs=book.rhs;events=len(book.events)
            output=quotient(model,kernel,case,y,contract['material_box']) if method=='general_quotient20' else reduced64(model,kernel,case,y,contract['material_box'],contract['initialization']['ReducedFullJ64_starts'])
            output.update(exposure='EXPOSED_REPRODUCTION_NOT_FRESH_BLIND',incremental_RHS=book.rhs-rhs,events=book.events[events:])
            write(destination/f'{case}_{method}.json',output);outputs.append(output)
            print(case,method,output['status'],'qualified roots',len(roots(output)),flush=True)
            if output['status']!='COMPLETE':stopped=True;break
        receipt=dict(method=method,accounting=book.summary(model.ns,model.nr),
                     wall_including_logging=time.perf_counter()-tick,model_counts=dict(model.counts))
        receipts.append(receipt);write(destination/f'{method}_accounting.json',receipt)
    summary=dict(status='STOPPED_ON_NUMERICAL_GATE' if stopped else 'COMPLETED_EXPOSED_REPRODUCTION',
        exposure='EXPOSED_REPRODUCTION_NOT_FRESH_BLIND',truth_access='NONE',
        new_actual_full_RHS=sum(r['accounting']['total_complex_rhs'] for r in receipts),
        wall_including_logging=time.perf_counter()-begin,
        qualified_root_counts=[dict(case=o['case'],method=o['method'],count=len(roots(o))) for o in outputs],
        certificate_flag='NOT_RIGOROUSLY_CERTIFIED')
    write(destination/'REPLAY_SUMMARY.json',summary)
    print(json.dumps(summary),flush=True)
    return 1 if stopped else 0


if __name__=='__main__':raise SystemExit(main())
