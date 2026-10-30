#!/usr/bin/env python3
"""Anonymous CSI batch with one paid, data-independent fixed-reference bank."""
from pathlib import Path
import argparse
import json
import os
import re
import resource
import sys
import time
sys.dont_write_bytecode=True
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[name]='1'
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from gsda.accounting import BudgetBook,BudgetExceeded,dump
from gsda.fem_model import FEMModel
from gsda.currents import FixedReference,CSIOptions,invert
from gsda.fwi import FWIOptions
from gsda.guard import OnlineGuard
from run_currents import parse_options

PHYSICS_KEYS=('frequencies','heldout_frequencies','h','targetstep','wall_refine',
              'wall_x','wall_y','target_x_edges','target_y_edges','target_degree','order',
              'tx','rx','source_indices','pml','quadrature_order')


def physics_signature(config):
    from gsda.fem_model import default_config
    defaults=default_config();defaults.update(order=3,target_degree=8,source_indices=[0,1,2,3])
    signature={key:config.get(key,defaults.get(key)) for key in PHYSICS_KEYS}
    signature['sigma_real_imag']=config.get('sigma_real_imag',config.get('sigma',.01))
    return signature


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,required=True)
    ap.add_argument('--options',type=Path,required=True)
    ap.add_argument('--online-root',type=Path,default=ROOT/'inputs/online')
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--cases',nargs='+',default=[f'c{i:03d}' for i in range(1,10)])
    ap.add_argument('--method',choices=('B1','B1_RESET','B2','B3','B4','B5'),default='B4')
    ap.add_argument('--job-prefix',default='gsda_batch')
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    if len(set(args.cases))!=len(args.cases) or any(re.fullmatch(r'c\d{3}',c) is None for c in args.cases):
        ap.error('cases must be distinct anonymous public cNNN directories')
    if (args.output/'BATCH_RECEIPT.json').exists():raise RuntimeError('refusing overwrite of an existing batch')
    guard=OnlineGuard(ROOT).install();batch_cpu=time.process_time();batch_wall=time.perf_counter()
    book=BudgetBook(ROOT,job=args.job_prefix+'_public_setup');pid=os.getpid()
    records=[];partition_receipts=[];setup_receipt=None;fatal=None;status='FAILED';active_partition=None

    def begin_partition(job):
        nonlocal active_partition
        book.job=job;book.cpu0=time.process_time();book.wall0=time.perf_counter()
        active_partition=dict(job=job,before=book.snapshot(pid=pid))

    def end_partition(partition_status,output,**metadata):
        nonlocal active_partition
        after=book.snapshot(pid=pid);before=active_partition['before']
        row=dict(job=book.job,pid=pid,status=partition_status,
            wrapper_process_cpu_seconds=time.process_time()-book.cpu0,
            wall_seconds=time.perf_counter()-book.wall0,
            peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),
            rhs=after['rhs']-before['rhs'],full_fe=after['full_fe']-before['full_fe'],
            output=str(output),partition_only=True,**metadata)
        dump(ROOT/'receipts'/f'{book.job}_{pid}.json',row)
        partition_receipts.append(row);active_partition=None
        return row

    try:
        config=json.loads(args.config.read_text());freeze=json.loads(args.options.read_text())
        csi,prior=parse_options(freeze);csi=dict(csi)
        csi['spectral_preconditioner']=args.method in ('B3','B4','B5')
        if args.method=='B1_RESET':csi['physical_reset_at']=csi.get('absorption_at',CSIOptions().absorption_at)
        else:csi['physical_reset_at']=0
        options=CSIOptions(**csi);FWIOptions(**prior)
        observation_config=None
        chosen_physics=physics_signature(config)
        for case in args.cases:
            if (args.output/case/'ONLINE_SEALED.json').exists():raise RuntimeError('sealed case output already exists: '+case)
            cfg=json.loads((args.online_root/case/'config.json').read_text());cfg.pop('case_id',None)
            measured_physics=physics_signature(cfg)
            differences={key:dict(chosen=chosen_physics[key],measurement=measured_physics[key])
                         for key in chosen_physics if chosen_physics[key]!=measured_physics[key]}
            if differences:raise ValueError('chosen shared physics differs from '+case+': '+json.dumps(differences))
            if observation_config is None:observation_config=cfg
            elif cfg!=observation_config:raise ValueError('public measurement configurations differ across cases')
        dump(args.output/'FROZEN_BATCH_CONTRACT.json',dict(cases=args.cases,method=args.method,
            csi_options=csi,fwi_options=prior,model_config=config,measurement_config=observation_config,
            asserted_physics_keys=list(chosen_physics),asserted_physics=chosen_physics,
            shared_only=['air/PML sparse L0','actual quadrature P/weights','GS0','public GD right modes Q/CQ','public GS Gram'],
            data_dependent_carry=False,case_starts='two fresh public material starts, current seeds, histories and endpoint selection per case',
            truth_access=False,ordinary_B5_same_public_cache_permission=True,
            standalone_cli='cold setup per invocation; this explicitly declared batch charges shared setup once'))
        begin_partition(args.job_prefix+'_public_setup');book.check(new_job=True)
        model=FEMModel(config,book,phase='R4_SHARED_PUBLIC_CSI')
        reference=FixedReference(model);bank=args.output/'public_reference'
        reference.save_public(bank)
        if options.spectral_preconditioner:reference.prepare_modes(options.spectral_rank,bank)
        setup_receipt=end_partition('PUBLIC_FIXED_REFERENCE_SETUP_COMPLETE',bank,
            data_dependent_arrays=False,shared_setup_charged_once=True,requested_case_count=len(args.cases))
        dump(args.output/'PUBLIC_OPERATOR_SETUP_RECEIPT.json',setup_receipt)
        for index,case in enumerate(args.cases):
            caseout=args.output/case;caseout.mkdir(parents=True,exist_ok=True)
            begin_partition(f'{args.job_prefix}_{case}')
            try:
                book.check(new_job=True)
                model.config['case_id']=case
                with np.load(args.online_root/case/'observation.npz',allow_pickle=False) as source:
                    if set(source.files)!={'y'}:raise ValueError('online observation must contain only y')
                    y=source['y'].copy()
                result=invert(model,y,caseout,csi,prior,
                    absorption=args.method in ('B2','B3','B4'),risk_controlled=args.method=='B4',
                    shared_reference=reference)
                del y
                row=dict(case=case,index=index,status=result['status'],objective=result['objective'],
                         original_air_reference_restored=result['shared_original_air_reference_restored'],
                         shared_public_reference=str(bank.resolve()))
                receipt=end_partition(result['status'],caseout,case=case,warm_case_only=True,
                    shared_setup_included=False,public_reference=str(bank.resolve()),
                    data_dependent_carry=False)
                dump(caseout/'RECEIPT.json',receipt);row['warm_cost']=receipt
            except BaseException as error:
                book.failure('R4_SHARED_PUBLIC_CSI',case,error)
                row=dict(case=case,index=index,status='FAILED',error=dict(type=type(error).__name__,message=str(error)))
                receipt=end_partition('FAILED',caseout,case=case,warm_case_only=True,shared_setup_included=False)
                dump(caseout/'RECEIPT.json',receipt);dump(caseout/'FAILURE.json',row);row['warm_cost']=receipt
                if isinstance(error,(BudgetExceeded,KeyboardInterrupt,SystemExit)):
                    records.append(row);raise
            records.append(row)
            dump(args.output/'BATCH_PROGRESS.json',dict(records=records,
                remaining=args.cases[len(records):],shared_setup=setup_receipt))
        status='COMPLETE_BOUNDED_BATCH' if all(r['status']=='COMPLETE_BOUNDED_ATTEMPT' for r in records) else 'PARTIAL_BOUNDED_BATCH'
    except BaseException as error:
        fatal=dict(type=type(error).__name__,message=str(error));book.failure('R4_SHARED_PUBLIC_CSI','batch',error)
        status='PARTIAL_BOUNDED_BATCH' if records else 'FAILED'
        if active_partition is not None:
            receipt=end_partition('FAILED',args.output,error=fatal)
            if setup_receipt is None:setup_receipt=receipt
        if isinstance(error,(KeyboardInterrupt,SystemExit)):raise
    finally:
        aggregate=dict(job=args.job_prefix,pid=pid,status=status,
            wrapper_process_cpu_seconds=time.process_time()-batch_cpu,wall_seconds=time.perf_counter()-batch_wall,
            peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),
            **book.snapshot(pid=pid),output=str(args.output),aggregate_only=True,
            accounting='ledger reservations authoritative; aggregate stored outside receipts to prevent duplicate CPU aggregation')
        aggregate['unpartitioned_process_cpu_seconds']=max(0.,aggregate['wrapper_process_cpu_seconds']-
            sum(r['wrapper_process_cpu_seconds'] for r in partition_receipts))
        dump(args.output/'BATCH_RECEIPT.json',aggregate)
        denominator=len(args.cases)
        amortization=None
        if setup_receipt is not None:
            shared={name:setup_receipt[name]/denominator for name in ('rhs','full_fe','wrapper_process_cpu_seconds','wall_seconds')}
            completed_count=sum(r['status']=='COMPLETE_BOUNDED_ATTEMPT' for r in records)
            completed_shared=None if not completed_count else {name:setup_receipt[name]/completed_count for name in shared}
            first=records[0].get('warm_cost') if records else None
            cold=None if first is None else {name:setup_receipt[name]+first[name] for name in shared}
            amortization=dict(requested_cases=denominator,attempted_cases=len(records),
                actual_completed_count=completed_count,
                setup_charged_once=setup_receipt,setup_per_completed_case=completed_shared,
                setup_per_requested_case_forecast_only=shared,
                cold_first_case_setup_plus_first_warm=cold,
                scope='public fixed-reference reuse, equally available to ordinary B5; not an absorption-specific advantage')
        dump(args.output/'BATCH_SUMMARY.json',dict(status=status,method=args.method,cases=args.cases,
            records=records,not_run=args.cases[len(records):],error=fatal,amortization=amortization,
            partition_receipts=partition_receipts,aggregate=aggregate,
            data_dependent_carry=False,truth_access=False))
        guard.save(args.output/'FILE_ACCESS_SUMMARY.json')
    print(json.dumps(dict(status=status,attempted=len(records),output=str(args.output),
                          cost=book.snapshot(pid=pid))),flush=True)
    if status!='COMPLETE_BOUNDED_BATCH':raise SystemExit(1)

if __name__=='__main__':main()
