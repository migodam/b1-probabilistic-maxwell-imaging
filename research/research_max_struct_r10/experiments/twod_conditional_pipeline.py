"""One frozen final4 conditional cohort using the existing finite CPU queue.

Only publicY/sigma/paid physical predictions trigger/choose attempts. No GT.
Shared prefix is charged once in production, repeated only in deployment views.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import argparse,json,sys,subprocess,time,shutil,copy
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]


def dump(p,d):
    p=Path(p);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n');tmp.replace(p)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--execute',action='store_true')
    ap.add_argument('--cpus',nargs='+',type=int,default=[8,10]);ap.add_argument('--vendor-path',type=Path)
    ap.add_argument('--resume-startup-guard',action='store_true',help='One explicit zero-physics path-guard startup repair; never retry a numerical solve')
    args=ap.parse_args()
    if not args.execute:raise ValueError('Explicit --execute required')
    import psutil
    process=psutil.Process();process.cpu_affinity(args.cpus)
    if os.name=='nt':process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    if args.vendor_path:sys.path.insert(0,str(args.vendor_path))
    from twod_queue import external_cpu
    out=ROOT/'delegated/twod/conditional_confirmation'
    if not args.resume_startup_guard:out.mkdir(parents=True,exist_ok=False)
    receipt_path=out/'RECEIPT.json';cfg=ROOT/'configs/twod_conditional_solver.json'
    prefix_manifest=ROOT/'configs/twod_conditional_prefix_tasks.json'
    start=time.perf_counter();receipt=dict(status='RUNNING_PREFIX4',started_utc=datetime.now(timezone.utc).isoformat(),
        truth_access='NONE',new_forward_calls_for_trigger_or_selection=0,queues=[],triggers=[],deployments=[],shared_prefix_charged_once=True)
    old_controller_CPU=0.
    if args.resume_startup_guard:
        old=json.loads(receipt_path.read_text())
        if old['status']!='ALL_ACTUAL_CONDITIONAL_PREDICTIONS_SEALED':raise ValueError('Only resume after original pipeline terminal')
        if (out/'PIPELINE_BEFORE_STARTUP_REPAIR.json').exists():raise FileExistsError('One controlled startup repair already used')
        dump(out/'PIPELINE_BEFORE_STARTUP_REPAIR.json',old)
        for name in ('TRIGGER_DECISIONS.json','FALLBACK_TASKS.json'):
            shutil.copy2(out/name,out/('BEFORE_STARTUP_REPAIR_'+name))
        old_controller_CPU=old['pipeline_controller_CPU_seconds']
        receipt.update(queues=old['queues'],status='CONTROLLED_ZERO_PHYSICS_STARTUP_REPAIR',
            original_pipeline_receipt='PIPELINE_BEFORE_STARTUP_REPAIR.json',
            qualification_note='One worker rejected output guard before model/forward; preserve its0.125s cost and retry only its previously unexecuted solve with unchanged numerical contract')
    dump(receipt_path,receipt)
    def run_queue(manifest,name):
        directory=out/name
        command=[sys.executable,str(ROOT/'experiments/twod_queue.py'),'--execute','--manifest',str(manifest),'--config',str(cfg),
            '--receipt-dir',str(directory),'--cpus',*[str(v) for v in args.cpus]]
        if args.vendor_path:command+=['--vendor-path',str(args.vendor_path)]
        with (out/(name+'.log')).open('x') as fh:
            p=subprocess.Popen(command,cwd=REPO,stdout=fh,stderr=subprocess.STDOUT);code=p.wait();cost=external_cpu(p)
        q=json.loads((directory/'RECEIPT.json').read_text())
        row=dict(name=name,command=command,returncode=code,external_queue_process_cost=cost,
            internal_controller_CPU_seconds=q['controller_CPU_seconds'],whole_inverse_children_CPU_seconds=q['completed_child_CPU_seconds'],
            jobs=q['jobs'],status=q['status'],receipt=str((directory/'RECEIPT.json').relative_to(ROOT)))
        receipt['queues'].append(row);dump(receipt_path,receipt)
        if q['status'] not in ('COMPLETED','COMPLETED_WITH_ERRORS'):raise RuntimeError('Queue did not reach all planned terminal outcomes')
        return q
    try:
        import numpy as np
        prefix=json.loads(prefix_manifest.read_text())
        if args.resume_startup_guard:
            q=json.loads((out/'prefix_queue/RECEIPT.json').read_text())
            repairs=[t for t in prefix['tasks'] if not (REPO/t['output_dir']/'RECEIPT.json').exists()]
            if len(repairs)!=1:raise ValueError('Expected exactly one unexecuted startup, never retry paid numerical runs')
            oldjob=next(j for j in q['jobs'] if j['task_id']==repairs[0]['task_id'])
            if oldjob.get('physical_model_cost') or 'All R10 outputs must stay inside' not in oldjob.get('error_log_tail',''):
                raise ValueError('Repair is restricted to the observed zero-physics output-guard error')
            (REPO/repairs[0]['output_dir']).parent.mkdir(parents=True,exist_ok=True)
            rp=out/'STARTUP_REPAIR_PREFIX_TASK.json';dump(rp,dict(input_base='.',tasks=repairs))
            rq=run_queue(rp,'startup_repair_queue')
            replacement=rq['jobs'][0]
            if replacement['status']!='FINISHED':raise RuntimeError('Controlled startup repair failed; preserve, do not retry again')
            collapsed=dict(replacement,child_process_CPU_seconds=replacement['child_process_CPU_seconds']+oldjob['child_process_CPU_seconds'],
                engineering_startup_CPU_seconds=oldjob['child_process_CPU_seconds'],
                deployment_cost_scope='Valid prefix plus preserved zero-physics failed startup; production rows retain both attempts separately')
            q['jobs']=[collapsed if j['task_id']==collapsed['task_id'] else j for j in q['jobs']]
        else:q=run_queue(prefix_manifest,'prefix_queue')
        prefix_jobs={j['task_id']:j for j in q['jobs']};fallbacks=[];prefix_rows={}
        for task in prefix['tasks']:
            case=task['case_id'];p=REPO/task['output_dir'];r=json.loads((p/'RECEIPT.json').read_text()) if (p/'RECEIPT.json').exists() else {}
            gate=dict(case_id=case,prefix_task_id=task['task_id'],prefix_status=r.get('status','MISSING'),trigger=False,
                strict_rule='||paidY_prediction-publicY||_F/sigma > public noise_radius')
            valid=r.get('status')=='SEALED_PUBLIC_PREDICTION' and (p/'prediction.npz').exists()
            if valid:
                obs=REPO/task['public_files'][-1]
                with np.load(obs,allow_pickle=False) as z:Y=z['Y'];sigma=float(z['sigma']);radius=float(z['noise_radius'])
                with np.load(p/'prediction.npz',allow_pickle=False) as z:prediction=z['Y_prediction']
                statistic=float(np.linalg.norm(prediction-Y)/sigma)
                gate.update(status='EVALUATED',whitened_port_residual=statistic,CN99_radius=radius,trigger=bool(statistic>radius),public_sigma=sigma)
            else:gate['status']='INVALID_PREFIX_NO_COMPUTABLE_TRIGGER'
            receipt['triggers'].append(gate);prefix_rows[case]=(task,r,prefix_jobs[task['task_id']])
            if gate['trigger']:
                for label,arm in [('G19_PGP','a_plus'),('DualCIE','cie_dual')]:
                    t=copy.deepcopy(task);ident=case+'__conditional_fallback_'+label
                    t.update(task_id=ident,arm=arm,output_dir='research/research_max_struct_r10/RESULTS/conditional_confirmation/fallback/'+ident)
                    if label=='G19_PGP':t['material_readout']=dict(module_file='research/research_max_struct_r10/src/lowrank_denominator.py',callable='matrix_material_readout',prepare='prepare_green_strong',label='conditional_G19_PGP')
                    fallbacks.append(t)
        dump(out/'TRIGGER_DECISIONS.json',dict(status='PREFIX4_COMPLETE_PUBLIC_DISCREPANCY_ONLY',rows=receipt['triggers'],truth_access='NONE',new_forward_calls=0))
        fallback_manifest=out/'FALLBACK_TASKS.json';dump(fallback_manifest,dict(input_base='.',tasks=fallbacks,status='DERIVED_ONLY_FROM_FROZEN_PUBLIC_TRIGGER'))
        fallback_jobs={j['task_id']:j for oldq in receipt['queues'] if oldq['name']=='fallback_queue' for j in oldq['jobs']}
        pending_fallbacks=[t for t in fallbacks if not (REPO/t['output_dir']/'RECEIPT.json').exists()]
        if pending_fallbacks:
            receipt['status']='RUNNING_TRIGGERED_FALLBACKS';dump(receipt_path,receipt)
            pending_manifest=out/'PENDING_FALLBACK_TASKS.json';dump(pending_manifest,dict(input_base='.',tasks=pending_fallbacks))
            fq=run_queue(pending_manifest,'startup_repair_fallback_queue' if args.resume_startup_guard else 'fallback_queue');fallback_jobs.update({j['task_id']:j for j in fq['jobs']})
        if args.resume_startup_guard:
            old_deployed=ROOT/'RESULTS/conditional_confirmation/deployed'
            if old_deployed.exists():old_deployed.rename(old_deployed.with_name('deployed_before_startup_repair'))
        for gate in receipt['triggers']:
            case=gate['case_id'];t,r,j=prefix_rows[case];prefix_valid=r.get('status')=='SEALED_PUBLIC_PREDICTION'
            for label in ['G19_PGP','DualCIE']:
                attempts=[dict(task_id=t['task_id'],Phi=r.get('physical_objective'),valid=prefix_valid,path=t['output_dir'],kind='prefix')]
                fallback_cpu=0.;fallback_status='NOT_TRIGGERED'
                if gate['trigger']:
                    ft=next(x for x in fallbacks if x['case_id']==case and x['task_id'].endswith(label));fp=REPO/ft['output_dir']
                    fr=json.loads((fp/'RECEIPT.json').read_text()) if (fp/'RECEIPT.json').exists() else {}
                    fallback_status=fr.get('status','MISSING');fallback_cpu=fallback_jobs[ft['task_id']].get('child_process_CPU_seconds')
                    attempts.append(dict(task_id=ft['task_id'],Phi=fr.get('physical_objective'),valid=fallback_status=='SEALED_PUBLIC_PREDICTION',path=ft['output_dir'],kind='fallback'))
                eligible=[a for a in attempts if a['valid'] and a['Phi'] is not None]
                chosen=min(eligible,key=lambda a:a['Phi']) if eligible else None
                d=ROOT/'RESULTS/conditional_confirmation/deployed'/(case+'__CIE_then_'+label);d.mkdir(parents=True,exist_ok=False)
                row=dict(case_id=case,method='CIE_then_'+label,trigger=gate['trigger'],fallback_status=fallback_status,
                    selected=chosen,selection='Minimum same originalPhi; ties prefix; no GT',attempts=attempts,
                    prefix_full_child_CPU_seconds=j.get('child_process_CPU_seconds'),fallback_full_child_CPU_seconds=fallback_cpu,
                    prefix_engineering_startup_CPU_seconds=j.get('engineering_startup_CPU_seconds',0.),
                    deployment_full_child_CPU_seconds=(j.get('child_process_CPU_seconds')+fallback_cpu) if j.get('child_process_CPU_seconds') is not None and fallback_cpu is not None else None,
                    deployment_cost_additive=False,scope='Deployment method view shares prefix with alternative method; do not sum these views as production cost',truth_access='NONE')
                if chosen:
                    shutil.copy2(REPO/chosen['path']/'prediction.npz',d/'prediction.npz')
                    shutil.copy2(REPO/chosen['path']/'prediction.json',d/'prediction.json')
                dump(d/'SELECTION.json',row);receipt['deployments'].append(row)
        actual_jobs=[j for q in receipt['queues'] for j in q['jobs']]
        startup=[dict(task_id=j['task_id'],GetProcessTimes_CPU_seconds=j['child_process_CPU_seconds'],physical_LU=0,physical_complex_RHS=0,
            source_log=j['log'],classification='ZERO_PHYSICS_PATH_GUARD_STARTUP_NOT_A_SECOND_CIE_SOLUTION') for j in actual_jobs if
            j['status']!='FINISHED' and 'All R10 outputs must stay inside' in j.get('error_log_tail','') and not j.get('physical_model_cost')]
        receipt.update(status='ALL_ACTUAL_CONDITIONAL_PREDICTIONS_SEALED',actual_inverse_runs=len(actual_jobs)-len(startup),actual_worker_process_attempts=len(actual_jobs),
            engineering_startup_failures=startup,prefix_runs=4,
            fallback_runs=len(fallbacks),triggered_cases=sum(x['trigger'] for x in receipt['triggers']),
            failed_inverse_runs=sum(j['status']!='FINISHED' and bool(j.get('physical_model_cost')) for j in actual_jobs),
            failed_worker_startups=sum(j['status']!='FINISHED' and not j.get('physical_model_cost') for j in actual_jobs),
            actual_whole_inverse_child_CPU_seconds=sum(j.get('child_process_CPU_seconds') or 0 for j in actual_jobs),
            actual_queue_process_CPU_seconds=sum(q['external_queue_process_cost']['child_process_CPU_seconds'] or 0 for q in receipt['queues']),
            GT_gate='Root separately authorizes offline scoring after all actual terminal outcomes and deployed outputs seal')
    except BaseException as exc:
        receipt.update(status='FAILED_PIPELINE_PRESERVED',error=repr(exc));raise
    finally:
        receipt.update(pipeline_controller_CPU_seconds=old_controller_CPU+time.process_time(),
            pipeline_controller_prior_segment_CPU_seconds=old_controller_CPU,
            wall_seconds=time.perf_counter()-start,finished_utc=datetime.now(timezone.utc).isoformat())
        dump(receipt_path,receipt)
    print(json.dumps({k:receipt[k] for k in ('status','actual_inverse_runs','triggered_cases','fallback_runs','failed_inverse_runs','actual_whole_inverse_child_CPU_seconds')}),flush=True)


if __name__=='__main__':main()
