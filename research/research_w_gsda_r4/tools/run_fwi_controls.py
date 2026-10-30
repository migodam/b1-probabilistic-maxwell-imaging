"""Execute the predeclared four FWI jobs only after explicit queue grant.

Whole-child CPU receipts are separate from historical Book CPU values. The
manifest supplies all numerical settings; this launcher never opens truth.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,os,resource,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from gsda.accounting import BudgetBook,dump


def utc_now():return datetime.now(timezone.utc).isoformat()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--manifest',type=Path,default=ROOT/'FWI_CONTROL_ABLATION_FREEZE.json')
    parser.add_argument('--execute',action='store_true',help='Only use after parent grants the heavy queue')
    args=parser.parse_args();manifest=json.loads(args.manifest.read_text());jobs=manifest['jobs']
    if not args.execute:
        print(json.dumps(dict(status='PREPARED_NOT_EXECUTED',jobs=len(jobs),new_physics=dict(full_fe=0,rhs=0),child_cpu_coverage='whole child via RUSAGE_CHILDREN around subprocess.run')));return
    records=[];end_to_end=[]
    for job in jobs:
        output=Path(job['output']);record=dict(case=job['case'],campaign=job['campaign'],status='NOT_RUN',argv=job['argv'])
        try:
            actual_options=json.loads((ROOT/job['options']).read_text())['fwi_options']
            expected_options=manifest['ablation_options'] if job['campaign']=='B0_TV0_ABLATION' else manifest['main_options']
            if actual_options!=expected_options:raise RuntimeError('numerical options changed after the predeclared freeze')
            book=BudgetBook(ROOT,job='FWI_control_launcher_preflight');book.check(new_job=True)
            if (output/'ONLINE_SEALED.json').exists():raise RuntimeError('refusing to rerun sealed job '+str(output))
            output.mkdir(parents=True,exist_ok=True)
            dump(output/'EXPERIMENT_SCOPE.json',dict(freeze=str(args.manifest),case=job['case'],campaign=job['campaign'],argv=job['argv'],scope=job['postseal_metadata']))
            before=resource.getrusage(resource.RUSAGE_CHILDREN);started_utc=utc_now();start_perf=time.perf_counter();start_unix=time.time()
            launch_error=None;returncode=None
            try:
                with (output/'stdout.log').open('w') as log:
                    result=subprocess.run(job['argv'],stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'})
                returncode=result.returncode
            except BaseException as error:
                launch_error=dict(type=type(error).__name__,message=str(error))
                if isinstance(error,(KeyboardInterrupt,SystemExit)):raise
            finally:
                elapsed=time.perf_counter()-start_perf;ended_utc=utc_now();after=resource.getrusage(resource.RUSAGE_CHILDREN)
                cpu_user=after.ru_utime-before.ru_utime;cpu_system=after.ru_stime-before.ru_stime
                child_receipt=dict(case=job['case'],campaign=job['campaign'],coverage='WHOLE_CHILD_PROCESS_INCLUDING_BOOTSTRAP_IMPORTS_SOLVES_AND_OUTPUTS',measurement='RUSAGE_CHILDREN user+system delta around blocking subprocess.run; this launcher has no parallel child',child_user_cpu_seconds=cpu_user,child_system_cpu_seconds=cpu_system,child_total_cpu_seconds=cpu_user+cpu_system,external_wall_seconds=elapsed,started_utc=started_utc,ended_utc=ended_utc,started_unix=start_unix,ended_unix=time.time(),argv=job['argv'],returncode=returncode,launch_error=launch_error,additive_with_authoritative_solver_receipt=False,do_not_double_add=True)
                old_path=output/'RECEIPT.json'
                if old_path.exists():
                    old=json.loads(old_path.read_text());old_cpu=old.get('wrapper_process_cpu_seconds')
                    child_receipt.update(child_pid=old.get('pid'),legacy_book_cpu_seconds=old_cpu,legacy_book_cpu_scope='from Book construction after imports through solver receipt generation',authoritative_solver_receipt=str(old_path))
                    if old_cpu is not None:child_receipt['observed_whole_child_minus_book_cpu_seconds']=cpu_user+cpu_system-old_cpu
                dump(output/'END_TO_END_CPU.json',child_receipt);end_to_end.append(child_receipt)
            if (output/'ONLINE_SEALED.json').exists():
                seal=json.loads((output/'ONLINE_SEALED.json').read_text())
                seal['inverter_method_before_scope_label']=seal.get('method')
                seal.update(job['postseal_metadata']);seal['scope_metadata_attached_unix']=time.time()
                dump(output/'ONLINE_SEALED.json',seal)
                record['online_status']=seal['status']
                record['objective']=seal.get('objective')
                record['selected_scaled_projected_gradient_inf']=seal.get('selected_scaled_projected_gradient_inf')
            record.update(returncode=returncode,status='PROCESS_COMPLETED' if returncode==0 else 'FAILED',launch_error=launch_error)
            record['missing_required_outputs']=[name for name in manifest['required_outputs'] if not (output/name).is_file() or (output/name).stat().st_size==0]
        except Exception as error:
            record['reason']=str(error)
        records.append(record)
        dump(ROOT/'FWI_CONTROL_ABLATION_EXECUTION.json',dict(freeze=str(args.manifest),attempts=records,remaining=[row['case'] for row in jobs[len(records):]],no_other_physics_appended=True))
        dump(ROOT/'FWI_CONTROL_ABLATION_END_TO_END_CPU.json',dict(coverage='whole child only, separate from existing solver receipts',receipts=end_to_end,total_child_cpu_seconds=sum(row['child_total_cpu_seconds'] for row in end_to_end),total_external_wall_seconds=sum(row['external_wall_seconds'] for row in end_to_end),additive_with_solver_receipts=False,do_not_double_add=True))
        print(json.dumps(record),flush=True)

if __name__=='__main__':main()
