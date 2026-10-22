"""Queue the approved conventional surrogate-force pilot after locked risk.

No action occurs on import. Wait for frozen_locked_risk/CLOSED COMPLETE and
the global paid-worker lock; run linear/material x four chains serially,
then audit each arm with the unchanged reference thresholds. Never retry or
overwrite failed/partial jobs. Parent retains final scientific review.
"""
from pathlib import Path
import fcntl
import json
import shutil
import signal
import subprocess
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,write
from p1_material.prior import MaterialPrior
from p3l.reference_audit import (audit_reference,RHAT_LIMIT,MATERIAL_ESS_MIN,MEAN_RMS_MCSE_MAX)

FOLDER=OUT/'schedules/surrogate_force_hmc'
DEPENDENCY=OUT/'schedules/frozen_locked_risk/CLOSED.json'
ARMS=('linear','material')
WARMUP=1024
OWNED_SOURCE=['tools/run_p3l_surrogate_sequence.py','tools/run_p3l_surrogate_hmc.py',
    'research/p3l_local_posterior_imaging/C/surrogate_code.py']
DEPENDENCY_SOURCE=['src/p3l/hmc.py','src/p3l/kernels.py','src/p3l/core.py',
    'src/p3l/reference.py','src/p3l/reference_audit.py','src/p3l/diagnostics.py',
    'src/p1_material/prior.py']
COST_KEYS=('RHS','states','LU','controlled_wall_seconds','physics_wall_seconds','GPU_seconds')


def relative(path):return Path(path).relative_to(OUT).as_posix()


def tasks():
    rows=[]
    for arm in ARMS:
        for chain in range(4):
            rows.append(dict(name=f'{arm}_c{chain}',kind='paid_chain',arm=arm,chain=chain,
                job=f'surrogate_force_{arm}_DEV030_c{chain}',
                command=['run_p3l_surrogate_hmc.py','--arm',arm,'--chain',str(chain)]))
        rows.append(dict(name=f'audit_{arm}',kind='analysis_only_audit',arm=arm,
            job=f'audit_surrogate_force_{arm}_DEV030',
            destination=f'C/references/DEV_030/surrogate_force_{arm}'))
    return rows


def ledger_offset():
    path=OUT/'RESOURCE_LEDGER.jsonl'
    return path.stat().st_size if path.exists() else 0


def job_costs(jobs,offset=0):
    totals={j:{k:0 for k in COST_KEYS} for j in jobs}
    path=OUT/'RESOURCE_LEDGER.jsonl'
    if not path.exists():return totals
    with path.open('r') as stream:
        stream.seek(offset)
        for line in stream:
            row=json.loads(line)
            if row.get('stage')!='C' or row.get('job') not in totals:continue
            for key in COST_KEYS:totals[row['job']][key]+=row.get(key,0)
    return totals


def summed_cost(done):
    return {key:sum(row.get('cost',{}).get(key,0) for row in done) for key in COST_KEYS}


def source_revision():
    version_path=OUT/'configs/LOCKED_CODE_VERSION.json'
    version=json.loads(version_path.read_text())
    baseline=version['git_commit']
    changed=subprocess.check_output(['git','diff',baseline,'--name-only','--',
        'src/p3l','src/p1_material'],cwd=ROOT,text=True).splitlines()
    if changed:raise RuntimeError('frozen src/p3l or original prior changed; preserve queue without execution')
    return dict(git_HEAD=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        frozen_dependency_commit=baseline,dependency_version_path=relative(version_path),
        frozen_dependency_diff=changed,
        source_status=subprocess.check_output(['git','status','--short','--',*OWNED_SOURCE,*DEPENDENCY_SOURCE],
            cwd=ROOT,text=True).splitlines(),
        identity_policy='Git revision plus exact source snapshots; no digest checks')


def freeze(manifest):
    protocol=FOLDER/'PROTOCOL.json'
    controls=OUT/'AUDIT/surrogate_force/RESULTS.json'
    result=json.loads(controls.read_text())
    if (result.get('status')!='PASS_PURE_SURROGATE_FORCE_AND_SAVED_ARRAY_CONTROLS' or
            result.get('new_physical_RHS')!=0 or result.get('GT_read') is not False):
        raise RuntimeError('approved pure surrogate-force controls absent or incompatible')
    revision=source_revision()
    if protocol.exists():
        previous=json.loads(protocol.read_text())
        if previous.get('tasks')!=manifest:raise RuntimeError('existing sequence manifest differs; no overwrite')
        verify_source()
        return previous
    project_paths={name:relative(OUT/name) for name in
        ['configs/PROTOCOL.json','configs/AUTHORIZATION.json','configs/TARGET.json',
         'configs/LOCKED_CODE_VERSION.json','THEORY_AND_ALGORITHM.md']}
    spec=dict(tasks=manifest,dependency=relative(DEPENDENCY),
        scientific_method='CONVENTIONAL_SURROGATE_FORCE_HMC_FULL_ENDPOINT_MH',
        novelty_claim=False,final_scientific_reviewer='parent Codex',
        control_paths=dict(results=relative(controls),
            source='AUDIT/surrogate_force/check_surrogate_force.py',
            explanation='AUDIT/surrogate_force/README.md'),
        project_protocol_paths=project_paths,source_revision=revision,
        source_files=OWNED_SOURCE+DEPENDENCY_SOURCE,
        interpreter=sys.executable,arms=list(ARMS),chains_per_arm=4,
        warmup=1024,production=4096,seed_rule='2026102341+chain, shared across arms',
        matched_mechanism_comparison='linear versus material only; identical new two-arm protocol',
        ordinary_HMC_external_context=dict(warmup=512,production=1024,lengths=[4,8,12],
            seed='SeedSequence([2026102290,chain])',
            interpretation='external cost/Pareto context, NOT a matched sole-force ablation'),
        superiority_claim='none from unequal lengths/warmup; compare statistical errors and actual total RHS/wall only',
        qualified_reference_guarantee=False,
        audit=dict(function='unchanged p3l.reference_audit.audit_reference',
            warmup=1024,Book_stage='C',analysis_only=True,Target_constructed=False,
            original_MaterialPrior=True,
            thresholds=dict(rhat_max_exclusive=RHAT_LIMIT,
                minimum_material_function_bulk_and_tail_ess=MATERIAL_ESS_MIN,
                full_grid_mean_rms_mcse_max_inclusive=MEAN_RMS_MCSE_MAX),
            qualification='preserve unchanged reference_qualified=False and crosscheck requirement'),
        numerical_events='preserve all saved divergence/large-energy-error flags and transition summaries',
        upper_RHS_per_arm=164384,actual_planned_upper_RHS_both_arms=328768,
        setup_source_reuse='same previously paid A mala-mode initialization as ordinary HMC; each chain rebuilds Full state paid C',
        within_job_setup_reuse='an unchanged paid current endpoint may reuse already paid J/T; ledger remains authoritative',
        ledger_cost_policy='sum actual events by exact C job identity; includes failed/rejected Full8 and J32 setup',
        wall_policy='actual per-job C Book controlled wall plus proxy CPU counters; idle queue wait disclosed separately',
        max_wall_per_chain=7200,max_wall_per_audit=600,
        wait_poll_seconds=30,automatic_retries=False,failed_or_partial_jobs='STOP_PRESERVE_NO_REPLAY',
        already_COMPLETE_jobs='reuse explicit receipts and saved outputs; never run twice',
        GT_read=False,training=False)
    write(protocol,spec)
    for rel in OWNED_SOURCE+DEPENDENCY_SOURCE:
        dest=FOLDER/'SOURCE'/rel;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(ROOT/rel,dest)
    for path in ['configs/PROTOCOL.json','configs/AUTHORIZATION.json','configs/TARGET.json',
        'configs/LOCKED_CODE_VERSION.json','AUDIT/surrogate_force/RESULTS.json',
        'AUDIT/surrogate_force/check_surrogate_force.py','AUDIT/surrogate_force/README.md']:
        dest=FOLDER/'FROZEN_INPUTS'/path;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(OUT/path,dest)
    return spec


def verify_source():
    source_revision()  # Compare immutable estimator/prior dependencies to the locked baseline.
    for rel in OWNED_SOURCE+DEPENDENCY_SOURCE:
        if (ROOT/rel).read_bytes()!=(FOLDER/'SOURCE'/rel).read_bytes():
            raise RuntimeError(f'queued source changed since freeze: {rel}; no automatic execution')


def progress(done,state,**fields):
    write(FOLDER/'PROGRESS.json',dict(done=done,state=state,actual_cost_so_far=summed_cost(done),**fields))


def wait_dependency(done):
    started=time.monotonic();progress(done,'WAITING_FOR_FROZEN_LOCKED_RISK',dependency=relative(DEPENDENCY))
    while not DEPENDENCY.exists():time.sleep(30)
    receipt=json.loads(DEPENDENCY.read_text())
    if receipt.get('status')!='COMPLETE':raise RuntimeError('frozen locked-risk sequence failed; surrogate pilot must not start')
    shutil.copyfile(DEPENDENCY,FOLDER/'LOCKED_RISK_COMPLETION.json')
    return time.monotonic()-started


def wait_paid_lock(done,name):
    started=time.monotonic()
    progress(done,'WAITING_FOR_GLOBAL_PAID_WORKER',next_task=name)
    with (OUT/'paid_worker.lock').open('a') as stream:
        while True:
            try:fcntl.flock(stream,fcntl.LOCK_EX|fcntl.LOCK_NB);break
            except BlockingIOError:time.sleep(30)
        # Child Book performs the authoritative atomic lock/reservation itself.
        # A competing claimant causes a preserved child failure, never a retry.
        fcntl.flock(stream,fcntl.LOCK_UN)
    return time.monotonic()-started


def complete_receipt(task):
    jobdir=OUT/'C/jobs'/task['job'];closed=jobdir/'CLOSED.json'
    if closed.exists():
        receipt=json.loads(closed.read_text())
        if receipt.get('status')!='COMPLETE' or receipt.get('error') is not None:
            raise RuntimeError(f'failed job retained, never replay: {task["job"]}')
        if task['kind']=='paid_chain' and not (jobdir/'CHAIN.npz').is_file():
            raise RuntimeError('completed chain receipt has no saved chain')
        if task['kind']=='analysis_only_audit' and not (OUT/task['destination']/'REFERENCE_STATS.json').is_file():
            raise RuntimeError('completed audit receipt has no saved statistics')
        return receipt
    if jobdir.exists() or (FOLDER/(task['name']+'.log')).exists():
        raise RuntimeError(f'partial/started job retained without replay: {task["job"]}')
    if task['kind']=='analysis_only_audit' and (OUT/task['destination']).exists():
        raise RuntimeError('audit destination exists without completed receipt; no overwrite')
    return None


def numerical_events(paths):
    rows=[]
    for path in paths:
        with np.load(path,allow_pickle=False) as archive:
            divergence=np.asarray(archive['divergences'])
            large=np.asarray(archive['large_energy_error'])
            rows.append(dict(chain=relative(path),
                all_draws=len(divergence),warmup_divergences=int(divergence[:WARMUP].sum()),
                production_divergences=int(divergence[WARMUP:].sum()),
                total_divergences=int(divergence.sum()),
                warmup_large_energy_errors=int(large[:WARMUP].sum()),
                production_large_energy_errors=int(large[WARMUP:].sum()),
                transition_summary=json.loads((path.parent/'TRANSITION_SUMMARY.json').read_text())))
    return dict(rows=rows,total_divergences=sum(r['total_divergences'] for r in rows),
        production_divergences=sum(r['production_divergences'] for r in rows),
        qualification='NUMERICAL_FLAGS_RETAINED_FOR_PARENT_REVIEW; NOT A MIXING_OR_INVARIANCE_CERTIFICATE')


def audit_arm(task):
    arm=task['arm']
    paths=[OUT/f'C/jobs/surrogate_force_{arm}_DEV030_c{c}/CHAIN.npz' for c in range(4)]
    for path in paths:
        closed=json.loads((path.parent/'CLOSED.json').read_text())
        if closed.get('status')!='COMPLETE' or closed.get('error') is not None:
            raise RuntimeError('all four Full-target chain jobs must complete before audit')
    destination=OUT/task['destination']
    if destination.exists():raise FileExistsError('no audit overwrite')
    book=Book('C',task['job'],max_wall=600,analysis_only=True)
    status='ERROR';error=None
    try:
        write(book.jobdir/'PROTOCOL.json',dict(arm=arm,paths=[relative(p) for p in paths],
            audit_function='p3l.reference_audit.audit_reference',warmup=1024,
            sequence_protocol=relative(FOLDER/'PROTOCOL.json'),analysis_only=True,
            Target_constructed=False,GT_read=False,thresholds_unchanged=True,
            scientific_method='CONVENTIONAL_SURROGATE_FORCE_HMC',novelty_claim=False))
        cfg=json.loads((OUT/'configs/TARGET.json').read_text())
        points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
        prior=MaterialPrior(cfg,points)
        stats=audit_reference(paths,prior,1024,destination)
        events=numerical_events(paths)
        write(destination/'NUMERICAL_EVENT_AUDIT.json',events)
        result=dict(arm=arm,audit_status=stats['status'],reference_qualified=stats['reference_qualified'],
            cross_mechanism_status=stats['cross_mechanism_status'],thresholds=stats['thresholds'],
            full_grid_mean_rms_mcse=stats['full_grid_mean_rms_mcse'],
            summaries=stats['summaries'],total_divergences=events['total_divergences'],
            production_divergences=events['production_divergences'],
            interpretation='UNCHANGED_WITHIN_MECHANISM_DIAGNOSTICS_ONLY; PARENT_FINAL_SCIENTIFIC_REVIEW_REQUIRED',
            novelty_claim=False,new_physical_RHS=0,GT_read=False)
        write(book.jobdir/'AUDIT_SUMMARY.json',result)
        status='COMPLETE';return result
    except BaseException as exc:error=repr(exc);raise
    finally:book.close(status,error)


def main():
    FOLDER.mkdir(parents=True,exist_ok=True)
    with (FOLDER/'sequence.lock').open('a') as sequence_lock:
        fcntl.flock(sequence_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if (FOLDER/'CLOSED.json').exists():raise RuntimeError('sequence already closed; no overwrite or replay')
        manifest=tasks();done=[];started=time.monotonic();queue_wait=0.;worker_wait=0.
        try:
            freeze(manifest)
            queue_wait=wait_dependency(done)
            # Existing COMPLETE jobs are reused only after the protected locked
            # risk sequence closes. Failure receipts are never reset or hidden.
            existing_cost=job_costs({task['job'] for task in manifest})
            for task in manifest:
                verify_source()
                receipt=complete_receipt(task)
                if receipt is not None:
                    done.append(dict(name=task['name'],job=task['job'],previously_completed=True,
                        closure=relative(OUT/'C/jobs'/task['job']/'CLOSED.json'),cost=existing_cost[task['job']]))
                    progress(done,'PREVIOUS_COMPLETE_RECEIPT_REUSED');continue
                worker_wait+=wait_paid_lock(done,task['name'])
                verify_source();offset=ledger_offset();result=None
                if task['kind']=='paid_chain':
                    log_path=FOLDER/(task['name']+'.log')
                    with log_path.open('x') as log:
                        child=subprocess.Popen([sys.executable,str(ROOT/'tools'/task['command'][0]),
                            *task['command'][1:]],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
                        try:
                            progress(done,'RUNNING_PAID_CHAIN',task=task['name'],PID=child.pid)
                            print('START',task['name'],child.pid,flush=True)
                            code=child.wait()
                        except BaseException:
                            if child.poll() is None:
                                child.send_signal(signal.SIGINT)
                                try:child.wait(timeout=60)
                                except subprocess.TimeoutExpired:
                                    child.terminate();child.wait(timeout=10)
                            raise
                    if code:raise RuntimeError(f'child failed ({code}); preserve all costs/states, stop: {task["name"]}')
                else:
                    progress(done,'RUNNING_ANALYSIS_ONLY_AUDIT',task=task['name'])
                    result=audit_arm(task)
                complete_receipt(task)
                done.append(dict(name=task['name'],job=task['job'],previously_completed=False,
                    closure=relative(OUT/'C/jobs'/task['job']/'CLOSED.json'),
                    cost=job_costs({task['job']},offset)[task['job']],audit=result))
                progress(done,'TASK_COMPLETE',task=task['name'])
                print('END',task['name'],flush=True)
            write(FOLDER/'CLOSED.json',dict(status='COMPLETE',done=done,
                actual_cost=summed_cost(done),planned_upper_RHS=328768,
                dependency_wait_wall_seconds=queue_wait,worker_wait_wall_seconds=worker_wait,
                sequence_elapsed_wall_seconds=time.monotonic()-started,
                scientific_status='UNQUALIFIED_PENDING_PARENT_SCIENTIFIC_REVIEW',
                matched_mechanism_comparison='linear versus material only',
                ordinary_HMC='external cost/Pareto context, not a matched sole-force ablation',
                method='CONVENTIONAL_SURROGATE_FORCE_HMC',novelty_claim=False,GT_read=False))
        except BaseException as exc:
            # Include failed or interrupted jobs' actual reservations and wall,
            # even when they never entered the successful-task list.
            all_costs=job_costs({task['job'] for task in manifest})
            write(FOLDER/'CLOSED.json',dict(status='ERROR',error=repr(exc),done=done,
                all_planned_job_costs_including_failures=all_costs,
                planned_upper_RHS=328768,dependency_wait_wall_seconds=queue_wait,
                worker_wait_wall_seconds=worker_wait,sequence_elapsed_wall_seconds=time.monotonic()-started,
                failure_policy='STOP_PRESERVE_NO_AUTOMATIC_RETRY',GT_read=False))
            raise


if __name__=='__main__':main()
