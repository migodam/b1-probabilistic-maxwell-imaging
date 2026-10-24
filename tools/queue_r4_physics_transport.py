"""Serial, explicitly released paired physics/supervised refinement controller.

Default --check is read-only/inert. Execution requires --execute and a parent-
reviewed registry with released:true. The individual trainer owns both locks
and its physics Book. Never put this controller inside run_r4_job.py.
"""
from __future__ import annotations

import argparse
import copy
import fcntl
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,write
from p3r4.train_physics_transport import identity,verify_identity,release_template,validate_release
from types import SimpleNamespace

FORMAT='p3r4_paired_physics_transport_campaign_v1'
CAMPAIGN='P3LR4-PHYSICS-TRANSPORT-PAIRED-V1'
SEEDS=[2026102501,2026102502,2026102503]
REGISTRY=OUT/'configs/PHYSICS_TRANSPORT_PAIRED_V1.json'
SOURCES=['src/p3r4/train_physics_transport.py','src/p3r4/density.py','src/p3r4/anchor_density.py',
         'src/p3r4/conditional_coordinates.py','src/p3r4/core.py','src/p3r4/train_density.py',
         'src/p3r4/train_anchor_density.py','tools/queue_r4_physics_transport.py']
GATE=dict(reference='physics seed1 initialization at update0',
    composite_improvement_fraction_minimum=.05,composite_denominator='max(abs(initial_composite),1e-12)',
    selected_meanPhi_over_initial_maximum=.8,proper_NLL_increase_maximum_nats=10.,
    all_conditions_required=True,paired_supervised_complete_required=True,
    selection='trainer-selected fixed validation composite only; no alternate checkpoint or DEV/locked result',
    replicas='if eligible, run BOTH original-source seed2 and seed3 pairs; no seed2 outcome gate',
    interpretation='preregistered validation continuation gate; not posterior qualification or independent holdout')


def load(path):return json.loads(Path(path).read_text())


def build_registry():
    """Create inert metadata; performs no queries, fit, or process launch."""
    definitions=[
        ('raw_spline','density','density_spline_32768','spline',32768,2048),
        ('raw_diag_gaussian','density','extended_diag_gaussian','diag_gaussian',32768,2048),
        ('anchor_spline_none','anchor_density','anchor_spline_none','spline',1024,256)]
    families=[]
    for name,folder,prefix,kind,ntrain,nval in definitions:
        data=OUT/('inputs/strong_context_1024_256' if folder=='anchor_density' else 'inputs/joint')
        train=data/('density_train.npz' if folder=='anchor_density' else 'train32768.npz')
        val=data/('density_validation.npz' if folder=='anchor_density' else 'validation.npz')
        sources=[]
        for seed in SEEDS:
            job=f'{prefix}_s{seed}'
            receipt=OUT/'I/jobs'/job/'CLOSED.json'
            if load(receipt)['status']!='COMPLETE':raise ValueError('source fit not CLOSED COMPLETE: '+job)
            result=OUT/'models'/folder/job/'RESULT.json'
            r=load(result)
            if r['kind']!=kind or r['seed']!=seed or (r['train_size'],r['validation_size'])!=(ntrain,nval):
                raise ValueError('original source checkpoint family/data/seed mismatch: '+job)
            if folder=='anchor_density' and r['anchor_contract']['coordinate_preconditioner']!='none':
                raise ValueError('anchor_none original required')
            sources.append(dict(seed=seed,checkpoint=identity(result.parent/'BEST.pt'),
                original_result=identity(result),original_receipt=identity(receipt),
                parameter_count=r['parameter_count'],original_selected_update=r['selected_update']))
        families.append(dict(id=name,kind=kind,train=identity(train),validation=identity(val),
            train_size=ntrain,validation_size=nval,
            information='raw512+explicit observation-derived anchor66; coordinate none' if folder=='anchor_density' else 'raw512 measurements only',
            sources=sources))
    template=release_template(families[0]['sources'][0]['checkpoint']['path'],families[0]['train']['path'],families[0]['validation']['path'],SEEDS)
    common={k:v for k,v in template.items() if k not in ('checkpoint','train','validation','seeds','released')}
    common['max_wall']=2400
    return dict(format=FORMAT,campaign_id=CAMPAIGN,released=False,
        authorization='parent review and explicit released:true required before --execute',
        common_protocol=common,families=families,replication_gate=GATE,
        seed1_order=[f['id'] for f in families],within_pair_order=['physics','supervised'],
        replication_order='after all seed1 pairs/gates, eligible families in declared order; seed2 pair then seed3 pair',
        maximum_physics_fits=9,maximum_supervised_fits=9,
        maximum_planned_RHS=9*(203264+11264),maximum_hard_RHS=9*(256000+20000),
        max_wall_per_fit_seconds=2400,maximum_fit_wall_seconds=18*2400,
        preflight_gate=identity(OUT/'AUDIT/PHYSICS_TRANSPORT_SAME_INPUT_AUDIT.json'),
        source_identities=[identity(ROOT/path) for path in SOURCES],
        execution_root=str(OUT/'results/physics_transport_paired_v1'),
        output_root=str(OUT/'models/physics_transport/paired_v1'),
        device='mps',fresh_holdout=False,no_locked_GT=True,
        target_rule='unchanged original Full target; no trainer-math or target changes',
        failure_rule='stop controller on failed/open/inconsistent job; never automatically replay or replace',
        resume_rule='explicit --resume may skip verified COMPLETE jobs only; never retry failed/open job',
        cost_rule='each trainer owns standalone physics Book and both locks; controller idle wall not added again to fit wall')


def trainer_release(registry,family,source):
    return dict(copy.deepcopy(registry['common_protocol']),released=True,
        checkpoint=copy.deepcopy(source['checkpoint']),train=copy.deepcopy(family['train']),
        validation=copy.deepcopy(family['validation']),seeds=[source['seed']])


def validate_registry(registry,*,execution=False):
    if registry.get('format')!=FORMAT or registry.get('campaign_id')!=CAMPAIGN:
        raise ValueError('unknown campaign registry')
    if execution and registry.get('released') is not True:
        raise ValueError('campaign inert: parent has not set released:true')
    if registry['replication_gate']!=GATE or registry['within_pair_order']!=['physics','supervised']:
        raise ValueError('frozen gate/pair order changed')
    if registry['seed1_order']!=['raw_spline','raw_diag_gaussian','anchor_spline_none']:
        raise ValueError('frozen family order changed')
    if registry['maximum_physics_fits']!=9 or registry['maximum_supervised_fits']!=9:
        raise ValueError('at most9 paired fits per arm')
    if registry['common_protocol']['max_wall']!=2400 or registry['max_wall_per_fit_seconds']!=2400:
        raise ValueError('2400s hard wall per job required')
    for source in registry['source_identities']:verify_identity(source)
    verify_identity(registry['preflight_gate'])
    preflight=load(registry['preflight_gate']['path'])
    if preflight.get('status')!='PASS_CORRECTED_REFERENCE' or preflight.get('target_unchanged') is not True:
        raise ValueError('same-input actual VJP reference gate not passed')
    if registry['device']!='mps':raise ValueError('frozen campaign device is MPS')
    if registry['common_protocol']['max_rhs']!={'physics':256000,'supervised':20000}:
        raise ValueError('frozen256000/20000 per-fit caps required')
    if registry['maximum_planned_RHS']!=1930752 or registry['maximum_hard_RHS']!=2484000:
        raise ValueError('campaign overall RHS ceiling changed')
    families=registry['families']
    if [f['id'] for f in families]!=registry['seed1_order']:
        raise ValueError('family registry/order mismatch')
    plans=[]
    for family in families:
        verify_identity(family['train']);verify_identity(family['validation'])
        if [s['seed'] for s in family['sources']]!=SEEDS:
            raise ValueError('original matching seeds1/2/3 required')
        for source in family['sources']:
            for k in ('checkpoint','original_result','original_receipt'):verify_identity(source[k])
            if load(source['original_receipt']['path'])['status']!='COMPLETE':
                raise ValueError('original source receipt not complete')
            for mode in registry['within_pair_order']:
                release=trainer_release(registry,family,source)
                plan=validate_release(release,SimpleNamespace(seed=source['seed'],mode=mode))
                plans.append(dict(family=family['id'],seed=source['seed'],mode=mode,
                    conditional_on='seed1 replication gate' if source['seed']!=SEEDS[0] else 'parent release',**plan))
    return plans


def replication_gate(physics_result,supervised_result):
    """Gate at the trainer-selected seed1 checkpoint, never reselect a row."""
    for r in (physics_result,supervised_result):
        if r.get('status')!='COMPLETE' or r.get('completed_updates')!=3000:
            raise ValueError('replication gate requires both complete3000-update arms')
    if physics_result.get('mode')!='physics' or supervised_result.get('mode')!='supervised':
        raise ValueError('paired arm mismatch')
    history=physics_result['history']
    if [h['update'] for h in history]!=list(range(0,3001,300)):
        raise ValueError('fixed update0/every300 roster incomplete')
    initial=history[0]
    selected=next((h for h in history if h['update']==physics_result['selected_update']),None)
    if selected is None:raise ValueError('selected update absent from fixed curve')
    keys=('validation_composite','validation_Phi','validation_full_NLL')
    if not all(math.isfinite(float(row[k])) for row in (initial,selected) for k in keys):
        raise ValueError('nonfinite gate input; no survivor-based selection')
    if selected['validation_composite']!=physics_result['selected_validation_metric']:
        raise ValueError('selected composite/result mismatch')
    improvement=(initial['validation_composite']-selected['validation_composite'])/max(abs(initial['validation_composite']),1e-12)
    conditions=dict(composite_improves_at_least5percent=improvement>=.05,
        meanPhi_at_most80percent_initial=selected['validation_Phi']<=.8*initial['validation_Phi'],
        proper_NLL_increases_at_most10nat=selected['validation_full_NLL']<=initial['validation_full_NLL']+10.)
    return dict(eligible=all(conditions.values()),conditions=conditions,selected_update=selected['update'],
        composite_improvement_fraction=improvement,initial={k:initial[k] for k in keys},
        selected={k:selected[k] for k in keys},replication_seeds=SEEDS[1:] if all(conditions.values()) else [],
        criterion='all3 preregistered validation conditions; physics vs its initialization; paired control retained',
        posterior_qualification='NOT_ESTABLISHED')


def fit_paths(registry,family,source,mode):
    name=f"{registry.get('job_prefix','transport')}_{family['id']}_{mode}_s{source['seed']}"
    return name,Path(registry['output_root'])/name,OUT/'I/jobs'/name


def verify_complete_fit(registry,family,source,mode,release):
    job,out,jobdir=fit_paths(registry,family,source,mode)
    closed=load(jobdir/'CLOSED.json');r=load(out/'RESULT.json')
    if closed['status']!='COMPLETE' or r['status']!='COMPLETE':raise RuntimeError('fit is not COMPLETE: '+job)
    if load(out/'FROZEN_RELEASE.json')!=release:raise RuntimeError('fit used a different release: '+job)
    expected=203264 if mode=='physics' else 11264
    if r['costs']['RHS']!=expected or r['completed_updates']!=3000 or r['mode']!=mode or r['seed']!=source['seed']:
        raise RuntimeError('completed fit violates fixed roster: '+job)
    return r


def wait_for_slots():
    """No held locks while waiting; the child atomically reacquires its locks."""
    while True:
        handles=[]
        try:
            for name in ('training.lock','physics.lock'):
                handle=(OUT/name).open('a');handles.append(handle)
                fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
            return
        except BlockingIOError:
            pass
        finally:
            for handle in handles:
                fcntl.flock(handle,fcntl.LOCK_UN);handle.close()
        time.sleep(15)


def run_campaign(registry,path,*,resume=False):
    validate_registry(registry,execution=True) # Reject before directory/process creation.
    if os.environ.get('R4_JOBDIR'):raise ValueError('do not nest controller in an analysis Book runner')
    registry_identity=identity(path)
    controller_lock=(OUT/'physics_transport_controller.lock').open('a')
    fcntl.flock(controller_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    root=Path(registry['execution_root']);root.mkdir(parents=True,exist_ok=True)
    state_path=root/'CONTROLLER_START.json'
    if state_path.exists() and not resume:raise FileExistsError('controller exists; explicit --resume required; no automatic replay')
    if not state_path.exists():write(state_path,dict(registry=registry_identity,pid=os.getpid(),time_ns=time.time_ns()))
    elif load(state_path)['registry']!=registry_identity:raise ValueError('registry changed after controller start')
    began,cpu=time.monotonic(),time.process_time()
    fits=[];gates={};status,error='COMPLETE',None
    try:
        def fit(family,source,mode):
            verify_identity(registry_identity)
            validate_registry(registry,execution=True)
            release=trainer_release(registry,family,source)
            job,out,jobdir=fit_paths(registry,family,source,mode)
            release_path=root/'releases'/(job+'.json')
            if release_path.exists():
                if load(release_path)!=release:raise RuntimeError('existing per-fit release changed')
            else:write(release_path,release)
            if jobdir.exists() or out.exists():
                if not resume:raise FileExistsError('existing fit requires explicit resume: '+job)
                r=verify_complete_fit(registry,family,source,mode,release)
            else:
                wait_for_slots()
                verify_identity(registry_identity)
                validate_registry(registry,execution=True)
                command=[sys.executable,'-m','p3r4.train_physics_transport','--release',str(release_path),
                    '--out',str(out),'--job',job,'--mode',mode,'--seed',str(source['seed']),'--device',registry['device']]
                env=os.environ.copy();env.update(PYTHONPATH=str(ROOT/'src'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',
                    VECLIB_MAXIMUM_THREADS='1',PYTORCH_MPS_HIGH_WATERMARK_RATIO='0.5',PYTORCH_MPS_LOW_WATERMARK_RATIO='0.4',
                    R4_ARRAY_CODEC='stored',PYTHONUNBUFFERED='1')
                log=root/'logs'/job;log.mkdir(parents=True)
                write(log/'COMMAND.json',dict(argv=command,registry=registry_identity,physical=True))
                with (log/'stdout.log').open('x') as stdout,(log/'stderr.log').open('x') as stderr:
                    child=subprocess.run(command,cwd=ROOT,env=env,stdout=stdout,stderr=stderr)
                if child.returncode:raise RuntimeError(f'{job} failed exit{child.returncode}; preserve receipt/logs; no replay')
                r=verify_complete_fit(registry,family,source,mode,release)
            fits.append(dict(job=job,family=family['id'],seed=source['seed'],mode=mode,
                result=identity(out/'RESULT.json'),receipt=identity(jobdir/'CLOSED.json'),costs=r['costs']))
            write(root/'PROGRESS.json',dict(status='RUNNING',completed_fits=fits,gates=gates))
            return r
        # Every seed1 family is tested in fixed order before any replication.
        for family in registry['families']:
            source=family['sources'][0]
            results={mode:fit(family,source,mode) for mode in registry['within_pair_order']}
            gate=replication_gate(results['physics'],results['supervised'])
            gates[family['id']]=gate
            gate_path=root/'gates'/(family['id']+'.json')
            if gate_path.exists() and load(gate_path)!=gate:raise RuntimeError('existing gate changed')
            write(gate_path,gate)
        for family in registry['families']:
            if gates[family['id']]['eligible']:
                for source in family['sources'][1:]:
                    for mode in registry['within_pair_order']:fit(family,source,mode)
        if sum(f['mode']=='physics' for f in fits)>9 or sum(f['mode']=='supervised' for f in fits)>9:
            raise AssertionError('campaign fit count exceeds frozen ceiling')
    except BaseException as exc:status,error='HOLD_FAILED_NO_AUTORETRY',repr(exc);raise
    finally:
        write(root/'CONTROLLER_RESULT.json',dict(status=status,error=error,completed_fits=fits,gates=gates,
            registry=registry_identity,elapsed_controller_wall_seconds=time.monotonic()-began,
            controller_process_CPU_seconds=time.process_time()-cpu,
            fits_physics_Books_own_paid_cost=True,no_double_billing_controller_wait=True,
            posterior_qualification='NOT_ESTABLISHED'))
        fcntl.flock(controller_lock,fcntl.LOCK_UN);controller_lock.close()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--registry',type=Path,default=REGISTRY)
    p.add_argument('--execute',action='store_true')
    p.add_argument('--resume',action='store_true')
    a=p.parse_args();registry=load(a.registry)
    plans=validate_registry(registry,execution=a.execute)
    if a.execute:run_campaign(registry,a.registry,resume=a.resume)
    else:print(json.dumps(dict(status='READY_BUT_INERT' if not registry['released'] else 'RELEASED_CHECK_ONLY',
        conditional_fit_plans=len(plans),no_process_launch=True,physics_RHS=0,optimizer_updates=0),indent=2))


if __name__=='__main__':main()
