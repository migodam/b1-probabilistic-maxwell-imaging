"""Optional frozen locked Full-target importance correction; never reads GT.

No action occurs on import. All128 online imaging parents must already have
completed under a COMPLETE closed job. The parent must explicitly enable
METHODS_FROZEN_FOR_LOCKED.json:locked_importance with enabled=true,
parents=[original128], draws_per_parent=128, sampling_seed=2026102331,
selected_model=<exact locked_imaging.selected_model>, defensive_weight=.1,
and max_wall_seconds in (0,7200]. No freeze is created or altered here.

Every draw is a fresh paid Full forward; no gradients, new context fits,
cached-physics-free evaluations, adaptive early stopping, or parent filtering.
Finite SNIS is biased and remains unqualified regardless of ESS.
"""
from pathlib import Path
from contextlib import nullcontext
import gc
import json
import shutil
import sys
import time
import numpy as np
from scipy.special import logsumexp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,Target,BudgetStop,write,save,from_probit
from p3l.learned_images import (validation_selection,load_selected_density,
    context_features,make_proposal)

JOB='locked_importance_frozen_128'
PARENTS=[f'SBI_locked_{i:04d}' for i in range(128)]
COUNT=128
SEED=2026102331
DEST=OUT/'D/locked_importance'
SUMMARY=OUT/'D/LOCKED_FULL_TARGET_IMPORTANCE.json'
QUALIFICATION='FINITE_SNIS_BIASED_UNQUALIFIED_POSTERIOR_NO_ESS_ONLY_CERTIFICATION'


def preflight():
    frozen_path=OUT/'METHODS_FROZEN_FOR_LOCKED.json'
    frozen=json.loads(frozen_path.read_text())
    spec=frozen.get('locked_importance',{})
    required=dict(enabled=True,parents=PARENTS,draws_per_parent=COUNT,
        sampling_seed=SEED,defensive_weight=.1)
    for name,value in required.items():
        if spec.get(name)!=value:raise RuntimeError(f'locked importance not explicitly frozen: {name}')
    if not isinstance(spec.get('max_wall_seconds'),(int,float)) or not 0<spec['max_wall_seconds']<=7200:
        raise RuntimeError('explicit bounded locked-importance wall cap required')
    selected=spec.get('selected_model')
    if not isinstance(selected,dict) or selected!=frozen.get('locked_imaging',{}).get('selected_model'):
        raise RuntimeError('importance and online imaging must use the identical frozen selected model')
    manifest=json.loads((OUT/'configs/PARENTS.json').read_text())
    if [r['parent'] for r in manifest if r['group']=='locked']!=PARENTS:
        raise RuntimeError('original frozen128 locked parent identities/order differ')
    closure_path=OUT/'D/jobs/locked_imaging_frozen_panel/CLOSED.json'
    closure=json.loads(closure_path.read_text())
    if closure.get('status')!='COMPLETE' or closure.get('error') is not None:
        raise RuntimeError('all locked online imaging must finish and close before importance')
    online_path=OUT/'D/LOCKED_ONLINE_IMAGING.json'
    online=json.loads(online_path.read_text())
    rows=online.get('rows',[])
    if (len(rows)!=128 or [r.get('parent') for r in rows]!=PARENTS or
            online.get('all128_retained') is not True or online.get('GT_read') is not False):
        raise RuntimeError('complete original128 truth-free online results required')
    completed={'ONLINE_ESTIMATES_RECORDED_FAILURES_RETAINED','PARENT_FAILURE_RETAINED'}
    if any(r.get('status') not in completed for r in rows):
        raise RuntimeError('a locked online parent remains unattempted or interrupted')
    # All15+simple selection depends only on completed train/validation reports.
    selection=validation_selection(OUT)
    if selection['selected']!=selected:
        raise RuntimeError('selected model changed since the locked method freeze')
    if DEST.exists() or SUMMARY.exists():
        raise FileExistsError('locked importance output exists; no overwrite or implicit resume')
    frame_path=OUT/'D/FIXED_FRAME.npz'
    with np.load(frame_path,allow_pickle=False) as data:
        frame={name:data[name].copy() for name in ('Vw','Vs')}
    F=np.column_stack((frame['Vw'],frame['Vs']))
    metadata=json.loads(frame_path.with_suffix('.json').read_text())
    if (frame['Vw'].shape!=(64,8) or frame['Vs'].shape!=(64,56) or
            not np.isfinite(F).all() or not np.allclose(F.T@F,np.eye(64),atol=1e-10,rtol=0) or
            metadata.get('training_only') is not True or metadata.get('GT_used') is not False):
        raise RuntimeError('frozen complete orthogonal training-only frame required')
    return spec,selected,frame,dict(freeze=frozen_path,online_closure=closure_path,online_results=online_path)


def ledger(book):
    record=book._latest()
    return dict(sequence=record['sequence'],**record['stages']['D'])


def difference(before,after):
    return {k:after[k]-before[k] for k in
        ('RHS','states','LU','controlled_wall_seconds','physics_wall_seconds','GPU_seconds')}


def freeze_sources(book,spec,selected,sources):
    if (book.jobdir/'PROTOCOL.json').exists():raise FileExistsError('job protocol exists; no replay')
    write(book.jobdir/'PROTOCOL.json',dict(locked_importance=spec,
        parents=PARENTS,draws_per_parent=128,sampling_seed=SEED,
        proposal_numpy_seed='[2026102331,parent_index,0]',
        numpy_density_seed='[2026102331,parent_index,1]',torch_seed='2026102331+parent_index',
        selected_model=selected,mixture_standard_probit_prior_weight=.1,
        target_log_kernel='-Phi(from_probit(x))-.5*x@x-33*log(2*pi)',
        likelihood_normalization='original exp(-Phi) kernel; observation-only constants omitted',
        logweights='log_target-logQ, exact normalized full66 defensive mixture logQ',
        forward_RHS_per_draw=8,per_parent_RHS=1024,total_fresh_states=16384,total_RHS_upper=131072,
        context='saved CHEAP16_CONTEXT from completed locked online imaging; no new fitting/J/cache physics',
        context_setup_cost='previous online D cost shared and disclosed; not hidden or billed again',
        every_proposal='fresh Target.evaluate Full state billed D; no cached-state likelihood reuse',
        stopped_early='budget or interruption only; never ESS/logweight/risk-based pruning',
        parent_or_draw_failures='all recorded; no replacement or filtering; incomplete parents have no final ESTIMATE',
        all128_parents_retained=True,GT_read=False,training=False,gradients=False,
        qualification=QUALIFICATION,
        output_contract=dict(path='D/locked_importance/{parent}/ESTIMATE.npz',
            keys=['weighted_material_mean','x','q','Phi','log_target','logQ','logweights','weights','chi','predictions'],
            companion='SUMMARY.json with ESS,max_weight,status'),
        global_paid_worker='shared Book one-worker lock',
        source_receipts={name:path.relative_to(OUT).as_posix() for name,path in sources.items()}))
    for rel in ['tools/run_p3l_locked_importance.py','src/p3l/learned_images.py',
        'src/p3l/learned_proposal.py','src/p3l/simple_density.py','src/p3l/neural.py',
        'src/p3l/moments.py','src/p3l/core.py','src/p1_material/prior.py']:
        dst=book.jobdir/'SOURCE'/rel;dst.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(ROOT/rel,dst)
    for name,path in sources.items():shutil.copyfile(path,book.jobdir/f'{name.upper()}_FROZEN.json')
    shutil.copyfile(OUT/'D/FIXED_FRAME.npz',book.jobdir/'FIXED_FRAME_FROZEN.npz')
    shutil.copyfile(OUT/'D/FIXED_FRAME.json',book.jobdir/'FIXED_FRAME_FROZEN.json')
    shutil.copyfile(OUT/'configs/TARGET.json',book.jobdir/'TARGET_FROZEN.json')
    shutil.copyfile(ROOT/'research/p1_material_target_20261018/inputs/geometry.json',book.jobdir/'GEOMETRY_FROZEN.json')
    logical_checkpoint=Path(selected['checkpoint_path'])
    if logical_checkpoint.is_absolute() or '..' in logical_checkpoint.parts:
        raise ValueError('selected checkpoint must be portable relative to campaign OUT')
    checkpoint=(OUT/logical_checkpoint).resolve()
    if not checkpoint.is_relative_to(OUT.resolve()):
        raise ValueError('selected checkpoint escaped campaign OUT')
    destination=book.jobdir/('MODEL_FROZEN'+checkpoint.suffix)
    shutil.copyfile(checkpoint,destination)
    # The strict loader validates the original relative path against completed
    # fit/NN identity. The frozen copy is retained as provenance; completed
    # selected checkpoints are immutable and are not retrained by this runner.
    loaded=load_selected_density(selected,out=OUT,rng=np.random.default_rng(SEED))
    if (loaded['model_type']!=selected['model_type'] or
            loaded['metadata'].get('seed')!=selected['seed'] or
            loaded['best_validation_raw_target_nll']!=selected['best_validation_raw_target_nll']):
        raise RuntimeError('restored normalized model differs from frozen selection')
    write(book.jobdir/'MODEL_RESTORATION.json',dict(original=selected,
        original_checkpoint_relative=logical_checkpoint.as_posix(),
        frozen_copy_checkpoint_relative=destination.relative_to(OUT).as_posix(),
        restoration_source='VALIDATED_IMMUTABLE_ORIGINAL_CHECKPOINT_RELATIVE_TO_OUT',
        frozen_copy_role='retained execution snapshot, not passed as a fabricated selected-fit path',
        restored_type=loaded['model_type'],metadata=loaded['metadata'],CPU=True))
    return loaded['model']


def context(target,parent,y,folder):
    source=OUT/f'D/locked_imaging/{parent}/CHEAP16_CONTEXT.npz'
    with np.load(source,allow_pickle=False) as archive:
        a={k:archive[k].copy() for k in ('scalar','image','q','prediction','J','residual')}
    shapes=dict(scalar=(1154,),image=(4,24,24),q=(66,),prediction=(32,8),J=(512,66),residual=(512,))
    for name,shape in shapes.items():
        if a[name].shape!=shape or not np.isfinite(a[name]).all():
            raise ValueError(f'invalid saved cheap16 context {name}')
    expected_scalar,expected_image=context_features(target.prior,a['q'],a['prediction'],a['J'],y,target.sigma)
    if (not np.allclose(a['scalar'],expected_scalar,rtol=1e-12,atol=1e-12) or
            not np.allclose(a['image'],expected_image,rtol=1e-12,atol=1e-12) or
            not np.allclose(a['residual'],target.white(a['prediction']-y),rtol=1e-12,atol=1e-12)):
        raise ValueError('saved cheap16 context differs from exact training feature API/observation')
    save(folder/'CONTEXT_FROZEN.npz',**a,observation=y)
    write(folder/'CONTEXT_SOURCE.json',dict(source=source.relative_to(OUT).as_posix(),
        setup_owner='completed D locked_imaging_frozen_panel',new_context_physics_RHS=0,
        previous_online_setup_cost_must_be_disclosed=True,GT_read=False))
    return a


def checkpoint_samples(folder,samples):
    if not samples:return None
    keys=('draw_index','x','q','Phi','log_target','logQ','logweights','chi','predictions','endpoint_key','endpoint_path')
    arrays={k:np.asarray([s[k] for s in samples]) for k in keys}
    save(folder/'SAMPLES.npz',**arrays,requested_draws=np.array(COUNT),
        qualification=np.array(QUALIFICATION))
    return arrays


def complete_estimate(folder,samples):
    arrays=checkpoint_samples(folder,samples)
    if arrays is None or len(samples)!=COUNT:raise ValueError('all128 completed draw slots required')
    lw=arrays['logweights'];logsum=float(logsumexp(lw))
    if not np.isfinite(logsum):raise FloatingPointError('no finite importance normalization; raw draws retained')
    # Normalize from centered weights directly: adding log(sum) to a huge
    # common log-energy offset can otherwise lose the normalization in float64.
    centered=np.exp(lw-np.max(lw));weights=centered/np.sum(centered)
    mean=np.einsum('i,ij->j',weights,arrays['chi'])
    second=np.einsum('i,ij->j',weights,np.abs(arrays['chi'])**2)
    variance=np.einsum('i,ij->j',weights,np.abs(arrays['chi']-mean)**2)
    ESS=float(1./np.sum(weights*weights));maximum=float(np.max(weights))
    save(folder/'ESTIMATE.npz',**arrays,weights=weights,weighted_material_mean=mean,
        weighted_material_variance=variance,weighted_material_second_abs_moment=second,
        ESS=np.array(ESS),max_weight=np.array(maximum),
        logmeanweight=np.array(logsum-np.log(COUNT)),
        qualification=np.array(QUALIFICATION),posterior_qualified=np.array(False))
    return dict(ESS=ESS,max_weight=maximum,logmeanweight=logsum-np.log(COUNT),
        status='COMPLETE_FINITE_SNIS_UNQUALIFIED_POSTERIOR',draws=COUNT,
        finite_SNIS_biased=True,ESS_alone_never_qualifies=True,posterior_qualified=False)


def run_parent(book,target,model,frame,parent,index,row):
    folder=DEST/parent;folder.mkdir(parents=True,exist_ok=False)
    before=ledger(book);started=time.perf_counter();samples=[];failures=[]
    proposal=None;point=None
    try:
        y=target.observation(parent);a=context(target,parent,y,folder)
        if hasattr(model,'rng'):
            model.rng=np.random.default_rng(np.random.SeedSequence([SEED,index,1]))
        proposal=make_proposal(model,a['scalar'],a['image'],frame,a['q'],a['J'],a['residual'])
        if proposal.defensive_weight!=.1:raise RuntimeError('exact frozen10% defensive mixture required')
        rng=np.random.default_rng(np.random.SeedSequence([SEED,index,0]))
        stream=nullcontext()
        if proposal._torch is not None:
            stream=proposal._torch.random.fork_rng(devices=[])
        with stream:
            if proposal._torch is not None:proposal._torch.manual_seed(SEED+index)
            for i in range(COUNT):
                draw_before=ledger(book);point=None
                try:
                    x=proposal.sample(1,rng=rng)[0];q=from_probit(x)
                    logQ=float(proposal.log_prob(x))
                    if not np.isfinite(logQ):raise FloatingPointError('finite normalized defensive mixture logQ required')
                    save(folder/'intents'/f'{i:04d}.npz',x=x,q=q,logQ=np.array(logQ),draw_index=np.array(i))
                    point=target.evaluate(q,f'LOCKED_FULL_IMPORTANCE_{parent}_{i:04d}')
                    Phi=target.phi(point,y)
                    log_target=-Phi-.5*float(x@x)-33*np.log(2*np.pi)
                    logweight=log_target-logQ
                    endpoint=(book.jobdir/'endpoints'/f'{point.key}.npz').relative_to(OUT).as_posix()
                    sample=dict(draw_index=i,x=x,q=q,Phi=Phi,log_target=log_target,logQ=logQ,
                        logweights=logweight,chi=point.state.chi.copy(),predictions=point.state.data.copy(),
                        endpoint_key=point.key,endpoint_path=endpoint)
                    # Every paid result is saved before aggregate arithmetic.
                    save(folder/'draws'/f'{i:04d}.npz',**sample,observation=y)
                    if not np.isfinite([Phi,log_target,logweight]).all():
                        raise FloatingPointError('nonfinite Full importance energy; paid output retained')
                    samples.append(sample);checkpoint_samples(folder,samples)
                except BudgetStop as exc:
                    failure=dict(draw_index=i,status='BUDGET_STOP_RETAINED',error=repr(exc))
                    failures.append(failure);write(folder/'failures'/f'{i:04d}.json',failure)
                    raise
                except Exception as exc:
                    failure=dict(draw_index=i,status='FAILED_FIXED_DRAW_RETAINED_NO_REPLACEMENT',error=repr(exc))
                    failures.append(failure);write(folder/'failures'/f'{i:04d}.json',failure)
                finally:
                    after=ledger(book)
                    write(folder/'costs'/f'{i:04d}.json',dict(draw_index=i,
                        ledger_sequence_before=draw_before['sequence'],ledger_sequence_after=after['sequence'],
                        cost=difference(draw_before,after)))
                    point=None
                write(folder/'CHECKPOINT.json',dict(attempted=i+1,requested=COUNT,completed=len(samples),
                    failures=failures,numpy_proposal_RNG=rng.bit_generator.state,
                    numpy_model_RNG=model.rng.bit_generator.state if hasattr(model,'rng') else None,
                    qualification=QUALIFICATION,cost=difference(before,ledger(book))))
                if proposal._torch is not None:
                    proposal._torch.save(proposal._torch.get_rng_state(),folder/'TORCH_RNG_STATE.pt')
        if failures:
            row.update(status='INCOMPLETE_FIXED_DRAW_FAILURES_RETAINED_NO_FINAL_ESTIMATE',
                completed=len(samples),failures=failures)
        else:row.update(complete_estimate(folder,samples))
    except BudgetStop:
        row.update(status='BUDGET_STOP_PARTIAL_DRAWS_RETAINED',completed=len(samples),failures=failures)
        raise
    except Exception as exc:
        row.update(status='PARENT_FAILURE_RETAINED_NO_FINAL_ESTIMATE',completed=len(samples),
            failures=failures,error=repr(exc))
    finally:
        point=proposal=None
        row.update(parent=parent,requested_draws=COUNT,qualification=QUALIFICATION,
            GT_read=False,cost=difference(before,ledger(book)),wall_seconds=time.perf_counter()-started)
        write(folder/'SUMMARY.json',row)
        if target.model is not None:target.model._transfers.clear()
        gc.collect()


def main():
    spec,selected,frame,sources=preflight()
    book=Book('D',JOB,max_wall=spec['max_wall_seconds']);status='ERROR';error=None
    rows=[dict(parent=p,status='NOT_RUN',requested_draws=COUNT) for p in PARENTS]
    try:
        model=freeze_sources(book,spec,selected,sources)
        DEST.mkdir(parents=True,exist_ok=False)
        target=Target(book)
        for index,row in enumerate(rows):
            try:run_parent(book,target,model,frame,row['parent'],index,row)
            finally:
                write(SUMMARY,dict(rows=rows,all128_retained=True,GT_read=False,
                    selected_model=selected,draws_per_parent=COUNT,seed=SEED,
                    total_RHS_upper=131072,qualification=QUALIFICATION,
                    no_parent_filtering=True,no_adaptive_early_stop=True,
                    context_setup='previous completed online-imaging D cost disclosed separately'))
            print(json.dumps(dict(parent=row['parent'],status=row['status'],
                ESS=row.get('ESS'),max_weight=row.get('max_weight'))),flush=True)
        status='COMPLETE'
    except BaseException as exc:error=repr(exc);raise
    finally:book.close(status,error)


if __name__=='__main__':main()
