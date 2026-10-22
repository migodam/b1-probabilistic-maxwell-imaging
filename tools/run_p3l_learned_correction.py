"""Fixed exposed-DEV learned independence-MH diagnostic; no training or truth IO.

All five fixed families and all four parents are preregistered together. A
checkpoint that is absent is explicitly NOT_RUN. This script never imports an
offline risk evaluator and does nothing on import.
"""
from pathlib import Path
import argparse
import json
import shutil
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,Target,write,save,from_probit,to_probit
from p3l.learned_proposal import LearnedFullProposal

PARENTS=('DEV_000','DEV_028','DEV_030','DEV_026')
FAMILIES=('mlp_diag','mlp_full','mlp_mdn','mlp_flow','cnn_full')
SEED=2026102270
TRANSITIONS=256
LABEL_COORDINATES='FIXED_GAUSSIAN_WEAK8_ANCHOR_RELATIVE_PLUS_EXACT_PROBIT_PROPERTY2'
JOB='learned_correction_exposed_DEV_s2026102270'


def context_arrays(target,point,y):
    """Exact scalar/image feature law used by run_p3l_contexts.py.

    The current diagnostic explicitly substitutes the saved A anchor for the
    training context's cheap 16-state estimate; this shift is disclosed.
    """
    J=target.jacobian(point)
    singular=np.linalg.svd(J[:,:64],compute_uv=False)
    residual=target.white(point.state.data-y)
    scalar=np.r_[target.white(y),residual,point.q,np.log1p(singular*singular)]
    derivative=target.prior.derivative(point.q)
    image=np.stack((point.state.chi.real,point.state.chi.imag,
        target.prior.A@point.q[:64],
        np.sqrt(np.sum(np.abs(derivative[:,:64])**2,axis=1)))).reshape(4,24,24)
    if scalar.shape!=(1154,) or image.shape!=(4,24,24):
        raise ValueError('frozen scalar1154/image4x24x24 feature contract')
    if not np.isfinite(scalar).all() or not np.isfinite(image).all():
        raise ValueError('context features must be finite')
    return scalar,image,J,residual


def finite_start(parent,anchor_q):
    """Deterministic saved-mode preference, then explicit inward anchor start."""
    for chain in range(4):
        path=OUT/f'A/jobs/reference_{parent}_mala_mode_c{chain}/PROBIT_MODE_INITIALIZATION.json'
        if not path.exists():continue
        document=json.loads(path.read_text())
        rows=[r for r in document['mode_history'] if r['endpoint']==document['mode_endpoint']]
        if len(rows)!=1:raise ValueError(f'ambiguous saved probit mode in {path}')
        x=np.asarray(rows[0]['x'],float)
        if x.shape!=(66,) or not np.isfinite(x).all():
            raise ValueError(f'nonfinite saved probit mode in {path}; no silent fallback')
        return x,dict(kind='SAVED_FINITE_PROBIT_MODE_NOT_POSTERIOR_DRAW',
            source=str(path),source_chain=chain,source_endpoint=document['mode_endpoint'])
    x=to_probit(anchor_q,initialization=True)
    return x,dict(kind='EXPLICIT_INWARD_ANCHOR_BURN_IN_START_NOT_POSTERIOR_DRAW',
        source=str(OUT/f'A/anchors/{parent}.npz'),inward_property_clip=[1e-4,1-1e-4],
        original_native_properties=np.asarray(anchor_q)[64:].tolist(),
        represented_native_properties=from_probit(x)[64:].tolist())


def load_frozen_density(path,family):
    """Restore normalized affine wrapper, never train or infer feature transforms."""
    from p3l.neural_training import load_density_checkpoint
    import torch
    restored=load_density_checkpoint(path,device='cpu')
    payload=restored['checkpoint'];model=restored['model'];meta=restored['metadata']
    expected=dict(scalar_dim=1154,image_shape=(4,24,24),output_dim=10)
    options=payload['model_options']
    if payload['model_type']!=family or int(meta['seed'])!=SEED:
        raise ValueError('frozen family/seed differs from checkpoint')
    for key,value in expected.items():
        actual=tuple(options[key]) if key=='image_shape' else options[key]
        if actual!=value:raise ValueError('checkpoint feature/label dimensions differ')
    if meta['label_coordinates']!=LABEL_COORDINATES:
        raise ValueError('checkpoint raw coordinates must be weak8 plus exact probit2')
    if meta['supervision']!='JOINT_SIMULATION':
        raise ValueError('this fixed comparison requires joint-simulation supervision')
    if meta.get('best_checkpoint_available') is not True:
        raise ValueError('checkpoint is not a completed available validation selection')
    if tuple(model.scalar_mean.shape)!=(1154,) or tuple(model.label_mean.shape)!=(10,):
        raise ValueError('frozen scalar/label affine buffer shape')
    for value in model.state_dict().values():
        if not torch.isfinite(value).all():raise ValueError('nonfinite checkpoint parameter/buffer')
    if torch.any(model.scalar_scale<=0) or torch.any(model.label_scale<=0):
        raise ValueError('positive train-only affine scales required')
    model.eval();model.requires_grad_(False)
    if any(p.device.type!='cpu' for p in model.parameters()):raise ValueError('CPU only')
    return model,dict(format=payload['format'],model_type=payload['model_type'],
        model_options=options,metadata=meta,best_update=payload['best_update'],
        best_validation_raw_target_nll=payload['best_validation_raw_target_nll'],
        online_transforms='checkpoint train-only affine scalar/label buffers; images unchanged',
        density_measure='da8 deta2, then normalized correlated b56 conditional; dz64 deta2')


def costs(book):
    latest=book._latest()
    return dict(sequence=latest['sequence'],**latest['stages']['D'])


def delta(before,after):
    return {key:after[key]-before[key] for key in
        ('RHS','states','LU','controlled_wall_seconds','physics_wall_seconds','GPU_seconds')}


def independent_log_ratio(logpi_proposed,logpi_current,logQ_proposed,logQ_current):
    """Both directions of the same normalized, frozen independence mixture."""
    return float(logpi_proposed-logpi_current+logQ_current-logQ_proposed)


def save_arrays(folder,chain,proposals):
    # These are owned checkpoints for this new job, preserving every repeat.
    save(folder/'CHAIN.npz',**{k:np.asarray(v) for k,v in chain.items()},
         warmup=np.array(0),target_family=np.array('FULL'),
         method=np.array('NORMALIZED_LEARNED_DEFENSIVE_INDEPENDENCE_FULL_MH'),
         probability_status=np.array('INVARIANT_KERNEL_NOT_A_MIXING_CERTIFICATE'))
    save(folder/'PROPOSALS.npz',**{k:np.asarray(v) for k,v in proposals.items()},
         distribution=np.array('FROZEN_LEARNED_FULL66_PROPOSAL_WITH_10_PERCENT_STANDARD_PRIOR'),
         calibrated_posterior=np.array(False))


def run_pair(book,target,parent,family,model,features,start_x,start_meta,index):
    import torch
    scalar,image,J,residual,q0,Vw,Vs=features
    folder=book.jobdir/'comparisons'/parent/family
    folder.mkdir(parents=True,exist_ok=False)
    pair_before=costs(book);started=time.perf_counter()
    # Separate acceptance stream gives identical fixed uniforms across all
    # families even when their proposal samplers consume different RNG draws.
    proposal_rng=np.random.default_rng(np.random.SeedSequence([SEED,index,0]))
    accept_rng=np.random.default_rng(np.random.SeedSequence([SEED,index,1]))
    uniforms=accept_rng.random(TRANSITIONS)
    save(folder/'RANDOMIZATION.npz',MH_uniforms=uniforms,
         proposal_numpy_seed=np.array([SEED,index,0]),
         MH_uniform_seed=np.array([SEED,index,1]),torch_seed=np.array(SEED+index))
    write(folder/'PROTOCOL.json',dict(parent=parent,model=family,training_seed=SEED,
        transitions=TRANSITIONS,initialization=start_meta,initial_x=start_x,
        proposal_numpy_seed=[SEED,index,0],MH_uniform_seed=[SEED,index,1],
        torch_seed=SEED+index,common_MH_uniforms_across_families=True,
        proposal_weight=.9,standard_probit_Gaussian66_prior_weight=.1,
        target='-Phi(from_probit(x))-.5*x@x with respect to dz64 deta2',
        MH_ratio='logpi(xprime)-logpi(x)+logQ(x)-logQ(xprime)',
        proposal_fixed=True,GT_used=False,offline_risk_evaluator_imported=False,
        production_adaptation=False,all_rejections_retained=True,
        support='all R66; all strong coordinates sampled, no clipping',
        no_qualification='256 diagnostic transitions are not a posterior mixing certificate'))
    chain={k:[] for k in ('q','x','Phi','log_density','prediction','endpoint_key',
        'endpoint_path','accepted','log_mh_ratio','logQ_current','log_proposal_forward',
        'log_proposal_reverse','proposal_endpoint_key','RHS_cumulative')}
    raw={k:[] for k in ('q','x','Phi','log_density','prediction','endpoint_key',
        'endpoint_path','logQ','accepted','log_mh_ratio')}
    proposal_material=[]
    with torch.random.fork_rng(devices=[]):
        torch.manual_seed(SEED+index)
        proposal=LearnedFullProposal(model,scalar,image,Vw,Vs,q0,J,residual,defensive_weight=.1)
        # A new paid Full forward is required even for a previously saved mode.
        current_x=start_x.copy();current=target.evaluate(from_probit(current_x),
            f'LEARNED_FULL_MH_INITIAL_{parent}_{family}')
        current_phi=target.phi(current,target.observation(parent))
        current_logpi=-current_phi-.5*float(current_x@current_x)
        current_logQ=float(proposal.log_prob(current_x))
        if not np.isfinite([current_logpi,current_logQ]).all():
            raise FloatingPointError('finite initial full density and normalized proposal required')
        save(folder/'INITIAL.npz',x=current_x,q=current.q,Phi=np.array(current_phi),
            log_density=np.array(current_logpi),logQ=np.array(current_logQ),
            endpoint_key=np.array(current.key),prediction=current.state.data)
        target.persist_cache(current,folder/'cache/CURRENT.npz')
        y=target.observation(parent)
        for i in range(TRANSITIONS):
            before=costs(book)
            candidate=proposal.sample(1,rng=proposal_rng)[0]
            logQ_proposed=float(proposal.log_prob(candidate))
            if not np.isfinite(logQ_proposed):
                raise FloatingPointError('defensive proposal density must be finite on a finite draw')
            # Persist the proposal and MH randomness before its paid query.
            save(folder/'intents'/f'{i:04d}.npz',x=candidate,q=from_probit(candidate),
                logQ=np.array(logQ_proposed),MH_uniform=np.array(uniforms[i]))
            proposed=target.evaluate(from_probit(candidate),f'LEARNED_FULL_MH_{parent}_{family}_{i:04d}')
            phi=target.phi(proposed,y);lp=-phi-.5*float(candidate@candidate)
            endpoint=str(book.jobdir/'endpoints'/f'{proposed.key}.npz')
            after_full=costs(book)
            # Persist Full output before any later density/acceptance action.
            save(folder/'raw_proposals'/f'{i:04d}.npz',x=candidate,q=proposed.q,
                Phi=np.array(phi),log_density=np.array(lp),logQ=np.array(logQ_proposed),
                prediction=proposed.state.data,chi=proposed.state.chi,
                endpoint_key=np.array(proposed.key),endpoint_path=np.array(endpoint))
            write(folder/'costs'/f'{i:04d}.json',dict(iteration=i,
                ledger_sequence_before=before['sequence'],ledger_sequence_after=after_full['sequence'],
                delta=delta(before,after_full),proposal_endpoint=endpoint,
                status='FULL_PROPOSAL_SAVED_BEFORE_MH_DECISION'))
            # Re-evaluate reverse density arithmetically to audit both explicit
            # normalized directions; no physics or changing proposal involved.
            reverse=float(proposal.log_prob(current_x))
            if not np.isclose(reverse,current_logQ,rtol=0,atol=1e-10):
                raise RuntimeError('frozen proposal density changed at current state')
            ratio=independent_log_ratio(lp,current_logpi,logQ_proposed,reverse)
            if np.isnan(ratio) or np.isposinf(lp):raise FloatingPointError('invalid Full MH energy ratio')
            accepted=bool(np.log(max(uniforms[i],np.finfo(float).tiny))<min(0.,ratio))
            if accepted:
                current_x=candidate.copy();current=proposed;current_phi=phi
                current_logpi=lp;current_logQ=logQ_proposed
            current_path=str(book.jobdir/'endpoints'/f'{current.key}.npz')
            entries=dict(q=current.q.copy(),x=current_x.copy(),Phi=current_phi,
                log_density=current_logpi,prediction=current.state.data.copy(),endpoint_key=current.key,
                endpoint_path=current_path,accepted=accepted,log_mh_ratio=ratio,
                logQ_current=current_logQ,log_proposal_forward=logQ_proposed,
                log_proposal_reverse=reverse,proposal_endpoint_key=proposed.key,
                RHS_cumulative=after_full['RHS']-pair_before['RHS'])
            for key,value in entries.items():chain[key].append(value)
            entries_raw=dict(q=proposed.q.copy(),x=candidate.copy(),Phi=phi,log_density=lp,
                prediction=proposed.state.data.copy(),endpoint_key=proposed.key,endpoint_path=endpoint,
                logQ=logQ_proposed,accepted=accepted,log_mh_ratio=ratio)
            for key,value in entries_raw.items():raw[key].append(value)
            proposal_material.append(proposed.state.chi.copy())
            save(folder/'transitions'/f'{i:04d}.npz',**entries,
                MH_uniform=np.array(uniforms[i]),proposal_x=candidate)
            save_arrays(folder,chain,raw)
            # Retain exactly the paid current cache, including on rejections.
            target.persist_cache(current,folder/'cache/CURRENT.npz')
            write(folder/'CHECKPOINT.json',dict(iteration=i+1,total=TRANSITIONS,
                accepted=int(sum(chain['accepted'])),current_endpoint=current_path,
                current_Phi=current_phi,numpy_proposal_RNG=proposal_rng.bit_generator.state,
                torch_RNG_state_saved_separately=True,
                cost=delta(pair_before,costs(book)),GT_used=False))
            torch.save(torch.get_rng_state(),folder/'TORCH_RNG_STATE.pt')
            if (i+1)%64==0:
                print(json.dumps(dict(parent=parent,model=family,iteration=i+1,
                    acceptance=float(np.mean(chain['accepted'])),Phi=current_phi)),flush=True)
        material=np.asarray(proposal_material)
        save(folder/'MODEL_PROPOSAL_MATERIAL_MOMENTS.npz',mean=material.mean(axis=0),
            variance=np.mean(np.abs(material-material.mean(axis=0))**2,axis=0),
            draws=np.array(TRANSITIONS),
            status=np.array('EMPIRICAL_UNWEIGHTED_DEFENSIVE_MODEL_PROPOSAL_DIAGNOSTIC'),
            calibrated_posterior=np.array(False))
    summary=dict(parent=parent,model=family,status='COMPLETE_EXPOSED_FULL_MH_DIAGNOSTIC',
        transitions=TRANSITIONS,acceptance=float(np.mean(chain['accepted'])),
        accepted=int(sum(chain['accepted'])),unique_chain_endpoints=len(set(chain['endpoint_key'])),
        cost=delta(pair_before,costs(book)),wall_seconds=time.perf_counter()-started,
        probability_status='NORMALIZED_FULL_TARGET_INVARIANT_KERNEL_NOT_A_MIXING_CERTIFICATE',
        proposal_moments='UNWEIGHTED_MODEL_PROPOSAL_DIAGNOSTICS_NOT_CALIBRATED_POSTERIOR',
        common_uniforms_across_models=True,all_rejections_retained=True,GT_used=False)
    write(folder/'SUMMARY.json',summary)
    del current,proposed,proposal
    target.backend()._transfers.clear()
    return summary


def run(max_wall=7200):
    import torch
    if max_wall<1 or max_wall>7200:raise ValueError('bounded max_wall must be within 1..7200 seconds')
    frame_path=OUT/'D/FIXED_FRAME.npz'
    with np.load(frame_path,allow_pickle=False) as frame:
        Vw=frame['Vw'].copy();Vs=frame['Vs'].copy()
    if Vw.shape!=(64,8) or Vs.shape!=(64,56) or not np.allclose(
            np.column_stack((Vw,Vs)).T@np.column_stack((Vw,Vs)),np.eye(64),atol=1e-10,rtol=0):
        raise ValueError('fixed complete orthogonal weak8/strong56 frame required')
    book=Book('D',JOB,max_wall=max_wall);status='ERROR';error=None
    previous_threads=torch.get_num_threads();torch.set_num_threads(1)
    summaries=[]
    try:
        if (book.jobdir/'PROTOCOL.json').exists():raise FileExistsError('no overwrite or automatic replay')
        checkpoint_sources={family:OUT/f'D/models/{family}/seed{SEED}/BEST_DENSITY.pt' for family in FAMILIES}
        write(book.jobdir/'PROTOCOL.json',dict(parents=PARENTS,families=FAMILIES,
            training_seed=SEED,transitions_per_pair=TRANSITIONS,device='cpu',threads=1,
            normalized_proposal=True,defensive_standard_prior_weight=.1,
            family_or_parent_selection='FIXED_EXPOSED_DEV_SET; NO_GT_SELECTION',
            checkpoint_sources=checkpoint_sources,missing_checkpoint='EXPLICIT_NOT_RUN_NO_TRAINING',
            max_wall_seconds=max_wall,expected_full_forward_states=4*5*(TRANSITIONS+1),
            conservative_RHS_bound_if_all_run=4*5*(TRANSITIONS+1)*8+4*32,
            online_context='SAVED_A_ANCHOR; exact run_p3l_contexts scalar/image feature definitions',
            context_distribution_shift='training used cheap16-state contexts; this uses paid A anchors',
            feature_order=dict(scalar=['white(y)','white(F(anchor)-y)','native_anchor_q66','log1p(sv(Jz)^2)'],
                image=['chi.real','chi.imag','A@z','sqrt(sum(abs(Dchi_dz)^2))']),
            scalar_dim=1154,image_shape=[4,24,24],label_coordinates=LABEL_COORDINATES,
            target='UNCHANGED_FULL_PROBIT_GAUSSIAN66_POSTERIOR',
            budget='D only; global Book one-worker lock; no other-stage transfer',
            GT_read=False,training=False,offline_risk_evaluator=False,
            scientific_status='EXPOSED_DIAGNOSTIC_NOT_REFERENCE_OR_CALIBRATED_POSTERIOR'))
        for rel in ['tools/run_p3l_learned_correction.py','tools/run_p3l_contexts.py',
            'src/p3l/learned_proposal.py','src/p3l/neural_training.py','src/p3l/neural.py',
            'src/p3l/core.py','src/p3l/moments.py','src/p1_material/prior.py']:
            dest=book.jobdir/'SOURCE'/rel;dest.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/rel,dest)
        shutil.copyfile(frame_path,book.jobdir/'FIXED_FRAME_FROZEN.npz')
        shutil.copyfile(OUT/'configs/TARGET.json',book.jobdir/'TARGET_FROZEN.json')
        models={}
        for family,source in checkpoint_sources.items():
            if not source.is_file():
                for parent in PARENTS:summaries.append(dict(parent=parent,model=family,
                    status='NOT_RUN_CHECKPOINT_ABSENT',checkpoint=str(source),new_physics_RHS=0))
                continue
            destination=book.jobdir/'CHECKPOINTS'/family/'BEST_DENSITY.pt'
            destination.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,destination)
            try:
                model,metadata=load_frozen_density(destination,family)
                models[family]=model
                write(destination.with_name('VALIDATED_METADATA.json'),metadata)
            except Exception as exc:
                write(destination.with_name('VALIDATION_FAILURE.json'),dict(error=repr(exc),status='NOT_RUN_INVALID_CHECKPOINT'))
                for parent in PARENTS:summaries.append(dict(parent=parent,model=family,
                    status='NOT_RUN_INVALID_CHECKPOINT',checkpoint=str(source),error=repr(exc),new_physics_RHS=0))
        write(book.jobdir/'CHECKPOINT_PREFLIGHT.json',dict(available_families=list(models),not_run=summaries))
        if models:
            target=Target(book)
            for index,parent in enumerate(PARENTS):
                y=target.observation(parent)
                anchor=target.load_cache(OUT/f'A/anchors/{parent}.npz',owner='saved paid y-only A anchor for learned correction online features')
                scalar,image,J,residual=context_arrays(target,anchor,y)
                start_x,start_meta=finite_start(parent,anchor.q)
                save(book.jobdir/'CONTEXTS'/f'{parent}.npz',scalar=scalar,image=image,
                    q=anchor.q,J=J,residual=residual,initial_x=start_x)
                write(book.jobdir/'CONTEXTS'/f'{parent}.json',dict(parent=parent,
                    source=str(OUT/f'A/anchors/{parent}.npz'),initialization=start_meta,
                    GT_used=False,context_distribution_shift='TRAIN_CHEAP16_TO_A_ANCHOR_DISCLOSED'))
                features=(scalar,image,J,residual,anchor.q.copy(),Vw,Vs)
                del anchor
                for family in FAMILIES:
                    if family not in models:continue
                    try:summary=run_pair(book,target,parent,family,models[family],features,start_x,start_meta,index)
                    except BaseException as exc:
                        write(book.jobdir/'PAIR_FAILURE.json',dict(parent=parent,model=family,error=repr(exc),
                            status='STOP_PRESERVE_ALL_PAID_ENDPOINTS_AND_LAST_CHECKPOINT',GT_used=False))
                        raise
                    summaries.append(summary)
                    write(book.jobdir/'PROGRESS.json',dict(comparisons=summaries))
                target.backend()._transfers.clear()
        write(book.jobdir/'SUMMARY.json',dict(status='COMPLETE_FIXED_EXPOSED_COMPARISON' if models else 'NOT_RUN_NO_VALID_CHECKPOINTS',
            comparisons=summaries,GT_read=False,training=False,device='cpu',
            posterior_certification=False,normalization='architecture-and-affine-exact; not numerical mixing evidence'))
        status='COMPLETE'
    except BaseException as exc:error=repr(exc);raise
    finally:
        torch.set_num_threads(previous_threads)
        book.close(status,error)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-wall',type=int,default=7200)
    args=parser.parse_args();run(args.max_wall)


if __name__=='__main__':main()
