"""Fitted simple Gaussian Full-MH comparison; no fitting, GT, or risk IO.

Use the two genuine train/validation-fitted controls already selected in
D/simple_density/RESULTS.json. Contexts are fresh matched coarse_TRF16 from
the original prior start, shared by both controls. No action on import.
"""
from pathlib import Path
import argparse
import json
import shutil
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,Target,write,save
from p3l.inversion import reconstruct
from p3l.simple_density import SimpleConditionalGaussian
from run_p3l_learned_correction import (PARENTS,SEED,TRANSITIONS,context_arrays,
    finite_start,costs,delta,run_pair)

JOB='simple_correction_exposed_DEV'
FAMILIES=('selected','unconditional_selected')
KINDS=dict(selected='FITTED_PCA_CONDITIONAL_GAUSSIAN',
    unconditional_selected='FITTED_UNCONDITIONAL_GAUSSIAN')


def selected_sources():
    path=OUT/'D/simple_density/RESULTS.json'
    result=json.loads(path.read_text())
    if result.get('status')!='FITTED_NORMALIZED_JOINT_SUPERVISION_CONTROLS':
        raise ValueError('completed fitted simple-density selection required')
    if result.get('train_parents')!=512 or result.get('validation_parents')!=128:
        raise ValueError('frozen train512/validation128 fitted controls required')
    if result.get('posterior_qualified') is not False or result.get('GT_direct_read') is not False:
        raise ValueError('simple-density selection provenance differs')
    rows={family:result[family] for family in FAMILIES}
    if int(rows['selected']['dimension'])<=0 or int(rows['unconditional_selected']['dimension'])!=0:
        raise ValueError('one selected conditional PCA model and one unconditional model required')
    sources={}
    for family,row in rows.items():
        source=(OUT/row['checkpoint']).resolve()
        if not source.is_relative_to((OUT/'D/simple_density/fits').resolve()):
            raise ValueError('selected fitted checkpoint must be inside simple_density/fits')
        if not source.is_file():raise FileNotFoundError(source)
        sources[family]=source
    return result,rows,sources


def load_simple(path,dimension,*,rng):
    model=SimpleConditionalGaussian.load(path,rng=rng)
    a=model.arrays
    shapes=dict(scalar_mean=(1154,),scalar_scale=(1154,),label_mean=(10,),
        label_scale=(10,),projection=(1154,dimension),weights=(dimension,10),covariance=(10,10))
    for name,shape in shapes.items():
        if a[name].shape!=shape or not np.isfinite(a[name]).all():
            raise ValueError(f'invalid fitted {name} shape/finite contract')
    if np.any(a['scalar_scale']<=0) or np.any(a['label_scale']<=0):
        raise ValueError('positive train-fitted raw-coordinate affine scales required')
    if not np.allclose(a['covariance'],a['covariance'].T,rtol=0,atol=1e-12):
        raise ValueError('symmetric normalized conditional Gaussian covariance required')
    return model


def run(max_wall=7200):
    import torch
    if not 1<=max_wall<=7200:raise ValueError('max wall must be within1..7200 seconds')
    result,rows,sources=selected_sources()
    frame_path=OUT/'D/FIXED_FRAME.npz'
    with np.load(frame_path,allow_pickle=False) as z:Vw=z['Vw'].copy();Vs=z['Vs'].copy()
    frame=np.column_stack((Vw,Vs))
    if Vw.shape!=(64,8) or Vs.shape!=(64,56) or not np.allclose(frame.T@frame,np.eye(64),atol=1e-10,rtol=0):
        raise ValueError('fixed complete orthogonal weak8/strong56 frame required')
    book=Book('D',JOB,max_wall=max_wall);status='ERROR';error=None
    previous_threads=torch.get_num_threads();torch.set_num_threads(1)
    summaries=[];context_summaries=[]
    try:
        if (book.jobdir/'PROTOCOL.json').exists():raise FileExistsError('no overwrite or automatic replay')
        write(book.jobdir/'PROTOCOL.json',dict(parents=PARENTS,families=FAMILIES,
            selected_fitted_controls=rows,selection_source=str(OUT/'D/simple_density/RESULTS.json'),
            selection_rule='reuse completed train/validation choices; no DEV or GT model selection',
            seed=SEED,device='cpu',transitions_per_pair=TRANSITIONS,
            defensive_standard_probit_prior_weight=.1,all_66_coordinates_retained=True,
            context=dict(method='coarse_TRF',max_evaluations=16,initial_point=None,
                native_original_prior_start=[0.]*64+[.5,.5],polish=False,
                shared_by_both_controls=True,new_states_and_J_and_cache_billed='D in this job'),
            initialization_preference=['completed B joint x_mode','old A probit mode','inward A-anchor fallback'],
            initialization_shared_setup='previous B/A cost disclosed; every control rebuilds its initial Full state paid D',
            MH_uniform_seed_per_parent='[2026102270,parent_index,1], shared across controls',
            proposal_numpy_seed_per_parent='[2026102270,parent_index,0]',
            model_numpy_seed_per_parent_family='[2026102270,parent_index,2,family_index]',
            MH_RHS_upper=16448,context_RHS_upper=2560,total_RHS_upper=19008,
            original_target_unchanged=True,production_adaptation=False,
            scientific_status='STRONG_FITTED_SIMPLE_PROPOSAL_BASELINE_NOT_POSTERIOR_CALIBRATION',
            model_proposal_moments='separate unweighted proposal diagnostics; not corrected-chain posterior means',
            GT_read=False,training_or_fitting=False,offline_risk_evaluator=False,
            max_wall_seconds=max_wall,global_paid_worker='shared Book one-worker lock',
            neural_status='NOT_RUN_VALIDATION_SELECTED_INITIALIZATION_NO_TRAINED_GAIN'))
        for rel in ['tools/run_p3l_simple_correction.py','tools/run_p3l_learned_correction.py',
            'tools/run_p3l_contexts.py','src/p3l/simple_density.py','src/p3l/learned_proposal.py',
            'src/p3l/inversion.py','src/p3l/core.py','src/p3l/moments.py','src/p1_material/prior.py']:
            destination=book.jobdir/'SOURCE'/rel;destination.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/rel,destination)
        shutil.copyfile(frame_path,book.jobdir/'FIXED_FRAME_FROZEN.npz')
        shutil.copyfile(OUT/'configs/TARGET.json',book.jobdir/'TARGET_FROZEN.json')
        shutil.copyfile(ROOT/'research/p1_material_target_20261018/inputs/geometry.json',book.jobdir/'GEOMETRY_FROZEN.json')
        write(book.jobdir/'FITTED_SELECTION_FROZEN.json',result)
        frozen={}
        for family in FAMILIES:
            destination=book.jobdir/'FITTED_MODELS'/f'{family}.npz'
            destination.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(sources[family],destination)
            # Validate all model arrays before any Target/context query.
            load_simple(destination,int(rows[family]['dimension']),rng=np.random.default_rng(SEED))
            frozen[family]=destination
        starts={parent:finite_start(parent) for parent in PARENTS}
        write(book.jobdir/'INITIALIZATIONS.json',dict(
            parents={p:dict(x=starts[p][0],metadata=starts[p][1]) for p in PARENTS},
            same_initial_material_for_both_controls=True,
            previous_setup_cost_shared_and_disclosed=True,all_new_initial_full_states_paid_D=True))
        target=Target(book)
        for index,parent in enumerate(PARENTS):
            y=target.observation(parent);before=costs(book)
            context,report=reconstruct(target,y,parent,method='coarse_TRF',
                max_evaluations=16,initial_point=None,polish=False)
            if report.get('failure') is not None:
                raise RuntimeError('context failure; preserve paid best cache and stop')
            scalar,image,J,residual=context_arrays(target,context,y)
            start_x,start_meta=starts[parent]
            save(book.jobdir/'CONTEXTS'/f'{parent}.npz',scalar=scalar,image=image,q=context.q,
                prediction=context.state.data,J=J,residual=residual,
                objective=np.array(report['objective']),Phi=np.array(report['phi']),
                endpoint_key=np.array(context.key),initial_x=start_x)
            target.persist_cache(context,book.jobdir/'CONTEXTS/cache'/f'{parent}.npz')
            setup=dict(parent=parent,source_job=book.job,report=report,cost=delta(before,costs(book)),
                context_law='MATCHED_TRAINING_COARSE_TRF16_ORIGINAL_PRIOR_START_NO_POLISH',
                shared_by=list(FAMILIES),no_per_model_double_billing=True,GT_read=False)
            write(book.jobdir/'CONTEXTS'/f'{parent}.json',setup);context_summaries.append(setup)
            features=(scalar,image,J,residual,context.q.copy(),Vw,Vs)
            del context
            for family_index,family in enumerate(FAMILIES):
                model_seed=[SEED,index,2,family_index]
                model=load_simple(frozen[family],int(rows[family]['dimension']),
                    rng=np.random.default_rng(np.random.SeedSequence(model_seed)))
                try:
                    summary=run_pair(book,target,parent,family,model,features,start_x,start_meta,index,
                        model_kind=KINDS[family],model_rng_seed=model_seed)
                except BaseException as exc:
                    write(book.jobdir/'PAIR_FAILURE.json',dict(parent=parent,model=family,
                        error=repr(exc),status='STOP_PRESERVE_PAID_ENDPOINTS_AND_LAST_CHAIN',GT_read=False))
                    raise
                summaries.append(summary)
                write(book.jobdir/'PROGRESS.json',dict(comparisons=summaries,shared_context_setup=context_summaries))
                del model
            target.backend()._transfers.clear()
        write(book.jobdir/'SUMMARY.json',dict(status='COMPLETE_FITTED_SIMPLE_FULL_MH_DIAGNOSTIC',
            comparisons=summaries,shared_context_setup=context_summaries,total_RHS_upper=19008,
            posterior_calibrated=False,GT_read=False,training_or_fitting=False,
            scientific_status='NORMALIZED_FULL_TARGET_INVARIANT_KERNEL_NOT_A_MIXING_CERTIFICATE'))
        status='COMPLETE'
    except BaseException as exc:error=repr(exc);raise
    finally:
        torch.set_num_threads(previous_threads);book.close(status,error)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-wall',type=int,default=7200)
    args=parser.parse_args();run(args.max_wall)


if __name__=='__main__':main()
