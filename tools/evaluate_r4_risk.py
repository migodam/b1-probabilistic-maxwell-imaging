#!/usr/bin/env python3
"""Freeze, online inference, then separate-process offline risk evaluation.

freeze --methods methods.json --out METHODS_FROZEN.json
infer --freeze METHODS_FROZEN.json --contexts contexts_dir_or_packed.npz --out ONLINE_DIR
score --freeze METHODS_FROZEN.json --online ONLINE_DIR --truth inputs/offline_only --out SCORE_DIR
"""
import argparse,json,os,sys,time,resource
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.evaluation import *


def load_json(path):return json.loads(Path(path).read_text())


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('command',choices=('freeze','infer','score'))
    ap.add_argument('--methods',type=Path);ap.add_argument('--freeze',type=Path);ap.add_argument('--contexts',type=Path)
    ap.add_argument('--online',type=Path);ap.add_argument('--truth',type=Path);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--device',choices=('cpu','mps'),default='cpu');a=ap.parse_args()
    wall,cpu=time.monotonic(),time.process_time()
    if a.command=='freeze':
        spec=load_json(a.methods)
        field_only=all(m['loader']=='field' for m in spec['methods'])
        if field_only:
            if spec.get('fallback')!='none':raise ValueError('measurement-only methods require explicit fallback=none')
        elif spec.get('fallback')!='anchor':raise ValueError('explicit fallback=anchor protocol required')
        ids=[m['id'] for m in spec['methods']]
        if len(ids)!=len(set(ids)) or any(not i.replace('_','').replace('-','').isalnum() for i in ids):raise ValueError('unique safe method IDs required')
        parents=spec['parents']
        if len(parents)!=len(set(parents)) or any('/' in p or '..' in p for p in parents):raise ValueError('unique safe parent IDs required')
        for m in spec['methods']:
            if m['loader'] not in ('r2','r4','anchorbin','anchor','field','fusion'):raise ValueError('loader must be r2/r4/anchorbin/anchor/field/fusion')
            if 'seed' not in m or 'family' not in m:raise ValueError('family and seed required')
            if m['loader']!='anchor':m['checkpoint_identity']=checkpoint_identity(m['checkpoint']);m['checkpoint']=m['checkpoint_identity']['path']
            if m['loader']=='fusion':
                from p3r4.risk_fusion import FusionPredictor,read_moments
                fusion=FusionPredictor(m['checkpoint'],'cpu')
                if not m.get('moments_directory'):raise ValueError('fusion requires frozen online moments_directory')
                m['moments_directory']=str(Path(m['moments_directory']).resolve())
                read_moments(m['moments_directory'],spec['parents'],fusion.p['moment_checkpoint_identity'],fusion.moment_provider,fusion.moment_K)
                m['moment_provider']=fusion.moment_provider;m['moment_K']=fusion.moment_K
                m['moment_gate_identity']=checkpoint_identity(Path(m['moments_directory'])/'MOMENTS_COMPLETE.json')
                del fusion
            if m['loader']=='field':
                from p3r4.field_moments import load_field_moments
                model=load_field_moments(m['checkpoint'],'cpu')
                if hasattr(model,'geometry'):g=model.geometry
                else:
                    from p3r4.probability_eval import decoder_from_geometry
                    cfg=load_json(ROOT/'research/p3l_r4_conditional_transport/configs/TARGET.json')
                    points=load_json(ROOT/'research/p1_material_target_20261018/inputs/geometry.json')['points']
                    g=decoder_from_geometry(cfg,points)
                    if not np.array_equal(model.reference.roi_indices,g['roi']):raise ValueError('historical field ROI mismatch')
                    if model.p.get('target_id') not in (None,cfg['target_id']):raise ValueError('historical field target mismatch')
                if 'field_geometry' not in spec:spec['field_geometry']={'config':g['config'],'roi':np.asarray(g['roi']).tolist(),'A':np.asarray(g['A']).tolist()}
                elif not np.array_equal(spec['field_geometry']['roi'],g['roi']):raise ValueError('field geometry mismatch')
                m.update(information_set=model.p.get('contract',{}).get('information_set','all raw512 measurements only; no anchor/context/J/residual inputs'),
                         fixed_action_ignores_y=model.p.get('architecture')=='train_mean',
                         moment_origin=model.p.get('moment_source','historical mean-only reference'),
                         fallback='none; field prediction failures never use anchor',
                         field_second_moments=bool(model.p.get('second_moments',False)),
                         moment_scope='pointwise joint Re/Im moments only; no spatial posterior or density')
                del model
        if not spec.get('comparisons'):raise ValueError('predeclare comparisons')
        for i,c in enumerate(spec['comparisons']):
            if c['action'] not in ('raw','clipped'):raise ValueError('comparison action raw/clipped')
            c.setdefault('id',f'comparison_{i:02d}')
        spec.update(status='FROZEN_BEFORE_LOCKED_TRUTH',freeze_time_ns=time.time_ns(),freeze_pid=os.getpid(),truth_read=False,physics_calls=0,
                    statistics=dict(draws=4000,seed=BOOTSTRAP_SEED,primary_CI=CI_LEVEL,unit='parent; average seed losses'))
        write_new_json(a.out,spec);print(json.dumps({'status':'FROZEN','path':str(a.out)}));return
    frozen=load_json(a.freeze)
    if frozen['status']!='FROZEN_BEFORE_LOCKED_TRUTH':raise ValueError('invalid freeze')
    if a.command=='infer':
        # No truth argument is accepted in the inference stage.
        if a.truth is not None:raise ValueError('truth path prohibited during online inference')
        from p3r4.risk import RiskActionPredictor,load_r2_baseline
        from p3r2.risk_models import observation_features
        from p3r2.fusion import passive_box
        from p3r4.field_moments import load_field_moments
        import torch
        torch.set_num_threads(1)
        a.out.mkdir(parents=True,exist_ok=True)
        if (a.out/'ONLINE_COMPLETE.json').exists():raise FileExistsError('online run already complete')
        geometry=None;models={};load_errors={}
        for m in frozen['methods']:
            if m['loader']=='anchor':models[m['id']]=None;continue
            if checkpoint_identity(m['checkpoint'])!=m['checkpoint_identity']:raise ValueError('checkpoint changed since freeze')
            try:
                if m['loader']=='anchorbin':
                    from p3r4.risk_anchorbin import AnchorBinPredictor
                    model=AnchorBinPredictor(m['checkpoint'],a.device)
                elif m['loader']=='fusion':
                    from p3r4.risk_fusion import FusionPredictor
                    if checkpoint_identity(Path(m['moments_directory'])/'MOMENTS_COMPLETE.json')!=m['moment_gate_identity']:raise ValueError('fusion moment gate changed')
                    model=FusionPredictor(m['checkpoint'],a.device)
                elif m['loader']=='field':model=load_field_moments(m['checkpoint'],a.device)
                else:
                    from p3r4.risk_anchorbin import load_risk_compatible
                    model=load_r2_baseline(m['checkpoint'],a.device) if m['loader']=='r2' else load_risk_compatible(m['checkpoint'],a.device)
                models[m['id']]=model;g=getattr(model,'geometry',frozen.get('field_geometry')) if m['loader']=='field' else model.p['geometry'] if m['loader']=='r2' else model.geometry
                if geometry is None:geometry=g
                elif not np.array_equal(geometry['roi'],g['roi']):raise ValueError('ROI mismatch')
            except Exception as e:models[m['id']]=None;load_errors[m['id']]=repr(e)
        if geometry is None:raise ValueError('at least one loadable model must supply frozen geometry')
        parents=frozen['parents'];anchors=[];q0=[];all_errors={m['id']:[] for m in frozen['methods']}
        predictions={m['id']:[] for m in frozen['methods']};flags={m['id']:[] for m in frozen['methods']}
        need_anchor=any(m['loader']!='field' for m in frozen['methods'])
        field_moments={m['id']:[] for m in frozen['methods'] if m['loader']=='field'}
        field_failures={key:[] for key in field_moments};field_warnings={key:[] for key in field_moments}
        for parent in parents:
            path=a.contexts/f'{parent}.npz' if a.contexts.is_dir() else a.contexts
            context_error=None
            if need_anchor:
                try:online=read_online_context(path,parent);anchor=online['chi0'];q=online['q0']
                except Exception as e:
                    context_error=repr(e);online=None;anchor,q=read_anchor_fallback(path,parent)
                anchors.append(anchor[0]);q0.append(q[0])
            else:online=anchor=q=None
            for m in frozen['methods']:
                error=None;fallback=False;moment_result=None
                try:
                    if m['loader']=='field':
                        if m['id'] in load_errors:raise ValueError(load_errors[m['id']])
                        raw_y=read_online_measurement(path,parent)
                        moment_result=models[m['id']].predict_moments(raw_y,batch_size=16)
                        prediction=moment_result['material_mean']
                    elif m['loader']=='anchor':prediction=anchor.copy()
                    else:
                        if context_error:raise ValueError(context_error)
                        if m['id'] in load_errors:raise ValueError(load_errors[m['id']])
                        model=models[m['id']]
                        if m['loader']=='fusion':
                            from p3r4.risk_fusion import read_moments
                            extra=read_moments(m['moments_directory'],[parent],model.p['moment_checkpoint_identity'],model.moment_provider,model.moment_K)
                            prediction=model.predict(observation_features(online['observation']),online,moments=extra,parents=[parent],moment_provider=model.moment_provider,moment_K=model.moment_K)
                        else:prediction=model.predict_full(online) if m['loader']=='r2' else model.predict(observation_features(online['observation']),online)
                    if prediction.shape!=(1,576) or not np.iscomplexobj(prediction) or not np.isfinite(prediction).all():raise ValueError('invalid material prediction')
                except Exception as e:
                    error=repr(e)
                    if m['loader']=='field':prediction=np.full((1,576),np.nan+1j*np.nan);moment_result=None
                    else:prediction=anchor.copy();fallback=True
                predictions[m['id']].append(prediction[0]);flags[m['id']].append(fallback)
                if error:all_errors[m['id']].append(dict(parent=parent,error=error,fallback='none' if m['loader']=='field' else 'anchor'))
                if m['loader']=='field':
                    field_moments[m['id']].append(moment_result);field_failures[m['id']].append(error is not None)
                    if moment_result is not None:field_warnings[m['id']].append(dict(parent=parent,status=moment_result['second_moment_status'],projection=moment_result.get('psd_projection_warning',{})))
        if need_anchor:np.savez_compressed(a.out/'online_anchor.npz',parent=np.array(parents),chi0=np.array(anchors),q0=np.array(q0))
        field_artifacts=[]
        for m in frozen['methods']:
            raw=np.array(predictions[m['id']]);clipped=passive_box(raw,geometry['config'],geometry['roi'])
            dest=a.out/(m['id']+'.npz')
            if dest.exists():raise FileExistsError(dest)
            arrays=dict(parent=np.array(parents),prediction_raw=raw,prediction_clipped=clipped,fallback_anchor=np.array(flags[m['id']],dtype=bool))
            if m['loader']=='field':
                records=field_moments[m['id']]
                keys=set().union(*(set(k for k,v in row.items() if isinstance(v,np.ndarray) and k!='material_mean') for row in records if row is not None))
                for key in keys:
                    example=next(row[key][0] for row in records if row is not None and key in row)
                    arrays['field__'+key]=np.stack([row[key][0] if row is not None and key in row else np.full(example.shape,np.nan) for row in records])
                arrays['prediction_failed']=np.array(field_failures[m['id']],bool)
            np.savez_compressed(dest,**arrays)
            if m['loader']=='field':field_artifacts.append(checkpoint_identity(dest))
        write_new_json(a.out/'ONLINE_COMPLETE.json',dict(status='ONLINE_COMPLETE',pid=os.getpid(),completed_time_ns=time.time_ns(),
            freeze=frozen,truth_read=False,physics_calls=0,parents=parents,geometry={'config':geometry['config'],'roi':np.asarray(geometry['roi']).tolist(),'A':np.asarray(geometry['A']).tolist()},
            failures=all_errors,fallback=('anchor; all parents retained' if not field_moments else 'anchor for context methods; none for field' if need_anchor else 'none; field failures retained'),wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,
            anchor_context_used=need_anchor,field_artifact_identities=field_artifacts,field_moment_warnings=field_warnings,
            field_information='raw512 only; field failures never use anchor; moments correspond to raw mean, not clipped action',
            peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024)))
        print(json.dumps({'status':'ONLINE_COMPLETE','parents':len(parents),'failures':{k:len(v) for k,v in all_errors.items()}}));return
    gate=load_json(a.online/'ONLINE_COMPLETE.json')
    if gate['status']!='ONLINE_COMPLETE' or gate['truth_read'] or gate['freeze']!=frozen:raise ValueError('missing/mismatched online gate')
    if gate['pid']==os.getpid():raise ValueError('offline scoring must be a separate process')
    if gate['completed_time_ns']<=frozen['freeze_time_ns']:raise ValueError('invalid freeze chronology')
    for m in frozen['methods']:
        if m['loader']=='field' and checkpoint_identity(m['checkpoint'])!=m['checkpoint_identity']:raise ValueError('field checkpoint changed after freeze')
    for identity in gate.get('field_artifact_identities',[]):
        if checkpoint_identity(identity['path'])!=identity:raise ValueError('sealed field prediction artifact changed')
    # All checks above precede the first locked truth access.
    parents=frozen['parents'];truth=[];truth_phase=[];g=gate['geometry'];A=np.array(g['A']);roi=np.array(g['roi']);cfg=g['config']
    for parent in parents:
        with np.load(a.truth/f'{parent}.npz',allow_pickle=False) as z:
            truth.append(z['chi']);q=z['q'];truth_phase.append((cfg['field_mean']+A[roi]@q[:64])/cfg['phase_transition_halfwidth'])
    truth=np.array(truth);truth_phase=np.array(truth_phase)
    anchor_phase=None
    if gate.get('anchor_context_used',True):
        with np.load(a.online/'online_anchor.npz',allow_pickle=False) as z:
            if z['parent'].astype(str).tolist()!=parents:raise ValueError('anchor identity mismatch')
            aq=z['q0'];anchor_phase=(cfg['field_mean']+aq[:,:64]@A[roi].T)/cfg['phase_transition_halfwidth']
    losses={};summary={};a.out.mkdir(parents=True,exist_ok=True)
    for m in frozen['methods']:
        with np.load(a.online/(m['id']+'.npz'),allow_pickle=False) as z:
            if z['parent'].astype(str).tolist()!=parents:raise ValueError('prediction parent mismatch')
            fallback=z['fallback_anchor'];losses[m['id']]={}
            perparent={'parent':np.array(parents),'fallback_anchor':fallback}
            for action in ('raw','clipped'):
                rows=material_losses(z['prediction_'+action],truth,roi,anchor_phase,truth_phase);losses[m['id']][action]=rows
                for k,v in rows.items():perparent[action+'__'+k]=v
            if m['loader']=='field':
                perparent['prediction_failed']=z['prediction_failed']
                if 'field__raw_second_moments' in z:
                    predicted=z['field__raw_second_moments']
                    labels=np.stack((truth.real**2,truth.imag**2,truth.real*truth.imag),axis=1)
                    perparent['field_second_moment_full_mse']=((predicted-labels)**2).mean((1,2))
                    for j,name in enumerate(('real_square','imag_square','real_imag')):
                        perparent['field_second_moment_'+name+'_mse']=((predicted[:,j]-labels[:,j])**2).mean(1)
        dest=a.out/(m['id']+'_per_parent.npz')
        if dest.exists():raise FileExistsError(dest)
        np.savez_compressed(dest,**perparent)
        summary[m['id']]={'family':m['family'],'seed':m['seed'],'fallback_count':int(fallback.sum()),'raw':{},'clipped':{}}
        if m['loader']=='field':
            summary[m['id']].update(information_set=m['information_set'],prediction_failures=int(perparent['prediction_failed'].sum()),
                moment_scope=m['moment_scope'],variance_calibration='NOT_ESTABLISHED',
                second_moment_scores={k:dict(mean=float(v.mean()) if np.isfinite(v).all() else None,finite_parents=int(np.isfinite(v).sum())) for k,v in perparent.items() if k.startswith('field_second_moment_')})
        for action in ('raw','clipped'):
            for metric,rows in losses[m['id']][action].items():
                finite=np.isfinite(rows)
                value=float(rows[finite].mean()) if finite.any() else None
                if m['loader']=='field' and metric in ('full_mse','full_mae','roi_mse','roi_mae') and not finite.all():value=None
                summary[m['id']][action][metric]={'mean':value,'parents_nonempty':int(finite.sum())}
    comparisons=[]
    for i,c in enumerate(frozen['comparisons']):
        # Older frozen DEV specifications omitted this artifact-only name.
        # Keep their frozen scientific comparisons unchanged on recovery.
        comparison_id=c.get('id',f'comparison_{i:02d}')
        candidate,cm=family_parent_mean(losses,frozen['methods'],c['candidate_family'],c['action'])
        reference,rm=family_parent_mean(losses,frozen['methods'],c['reference_family'],c['action'])
        if np.isfinite(candidate).all() and np.isfinite(reference).all():
            row={**c,**paired_bootstrap(candidate,reference), 'candidate_seed_ids':[m['seed'] for m in cm],'reference_seed_ids':[m['seed'] for m in rm]}
        else:row={**c,'status':'NOT_EVALUATED_INCOMPLETE_PREDICTIONS','candidate_seed_ids':[m['seed'] for m in cm],'reference_seed_ids':[m['seed'] for m in rm]}
        comparisons.append(row)
        np.savez_compressed(a.out/(comparison_id+'_paired_parent_losses.npz'),parent=np.array(parents),candidate_parent_mean_loss=candidate,reference_parent_mean_loss=reference)
    write_new_json(a.out/'RESULTS.json',dict(status='OFFLINE_SCORE_COMPLETE',n_parent=len(parents),methods=summary,comparisons=comparisons,
        interpretation='finite frozen method risk comparison; no posterior or Maxwell compatibility qualification',
        all_failures=gate['failures'],physics_calls=0,training_updates=0,truth_access='separate process after ONLINE_COMPLETE',
        scoring_pid=os.getpid(),online_pid=gate['pid'],wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,
        peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024)))
    print(json.dumps({'status':'OFFLINE_SCORE_COMPLETE','parents':len(parents),'comparisons':comparisons}))

if __name__=='__main__':main()
