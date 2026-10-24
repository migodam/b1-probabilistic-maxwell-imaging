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
        if spec.get('fallback')!='anchor':raise ValueError('explicit fallback=anchor protocol required')
        ids=[m['id'] for m in spec['methods']]
        if len(ids)!=len(set(ids)) or any(not i.replace('_','').replace('-','').isalnum() for i in ids):raise ValueError('unique safe method IDs required')
        parents=spec['parents']
        if len(parents)!=len(set(parents)) or any('/' in p or '..' in p for p in parents):raise ValueError('unique safe parent IDs required')
        for m in spec['methods']:
            if m['loader'] not in ('r2','r4','anchor'):raise ValueError('loader must be r2/r4/anchor')
            if 'seed' not in m or 'family' not in m:raise ValueError('family and seed required')
            if m['loader']!='anchor':m['checkpoint_identity']=checkpoint_identity(m['checkpoint']);m['checkpoint']=m['checkpoint_identity']['path']
        if not spec.get('comparisons'):raise ValueError('predeclare comparisons')
        for c in spec['comparisons']:
            if c['action'] not in ('raw','clipped'):raise ValueError('comparison action raw/clipped')
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
        import torch
        torch.set_num_threads(1)
        a.out.mkdir(parents=True,exist_ok=True)
        if (a.out/'ONLINE_COMPLETE.json').exists():raise FileExistsError('online run already complete')
        geometry=None;models={};load_errors={}
        for m in frozen['methods']:
            if m['loader']=='anchor':models[m['id']]=None;continue
            if checkpoint_identity(m['checkpoint'])!=m['checkpoint_identity']:raise ValueError('checkpoint changed since freeze')
            try:
                model=load_r2_baseline(m['checkpoint'],a.device) if m['loader']=='r2' else RiskActionPredictor(m['checkpoint'],a.device)
                models[m['id']]=model;g=model.p['geometry'] if m['loader']=='r2' else model.geometry
                if geometry is None:geometry=g
                elif not np.array_equal(geometry['roi'],g['roi']):raise ValueError('ROI mismatch')
            except Exception as e:models[m['id']]=None;load_errors[m['id']]=repr(e)
        if geometry is None:raise ValueError('at least one loadable model must supply frozen geometry')
        parents=frozen['parents'];anchors=[];q0=[];all_errors={m['id']:[] for m in frozen['methods']}
        predictions={m['id']:[] for m in frozen['methods']};flags={m['id']:[] for m in frozen['methods']}
        for parent in parents:
            path=a.contexts/f'{parent}.npz' if a.contexts.is_dir() else a.contexts
            context_error=None
            try:online=read_online_context(path,parent);anchor=online['chi0'];q=online['q0']
            except Exception as e:
                context_error=repr(e);online=None;anchor,q=read_anchor_fallback(path,parent)
            anchors.append(anchor[0]);q0.append(q[0])
            for m in frozen['methods']:
                error=None;fallback=False
                try:
                    if m['loader']=='anchor':prediction=anchor.copy()
                    else:
                        if context_error:raise ValueError(context_error)
                        if m['id'] in load_errors:raise ValueError(load_errors[m['id']])
                        model=models[m['id']]
                        prediction=model.predict_full(online) if m['loader']=='r2' else model.predict(observation_features(online['observation']),online)
                    if prediction.shape!=(1,576) or not np.iscomplexobj(prediction) or not np.isfinite(prediction).all():raise ValueError('invalid material prediction')
                except Exception as e:error=repr(e);prediction=anchor.copy();fallback=True
                predictions[m['id']].append(prediction[0]);flags[m['id']].append(fallback)
                if error:all_errors[m['id']].append(dict(parent=parent,error=error,fallback='anchor'))
        np.savez_compressed(a.out/'online_anchor.npz',parent=np.array(parents),chi0=np.array(anchors),q0=np.array(q0))
        for m in frozen['methods']:
            raw=np.array(predictions[m['id']]);clipped=passive_box(raw,geometry['config'],geometry['roi'])
            dest=a.out/(m['id']+'.npz')
            if dest.exists():raise FileExistsError(dest)
            np.savez_compressed(dest,parent=np.array(parents),prediction_raw=raw,prediction_clipped=clipped,fallback_anchor=np.array(flags[m['id']],dtype=bool))
        write_new_json(a.out/'ONLINE_COMPLETE.json',dict(status='ONLINE_COMPLETE',pid=os.getpid(),completed_time_ns=time.time_ns(),
            freeze=frozen,truth_read=False,physics_calls=0,parents=parents,geometry={'config':geometry['config'],'roi':np.asarray(geometry['roi']).tolist(),'A':np.asarray(geometry['A']).tolist()},
            failures=all_errors,fallback='anchor; all parents retained',wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,
            peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024)))
        print(json.dumps({'status':'ONLINE_COMPLETE','parents':len(parents),'failures':{k:len(v) for k,v in all_errors.items()}}));return
    gate=load_json(a.online/'ONLINE_COMPLETE.json')
    if gate['status']!='ONLINE_COMPLETE' or gate['truth_read'] or gate['freeze']!=frozen:raise ValueError('missing/mismatched online gate')
    if gate['pid']==os.getpid():raise ValueError('offline scoring must be a separate process')
    if gate['completed_time_ns']<=frozen['freeze_time_ns']:raise ValueError('invalid freeze chronology')
    # All checks above precede the first locked truth access.
    parents=frozen['parents'];truth=[];truth_phase=[];g=gate['geometry'];A=np.array(g['A']);roi=np.array(g['roi']);cfg=g['config']
    for parent in parents:
        with np.load(a.truth/f'{parent}.npz',allow_pickle=False) as z:
            truth.append(z['chi']);q=z['q'];truth_phase.append((cfg['field_mean']+A[roi]@q[:64])/cfg['phase_transition_halfwidth'])
    truth=np.array(truth);truth_phase=np.array(truth_phase)
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
        dest=a.out/(m['id']+'_per_parent.npz')
        if dest.exists():raise FileExistsError(dest)
        np.savez_compressed(dest,**perparent)
        summary[m['id']]={'family':m['family'],'seed':m['seed'],'fallback_count':int(fallback.sum()),'raw':{},'clipped':{}}
        for action in ('raw','clipped'):
            for metric,rows in losses[m['id']][action].items():
                finite=np.isfinite(rows)
                summary[m['id']][action][metric]={'mean':float(rows[finite].mean()) if finite.any() else None,'parents_nonempty':int(finite.sum())}
    comparisons=[]
    for c in frozen['comparisons']:
        candidate,cm=family_parent_mean(losses,frozen['methods'],c['candidate_family'],c['action'])
        reference,rm=family_parent_mean(losses,frozen['methods'],c['reference_family'],c['action'])
        row={**c,**paired_bootstrap(candidate,reference), 'candidate_seed_ids':[m['seed'] for m in cm],'reference_seed_ids':[m['seed'] for m in rm]}
        comparisons.append(row)
        np.savez_compressed(a.out/(c['id']+'_paired_parent_losses.npz'),parent=np.array(parents),candidate_parent_mean_loss=candidate,reference_parent_mean_loss=reference)
    write_new_json(a.out/'RESULTS.json',dict(status='OFFLINE_SCORE_COMPLETE',n_parent=len(parents),methods=summary,comparisons=comparisons,
        interpretation='finite frozen method risk comparison; no posterior or Maxwell compatibility qualification',
        all_failures=gate['failures'],physics_calls=0,training_updates=0,truth_access='separate process after ONLINE_COMPLETE',
        scoring_pid=os.getpid(),online_pid=gate['pid'],wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,
        peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024)))
    print(json.dumps({'status':'OFFLINE_SCORE_COMPLETE','parents':len(parents),'comparisons':comparisons}))

if __name__=='__main__':main()
