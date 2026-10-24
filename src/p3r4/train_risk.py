"""Bounded zero-physics risk training; run as python -m p3r4.train_risk.

Accepts R4 x/y/material or original R2 joint datasets. Optional separate anchor
contexts restrict both splits to their matching independent parent subsets.
"""
import argparse,copy,json,time,resource,sys
from pathlib import Path
import numpy as np
import torch
from scipy.special import ndtri
from .risk import KINDS,LocalRisk,frozen_geometry,feature_image,phase_labels,roi_components,full_material,load_r2_baseline
from p3r2.risk_models import observation_features
from p3r2.residual_risk import read_anchor,fit_normalization,online_features


def load_joint(path):
    if 'locked' in str(path).lower(): raise ValueError('locked data prohibited')
    with np.load(path,allow_pickle=False) as z:
        parent=z['parent'].astype(str)
        if 'material' in z:
            material=z['material'];y=z['y'];x=z['x'];s=z['s'] if 's' in z else None
        else:
            material=z['chi'];y=observation_features(z['observation']);q=z['q'];x=np.c_[q[:,:64],ndtri(q[:,64:])];s=None
    n=len(parent)
    if x.shape!=(n,66) or y.shape!=(n,512) or material.shape!=(n,576) or not np.iscomplexobj(material): raise ValueError('joint shapes')
    if len(set(parent))!=n or not all(np.isfinite(a).all() for a in (x,y,material)): raise ValueError('unique finite joint parents required')
    return dict(x=x.astype('float32'),y=y.astype('float32'),material=material,parent=parent,s=s)


def align(data,anchor):
    lookup={p:i for i,p in enumerate(data['parent'])}
    try: ids=np.array([lookup[p] for p in anchor['parent']])
    except KeyError as e: raise ValueError('anchor parent absent from joint split') from e
    return {k:(v[ids] if v is not None else None) for k,v in data.items()}


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('train','validation','out'): p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--kind',choices=KINDS,required=True);p.add_argument('--seed',type=int,default=11)
    p.add_argument('--updates',type=int,default=400);p.add_argument('--device',choices=('cpu','mps'),default='cpu')
    p.add_argument('--context-mode',choices=('spectral','decoder_only','random'),default='spectral')
    p.add_argument('--refit-feature-normalization',action='store_true',help='TRAIN-anchor-only refit for a new context population; historical default unchanged')
    p.add_argument('--batch-size',type=int,default=64);p.add_argument('--learning-rate',type=float,default=5e-4)
    p.add_argument('--train-context','--anchor-train',type=Path);p.add_argument('--validation-context','--anchor-validation',type=Path)
    default=Path(__file__).resolve().parents[3]/'b1-p3l-r2-risk-posterior/research/p3l_r2_risk_posterior_exploration/A/residual_training/n512/local_geometry_cnn_s33/BEST.pt'
    p.add_argument('--geometry-checkpoint',type=Path,default=default)
    a=p.parse_args(argv)
    if a.updates<0 or a.batch_size<1: p.error('nonnegative updates and positive batch required')
    if bool(a.train_context)!=bool(a.validation_context):p.error('both contexts required')
    wall,cpu=time.monotonic(),time.process_time();torch.set_num_threads(1);torch.manual_seed(a.seed);rng=np.random.default_rng(a.seed)
    tr,va=load_joint(a.train),load_joint(a.validation)
    if set(tr['parent'])&set(va['parent']):raise ValueError('parent overlap')
    payload=torch.load(a.geometry_checkpoint,map_location='cpu',weights_only=False);g=payload['geometry']
    tc=vc=None
    if a.train_context:
        tc,vc=read_anchor(a.train_context),read_anchor(a.validation_context);tr,va=align(tr,tc),align(va,vc)
    feature_norm=payload['normalization']
    if a.context_mode=='random':
        if tc is None:raise ValueError('random frame control requires genuine anchor context')
        g=copy.deepcopy(g);frame,_=np.linalg.qr(np.random.default_rng(2026102504).normal(size=(64,64)))
        g['weak']=frame[:,:g['weak'].shape[1]].astype('float32');g['strong']=frame[:,g['weak'].shape[1]:g['weak'].shape[1]+g['strong'].shape[1]].astype('float32')
        feature_norm=fit_normalization(online_features(tc,g))
    if a.refit_feature_normalization:
        if tc is None:raise ValueError('TRAIN-only context refit requires actual anchor contexts')
        feature_norm=fit_normalization(online_features(tc,g))
    ti=feature_image(tr['y'],g,tc,feature_norm);vi=feature_image(va['y'],g,vc,feature_norm)
    if a.context_mode=='decoder_only':
        ti[:,2:6]=0;vi[:,2:6]=0
    if tc is None:
        mean=ti.mean((0,2,3),keepdims=True);scale=np.maximum(ti.std((0,2,3),keepdims=True),1e-8)
        ti=(ti-mean)/scale;vi=(vi-mean)/scale
    else:mean=scale=None
    target=roi_components(tr['material'],g['roi']);vt=roi_components(va['material'],g['roi'])
    ta=roi_components(tc['chi0'],g['roi']) if tc is not None else np.zeros_like(target)
    vanc=roi_components(vc['chi0'],g['roi']) if vc is not None else np.zeros_like(vt)
    out=a.out;out.mkdir(parents=True,exist_ok=True)
    common=dict(kind=a.kind,seed=a.seed,train_size=len(ti),validation_size=len(vi),physics_calls=0,
                information_set='R2 online anchor local image and inferred u2' if tc is not None else 'source-wise mean packed observation8; insufficient for full-y claim',
                conditional_density='NOT_CLAIMED',context_mode=a.context_mode,feature_normalization=feature_norm,geometry=g,feature_mean=mean,feature_scale=scale,
                feature_normalization_provenance='new TRAIN-only contexts' if a.refit_feature_normalization else 'frozen R2 TRAIN; random frame refit only',
                geometry_checkpoint=str(a.geometry_checkpoint.resolve()),training_parents=tr['parent'])
    if a.kind=='r2_cnn':
        if vc is None or a.context_mode!='spectral':raise ValueError('frozen R2 baseline requires original spectral anchor context')
        baseline=load_r2_baseline(a.geometry_checkpoint,a.device)
        raw=baseline.predict_full(vc);vp=roi_components(raw,g['roi'])
        count=sum(z.numel() for z in baseline.model.parameters());history=[];best_update=payload['selected_update'] if 'selected_update' in payload else payload.get('best_update',0)
        common['historical_state_dict']=payload['state_dict']
        common['historical_payload']=payload
        common['baseline_contract']='exact frozen R2 raw prediction; historical validation-selected projection reported separately'
    elif a.kind in ('affine','sharpen','spatial'):
        if tc is None:raise ValueError('classical controls require genuine anchor context')
        from scipy.ndimage import uniform_filter
        cfg=g['config'];host=np.array(cfg['host_chi'])[None,:,None]
        if a.kind=='affine':
            # TRAIN-only shared slope/intercept for each complex component.
            coef=[];pred=[]
            for c in range(2):
                design=np.c_[ta[:,c].ravel(),np.ones(ta[:,c].size)]
                beta=np.linalg.lstsq(design,target[:,c].ravel(),rcond=None)[0];coef.append(beta.tolist());pred.append(vanc[:,c]*beta[0]+beta[1])
            vp=np.stack(pred,1);common['coefficient']=coef
        elif a.kind=='spatial':
            full=vc['chi0'];im=np.stack((full.real,full.imag),1).reshape(-1,2,24,24)
            vp=uniform_filter(im,size=(1,1,3,3),mode='nearest').reshape(-1,2,576)[:,:,g['roi']]
        else:
            # Prespecified small grid selected by validation, no test access.
            delta0=np.stack((cfg['inclusion_increment_real_interval'][0]+np.diff(cfg['inclusion_increment_real_interval'])[0]*vc['q0'][:,64],cfg['inclusion_increment_imag_interval'][0]+np.diff(cfg['inclusion_increment_imag_interval'])[0]*vc['q0'][:,65]),1)[:,:,None]
            fraction=np.clip((vanc-host)/np.maximum(delta0,1e-8),0,1);best=float('inf')
            for slope in (0.5,0.75,1.,1.25,1.5,2.,3.):
                candidate=host+delta0*np.clip(.5+slope*(fraction-.5),0,1);score=np.mean((candidate-vt)**2)
                if score<best:best=score;vp=candidate;common['slope']=slope
        count=0;history=[];best_update=0
    else:
        model=LocalRisk(a.kind,g).to(a.device);count=sum(z.numel() for z in model.parameters())
        if count>30000:raise ValueError('parameter cap')
        optimizer=torch.optim.AdamW(model.parameters(),lr=a.learning_rate,weight_decay=1e-3)
        thresholds=model.thresholds.detach().cpu().numpy() if model.k else None
        # Build only batch labels to bound memory for32768 x16 x256.
        def prediction(image,base):
            return model(torch.tensor(image,device=a.device),torch.tensor(base,device=a.device) if not model.k else None)
        def evaluate():
            model.eval();total=0.;rows=[]
            with torch.inference_mode():
                for start in range(0,len(vi),a.batch_size):
                    pred=prediction(vi[start:start+a.batch_size],vanc[start:start+a.batch_size]).cpu().numpy();rows.append(pred);total+=float(np.sum((pred-vt[start:start+a.batch_size])**2))
            return total/vt.size,np.concatenate(rows)
        score,vp=evaluate();best=score;state=copy.deepcopy(model.state_dict());best_update=0
        history=[dict(update=0,validation_component_mse=score)]
        for update in range(1,a.updates+1):
            ids=rng.integers(len(ti),size=min(a.batch_size,len(ti)));model.train();optimizer.zero_grad(set_to_none=True)
            if model.k:
                _,mom=model(torch.tensor(ti[ids],device=a.device),return_moments=True)
                labels=phase_labels(tr['x'][ids],g,thresholds,None if tr['s'] is None else tr['s'][ids])
                # Actual joint complex squared moment loss; weights integrate h'.
                loss=(((mom-torch.tensor(labels,device=a.device))**2)*model.weights[None,None,:,None]).sum(2).mean()
            else:loss=((prediction(ti[ids],ta[ids])-torch.tensor(target[ids],device=a.device))**2).mean()
            if not torch.isfinite(loss):raise RuntimeError('nonfinite loss')
            loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5.);optimizer.step()
            if update%max(1,min(100,a.updates))==0 or update==a.updates:
                score,current=evaluate();history.append(dict(update=update,validation_component_mse=score))
                if score<best:best=score;state=copy.deepcopy(model.state_dict());best_update=update;vp=current
        common['state_dict']={k:v.cpu() for k,v in state.items()}
    prediction_full=full_material(vp,g)
    rows=np.mean(np.abs(prediction_full-va['material'])**2,1)
    host=complex(*g['config']['host_chi']);upper=host+complex(g['config']['inclusion_increment_real_interval'][1],g['config']['inclusion_increment_imag_interval'][1])
    clipped=np.clip(prediction_full.real,host.real,upper.real)+1j*np.clip(prediction_full.imag,host.imag,upper.imag)
    clipped_rows=np.mean(np.abs(clipped-va['material'])**2,1)
    common['selected_update']=best_update;torch.save(common,out/'BEST.pt')
    np.savez_compressed(out/'validation.npz',parent=va['parent'],prediction=prediction_full,risk=rows,clipped_prediction=clipped,clipped_risk=clipped_rows)
    result=dict(kind=a.kind,seed=a.seed,parameter_count=count,updates_requested=a.updates,selected_update=best_update,
                train_size=len(ti),validation_size=len(vi),validation_full_complex_mse=float(rows.mean()),validation_clipped_full_complex_mse=float(clipped_rows.mean()),
                wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,
                peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024),
                physics_calls=0,history=history,information_set=common['information_set'],density='NOT_CLAIMED')
    (out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
