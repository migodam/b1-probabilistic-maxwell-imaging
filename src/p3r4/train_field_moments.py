"""One bounded raw-measurement material mean/moment fit, with zero physics.

TRAIN/validation: y[N,512]+material[N,576], or original complex
observation[N,32,8]+chi[N,576]. x[N,66] may reconstruct labels through the
unchanged pure P1 decoder only when material labels are absent. No anchor,
Jacobian, residual or GT-derived feature is ever a network input.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import resource
import signal
import sys
import time

import numpy as np
import torch

if __package__ in (None,''):
    sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
    from src.p3r4.field_moments import *
    from src.p3r4.probability_eval import decoder_from_geometry,decode_gaussianized
    from src.p3r4.evaluation import checkpoint_identity
else:
    from .field_moments import *
    from .probability_eval import decoder_from_geometry,decode_gaussianized
    from .evaluation import checkpoint_identity


def write_json(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def load_joint(path,geometry):
    path=Path(path)
    if 'locked' in str(path).lower() or 'offline_only' in path.resolve().parts:
        raise ValueError('locked/offline_only data prohibited in supervised training')
    with np.load(path,allow_pickle=False) as z:
        y=packed_observation(z['y'] if 'y' in z else z['observation'])
        parent_key='parent' if 'parent' in z else 'parents'
        parents=z[parent_key].astype(str)
        if 'material' in z:
            material=z['material'];source='stored original material labels'
        elif 'chi' in z:
            material=z['chi'];source='stored original chi labels'
        else:
            x=z['x']
            if x.shape!=(len(y),66):raise ValueError('Gaussianized x labels required for pure decoder reconstruction')
            material=np.empty((len(y),576),np.complex128)
            for start in range(0,len(y),128):
                material[start:start+128]=decode_gaussianized(x[start:start+128],geometry)[0]
            source='reconstructed from supplied Gaussianized labels through original pure P1 decoder; not new simulation'
    if len(y)<2 or parents.shape!=(len(y),) or len(set(parents))!=len(parents):
        raise ValueError('at least2 independent unique parents required')
    if material.shape!=(len(y),576) or not np.iscomplexobj(material) or not np.isfinite(material).all():
        raise ValueError('finite full576 complex material labels required')
    return dict(y=y,material=material,parent=parents,label_source=source)


def fit_statistics(data,roi,second_moments=False,batch_size=256):
    y=data['y'];mean_y=y.mean(0,dtype=np.float64);scale_y=y.std(0,dtype=np.float64)
    scale_y=np.where(scale_y<1e-12,1.,scale_y)
    channels=5 if second_moments else 2
    total=np.zeros((channels,len(roi)),np.float64);squares=np.zeros_like(total)
    for start in range(0,len(y),batch_size):
        label=moment_labels(data['material'][start:start+batch_size],roi,second_moments).astype(np.float64)
        total+=label.sum(0);squares+=(label*label).sum(0)
    mean=total/len(y);variance=np.maximum(squares/len(y)-mean*mean,0.)
    # A common Re/Im scale keeps the mean term proportional to the physical
    # full-grid complex MSE. Second-moment channels have separate TRAIN scales.
    mean_scale=np.sqrt(variance[:2].mean())
    mean_scale=mean_scale if mean_scale>=1e-8 else 1.
    scale=np.ones((channels,1))*mean_scale
    if second_moments:
        moment_scale=np.sqrt(variance[2:].mean(1))
        scale[2:,0]=np.where(moment_scale>=1e-8,moment_scale,1.)
    if not all(np.isfinite(a).all() for a in (mean_y,scale_y,mean,scale)):
        raise ValueError('nonfinite TRAIN statistics')
    return dict(y_mean=mean_y,y_scale=scale_y,target_mean=mean,target_scale=scale)


def run(args):
    wall,cpu=time.monotonic(),time.process_time()
    torch.set_num_threads(1);torch.manual_seed(args.seed)
    rng=np.random.default_rng(args.seed)
    out=args.out
    if out.exists() and any(out.iterdir()):raise FileExistsError('fresh output destination required')
    out.mkdir(parents=True,exist_ok=True)
    cfg=json.loads(args.target.read_text());points=json.loads(args.geometry.read_text())['points']
    geometry=decoder_from_geometry(cfg,points);roi=np.asarray(geometry['roi'],int)
    tr,va=load_joint(args.train,geometry),load_joint(args.validation,geometry)
    if set(tr['parent'])&set(va['parent']):raise ValueError('TRAIN/validation parent overlap')
    if args.train_size:
        if not 2<=args.train_size<=len(tr['y']):raise ValueError('2<=train-size<=available required')
        tr={k:(v[:args.train_size] if isinstance(v,np.ndarray) else v) for k,v in tr.items()}
    moments=args.second_moments or args.kind.endswith('_moments')
    architecture=args.kind.replace('_moments','')
    if architecture=='r2_reference' and (args.reference_checkpoint is None or moments):
        raise ValueError('reference requires an old R2 MLP checkpoint and cannot invent second moments')
    stats=fit_statistics(tr,roi,moments,args.batch_size)
    device=torch.device(args.device)
    if args.device=='mps' and not torch.backends.mps.is_available():raise RuntimeError('MPS requested but unavailable')
    model=None;reference=None;reference_payload=None;reference_identity=None
    if architecture=='r2_reference':
        reference_identity=checkpoint_identity(args.reference_checkpoint)
        reference_payload=torch.load(args.reference_checkpoint,map_location='cpu',weights_only=False)
        reference=HistoricalR2Reference(reference_payload,args.device,reference_identity)
        if not np.array_equal(reference.reference.roi_indices,roi):raise ValueError('historical reference ROI mismatch')
        count=sum(p.numel() for p in reference.model.parameters())
    elif architecture=='train_mean':count=0
    else:
        model=FieldMomentRegressor(architecture,len(roi),moments).to(device)
        count=sum(p.numel() for p in model.parameters())
    if count>args.capacity_limit:raise ValueError('actual parameter count exceeds declared cap')
    warning=dict(moment_labels='same joint parent Re,Im,Re²,Im²,ReIm; cross moment is never product of predicted means',
        label_bias='raw moment labels target conditional expectations; no per-parent posterior variance labels exist from one joint draw',
        covariance_bias='predicted second moment minus predicted mean product is a nonlinear plug-in estimate; finite data/optimization can bias it',
        PSD_projection='per-pixel eigenvalue clipping changes the raw covariance estimate; it does not change the raw mean or prove calibration',
        scope='no cross-pixel covariance, normalized material density, ambient Gaussian posterior, joint material sampler, or credible interval qualification',
        mean_action='raw conditional-risk estimate may be outside a single decoder image; optional convex rectangle clipping remains a separate action')
    contract=dict(format=REFERENCE_FORMAT if reference_payload is not None else FORMAT,kind=args.kind,architecture=architecture,seed=args.seed,
        second_moments=moments,train=str(args.train.resolve()),validation=str(args.validation.resolve()),
        train_size=len(tr['y']),validation_size=len(va['y']),train_label_source=tr['label_source'],validation_label_source=va['label_source'],
        information_set='fixed empirical TRAIN mean, no online observation dependence' if architecture=='train_mean' else 'all raw512 source-major measurements only; no anchor/context/J/residual inputs',
        packing='8 source blocks, each Re32 then Im32',parameter_count=count,capacity_limit=args.capacity_limit,
        encoder_width=64 if architecture.endswith('small') else 128 if architecture.endswith('medium') else None,
        CNN_control='measurement2x32x8 convolutions plus raw512-to32 skip; same raw information, comparable reported parameter scale',
        criterion='validation full576 complex mean squared error, raw mean; update0 TRAIN mean eligible',
        mean_loss='mean of normalized squared Re/Im ROI errors with one common TRAIN scale, proportional to full576 complex MSE',
        moment_loss='mean separately TRAIN-scaled squared raw second-moment errors; fixed declared auxiliary weight',
        moment_loss_weight=args.moment_loss_weight,learning_rate=args.lr,updates_requested=args.updates,
        batch_size=args.batch_size,eval_every=args.eval_every,patience=args.patience,
        target_id=cfg['target_id'],outside_ROI='fixed original host, variance0',
        statistics='TRAIN only; input column scales, spatial target means, common Re/Im target scale, per-second-moment-channel scales',
        reference_identity=reference_identity,reference_historical_seed=reference_payload.get('seed') if reference_payload else None,
        reference_normalization='historical checkpoint used verbatim; new TRAIN statistics are unused diagnostics for reference' if reference_payload else None,
        reference_run_seed_ignored=bool(reference_payload is not None),
        posterior='NOT_CLAIMED',warnings=warning,physics_calls=0,device=args.device)
    write_json(out/'RUN_CONTRACT.json',contract)
    np.savez_compressed(out/'TRAIN_STATISTICS.npz',**stats)
    common=dict(format=FORMAT,kind=args.kind,architecture=architecture,seed=args.seed,
                second_moments=moments,geometry=geometry,**stats,contract=contract,warnings=warning,
                moment_source='EMPIRICAL_TRAIN_MOMENTS' if architecture=='train_mean' else 'SUPERVISED_RAW_MOMENT_REGRESSION')
    def payload(state,selected_update):
        if reference_payload is not None:
            return dict(format=REFERENCE_FORMAT,reference_payload=reference_payload,reference_identity=reference_identity,
                        geometry=geometry,contract=contract,selected_update=0)
        return dict(**common,state_dict=state,selected_update=selected_update)
    def inference(current_payload):
        return HistoricalR2Reference(reference_payload,args.device,reference_identity) if reference_payload is not None else FieldMomentPredictor(current_payload,args.device)
    def validation():
        if model is not None:model.eval()
        losses=[];mean_losses=[];moment_losses=[]
        for start in range(0,len(va['y']),args.validation_batch_size):
            end=start+args.validation_batch_size
            label=moment_labels(va['material'][start:end],roi,moments).astype(np.float64)
            if reference is not None:
                prediction=reference.predict(va['y'][start:end],args.validation_batch_size)
                predicted=np.stack((prediction[:,roi].real,prediction[:,roi].imag),1)
            elif model is None:
                predicted=np.broadcast_to(stats['target_mean'],label.shape)
                prediction=material_from_mean(predicted[:,:2],roi,cfg['host_chi'])
            else:
                yn=((va['y'][start:end]-stats['y_mean'])/stats['y_scale']).astype(np.float32)
                with torch.no_grad():
                    z=model(torch.from_numpy(yn).to(device)).cpu().numpy()
                predicted=z*stats['target_scale']+stats['target_mean']
                prediction=material_from_mean(predicted[:,:2],roi,cfg['host_chi'])
            loss=np.mean(abs(prediction-va['material'][start:end])**2,1)
            if not np.isfinite(loss).all():raise FloatingPointError('nonfinite full-grid validation MSE')
            losses.extend(loss)
            mean_losses.extend(((predicted[:,:2]-label[:,:2])**2).mean((1,2)))
            if moments:moment_losses.extend(((predicted[:,2:]-label[:,2:])**2).mean((1,2)))
        return np.asarray(losses),float(np.mean(mean_losses)),float(np.mean(moment_losses)) if moments else None
    optimizer=torch.optim.AdamW(model.parameters(),lr=args.lr,weight_decay=args.weight_decay) if model is not None else None
    selected=0;best=float('inf');best_state=None;history=[];update=0;training_loss=None;stale=0
    stop='UPDATE_LIMIT' if model is not None else 'FROZEN_REFERENCE' if reference is not None else 'TRAIN_MEAN_BASELINE'
    requested=[False];old={}
    for sig in (signal.SIGTERM,signal.SIGINT):old[sig]=signal.signal(sig,lambda *_:requested.__setitem__(0,True))
    def event():
        nonlocal best,selected,best_state,stale
        rows,component_mse,moment_mse=validation();value=float(rows.mean())
        row=dict(update=update,train_minibatch_normalized_loss=training_loss,
                 validation_full_complex_mse=value,validation_ROI_mean_component_mse=component_mse,
                 validation_raw_second_moment_mse=moment_mse,learning_rate=args.lr,wall_seconds=time.monotonic()-wall)
        history.append(row)
        with (out/'loss_curve.jsonl').open('a') as file:file.write(json.dumps(row,allow_nan=False)+'\n')
        state={k:v.detach().cpu().clone() for k,v in model.state_dict().items()} if model is not None else None
        if value<best:
            best,selected,best_state,stale=value,update,state,0
            torch.save(payload(best_state,selected),out/'BEST.pt')
            np.savez_compressed(out/'validation_metric.npz',parent=va['parent'],full_complex_mse=rows)
        else:stale+=1
        torch.save(dict(**payload(state,update),optimizer=optimizer.state_dict() if optimizer else None,
                        numpy_rng_state=rng.bit_generator.state,torch_rng_state=torch.get_rng_state()),out/'LAST.pt')
        print(json.dumps(dict(event='validation',**row,selected_update=selected)),flush=True)
    try:
        event()
        if model is not None:
            for update in range(1,args.updates+1):
                if requested[0] or (args.max_wall and time.monotonic()-wall>=args.max_wall):
                    update-=1;stop='SIGNAL' if requested[0] else 'WALL_LIMIT';break
                ids=rng.integers(0,len(tr['y']),min(args.batch_size,len(tr['y'])))
                yb=((tr['y'][ids]-stats['y_mean'])/stats['y_scale']).astype(np.float32)
                target=((moment_labels(tr['material'][ids],roi,moments)-stats['target_mean'])/stats['target_scale']).astype(np.float32)
                model.train();optimizer.zero_grad(set_to_none=True)
                pred=model(torch.from_numpy(yb).to(device));label=torch.from_numpy(target).to(device)
                loss=(pred[:,:2]-label[:,:2]).square().mean()
                if moments:loss=loss+args.moment_loss_weight*(pred[:,2:]-label[:,2:]).square().mean()
                if not bool(torch.isfinite(loss)):raise FloatingPointError('nonfinite training loss')
                loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5.,error_if_nonfinite=True);optimizer.step()
                training_loss=float(loss.detach().cpu())
                if update%args.eval_every==0 or update==args.updates:
                    event()
                    if args.patience and stale>=args.patience:stop='VALIDATION_PATIENCE';break
            if history[-1]['update']!=update:event()
    finally:
        for sig,handler in old.items():signal.signal(sig,handler)
    if model is not None:model.load_state_dict(best_state)
    predictor=inference(payload(best_state,selected))
    moments_out=predictor.predict_moments(va['y'],args.validation_batch_size)
    mean=moments_out['material_mean'];host=complex(*cfg['host_chi'])
    upper=host+complex(cfg['inclusion_increment_real_interval'][1],cfg['inclusion_increment_imag_interval'][1])
    clipped=np.clip(mean.real,host.real,upper.real)+1j*np.clip(mean.imag,host.imag,upper.imag)
    arrays=dict(parent=va['parent'],material_mean=mean,full_complex_mse=np.mean(abs(mean-va['material'])**2,1),
                clipped_material_mean=clipped,clipped_full_complex_mse=np.mean(abs(clipped-va['material'])**2,1))
    for key,value in moments_out.items():
        if isinstance(value,np.ndarray) and key!='material_mean':arrays[key]=value
    np.savez_compressed(out/'validation_predictions.npz',**arrays)
    projection=moments_out.get('psd_projection_warning',{})
    write_json(out/'COVARIANCE_WARNINGS.json',dict(**warning,variance_status=moments_out['second_moment_status'],
        PSD_projection_diagnostics=projection,
        mean_outside_legal_rectangle_count=int(((mean[:,roi].real<host.real)|(mean[:,roi].real>upper.real)|(mean[:,roi].imag<host.imag)|(mean[:,roi].imag>upper.imag)).sum())))
    np.savez_compressed(out/'loss_curve.npz',update=np.array([r['update'] for r in history]),validation_full_complex_mse=np.array([r['validation_full_complex_mse'] for r in history]))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7,4))
    ax.plot([r['update'] for r in history],[r['validation_full_complex_mse'] for r in history],marker='o')
    ax.axvline(selected,ls='--',color='gray')
    ax.set(xlabel='Optimizer update (0 eligible)',ylabel='Full576 complex material MSE',title=f'Raw512 {args.kind}; parameters {count:,}')
    fig.tight_layout();fig.savefig(out/'LOSS_CURVE.png',dpi=150);plt.close(fig)
    result=dict(**contract,selected_update=selected,selected_initialization=selected==0,completed_updates=update,
        selected_validation_full_complex_mse=best,validation_clipped_full_complex_mse=float(arrays['clipped_full_complex_mse'].mean()),
        stop_reason=stop,history=history,covariance_projection_diagnostics=projection,
        second_moment_status=moments_out['second_moment_status'],wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024))
    write_json(out/'RESULT.json',result)
    return result


def main(argv=None):
    root=Path(__file__).resolve().parents[2]
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('train','validation','out'):parser.add_argument('--'+name,required=True,type=Path)
    parser.add_argument('--kind',choices=KINDS,required=True);parser.add_argument('--seed',type=int,required=True)
    parser.add_argument('--updates',type=int,required=True);parser.add_argument('--device',choices=('cpu','mps'),default='cpu')
    parser.add_argument('--second-moments',action='store_true');parser.add_argument('--moment-loss-weight',type=float,default=1.)
    parser.add_argument('--target',type=Path,default=root/'research/p3l_r4_conditional_transport/configs/TARGET.json')
    parser.add_argument('--geometry',type=Path,default=root/'research/p1_material_target_20261018/inputs/geometry.json')
    parser.add_argument('--reference-checkpoint',type=Path);parser.add_argument('--train-size',type=int)
    parser.add_argument('--batch-size',type=int,default=128);parser.add_argument('--validation-batch-size',type=int,default=128)
    parser.add_argument('--eval-every',type=int,default=100);parser.add_argument('--patience',type=int,default=12)
    parser.add_argument('--lr','--learning-rate',type=float,default=5e-4);parser.add_argument('--weight-decay',type=float,default=1e-4)
    parser.add_argument('--capacity-limit',type=int,default=1_000_000);parser.add_argument('--max-wall','--max-seconds',type=float,default=0.)
    args=parser.parse_args(argv)
    if args.updates<0 or args.patience<0 or min(args.batch_size,args.validation_batch_size,args.eval_every,args.capacity_limit)<1:
        parser.error('nonnegative updates/patience and positive batches/evaluation interval/cap required')
    if not np.isfinite([args.lr,args.weight_decay,args.max_wall,args.moment_loss_weight]).all() or args.lr<=0 or args.weight_decay<0 or args.max_wall<0 or args.moment_loss_weight<=0:
        parser.error('finite positive learning rate/moment weight; nonnegative decay/wall cap required')
    run(args)


if __name__=='__main__':main()
