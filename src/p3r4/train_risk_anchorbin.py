"""Revision2: same exact-bin loss/optimizer/loop with explicit anchor mass tilt.

Parent schedules training; no inference, data generation, or physics is launched
beyond this one explicitly requested fit. Old trainers/checkpoints are untouched.
"""
import argparse,copy,json,time,resource,sys
from pathlib import Path
import numpy as np
import torch
from .risk_anchorbin import AnchorBinRisk,FORMAT,EPSILON,FRACTION_TOLERANCE
from .risk_bin import bin_phase_labels,MOMENT_RULE
from .risk import feature_image,roi_components,full_material
from .train_risk import load_joint,align
from p3r2.residual_risk import read_anchor,fit_normalization,online_features


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('train','validation','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--anchor-train','--train-context',dest='train_context',type=Path,required=True)
    p.add_argument('--anchor-validation','--validation-context',dest='validation_context',type=Path,required=True)
    p.add_argument('--kind',choices=('phase8','phase16'),required=True);p.add_argument('--seed',type=int,default=2026102501)
    p.add_argument('--updates',type=int,default=3000);p.add_argument('--batch-size',type=int,default=64);p.add_argument('--learning-rate',type=float,default=5e-4)
    p.add_argument('--device',choices=('cpu','mps'),default='cpu');p.add_argument('--context-mode',choices=('spectral','decoder_only','random'),default='spectral')
    default=Path(__file__).resolve().parents[3]/'b1-p3l-r2-risk-posterior/research/p3l_r2_risk_posterior_exploration/A/residual_training/n512/local_geometry_cnn_s33/BEST.pt'
    p.add_argument('--geometry-checkpoint',type=Path,default=default);a=p.parse_args(argv)
    if not 0<=a.updates<=3000 or a.batch_size!=64 or a.learning_rate!=5e-4:raise ValueError('frozen max3000updates,batch64,lr5e-4; update0 allowed for interface checks')
    if a.out.exists():raise FileExistsError('new revision output only')
    wall,cpu=time.monotonic(),time.process_time();torch.set_num_threads(1);torch.manual_seed(a.seed);rng=np.random.default_rng(a.seed)
    tr,va=load_joint(a.train),load_joint(a.validation)
    if set(tr['parent'])&set(va['parent']):raise ValueError('parent overlap')
    tc,vc=read_anchor(a.train_context),read_anchor(a.validation_context);tr,va=align(tr,tc),align(va,vc)
    saved=torch.load(a.geometry_checkpoint,map_location='cpu',weights_only=False);g=saved['geometry'];normalization=saved['normalization']
    if a.context_mode=='random':
        g=copy.deepcopy(g);frame,_=np.linalg.qr(np.random.default_rng(2026102504).normal(size=(64,64)));w=g['weak'].shape[1];s=g['strong'].shape[1]
        g['weak']=frame[:,:w].astype('float32');g['strong']=frame[:,w:w+s].astype('float32');normalization=fit_normalization(online_features(tc,g))
    ti=feature_image(tr['y'],g,tc,normalization);vi=feature_image(va['y'],g,vc,normalization)
    if a.context_mode=='decoder_only':ti[:,2:6]=0;vi[:,2:6]=0
    target=roi_components(tr['material'],g['roi']);vt=roi_components(va['material'],g['roi'])
    model=AnchorBinRisk(a.kind,g).to(a.device);count=sum(v.numel() for v in model.parameters());assert count==(6114 if a.kind=='phase8' else 6386)
    optimizer=torch.optim.AdamW(model.parameters(),lr=a.learning_rate,weight_decay=1e-3)
    thresholds=model.thresholds.detach().cpu().numpy()
    def prediction(image,online,ids,return_moments=False):
        return model(torch.tensor(image[ids],device=a.device),torch.tensor(online['q0'][ids],dtype=torch.float32,device=a.device),torch.tensor(np.stack((online['chi0'][ids].real,online['chi0'][ids].imag),1),dtype=torch.float32,device=a.device),return_moments=return_moments)
    def evaluate():
        model.eval();total=0.;rows=[]
        with torch.inference_mode():
            for start in range(0,len(vi),a.batch_size):
                ids=slice(start,start+a.batch_size);pred=prediction(vi,vc,ids).cpu().numpy();rows.append(pred);total+=float(np.sum((pred-vt[ids])**2))
        return total/vt.size,np.concatenate(rows)
    score,vp=evaluate();best=score;state=copy.deepcopy(model.state_dict());best_update=0;history=[dict(update=0,validation_component_mse=score)]
    update0=full_material(vp,g);update0_anchor_rms=float(np.sqrt(np.mean(abs(update0-vc['chi0'])**2)))
    for update in range(1,a.updates+1):
        ids=rng.integers(len(ti),size=min(a.batch_size,len(ti)));model.train();optimizer.zero_grad(set_to_none=True)
        _,mom=prediction(ti,tc,ids,True);labels=bin_phase_labels(tr['x'][ids],g,thresholds,None if tr['s'] is None else tr['s'][ids])
        loss=(((mom-torch.tensor(labels,device=a.device))**2)*model.weights[None,None,:,None]).sum(2).mean()
        if not torch.isfinite(loss):raise RuntimeError('nonfinite loss')
        loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5.);optimizer.step()
        if update%max(1,min(100,a.updates))==0 or update==a.updates:
            score,current=evaluate();history.append(dict(update=update,validation_component_mse=score))
            if score<best:best=score;state=copy.deepcopy(model.state_dict());best_update=update;vp=current
    a.out.mkdir(parents=True)
    common=dict(format=FORMAT,revision='FINAL_PHASE_REVISION2_ANCHOR_INITIALIZATION',kind=a.kind,seed=a.seed,moment_rule=MOMENT_RULE,
                anchor_epsilon=EPSILON,anchor_fraction_tolerance=FRACTION_TOLERANCE,geometry=g,feature_normalization=normalization,context_mode=a.context_mode,
                state_dict={k:v.cpu() for k,v in state.items()},training_parents=tr['parent'],selected_update=best_update,
                online_contract='explicit native q0 and raw complex chi0; no GT or normalized feature ratios',
                initialization='log((1-eps)*anchor joint-bin masses+eps/(K+1)) plus zero initialized full mass correction',
                stopping_rule='final substantive phase revision; if validation gain absent stop phase family; parent owns decision',
                physics_calls=0,conditional_density='NOT_CLAIMED')
    torch.save(common,a.out/'BEST.pt');pred=full_material(vp,g);host=complex(*g['config']['host_chi']);upper=host+complex(g['config']['inclusion_increment_real_interval'][1],g['config']['inclusion_increment_imag_interval'][1])
    clipped=np.clip(pred.real,host.real,upper.real)+1j*np.clip(pred.imag,host.imag,upper.imag)
    risk=np.mean(abs(pred-va['material'])**2,1);crisk=np.mean(abs(clipped-va['material'])**2,1)
    np.savez_compressed(a.out/'validation.npz',parent=va['parent'],prediction=pred,risk=risk,clipped_prediction=clipped,clipped_risk=crisk)
    result=dict(format=FORMAT,kind=a.kind,seed=a.seed,parameter_count=count,updates_requested=a.updates,selected_update=best_update,
        train_size=len(ti),validation_size=len(vi),validation_full_complex_mse=float(risk.mean()),validation_clipped_full_complex_mse=float(crisk.mean()),
        update0_anchor_full_RMS=update0_anchor_rms,anchor_epsilon=EPSILON,moment_rule=MOMENT_RULE,
        learning_rate=a.learning_rate,batch_size=a.batch_size,weight_decay=1e-3,gradient_clip_norm=5.,history=history,
        wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024),physics_calls=0)
    (a.out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
