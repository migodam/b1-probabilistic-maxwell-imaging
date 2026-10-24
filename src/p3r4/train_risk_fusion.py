"""Six parent-scheduled fits: two frozen moment contexts x three fixed seeds."""
import argparse,json,copy,time,resource,sys
from pathlib import Path
import numpy as np,torch
from .risk_fusion import *
from .risk import feature_image,roi_components,full_material
from .train_risk import align
from p3r2.residual_risk import read_anchor,fit_normalization,online_features
from p3r2.risk_models import observation_features


def material_data(path):
    if 'locked' in str(path).lower():raise ValueError('locked training labels prohibited')
    with np.load(path,allow_pickle=False) as z:
        parent=z['parent'].astype(str);y=z['y'] if 'y' in z else observation_features(z['observation']);material=z['material'] if 'material' in z else z['chi']
    if material.shape!=(len(parent),576) or y.shape!=(len(parent),512) or len(set(parent))!=len(parent):raise ValueError('joint material/y/independentparents')
    return dict(parent=parent,y=y,material=material)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('train','validation','out','moments-train','moments-validation','anchor-train','anchor-validation'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--control',choices=('moments','permuted'),required=True);p.add_argument('--seed',type=int,choices=(2026102501,2026102502,2026102503),default=2026102501)
    p.add_argument('--updates',type=int,default=4000);p.add_argument('--device',choices=('cpu','mps'),default='cpu')
    p.add_argument('--refit-feature-normalization',action='store_true')
    default=Path(__file__).resolve().parents[3]/'b1-p3l-r2-risk-posterior/research/p3l_r2_risk_posterior_exploration/A/residual_training/n512/local_geometry_cnn_s33/BEST.pt';p.add_argument('--geometry-checkpoint',type=Path,default=default);a=p.parse_args(argv)
    if a.out.exists() or not 0<=a.updates<=4000:raise ValueError('freshoutput,max4000updates')
    wall,cpu=time.monotonic(),time.process_time();torch.set_num_threads(1);torch.manual_seed(a.seed);rng=np.random.default_rng(a.seed)
    tr,va=material_data(a.train),material_data(a.validation)
    if set(tr['parent'])&set(va['parent']):raise ValueError('foldoverlap')
    tc,vc=read_anchor(a.anchor_train),read_anchor(a.anchor_validation);tr,va=align(tr,tc),align(va,vc)
    saved=torch.load(a.geometry_checkpoint,map_location='cpu',weights_only=False);g=saved['geometry'];normalization=saved['normalization']
    if a.refit_feature_normalization:normalization=fit_normalization(online_features(tc,g))
    gate_t=json.loads((a.moments_train/'MOMENTS_COMPLETE.json').read_text());gate_v=json.loads((a.moments_validation/'MOMENTS_COMPLETE.json').read_text());identity=gate_t['checkpoint_identity']
    provider,moment_count,identity=matching_moment_gates(gate_t,gate_v)
    if provider!='spline' and not a.refit_feature_normalization:raise ValueError('strong provider controls require --refit-feature-normalization, TRAIN-only online scaling')
    if provider!='spline' and (len(tc['parent']),len(vc['parent']))!=(1024,256):raise ValueError('strong provider controls freeze exactly1024TRAIN/256validation contexts')
    raw_t=read_moments(a.moments_train,tr['parent'],identity,provider,moment_count);raw_v=read_moments(a.moments_validation,va['parent'],identity,provider,moment_count)
    mean=raw_t.mean((0,2,3),keepdims=True,dtype=np.float64).astype('float32');scale=raw_t.std((0,2,3),keepdims=True,dtype=np.float64).astype('float32');scale[scale<1e-8]=1
    mt=transform_moments(raw_t,mean,scale,tr['parent'],a.control);mv=transform_moments(raw_v,mean,scale,va['parent'],a.control)
    ti=feature_image(tr['y'],g,tc,normalization);vi=feature_image(va['y'],g,vc,normalization);ta=roi_components(tc['chi0'],g['roi']);vanc=roi_components(vc['chi0'],g['roi'])
    target=roi_components(tr['material'],g['roi']);vt=roi_components(va['material'],g['roi']);model=FusionRisk(g).to(a.device);count=sum(z.numel() for z in model.parameters());assert count==6418
    optimizer=torch.optim.AdamW(model.parameters(),lr=5e-4,weight_decay=1e-3)
    def prediction(image,moments,base,ids):return model(torch.tensor(image[ids],device=a.device),torch.tensor(moments[ids],device=a.device),torch.tensor(base[ids],device=a.device))
    def evaluate():
        model.eval();rows=[];total=0.
        with torch.inference_mode():
            for start in range(0,len(vi),64):
                sl=slice(start,start+64);pred=prediction(vi,mv,vanc,sl).cpu().numpy();rows.append(pred);total+=float(np.sum((pred-vt[sl])**2))
        return total/vt.size,np.concatenate(rows)
    best,vp=evaluate();state=copy.deepcopy(model.state_dict());best_update=0;history=[dict(update=0,validation_component_mse=best)]
    for update in range(1,a.updates+1):
        ids=rng.integers(len(ti),size=min(64,len(ti)));model.train();optimizer.zero_grad(set_to_none=True)
        loss=((prediction(ti,mt,ta,ids)-torch.tensor(target[ids],device=a.device))**2).mean()
        if not torch.isfinite(loss):raise RuntimeError('nonfinite material MSE')
        loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5.);optimizer.step()
        if update%100==0 or update==a.updates:
            score,current=evaluate();history.append(dict(update=update,validation_component_mse=score))
            if score<best:best=score;best_update=update;state=copy.deepcopy(model.state_dict());vp=current
    a.out.mkdir(parents=True)
    payload=dict(format=FORMAT,kind='fusion',control=a.control,seed=a.seed,state_dict={k:v.cpu() for k,v in state.items()},geometry=g,feature_normalization=normalization,
        moment_mean=mean,moment_scale=scale,moment_checkpoint_identity=identity,moment_K=moment_count,moment_provider=provider,moment_contract_version=2,permutation_seed=PERMUTATION_SEED,selected_update=best_update,
        training_parents=tr['parent'],moment_semantics='frozen provider moments as coarse context; model pushforward for explicit densities, deterministic learned pointwise joint ReIm moments forfield; not true posterior/no independence claim',
        optimizer=dict(lr=5e-4,weight_decay=1e-3,batch=64,clip=5,updates=a.updates),no_GT_features=True,physics_calls=0)
    torch.save(payload,a.out/'BEST.pt');pred=full_material(vp,g)
    from p3r2.fusion import passive_box
    clipped=passive_box(pred,g['config'],g['roi']);risks=np.mean(abs(pred-va['material'])**2,1);crisk=np.mean(abs(clipped-va['material'])**2,1)
    np.savez_compressed(a.out/'validation.npz',parent=va['parent'],prediction=pred,clipped_prediction=clipped,risk=risks,clipped_risk=crisk)
    result=dict(format=FORMAT,moment_provider=provider,moment_K=moment_count,refit_feature_normalization=a.refit_feature_normalization,control=a.control,seed=a.seed,parameter_count=count,train_size=len(ti),validation_size=len(vi),updates_requested=a.updates,selected_update=best_update,
        validation_full_complex_mse=float(risks.mean()),validation_clipped_full_complex_mse=float(crisk.mean()),history=history,
        moment_scope=payload['moment_semantics'],stop_rule='stop if development no gain; no locked tuning',physics_calls=0,
        wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024))
    (a.out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
