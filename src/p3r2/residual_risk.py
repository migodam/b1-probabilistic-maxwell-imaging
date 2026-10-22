"""Anchor-aware material-risk estimators; training-only fixed geometry frame.

No physical forward or adjoint is called. Inference consumes saved online anchor
observations, predictions/residual, q0, chi0 and J, never q_GT or chi_GT.
"""
import copy
import json
import time
from pathlib import Path
import numpy as np
import torch
from torch import nn
from p1_material.prior import MaterialPrior, saturation_prime
from .risk_models import observation_features
from .risk_training import material_targets, risk_rows, risk_summary, write_json

MODELS = ('anchor', 'residual_ridge', 'residual_mlp', 'geometry_mlp', 'physics_cnn',
          'weak_decoder', 'weak_strong_decoder', 'local_anchor_cnn', 'local_geometry_cnn')


def read_anchor(path, supervised=False):
    if supervised and 'locked' in str(path).lower():
        raise ValueError('locked truth is prohibited in supervised training')
    online = ('observation','prediction','residual','q0','chi0','J','parent')
    with np.load(path, allow_pickle=False) as z:
        data = {key:z[key] for key in online if key in z}
        if 'parent' not in data and 'parents' in z:
            data['parent'] = z['parents']
        if supervised:
            data['chi_GT'] = z['chi_GT']
    if 'residual' not in data and 'prediction' in data:
        data['residual'] = data['observation']-data['prediction']
    for key in ('observation','residual','q0','chi0','J','parent'):
        if key not in data:
            raise ValueError(f'missing anchor field {key}')
    n=len(data['q0'])
    shapes={'observation':(n,32,8),'residual':(n,32,8),'q0':(n,66),
            'chi0':(n,576),'J':(n,512,66)}
    for key, shape in shapes.items():
        if data[key].shape!=shape or not np.isfinite(data[key]).all():
            raise ValueError(f'invalid {key}: expected finite {shape}')
    data['parent']=data['parent'].astype(str)
    if len(set(data['parent']))!=n:
        raise ValueError('independent unique parents required')
    if supervised and (data['chi_GT'].shape!=(n,576) or not np.isfinite(data['chi_GT']).all()):
        raise ValueError('invalid material supervision')
    return data


def subset(data, n):
    return {k:v[:n] for k,v in data.items()}


def fit_geometry(train, config, points, weak_rank=16, strong_rank=12, j_is_whitened=True):
    """Frame uses train J only. Singular directions are z64, with u2 separate."""
    prior=MaterialPrior(config,points)
    sigma=float(config['noise']['sigma'])
    scale=1. if j_is_whitened else np.sqrt(2.)/sigma
    gram=np.zeros((64,64),dtype=np.float64)
    for begin in range(0,len(train['J']),32):
        j=np.asarray(train['J'][begin:begin+32,:,:64],dtype=np.float64)*scale
        gram+=np.einsum('nri,nrj->ij',j,j,optimize=True)
    gram/=len(train['J'])
    val,vec=np.linalg.eigh(gram)
    # Canonicalize signs for reproducible saved frames across eigensolver runs.
    vec*=np.where(vec[np.argmax(np.abs(vec),axis=0),np.arange(64)]<0,-1.,1.)[None,:]
    return dict(A=prior.A.astype('float32'), roi=np.flatnonzero(prior.mask),
                weak=vec[:,:weak_rank].astype('float32'), strong=vec[:,-strong_rank:].astype('float32'),
                eigenvalues=val, config=config, j_is_whitened=j_is_whitened,
                frame_fit='train J only; ascending eigenspaces of E[Jz.T Jz]; no material or latent GT')


def online_features(data, geometry):
    """All inputs are online fields; chi_GT/q_GT are intentionally inaccessible."""
    obs=observation_features(data['observation'])
    residual=observation_features(data['residual'])
    q=np.asarray(data['q0'],dtype='float32')
    cfg=geometry['config']; sigma=float(cfg['noise']['sigma'])
    rw=residual*(np.sqrt(2.)/sigma)
    base=np.concatenate((obs,rw,q),axis=1)
    j=np.asarray(data['J'],dtype='float32')
    if not geometry['j_is_whitened']:
        j=j*(np.sqrt(2.)/sigma)
    grad=np.einsum('nri,nr->ni',j,rw,optimize=True)
    # Per-column norm preserves the original coordinate convention. A complete
    # per-anchor SVD would be expensive and its sign arbitrary; no SVD is used.
    norms=np.log1p(np.sum(j*j,axis=1))
    weak=np.asarray(geometry['weak']); strong=np.asarray(geometry['strong'])
    gweak=grad[:,:64]@weak; gstrong=grad[:,:64]@strong
    extra=np.concatenate((grad,norms,gweak,gstrong),axis=1)
    images=[]
    A=np.asarray(geometry['A'],dtype='float32'); roi=geometry['roi']
    mask=np.zeros(576,dtype='float32'); mask[roi]=1
    for n in range(len(q)):
        f=cfg['field_mean']+A@q[n,:64]
        df=saturation_prime(f/cfg['phase_transition_halfwidth'])/cfg['phase_transition_halfwidth']*mask
        delta=complex(cfg['inclusion_increment_real_interval'][0]+
                      np.diff(cfg['inclusion_increment_real_interval'])[0]*q[n,64],
                      cfg['inclusion_increment_imag_interval'][0]+
                      np.diff(cfg['inclusion_increment_imag_interval'])[0]*q[n,65])
        dw=(df[:,None]*(A@weak))*abs(delta)
        ds=(df[:,None]*(A@strong))*abs(delta)
        # Normalize the latent physical-gradient direction before mapping it to
        # material so sensor units do not overwhelm the image channels.
        direction=grad[n,:64]/max(np.linalg.norm(grad[n,:64]),1.)
        back=delta*df*(A@direction)
        anchor=data['chi0'][n]
        images.append(np.stack((anchor.real,anchor.imag,np.sqrt(np.sum(dw*dw,axis=1)),
                                np.sqrt(np.sum(ds*ds,axis=1)),back.real,back.imag)).reshape(6,24,24))
    return dict(base=base.astype('float32'), geometry=np.concatenate((base,extra),axis=1).astype('float32'),
                image=np.asarray(images,dtype='float32'))


def fit_normalization(features):
    result={}
    for key,arr in features.items():
        axes=(0,2,3) if key=='image' else (0,)
        mean=arr.mean(axis=axes,keepdims=True,dtype=np.float64).astype('float32')
        scale=arr.std(axis=axes,keepdims=True,dtype=np.float64).astype('float32')
        scale[scale<1e-10]=1
        result[key+'_mean']=mean; result[key+'_scale']=scale
    return result


def normalize(features, stats):
    return {key:((a-stats[key+'_mean'])/stats[key+'_scale']).astype('float32') for key,a in features.items()}


class ResidualNet(nn.Module):
    def __init__(self, kind, dimensions, geometry, residual_scale, width=None):
        super().__init__(); self.kind=kind
        self.residual_scale=float(residual_scale)
        self.config=geometry['config']
        self.register_buffer('A',torch.as_tensor(geometry['A']))
        self.register_buffer('roi',torch.as_tensor(geometry['roi'],dtype=torch.long))
        self.register_buffer('weak',torch.as_tensor(geometry['weak']))
        self.register_buffer('strong',torch.as_tensor(geometry['strong']))
        output=2*len(geometry['roi'])
        if kind in ('local_anchor_cnn','local_geometry_cnn'):
            # Shared local filters impose a translation-shared material-risk bias.
            # Anchor-only: chi0 Re/Im and broadcast u2. Geometry: those inputs plus
            # weak/strong sensitivity and complex residual-backprojection images.
            channels=4 if kind=='local_anchor_cnn' else 8
            self.local=nn.Sequential(nn.Conv2d(channels,16,3,padding=1),nn.GELU(),
                       nn.Conv2d(16,16,3,padding=1),nn.GELU(),
                       nn.Conv2d(16,16,3,padding=1),nn.GELU())
            self.head=nn.Conv2d(16,2,1)
        elif kind=='physics_cnn':
            self.image=nn.Sequential(nn.Conv2d(6,16,3,padding=1),nn.GELU(),
                       nn.Conv2d(16,32,3,stride=2,padding=1),nn.GELU(),
                       nn.Conv2d(32,32,3,stride=2,padding=1),nn.GELU(),
                       nn.AdaptiveAvgPool2d((4,4)),nn.Flatten())
            self.scalar=nn.Sequential(nn.Linear(dimensions['geometry'],128),nn.GELU())
            self.head=nn.Sequential(nn.Linear(640,128),nn.GELU(),nn.Linear(128,output))
        else:
            ins=dimensions['base'] if kind=='residual_mlp' else dimensions['geometry']
            if kind=='geometry_mlp':
                ins+=6*24*24
            if kind=='weak_decoder': output=self.weak.shape[1]+2
            if kind=='weak_strong_decoder': output=self.weak.shape[1]+self.strong.shape[1]+2
            self.head=nn.Sequential(nn.Linear(ins,width or 128),nn.GELU(),nn.Dropout(.05),
                                   nn.Linear(width or 128,128),nn.GELU(),nn.Linear(128,output))
        last=self.head if isinstance(self.head,nn.Conv2d) else self.head[-1]
        nn.init.zeros_(last.weight); nn.init.zeros_(last.bias)

    def decode(self,q):
        cfg=self.config
        field=cfg['field_mean']+q[:,:64]@self.A.T
        t=torch.clamp((field[:,self.roi]/cfg['phase_transition_halfwidth']+1)/2,0,1)
        fraction=t*t*t*(10+t*(-15+6*t))
        real=cfg['inclusion_increment_real_interval']; imag=cfg['inclusion_increment_imag_interval']
        dr=real[0]+(real[1]-real[0])*q[:,64:65]
        di=imag[0]+(imag[1]-imag[0])*q[:,65:66]
        return torch.cat((cfg['host_chi'][0]+dr*fraction,cfg['host_chi'][1]+di*fraction),dim=1)

    def forward(self,base,geometry,image,q0,anchor):
        if self.kind in ('local_anchor_cnn','local_geometry_cnn'):
            material=image[:,:2] if self.kind=='local_anchor_cnn' else image
            # The last two base features are u2, standardized from train only.
            properties=base[:,-2:,None,None].expand(-1,-1,24,24)
            local_image=torch.cat((material,properties),dim=1)
            grid=self.head(self.local(local_image)).flatten(2)
            raw=grid[:,:,self.roi].reshape(len(base),-1)
        elif self.kind=='physics_cnn':
            raw=self.head(torch.cat((self.scalar(geometry),self.image(image)),dim=1))
        elif self.kind=='residual_mlp': raw=self.head(base)
        elif self.kind=='geometry_mlp': raw=self.head(torch.cat((geometry,image.flatten(1)),dim=1))
        else: raw=self.head(geometry)
        if 'decoder' not in self.kind:
            return anchor+raw*self.residual_scale
        w=self.weak.shape[1]
        dz=3*torch.tanh(raw[:,:w])@self.weak.T
        offset=w
        if self.kind=='weak_strong_decoder':
            s=self.strong.shape[1]
            # Smoothly bounded smaller strong correction; factorized weak/strong
            # heads share context but do not assert probabilistic independence.
            dz=dz+.75*torch.tanh(raw[:,w:w+s])@self.strong.T
            offset+=s
        odds_factor=torch.exp(2*torch.tanh(raw[:,offset:offset+2]))
        u=q0[:,64:]*odds_factor/(1-q0[:,64:]+q0[:,64:]*odds_factor)
        q=torch.cat((q0[:,:64]+dz,u),dim=1)
        return self.decode(q)


def parameter_count(net):
    return sum(p.numel() for p in net.parameters())


def matched_width(input_dim,target_parameters,output_dim=512):
    # Two hidden layers W,128, with fixed output; integer rounding is recorded.
    return max(8,int(round((target_parameters-128-(128+1)*output_dim)/(input_dim+1+128))))


class ResidualPredictor:
    def __init__(self,payload,device='cpu'):
        self.p=payload; self.device=torch.device(device)
        self.model=None
        if payload['kind'] not in ('anchor','residual_ridge'):
            self.model=ResidualNet(payload['kind'],payload['dimensions'],payload['geometry'],
                                   payload['residual_scale'],payload['width']).to(self.device)
            self.model.load_state_dict(payload['state_dict']); self.model.eval()

    def predict(self,online,batch_size=128):
        # Construct a strict online-only dictionary even if caller supplies labels.
        fields={k:online[k] for k in ('observation','residual','q0','chi0','J')}
        roi=self.p['geometry']['roi']; predictions=[]
        with torch.inference_mode():
            for begin in range(0,len(fields['q0']),batch_size):
                batch={k:v[begin:begin+batch_size] for k,v in fields.items()}
                anchor=material_targets(batch['chi0'],roi)
                if self.p['kind']=='anchor': y=anchor
                else:
                    f=normalize(online_features(batch,self.p['geometry']),self.p['normalization'])
                    if self.p['kind']=='residual_ridge': y=anchor+f['base']@self.p['coef']
                    else:
                        args=[torch.from_numpy(f[k]).to(self.device) for k in ('base','geometry','image')]
                        args += [torch.as_tensor(batch['q0'],dtype=torch.float32,device=self.device),
                                 torch.from_numpy(anchor).to(self.device)]
                        y=self.model(*args).cpu().numpy()
                c=len(roi); predictions.append(y[:,:c]+1j*y[:,c:])
        return np.concatenate(predictions)

    def predict_full(self,online,batch_size=128):
        out=np.asarray(online['chi0']).copy()
        out[:,self.p['geometry']['roi']]=self.predict(online,batch_size)
        return out


def load_residual_predictor(path,device='cpu'):
    return ResidualPredictor(torch.load(path,map_location='cpu',weights_only=False),device)


def train_residual(kind,seed,train,validation,geometry,output,epochs=100,batch_size=128,
                   patience=15,max_seconds=600,learning_rate=5e-4,weight_decay=1e-3,device='cpu'):
    torch.set_num_threads(1)
    if torch.get_num_interop_threads()!=1: torch.set_num_interop_threads(1)
    torch.manual_seed(seed); rng=np.random.default_rng(seed)
    out=Path(output); out.mkdir(parents=True,exist_ok=True)
    (out/'epochs.jsonl').write_text('')
    wall,cpu=time.monotonic(),time.process_time()
    tf,vf=online_features(train,geometry),online_features(validation,geometry)
    stats=fit_normalization(tf); tf,vf=normalize(tf,stats),normalize(vf,stats)
    roi=geometry['roi']; y=material_targets(train['chi_GT'],roi); vy=material_targets(validation['chi_GT'],roi)
    anchor=material_targets(train['chi0'],roi); va=material_targets(validation['chi0'],roi)
    residual_scale=float(max(np.sqrt(np.mean((y-anchor)**2)),1e-8))
    dimensions={k:tf[k].shape[1] for k in ('base','geometry')}
    comparator=ResidualNet('physics_cnn',dimensions,geometry,residual_scale)
    target_count=parameter_count(comparator); del comparator
    ins=dimensions['base'] if kind=='residual_mlp' else dimensions['geometry']+6*24*24
    width=matched_width(ins,target_count) if kind in ('residual_mlp','geometry_mlp') else 128
    payload=dict(format='p3r2_anchor_material_risk_v1',kind=kind,seed=seed,geometry=geometry,
                 normalization=stats,dimensions=dimensions,width=width,residual_scale=residual_scale,
                 training_parents=train['parent'],selection='validation complex ROI MSE; update zero included',
                 estimator_type='legal P1 decoder point' if 'decoder' in kind else 'unrestricted material mean point',
                 calibrated_density='NOT_CLAIMED')
    history=[]; best_update=0; stop='BASELINE'; model=None
    if kind=='residual_ridge':
        x=tf['base'].astype('float64'); vx=vf['base']; target=(y-anchor).astype('float64')
        # Dual ridge is bounded by the number of parents, suitable for n=512.
        gram=x@x.T/len(x); eig,vec=np.linalg.eigh(gram); rotated=vec.T@target
        best=float(np.mean((va-vy)**2)); payload['coef']=np.zeros((x.shape[1],y.shape[1]),dtype='float32')
        payload['alpha']=None
        history.append(dict(alpha='anchor_zero_correction',validation_complex_roi_mse=2*best))
        for alpha in (.001,.01,.1,1.,10.,100.):
            coef=(x.T@(vec@(rotated/(eig[:,None]+alpha)))/len(x)).astype('float32')
            score=float(np.mean((va+vx@coef-vy)**2))
            history.append(dict(alpha=alpha,validation_complex_roi_mse=2*score))
            if score<best: best=score; payload['coef']=coef.copy(); payload['alpha']=alpha
        stop='VALIDATION_RIDGE_ALPHA_WITH_ANCHOR'
    elif kind!='anchor':
        dev=torch.device(device)
        model=ResidualNet(kind,dimensions,geometry,residual_scale,width).to(dev)
        optimizer=torch.optim.AdamW(model.parameters(),lr=learning_rate,weight_decay=weight_decay)
        scheduler=torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer,factor=.5,patience=5,min_lr=1e-6)
        def batch(features,data,a,indices):
            return [torch.from_numpy(features[k][indices]).to(dev) for k in ('base','geometry','image')]+[
                torch.as_tensor(data['q0'][indices],dtype=torch.float32,device=dev),torch.from_numpy(a[indices]).to(dev)]
        def evaluate(features,data,a,target):
            model.eval(); total=0
            with torch.inference_mode():
                for begin in range(0,len(a),batch_size):
                    ix=slice(begin,begin+batch_size)
                    pred=model(*batch(features,data,a,ix))
                    total+=float(torch.sum((pred-torch.from_numpy(target[ix]).to(dev))**2).cpu())
            return total/target.size
        best=evaluate(vf,validation,va,vy); state=copy.deepcopy(model.state_dict())
        history.append(dict(epoch=0,update=0,validation_complex_roi_mse=2*best,
                            train_complex_roi_mse=2*evaluate(tf,train,anchor,y),wall_seconds=time.monotonic()-wall))
        with (out/'epochs.jsonl').open('a') as f: f.write(json.dumps(history[-1])+'\n')
        stale=0; update=0; stop='EPOCH_LIMIT'
        for epoch in range(1,epochs+1):
            if time.monotonic()-wall>=max_seconds: stop='TIME_LIMIT'; break
            model.train(); order=rng.permutation(len(y)); total=0.; seen=0; interrupted=False
            for begin in range(0,len(y),batch_size):
                if time.monotonic()-wall>=max_seconds: interrupted=True; break
                ix=order[begin:begin+batch_size]
                optimizer.zero_grad(set_to_none=True)
                pred=model(*batch(tf,train,anchor,ix))
                loss=nn.functional.mse_loss(pred,torch.from_numpy(y[ix]).to(dev))/(residual_scale**2)
                if not torch.isfinite(loss): raise RuntimeError('nonfinite material loss')
                loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(),5.); optimizer.step()
                total+=float(loss.detach().cpu())*len(ix); seen+=len(ix); update+=1
            val=evaluate(vf,validation,va,vy); scheduler.step(val)
            row=dict(epoch=epoch,update=update,completed_epoch=not interrupted,
                     validation_complex_roi_mse=2*val,train_minibatch_complex_roi_mse=2*total/max(seen,1)*residual_scale**2,
                     learning_rate=optimizer.param_groups[0]['lr'],wall_seconds=time.monotonic()-wall)
            history.append(row)
            with (out/'epochs.jsonl').open('a') as f: f.write(json.dumps(row)+'\n')
            if val<best: best=val; state=copy.deepcopy(model.state_dict()); best_update=update; stale=0
            else: stale+=1
            if interrupted: stop='TIME_LIMIT'; break
            if stale>=patience: stop='VALIDATION_PATIENCE'; break
        payload['state_dict']={k:v.detach().cpu() for k,v in state.items()}
    payload['best_update']=best_update
    torch.save(payload,out/'BEST.pt')
    predictor=ResidualPredictor(payload,device)
    pred=predictor.predict(validation,batch_size); train_pred=predictor.predict(train,batch_size)
    risks=risk_rows(pred,validation['chi_GT'][:,roi]); tr=risk_rows(train_pred,train['chi_GT'][:,roi])
    anchor_risks=risk_rows(validation['chi0'][:,roi],validation['chi_GT'][:,roi])
    np.savez_compressed(out/'validation_risks.npz',parent=validation['parent'],risk=risks,
                        anchor_risk=anchor_risks,prediction_roi=pred,roi_indices=roi)
    count=parameter_count(model) if model else (payload['coef'].size if kind=='residual_ridge' else 0)
    summary=dict(kind=kind,seed=seed,train_size=len(y),validation_size=len(vy),
                 validation=risk_summary(risks),train=risk_summary(tr),anchor_validation=risk_summary(anchor_risks),
                 parameter_count=count,cnn_comparator_parameters=target_count,
                 parameter_count_relative_to_cnn=float(count/target_count),selected_update=best_update,
                 selected_initialization=bool(model is not None and best_update==0),stop_reason=stop,
                 cpu_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall,device=device,threads=1,
                 batch_size=batch_size,epochs_limit=epochs,max_seconds=max_seconds,learning_rate=learning_rate,
                 weight_decay=weight_decay,ridge_alpha=payload.get('alpha'),weak_rank=geometry['weak'].shape[1],
                 strong_rank=geometry['strong'].shape[1],j_is_whitened=geometry['j_is_whitened'],
                 input_contract=('chi0 real/imag + broadcast u2; train-only normalization' if kind=='local_anchor_cnn' else
                     'chi0 real/imag + broadcast u2 + normalized weak/strong sensitivity and complex residual-backprojection images'
                     if kind=='local_geometry_cnn' else 'see model kind and checkpoint feature normalization'),
                 architecture=('3 shared 3x3 convolution blocks, 16 channels; zero-init 1x1 two-component material residual head'
                     if kind in ('local_anchor_cnn','local_geometry_cnn') else kind),
                 estimator_type=payload['estimator_type'],calibrated_density='NOT_CLAIMED',history=history)
    write_json(out/'RESULT.json',summary)
    return summary
