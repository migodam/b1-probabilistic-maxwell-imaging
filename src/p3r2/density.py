"""Proper conditional density in ten fixed prior coordinates.

The learned variable is (z[:8], Phi^-1(u2)). Full66 proposals multiply this
density by the independent Gaussian density in the complementary z56, with
the exact probit Jacobian. These are proposal densities, not validated
posterior references. No local rotating coordinate charts are used.
"""
import copy
import time
from pathlib import Path
import numpy as np
import torch
from torch import nn
from scipy.special import ndtr, ndtri
from .risk_models import observation_features
from .risk_training import write_json

KINDS=('unit_prior','diag','full','mdn','flow')
LOG2PI=float(np.log(2*np.pi))


def load_density_data(path):
    if 'locked' in str(path).lower(): raise ValueError('locked truth forbidden in density training')
    with np.load(path,allow_pickle=False) as z:
        obs=z['observation']; q=z['q']; parents=z['parent'] if 'parent' in z else z['parents']
    if q.shape!=(len(obs),66) or not np.isfinite(q).all() or np.any(q[:,64:]<=0) or np.any(q[:,64:]>=1):
        raise ValueError('finite q66 and interior u2 required for probit labels')
    parents=parents.astype(str)
    if len(set(parents))!=len(obs): raise ValueError('unique independent parents required')
    observation_features(obs)  # Shape and complex packing validation.
    return obs,q,parents


def coordinate_frame(basis=None):
    if basis is None: return np.eye(64,dtype=np.float64)
    basis=np.asarray(basis,dtype=np.float64)
    if basis.shape!=(64,8) or np.max(np.abs(basis.T@basis-np.eye(8)))>1e-6:
        raise ValueError('optional fixed train-only frame must have orthonormal columns [64,8]')
    # Keep the supplied first eight vectors, complete the orthogonal complement.
    frame=np.linalg.qr(basis,mode='complete')[0]
    signs=np.sign(np.sum(frame[:,:8]*basis,axis=0)); frame[:,:8]*=signs
    return frame


def labels_from_q(q,frame):
    q=np.asarray(q)
    return np.column_stack((q[:,:64]@frame[:,:8],ndtri(q[:,64:]))).astype('float32')


def fit_context(observation,components=32):
    raw=observation_features(observation)
    mean=raw.mean(0,dtype=np.float64).astype('float32')
    scale=raw.std(0,dtype=np.float64).astype('float32'); scale[scale<1e-12]=1
    x=((raw-mean)/scale).astype('float64')
    covariance=x.T@x/len(x)
    eig,vec=np.linalg.eigh(covariance)
    rank=min(components,len(x)-1,x.shape[1])
    vectors=vec[:,-rank:][:,::-1]
    vectors*=np.where(vectors[np.argmax(np.abs(vectors),axis=0),np.arange(rank)]<0,-1.,1.)
    pca_scale=np.sqrt(np.maximum(eig[-rank:][::-1],1e-8))
    return dict(mean=mean,scale=scale,vectors=vectors.astype('float32'),
                pca_scale=pca_scale.astype('float32'),components=rank,
                fit='train observations only; per-feature scale followed by PCA32 and score scaling')


def context_features(observation,stats):
    x=(observation_features(observation)-stats['mean'])/stats['scale']
    return ((x@stats['vectors'])/stats['pca_scale']).astype('float32')


class Coupling(nn.Module):
    def __init__(self,context_dim,flip):
        super().__init__(); self.flip=flip
        self.net=nn.Sequential(nn.Linear(context_dim+5,96),nn.GELU(),nn.Linear(96,96),nn.GELU(),nn.Linear(96,10))
        nn.init.zeros_(self.net[-1].weight); nn.init.zeros_(self.net[-1].bias)

    def transform(self,y,context,inverse=False):
        if self.flip: y=torch.cat((y[:,5:],y[:,:5]),dim=1)
        a,b=y[:,:5],y[:,5:]
        s,t=self.net(torch.cat((context,a),dim=1)).chunk(2,dim=1)
        s=.5*torch.tanh(s)
        b=(b-t)*torch.exp(-s) if inverse else b*torch.exp(s)+t
        result=torch.cat((a,b),dim=1)
        if self.flip: result=torch.cat((result[:,5:],result[:,:5]),dim=1)
        return result,(-s.sum(1) if inverse else s.sum(1))


class ConditionalDensity(nn.Module):
    dimension=10
    def __init__(self,kind,context_dim,mixtures=4):
        super().__init__(); self.kind=kind; self.mixtures=mixtures
        if kind=='unit_prior': return
        self.context=nn.Sequential(nn.Linear(context_dim,96),nn.GELU(),nn.Linear(96,96),nn.GELU())
        if kind=='flow':
            self.layers=nn.ModuleList([Coupling(96,False),Coupling(96,True)])
        else:
            out={'diag':20,'full':65,'mdn':mixtures*21}[kind]
            self.head=nn.Linear(96,out)
            nn.init.zeros_(self.head.weight); nn.init.zeros_(self.head.bias)
            if kind=='mdn':
                # Break mixture symmetry while retaining a nearly unit-prior start.
                with torch.no_grad():
                    for k in range(mixtures):
                        self.head.bias[mixtures+k*20]=(k-(mixtures-1)/2)*.05
        if kind=='full':
            i,j=torch.tril_indices(10,10)
            self.register_buffer('tril_i',i); self.register_buffer('tril_j',j)

    def parameters_given(self,context):
        raw=self.head(self.context(context))
        if self.kind=='diag':
            mu,logsd=raw.chunk(2,dim=1)
            return mu,torch.clamp(logsd,-3,2)
        if self.kind=='mdn':
            logits=raw[:,:self.mixtures]
            components=raw[:,self.mixtures:].reshape(-1,self.mixtures,20)
            return logits,components[:,:,:10],torch.clamp(components[:,:,10:],-3,2)
        mu=raw[:,:10]; value=raw[:,10:]
        L=raw.new_zeros((len(raw),10,10))
        L[:,self.tril_i,self.tril_j]=.5*torch.tanh(value)
        diagonal=self.tril_i==self.tril_j
        L[:,self.tril_i[diagonal],self.tril_j[diagonal]]=torch.exp(torch.clamp(value[:,diagonal],-3,2))
        return mu,L

    def log_prob(self,context,y):
        if self.kind=='unit_prior': return -.5*(y*y+LOG2PI).sum(1)
        if self.kind=='flow':
            h=self.context(context); logdet=y.new_zeros(len(y)); z=y
            for layer in reversed(self.layers):
                z,ld=layer.transform(z,h,True); logdet+=ld
            return -.5*(z*z+LOG2PI).sum(1)+logdet
        if self.kind=='diag':
            mu,logsd=self.parameters_given(context)
            return (-.5*((y-mu)*torch.exp(-logsd))**2-logsd-.5*LOG2PI).sum(1)
        if self.kind=='mdn':
            logits,mu,logsd=self.parameters_given(context)
            components=(-.5*((y[:,None,:]-mu)*torch.exp(-logsd))**2-logsd-.5*LOG2PI).sum(2)
            return torch.logsumexp(torch.log_softmax(logits,dim=1)+components,dim=1)
        mu,L=self.parameters_given(context)
        standardized=torch.linalg.solve_triangular(L,(y-mu).unsqueeze(-1),upper=False).squeeze(-1)
        return -.5*(standardized**2+LOG2PI).sum(1)-torch.log(torch.diagonal(L,dim1=1,dim2=2)).sum(1)

    def sample(self,context,generator):
        noise=torch.randn((len(context),10),device=context.device,dtype=context.dtype,generator=generator)
        if self.kind=='unit_prior': return noise
        if self.kind=='flow':
            h=self.context(context); y=noise
            for layer in self.layers: y,_=layer.transform(y,h,False)
            return y
        if self.kind=='diag':
            mu,logsd=self.parameters_given(context); return mu+torch.exp(logsd)*noise
        if self.kind=='mdn':
            logits,mu,logsd=self.parameters_given(context)
            k=torch.multinomial(torch.softmax(logits,dim=1),1,generator=generator).squeeze(1)
            row=torch.arange(len(context),device=context.device)
            return mu[row,k]+torch.exp(logsd[row,k])*noise
        mu,L=self.parameters_given(context)
        return mu+(L@noise.unsqueeze(-1)).squeeze(-1)


class FixedCoordinateDensity:
    def __init__(self,payload,device='cpu'):
        self.p=payload; self.device=torch.device(device)
        self.frame=np.asarray(payload['frame'],dtype=np.float64)
        self.model=ConditionalDensity(payload['kind'],payload['context']['components']).to(self.device)
        self.model.load_state_dict(payload['state_dict']); self.model.eval()

    def log_prob_labels(self,observation,labels,batch_size=512):
        context=context_features(observation,self.p['context']); labels=np.asarray(labels,dtype='float32')
        if labels.shape!=(len(context),10): raise ValueError('aligned labels [N,10] required')
        result=[]
        with torch.inference_mode():
            for start in range(0,len(labels),batch_size):
                result.append(self.model.log_prob(torch.from_numpy(context[start:start+batch_size]).to(self.device),
                              torch.from_numpy(labels[start:start+batch_size]).to(self.device)).cpu().numpy())
        return np.concatenate(result).astype('float64')

    def log_prob(self,observation,q,batch_size=512):
        q=np.asarray(q,dtype=np.float64)
        if q.shape!=(len(observation),66): raise ValueError('aligned q[N,66] required')
        legal=np.isfinite(q).all(1)&(q[:,64:]>0).all(1)&(q[:,64:]<1).all(1)
        result=np.full(len(q),-np.inf)
        if np.any(legal):
            valid=q[legal]; labels=labels_from_q(valid,self.frame)
            lp=self.log_prob_labels(np.asarray(observation)[legal],labels,batch_size)
            rest=valid[:,:64]@self.frame[:,8:]
            eta=ndtri(valid[:,64:])
            lp+=-.5*np.sum(rest*rest+LOG2PI,axis=1)
            lp+=.5*np.sum(eta*eta+LOG2PI,axis=1)  # log |d eta / d u|.
            result[legal]=lp
        return result

    log_prob_q=log_prob

    def sample_labels(self,observation,n_samples=1,seed=0,batch_size=512):
        context=context_features(observation,self.p['context'])
        repeated=np.repeat(context,n_samples,axis=0)
        gen=torch.Generator(device=self.device); gen.manual_seed(seed)
        values=[]
        with torch.inference_mode():
            for start in range(0,len(repeated),batch_size):
                c=torch.from_numpy(repeated[start:start+batch_size]).to(self.device)
                values.append(self.model.sample(c,gen).cpu().numpy())
        return np.concatenate(values).reshape(len(context),n_samples,10).astype('float64')

    def sample(self,observation,n_samples=1,seed=0,batch_size=512):
        labels=self.sample_labels(observation,n_samples,seed,batch_size)
        rng=np.random.default_rng(seed+104729)
        z=labels[:,:,:8]@self.frame[:,:8].T
        z+=rng.normal(size=(len(observation),n_samples,56))@self.frame[:,8:].T
        u=ndtr(labels[:,:,8:])
        if np.any((u<=0)|(u>=1)):
            raise FloatingPointError('probit sample reached floating-point boundary; no silent clipping or atom introduced')
        return np.concatenate((z,u),axis=2)


def load_density(path,device='cpu'):
    return FixedCoordinateDensity(torch.load(path,map_location='cpu',weights_only=False),device)


def train_density(kind,seed,train_obs,train_q,train_parents,val_obs,val_q,val_parents,
                  output,context,frame,epochs=100,batch_size=256,patience=12,max_seconds=300,
                  learning_rate=2e-4,weight_decay=1e-4,device='cpu'):
    torch.set_num_threads(1)
    if torch.get_num_interop_threads()!=1: torch.set_num_interop_threads(1)
    torch.manual_seed(seed); rng=np.random.default_rng(seed)
    out=Path(output); out.mkdir(parents=True,exist_ok=True); (out/'epochs.jsonl').write_text('')
    wall,cpu=time.monotonic(),time.process_time(); dev=torch.device(device)
    x,vx=context_features(train_obs,context),context_features(val_obs,context)
    y,vy=labels_from_q(train_q,frame),labels_from_q(val_q,frame)
    model=ConditionalDensity(kind,context['components']).to(dev)
    def rows(features,labels):
        model.eval(); result=[]
        with torch.inference_mode():
            for begin in range(0,len(labels),batch_size):
                c=torch.from_numpy(features[begin:begin+batch_size]).to(dev)
                l=torch.from_numpy(labels[begin:begin+batch_size]).to(dev)
                result.append((-model.log_prob(c,l)).cpu().numpy())
        return np.concatenate(result).astype('float64')
    initial=rows(vx,vy); best=float(initial.mean()); state=copy.deepcopy(model.state_dict())
    history=[dict(epoch=0,update=0,validation_nll=best,train_nll=float(rows(x,y).mean()),wall_seconds=time.monotonic()-wall)]
    with (out/'epochs.jsonl').open('a') as f: f.write(__import__('json').dumps(history[-1])+'\n')
    update=0; best_update=0; stale=0; stop='UNIT_PRIOR_BASELINE'
    if kind!='unit_prior':
        optimizer=torch.optim.AdamW(model.parameters(),lr=learning_rate,weight_decay=weight_decay)
        scheduler=torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer,factor=.5,patience=4,min_lr=1e-6)
        stop='EPOCH_LIMIT'
        for epoch in range(1,epochs+1):
            if time.monotonic()-wall>=max_seconds: stop='TIME_LIMIT'; break
            model.train(); order=rng.permutation(len(y)); total=0.; seen=0; interrupted=False
            for begin in range(0,len(y),batch_size):
                if time.monotonic()-wall>=max_seconds: interrupted=True; break
                ii=order[begin:begin+batch_size]; c=torch.from_numpy(x[ii]).to(dev); target=torch.from_numpy(y[ii]).to(dev)
                optimizer.zero_grad(set_to_none=True); loss=-model.log_prob(c,target).mean()
                if not torch.isfinite(loss): raise RuntimeError('nonfinite density NLL')
                loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(),5.); optimizer.step()
                total+=float(loss.detach().cpu())*len(ii); seen+=len(ii); update+=1
            val=float(rows(vx,vy).mean()); scheduler.step(val)
            row=dict(epoch=epoch,update=update,completed_epoch=not interrupted,validation_nll=val,
                     train_minibatch_nll=total/max(seen,1),learning_rate=optimizer.param_groups[0]['lr'],
                     wall_seconds=time.monotonic()-wall)
            history.append(row)
            with (out/'epochs.jsonl').open('a') as f: f.write(__import__('json').dumps(row)+'\n')
            if val<best: best=val; state=copy.deepcopy(model.state_dict()); best_update=update; stale=0
            else: stale+=1
            if interrupted: stop='TIME_LIMIT'; break
            if stale>=patience: stop='VALIDATION_PATIENCE'; break
    model.load_state_dict(state)
    val_rows=rows(vx,vy); train_rows=rows(x,y)
    prior_rows=.5*np.sum(vy.astype('float64')**2+LOG2PI,axis=1)
    payload=dict(format='p3r2_fixed10_density_v1',kind=kind,seed=seed,context=context,frame=frame,
                 state_dict={k:v.detach().cpu() for k,v in state.items()},best_update=best_update,
                 train_parent=train_parents,label_definition='first8 fixed orthonormal z64 coordinates + probit(u2)',
                 label_normalization='fixed standard-normal prior reference: mean0 scale1, no fitted label rescaling',
                 reference_measure='dz8 d eta2; full66 is dz64 du2 with probit Jacobian and prior z56',
                 scientific_status='learned proposal density; posterior fidelity NOT_ESTABLISHED')
    torch.save(payload,out/'BEST.pt')
    np.savez_compressed(out/'validation_nll.npz',parent=val_parents,nll=val_rows,unit_prior_nll=prior_rows)
    result=dict(kind=kind,seed=seed,train_size=len(y),validation_size=len(vy),selected_update=best_update,
                selected_initialization=bool(kind!='unit_prior' and best_update==0),stop_reason=stop,
                validation_nll=float(val_rows.mean()),validation_nll_sem=float(val_rows.std(ddof=1)/np.sqrt(len(val_rows))),
                unit_prior_validation_nll=float(prior_rows.mean()),train_nll=float(train_rows.mean()),
                parameter_count=sum(p.numel() for p in model.parameters()),context_components=context['components'],
                dimension=10,full66_proposal='learned10 times independent prior56 with exact probit Jacobian',
                posterior_fidelity='NOT_ESTABLISHED',device=device,threads=1,cpu_seconds=time.process_time()-cpu,
                wall_seconds=time.monotonic()-wall,epochs_limit=epochs,batch_size=batch_size,max_seconds=max_seconds,
                learning_rate=learning_rate,weight_decay=weight_decay,history=history)
    write_json(out/'RESULT.json',result)
    return result
