"""Small proposal learning from paid discovery trajectories, not posterior distillation."""
import json,time,copy
import numpy as np
from .core import OUT,write,save,BlockStop
from .explore import geometry,direction,constrained_fit
from p1_material.proposals import frozen_latent_step

def rows(stage):
    p=OUT/f'{stage}/candidates.jsonl'
    return [json.loads(x) for x in p.read_text().splitlines()] if p.exists() else []

def feature(t,p,y,s,g):
    rng=np.random.default_rng(2026101960);proj=np.linalg.qr(rng.normal(size=(512,32)))[0]
    r=t.white(p.state.data-y);eig=np.log1p(np.quantile(g['lam']*s*s,np.linspace(0,1,16)))
    return np.r_[p.q,proj.T@r/24,eig,s/.02].astype(np.float32)

def teachers(t,split,anchors):
    parents=json.loads((OUT/'configs/PARENTS.json').read_text());ids={r['parent'] for r in parents if r['split']==split};data=[];seen=set();context={}
    for row in rows('B')+rows('C'):
        parent=row['parent']
        if parent not in ids or not row['data_feasible']:continue
        if parent not in context:
            p=t.load(anchors[parent]);y=t.observation(parent);g=geometry(t,p,.02);context[parent]=(p,y,g)
        anchor,y,g=context[parent]
        path=next((OUT/s/'endpoints'/f'{row["key"]}.npz' for s in ['B','C'] if (OUT/s/'endpoints'/f'{row["key"]}.npz').exists()),None)
        if path is None:continue
        with np.load(path) as a:q=a['q'].copy();chi=a['chi'].copy()
        distance=t.prior.material_rms(chi,anchor.state.chi);dphi=row['Phi']-t.phi(anchor,y)
        if distance<.0196 or dphi>2:continue
        k=(parent,tuple(np.round(q,8)))
        if k in seen:continue
        seen.add(k);s=.02 if distance<.03 else .04
        data.append((feature(t,anchor,y,s,g),q[:64]-anchor.q[:64],parent))
    return data

def build_model(kind,nfeatures):
    import torch
    from torch import nn
    class Model(nn.Module):
        def __init__(self):
            super().__init__();self.kind=kind;self.encoder=nn.Sequential(nn.Linear(nfeatures,96),nn.Tanh(),nn.Linear(96,96),nn.Tanh())
            self.gaussian=nn.Linear(96,128)
            self.couplings=nn.ModuleList([nn.Sequential(nn.Linear(128,96),nn.Tanh(),nn.Linear(96,64)) for _ in range(4)]) if kind=='flow' else nn.ModuleList()
            self.register_buffer('xmean',torch.zeros(nfeatures));self.register_buffer('xscale',torch.ones(nfeatures));self.register_buffer('dmean',torch.zeros(64));self.register_buffer('dscale',torch.ones(64))
        def context(self,x):return self.encoder((x-self.xmean)/self.xscale)
        def parameters_for(self,c):
            a=self.gaussian(c);return a[...,:64],-1+3*torch.tanh(a[...,64:]/3)
        def reverse(self,y,c):
            ld=torch.zeros(y.shape[:-1],device=y.device)
            for i in range(len(self.couplings)-1,-1,-1):
                flip=i%2;z=y.flip(-1) if flip else y;a,b=z[...,:32],z[...,32:];st=self.couplings[i](torch.cat([a,c],-1));sc=.5*torch.tanh(st[...,:32]);b=(b-st[...,32:])*torch.exp(-sc);y=torch.cat([a,b],-1);y=y.flip(-1) if flip else y;ld-=sc.sum(-1)
            return y,ld
        def log_prob(self,d,x):
            c=self.context(x);y=(d-self.dmean)/self.dscale;ld=-torch.log(self.dscale).sum()
            if self.kind=='flow':y,l=self.reverse(y,c);ld=ld+l
            mu,ls=self.parameters_for(c);z=(y-mu)*torch.exp(-ls)
            return -.5*(z*z+np.log(2*np.pi)).sum(-1)-ls.sum(-1)+ld
        def sample(self,x):
            c=self.context(x);mu,ls=self.parameters_for(c);y=mu+torch.exp(ls)*torch.randn_like(mu)
            for i,net in enumerate(self.couplings):
                flip=i%2;z=y.flip(-1) if flip else y;a,b=z[...,:32],z[...,32:];st=net(torch.cat([a,c],-1));sc=.5*torch.tanh(st[...,:32]);y=torch.cat([a,b*torch.exp(sc)+st[...,32:]],-1);y=y.flip(-1) if flip else y
            return self.dmean+self.dscale*y
    return Model()

def train(t,anchors):
    import torch
    torch.set_num_threads(2);tr=teachers(t,'train',anchors);va=teachers(t,'validation',anchors)
    write(OUT/'D/TEACHER_PROVENANCE.json',dict(train_examples=len(tr),validation_examples=len(va),train_parents=sorted({x[2] for x in tr}),teacher='full-validated task-scale likelihood-competitive discovery endpoints; NOT posterior draws',heldout_examples=0,teacher_physics_cost_in_B_C=True,new_teacher_physics=0,learned_coordinates='64 Gaussian updates; uniform properties use shared reflected proposal',GT_supervision=False))
    if not tr:
        write(OUT/'D/TRAINING.json',dict(status='NOT_RUN_NO_NONZERO_TEACHER',GPU_seconds=0));return {},None
    X=torch.tensor(np.array([a[0] for a in tr]));Y=torch.tensor(np.array([a[1] for a in tr]),dtype=torch.float32)
    VX=torch.tensor(np.array([a[0] for a in va])) if va else X;VY=torch.tensor(np.array([a[1] for a in va]),dtype=torch.float32) if va else Y
    models={};log=[]
    for kind in ['gaussian','flow']:
        for seed in [2026101970,2026101971,2026101972]:
            start=time.perf_counter();torch.manual_seed(seed);model=build_model(kind,X.shape[1]);model.xmean.copy_(X.mean(0));model.xscale.copy_(X.std(0,unbiased=False).clamp_min(.05));model.dmean.copy_(Y.mean(0));model.dscale.copy_(Y.std(0,unbiased=False).clamp_min(.01))
            opt=torch.optim.Adam(model.parameters(),lr=1e-3);best=None;score=float('inf');history=[]
            for step in range(3000):
                if time.perf_counter()-start>600:break
                index=torch.randint(len(X),(128,));opt.zero_grad();loss=-model.log_prob(Y[index],X[index]).mean()
                if not torch.isfinite(loss):raise RuntimeError('nonfinite normalized proposal loss')
                loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),10.);opt.step()
                if step%100==0 or step==2999:
                    with torch.no_grad():v=float(-model.log_prob(VY,VX).mean())
                    history.append(dict(update=step,train_NLL=float(loss.detach()),validation_NLL=v))
                    if v<score:score=v;best=copy.deepcopy(model.state_dict())
            model.load_state_dict(best);model.eval();path=OUT/f'D/models/{kind}_{seed}.pt';path.parent.mkdir(parents=True,exist_ok=True);torch.save(dict(kind=kind,nfeatures=X.shape[1],state_dict=best,normalization='exact Gaussian / affine-coupling log density in delta-z',target_id=t.target_id),path)
            models[kind,seed]=model;row=dict(kind=kind,seed=seed,updates=step+1,wall_seconds=time.perf_counter()-start,parameters=sum(p.numel() for p in model.parameters()),validation_NLL=score,history=history,device='CPU',GPU_seconds=0,Maxwell_calls_in_training=0);log.append(row)
            t.book.append('NN_TRAINING',**row);print(json.dumps({k:v for k,v in row.items() if k!='history'}),flush=True)
    write(OUT/'D/TRAINING.json',dict(status='COMPLETE',runs=log,posterior_distillation=False,train_parent_limitations='12 train parents; only subset has qualified teachers; no broad generalization claim',GPU_seconds=0))
    return models,X.shape[1]

def learned_trial(t,parent,anchor,method,seed,s,model=None,trial=0):
    import torch
    rng=np.random.default_rng(np.random.SeedSequence([seed,trial,int(s*10000)]));torch.manual_seed(seed+trial)
    t.parent=parent;y=t.observation(parent);g=geometry(t,anchor,s);before=t.book.used['D'].copy();start=time.perf_counter()
    if model is None:
        v,w=direction(t,anchor,g,'soft_dynamic',rng)
        raw=anchor.q.copy()
        if v is not None:raw[:64]+=s*v[:64]
    else:
        with torch.no_grad():dz=model.sample(torch.tensor(feature(t,anchor,y,s,g))[None])[0].numpy().astype(float)
        raw=anchor.q.copy();raw[:64]+=dz
    # Identical independent property draw for every proposal family.
    prng=np.random.default_rng(np.random.SeedSequence([seed,trial,int(s*10000),991]))
    z=np.mod(raw[64:]+.04*prng.normal(size=2),2);raw[64:]=np.where(z<=1,z,2-z)
    from .explore import feasible_direction
    v,w=feasible_direction(t,anchor,raw-anchor.q)
    result=dict(parent=parent,method=method,seed=seed,trial=trial,radius=s)
    with t.book.scope(2,56):
        if v is None:result['status']='NO_MATERIAL_DIRECTION';return result
        q,matched=frozen_latent_step(t.prior,anchor.q,v,s)
        if matched:
            q,fit=constrained_fit(t,anchor,anchor,w,s,q,objective='data');matched=fit['radius_matched']
        if not matched:result['status']='DECODER_DISTANCE_FAILURE';return result
        p=t.evaluate(q,f'{parent}_{method}_{seed}_{trial}_predictor');result['predictor']=t.record(p,parent,method,anchor,seed=seed,trial=trial,role='learned_or_matched_control_predictor',target_RMS=s)
        try:
            t.jacobian(p);q2,fit=constrained_fit(t,p,anchor,w,s,p.q,objective='data')
            if fit['radius_matched'] and np.linalg.norm(q2-p.q)>1e-10:
                c=t.evaluate(q2,f'{parent}_{method}_{seed}_{trial}_corrector');result['corrector']=t.record(c,parent,method,anchor,seed=seed,trial=trial,role='common_FullJ_corrector',target_RMS=s)
        except BlockStop:result['block_limit']=True
    result.update(states=t.book.used['D']['STATE']-before['STATE'],RHS=t.book.used['D']['RHS']-before['RHS'],wall_seconds=time.perf_counter()-start)
    return result

def stage_D(t):
    parents=json.loads((OUT/'configs/PARENTS.json').read_text());anchors=json.loads((OUT/'A/ANCHORS.json').read_text());models,nf=train(t,anchors)
    if not models:
        write(OUT/'D/SUMMARY.json',dict(status='NOT_RUN_NO_TEACHERS',new_physics=0));return
    write(OUT/'D/MODELS_FROZEN_BEFORE_HELDOUT.json',dict(seeds=[2026101970,2026101971,2026101972],training_only='train',checkpoint='stored validation NLL; no physical teacher or GT in epochs',heldout_GT_opened=False))
    runs=[]
    for r in parents:
        if r['split']=='train':continue
        p=t.load(anchors[r['parent']]);y=t.observation(r['parent'])
        if not t.qualified(p,y):runs.append(dict(parent=r['parent'],status='INVALID_BASE'));continue
        seeds=[2026101970,2026101971] if r['split']=='validation' else [2026101970,2026101971,2026101972]
        radii=[.02,.04] if r['split']=='validation' else [.02,.04,.02,.04]
        for seed in seeds:
            for kind in ['nonlearned','gaussian','flow']:
                for j,s in enumerate(radii):runs.append(learned_trial(t,r['parent'],p,kind,seed,s,models.get((kind,seed)),j))
        print(json.dumps(dict(stage='D',parent=r['parent'],consumed=dict(t.book.used['D']))),flush=True)
    write(OUT/'D/SUMMARY.json',dict(runs=runs,normalized_density='raw delta-z proposal only; retraction/corrector output candidate-only',posterior_mass=None,GT_online=False))
