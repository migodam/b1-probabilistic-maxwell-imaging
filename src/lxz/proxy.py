"""Training-free inference for the EXISTING Joint P2 proxy target.

This module never imports a Maxwell forward solver.  Geometry is frozen before
sampling. Rejected states remain in the chain; correctors are not sampling maps.
"""
from pathlib import Path
import argparse, json, time, sys
import numpy as np
from scipy.special import ndtr, ndtri, logsumexp, expit
from b1_joint.material import BOUNDS, WIDTH, geometry, event_ids, prior_sample, prior_logpdf
from b1_joint.response import Response
from b1_joint.likelihood import log_likelihood
from b1_joint.proposal_correction import measurement_condition, load_conditioned_proposal

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'research/research_lxz_20261016'
PAID=ROOT/'research/joint_p2/results'
IDENTITY='joint-p2-radial-v1:C:tensor-Chebyshev-five-per-physical-axis-radial-v1'

def plain(x):
    if isinstance(x,dict):return {k:plain(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [plain(v) for v in x]
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    return x

def dump(path,x):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(plain(x),ensure_ascii=False,indent=2,allow_nan=False)+'\n')

class MaterialPrior:
    """Original coefficient-uniform measure; Gaussian coordinates are bijective.

    Decoder non-injectivity is retained, including inactive coefficient volume.
    """
    @staticmethod
    def from_u(u):
        x=BOUNDS[:,0]+WIDTH*ndtr(u)
        if np.any((x<=BOUNDS[:,0])|(x>=BOUNDS[:,1])):
            raise FloatingPointError('rounded prior coordinate; no clipping')
        return x
    @staticmethod
    def to_u(x):return ndtri((np.asarray(x)-BOUNDS[:,0])/WIDTH)

class MaxwellTarget:
    """Explicitly proxy-only; counts each attempted likelihood, including setup."""
    identity=IDENTITY
    def __init__(self,response,y,sigma,limit=4096):
        self.response=response;self.y=np.asarray(y);self.sigma=np.asarray(sigma)
        self.calls=0;self.failures=0;self.limit=limit
    def __call__(self,x):
        a=np.atleast_2d(x);n=len(a)
        if self.calls+n>self.limit:raise RuntimeError('proxy likelihood budget exhausted')
        self.calls+=n
        try:
            f=self.response.predict(a,(1.6,))
            ll=log_likelihood(f,self.y,self.sigma,.15)
        except Exception:
            self.failures+=n;raise
        return ll
    def frozen_geometry(self,x,nuisance_projected=False):
        # 8 paid proxy likelihoods also produce the actual response derivatives.
        u=MaterialPrior.to_u(x);h=1e-3
        rows=np.array([u+s*h*np.eye(4)[j] for j in range(4) for s in [-1,1]])
        xx=MaterialPrior.from_u(rows)
        if self.calls+8>self.limit:raise RuntimeError('geometry setup budget')
        self.calls+=8
        f=self.response.predict(xx,(1.6,))
        # All response evaluations count even though only derivatives are needed.
        log_likelihood(f,self.y,self.sigma,.15)
        A=np.column_stack([(f[2*j+1]-f[2*j])/(2*h)/self.sigma for j in range(4)])
        H=A.T@A
        if nuisance_projected:
            # Local Gaussian gain geometry is a proposal metric ONLY. The target
            # always uses the exact bounded-gain marginal with its determinant.
            a=np.mean(f,axis=0)/self.sigma
            H-=np.outer(A.T@a,A.T@a)/(a@a+1/.15**2)
        vals,V=np.linalg.eigh((H+H.T)/2)
        rates=np.maximum(.025,.45/np.sqrt(1+np.maximum(vals,0)))
        return (V*np.sqrt(1-rates*rates))@V.T,(V*rates)@V.T,vals

class FrozenDensity:
    """q_x=.9 normalized C1 MDN + .1 original prior, evaluable in dx."""
    def __init__(self,path,y,sigma):
        start=time.perf_counter()
        self.p,self.receipt,self.api=load_conditioned_proposal(path,measurement_condition(y,sigma))
        self.cold_seconds=time.perf_counter()-start
        self.mu=self.p['means'][0];self.chol=self.p['cholesky'][0]
        self.weights=self.p['weights'][0]
        self.inverse=np.linalg.inv(self.chol)
        self.constant=np.log(self.weights)-np.log(np.diagonal(self.chol,axis1=-2,axis2=-1)).sum(1)-2*np.log(2*np.pi)
    def draw(self,rng,n):
        prior=rng.random(n)<.1
        ids=rng.choice(len(self.weights),n,p=self.weights)
        z=self.mu[ids]+np.einsum('nij,nj->ni',self.chol[ids],rng.normal(size=(n,4)))
        x=BOUNDS[:,0]+WIDTH*expit(z)
        x[prior]=prior_sample(rng,int(prior.sum()))
        if np.any((x<=BOUNDS[:,0])|(x>=BOUNDS[:,1])):raise FloatingPointError('rounded MDN endpoint')
        return x,self.logpdf(x)
    def logpdf(self,x):
        x=np.atleast_2d(x);low=x-BOUNDS[:,0];high=BOUNDS[:,1]-x
        z=np.log(low)-np.log(high)
        r=np.einsum('kij,nkj->nki',self.inverse,z[:,None,:]-self.mu)
        mdn=logsumexp(self.constant-.5*np.sum(r*r,axis=-1),axis=1)+(np.log(WIDTH)-np.log(low)-np.log(high)).sum(1)
        return np.logaddexp(np.log(.9)+mdn,np.log(.1)+prior_logpdf(x))

def ess_1d(x):
    x=np.asarray(x,float);n=len(x);x=x-x.mean();v=x@x
    if v<=1e-28:return None
    size=1<<(2*n-1).bit_length()
    f=np.fft.rfft(x,size);ac=np.fft.irfft(f*np.conj(f),size)[:n]/v
    pairs=ac[1:1+2*((n-1)//2)].reshape(-1,2).sum(1)
    stop=np.where(pairs<=0)[0];m=int(stop[0]) if len(stop) else len(pairs)
    return float(min(n,n/max(1.,1+2*pairs[:m].sum())))

def metrics(x,w,ref,event_ref,cdf_ref):
    g=np.column_stack((geometry(x),x[:,2:]));e=event_ids(x)
    p=np.bincount(e,weights=w,minlength=12)
    errs=[]
    for j in range(4):
        grid=np.array(cdf_ref['CDF_grid'][j]);expected=np.array(cdf_ref['CDF'][j])
        estimate=np.array([w[g[:,j]<=t].sum() for t in grid])
        errs.append(float(abs(estimate-expected).max()))
    mean=w@g;sd=np.sqrt(w@((g-mean)**2));important=event_ref>=.05
    return dict(event_mass=p,event_MAE=float(abs(p-event_ref).mean()),event_max_error=float(abs(p-event_ref).max()),
        significant_event_recall=float(np.mean(p[important]>=.005)),false_exclusions=int(np.sum(important&(p<.005))),
        false_confidence=float(p.max()-event_ref.max()),CDF_max_error=errs,mean=mean,SD=sd,
        mean_absolute_fraction=float(np.max(abs(mean-np.array(ref['mean']))/np.array([.29256,.11702,.6,.12]))),
        unique_material_geometry=int(len(np.unique(np.round(g,12),axis=0))),
        probability_status='finite proxy estimate; no guarantee of all-mode discovery')

def local_global(method,target,initial,initial_ll,q,seed,round_id):
    rng=np.random.default_rng(seed);start=time.perf_counter()
    x=initial[np.argmax(initial_ll)].copy();ll=float(max(initial_ll));u=MaterialPrior.to_u(x)
    R=np.eye(4)*np.sqrt(1-.35**2);B=np.eye(4)*.35;spectrum=None
    if method=='geometry':R,B,spectrum=target.frozen_geometry(x,round_id==2)
    setup=target.calls
    n=4096-target.calls
    global_x,global_logq=q.draw(np.random.default_rng(seed+991),n)
    # Proposal acquisition is separately timed but all unused proposals are
    # charged as setup, avoiding a free pre-generated bank.
    cold_setup=time.perf_counter()-start
    logq=float(q.logpdf(x)[0]);samples=[];lls=[];accept=[];branches=[];calls=[];offsets=[];snapshots=[]
    for i in range(n):
        glob=method in ('global','prior') or (method=='hybrid' and rng.random()<.25)
        if glob:
            if method=='prior':v=prior_sample(rng,1)[0];offset=0.
            else:v=global_x[i];offset=logq-float(global_logq[i])
        else:
            v=MaterialPrior.from_u(R@u+B@rng.normal(size=4));offset=0.
        lp=float(target(v)[0]);accepted=np.log(rng.random())<min(0.,lp-ll+offset)
        if accepted:
            x=v;u=MaterialPrior.to_u(x);ll=lp;logq=float(q.logpdf(x)[0])
        samples.append(x.copy());lls.append(ll);accept.append(accepted);branches.append(glob);calls.append(target.calls);offsets.append(offset)
        if target.calls in (1024,4096):snapshots.append((target.calls,len(samples)))
    arr=np.array(samples)
    meta=dict(method=method,seed=seed,round=round_id,proxy_calls=target.calls,setup_likelihood_calls=setup,
        setup_seconds=cold_setup,warm_search_seconds=time.perf_counter()-start-cold_setup,
        total_seconds=time.perf_counter()-start,metric_eigenvalues=spectrum,acquisition_scope='proposal pre-generation + exact q; no teacher fitting',
        accepted=int(sum(accept)),global_trials=int(sum(branches)),global_accepted=int(sum(a and b for a,b in zip(accept,branches))))
    arrays=dict(latent=arr,log_likelihood=np.array(lls),accepted=np.array(accept),global_move=np.array(branches),
        proposal_log_offset=np.array(offsets),likelihood_calls=np.array(calls))
    return meta,arrays,snapshots

def tempered(target,initial,initial_ll,seed,round_id):
    rng=np.random.default_rng(seed);start=time.perf_counter();x=initial.copy();ll=initial_ll.copy();n=len(x)
    logw=np.full(n,-np.log(n));beta=0.;trace=[];snapshots=[];store={};accepted=0
    # 64 initial + 63 complete 64-particle sweeps = 4096 evaluations.
    for sweep in range(63):
        bn=[.10,.20,.35,.55,.80,1.][min(sweep,5)]
        logw+=(bn-beta)*ll;logw-=logsumexp(logw);w=np.exp(logw)
        ess=1/(w@w)
        if ess<.8*n:
            ids=np.searchsorted(np.cumsum(w),(np.arange(n)+rng.random())/n)
            x=x[ids].copy();ll=ll[ids].copy();logw.fill(-np.log(n))
        u=MaterialPrior.to_u(x);proposal=MaterialPrior.from_u(np.sqrt(1-.35**2)*u+.35*rng.normal(size=u.shape))
        lp=target(proposal);ok=np.log(rng.random(n))<np.minimum(0,bn*(lp-ll))
        x[ok]=proposal[ok];ll[ok]=lp[ok];accepted+=int(ok.sum());beta=bn
        trace.append(dict(sweep=sweep,beta=beta,pre_resampling_ESS=ess,accepted=int(ok.sum()),calls=target.calls))
        if target.calls in (1024,4096):
            store[str(target.calls)+'_latent']=x.copy();store[str(target.calls)+'_weights']=np.exp(logw);snapshots.append((target.calls,0))
    meta=dict(method='tempered',seed=seed,round=round_id,proxy_calls=target.calls,setup_likelihood_calls=64,
        total_seconds=time.perf_counter()-start,setup_seconds=0.,warm_search_seconds=time.perf_counter()-start,
        accepted=accepted,global_trials=0,global_accepted=0,beta=beta,trace=trace)
    return meta,store,snapshots

class LocalGlobalSampler:
    """The inference boundary: fixed proposals, fixed material, metered target.

    Discovery correctors are intentionally absent. The output retains rejected
    states; archive images or event labels never alter the transition density.
    """
    def __init__(self,method,proposal,seed,round_id=1):
        self.method,self.proposal,self.seed,self.round_id=method,proposal,seed,round_id
    def run(self,target,initial,initial_ll):
        if self.method=='tempered':return tempered(target,initial,initial_ll,self.seed,self.round_id)
        return local_global(self.method,target,initial,initial_ll,self.proposal,self.seed,self.round_id)

def controls():
    """Deterministic finite-state probability mechanism checks, no physics."""
    p=np.array([.1,.15,.2,.25,.3]);n=len(p)
    # Symmetric base proposals are uniform-prior reversible. Selector ratio is
    # needed even though each individual MH kernel preserves the same target.
    Q=[np.full((n,n),1/n),.7*np.eye(n)+.3*np.full((n,n),1/n)]
    selectors=np.column_stack((np.linspace(.1,.9,n),np.linspace(.9,.1,n)))
    def kernel(correct):
        K=np.zeros((n,n))
        for i in range(n):
            for a,q in enumerate(Q):
                for j in range(n):
                    ratio=p[j]/p[i]*(selectors[j,a]/selectors[i,a] if correct else 1)
                    K[i,j]+=selectors[i,a]*q[i,j]*min(1,ratio)
        K[np.diag_indices(n)]+=1-K.sum(1);return K
    good,bad=kernel(True),kernel(False)
    q=np.array([.32,.21,.12,.18,.17]);K=np.zeros((n,n));wrong=K.copy()
    for i in range(n):
        for j in range(n):
            K[i,j]=q[j]*min(1,p[j]*q[i]/(p[i]*q[j]));wrong[i,j]=q[j]*min(1,p[j]/p[i])
    for a in [K,wrong]:a[np.diag_indices(n)]+=1-a.sum(1)
    result=dict(selector_correct_invariance=float(abs(p@good-p).max()),selector_uncorrected_bias=float(abs(p@bad-p).max()),
        global_MH_correct_invariance=float(abs(p@K-p).max()),global_omitted_q_bias=float(abs(p@wrong-p).max()),
        high_ESS_missed_half_mass=dict(ESS=10000,represented_true_mass=.5,event_error=.5),
        prior_coordinate_CDF_bijection='x=lower+width*Phi(u); prior density cancels only for prior-reversible local kernel',
        corrector_in_sampling=False,new_physics=0)
    assert result['selector_correct_invariance']<1e-14 and result['global_MH_correct_invariance']<1e-14
    assert result['selector_uncorrected_bias']>.001 and result['global_omitted_q_bias']>.001
    dump(OUT/'results/PROBABILITY_CONTROLS.json',result)

def run(round_id=1,repair=False):
    protocol=json.loads((OUT/'configs/PROTOCOL.json').read_text());cases=protocol['E1']['observations'];seeds=protocol['E1']['seeds']
    dest=OUT/'results'/f'round{round_id}{"_geometry_repair" if repair else ""}';dest.mkdir(parents=True,exist_ok=True)
    if (dest/'COMPLETE.json').exists():raise RuntimeError('completed campaign cannot silently repeat')
    response=Response(PAID/'COMMON_RESPONSE_BANK.npz');data=json.loads((PAID/'DEVELOPMENT_CASES.json').read_text())['cases']
    reference=json.loads((PAID/'ADAPTIVE_PHYSICAL_REFERENCE.json').read_text())['cases']
    continuous=json.loads((PAID/'PILOT_REFERENCE.json').read_text())['cases'];records=[];total=0;t0=time.perf_counter()
    c1=Path(protocol['C1_root']);sys.path.insert(0,str(c1/'p2'))
    for name in cases:
        item=next(c for c in data if c['id']==name);y=np.array(item['y'])[:32];sigma=np.array(item['sigma'])[:32]
        loadstart=time.perf_counter();q=FrozenDensity(c1/'p2/results/models/full_distill-26101321.pt',y,sigma)
        dump(dest/f'{name}_PROPOSAL.json',dict(receipt=q.receipt,cold_model_seconds=time.perf_counter()-loadstart,parameters=q.p))
        for seed in seeds:
            initial=prior_sample(np.random.default_rng(seed),64)
            methods=protocol['E1']['methods'] if round_id==1 and not repair else ['geometry']
            # Rotate order; each arm repeats and pays its own 64-point setup.
            shift=seed%len(methods);methods=methods[shift:]+methods[:shift]
            for method in methods:
                target=MaxwellTarget(response,y,sigma);st=time.perf_counter();initial_ll=target(initial)
                meta,arrays,snaps=LocalGlobalSampler(method,q,seed,round_id).run(target,initial,initial_ll)
                meta['inclusive_seconds']=time.perf_counter()-st
                total+=target.calls;np.savez_compressed(dest/f'{name}_{method}_{seed}.npz',**arrays)
                for calls,count in snaps:
                    if method=='tempered':x=arrays[str(calls)+'_latent'];w=arrays[str(calls)+'_weights'];burn=0
                    else:
                        burn=count//4;x=arrays['latent'][burn:count];w=np.full(len(x),1/len(x))
                    r=metrics(x,w,reference[name],np.array(reference[name]['event_mass']),continuous[name]['C']['continuous'][0])
                    e=event_ids(x)
                    r.update(meta);r.update(observation_id=name,checkpoint=calls,burnin_states=burn,
                        ESS_importance=float(1/(w@w)) if method=='tempered' else None,
                        ESS_material=[ess_1d(z) for z in np.column_stack((geometry(x),x[:,2:])).T] if method!='tempered' else None,
                        ESS_event=[ess_1d(e==j) for j in range(12)] if method!='tempered' else None)
                    if method!='tempered':
                        g=arrays['global_move'][:count];a=arrays['accepted'][:count]
                        r['acceptance_checkpoint']=float(a.mean());r['global_acceptance_checkpoint']=float(a[g].mean()) if g.any() else None
                    records.append(r)
                dump(dest/'METRICS.json',dict(rows=records,target=IDENTITY,scope='PROXY_TARGET_MECHANISM_ONLY'))
                with (OUT/'COST_LEDGER.jsonl').open('a') as f:
                    f.write(json.dumps(plain(dict(event='E1_CHAIN_COMPLETE',round=round_id,repair=repair,observation=name,method=method,seed=seed,
                        proxy_likelihood_attempts=target.calls,failures=target.failures,new_full_states=0,new_RHS=0,
                        wall_seconds=meta['inclusive_seconds'])))+'\n')
                print(json.dumps(dict(round=round_id,observation=name,method=method,seed=seed,calls=target.calls,event_MAE=records[-1]['event_MAE'])),flush=True)
    dump(dest/'COMPLETE.json',dict(proxy_calls=total,wall_seconds=time.perf_counter()-t0,full_states=0,RHS=0))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('task',choices=['controls','run']);p.add_argument('--round',type=int,default=1);p.add_argument('--repair',action='store_true');args=p.parse_args()
    controls() if args.task=='controls' else run(args.round,args.repair)
