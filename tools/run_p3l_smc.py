"""Independent full-target tempering; every new likelihood is paid."""
from pathlib import Path
import sys,argparse,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,Target,write,save,probit_derivative
from p3l.reference import Posterior,FrozenMaterialMobility
from p3l.tempered import tempered_smc

def main():
    a=argparse.ArgumentParser();a.add_argument('parent');a.add_argument('--replicate',type=int,default=0);a.add_argument('--particles',type=int,default=128);a.add_argument('--mutation',choices=['pcn','mala'],default='mala');a.add_argument('--memory-resume',action='store_true');args=a.parse_args()
    suffix='_memory_resume1' if args.memory_resume else ''
    b=Book('A',f'smc_{args.parent}_{args.mutation}_r{args.replicate}'+suffix,max_wall=7200);status='COMPLETE';error=None
    try:
        t=Target(b);y=t.observation(args.parent);post=Posterior(t,y,probit=True)
        anchor=t.load_cache(OUT/f'A/anchors/{args.parent}.npz',owner='P3L paid deterministic anchor; proposal setup only')
        mobility=FrozenMaterialMobility(t,anchor) if args.mutation=='mala' else None
        rng=np.random.default_rng(np.random.SeedSequence([2026102281,int(args.parent.split('_')[-1]),args.replicate]))
        saved_particles=None;paid_ll={}
        if args.memory_resume:
            origin=OUT/f'A/jobs/smc_{args.parent}_{args.mutation}_r{args.replicate}'
            if json.loads((origin/'CLOSED.json').read_text())['status']!='INTERRUPTED_MEMORY_PRESSURE':raise RuntimeError('only documented interrupted SMC can resume')
            if len(list((origin/'particles').glob('STAGE_*.npz')))!=1:raise RuntimeError('SMC progressed beyond initialized prior; require separately reviewed later-stage resume')
            with np.load(origin/'particles/STAGE_000.npz') as z:saved_particles=z['x'].copy();oldll=z['loglikelihood'].copy()
            if saved_particles.shape!=(args.particles,66):raise RuntimeError('prior particle shape changed')
            rng.bit_generator.state=json.loads((origin/'CHECKPOINT.json').read_text())['RNG']
            paid_ll={x.tobytes():float(v) for x,v in zip(saved_particles,oldll)}
            b._update('PAID_PRIOR_PARTICLE_CACHE_REUSE',original_job=origin.name,paid_particles=args.particles,old_partial_mutation_costs_preserved=True,RNG_resume='initial prior checkpoint; unfinished mutation is retried and billed')
        def ll(x):
            k=np.asarray(x).tobytes()
            if k in paid_ll:b.counts['paid_initial_likelihood_reuses']+=1;return paid_ll[k]
            return -t.phi(post.point(x),y)
        def gl(x):return -t.phi_gradient(post.point(x),y)*probit_derivative(x)
        def precision(x,beta):return np.eye(66)+beta*(mobility(x)-np.eye(66))
        events=(b.jobdir/'EVENTS.jsonl').open('x')
        def callback(v):
            if v['event'] not in ['INITIALIZED','STAGE_COMPLETE','FINISHED','TEMPERING_FAILED']:return
            if 'particles' in v:
                stage=v['stage'];save(b.jobdir/f'particles/STAGE_{stage:03d}.npz',x=v['particles'],logweights=v['logweights'],loglikelihood=v['loglikelihood'],logprior=v['logprior'],ancestor_ids=v.get('ancestor_ids',v.get('history',{}).get('ancestor_ids',np.arange(args.particles))),beta=v['beta'])
                v={k:a for k,a in v.items() if k not in ['particles','logweights','loglikelihood','logprior','ancestor_ids']}
            from p3l.core import plain
            events.write(json.dumps(plain(v),allow_nan=False)+'\n');events.flush()
            write(b.jobdir/'CHECKPOINT.json',dict(event=v,RNG=rng.bit_generator.state))
            if v['event']=='STAGE_COMPLETE':print(json.dumps(dict(parent=args.parent,replicate=args.replicate,beta=v['beta'],stage=v['stage'],ancestors=v['history']['distinct_prior_ancestors'],acceptance=v['history']['mutation']['acceptance'])),flush=True)
        result=tempered_smc(lambda rng,n:rng.normal(size=(n,66)) if saved_particles is None else saved_particles.copy(),lambda x:-.5*float(x@x)-33*np.log(2*np.pi),ll,rng=rng,N=args.particles,mutation_steps=2,mutation=args.mutation,gaussian_dimension=66,gradloglikelihood=gl if mobility else None,precision=precision if mobility else None,step_size=.01,proposal_scale=.02,cess_target=.8,resample_threshold=.5,max_stages=128,callback=callback)
        events.close()
        save(b.jobdir/'PARTICLES.npz',**{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
        summary={k:v for k,v in result.items() if not isinstance(v,np.ndarray)}
        summary.update(target_family='FULL',coordinate='exact probit Gaussian66 pushforward of original Gaussian64 x uniform2',GT_used=False,loglikelihood_constant_omitted=-256*np.log(np.pi*t.sigma**2),no_reference_qualification_from_weight_ESS=True)
        write(b.jobdir/'RESULT.json',summary)
    except BaseException as e:status='ERROR';error=repr(e);raise
    finally:b.close(status,error)

if __name__=='__main__':main()
