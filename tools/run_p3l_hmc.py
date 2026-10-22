"""Independent full-target HMC; serial caller-owned physical accounting."""
from pathlib import Path
import sys,argparse,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,Target,write,save,from_probit
from p3l.reference import Posterior,metric
from p3l.kernels import ChainState
from p3l.hmc import PreconditionedHMC

def main():
    a=argparse.ArgumentParser();a.add_argument('--chain',type=int,required=True);a.add_argument('--warmup',type=int,default=512);a.add_argument('--production',type=int,default=1024);a.add_argument('--version',default='pilot');args=a.parse_args()
    parent='DEV_030';job=f'hmc_{parent}_{args.version}_c{args.chain}'
    b=Book('A',job,max_wall=7200);status='COMPLETE';error=None
    try:
        t=Target(b);y=t.observation(parent);post=Posterior(t,y,probit=True)
        source=OUT/f'A/jobs/reference_{parent}_mala_mode_c{args.chain}/CHAIN.npz'
        with np.load(source,allow_pickle=False) as z:x=z['x'][0].copy()
        # Paid common initialization is reused, not called a posterior draw.
        # Rebuilding its physics state costs a fresh forward and receiver-J.
        p=t.evaluate(from_probit(x),'HMC_state_rebuild_from_paid_y_only_initialization')
        post.cache[x.tobytes()]=p;post.pin(x);C=metric(t,p,x)
        current=ChainState(x,post(x),post.gradient(x))
        rng=np.random.default_rng(np.random.SeedSequence([2026102290,args.chain]))
        h=.05;mu=np.log(10*h);t_da=0;hbar=0.;logbar=np.log(h)
        warmup,production=args.warmup,args.production;total=warmup+production
        refresh=tuple(int(warmup*f) for f in [.25,.5,.75]);history=[];qs=[];xs=[];phis=[];lps=[];keys=[]
        write(b.jobdir/'PROTOCOL.json',dict(parent=parent,source_initialization=str(source),source_index=0,source_is_not_teacher=True,initialization_setup_must_be_disclosed=True,GT_used=False,warmup=warmup,production=production,mass_refresh_warmup=refresh,trajectory_length_choices=[4,8,12],production_parameters_fixed=True,coordinate='x=(z,eta), u=normal_CDF(eta)',original_prior_unchanged=True))
        for i in range(total):
            if i in refresh and i<warmup:
                C=metric(t,post.point(current.q),current.q)
                save(b.jobdir/f'MASS_WARMUP_{i}.npz',covariance=C,x=current.q)
                h=.05;mu=np.log(10*h);t_da=0;hbar=0.;logbar=np.log(h)
            length=int(rng.choice([4,8,12]));kernel=PreconditionedHMC(post,post.gradient,h,C,length)
            tr=kernel.step(current,rng);current=tr.state;post.pin(current.q)
            m=tr.metadata
            if i<warmup:
                ratio=m['log_mh_ratio'];alpha=0. if not np.isfinite(ratio) else np.exp(min(0.,ratio))
                t_da+=1;eta=1/(t_da+10);hbar=(1-eta)*hbar+eta*(.8-alpha)
                logh=mu-np.sqrt(t_da)/.05*hbar;weight=t_da**(-.75);logbar=weight*logh+(1-weight)*logbar
                h=float(np.clip(np.exp(np.clip(logbar if i==warmup-1 else logh,-20,0)),1e-6,1.))
            p=post.point(current.q);qs.append(from_probit(current.q));xs.append(current.q.copy());phis.append(t.phi(p,y));lps.append(current.log_density);keys.append(p.key);history.append(m)
            if (i+1)%128==0 or i+1==total:
                save(b.jobdir/'CHAIN.npz',q=np.asarray(qs),x=np.asarray(xs),Phi=np.asarray(phis),log_density=np.asarray(lps),accepted=np.array([v['accepted'] for v in history]),log_mh_ratio=np.array([v['log_mh_ratio'] for v in history]),log_hamiltonian_error=np.array([v['log_hamiltonian_error'] for v in history]),leapfrog_steps=np.array([v['leapfrog_steps'] for v in history]),gradient_evaluations=np.array([v['gradient_evaluations'] for v in history]),divergences=np.array([v['divergences'] for v in history]),large_energy_error=np.array([v.get('large_energy_error',False) for v in history]),endpoint_key=np.asarray(keys),warmup=np.array(warmup),target_family=np.array('FULL'))
                t.persist_cache(p,b.jobdir/'cache/CURRENT.npz')
                write(b.jobdir/'CHECKPOINT.json',dict(iteration=i+1,total=total,RNG=rng.bit_generator.state,step_size=h,production_parameters_frozen=i>=warmup,all_rejections_retained=True,acceptance_last128=float(np.mean([v['accepted'] for v in history[-128:]])),Phi=phis[-1]))
                print(json.dumps(dict(chain=args.chain,iteration=i+1,total=total,Phi=phis[-1],step_size=h,acceptance_last128=float(np.mean([v['accepted'] for v in history[-128:]])))),flush=True)
        save(b.jobdir/'MASS.npz',covariance=C)
        write(b.jobdir/'TRANSITION_SUMMARY.json',dict(target_family='FULL',method='fixed_dense_mass_HMC',production_acceptance=float(np.mean([v['accepted'] for v in history[warmup:]])),gradient_invocations=sum(v['gradient_evaluations'] for v in history),completed_leapfrog_steps=sum(v['leapfrog_steps'] for v in history),numerical_divergences=sum(v['divergences'] for v in history),large_energy_errors=sum(v.get('large_energy_error',False) for v in history),step_size=h,probability_status='INVARIANT_KERNEL_NOT_A_MIXING_OR_GLOBAL_REFERENCE_CERTIFICATE',GT_used=False))
    except BaseException as e:status='ERROR';error=repr(e);raise
    finally:b.close(status,error)

if __name__=='__main__':main()
