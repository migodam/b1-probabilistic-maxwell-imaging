"""Conventional surrogate-force HMC, exact Full endpoint Hamiltonian MH.

Explicit --arm and --chain select one serial Book-C job. No action on import.
No internal Full gradients, training, GT, or frozen src/p3l modifications.
"""
from pathlib import Path
import argparse
import importlib.util
import json
import shutil
import sys
import time
import numpy as np
from scipy.linalg import cho_factor,cho_solve

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,Target,write,save,from_probit,probit_derivative
from p3l.reference import Posterior
from p3l.hmc import PreconditionedHMC
from p3l.kernels import ChainState

MODULE_PATH=OUT/'C/surrogate_code.py'
module_spec=importlib.util.spec_from_file_location('p3l_C_frozen_surrogate',MODULE_PATH)
proxy=importlib.util.module_from_spec(module_spec);module_spec.loader.exec_module(proxy)
PARENT='DEV_030';WARMUP=1024;PRODUCTION=4096;REFRESH=(256,512,768)


def clean(value):
    if isinstance(value,dict):return {k:clean(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [clean(v) for v in value]
    if isinstance(value,(float,np.floating)) and not np.isfinite(value):return None
    if isinstance(value,np.ndarray):return clean(value.tolist())
    if isinstance(value,np.generic):return value.item()
    return value


def stage(book):
    d=book._latest();return dict(sequence=d['sequence'],**d['stages']['C'])


def cost(before,after):
    return {k:after[k]-before[k] for k in
        ('RHS','states','LU','controlled_wall_seconds','physics_wall_seconds','GPU_seconds')}


class CountedForce:
    def __init__(self):self.force=None;self.calls=0;self.seconds=0.
    def __call__(self,x):
        self.calls+=1;started=time.perf_counter()
        try:
            with np.errstate(over='ignore',invalid='ignore'):
                if not np.isfinite(x).all():return np.full(66,np.nan)
                return self.force.gradient(x)
        finally:self.seconds+=time.perf_counter()-started


class FullEndpoint:
    """Record the paid proposal endpoint; proxy diagnostics occur after MH."""
    def __init__(self,post):self.post=post;self.reset()
    def reset(self):self.x=None;self.point=None;self.log_density=None
    def __call__(self,x):
        value=self.post(x)
        self.x=np.asarray(x).copy();self.point=self.post.point(x);self.log_density=value
        return value


def setup(target,p,x,y,arm,folder):
    # This is the only receiver-J/T action; paid by Target if not cached on the
    # exact already-paid current Full state. No source-adjoint action exists.
    Jq=target.jacobian(p);Jx=Jq*probit_derivative(x)[None,:]
    H=Jx.T@Jx+np.eye(66)
    C=cho_solve(cho_factor(H,lower=True),np.eye(66));C=.5*(C+C.T)
    arrays=dict(x0=x,q0=p.q,chi0=p.state.chi,F0=p.state.data,Jx=Jx,H=H,
        covariance=C,y=y,sigma=np.array(target.sigma),endpoint_key=np.array(p.key))
    if arm=='material':
        Cchi=np.einsum('qp,ps->qsp',p.T,p.state.field).reshape(256,576)
        force=proxy.MaterialResponseForce(target.prior,x,p.state.chi,p.state.data,Cchi,y,target.sigma)
        arrays['Cchi']=Cchi
    else:force=proxy.LinearResponseForce(x,p.state.data,Jx,y,target.sigma)
    save(folder,**arrays)
    return force,C


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arm',choices=['linear','material'],required=True)
    parser.add_argument('--chain',type=int,choices=range(4),required=True)
    args=parser.parse_args();arm=args.arm;chain=args.chain
    source=OUT/f'A/jobs/reference_{PARENT}_mala_mode_c{chain}/CHAIN.npz'
    with np.load(source,allow_pickle=False) as archive:x0=archive['x'][0].copy()
    if x0.shape!=(66,) or not np.isfinite(x0).all():raise ValueError('finite common physical probit66 start required')
    book=Book('C',f'surrogate_force_{arm}_DEV030_c{chain}',max_wall=7200)
    status='ERROR';error=None;before=stage(book);history=[];force=None
    try:
        if (book.jobdir/'PROTOCOL.json').exists():raise FileExistsError('no overwrite or automatic replay')
        write(book.jobdir/'PROTOCOL.json',dict(parent=PARENT,arm=arm,chain=chain,
            method='CONVENTIONAL_FROZEN_SURROGATE_FORCE_HMC_WITH_FULL_ENDPOINT_HAMILTONIAN_MH',
            novelty_claim=False,scientific_mixing_certificate=False,
            dependency_baseline='ca2ec9, parent-declared frozen src/p3l; unchanged by this implementation',
            source_initialization=source.relative_to(OUT).as_posix(),source_index=0,initial_x=x0,
            common_initialization='original mala_mode first row, matching ordinary HMC pilot',
            matched_mechanism_comparison='linear versus material arms share1024warm/4096prod,16/32/64 lengths,seed2341+chain',
            ordinary_HMC_comparison='external cost/Pareto context only; not a matched sole-force ablation',
            ordinary_HMC_protocol=dict(warmup=512,production=1024,lengths=[4,8,12],
                seed='SeedSequence([2026102290,chain])'),
            comparison_claim_limit='no superiority from unequal warmup/lengths; compare statistical errors with actual total RHS and wall; no qualified-reference guarantee',
            previous_initialization_setup='paid A y-only setup shared and disclosed; initial Full rebuild paid C',
            warmup=WARMUP,production=PRODUCTION,trajectory_lengths=[16,32,64],
            seed=2026102341+chain,initial_step_size=.05,target_acceptance=.8,
            refresh=REFRESH,refresh_actions='true current Full endpoint J/T and GN inverse mass; force/current-gradient/dual-averaging reset',
            production_policy='surrogate, mass, and step size frozen',
            full_target='-Phi(from_probit(x))-.5*x@x in physical probit66 coordinates',
            material_response='Fhat=F0+(Cchi@(chi(q)-chi0)).reshape(32,8)',
            Cchi='einsum(qp,ps->qsp,T,field).reshape(256,576)',
            linear_force='-Jx.T@(white(F0-y)+Jx@(x-x0))-x',
            material_force='-2/sigma^2 Re(Dchi.conj().T@Cchi.conj().T@vec(Fhat-y))*probit_derivative(x)-x',
            internal_force='frozen proxy arithmetic only; no Full gradients or physics',
            endpoint='exact Full Hamiltonian MH including true likelihood and full Gaussian66 prior',
            proxy_error='computed after MH decision at already-paid Full proposal endpoint; no correction/selection',
            all_rejection_repeats_retained=True,checkpoint_every=128,max_wall_seconds=7200,
            RHS_upper_this_chain=(WARMUP+PRODUCTION+1)*8+4*32,
            RHS_upper_four_chains_per_arm=164384,RHS_upper_both_arms=328768,
            GT_read=False,training=False,serial_paid_worker=True))
        for rel in ['tools/run_p3l_surrogate_hmc.py','src/p3l/hmc.py','src/p3l/kernels.py',
            'src/p3l/reference.py','src/p3l/core.py','src/p1_material/prior.py']:
            dest=book.jobdir/'SOURCE'/rel;dest.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/rel,dest)
        shutil.copyfile(MODULE_PATH,book.jobdir/'SOURCE/surrogate_code.py')
        shutil.copyfile(OUT/'configs/TARGET.json',book.jobdir/'TARGET_FROZEN.json')
        shutil.copyfile(ROOT/'research/p1_material_target_20261018/inputs/geometry.json',book.jobdir/'GEOMETRY_FROZEN.json')
        target=Target(book);y=target.observation(PARENT);post=Posterior(target,y,probit=True)
        p=target.evaluate(from_probit(x0),'SURROGATE_HMC_FRESH_COMMON_INITIAL_STATE')
        post.cache[x0.tobytes()]=p;post.pin(x0)
        force=CountedForce();force.force,C=setup(target,p,x0,y,arm,book.jobdir/'SURROGATE_INITIAL.npz')
        current=ChainState(x0,post(x0),force(x0));full=FullEndpoint(post)
        rng=np.random.default_rng(2026102341+chain)
        h=.05;mu=np.log(10*h);da_t=0;hbar=0.;logbar=np.log(h)
        arrays={k:[] for k in ('q','x','Phi','log_density','endpoint_key','accepted','log_mh_ratio',
            'log_hamiltonian_error','requested_leapfrog_steps','leapfrog_steps','gradient_evaluations',
            'target_evaluations','divergences','large_energy_error','force_calls','proxy_Phi_error',
            'proxy_white_response_RMS','proposal_endpoint_key')}
        setup_index=0
        for i in range(WARMUP+PRODUCTION):
            if i in REFRESH:
                p=post.point(current.q)
                force.force,C=setup(target,p,current.q,y,arm,book.jobdir/f'SURROGATE_WARMUP_{i}.npz')
                # The cached gradient belongs to the previous frozen force.
                current=ChainState(current.q,current.log_density,force(current.q))
                setup_index+=1;h=.05;mu=np.log(10*h);da_t=0;hbar=0.;logbar=np.log(h)
            full.reset();calls_before=force.calls;force_seconds_before=force.seconds;draw_before=stage(book)
            length=int(rng.choice([16,32,64]))
            transition=PreconditionedHMC(full,force,h,C,length).step(current,rng)
            current=transition.state;post.pin(current.q);m=transition.metadata
            # All Full/proxy discrepancies are measured after the decision.
            discrepancy=None
            if full.point is not None:
                true_Phi=target.phi(full.point,y)
                discrepancy=proxy.response_error(force.force,full.x,full.point.state.data,true_Phi)
                save(book.jobdir/'PROPOSALS'/f'{i:05d}.npz',x=full.x,q=full.point.q,
                    Phi=np.array(true_Phi),log_density=np.array(full.log_density),
                    accepted=np.array(m['accepted']),endpoint_key=np.array(full.point.key),
                    surrogate_setup_index=np.array(setup_index),**discrepancy)
            if i<WARMUP:
                ratio=m['log_mh_ratio'];alpha=0. if not np.isfinite(ratio) else np.exp(min(0.,ratio))
                da_t+=1;eta=1/(da_t+10);hbar=(1-eta)*hbar+eta*(.8-alpha)
                logh=mu-np.sqrt(da_t)/.05*hbar;weight=da_t**(-.75)
                logbar=weight*logh+(1-weight)*logbar
                h=float(np.clip(np.exp(np.clip(logbar if i==WARMUP-1 else logh,-20,0)),1e-6,1.))
            p=post.point(current.q);Phi=target.phi(p,y)
            entries=dict(q=from_probit(current.q),x=current.q.copy(),Phi=Phi,
                log_density=current.log_density,endpoint_key=p.key,force_calls=force.calls-calls_before,
                proxy_Phi_error=np.nan if discrepancy is None else discrepancy['Phi_error'],
                proxy_white_response_RMS=np.nan if discrepancy is None else discrepancy['white_response_RMS'],
                proposal_endpoint_key='' if full.point is None else full.point.key)
            for k in arrays:
                arrays[k].append(entries[k] if k in entries else m[k])
            record=dict(iteration=i+1,phase='warmup' if i<WARMUP else 'production',
                transition=m,force_calls=force.calls-calls_before,total_force_calls=force.calls,
                proxy_force_seconds=force.seconds-force_seconds_before,
                surrogate_setup_index=setup_index,proxy_error_after_decision=discrepancy,
                cost=cost(draw_before,stage(book)),current_endpoint=p.key)
            history.append(record)
            path=book.jobdir/'TRANSITIONS.jsonl'
            with path.open('a') as stream:stream.write(json.dumps(clean(record),allow_nan=False)+'\n')
            if (i+1)%128==0:
                save(book.jobdir/'CHAIN.npz',**{k:np.asarray(v) for k,v in arrays.items()},
                    warmup=np.array(WARMUP),target_family=np.array('FULL'),arm=np.array(arm),
                    method=np.array('CONVENTIONAL_SURROGATE_FORCE_FULL_ENDPOINT_HMC'))
                target.persist_cache(p,book.jobdir/'cache/CURRENT.npz')
                write(book.jobdir/'CHECKPOINT.json',dict(iteration=i+1,total=WARMUP+PRODUCTION,
                    RNG=rng.bit_generator.state,step_size=h,current_Phi=Phi,
                    production_parameters_frozen=i>=WARMUP,total_proxy_force_calls=force.calls,
                    proxy_force_seconds=force.seconds,elapsed_wall_seconds=time.perf_counter()-book.start,
                    Full_gradient_calls=0,cost=cost(before,stage(book)),
                    acceptance_last128=float(np.mean(arrays['accepted'][-128:]))))
                print(json.dumps(dict(arm=arm,chain=chain,iteration=i+1,total=WARMUP+PRODUCTION,
                    Phi=Phi,step_size=h,acceptance_last128=float(np.mean(arrays['accepted'][-128:])))),flush=True)
        save(book.jobdir/'MASS.npz',covariance=C)
        write(book.jobdir/'TRANSITION_SUMMARY.json',dict(arm=arm,chain=chain,parent=PARENT,
            production_acceptance=float(np.mean(arrays['accepted'][WARMUP:])),
            total_proxy_force_calls=force.calls,Full_gradient_calls=0,
            proxy_force_seconds=force.seconds,elapsed_wall_seconds=time.perf_counter()-book.start,
            numerical_divergences=int(sum(arrays['divergences'])),
            large_energy_errors=int(sum(arrays['large_energy_error'])),
            cost=cost(before,stage(book)),step_size=h,GT_read=False,
            method='CONVENTIONAL_SURROGATE_FORCE_HMC',novelty_claim=False,
            probability_status='FULL_ENDPOINT_INVARIANT_KERNEL_NOT_A_MIXING_CERTIFICATE'))
        status='COMPLETE'
    except BaseException as exc:error=repr(exc);raise
    finally:
        if force is not None:
            book.counts['proxy_force_calls']=force.calls
            book.counts['proxy_force_seconds']=force.seconds
        book.close(status,error)


if __name__=='__main__':main()
