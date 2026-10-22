"""Parent-queued frozen versus true Full-J metric MALA at explicitly matched RHS.

No solver work at import. Production epsilon is fixed. Dynamic endpoints cost
8 forward +32 receiver-J RHS; frozen endpoints cost8 forward +8 source-adjoint.
The same finite paid probit mode initializes both arms; it is not a posterior
sample or evidence of dispersed initialization.
"""
from pathlib import Path
import argparse,json,os,sys,time
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(name,'1')
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'tools'))
from p3r2.core import OUT,Book,Target,BudgetStop,write,save,to_probit,from_probit,probit_derivative
from p3r2.metric_mala import MetricState,FullMetricMALA,verify_probability_arrays
from p3l.diagnostics import material_diagnostics
from run_p3r2_posterior import finite,freeze


class PilotBook(Book):
    """Count actual atomic Book reservations, including failed physical attempts."""
    def __init__(self,*args,pilot_rhs_cap=262144,**kwargs):
        self.pilot_rhs_cap=int(pilot_rhs_cap);self.pilot_rhs=0;self.cost_owner='setup';self.costs={}
        super().__init__(*args,**kwargs)
    def _update(self,event,**fields):
        rhs=int(fields.get('rhs',0))
        if self.pilot_rhs+rhs>self.pilot_rhs_cap: raise BudgetStop('metric-MALA independent pilot RHS cap')
        result=super()._update(event,pilot_rhs_before=self.pilot_rhs,pilot_rhs_cap=self.pilot_rhs_cap,**fields)
        self.pilot_rhs+=rhs;self.costs[self.cost_owner]=self.costs.get(self.cost_owner,0)+rhs
        return result


def material_summary(chi):
    d=material_diagnostics(chi);ess=d['mean_ess'];bulk=d['bulk_ess'];tail=d['tail_ess'];rh=d['rhat'];mc=d['mean_mcse']
    f=lambda a,fun:float(fun(a[np.isfinite(a)])) if np.isfinite(a).any() else None
    return d,dict(mean_ESS_minimum_finite=f(ess,np.min),bulk_ESS_minimum_finite=f(bulk,np.min),
        tail_ESS_minimum_finite=f(tail,np.min),rhat_maximum_finite=f(rh,np.max),
        undefined_mean_ESS_functions=int((~np.isfinite(ess)).sum()),undefined_bulk_ESS_functions=int((~np.isfinite(bulk)).sum()),
        undefined_tail_ESS_functions=int((~np.isfinite(tail)).sum()),undefined_rhat_functions=int((~np.isfinite(rh)).sum()),
        complex_material_mean_RMS_MCSE=float(np.sqrt(np.mean(np.sum(mc**2,axis=-1)))),
        constant_unassessed_functions=int((d['status']=='CONSTANT_UNASSESSED').sum()),
        chains_stuck_functions=int((d['status']=='CHAINS_STUCK_AT_DIFFERENT_VALUES').sum()),
        production_draws_per_chain=chi.shape[1],chains=chi.shape[0],
        diagnostic_scope='Asymptotic material-function diagnostics; constants/undefined functions remain unassessed; no global mixing certificate')


def schedule(args):
    dynamic=dict(warmup=args.warmup,production=args.production,RHS_per_transition=40)
    if args.comparison=='equal-rhs':
        if (args.warmup*40)%16 or (args.production*40)%16: raise ValueError('equal-RHS schedule needs even dynamic warmup/production counts')
        frozen=dict(warmup=args.warmup*40//16,production=args.production*40//16,RHS_per_transition=16)
    else: frozen=dict(warmup=args.warmup,production=args.production,RHS_per_transition=16)
    return dict(dynamic=dynamic,frozen=frozen)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job');parser.add_argument('--initial-cache',type=Path)
    parser.add_argument('--verify-only',action='store_true');parser.add_argument('--validation-output',type=Path,default=OUT/'B/metric_mala_array_validation/RESULT.json')
    parser.add_argument('--parent',default='DEV_030');parser.add_argument('--arm',choices=('both','dynamic','frozen'),default='both')
    parser.add_argument('--comparison',choices=('equal-rhs','equal-steps'),default='equal-rhs')
    parser.add_argument('--chains',type=int,default=4);parser.add_argument('--warmup',type=int,default=128);parser.add_argument('--production',type=int,default=512)
    parser.add_argument('--eps',type=float,default=.05);parser.add_argument('--eps-dynamic',type=float);parser.add_argument('--eps-frozen',type=float)
    parser.add_argument('--rhs-cap',type=int,default=262144);parser.add_argument('--max-wall',type=int,default=7200)
    parser.add_argument('--seed',type=int,default=2026102391);parser.add_argument('--save-full-j',action='store_true')
    args=parser.parse_args()
    if args.verify_only:
        result=verify_probability_arrays();freeze(args.validation_output,result);print(json.dumps(finite(result),indent=2));return
    if not args.job or args.initial_cache is None: parser.error('--job and --initial-cache are required for a paid pilot')
    if args.rhs_cap>262144 or args.chains<2 or args.production<4 or args.warmup<0: raise ValueError('bounded multi-chain pilot required')
    arms=('dynamic','frozen') if args.arm=='both' else (args.arm,);plans=schedule(args)
    with np.load(args.initial_cache,allow_pickle=False) as z:
        q0=z['q'].copy();hasT='receiver_T' in z.files and bool(z['receiver_T'].size)
        if np.any(q0[64:]<=0) or np.any(q0[64:]>=1): raise ValueError('finite paid probit mode required; no implicit inward reinitialization')
    maximum_rhs=sum(args.chains*(plans[a]['warmup']+plans[a]['production'])*plans[a]['RHS_per_transition'] for a in arms)+(0 if hasT else 32)
    if maximum_rhs>args.rhs_cap: raise ValueError(f'planned total {maximum_rhs} RHS exceeds pilot cap {args.rhs_cap}')
    eps={a:float(args.eps if getattr(args,'eps_'+a) is None else getattr(args,'eps_'+a)) for a in ('dynamic','frozen')}
    if any(not np.isfinite(v) or v<=0 for v in eps.values()): raise ValueError('positive frozen epsilon required')
    book=PilotBook('B',args.job,max_wall=args.max_wall,pilot_rhs_cap=args.rhs_cap)
    status='COMPLETE';error=None;stream=None;armresults={};allchi={};start=time.perf_counter()
    try:
        dest=book.jobdir;target=Target(book);y=target.observation(args.parent)
        protocol=dict(job=args.job,target_family='FULL66_ORIGINAL_PROBIT_GAUSSIAN_PRIOR',target_id=target.cfg['target_id'],parent=args.parent,
            initial_cache=str(args.initial_cache.resolve()),initialization='Same already-paid finite probit mode for every chain and arm; no new optimization',
            shared_initialization_is_not_overdispersed=True,initialization_is_not_posterior_draws=True,GT_used=False,source_isolated_online_observation=True,
            dynamic_metric='H(x)=I+Jx(x).T@Jx(x), Jx=TRUE FullMaxwell receiver-J times exact probit derivative; C=H^-1',
            frozen_metric='H at the common paid initial anchor; subsequent endpoint gradient from true source adjoint, no endpoint J',
            drift='-0.5*eps^2*C(x)*grad U(x)',noise_covariance='eps^2*C(x)',
            exact_MH='log pi(y)-log pi(x)+log N(x;mu(y),eps²C(y))-log N(y;mu(x),eps²C(x)); both log determinants included',
            Gamma_term='Not added; exact MH corrects the explicitly evaluated Gaussian proposal',
            target_energy='U=Phi(FullMaxwell(from_probit(x)),y)+0.5*x@x; original exact probit pushforward prior',
            arms=arms,comparison=args.comparison,schedules=plans,eps=eps,eps_fixed_during_warmup_and_production=True,
            new_setup_RHS_maximum=0 if hasT else 32,planned_total_new_RHS_maximum=maximum_rhs,pilot_RHS_cap=args.rhs_cap,
            dynamic_RHS_per_endpoint=40,frozen_RHS_per_endpoint=16,costs='Actual atomic Book reservations; setup and warmup separately disclosed; all rejected endpoints billed',
            coupled_rng='Identical SeedSequence([seed,chain]) for both arms; same initial normal/uniform stream, different trajectories and draw counts',
            seed=args.seed,chains=args.chains,save_full_J=args.save_full_j,posterior_qualification=False,reference_qualification=False,
            failure_policy='A nonfinite state/gradient/metric stops as unqualified numerical failure; no silent clipping, metric flooring, or free rejection')
        freeze(dest/'PROTOCOL.json',protocol)
        anchor=target.load_cache(args.initial_cache,owner='Historical P3L paid DEV030 finite probit mode; online-y only, old ledger immutable')
        x0=to_probit(anchor.q);J0=target.jacobian(anchor)*probit_derivative(x0)[None,:]
        H0=np.eye(66)+J0.T@J0;H0=(H0+H0.T)/2
        gradient0=x0+target.phi_gradient(anchor,y)*probit_derivative(x0);Phi0=target.phi(anchor,y)
        field0=target.prior.decode(anchor.q)[1]['field'];stratum0=np.where(field0<=-target.prior.temperature,-1,np.where(field0>=target.prior.temperature,1,0)).astype(np.int8)
        initial=MetricState(x0,Phi0+.5*float(x0@x0),gradient0,H0,
            dict(point=anchor,Phi=Phi0,chi=anchor.state.chi.copy(),stratum=stratum0,endpoint_key=anchor.key))
        save(dest/'INITIAL_STATE_AND_FROZEN_METRIC.npz',x=x0,q=anchor.q,Jx=J0,H=H0,gradientU=gradient0,Phi=Phi0,
             stratum=stratum0,material=anchor.state.chi,probit_roundtrip_native_q_error=np.linalg.norm(from_probit(x0)-anchor.q))
        write(dest/'SETUP_RECEIPT.json',dict(actual_new_RHS=book.pilot_rhs,paid_receiver_T=hasT,load_cache_checks='Target.load_cache verifies unchanged decoder and Maxwell equations',initial_Phi=Phi0,initial_U=initial.energy))
        stream=(dest/'TRANSITIONS.jsonl').open('x')
        for arm in arms:
            plan=plans[arm];warmup=plan['warmup'];production=plan['production'];cost=plan['RHS_per_transition'];armdir=dest/arm;armdir.mkdir()
            arm_start=book.pilot_rhs;chains=[];rows=[];production_chi=[]
            def evaluate(x):
                if book.pilot_rhs+cost>args.rhs_cap: raise BudgetStop('stop before complete metric-MALA endpoint; local RHS cap')
                q=from_probit(x);p=target.evaluate(q,purpose='P3R2_'+arm+'_TRUE_FULL_METRIC_MALA')
                if arm=='dynamic':
                    Jx=target.jacobian(p)*probit_derivative(x)[None,:];H=np.eye(66)+Jx.T@Jx;H=(H+H.T)/2
                else: Jx=None;H=H0
                gradient=x+target.phi_gradient(p,y)*probit_derivative(x);Phi=target.phi(p,y)
                _,decoded=target.prior.decode(q);field=decoded['field'];stratum=np.where(field<=-target.prior.temperature,-1,np.where(field>=target.prior.temperature,1,0)).astype(np.int8)
                payload=dict(point=p,Phi=Phi,chi=p.state.chi.copy(),stratum=stratum,endpoint_key=p.key)
                state=MetricState(x,Phi+.5*float(x@x),gradient,H,payload)
                arrays=dict(x=x,q=q,Phi=Phi,U=state.energy,gradientU=gradient,stratum=stratum,
                    precision=H if arm=='dynamic' else np.zeros((0,0)),logdet_precision=state.logdet,
                    receiver_Jx_Frobenius=np.nan if Jx is None else np.linalg.norm(Jx))
                if Jx is not None and args.save_full_j: arrays['Jx']=Jx
                save(armdir/'endpoint_geometry'/f'{p.key}.npz',**arrays)
                # Current/proposed MetricState retains exactly the Point it
                # needs; no historical native LU cache is accumulated.
                if target.model is not None: target.model._transfers.clear()
                return state
            kernel=FullMetricMALA(evaluate,eps[arm])
            for chain in range(args.chains):
                rng=np.random.default_rng(np.random.SeedSequence([args.seed,chain]));current=initial
                states=[];phis=[];energy=[];accepted=[];keys=[];metadata=[];materials=[];charged=[];proposed_strata=[]
                chain_start=book.pilot_rhs
                for iteration in range(warmup+production):
                    phase='warmup' if iteration<warmup else 'production';book.cost_owner=f'{arm}/c{chain}/{phase}'
                    before=book.pilot_rhs;old=current
                    current,proposed,row=kernel.step(current,rng);rhs=book.pilot_rhs-before
                    if rhs!=cost: raise RuntimeError(f'{arm} endpoint actual RHS {rhs} != declared {cost}')
                    oldp=old.payload;newp=proposed.payload
                    RMS=float(np.linalg.norm(newp['chi']-oldp['chi'])/24)
                    strata_changed=(oldp['stratum']!=newp['stratum'])&target.prior.mask
                    row.update(arm=arm,chain=chain,iteration=iteration,phase=phase,proposal_endpoint=newp['endpoint_key'],state_endpoint=current.payload['endpoint_key'],
                        Phi_current=oldp['Phi'],Phi_proposed=newp['Phi'],Phi_retained=current.payload['Phi'],
                        delta_Phi=newp['Phi']-oldp['Phi'],delta_prior_energy=.5*float(proposed.x@proposed.x-old.x@old.x),
                        proposed_material_RMS=RMS,accepted_material_RMS=RMS if row['accepted'] else 0.,
                        decoder_stratum_changed_ROI_cells=int(strata_changed.sum()),
                        proposed_host_plateau_ROI_fraction=float(np.mean(newp['stratum'][target.prior.mask]==-1)),
                        proposed_pure_plateau_ROI_fraction=float(np.mean(newp['stratum'][target.prior.mask]==1)),
                        proposed_transition_ROI_fraction=float(np.mean(newp['stratum'][target.prior.mask]==0)),
                        actual_endpoint_RHS=rhs,pilot_cumulative_RHS=book.pilot_rhs,
                        arm_chain_cumulative_RHS=book.pilot_rhs-chain_start,
                        arm_chain_production_RHS=max(0,(iteration-warmup+1))*cost)
                    stream.write(json.dumps(finite(row),allow_nan=False,separators=(',',':'))+'\n');stream.flush()
                    states.append(current.x.copy());phis.append(current.payload['Phi']);energy.append(current.energy);accepted.append(row['accepted']);keys.append(current.payload['endpoint_key'])
                    metadata.append(row);materials.append(current.payload['chi'].copy());charged.append(rhs);proposed_strata.append(newp['stratum'].copy())
                    if (iteration+1)%64==0 or iteration+1==warmup+production:
                        chaindir=armdir/f'chain_{chain}';chaindir.mkdir(exist_ok=True)
                        save(chaindir/'CHAIN.npz',x=np.asarray(states),Phi=np.asarray(phis),U=np.asarray(energy),accepted=np.asarray(accepted),endpoint_key=np.asarray(keys),
                            actual_RHS_per_step=np.asarray(charged),proposed_stratum=np.asarray(proposed_strata),warmup=warmup,
                            proposed_material_RMS=np.array([r['proposed_material_RMS'] for r in metadata]),
                            accepted_material_RMS=np.array([r['accepted_material_RMS'] for r in metadata]),
                            forward_quadratic=np.array([r['forward_quadratic'] for r in metadata]),reverse_quadratic=np.array([r['reverse_quadratic'] for r in metadata]),
                            log_mh_ratio=np.array([r['log_mh_ratio'] for r in metadata]))
                        target.persist_cache(current.payload['point'],chaindir/'CURRENT_PAID_CACHE.npz')
                        write(chaindir/'CHECKPOINT.json',dict(iteration=iteration+1,phase=phase,warmup=warmup,production=production,
                            RNG=rng.bit_generator.state,actual_RHS=book.pilot_rhs-chain_start,eps=eps[arm],production_eps_fixed=True,
                            every_rejection_retained=True,current_Phi=current.payload['Phi'],posterior_qualification=False))
                        print(json.dumps(dict(arm=arm,chain=chain,iteration=iteration+1,phase=phase,
                            actual_chain_RHS=book.pilot_rhs-chain_start,Phi=current.payload['Phi'],acceptance_last64=float(np.mean(accepted[-64:])))),flush=True)
                    del old,proposed
                chi=np.asarray(materials[warmup:]);production_chi.append(chi)
                chainrow=dict(chain=chain,actual_RHS=book.pilot_rhs-chain_start,warmup_RHS=warmup*cost,production_RHS=production*cost,
                    production_acceptance=float(np.mean(accepted[warmup:])),warmup=warmup,production=production,
                    proposed_material_RMS_median=float(np.median([r['proposed_material_RMS'] for r in metadata[warmup:]])),
                    accepted_material_RMS_mean=float(np.mean([r['accepted_material_RMS'] for r in metadata[warmup:]])))
                rows.append(chainrow);chains.append(np.asarray(states));del current,materials
            chi=np.asarray(production_chi);allchi[arm]=chi
            save(armdir/'PRODUCTION_MATERIAL.npz',chi=chi)
            diagnostics,summary=material_summary(chi);save(armdir/'MATERIAL_DIAGNOSTICS.npz',**{k:v for k,v in diagnostics.items() if isinstance(v,np.ndarray)})
            result=dict(arm=arm,actual_RHS=book.pilot_rhs-arm_start,chains=rows,material_functions=summary,
                production_RHS_total=args.chains*production*cost,warmup_RHS_total=args.chains*warmup*cost,
                posterior_qualification=False,reference_qualification=False,
                route='TRUE_STATE_DEPENDENT_FULL_J_GN' if arm=='dynamic' else 'FROZEN_ANCHOR_GN_SOURCE_ADJOINT_GRADIENT')
            armresults[arm]=result;write(armdir/'RESULT.json',finite(result))
        checkpoints=[]
        if len(arms)==2:
            maximum=min(allchi[a].shape[1]*plans[a]['RHS_per_transition'] for a in arms)
            # Multiples of lcm(16,40)=80 give exactly matching paid production
            # work per chain. Warmup and setup remain explicitly separate.
            budgets=sorted(set(maximum//(4*80)*80*i for i in (1,2,3,4)))
            for budget in budgets:
                if budget<=0: continue
                comparison=dict(production_RHS_per_chain=budget,arms={})
                for arm in arms:
                    n=budget//plans[arm]['RHS_per_transition']
                    if n<4: continue
                    _,ds=material_summary(allchi[arm][:,:n]);comparison['arms'][arm]=dict(draws_per_chain=n,actual_production_RHS_per_chain=n*plans[arm]['RHS_per_transition'],material_functions=ds)
                checkpoints.append(comparison)
            dynamic=allchi['dynamic'].mean(axis=(0,1));frozen=allchi['frozen'].mean(axis=(0,1))
            difference=float(np.linalg.norm(dynamic-frozen)/24)
        else: difference=None
        summary=dict(protocol=protocol,status='FULL_METRIC_MALA_COMPARISON_COMPLETE_UNQUALIFIED',actual_new_RHS=book.pilot_rhs,
            cost_receipt_by_arm_chain_phase=book.costs,arms=armresults,matched_production_RHS_checkpoints=checkpoints,
            material_mean_RMS_difference_between_arms=difference,elapsed_seconds=time.perf_counter()-start,
            posterior_qualification=False,reference_qualification=False,
            scope='True Full-J metric versus frozen metric; acceptance alone is not a mixing result. Shared mode starts weaken convergence diagnostics.')
        write(dest/'RESULT.json',finite(summary));print(json.dumps(finite(summary)),flush=True)
    except BaseException as exc:
        status='ERROR';error=repr(exc);write(book.jobdir/'ERROR.json',dict(error=error,actual_new_RHS=book.pilot_rhs,costs=book.costs,all_completed_endpoints_and_transition_receipts_preserved=True));raise
    finally:
        if stream is not None: stream.close()
        book.cost_owner='closure';book.close(status,error)


if __name__=='__main__': main()
