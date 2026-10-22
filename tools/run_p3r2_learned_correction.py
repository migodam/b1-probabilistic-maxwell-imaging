"""D2: exact FULL-target corrections for frozen learned proposal densities.

One explicit parent-queued job compares unit-prior independence, learned full
independence, and alternating learned-active/pCN-rest updates. Every proposal
uses the source-isolated online observation and the original native target.
No training arrays, truth labels, gradients or Jacobian queries are accessed.
"""
from pathlib import Path
import argparse,json,os,shutil,sys,time
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(name,'1')
import numpy as np
from scipy.special import ndtr,ndtri
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'tools'))
from p3r2.core import OUT,Book,Target,BudgetStop,write,save
from p3l.diagnostics import material_diagnostics
from run_p3r2_posterior import finite,freeze
LOG2PI=float(np.log(2*np.pi))
ARMS=('unit_prior_full','learned_full','learned_blocked')


class CorrectionBook(Book):
    def __init__(self,*args,pilot_rhs_cap=131072,**kwargs):
        self.pilot_rhs_cap=int(pilot_rhs_cap);self.pilot_rhs=0;self.cost_owner='setup';self.costs={}
        super().__init__(*args,**kwargs)
    def _update(self,event,**fields):
        rhs=int(fields.get('rhs',0))
        if self.pilot_rhs+rhs>self.pilot_rhs_cap: raise BudgetStop('D2 independent pilot RHS cap')
        result=super()._update(event,pilot_rhs_before=self.pilot_rhs,pilot_rhs_cap=self.pilot_rhs_cap,**fields)
        self.pilot_rhs+=rhs;self.costs[self.cost_owner]=self.costs.get(self.cost_owner,0)+rhs
        return result


def native_logprior(q): return -.5*float(np.asarray(q)[:64]@np.asarray(q)[:64])-32*LOG2PI


def labels_from_state(q,frame):
    q=np.asarray(q,float)
    if q.shape!=(66,) or not np.isfinite(q).all() or np.any(q[64:]<=0) or np.any(q[64:]>=1):
        raise ValueError('finite interior native state required; no boundary clipping')
    return np.r_[q[:64]@frame[:,:8],ndtri(q[64:])]


def block_logprob(density,y,q):
    labels=labels_from_state(q,density.frame)
    return float(density.log_prob_labels(y[None],labels[None])[0]+.5*np.sum(labels[8:]**2+LOG2PI))


def active_proposal(density,y,current,rng):
    seed=int(rng.integers(0,2**31-1));labels=density.sample_labels(y[None],1,seed)[0,0]
    coordinates=current[:64]@density.frame;coordinates[:8]=labels[:8]
    q=np.r_[coordinates@density.frame.T,ndtr(labels[8:])]
    labels_from_state(q,density.frame)
    return q,seed


def rest_proposal(current,frame,rng,beta):
    coordinates=current[:64]@frame;before=coordinates[8:].copy()
    coordinates[8:]=np.sqrt(1-beta**2)*before+beta*rng.normal(size=56)
    q=current.copy();q[:64]=coordinates@frame.T
    correction=.5*float(coordinates[8:]@coordinates[8:]-before@before)
    active_error=float(np.linalg.norm((q[:64]-current[:64])@frame[:,:8]))
    return q,correction,active_error


def _single_ess(values):
    values=np.asarray(values,float);n=len(values)
    if np.all(values==values[0]): return np.nan
    centered=values-values.mean();size=1<<(2*n-1).bit_length()
    f=np.fft.rfft(centered,n=size);cov=np.fft.irfft(f*np.conj(f),n=size)[:n]/n
    if cov[0]<=0: return np.nan
    rho=cov/cov[0];pairs=[]
    for index in range(0,n-1,2):
        value=float(rho[index]+rho[index+1])
        if value<0: break
        pairs.append(value if not pairs else min(value,pairs[-1]))
    tau=max(-1+2*sum(pairs),1/np.log10(n))
    return float(min(n/tau,n*np.log10(n)))


def summarize_material(chi):
    """Parent-specific functions only; different posteriors are never pooled."""
    chi=np.asarray(chi)
    if chi.shape[0]>=2:
        diagnostics=material_diagnostics(chi);ess=diagnostics['mean_ess'];mcse=diagnostics['mean_mcse'];rhat=diagnostics['rhat']
        constant=diagnostics['status']=='CONSTANT_UNASSESSED'
        rhat_status='MULTICHAIN_DIAGNOSTIC_WITH_SHARED_MODE_INITIALIZATION'
    else:
        a=np.stack((chi[0].real,chi[0].imag),axis=-1);flat=a.reshape(len(a),-1)
        ess=np.array([_single_ess(flat[:,k]) for k in range(flat.shape[1])]).reshape(576,2)
        sd=a.std(axis=0,ddof=1);constant=np.all(a==a[0],axis=0)
        with np.errstate(invalid='ignore',divide='ignore'): mcse=sd/np.sqrt(ess)
        mcse=np.where(constant,0.,mcse);rhat=np.full((576,2),np.nan)
        diagnostics=dict(mean=a.mean(axis=0),variance=sd**2,mean_ess=ess,mean_mcse=mcse,rhat=rhat,
                         status=np.where(constant,'CONSTANT_UNASSESSED','SINGLE_CHAIN_ASYMPTOTIC_DIAGNOSTIC_ONLY'))
        rhat_status='NOT_COMPUTED_SINGLE_CHAIN'
    finite_ess=ess[np.isfinite(ess)];finite_rhat=rhat[np.isfinite(rhat)]
    summary=dict(chains=chi.shape[0],retained_sweep_endpoints_per_chain=chi.shape[1],
        material_mean_ESS_minimum_finite=float(finite_ess.min()) if finite_ess.size else None,
        material_mean_ESS_median_finite=float(np.median(finite_ess)) if finite_ess.size else None,
        undefined_material_ESS_functions=int((~np.isfinite(ess)).sum()),constant_unassessed_material_functions=int(constant.sum()),
        complex_material_mean_RMS_MCSE=float(np.sqrt(np.mean(np.sum(mcse**2,axis=-1)))),
        rhat_status=rhat_status,rhat_maximum_finite=float(finite_rhat.max()) if finite_rhat.size else None,
        undefined_rhat_functions=int((~np.isfinite(rhat)).sum()),posterior_qualification=False,
        diagnostic_scope='Asymptotic material-function diagnostics; same paid-mode initialization is not overdispersion, ESS/MCSE cannot certify exploration')
    return diagnostics,summary


def probability_controls(seed=2026102394):
    rng=np.random.default_rng(seed);frame=np.linalg.qr(rng.normal(size=(64,64)))[0]
    full=[];blocked=[];pcn=[];native_jacobian=[];rest_errors=[]
    for _ in range(120):
        current=np.r_[rng.normal(size=64),rng.uniform(.02,.98,size=2)]
        proposed=np.r_[rng.normal(size=64),rng.uniform(.02,.98,size=2)]
        Phi=lambda q:float(.1*np.sum(q[:3]**2)+.7*q[64]+.3*q[65]**2)
        logtarget=lambda q:-Phi(q)-.5*float(q[:64]@q[:64])
        ratio=logtarget(proposed)-logtarget(current)+native_logprior(current)-native_logprior(proposed)
        full.append(abs(ratio-(-Phi(proposed)+Phi(current))))
        z=current[:64]@frame;lp=labels_from_state(current,frame);newlabels=rng.normal(size=10)
        z[:8]=newlabels[:8];q=np.r_[z@frame.T,ndtr(newlabels[8:])]
        def unit_block(a):
            label=labels_from_state(a,frame)
            return -.5*float(np.sum(label**2+LOG2PI))+.5*float(np.sum(label[8:]**2+LOG2PI))
        ratio=logtarget(q)-logtarget(current)+unit_block(current)-unit_block(q)
        blocked.append(abs(ratio-(-Phi(q)+Phi(current))))
        native_jacobian.append(abs(unit_block(current)-(-.5*float(lp[:8]@lp[:8])-4*LOG2PI)))
        q,correction,err=rest_proposal(current,frame,rng,.15)
        ratio=logtarget(q)-logtarget(current)+correction
        pcn.append(abs(ratio-(-Phi(q)+Phi(current))));rest_errors.append(err)
    return dict(status='PURE_ARRAY_UNIT_PRIOR_CANCELLATION_CONTROL_COMPLETE',new_physics_RHS=0,seed=seed,
        full_independence_prior_cancellation_max_error=float(max(full)),
        active_block_unit_prior_cancellation_max_error=float(max(blocked)),
        active_native_probit_Jacobian_cancellation_max_error=float(max(native_jacobian)),
        rest_pCN_prior_cancellation_max_error=float(max(pcn)),rest_pCN_active_coordinate_roundoff_max=float(max(rest_errors)),
        checks=dict(full_cancellation=max(full)<1e-10,active_cancellation=max(blocked)<1e-10,
                    probit_Jacobian=max(native_jacobian)<1e-10,rest_pCN_cancellation=max(pcn)<1e-10,
                    complementary_block_preserved=max(rest_errors)<1e-10),
        scope='Algebraic unit-prior controls only; no trained density quality, physics or posterior mixing qualification')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job');parser.add_argument('--selected',type=Path,default=OUT/'A/density_training/SELECTED.json')
    parser.add_argument('--parents',nargs='+',default=['DEV_000','DEV_026','DEV_028','DEV_030'])
    parser.add_argument('--cache-root',type=Path,default=ROOT.parent/'b1-p3l-local-posterior/research/p3l_local_posterior_imaging/B/joint_gaussian')
    parser.add_argument('--arms',nargs='+',choices=ARMS,default=list(ARMS));parser.add_argument('--chains',type=int,default=1)
    parser.add_argument('--warmup',type=int,default=128);parser.add_argument('--production',type=int,default=512)
    parser.add_argument('--pcn-beta',type=float,default=.15);parser.add_argument('--rhs-cap',type=int,default=131072)
    parser.add_argument('--seed',type=int,default=2026102395);parser.add_argument('--max-wall',type=int,default=7200)
    parser.add_argument('--verify-only',action='store_true');parser.add_argument('--validation-output',type=Path,default=OUT/'D/learned_correction_array_validation/RESULT.json')
    args=parser.parse_args()
    if args.verify_only:
        result=probability_controls();freeze(args.validation_output,result);print(json.dumps(finite(result),indent=2));return
    if not args.job: parser.error('--job is required for a paid pilot')
    if args.rhs_cap>131072 or args.chains<1 or args.production<8 or args.production%2 or args.warmup<0 or args.warmup%2 or not 0<args.pcn_beta<=1:
        raise ValueError('positive bounded pilot with even warmup/production transition counts required')
    if len(set(args.parents))!=len(args.parents) or len(set(args.arms))!=len(args.arms): raise ValueError('unique parents and arms required')
    maximum_rhs=len(args.parents)*len(args.arms)*args.chains*(args.warmup+args.production)*8
    if maximum_rhs>args.rhs_cap: raise ValueError(f'planned {maximum_rhs} RHS exceeds local pilot cap')
    selected=json.loads(args.selected.read_text());checkpoint=Path(selected['checkpoint'])
    if not checkpoint.is_absolute(): checkpoint=args.selected.parent/checkpoint
    if not checkpoint.is_file(): raise ValueError('frozen selected checkpoint missing; no implicit training or model substitution')
    paths={p:args.cache_root/p/'RAW_MODE_CACHE.npz' for p in args.parents}
    for p,path in paths.items():
        with np.load(path,allow_pickle=False) as z:
            q=z['q']
            if q.shape!=(66,) or np.any(q[64:]<=0) or np.any(q[64:]>=1): raise ValueError('paid finite probit mode required for '+p)
    book=CorrectionBook('D',args.job,max_wall=args.max_wall,pilot_rhs_cap=args.rhs_cap)
    status='COMPLETE';error=None;stream=None;start=time.perf_counter();results=[]
    try:
        dest=book.jobdir
        shutil.copy2(checkpoint,dest/'FROZEN_PROPOSAL.pt');write(dest/'SELECTED_SNAPSHOT.json',selected)
        from p3r2.density import load_density
        import torch
        torch.set_num_threads(1)
        if torch.get_num_interop_threads()!=1: torch.set_num_interop_threads(1)
        density=load_density(dest/'FROZEN_PROPOSAL.pt',device='cpu');frame=density.frame
        if frame.shape!=(64,64) or np.linalg.norm(frame.T@frame-np.eye(64))>1e-8: raise ValueError('frozen complete orthonormal Gaussian frame required')
        selection_result=selected.get('result',{});kind=density.p['kind'];update=int(density.p.get('best_update',0))
        flags=dict(kind=kind,selected_update=update,selected_initialization=selection_result.get('selected_initialization'),
            selected_trained_update_positive=bool(kind!='unit_prior' and update>0),unit_prior_selected=kind=='unit_prior',
            density_is_only_a_proposal=True,posterior_fidelity='NOT_ESTABLISHED',selection_criterion=selected.get('criterion'))
        target=Target(book)
        protocol=dict(job=args.job,stage='D',parents=args.parents,arms=args.arms,chains=args.chains,
            target_family='FULL66_ORIGINAL_NATIVE_GAUSSIAN64_UNIFORM2',target_id=target.cfg['target_id'],
            target_log_density='-Phi(FullMaxwell(q),online_y)-0.5*z64@z64 with respect to dz64 du2',
            GT_used_online=False,training_data_read=False,online_inputs='Source-isolated observation only plus already-paid y-only initialization caches',
            selected_json=str(args.selected.resolve()),selected_checkpoint=str(checkpoint.resolve()),frozen_checkpoint_copy=str(dest/'FROZEN_PROPOSAL.pt'),
            selected_density_flags=flags,frame='Frozen checkpoint first8 Gaussian directions and orthogonal prior56 complement; no state-dependent frame',
            learned_full_density='learned10(labels8,eta2|y) * normal_prior56 * abs(deta/du); normalized native66 proposal',
            learned_active_density='log_prob_labels(y,[frame8.T z,eta]) + 0.5*sum(eta²+log(2pi)); complement unchanged',
            prior_baseline='Normalized independent N(0,I64) x Uniform2; exact prior cancellation gives likelihood-only MH',
            learned_full_MH='Full native target difference + full native proposal logq(current)-logq(proposed)',
            active_MH='Full native66 target difference + native active-block logq(current)-logq(proposed)',
            rest_pCN_MH='Gaussian56 pCN in frozen complement; active8 and properties fixed; prior cancellation leaves -deltaPhi',
            block_schedule='Alternating learned active then pCN rest, fixed throughout warmup and production; no state-dependent mixture weights',
            transition_counting='Each paid block counts once. All arms have equal paid transitions/RHS. Blocked arm reports full sweep endpoints for ESS; all intermediate/rejected states retained.',
            warmup_paid_transitions=args.warmup,production_paid_transitions=args.production,
            blocked_warmup_sweeps=args.warmup//2,blocked_production_sweeps=args.production//2,
            pcn_beta=args.pcn_beta,seed=args.seed,maximum_new_RHS=maximum_rhs,pilot_RHS_cap=args.rhs_cap,
            actual_cost='8 FullMaxwell RHS per transition; no gradients, receiver-J or unbilled solver queries',
            initialization_cache_paths={k:str(v.resolve()) for k,v in paths.items()},
            parameters_and_density_frozen_before_production=True,same_mode_starts_not_overdispersed=True,
            full_target_exact_invariance_is_not_finite_run_mixing=True,reference_qualification=False,posterior_qualification=False)
        freeze(dest/'PROTOCOL.json',protocol);save(dest/'FIXED_COORDINATE_FRAME.npz',frame=frame)
        controls=probability_controls(args.seed);write(dest/'UNIT_PRIOR_ARRAY_CONTROLS.json',controls)
        if not all(controls['checks'].values()): raise RuntimeError('unit-prior probability control failed')
        stream=(dest/'TRANSITIONS.jsonl').open('x')
        for parent in args.parents:
            y=target.observation(parent)
            book.cost_owner=f'{parent}/initialization';anchor=target.load_cache(paths[parent],owner='Historical P3L paid y-only finite probit mode, old ledger immutable')
            initial=dict(q=anchor.q.copy(),Phi=target.phi(anchor,y),chi=anchor.state.chi.copy(),key=anchor.key)
            initial['prior_energy']=.5*float(initial['q'][:64]@initial['q'][:64]);initial['logtarget']=-initial['Phi']-initial['prior_energy']
            del anchor
            if target.model is not None: target.model._transfers.clear()
            # Verify the learned full proposal decomposes into the exact
            # active native density and independent Gaussian complement.
            lpfull=float(density.log_prob(y[None],initial['q'][None])[0]);lpblock=block_logprob(density,y,initial['q'])
            rest=initial['q'][:64]@frame[:,8:];decomposition_error=abs(lpfull-(lpblock-.5*np.sum(rest**2+LOG2PI)))
            if not np.isfinite(lpfull) or not np.isfinite(lpblock) or decomposition_error>1e-8: raise RuntimeError('learned normalized native density decomposition mismatch')
            write(dest/parent/'DENSITY_DECOMPOSITION_CHECK.json',dict(full_logq=lpfull,active_native_logq=lpblock,absolute_error=decomposition_error,GT_used=False,selected_density_flags=flags))
            for arm in args.arms:
                armdir=dest/parent/arm;armdir.mkdir(parents=True,exist_ok=True);chi_chains=[];chain_rows=[];arm_start=book.pilot_rhs
                for chain in range(args.chains):
                    rng=np.random.default_rng(np.random.SeedSequence([args.seed,int(parent.split('_')[-1]),chain,ARMS.index(arm)]))
                    current={k:(v.copy() if isinstance(v,np.ndarray) else v) for k,v in initial.items()}
                    qs=[];phis=[];priors=[];logtargets=[];accepted=[];blocks=[];keys=[];materials=[];RMS=[];ratios=[];chain_start=book.pilot_rhs
                    for transition in range(args.warmup+args.production):
                        phase='warmup' if transition<args.warmup else 'production';book.cost_owner=f'{parent}/{arm}/c{chain}/{phase}'
                        old=current;seed=None;correction=None;complement_error=0.;active_error=0.
                        if arm=='unit_prior_full':
                            block='full_prior';q=np.r_[rng.normal(size=64),rng.uniform(size=2)]
                            logq_current=native_logprior(old['q']);logq_proposed=native_logprior(q)
                        elif arm=='learned_full':
                            block='full_learned';seed=int(rng.integers(0,2**31-1));q=density.sample(y[None],1,seed)[0,0]
                            logq_current=float(density.log_prob(y[None],old['q'][None])[0]);logq_proposed=float(density.log_prob(y[None],q[None])[0])
                        elif transition%2==0:
                            block='learned_active8_probit2';q,seed=active_proposal(density,y,old['q'],rng)
                            logq_current=block_logprob(density,y,old['q']);logq_proposed=block_logprob(density,y,q)
                            complement_error=float(np.linalg.norm((q[:64]-old['q'][:64])@frame[:,8:]))
                        else:
                            block='pcn_prior_complement56';q,correction,active_error=rest_proposal(old['q'],frame,rng,args.pcn_beta)
                            logq_current=logq_proposed=None
                        labels_from_state(q,frame)
                        if complement_error>1e-10 or active_error>1e-10: raise RuntimeError('fixed block preservation failed')
                        if logq_current is not None and not np.isfinite(logq_current+logq_proposed): raise FloatingPointError('nonfinite learned native proposal density')
                        if book.pilot_rhs+8>args.rhs_cap: raise BudgetStop('stop before complete D2 endpoint due local RHS cap')
                        before=book.pilot_rhs;p=target.evaluate(q,purpose='D2_'+arm+'_'+block)
                        proposed=dict(q=q.copy(),Phi=target.phi(p,y),chi=p.state.chi.copy(),key=p.key)
                        proposed['prior_energy']=.5*float(q[:64]@q[:64]);proposed['logtarget']=-proposed['Phi']-proposed['prior_energy']
                        del p
                        if target.model is not None: target.model._transfers.clear()
                        actual=book.pilot_rhs-before
                        if actual!=8: raise RuntimeError('D2 full endpoint actual RHS differs from declared8')
                        if not np.isfinite(proposed['logtarget']): raise FloatingPointError('nonfinite exact Full target energy')
                        target_ratio=proposed['logtarget']-old['logtarget']
                        if block=='pcn_prior_complement56':
                            ratio=-proposed['Phi']+old['Phi'];proposal_ratio=correction
                            cancellation_error=abs((target_ratio+correction)-ratio)
                        else:
                            proposal_ratio=logq_current-logq_proposed;ratio=target_ratio+proposal_ratio
                            cancellation_error=abs(ratio-(-proposed['Phi']+old['Phi'])) if arm=='unit_prior_full' else None
                        uniform=float(rng.random());yes=bool(np.log(max(uniform,np.finfo(float).tiny))<min(0.,ratio))
                        current=proposed if yes else old;material_RMS=float(np.linalg.norm(proposed['chi']-old['chi'])/24)
                        row=dict(parent=parent,arm=arm,chain=chain,transition=transition,phase=phase,block=block,
                            proposal_seed=seed,accepted=yes,repeated_state=not yes,current_q=old['q'],proposed_q=q,retained_q=current['q'],
                            current_Phi=old['Phi'],proposed_Phi=proposed['Phi'],retained_Phi=current['Phi'],
                            current_prior_energy=old['prior_energy'],proposed_prior_energy=proposed['prior_energy'],
                            current_logtarget=old['logtarget'],proposed_logtarget=proposed['logtarget'],
                            log_proposal_current=logq_current,log_proposal_proposed=logq_proposed,
                            log_target_ratio=target_ratio,log_proposal_correction=proposal_ratio,log_mh_ratio=ratio,
                            prior_cancellation_error=cancellation_error,unchanged_complement_roundoff=complement_error,
                            unchanged_active_roundoff=active_error,proposal_material_RMS=material_RMS,accepted_material_RMS=material_RMS if yes else 0.,
                            current_endpoint=old['key'],proposed_endpoint=proposed['key'],retained_endpoint=current['key'],
                            actual_endpoint_RHS=actual,pilot_cumulative_RHS=book.pilot_rhs,chain_cumulative_RHS=book.pilot_rhs-chain_start,
                            selected_density_flags=flags,GT_used=False,posterior_qualification=False)
                        stream.write(json.dumps(finite(row),allow_nan=False,separators=(',',':'))+'\n');stream.flush()
                        qs.append(current['q'].copy());phis.append(current['Phi']);priors.append(current['prior_energy']);logtargets.append(current['logtarget'])
                        accepted.append(yes);blocks.append(block);keys.append(current['key']);materials.append(current['chi'].copy());RMS.append(material_RMS);ratios.append(ratio)
                        if (transition+1)%64==0 or transition+1==args.warmup+args.production:
                            chaindir=armdir/f'chain_{chain}';chaindir.mkdir(exist_ok=True)
                            save(chaindir/'CHAIN.npz',q=np.asarray(qs),Phi=np.asarray(phis),prior_energy=np.asarray(priors),logtarget=np.asarray(logtargets),
                                accepted=np.asarray(accepted),block=np.asarray(blocks),endpoint_key=np.asarray(keys),proposal_material_RMS=np.asarray(RMS),log_mh_ratio=np.asarray(ratios),
                                actual_RHS_per_transition=np.full(len(qs),8),warmup=args.warmup)
                            write(chaindir/'CHECKPOINT.json',dict(transition=transition+1,phase=phase,actual_RHS=book.pilot_rhs-chain_start,
                                RNG=rng.bit_generator.state,all_rejections_and_block_states_retained=True,parameters_frozen=True,posterior_qualification=False))
                            print(json.dumps(dict(parent=parent,arm=arm,chain=chain,transition=transition+1,actual_chain_RHS=book.pilot_rhs-chain_start,
                                Phi=current['Phi'],acceptance_last64=float(np.mean(accepted[-64:])))),flush=True)
                    endpoint_indices=np.arange(args.warmup+1,len(materials),2) if arm=='learned_blocked' else np.arange(args.warmup,len(materials))
                    chi_chains.append(np.asarray(materials)[endpoint_indices]);blockstats={}
                    for block in sorted(set(blocks)):
                        mask=np.array([i>=args.warmup and b==block for i,b in enumerate(blocks)])
                        blockstats[block]=dict(production_moves=int(mask.sum()),production_acceptance=float(np.mean(np.asarray(accepted)[mask])),
                            proposal_material_RMS_median=float(np.median(np.asarray(RMS)[mask])),accepted_material_RMS_mean=float(np.mean(np.asarray(RMS)[mask]*np.asarray(accepted)[mask])))
                    chain_rows.append(dict(chain=chain,actual_RHS=book.pilot_rhs-chain_start,warmup_RHS=args.warmup*8,production_RHS=args.production*8,
                        production_transitions=args.production,retained_production_sweep_endpoints=len(endpoint_indices),blocks=blockstats))
                chi=np.asarray(chi_chains);save(armdir/'PRODUCTION_MATERIAL.npz',chi=chi)
                diagnostics,summary=summarize_material(chi);save(armdir/'MATERIAL_DIAGNOSTICS.npz',**{k:v for k,v in diagnostics.items() if isinstance(v,np.ndarray)})
                result=dict(parent=parent,arm=arm,actual_RHS=book.pilot_rhs-arm_start,chains=chain_rows,material_functions=summary,
                    mean_material_RMS_from_initialization=float(np.linalg.norm(chi.mean(axis=(0,1))-initial['chi'])/24),
                    selected_density_flags=flags,proper_full66_proposal=arm!='learned_blocked',
                    proper_native_active10_density_with_prior56_reversible_block=arm=='learned_blocked',
                    reference_qualification=False,posterior_qualification=False,claim='Exact correction mechanism with finite-run diagnostics only')
                write(armdir/'RESULT.json',finite(result));results.append(result)
        summary=dict(protocol=protocol,status='D2_LEARNED_CORRECTION_COMPARISON_COMPLETE_UNQUALIFIED',actual_new_RHS=book.pilot_rhs,
            cost_receipt_by_parent_arm_chain_phase=book.costs,results=results,elapsed_seconds=time.perf_counter()-start,
            posterior_qualification=False,reference_qualification=False,
            interpretation='Learned10 times independent prior56 is a proposal. Full likelihood corrections establish invariant kernels, not finite-run posterior calibration.')
        write(dest/'RESULT.json',finite(summary));print(json.dumps(finite(summary)),flush=True)
    except BaseException as exc:
        status='ERROR';error=repr(exc);write(book.jobdir/'ERROR.json',dict(error=error,actual_new_RHS=book.pilot_rhs,costs=book.costs,all_completed_endpoints_and_transition_receipts_preserved=True));raise
    finally:
        if stream is not None: stream.close()
        book.cost_owner='closure';book.close(status,error)


if __name__=='__main__': main()
