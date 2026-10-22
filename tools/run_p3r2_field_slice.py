"""Conditional Gaussian-field2 SLICE: FullMaxwell GH versus exact pCN-MH.

The frozen two-dimensional span follows posterior material risk at a paid
anchor. Absolute Gaussian2 coordinates vary; the orthogonal Gaussian62 and
both native properties remain fixed. This is never a FULL marginal/reference.
No physics runs at import; the parent explicitly queues a registered job.
"""
from pathlib import Path
import argparse,gc,json,os,subprocess,sys,time
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(name,'1')
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.special import logsumexp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'tools'))
from p3r2.core import OUT,Book,Target,BudgetStop,write,save
from p3l.diagnostics import material_diagnostics
from run_p3r2_posterior import finite,freeze
LOG2PI=float(np.log(2*np.pi))
DEFAULT_ORDERS=(9,17,33,65,97)


class SliceBook(Book):
    """Guard each real Book reservation before the corresponding physical call."""
    def __init__(self,*args,pilot_rhs_cap=196608,**kwargs):
        self.pilot_rhs_cap=int(pilot_rhs_cap);self.pilot_rhs=0;self.cost_owner='setup';self.costs={}
        super().__init__(*args,**kwargs)
    def _update(self,event,**fields):
        rhs=int(fields.get('rhs',0))
        if self.pilot_rhs+rhs>self.pilot_rhs_cap: raise BudgetStop('field-SLICE independent pilot RHS cap')
        result=super()._update(event,pilot_rhs_before=self.pilot_rhs,pilot_rhs_cap=self.pilot_rhs_cap,**fields)
        self.pilot_rhs+=rhs;self.costs[self.cost_owner]=self.costs.get(self.cost_owner,0)+rhs
        return result


def risk_span(Jz,Dz):
    """Array-only leading2 C^.5 M C^.5 risk span, then Euclidean QR."""
    Jz=np.asarray(Jz,float);Dz=np.asarray(Dz)
    if Jz.shape[1]!=64 or Dz.shape!=(576,64): raise ValueError('frozen P1 Jz and material Dz layout required')
    H=np.eye(64)+Jz.T@Jz;H=(H+H.T)/2
    values,vectors=np.linalg.eigh(H)
    if values.min()<=0: raise FloatingPointError('prior+GN precision is not SPD')
    Croot=(vectors*(1/np.sqrt(values)))@vectors.T
    M=np.real(Dz.conj().T@Dz)/576;M=(M+M.T)/2
    risk=Croot@M@Croot;risk=(risk+risk.T)/2
    risk_values,risk_vectors=np.linalg.eigh(risk)
    mapped=Croot@risk_vectors[:,-2:][:,::-1];V,R=np.linalg.qr(mapped,mode='reduced')
    # Freeze a reproducible sign gauge; the selected orthogonal span is the
    # scientific object and a sign change only reverses a displayed CDF axis.
    for axis in range(2):
        pivot=int(np.argmax(abs(V[:,axis])))
        if V[pivot,axis]<0: V[:,axis]*=-1
    orth_error=float(np.linalg.norm(V.T@V-np.eye(2)))
    if orth_error>1e-11: raise FloatingPointError('field2 span failed orthogonality')
    return V,dict(H=H,Croot=Croot,material_Gram=M,risk=risk,risk_eigenvalues_descending=risk_values[::-1],
                  leading_mapped_span=mapped,orthogonality_error=orth_error,
                  selection_boundary_gap=float(risk_values[-2]-risk_values[-3]),
                  leading_selection_gap=float(risk_values[-1]-risk_values[-2]))


def gh_rule(order):
    nodes,weights=hermgauss(order);axis=np.sqrt(2)*nodes;one_logweights=np.log(weights)-.5*np.log(np.pi)
    a,b=np.meshgrid(axis,axis,indexing='ij');wa,wb=np.meshgrid(one_logweights,one_logweights,indexing='ij')
    return np.column_stack((a.ravel(),b.ravel())),(wa+wb).ravel(),axis


def empirical_marginal_cdf(coordinates,weights,grid):
    coordinates=np.asarray(coordinates);weights=np.asarray(weights,float);weights/=weights.sum()
    return np.array([[weights[coordinates[:,axis]<=threshold].sum() for threshold in grid] for axis in range(2)])


def stratum(prior,q):
    _,decoded=prior.decode(q);field=decoded['field']
    value=np.where(field<=-prior.temperature,-1,np.where(field>=prior.temperature,1,0)).astype(np.int8)
    value[~prior.mask]=-1  # Outside the ROI the decoder is exactly host.
    return value,decoded['fraction']


def quadrature_moments(prior,q0,rest,V,a,Phi,chi,strata,log_prior_weights,cdf_grid):
    logZ=float(logsumexp(-Phi+log_prior_weights));weights=np.exp(-Phi+log_prior_weights-logZ)
    mean=weights@chi;variance=np.einsum('i,ij->j',weights,np.abs(chi-mean)**2)
    a_mean=weights@a;mean_q=q0.copy();mean_q[:64]=rest+V@a_mean
    decoder_of_mean,decoded=prior.decode(mean_q);gap=mean-decoder_of_mean
    delta=decoded['inclusion_increment'];fraction_mean=np.real((mean-prior.host)/delta)
    centered=a-a_mean;cov=centered.T@(weights[:,None]*centered)
    cdf=empirical_marginal_cdf(a,weights,cdf_grid)
    return dict(weights=weights,mean=mean,variance=variance,coordinate_mean=a_mean,coordinate_covariance=cov,
        log_integral=logZ,cdf_grid=cdf_grid,quadrature_empirical_marginal_CDF=cdf,
        conditional_latent_mean_q=mean_q,decoder_at_conditional_latent_mean=decoder_of_mean,
        decoder_Jensen_gap=gap,decoder_Jensen_gap_material_RMS=float(np.linalg.norm(gap)/24),
        host_plateau_probability_per_cell=weights@(strata==-1),pure_plateau_probability_per_cell=weights@(strata==1),
        transition_probability_per_cell=weights@(strata==0),
        decoder_fraction_at_latent_mean=decoded['fraction'],conditional_decoder_fraction_mean=fraction_mean,
        decoder_fraction_Jensen_gap=fraction_mean-decoded['fraction'])


def array_controls():
    # A polynomial rule check catches the sqrt2 and1/pi Gaussian-measure
    # factors; pCN balance is checked independently from the likelihood ratio.
    a,lw,_=gh_rule(9);w=np.exp(lw);rng=np.random.default_rng(2026102396)
    ratio_errors=[];beta=.2;rho=np.sqrt(1-beta**2)
    for _ in range(100):
        x=rng.normal(size=2);y=rho*x+beta*rng.normal(size=2)
        forward=-np.log(2*np.pi*beta**2)-np.sum((y-rho*x)**2)/(2*beta**2)
        reverse=-np.log(2*np.pi*beta**2)-np.sum((x-rho*y)**2)/(2*beta**2)
        ratio_errors.append(abs(reverse-forward-.5*(y@y-x@x)))
    result=dict(status='PURE_FIELD_SLICE_MEASURE_CONTROLS_COMPLETE',new_physics_RHS=0,
        standard_normal_prior_integral=float(w.sum()),standard_normal_coordinate_mean=w@a,
        standard_normal_second_moments=w@(a*a),pCN_prior_ratio_max_abs_error=float(max(ratio_errors)),
        checks=dict(prior_mass=abs(w.sum()-1)<1e-12,prior_mean=np.max(abs(w@a))<1e-12,
                    prior_second_moments=np.max(abs(w@(a*a)-1))<1e-12,pCN_ratio=max(ratio_errors)<1e-11),
        scope='Pure Gaussian measure and reversible proposal controls; no FullMaxwell or posterior qualification')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job');parser.add_argument('--parent',default='DEV_030',choices=['DEV_030'])
    parser.add_argument('--anchor',type=Path,default=ROOT.parent/'b1-p3l-local-posterior/research/p3l_local_posterior_imaging/A/anchors/DEV_030.npz')
    parser.add_argument('--orders',default=','.join(map(str,DEFAULT_ORDERS)));parser.add_argument('--refinement-tolerance',type=float,default=1e-6)
    parser.add_argument('--chains',type=int,default=4);parser.add_argument('--warmup',type=int,default=256);parser.add_argument('--production',type=int,default=1024)
    parser.add_argument('--beta',type=float,default=.2);parser.add_argument('--rhs-cap',type=int,default=196608)
    parser.add_argument('--max-wall',type=int,default=3600);parser.add_argument('--seed',type=int,default=2026102397)
    parser.add_argument('--verify-only',action='store_true');parser.add_argument('--validation-output',type=Path,default=OUT/'B/field_slice_array_validation/RESULT.json')
    args=parser.parse_args()
    if args.verify_only:
        result=array_controls();freeze(args.validation_output,result);print(json.dumps(finite(result),indent=2));return
    if not args.job: parser.error('--job required for a paid pilot')
    orders=tuple(int(x) for x in args.orders.split(','))
    if any(x<3 for x in orders) or any(b<=a for a,b in zip(orders,orders[1:])): raise ValueError('strictly increasing GH orders>=3 required')
    if args.rhs_cap>196608 or args.max_wall>3600 or args.chains<2 or args.production<4 or args.warmup<0 or not 0<args.beta<=1:
        raise ValueError('bounded declared field-SLICE pilot required')
    maximum_states=sum(n*n for n in orders)+args.chains*(1+args.warmup+args.production)
    if maximum_states*8>args.rhs_cap: raise ValueError('maximum frozen recipe exceeds independent RHS cap')
    with np.load(args.anchor,allow_pickle=False) as z:
        if 'receiver_T' not in z.files or not z['receiver_T'].size: raise ValueError('already-paid receiver T required; no unregistered J solves')
    book=SliceBook('B',args.job,max_wall=args.max_wall,pilot_rhs_cap=args.rhs_cap)
    status='COMPLETE';error=None;stream=None;start=time.perf_counter();rows=[]
    try:
        dest=book.jobdir
        git_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
        own_status=subprocess.check_output(['git','status','--porcelain','--','tools/run_p3r2_field_slice.py'],cwd=ROOT,text=True).strip()
        source_files=['tools/run_p3r2_field_slice.py','src/p3r2/core.py','src/p3l/diagnostics.py','tools/run_p3r2_posterior.py','src/p1_material/prior.py']
        for name in source_files:
            copied=dest/'source'/name;copied.parent.mkdir(parents=True,exist_ok=True);copied.write_bytes((ROOT/name).read_bytes())
        write(dest/'SOURCE_PROVENANCE.json',dict(git_revision=git_revision,executed_script_worktree_status=own_status,source_files=source_files,
            exact_executed_script_copy=str(dest/'source/tools/run_p3r2_field_slice.py'),
            source_of_truth='Exact copied script bytes; a dirty script is not claimed identical to Git HEAD',
            python=sys.version,numpy=np.__version__,no_hash_check=True))
        target=Target(book);y=target.observation(args.parent)
        anchor=target.load_cache(args.anchor,owner='Historical P3L paid A DEV030 anchor; original observation, old ledger immutable')
        q0=anchor.q.copy();Jz=target.jacobian(anchor)[:,:64].copy();Dz=target.prior.derivative(q0)[:,:64]
        V,geometry=risk_span(Jz,Dz);a0=V.T@q0[:64];rest=q0[:64]-V@a0
        if np.linalg.norm(V.T@rest)>1e-11 or np.linalg.norm(rest+V@a0-q0[:64])>1e-11: raise RuntimeError('absolute orthogonal conditional coordinates')
        cdf_grid=np.linspace(-6.,6.,241)
        protocol=dict(job=args.job,parent=args.parent,stage='B',target_id=target.cfg['target_id'],
            git_revision=git_revision,executed_script_worktree_status=own_status,executed_source_snapshot='source/tools/run_p3r2_field_slice.py',
            target_family='SLICE_GAUSSIAN_FIELD2_FIXED_ORTHOGONAL62_AND_PROPERTIES',original_FULL_posterior=False,original_P1_target_unchanged=True,
            anchor=str(args.anchor.resolve()),q0=q0,V=V,rest=rest,a0=a0,
            span='C=(I+Jz.TJz)^-1; leading2 eigvectors of C^.5 Re(Dz.H Dz/576) C^.5; Euclidean QR(C^.5 eigvectors)',
            material_weight='Complex full-grid material RMS: Re(Dz.H Dz)/576',
            coordinates='Absolute a=V.Tz ~ normalized N(0,I2), rest=z0-VV.Tz0 fixed; z=rest+Va; native properties fixed at q0.u',
            density='exp(-Phi(FullMaxwell(rest+Va,q0.u),y))*standardNormal2(a); no Gram determinant or local shifted prior',
            likelihood_normalization='Constants independent of a omitted; logZ is the conditional Gaussian2 integral of exp(-Phi)',
            quadrature='Product Gauss-Hermite for normalized standard Gaussian2: a=sqrt2*hermite_nodes, prior weights=hermite_wi*hermite_wj/pi',
            orders=orders,stopping_rule='Two consecutive material-mean RMS and absolute logZ differences <= frozen tolerance',
            tolerance=args.refinement_tolerance,CDF_grid=cdf_grid,
            CDF_scope='Weighted empirical GH-node marginal CDFs at fixed grid, refinement diagnostic only; no continuum certification or credible interval claims',
            chains=args.chains,warmup=args.warmup,production=args.production,pcn_beta=args.beta,
            chain_initialization='Independently seeded posterior-weighted final-GH node draws; each initial state freshly reevaluated and billed8 RHS; no GT',
            pCN='a_new=sqrt(1-beta²)*a+beta*normal2; original conditional Gaussian2 prior cancels, accept min(1,exp(-deltaPhi))',
            Jensen_gap='Conditional E[decoder(q)]-decoder(conditional E[q]), exact decoder array arithmetic; not a FULL posterior moment',
            source_isolated_online_observation=True,GT_used=False,training_arrays_read=False,
            maximum_new_states=maximum_states,maximum_new_RHS=maximum_states*8,pilot_RHS_cap=args.rhs_cap,max_wall_seconds=args.max_wall,
            memory_policy='One physical Point at a time; raw arrays capped at400MiB per GH order; no retained historical LU',
            reference_qualification=False,posterior_qualification=False,seed=args.seed)
        freeze(dest/'PROTOCOL.json',protocol)
        save(dest/'FRAME.npz',q0=q0,V=V,rest=rest,a0=a0,Jz=Jz,Dz=Dz,**geometry)
        controls=array_controls();write(dest/'GAUSSIAN_MEASURE_CONTROLS.json',controls)
        if not all(controls['checks'].values()): raise RuntimeError('Gaussian2 measure control failed')
        del anchor,Jz,Dz
        if target.model is not None: target.model._transfers.clear()
        gc.collect()
        def compose(a):
            q=q0.copy();q[:64]=rest+V@np.asarray(a,float);return q
        def evaluate(a,purpose):
            if book.pilot_rhs+8>args.rhs_cap: raise BudgetStop('stop before next Full solve due field-SLICE pilot cap')
            q=compose(a);before=book.pilot_rhs;p=target.evaluate(q,purpose=purpose)
            value=dict(a=np.asarray(a,float).copy(),q=q,Phi=target.phi(p,y),chi=p.state.chi.copy(),key=p.key)
            value['stratum'],value['fraction']=stratum(target.prior,q)
            if not np.isfinite(value['Phi']): raise FloatingPointError('nonfinite exact field-SLICE Phi')
            del p
            if target.model is not None: target.model._transfers.clear()
            if book.pilot_rhs-before!=8: raise RuntimeError('field-SLICE exact Full endpoint RHS differs from8')
            return value
        previous_mean=None;previous_logZ=None;previous_CDF=None;stable=0;reference=None
        for order in orders:
            book.cost_owner=f'quadrature/order_{order}';a,logw,axis=gh_rule(order);n=len(a)
            estimated_bytes=n*(66*8+576*16+576+8*8)
            if estimated_bytes>400*1024**2: raise MemoryError('GH raw arrays exceed frozen400MiB policy')
            Phi=np.empty(n);chi=np.empty((n,576),complex);strata=np.empty((n,576),np.int8);keys=np.empty(n,dtype='U96')
            for i,coordinate in enumerate(a):
                p=evaluate(coordinate,f'field2_GH_order{order}')
                Phi[i]=p['Phi'];chi[i]=p['chi'];strata[i]=p['stratum'];keys[i]=p['key'];del p
                if (i+1)%256==0 or i+1==n:
                    write(dest/'GH_NODE_PROGRESS.json',dict(order=order,completed_nodes=i+1,total_nodes=n,actual_RHS=book.pilot_rhs))
            moments=quadrature_moments(target.prior,q0,rest,V,a,Phi,chi,strata,logw,cdf_grid)
            save(dest/f'GH_{order}_RAW.npz',a=a,Phi=Phi,chi=chi,stratum=strata,endpoint_key=keys,
                log_standard_normal_quadrature_weights=logw,normalized_conditional_posterior_weights=moments['weights'],axis=axis)
            save(dest/f'GH_{order}_MOMENTS.npz',**moments)
            rms=None if previous_mean is None else float(np.linalg.norm(moments['mean']-previous_mean)/24)
            logdiff=None if previous_logZ is None else float(abs(moments['log_integral']-previous_logZ))
            cdfdiff=None if previous_CDF is None else np.max(abs(moments['quadrature_empirical_marginal_CDF']-previous_CDF),axis=1)
            stable=stable+1 if previous_mean is not None and rms<=args.refinement_tolerance and logdiff<=args.refinement_tolerance else 0
            row=dict(order=order,new_Full_states=n,new_RHS=n*8,logZ=moments['log_integral'],
                material_mean_successive_RMS=rms,absolute_logZ_successive_difference=logdiff,
                empirical_marginal_CDF_successive_max_difference=cdfdiff,stable_refinements=stable,
                decoder_Jensen_gap_material_RMS=moments['decoder_Jensen_gap_material_RMS'],
                maximum_normalized_node_weight=float(moments['weights'].max()),reference_qualification=False)
            rows.append(row);write(dest/'QUADRATURE_PROGRESS.json',finite(dict(rows=rows,stopping_rule_met=stable>=2,actual_RHS=book.pilot_rhs,
                scope='Conditional field2 numerical refinement; empirical node CDF, not continuum/FULL certification')))
            print(json.dumps(finite(row)),flush=True)
            previous_mean=moments['mean'].copy();previous_logZ=moments['log_integral'];previous_CDF=moments['quadrature_empirical_marginal_CDF'].copy()
            # Only final nodes/Phi/weights and small moments survive to the
            # next order. Material raw arrays are immediately persisted.
            reference=dict(a=a.copy(),Phi=Phi.copy(),weights=moments['weights'].copy(),moments=moments,order=order)
            del chi,strata,Phi,keys,logw;gc.collect()
            if stable>=2: break
        write(dest/'CHAIN_PROTOCOL.json',dict(frozen_before_first_chain=True,final_GH_order=reference['order'],
            initialization='Independent posterior-weighted GH-node indices; starting points not exact posterior draws',
            every_initial_reevaluation_and_proposal_billed_RHS=8,chains=args.chains,warmup=args.warmup,production=args.production,
            beta=args.beta,full_target=False,conditional_slice_only=True))
        stream=(dest/'CHAIN_TRANSITIONS.jsonl').open('x');rho=np.sqrt(1-args.beta**2)
        chain_a=[];chain_chi=[];chain_Phi=[];chain_accept=[];chain_rows=[]
        for chain in range(args.chains):
            rng=np.random.default_rng(np.random.SeedSequence([args.seed,chain]));node=int(rng.choice(len(reference['a']),p=reference['weights']))
            book.cost_owner=f'chain_{chain}/initialization';current=evaluate(reference['a'][node],f'field2_pCN_chain{chain}_initial_node_reevaluation')
            chaindir=dest/'chains'/f'chain_{chain}';chaindir.mkdir(parents=True)
            initial_error=float(current['Phi']-reference['Phi'][node]);write(chaindir/'INITIALIZATION.json',dict(GH_order=reference['order'],GH_node_index=node,
                initial_a=current['a'],GH_normalized_weight=reference['weights'][node],fresh_endpoint=current['key'],fresh_RHS=8,
                paid_GH_Phi=reference['Phi'][node],reevaluated_Phi=current['Phi'],reevaluation_difference=initial_error,
                initialization_is_discrete_node_draw_not_continuum_posterior=True,GT_used=False))
            coordinates=[];materials=[];phis=[];accepted=[];proposed_RMS=[];endpoint_keys=[];log_ratios=[]
            for iteration in range(args.warmup+args.production):
                phase='warmup' if iteration<args.warmup else 'production';book.cost_owner=f'chain_{chain}/{phase}'
                old=current;proposal_a=rho*old['a']+args.beta*rng.normal(size=2)
                proposed=evaluate(proposal_a,f'field2_pCN_chain{chain}_{phase}')
                ratio=-proposed['Phi']+old['Phi'];logu=float(np.log(max(rng.random(),np.finfo(float).tiny)));yes=bool(logu<min(0.,ratio))
                forward=-np.log(2*np.pi*args.beta**2)-float(np.sum((proposal_a-rho*old['a'])**2))/(2*args.beta**2)
                reverse=-np.log(2*np.pi*args.beta**2)-float(np.sum((old['a']-rho*proposal_a)**2))/(2*args.beta**2)
                prior_delta=.5*float(proposal_a@proposal_a-old['a']@old['a'])
                current=proposed if yes else old;material_RMS=float(np.linalg.norm(proposed['chi']-old['chi'])/24)
                row=dict(chain=chain,iteration=iteration,phase=phase,accepted=yes,repeated_state=not yes,
                    current_a=old['a'],proposed_a=proposal_a,retained_a=current['a'],
                    current_Phi=old['Phi'],proposed_Phi=proposed['Phi'],retained_Phi=current['Phi'],
                    current_logtarget=-old['Phi']-.5*float(old['a']@old['a'])-LOG2PI,
                    proposed_logtarget=-proposed['Phi']-.5*float(proposal_a@proposal_a)-LOG2PI,
                    log_proposal_forward=forward,log_proposal_reverse=reverse,log_proposal_correction=reverse-forward,
                    prior_cancellation_error=abs(reverse-forward-prior_delta),log_mh_ratio=ratio,log_uniform=logu,
                    proposed_material_RMS=material_RMS,accepted_material_RMS=material_RMS if yes else 0.,
                    decoder_stratum_changed_ROI_cells=int(np.sum((proposed['stratum']!=old['stratum'])&target.prior.mask)),
                    current_endpoint=old['key'],proposed_endpoint=proposed['key'],retained_endpoint=current['key'],
                    actual_endpoint_RHS=8,pilot_cumulative_RHS=book.pilot_rhs,target_family='SLICE_GAUSSIAN_FIELD2',GT_used=False)
                stream.write(json.dumps(finite(row),allow_nan=False,separators=(',',':'))+'\n');stream.flush()
                coordinates.append(current['a'].copy());materials.append(current['chi'].copy());phis.append(current['Phi']);accepted.append(yes)
                proposed_RMS.append(material_RMS);endpoint_keys.append(current['key']);log_ratios.append(ratio)
                if (iteration+1)%128==0 or iteration+1==args.warmup+args.production:
                    save(chaindir/'CHAIN.npz',a=np.asarray(coordinates),chi=np.asarray(materials),Phi=np.asarray(phis),accepted=np.asarray(accepted),
                        proposed_material_RMS=np.asarray(proposed_RMS),endpoint_key=np.asarray(endpoint_keys),log_mh_ratio=np.asarray(log_ratios),warmup=args.warmup)
                    write(chaindir/'CHECKPOINT.json',dict(iteration=iteration+1,RNG=rng.bit_generator.state,beta=args.beta,parameters_frozen=True,
                        current_Phi=current['Phi'],actual_pilot_RHS=book.pilot_rhs,all_rejections_retained=True))
                    print(json.dumps(dict(chain=chain,iteration=iteration+1,phase=phase,Phi=current['Phi'],
                        acceptance_last128=float(np.mean(accepted[-128:])),pilot_RHS=book.pilot_rhs)),flush=True)
                del old,proposed
            chain_a.append(np.asarray(coordinates[args.warmup:]));chain_chi.append(np.asarray(materials[args.warmup:]));chain_Phi.append(np.asarray(phis[args.warmup:]));chain_accept.append(np.asarray(accepted[args.warmup:]))
            chain_rows.append(dict(chain=chain,initial_GH_node_index=node,production_acceptance=float(np.mean(accepted[args.warmup:])),
                warmup_RHS=args.warmup*8,production_RHS=args.production*8,initial_reevaluation_RHS=8))
        a=np.asarray(chain_a);chi=np.asarray(chain_chi);Phi=np.asarray(chain_Phi);accept=np.asarray(chain_accept)
        save(dest/'PRODUCTION_CHAINS.npz',a=a,chi=chi,Phi=Phi,accepted=accept)
        diagnostics=material_diagnostics(chi);coordinate_diagnostics=material_diagnostics(a)
        save(dest/'MATERIAL_DIAGNOSTICS.npz',**{k:v for k,v in diagnostics.items() if isinstance(v,np.ndarray)})
        save(dest/'COORDINATE_DIAGNOSTICS.npz',**{k:v for k,v in coordinate_diagnostics.items() if isinstance(v,np.ndarray)})
        mean=chi.mean(axis=(0,1));variance=np.mean(abs(chi-mean)**2,axis=(0,1));amean=a.mean(axis=(0,1))
        meanq=compose(amean);decoder_mean=target.prior.decode(meanq)[0];gap=mean-decoder_mean
        cdf=empirical_marginal_cdf(a.reshape(-1,2),np.ones(a.shape[0]*a.shape[1]),cdf_grid)
        save(dest/'CHAIN_MOMENTS_AND_JENSEN.npz',mean=mean,variance=variance,coordinate_mean=amean,
            empirical_marginal_CDF=cdf,cdf_grid=cdf_grid,conditional_latent_mean_q=meanq,
            decoder_at_conditional_latent_mean=decoder_mean,decoder_Jensen_gap=gap)
        ref=reference['moments'];ess=diagnostics['mean_ess'];rhat=diagnostics['rhat'];mcse=diagnostics['mean_mcse']
        fe=ess[np.isfinite(ess)];fr=rhat[np.isfinite(rhat)]
        result=dict(protocol=protocol,status='FIELD2_CONDITIONAL_SLICE_PILOT_COMPLETE_UNQUALIFIED',actual_new_RHS=book.pilot_rhs,
            cost_receipt=book.costs,quadrature_rows=rows,quadrature_stopping_rule_met=stable>=2,final_GH_order=reference['order'],
            chains=chain_rows,coordinate_mean=amean,GH_coordinate_mean=ref['coordinate_mean'],coordinate_mean_error=amean-ref['coordinate_mean'],
            coordinate_mean_ESS=coordinate_diagnostics['mean_ess'],coordinate_rhat=coordinate_diagnostics['rhat'],
            empirical_marginal_CDF_max_difference_chain_vs_GH=np.max(abs(cdf-ref['quadrature_empirical_marginal_CDF']),axis=1),
            material_mean_RMS_error_chain_vs_GH=float(np.linalg.norm(mean-ref['mean'])/24),
            material_variance_RMS_error_chain_vs_GH=float(np.linalg.norm(variance-ref['variance'])/24),
            material_mean_ESS_minimum_finite=float(fe.min()) if fe.size else None,
            material_rhat_maximum_finite=float(fr.max()) if fr.size else None,
            undefined_material_ESS_functions=int((~np.isfinite(ess)).sum()),undefined_material_rhat_functions=int((~np.isfinite(rhat)).sum()),
            complex_material_mean_RMS_MCSE=float(np.sqrt(np.mean(np.sum(mcse**2,axis=-1)))),
            GH_decoder_Jensen_gap_material_RMS=ref['decoder_Jensen_gap_material_RMS'],
            chain_decoder_Jensen_gap_material_RMS=float(np.linalg.norm(gap)/24),
            Jensen_gap_RMS_difference_chain_vs_GH=float(np.linalg.norm(gap-ref['decoder_Jensen_gap'])/24),
            elapsed_seconds=time.perf_counter()-start,original_FULL_posterior=False,posterior_qualification=False,reference_qualification=False,
            interpretation='Conditional field2 nonlinearity/saturation mechanism only; GH refinement and pCN agreement do not establish a continuum reference or FULL posterior calibration')
        write(dest/'RESULT.json',finite(result));print(json.dumps(finite(result)),flush=True)
    except BaseException as exc:
        status='ERROR';error=repr(exc);write(book.jobdir/'ERROR.json',dict(error=error,actual_new_RHS=book.pilot_rhs,cost_receipt=book.costs,
            completed_quadrature_rows=rows,all_completed_endpoints_immediately_persisted=True,no_automatic_replay=True));raise
    finally:
        if stream is not None: stream.close()
        book.cost_owner='closure';book.close(status,error)


if __name__=='__main__': main()
