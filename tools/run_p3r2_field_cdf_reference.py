"""Same-target field2 reference: Gaussian CDF coordinates and patch-integrated CDFs.

Absolute a~N(0,I2) is transformed to t=NormalCDF(a) on the entire unit square.
The Gaussian prior cancels the change-of-variable Jacobian exactly. Actual
FullMaxwell likelihoods, never a surrogate, are integrated by partitioned GL.
This is a conditional field2 numerical reference, never a FULL posterior.
"""
from pathlib import Path
import argparse,gc,json,os,subprocess,sys,time
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(name,'1')
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import ndtr,ndtri,logsumexp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'tools'))
from p3r2.core import OUT,Book,Target,BudgetStop,write,save
from run_p3r2_posterior import finite,freeze


class ReferenceBook(Book):
    def __init__(self,*args,pilot_rhs_cap=500000,**kwargs):
        self.pilot_rhs_cap=int(pilot_rhs_cap);self.pilot_rhs=0;self.cost_owner='setup';self.costs={}
        super().__init__(*args,**kwargs)
    def _update(self,event,**fields):
        rhs=int(fields.get('rhs',0))
        if self.pilot_rhs+rhs>self.pilot_rhs_cap: raise BudgetStop('independent field-CDF reference RHS cap')
        result=super()._update(event,pilot_rhs_before=self.pilot_rhs,pilot_rhs_cap=self.pilot_rhs_cap,**fields)
        self.pilot_rhs+=rhs;self.costs[self.cost_owner]=self.costs.get(self.cost_owner,0)+rhs
        return result


def freeze_cuts(mean,covariance,a0,maximum_node):
    sd=np.sqrt(np.diag(covariance));cuts=[]
    if not np.isfinite(sd).all() or np.any(sd<=0): raise ValueError('positive finite GH coordinate SD required')
    for axis in range(2):
        values=[0.,1.]
        values.extend(ndtr(mean[axis]+np.array([-2.,-1.,0.,1.,2.])*sd[axis]).tolist())
        values.extend([float(ndtr(a0[axis])),float(ndtr(maximum_node[axis]))])
        cut=np.unique(np.asarray(values,float))
        if len(cut)>9 or cut[0]!=0 or cut[-1]!=1 or np.any(np.diff(cut)<=0): raise ValueError('entire unit axis with at most9 cuts required')
        cuts.append(cut)
    return cuts,sd


def patch_rule(order,bounds):
    nodes,weights=leggauss(order);axes=[];logweights=[]
    for low,high in bounds:
        if not low<high: raise ValueError('positive patch width required')
        axis=(low+high)/2+(high-low)/2*nodes
        if np.any(axis<=0) or np.any(axis>=1): raise FloatingPointError('patch nodes not representable in open unit square; no clipping')
        axes.append(axis);logweights.append(np.log(weights)+np.log((high-low)/2))
    a,b=np.meshgrid(*axes,indexing='ij');wa,wb=np.meshgrid(*logweights,indexing='ij')
    t=np.column_stack((a.ravel(),b.ravel()))
    return t,ndtri(t),(wa+wb).ravel()


def summarize_patch(a,Phi,chi,strata,logweights):
    logmass=float(logsumexp(-Phi+logweights));w=np.exp(-Phi+logweights-logmass)
    mean=w@chi;amean=w@a;centered=a-amean
    return dict(logmass=logmass,mean=mean,variance=np.einsum('i,ij->j',w,abs(chi-mean)**2),
        coordinate_mean=amean,coordinate_covariance=centered.T@(w[:,None]*centered),
        host_plateau_probability_per_cell=w@(strata==-1),pure_plateau_probability_per_cell=w@(strata==1),
        transition_probability_per_cell=w@(strata==0))


def combine_patches(patches,shape):
    lm=np.array([p['logmass'] for p in patches]);logZ=float(logsumexp(lm));mass=np.exp(lm-logZ)
    means=np.asarray([p['mean'] for p in patches]);mean=mass@means
    variance=np.einsum('i,ij->j',mass,np.asarray([p['variance'] for p in patches])+abs(means-mean)**2)
    am=np.asarray([p['coordinate_mean'] for p in patches]);amean=mass@am;acov=np.zeros((2,2))
    for weight,p in zip(mass,patches):
        delta=p['coordinate_mean']-amean;acov+=weight*(p['coordinate_covariance']+np.outer(delta,delta))
    grid=lm.reshape(shape)
    CDF0=np.array([np.exp(logsumexp(grid[:index,:])-logZ) for index in range(1,shape[0])])
    CDF1=np.array([np.exp(logsumexp(grid[:,:index])-logZ) for index in range(1,shape[1])])
    return dict(log_integral=logZ,patch_log_integrals=lm.reshape(shape),normalized_patch_masses=mass.reshape(shape),
        mean=mean,variance=variance,coordinate_mean=amean,coordinate_covariance=acov,
        directly_integrated_marginal_CDF_axis0=CDF0,directly_integrated_marginal_CDF_axis1=CDF1,
        host_plateau_probability_per_cell=mass@np.asarray([p['host_plateau_probability_per_cell'] for p in patches]),
        pure_plateau_probability_per_cell=mass@np.asarray([p['pure_plateau_probability_per_cell'] for p in patches]),
        transition_probability_per_cell=mass@np.asarray([p['transition_probability_per_cell'] for p in patches]))


def node_CDF(a,w,cuts):
    a=np.asarray(a);w=np.asarray(w,float);w=w/w.sum()
    return [np.array([w[a[:,axis]<=ndtri(c)].sum() for c in cuts[axis][1:-1]]) for axis in range(2)]


def pure_controls():
    cuts=[np.array([0.,.2,.7,1.]),np.array([0.,.4,.9,1.])];logmass=[];prior_mass=0.
    for x0,x1 in zip(cuts[0][:-1],cuts[0][1:]):
        for y0,y1 in zip(cuts[1][:-1],cuts[1][1:]):
            t,a,lw=patch_rule(8,((x0,x1),(y0,y1)))
            prior_mass+=float(np.exp(lw).sum());f=1+2*t[:,0]+3*t[:,1]
            logmass.append(float(logsumexp(np.log(f)+lw)))
    lm=np.array(logmass).reshape(3,3);logZ=float(logsumexp(lm));Z=np.exp(logZ)
    c0=np.array([np.exp(logsumexp(lm[:i,:])-logZ) for i in (1,2)])
    c1=np.array([np.exp(logsumexp(lm[:,:i])-logZ) for i in (1,2)])
    expected0=(2.5*cuts[0][1:-1]+cuts[0][1:-1]**2)/3.5
    expected1=(2*cuts[1][1:-1]+1.5*cuts[1][1:-1]**2)/3.5
    error=float(max(np.max(abs(c0-expected0)),np.max(abs(c1-expected1))))
    return dict(status='PURE_FULL_SQUARE_AND_DIRECT_CDF_CONTROLS_COMPLETE',new_physics_RHS=0,
        transformed_standard_Gaussian_prior_mass=prior_mass,polynomial_integral=Z,expected_polynomial_integral=3.5,
        direct_CDF_maximum_abs_error=error,
        checks=dict(full_prior_mass=abs(prior_mass-1)<1e-12,polynomial_integral=abs(Z-3.5)<1e-12,direct_patch_CDF=error<1e-12),
        scope='Array-only measure and whole-patch CDF control; no physical reference qualification')


def compare_previous(source,moments,cuts,q0,rest,V,prior):
    comparison=dict(source_job=str(source),scope='Same conditional field2 target; old files read only; no GT')
    with np.load(source/'GH_97_MOMENTS.npz',allow_pickle=False) as z: old={k:z[k].copy() for k in z.files}
    with np.load(source/'GH_97_RAW.npz',allow_pickle=False) as z:
        old_a=z['a'].copy();old_w=z['normalized_conditional_posterior_weights'].copy()
    oldcdf=node_CDF(old_a,old_w,cuts)
    comparison['GH97']=dict(logZ=float(old['log_integral']),new_minus_old_logZ=float(moments['log_integral']-old['log_integral']),
        material_mean_RMS_difference=float(np.linalg.norm(moments['mean']-old['mean'])/24),
        material_variance_RMS_difference=float(np.linalg.norm(moments['variance']-old['variance'])/24),
        coordinate_mean_difference=moments['coordinate_mean']-old['coordinate_mean'],
        GH_node_empirical_CDF_axis0=oldcdf[0],GH_node_empirical_CDF_axis1=oldcdf[1],
        max_abs_direct_integral_minus_GH_node_CDF=np.array([np.max(abs(moments['directly_integrated_marginal_CDF_axis0']-oldcdf[0])),np.max(abs(moments['directly_integrated_marginal_CDF_axis1']-oldcdf[1]))]),
        old_decoder_Jensen_gap_RMS=float(old['decoder_Jensen_gap_material_RMS']),
        Jensen_gap_material_RMS_difference=float(np.linalg.norm(moments['decoder_Jensen_gap']-old['decoder_Jensen_gap'])/24),
        caution='Old GH CDF is a weighted discrete node CDF; the new CDF is directly integrated at cut values')
    chain_path=source/'PRODUCTION_CHAINS.npz'
    if chain_path.exists():
        with np.load(chain_path,allow_pickle=False) as z: a=z['a'].copy();chi=z['chi'].copy();accepted=z['accepted'].copy()
        mean=chi.mean(axis=(0,1));variance=np.mean(abs(chi-mean)**2,axis=(0,1));amean=a.mean(axis=(0,1));meanq=q0.copy();meanq[:64]=rest+V@amean
        chain_gap=mean-prior.decode(meanq)[0];cdf=node_CDF(a.reshape(-1,2),np.ones(a.shape[0]*a.shape[1]),cuts)
        perchain=[]
        for index in range(len(a)):
            ccdf=node_CDF(a[index],np.ones(len(a[index])),cuts)
            perchain.append(dict(chain=index,material_mean_RMS_difference=float(np.linalg.norm(moments['mean']-chi[index].mean(axis=0))/24),
                CDF_axis0=ccdf[0],CDF_axis1=ccdf[1],
                max_abs_CDF_difference=np.array([np.max(abs(moments['directly_integrated_marginal_CDF_axis0']-ccdf[0])),np.max(abs(moments['directly_integrated_marginal_CDF_axis1']-ccdf[1]))])))
        comparison['existing_pCN_chains']=dict(chains=len(a),draws_per_chain=a.shape[1],acceptance_per_chain=accepted.mean(axis=1),
            material_mean_RMS_difference=float(np.linalg.norm(moments['mean']-mean)/24),
            material_variance_RMS_difference=float(np.linalg.norm(moments['variance']-variance)/24),
            coordinate_mean_difference=moments['coordinate_mean']-amean,
            pooled_empirical_CDF_axis0=cdf[0],pooled_empirical_CDF_axis1=cdf[1],
            max_abs_direct_CDF_minus_pooled_chain_ECDF=np.array([np.max(abs(moments['directly_integrated_marginal_CDF_axis0']-cdf[0])),np.max(abs(moments['directly_integrated_marginal_CDF_axis1']-cdf[1]))]),
            chain_decoder_Jensen_gap_material_RMS=float(np.linalg.norm(chain_gap)/24),
            Jensen_gap_material_RMS_difference=float(np.linalg.norm(moments['decoder_Jensen_gap']-chain_gap)/24),per_chain=perchain,
            caution='Existing pCN draws have finite mixing error; comparison does not certify either mechanism')
        diag_path=source/'MATERIAL_DIAGNOSTICS.npz'
        if diag_path.exists():
            with np.load(diag_path,allow_pickle=False) as z:
                mcse=z['mean_mcse'];ess=z['mean_ess'];rhat=z['rhat']
            comparison['existing_pCN_chains']['stored_complex_material_mean_RMS_MCSE']=float(np.sqrt(np.mean(np.sum(mcse**2,axis=-1))))
            comparison['existing_pCN_chains']['undefined_material_ESS_functions']=int((~np.isfinite(ess)).sum())
            comparison['existing_pCN_chains']['undefined_material_rhat_functions']=int((~np.isfinite(rhat)).sum())
    else: comparison['existing_pCN_chains']=dict(status='NOT_AVAILABLE_AT_COMPARISON_TIME',new_MCMC_run=False)
    return comparison


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job');parser.add_argument('--source-job',type=Path,default=OUT/'B/jobs/B_field2_slice_DEV030_r0')
    parser.add_argument('--orders',default='8,12,20');parser.add_argument('--rhs-cap',type=int,default=500000);parser.add_argument('--max-wall',type=int,default=3600)
    parser.add_argument('--mean-tolerance',type=float,default=1e-6);parser.add_argument('--logZ-tolerance',type=float,default=1e-6);parser.add_argument('--CDF-tolerance',type=float,default=.01)
    parser.add_argument('--verify-only',action='store_true');parser.add_argument('--validation-output',type=Path,default=OUT/'B/field_cdf_reference_array_validation/RESULT.json')
    args=parser.parse_args()
    if args.verify_only:
        result=pure_controls();freeze(args.validation_output,result);print(json.dumps(finite(result),indent=2));return
    if not args.job: parser.error('--job required for a paid reference')
    orders=tuple(int(x) for x in args.orders.split(','));source=args.source_job.resolve()
    if len(orders)<3 or any(n<2 for n in orders) or any(b<=a for a,b in zip(orders,orders[1:])): raise ValueError('three increasing GL orders required')
    if not 0<args.rhs_cap<=500000 or not 0<args.max_wall<=3600 or min(args.mean_tolerance,args.logZ_tolerance,args.CDF_tolerance)<=0: raise ValueError('bounded declared reference required')
    oldprotocol=json.loads((source/'PROTOCOL.json').read_text())
    if oldprotocol['parent']!='DEV_030' or oldprotocol['target_family']!='SLICE_GAUSSIAN_FIELD2_FIXED_ORTHOGONAL62_AND_PROPERTIES': raise ValueError('declared DEV030 field2 source required')
    with np.load(source/'FRAME.npz',allow_pickle=False) as z: q0=z['q0'].copy();V=z['V'].copy();rest=z['rest'].copy();a0=z['a0'].copy()
    with np.load(source/'GH_97_MOMENTS.npz',allow_pickle=False) as z:
        oldmean=z['coordinate_mean'].copy();oldcov=z['coordinate_covariance'].copy();oldweights=z['weights'].copy()
    with np.load(source/'GH_97_RAW.npz',allow_pickle=False) as z:
        maxindex=int(np.argmax(oldweights));maxnode=z['a'][maxindex].copy()
    if not np.array_equal(q0,np.asarray(oldprotocol['q0'])) or not np.array_equal(V,np.asarray(oldprotocol['V'])) or not np.array_equal(rest,np.asarray(oldprotocol['rest'])):
        raise ValueError('paid FRAME does not match the source frozen protocol')
    if np.linalg.norm(V.T@V-np.eye(2))>1e-10 or np.linalg.norm(V.T@rest)>1e-10 or np.linalg.norm(rest+V@a0-q0[:64])>1e-10: raise ValueError('frozen same-target orthogonal coordinates')
    cuts,sd=freeze_cuts(oldmean,oldcov,a0,maxnode);shape=(len(cuts[0])-1,len(cuts[1])-1)
    patches=[(i,j,((float(x0),float(x1)),(float(y0),float(y1)))) for i,(x0,x1) in enumerate(zip(cuts[0][:-1],cuts[0][1:])) for j,(y0,y1) in enumerate(zip(cuts[1][:-1],cuts[1][1:]))]
    maximum_states=len(patches)*sum(n*n for n in orders);maximum_rhs=maximum_states*8
    if maximum_rhs>args.rhs_cap: raise ValueError(f'worst-case patch recipe needs{maximum_rhs}RHS above cap{args.rhs_cap}')
    book=ReferenceBook('B',args.job,max_wall=args.max_wall,pilot_rhs_cap=args.rhs_cap)
    status='COMPLETE';error=None;events=None;rows=[];start=time.perf_counter()
    try:
        dest=book.jobdir;git_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
        source_status=subprocess.check_output(['git','status','--porcelain','--','tools/run_p3r2_field_cdf_reference.py'],cwd=ROOT,text=True).strip()
        files=['tools/run_p3r2_field_cdf_reference.py','src/p3r2/core.py','src/p1_material/prior.py','tools/run_p3r2_posterior.py']
        for name in files:
            path=dest/'source'/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes((ROOT/name).read_bytes())
        write(dest/'SOURCE_PROVENANCE.json',dict(git_revision=git_revision,executed_script_worktree_status=source_status,
            exact_executed_script='source/tools/run_p3r2_field_cdf_reference.py',source_files=files,no_hash_check=True,
            interpretation='Copied source bytes are authoritative; Git HEAD alone does not certify a dirty script',python=sys.version,numpy=np.__version__))
        # Copy the paid frame and selection moments, preserving the old job.
        for name in ('FRAME.npz','GH_97_MOMENTS.npz','PROTOCOL.json'):
            copied=dest/'paid_source_snapshot'/name;copied.parent.mkdir(parents=True,exist_ok=True);copied.write_bytes((source/name).read_bytes())
        target=Target(book);y=target.observation('DEV_030')
        if target.cfg['target_id']!=oldprotocol['target_id']: raise ValueError('original target identity changed')
        protocol=dict(job=args.job,stage='B',parent='DEV_030',target_id=target.cfg['target_id'],source_job=str(source),
            target_family='SAME_SLICE_GAUSSIAN_FIELD2_FIXED_ORTHOGONAL62_AND_PROPERTIES',original_FULL_posterior=False,
            q0=q0,V=V,rest=rest,a0=a0,source_GH97_coordinate_mean=oldmean,source_GH97_coordinate_SD=sd,
            source_GH97_max_weight_node_index=maxindex,source_GH97_max_weight_node=maxnode,
            transform='t=NormalCDF(a), a=NormalInverseCDF(t); normalized Gaussian2 prior exactly cancels change-of-variable Jacobian',
            integration_target='exp(-Phi(FullMaxwell(rest+V*NormalInverseCDF(t),q0.u),online_y)) on entire(0,1)^2',
            domain='Entire unit square, including outer prior tails; no finite-a truncation, no added prior/Jacobian/Gram factor',
            cut_rule='0,1, NormalCDF(GH97mean+k*GH97SD) for k=-2,-1,0,1,2, NormalCDF(a0), NormalCDF(maxweightGH97node); exact dedup',
            cuts_t_axis0=cuts[0],cuts_t_axis1=cuts[1],patch_shape=shape,patch_count=len(patches),orders=orders,
            CDF_method='At interior cut values, sum whole-patch likelihood integrals and normalize by total integral; direct integrals, not node ECDF',
            qualification_rule='Two consecutive refinements meet material-meanRMS<=meanTol,abslogZdiff<=logZTol,and both direct-CDFmaxdiff<=CDFTol',
            material_mean_tolerance=args.mean_tolerance,logZ_tolerance=args.logZ_tolerance,CDF_tolerance=args.CDF_tolerance,
            qualification_scope='Empirical conditional-SLICE numerical refinement only; no continuum/global certification or FULL qualification',
            new_MCMC=False,source_chains_read_only=True,GT_used=False,online_observation_source_isolated=True,
            maximum_new_states=maximum_states,maximum_new_RHS=maximum_rhs,pilot_RHS_cap=args.rhs_cap,max_wall_seconds=args.max_wall,
            memory_policy='One patch at a time; all full endpoints immediately persisted, per-node receipts flushed, complete patch raw arrays saved',
            git_revision=git_revision,executed_script_worktree_status=source_status,executed_source_snapshot='source/tools/run_p3r2_field_cdf_reference.py')
        freeze(dest/'PROTOCOL.json',protocol)
        save(dest/'FROZEN_PARTITION.npz',cut_t_axis0=cuts[0],cut_t_axis1=cuts[1],CDF_t_cut_axis0=cuts[0][1:-1],CDF_t_cut_axis1=cuts[1][1:-1],
            CDF_a_cut_axis0=ndtri(cuts[0][1:-1]),CDF_a_cut_axis1=ndtri(cuts[1][1:-1]),source_mean=oldmean,source_SD=sd,a0=a0,maximum_node=maxnode)
        controls=pure_controls();write(dest/'UNIT_SQUARE_DIRECT_CDF_CONTROLS.json',controls)
        if not all(controls['checks'].values()): raise RuntimeError('full-prior/direct-CDF control failed')
        book._update('PAID_FIELD_FRAME_REUSE',source_job=str(source),source_FRAME=str(source/'FRAME.npz'),source_GH97=str(source/'GH_97_MOMENTS.npz'),new_setup_RHS=0)
        events=(dest/'NODE_EVENTS.jsonl').open('x');previous=None;stable=0
        for order in orders:
            patch_summaries=[];order_start=book.pilot_rhs
            for patch_index,(ix,iy,bounds) in enumerate(patches):
                book.cost_owner=f'order_{order}/patch_{patch_index:03d}';patch_start=book.pilot_rhs
                t,a,lw=patch_rule(order,bounds);n=len(t);Phi=np.empty(n);chi=np.empty((n,576),complex);strata=np.empty((n,576),np.int8);keys=np.empty(n,dtype='U128')
                patchdir=dest/f'GL_{order}'/f'patch_{patch_index:03d}';patchdir.mkdir(parents=True)
                save(patchdir/'PLAN.npz',t=t,a=a,log_prior_uniform_quadrature_weights=lw,bounds=np.asarray(bounds),axis_patch_index=np.array([ix,iy]))
                for node in range(n):
                    if book.pilot_rhs+8>args.rhs_cap: raise BudgetStop('stop before next Full solve due independent field-CDF cap')
                    q=q0.copy();q[:64]=rest+V@a[node];before=book.pilot_rhs
                    p=target.evaluate(q,purpose=f'field2_CDF_GL{order}_patch{patch_index:03d}_node{node:04d}')
                    Phi[node]=target.phi(p,y);chi[node]=p.state.chi;keys[node]=p.key
                    _,decoded=target.prior.decode(q);field=decoded['field'];s=np.where(field<=-target.prior.temperature,-1,np.where(field>=target.prior.temperature,1,0)).astype(np.int8);s[~target.prior.mask]=-1;strata[node]=s
                    del p
                    if target.model is not None: target.model._transfers.clear()
                    if book.pilot_rhs-before!=8: raise RuntimeError('exact Full endpoint cost differs from8 RHS')
                    if not np.isfinite(Phi[node]): raise FloatingPointError('nonfinite exact Full likelihood')
                    receipt=dict(order=order,patch=patch_index,patch_axis_index=[ix,iy],node=node,t=t[node],a=a[node],Phi=Phi[node],endpoint=keys[node],
                        log_quadrature_weight=lw[node],actual_endpoint_RHS=8,pilot_cumulative_RHS=book.pilot_rhs)
                    events.write(json.dumps(finite(receipt),allow_nan=False,separators=(',',':'))+'\n');events.flush()
                summary=summarize_patch(a,Phi,chi,strata,lw);patch_summaries.append(summary)
                save(patchdir/'RAW.npz',t=t,a=a,Phi=Phi,chi=chi,stratum=strata,endpoint_key=keys,log_prior_uniform_quadrature_weights=lw)
                save(patchdir/'MOMENTS.npz',**summary)
                write(patchdir/'COST.json',dict(actual_RHS=book.pilot_rhs-patch_start,completed_nodes=n,new_solver_gradients=0,new_receiver_J=0))
                write(dest/'PATCH_PROGRESS.json',dict(order=order,patch_completed=patch_index+1,patch_count=len(patches),actual_new_RHS=book.pilot_rhs))
                del t,a,lw,Phi,chi,strata,keys;gc.collect()
            moments=combine_patches(patch_summaries,shape);meanq=q0.copy();meanq[:64]=rest+V@moments['coordinate_mean'];decoder,decoded=target.prior.decode(meanq)
            gap=moments['mean']-decoder;fraction_mean=np.real((moments['mean']-target.prior.host)/decoded['inclusion_increment'])
            moments.update(conditional_latent_mean_q=meanq,decoder_at_conditional_latent_mean=decoder,decoder_Jensen_gap=gap,
                decoder_Jensen_gap_material_RMS=float(np.linalg.norm(gap)/24),conditional_decoder_fraction_mean=fraction_mean,
                decoder_fraction_at_latent_mean=decoded['fraction'],decoder_fraction_Jensen_gap=fraction_mean-decoded['fraction'],
                CDF_t_cut_axis0=cuts[0][1:-1],CDF_t_cut_axis1=cuts[1][1:-1],CDF_a_cut_axis0=ndtri(cuts[0][1:-1]),CDF_a_cut_axis1=ndtri(cuts[1][1:-1]))
            save(dest/f'GL_{order}_MOMENTS.npz',**moments)
            rms=None if previous is None else float(np.linalg.norm(moments['mean']-previous['mean'])/24)
            logdiff=None if previous is None else float(abs(moments['log_integral']-previous['log_integral']))
            cdfdiff=None if previous is None else np.array([np.max(abs(moments['directly_integrated_marginal_CDF_axis0']-previous['directly_integrated_marginal_CDF_axis0'])),np.max(abs(moments['directly_integrated_marginal_CDF_axis1']-previous['directly_integrated_marginal_CDF_axis1']))])
            met=previous is not None and rms<=args.mean_tolerance and logdiff<=args.logZ_tolerance and np.all(cdfdiff<=args.CDF_tolerance)
            stable=stable+1 if met else 0
            row=dict(order=order,new_nodes=len(patches)*order*order,actual_order_RHS=book.pilot_rhs-order_start,logZ=moments['log_integral'],
                material_mean_successive_RMS=rms,absolute_logZ_successive_difference=logdiff,direct_marginal_CDF_successive_max_difference=cdfdiff,
                material_variance_successive_RMS=None if previous is None else float(np.linalg.norm(moments['variance']-previous['variance'])/24),
                Jensen_gap_successive_material_RMS=None if previous is None else float(np.linalg.norm(moments['decoder_Jensen_gap']-previous['decoder_Jensen_gap'])/24),
                coordinate_mean_successive_norm=None if previous is None else float(np.linalg.norm(moments['coordinate_mean']-previous['coordinate_mean'])),
                all_refinement_criteria_met=bool(met),consecutive_successful_refinements=stable,
                decoder_Jensen_gap_material_RMS=moments['decoder_Jensen_gap_material_RMS'],
                maximum_patch_posterior_mass=float(moments['normalized_patch_masses'].max()))
            rows.append(row);write(dest/'REFINEMENT_PROGRESS.json',finite(dict(rows=rows,empirical_conditional_refinement_qualified=stable>=2,
                continuum_certification=False,FULL_reference_qualification=False,actual_new_RHS=book.pilot_rhs)))
            print(json.dumps(finite(row)),flush=True);previous=moments
            if stable>=2: break
        comparison=compare_previous(source,moments,cuts,q0,rest,V,target.prior);write(dest/'INDEPENDENT_SAME_TARGET_COMPARISON.json',finite(comparison))
        result=dict(protocol=protocol,status='SAME_TARGET_FIELD2_DIRECT_CDF_REFERENCE_COMPLETE',actual_new_RHS=book.pilot_rhs,
            cost_receipt=book.costs,rows=rows,final_order=order,empirical_conditional_refinement_qualified=stable>=2,
            continuum_certification=False,global_certification=False,FULL_reference_qualification=False,original_FULL_posterior=False,
            final_logZ=moments['log_integral'],final_decoder_Jensen_gap_material_RMS=moments['decoder_Jensen_gap_material_RMS'],
            independent_same_target_comparison=comparison,elapsed_seconds=time.perf_counter()-start,
            conclusion='Conditional field2 numerical refinement and nonlinear decoder Jensen diagnostic only; no new MCMC, no FULL calibration claim')
        write(dest/'RESULT.json',finite(result));print(json.dumps(finite(result)),flush=True)
    except BaseException as exc:
        status='ERROR';error=repr(exc);write(book.jobdir/'ERROR.json',dict(error=error,actual_new_RHS=book.pilot_rhs,cost_receipt=book.costs,
            completed_refinement_rows=rows,completed_full_endpoints_and_node_receipts_preserved=True,no_automatic_replay=True));raise
    finally:
        if events is not None: events.close()
        book.cost_owner='closure';book.close(status,error)


if __name__=='__main__': main()
