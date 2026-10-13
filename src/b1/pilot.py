"""Fresh synthetic vector-Maxwell evidence pilot. No historical arrays read."""
from __future__ import annotations
import json,time
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from scipy.stats import beta
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import minimize
from .physics import Forward,angular_data
from .budget import Ledger
from .inference import infer,quadrature_nodes,refinement,top_k,field_summaries,log_likelihood,weighted_quantiles


def dump(p,d):
    def convert(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        raise TypeError(type(v).__name__)
    Path(p).parent.mkdir(parents=True,exist_ok=True)
    Path(p).write_text(json.dumps(d,default=convert,indent=2,allow_nan=False)+'\n')

def compact(r):return {k:v for k,v in r.items() if k!='weights'}


def make_bank(root,c,f,order):
    out=root/f'results/banks/order{order}.npz'
    if out.exists():return dict(np.load(out))
    start=time.perf_counter();th,w=quadrature_nodes(order,c);H=[];T=[];W=[];F=[];D=[];actual_before=f.ledger.count()
    for h in range(3):
        for theta,weight in zip(th,w):
            data,a,b,_=f.run(h,*theta,kind=f'QUADRATURE_ORDER_{order}')
            H.append(h);T.append(theta);W.append(weight);F.append(data)
            D.append(angular_data(a[:1],b[:1],c['angles_degrees'],c['wavenumber']))
    bank=dict(H=np.array(H),theta=np.asarray(T),weights=np.array(W),F=np.array(F),F_dipole=np.array(D),
              bounds=np.array([c['material_prior']['eps_real'],c['material_prior']['eps_imag']]))
    out.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(out,**bank)
    dump(out.with_suffix('.json'),dict(order=order,states=len(F),new_full_calls=f.ledger.count()-actual_before,wall_seconds=time.perf_counter()-start,prior='continuous uniform product, quadrature weights normalized within branch'))
    return bank


def save_case(root,name,c,bank,data,sigma,pi,result,reference,*,gain_sd=0.):
    p=root/'results/cases'/name;p.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(p/'posterior.npz',data=data,noise_real_sd=sigma,topology_prior=pi,gain_sd=gain_sd,posterior_weights=result['weights'],theta=bank['theta'],H=bank['H'],forward=bank['F'])
    fields=field_summaries(result,bank,c);np.savez_compressed(p/'material_fields.npz',**{k:v for k,v in fields.items() if isinstance(v,np.ndarray)})
    dump(p/'probabilities.json',dict(case=name,target='3 fixed geometry branches, uniform continuous complex permittivity prior, real Gaussian scattering noise',
                                    likelihood=dict(real_noise_sd=sigma,gain_sd=gain_sd,covariance='sigma^2 I + gain_sd^2 F F^T' if gain_sd else 'sigma^2 I'),result=compact(result),reference=reference,top_k=top_k(result,bank,c['top_k_display']),
                                    cavity_location_scope=fields['cavity_center_probability'],conditional_quantile_method='piecewise constant on normalized prior-mass quadrature bins; no rigorous interval certificate'))


def analytic_proposal(coarse,r,h,c):
    ids=np.flatnonzero(coarse['H']==h);th=coarse['theta'][ids]; lw=r['weights'][ids];lw/=lw.sum()
    lo=np.array([c['material_prior']['eps_real'][0],c['material_prior']['eps_imag'][0]]);hi=np.array([c['material_prior']['eps_real'][1],c['material_prior']['eps_imag'][1]])
    z=(th-lo)/(hi-lo);mu=lw@z;var=lw@((z-mu)**2)
    strength=np.maximum(mu*(1-mu)/np.maximum(var,1e-8)-1,1.)
    aa=np.clip(mu*strength,.5,30.);bb=np.clip((1-mu)*strength,.5,30.)
    return aa,bb,lo,hi


def importance(root,c,f,coarse,fine,cases):
    p=root/'results/importance';p.mkdir(parents=True,exist_ok=True);rng=np.random.default_rng(c['seed']+200);case=cases['competitive'];y=case['data'];sigma=case['sigma'];pi=np.array(case['prior'])
    cr=infer(coarse,y,sigma,pi);before=f.ledger.count();start=time.perf_counter();allrecords=[];zh=[];bootstrap=[];diagnostics=[]
    for h in range(3):
        aa,bb,lo,hi=analytic_proposal(coarse,cr,h,c);points=[];qs=[]
        for j in range(16):
            z=rng.uniform(size=2) if rng.uniform()<.5 else rng.beta(aa,bb)
            q=.5+.5*np.prod(beta.pdf(z,aa,bb));theta=lo+z*(hi-lo)
            response=f.run(h,*theta,kind='DEFENSIVE_IMPORTANCE')[0]
            ll,_=log_likelihood(response[None],y,sigma);logweight=float(ll[0]-np.log(q))
            allrecords.append(dict(H=h,theta=theta,proposal_density_z=float(q),logweight=logweight));points.append(logweight);qs.append(q)
        points=np.array(points);lh=logsumexp(points)-np.log(len(points));zh.append(lh)
        w=np.exp(points-logsumexp(points));diagnostics.append(dict(H=h,ESS=float(1/(w@w)),max_normalized_weight=float(w.max()),beta_a=aa,beta_b=bb,defensive_prior_weight=.5))
        idx=rng.integers(0,len(points),size=(2000,len(points)));bootstrap.append(logsumexp(points[idx],axis=1)-np.log(len(points)))
    zh=np.array(zh);post=np.exp(np.log(pi)+zh-logsumexp(np.log(pi)+zh));bs=np.asarray(bootstrap).T+np.log(pi);bp=np.exp(bs-logsumexp(bs,axis=1)[:,None])
    ref=infer(fine,y,sigma,pi)
    result=dict(case='competitive',target='same continuous full vector-Mie posterior as quadrature',proposal='stratified 16 per topology, q_H(z)=0.5 uniform +0.5 product Beta fitted only to order3 posterior',
                correction='proper continuous prior density divided by q_H; each branch evidence is mean L/q, followed by declared topology prior',
                topology_probabilities=post,log_evidence=zh,probability_error_vs_order7=float(max(abs(post-ref['topology_probabilities']))),
                bootstrap_95_intervals=np.quantile(bp,[.025,.975],axis=0).T,diagnostics=diagnostics,new_full_calls=f.ledger.count()-before,
                coarse_setup_full_calls=27,wall_seconds=time.perf_counter()-start,finite_N_normalized_ratio_biased=True,
                missed_mode_warning='defensive support covers all 3 branches, but finite 16 samples/branch can miss within-branch modes; bootstrap does not certify absent modes')
    dump(p/'RESULT.json',result);dump(p/'SAMPLES.json',allrecords);return result


def comparisons(root,c,f,bank,cases):
    result={};before=f.ledger.count();start=time.perf_counter()
    for name in ['distinguishable','competitive','representation_mismatch']:
        case=cases[name];y=case['data'];sigma=case['sigma'];pi=np.array(case['prior']);logl,mah=log_likelihood(bank['F'],y,sigma)
        minima=np.array([mah[bank['H']==h].min() for h in range(3)])
        wrong=np.exp(-.5*minima+np.log(pi)-logsumexp(-.5*minima+np.log(pi)))
        Hmap=int(bank['H'][np.argmax(logl+np.log(pi[bank['H']]))]);opts=[]
        for h in range(3):
            ids=np.flatnonzero(bank['H']==h);th=bank['theta'][ids];xx=np.unique(th[:,0]);yy=np.unique(th[:,1]);vals=bank['F'][ids].reshape(len(xx),len(yy),-1)
            spl=[RectBivariateSpline(xx,yy,vals[:,:,j],kx=3,ky=3,s=0) for j in range(vals.shape[2])]
            def objective(theta):
                v=np.array([s(*theta)[0,0] for s in spl]);return float(-log_likelihood(v[None],y,sigma)[0][0])
            starts=[th[np.argmax(logl[ids])],[2.5,.1],[2.25,.06],[2.75,.14]]
            solutions=[minimize(objective,s,method='L-BFGS-B',bounds=[(xx[0],xx[-1]),(yy[0],yy[-1])]) for s in starts]
            selected=min(solutions,key=lambda r:r.fun);full=f.run(h,*selected.x,kind='MULTISTART_SELECTED_FULL_VERIFICATION')[0];verified=log_likelihood(full[None],y,sigma)[0][0]
            opts.append(dict(H=h,theta=selected.x,surrogate_log_likelihood=float(-selected.fun),verified_full_log_likelihood=float(verified),successful_starts=sum(r.success for r in solutions),starts=len(starts),probability=None))
        result[name]=dict(single_MAP_topology=Hmap,single_MAP_is_point_decision=True,grid_minimum_whitened_residual=minima,
                          residual_peak_softmax_HEURISTIC_NOT_POSTERIOR=wrong,
                          multiple_start_optimization=opts,optimization_scope='four starts per branch on frozen cubic interpolant within order7 node box; selected endpoints checked by full Maxwell. Not a global/full optimizer or posterior estimator')
    result['resources']=dict(new_physics=f.ledger.count()-before,wall_seconds=time.perf_counter()-start)
    dump(root/'results/POINT_ESTIMATOR_COMPARISONS.json',result);return result


def model_mechanism(root,c,f,bank,cases):
    before=f.ledger.count();start=time.perf_counter();linear=np.zeros_like(bank['F']);diplinear=linear.copy();ranks=[];centers=[];fullJ=[];cheapJ=[]
    th0=np.array([2.5,.1]);step=1e-3
    for h in range(3):
        center=f.run(h,*th0,kind='LINEAR_ANCHOR')[0];JF=[];JD=[]
        for j in range(2):
            plus=th0.copy();minus=th0.copy();plus[j]+=step;minus[j]-=step
            p=f.run(h,*plus,kind='PAID_LINEAR_DERIVATIVE');m=f.run(h,*minus,kind='PAID_LINEAR_DERIVATIVE')
            JF.append((p[0]-m[0])/(2*step))
            pd=angular_data(p[1][:1],p[2][:1],c['angles_degrees'],c['wavenumber']);md=angular_data(m[1][:1],m[2][:1],c['angles_degrees'],c['wavenumber']);JD.append((pd-md)/(2*step))
        JF=np.array(JF).T;JD=np.array(JD).T;ids=bank['H']==h;s=bank['theta'][ids]-th0
        linear[ids]=center+s@JF.T;diplinear[ids]=center+s@JD.T;centers.append(center);fullJ.append(JF);cheapJ.append(JD);ranks.append(int(np.linalg.matrix_rank(JF-JD)))
    report={}
    for name,case in cases.items():
        y=case['data'];sigma=case['sigma'];pi=case['prior'];actual=infer(bank,y,sigma,pi);bl=dict(bank,F=linear);bc=dict(bank,F=diplinear);ld=dict(bank,F=bank['F_dipole'])
        rl=infer(bl,y,sigma,pi);rc=infer(bc,y,sigma,pi);rd=infer(ld,y,sigma,pi)
        h_all=np.zeros(len(bank['H']));error=[]
        for h in range(3):
            ids=bank['H']==h;s=bank['theta'][ids]-th0;JF=fullJ[h]/sigma;JD=cheapJ[h]/sigma;r=(centers[h]-y)/sigma
            U,sv,Vh=np.linalg.svd(JF-JD,full_matrices=False);Y=U*sv;N=Vh
            aa=s@N.T;bb=s@JD.T@Y;hh=np.sum((r@Y+bb)*aa,axis=1)+.5*np.sum((aa@Y.T)**2,axis=1);h_all[ids]=hh
            p0=.5*np.sum((r+s@JD.T)**2,axis=1);p1=.5*np.sum((r+s@JF.T)**2,axis=1);error.append(float(np.max(abs(p1-p0-hh))))
        lf,_=log_likelihood(bank['F'],y,sigma);la,_=log_likelihood(diplinear,y,sigma);delta=la-lf
        # Full nonlinear potential change = paired affine h + remainder.
        remainder=delta-h_all
        corrected=infer(dict(bank,F=bank['F']),y,sigma,pi)
        log_reweight_branch=[]
        for h in range(3):
            ids=bank['H']==h;logq=np.log(bank['weights'][ids])+la[ids]-rc['log_evidence'][h]
            log_reweight_branch.append(float(logsumexp(logq-delta[ids])))
        predicted=np.array(rc['log_evidence'])+log_reweight_branch
        report[name]=dict(full_topology=actual['topology_probabilities'],single_anchor_full_derivative_linear=rl['topology_probabilities'],single_anchor_dipole_derivative_calibrated_linear=rc['topology_probabilities'],truncated_dipole_model=rd['topology_probabilities'],
                          affine_signature_identity_max_error=max(error),nonlinear_remainder_abs_max=float(max(abs(remainder))),
                          nonlinear_remainder_weighted_RMS=float(np.sqrt(actual['weights']@(remainder*remainder))),
                          full_vs_calibrated_linear_probability_change=float(max(abs(np.array(actual['topology_probabilities'])-rc['topology_probabilities']))),
                          exact_branch_reweight_log_evidence_error=float(max(abs(predicted-actual['log_evidence']))),
                          log_branch_correction=log_reweight_branch,
                          prior_support_remainder_certified=False)
    report['scope']=dict(material_dimension=2,affine_defect_ranks=ranks,signature_dimension_up_to=4,
                         A17_native_Galerkin_exchange_used=False,source='generic SVD factors of two actual Mie affine Jacobians; not A17 current feedback',
                         efficiency='no dimension reduction: signature dimension >= material dimension; remainder uses paid offline full bank',
                         costs=dict(new_physics=f.ledger.count()-before,wall_seconds=time.perf_counter()-start,derivative_material_frequency_calls=12,anchors=3,anchors_cache_hit=True))
    dump(root/'results/MODEL_CORRECTION_MECHANISM.json',report);np.savez_compressed(root/'results/LINEAR_MODELS.npz',F_linear=linear,F_dipole_linear=diplinear,J_full=np.array(fullJ),J_dipole=np.array(cheapJ),centers=np.array(centers),theta0=th0)
    return report


def calibration(root,c,f,coarse,fine):
    start=time.perf_counter();before=f.ledger.count();rng=np.random.default_rng(c['seed']+400);records=[];allp=[];allh=[];ci=[];references=[];predicted=[];captures=[]
    lo=np.array([c['material_prior']['eps_real'][0],c['material_prior']['eps_imag'][0]]);hi=np.array([c['material_prior']['eps_real'][1],c['material_prior']['eps_imag'][1]])
    # 24 independent continuous-material parents, stratified equal topology.
    # 30 independent noise replicates per parent; uncertainty must cluster by parent.
    design=[dict(id=i,H=i//8,theta=(lo+rng.uniform(size=2)*(hi-lo)).tolist()) for i in range(24)]
    dump(root/'results/CALIBRATION_GENERATION_MANIFEST.json',dict(seed=c['seed']+400,parents=design,noise_replicates_per_parent=30,scope='fresh continuous prior generative simulation; development model checks, not locked or external validation'))
    refRMS=np.sqrt(np.mean(fine['F'][fine['H']==0]**2));sigma=.65*refRMS
    for item in design:
        h=item['H'];theta=np.array(item['theta']);response=f.run(h,*theta,kind='FRESH_CONTINUOUS_CALIBRATION_PARENT')[0];batch=[]
        for rep in range(30):
            y=response+sigma*rng.normal(size=len(response));r=infer(fine,y,sigma,[1/3]*3);q=infer(coarse,y,sigma,[1/3]*3);rc=refinement(q,r,c);prob=np.array(r['topology_probabilities']);one=np.eye(3)[h];ordr=np.argsort(-prob);cum=np.cumsum(prob[ordr]);credible=ordr[:np.searchsorted(cum,.9)+1]
            intervals=[]
            for j in range(2):
                intervals.append(weighted_quantiles(fine['theta'][:,j],r['weights'],fine['weights'],fine['bounds'][j]))
            param_cov=[intervals[j][0]<=theta[j]<=intervals[j][2] for j in range(2)]
            B=float(np.sum((prob-one)**2));LS=float(-np.log(max(prob[h],np.finfo(float).tiny)))
            batch.append(dict(rep=rep,probabilities=prob.tolist(),brier=B,log_score=LS,topology_credible_90_contains_truth=bool(h in credible),material_90_contains_truth=param_cov,
                              top2_contains_truth=bool(h in ordr[:2]),top2_mass=float(prob[ordr[:2]].sum()),refinement=rc))
            allp.append(prob);allh.append(h);ci.append(h in credible);references.append(rc['status']=='NUMERICAL_REFERENCE_PASS');predicted.append(prob[ordr[:2]].sum());captures.append(h in ordr[:2])
        records.append(dict(**item,replicates=batch,mean_brier=float(np.mean([r['brier'] for r in batch])),mean_log_score=float(np.mean([r['log_score'] for r in batch]))))
    probs=np.asarray(allp);truth=np.array(allh);brier=np.array([r['mean_brier'] for r in records]);logscore=np.array([r['mean_log_score'] for r in records]);bins=[]
    for h in range(3):
        for low,high in zip(np.linspace(0,1,6)[:-1],np.linspace(0,1,6)[1:]):
            idx=(probs[:,h]>=low)&(probs[:,h]<(high if high<1 else 1+1e-14))
            if idx.any():bins.append(dict(H=h,low=float(low),high=float(high),count=int(idx.sum()),mean_probability=float(probs[idx,h].mean()),observed_frequency=float(np.mean(truth[idx]==h))))
    bootrng=np.random.default_rng(c['seed']+401); resamples=bootrng.integers(0,8,size=(2000,3,8));parent_brier=brier.reshape(3,8);parent_log=logscore.reshape(3,8)
    bootb=np.array([np.mean([parent_brier[h][ix[h]].mean() for h in range(3)]) for ix in resamples]);bootl=np.array([np.mean([parent_log[h][ix[h]].mean() for h in range(3)]) for ix in resamples])
    summary=dict(independent_continuous_material_parents=24,noise_replicates=720,Brier_multiclass=float(brier.mean()),Brier_cluster_bootstrap_95=np.quantile(bootb,[.025,.975]),log_score=float(logscore.mean()),log_score_cluster_bootstrap_95=np.quantile(bootl,[.025,.975]),
                 topology_credible_90_coverage=float(np.mean(ci)),top2_empirical_coverage=float(np.mean(captures)),top2_mean_reported_mass=float(np.mean(predicted)),reference_refinement_pass_fraction=float(np.mean(references)),
                 marginal_material_90_coverage=np.mean([r['material_90_contains_truth'] for parent in records for r in parent['replicates']],axis=0),
                 calibration_bins=bins,full_calls=f.ledger.count()-before,wall_seconds=time.perf_counter()-start,
                 calibration_status='EXPLORATORY_CONTINUOUS_PRIOR_MODEL_CHECK; not real-world calibration, reference convergence separately reported',
                 uncertainty='24 independent material clusters; 720 noise realizations are not 720 independent material objects; cluster CI excludes solver/representation bias')
    dump(root/'results/CALIBRATION_RESULTS.json',summary);dump(root/'results/CALIBRATION_RAW.json',records);return summary


def run(root,c):
    root=Path(root);out=root/'results/PRIMARY_SUMMARY.json'
    if out.exists():raise RuntimeError('primary experiment already completed; use cached summarize, never overwrite ledger/results')
    audit=json.loads((root/'results/FORWARD_AUDIT.json').read_text())
    if audit['status']!='PASS':raise RuntimeError('audit gate closed')
    t=time.perf_counter();led=Ledger(root/'results',c['physics_new_cap']);f=Forward(c,led,root/'results/forward_cache');banks={n:make_bank(root,c,f,n) for n in c['quadrature_orders']};fine=banks[7];coarse=banks[5]
    # Fixed experiment scale from one declared calibration material, not from each candidate.
    ref=f.run(0,2.5,.1,kind='EXPERIMENT_NOISE_SCALE')[0];rms=np.sqrt(np.mean(ref*ref));rng=np.random.default_rng(c['seed']);cases={};summary={}
    true_design={'distinguishable':(0,2.47,.09),'competitive':(1,2.47,.09),'information_limited':(2,2.47,.09),'gain_mismatch':(1,2.47,.09),'representation_mismatch':(1,2.47,.09)}
    for name,(h,er,ei) in true_design.items():
        if name=='representation_mismatch': truth=f.run(h,er,ei,kind='OUTSIDE_TEMPLATE_NEGATIVE_CONTROL',radii=[.68,1.],eps_layers=[1.,complex(er,ei)])[0]
        else:truth=f.run(h,er,ei,kind='FRESH_OBSERVATION_TRUTH')[0]
        if name=='gain_mismatch':truth=1.30*truth
        sigma=c['case_noise_relative_real_rms'][name]*rms;y=truth+sigma*rng.normal(size=len(truth));pi=c['topology_prior'];case=dict(data=y,sigma=sigma,prior=pi,truth_H=h,truth_theta=[er,ei],design=name);cases[name]=case
        r=infer(fine,y,sigma,pi);q=infer(coarse,y,sigma,pi);q3=infer(banks[3],y,sigma,pi);rf=refinement(q,r,c)
        save_case(root,name,c,fine,y,sigma,pi,r,rf)
        summary[name]=dict(result=compact(r),refinement_5_vs_7=rf,refinement_3_vs_5=refinement(q3,q,c),truth=dict(H=h,eps_real=er,eps_imag=ei,gain=1.30 if name=='gain_mismatch' else 1.,within_inference_support=name!='representation_mismatch'),noise_real_sd=float(sigma))
        if name=='gain_mismatch':
            fixed=infer(fine,y,sigma,pi,gain_sd=.15);fixed5=infer(coarse,y,sigma,pi,gain_sd=.15)
            save_case(root,'gain_marginalized',c,fine,y,sigma,pi,fixed,refinement(fixed5,fixed,c),gain_sd=.15);summary['gain_marginalized']=dict(result=compact(fixed),refinement_5_vs_7=refinement(fixed5,fixed,c),nuisance_prior='g~N(1,0.15^2), determinant included')
    # Early complete Top-K files precede optional sampler/mechanism/calibration work.
    dump(root/'results/INITIAL_IMAGING_SUMMARY.json',summary)
    sensitivity={};case=cases['information_limited']
    for name,pi in [('solid_favored',[.8,.15,.05]),('double_void_favored',[.05,.15,.8])]:
        r=infer(fine,case['data'],case['sigma'],pi);q=infer(coarse,case['data'],case['sigma'],pi)
        sensitivity[name]=dict(prior=pi,topology=r['topology_probabilities'],reference=refinement(q,r,c),same_conditional_material_prior=True)
    # Material-prior sensitivity: independent Beta(2,2) real / Beta(2,5) imaginary.
    z=(fine['theta']-fine['bounds'][:,0])/(fine['bounds'][:,1]-fine['bounds'][:,0]);bw=fine['weights']*beta.pdf(z[:,0],2,2)*beta.pdf(z[:,1],2,5)
    bbank=dict(fine,weights=bw);changed=infer(bbank,cases['competitive']['data'],cases['competitive']['sigma'],c['topology_prior'])
    sensitivity['material_beta_prior']=dict(type='proper Beta(2,2) real and Beta(2,5) imag on original support, exact analytic normalizations retained',topology=changed['topology_probabilities'],quadrature_prior_integrals=[float(bw[fine['H']==h].sum()) for h in range(3)],reference_status='order7 only; not independently refined',conditional_credible_intervals='not interpreted for this changed prior')
    dump(root/'results/PRIOR_SENSITIVITY.json',sensitivity)
    IS=importance(root,c,f,banks[3],fine,cases);point=comparisons(root,c,f,fine,cases);mechanism=model_mechanism(root,c,f,fine,cases);cal=calibration(root,c,f,coarse,fine)
    result=dict(protocol='B1-probabilistic-Maxwell-v1',cases=summary,prior_sensitivity=sensitivity,importance=IS,calibration=cal,resources=led.summary(),wall_seconds=time.perf_counter()-t,
                scientific_status='working vector-Mie continuous-material Top-K prototype; reference/paper gates assessed separately',locked_data_opened=False,new_nn_training=False)
    dump(out,result);dump(root/'results/DEPLOYMENT_CONTRACT.json',{k:c[k] for k in ['outer_radius','wavenumber','nmax','angles_degrees','topologies','material_prior','geometry']});dump(root/'results/RESOURCE_SUMMARY.json',led.summary());print(json.dumps(dict(probabilities={n:r['result']['topology_probabilities'] for n,r in summary.items()},references={n:r['refinement_5_vs_7']['status'] for n,r in summary.items()},resources=led.summary(),wall_seconds=result['wall_seconds'])))
