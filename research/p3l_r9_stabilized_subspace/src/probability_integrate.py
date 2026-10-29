"""R9 two-factor numerical R6 integration and whole-support residual audit.

No strict fields are populated: ordinary float Schur/triangular majorants lack
outward coefficient arithmetic. Positive numerical evidence and full-support
exclusion form an explicitly numerical reference qualification.
"""
from probability_transfer import *
from scipy.optimize import least_squares
from scipy.special import logsumexp


def summary(V,w,labels,phase,phis,shift):
    w=np.asarray(w);V=np.asarray(V);phis=np.asarray(phis);a=w*np.exp(shift-phis);Z=float(a.sum());p=a/Z
    mu=np.sum(p*V);variance=float(np.sum(p*abs(V-mu)**2));labels=np.asarray(labels);phase=np.asarray(phase)
    return dict(mean_increment=[mu.real,mu.imag],mean_real_ROI_material=HOST.real+mu.real,material_ROI_SD=float(np.sqrt(variance)),material_variance=variance,source_layer_mass=p@labels,latent_phase_gt_half=float(p@phase),material_real_increment_gt_002=float(p@(V.real>.02)),material_nonhost_event=float(p@(V!=0)),ROI256_hard_event=float(p@(HOST.real+V.real>=.75)),scaled_evidence=Z,log_evidence=float(np.log(Z)-shift),node_count=len(V),strict_evidence_interval=None,strict_mean_interval=None,strict_event_interval=None,qualification='NUMERICAL_INTEGRAL')


def locate(ph,name):
    y=ph.obs[name]
    def residual(x):return pack(ph.reduced(complex(*x))-y)*np.sqrt(2)/SIGMA
    with ph.book.span('fresh_profile_mode_and_GN'):
        fit=least_squares(residual,[.9,.09],bounds=([1e-8,1e-8],[RHI,IHI]),xtol=1e-12,ftol=1e-12,gtol=1e-9,max_nfev=100)
        center=fit.x;C=np.linalg.inv(fit.jac.T@fit.jac);chol=np.linalg.cholesky(C);phi=.5*fit.fun@fit.fun
    return dict(center=center,C=C,chol=chol,phi=phi,fit_nfev=fit.nfev,fit_status=int(fit.status),legal=legal(complex(*center)))


def gh(ph,name,mode,order,kind='reduced',flat=False):
    y=ph.obs[name];center=mode['center'];chol=mode['chol'];C=mode['C'];shift=mode['phi'];gh,w=hermgauss(order)
    V=[];weights=[];labels=[];phase=[];phis=[];maximum_residual=0
    before=dict(ph.book.counts);t=time.perf_counter();cpu=time.process_time()
    with ph.book.span('GH_prior_and_total'):
        for i,x in enumerate(gh):
            for j,z in enumerate(gh):
                standardized=np.sqrt(2)*np.array([x,z]);v=complex(*(center+chol@standardized))
                density,bran,event=prior(v,96)
                if density<=0:continue
                if kind=='direct':f,res=ph.direct(v);maximum_residual=max(maximum_residual,res)
                else:f=ph.reduced(v)
                phi=float(np.sum(abs(f-y)**2)/SIGMA**2)
                logq=-np.log(2*np.pi)-.5*np.linalg.slogdet(C)[1]-.5*standardized@standardized
                weights.append(w[i]*w[j]/np.pi*np.exp(-logq)*(1 if flat else density))
                V.append(v);labels.append(bran/density);phase.append(event/density);phis.append(phi)
    result=summary(V,weights,labels,phase,phis,shift)
    result.update(method='likelihood_GN_adapted_GH_flat_material_WRONG' if flat else 'likelihood_GN_adapted_GH_correct_quotient_prior',physics=kind,order=order,wall_seconds=time.perf_counter()-t,CPU_seconds=time.process_time()-cpu,max_direct_backward_residual=maximum_residual,cost_counts_delta={k:ph.book.counts[k]-before.get(k,0) for k in ph.book.counts})
    stem=f'{name}_{kind}_GH{order}'+('_flat_wrong' if flat else '')
    np.savez_compressed(RES/(stem+'.npz'),V=V,prior_importance_weights=weights,source_conditional=labels,phase_conditional=phase,phi=phis)
    write(RES/(stem+'.json'),result);return result


def source_integral(ph,name,mode,nouter,ninner):
    """Independent original (z0,u64,u65) integration, split host/pure/transition.

    A retained material rectangle determines property integration limits. The
    complement is handled by whole-support residual audit, never deleted.
    No Q^-2 factor occurs here: integrating the original source prior directly.
    """
    before=dict(ph.book.counts);t=time.perf_counter();cpu=time.process_time();y=ph.obs[name];center=mode['center'];sd=np.sqrt(np.diag(mode['C']));lo=center-12*sd;hi=center+12*sd
    qlo=max(lo[0]/RHI,lo[1]/IHI,0);qhi=min(hi[0]/RLO,hi[1]/ILO,1.)
    zl,zh=(2*betaincinv(3,3,np.array([qlo,qhi]))-1)/1.3
    # Cut where a property rectangle face switches; seam/event cuts explicit.
    cuts=[zl,zh,0.]
    for q in [lo[0]/RLO,lo[0]/RHI,hi[0]/RLO,hi[0]/RHI,lo[1]/ILO,lo[1]/IHI,hi[1]/ILO,hi[1]/IHI]:
        if qlo<q<qhi:cuts.append(float((2*betaincinv(3,3,q)-1)/1.3))
    cuts=sorted(set(z for z in cuts if zl<=z<=zh))
    V=[];weights=[];labels=[];phase=[];phis=[]
    with ph.book.span('original_source_conditional_integration'):
        def properties(q,zw,layer,zphase):
            rlo=max(ULO,lo[0]/q-.8);rhi=min(UHI,hi[0]/q-.8)
            ilo=max(ULO,(lo[1]/q-.03)/.15);ihi=min(UHI,(hi[1]/q-.03)/.15)
            if rhi<=rlo or ihi<=ilo:return
            rr,wr=gl(ninner,rlo,rhi);ii,wi=gl(ninner,ilo,ihi)
            for r,a in zip(rr,wr):
                for i,b in zip(ii,wi):
                    v=q*complex(.8+r,.03+.15*i);f=ph.reduced(v)
                    V.append(v);weights.append(zw*a*b/(UHI-ULO)**2);labels.append(np.eye(3)[layer]);phase.append(float(zphase>0));phis.append(float(np.sum(abs(f-y)**2)/SIGMA**2))
        for left,right in zip(cuts[:-1],cuts[1:]):
            zz,ww=gl(nouter,left,right)
            for z,w in zip(zz,ww):properties(float(sat(1.3*z)),float(w*norm.pdf(z)/Z0),1,z)
        if lo[0]<RHI and hi[0]>RLO and lo[1]<IHI and hi[1]>ILO:properties(1.,PH,2,1.)
        # host atom explicitly retained; underflow not interpreted as strict zero.
        V.append(0j);weights.append(PH);labels.append(np.eye(3)[0]);phase.append(0.);phis.append(float(np.sum(abs(ph.Fh-y)**2)/SIGMA**2))
    result=summary(V,weights,labels,phase,phis,mode['phi']);result.update(method='original_latent_conditional_GL',outer_order=nouter,inner_order=ninner,wall_seconds=time.perf_counter()-t,CPU_seconds=time.process_time()-cpu,retained_rectangle=[lo,hi],unintegrated_complement='see whole_support_audit.json; not silently renormalized',cost_counts_delta={k:ph.book.counts[k]-before.get(k,0) for k in ph.book.counts})
    write(RES/f'{name}_source_GL{nouter}_{ninner}.json',result);return result


def global_audit(ph,name,mode,maxcells=150000):
    """Exhaustive rectangular cover with exact-algebra/float majorant.

    All source strata map inside the full rectangle. Numerically excluded
    complement has L/e^-phi_mode <= exp(-50), so its PRIOR integrated evidence
    <= exp(-50) *1. This avoids summing a fictitious unit mass for each cell.
    The source-coordinate retained integration supplies positive denominator.
    """
    start=time.perf_counter();cpu=time.process_time();y=ph.obs[name];sd=np.sqrt(np.diag(mode['C']));lo=mode['center']-12*sd;hi=mode['center']+12*sd
    assert 0<lo[0]<hi[0]<RHI and 0<lo[1]<hi[1]<IHI, 'retained cover must lie strictly within outer support rectangle'
    xcuts=[0.,float(lo[0]),float(hi[0]),RHI];ycuts=[0.,float(lo[1]),float(hi[1]),IHI]
    todo=[];retained=[]
    for i in range(3):
        for j in range(3):
            box=[xcuts[i],xcuts[i+1],ycuts[j],ycuts[j+1]]
            if i==j==1:retained.append(box)
            else:todo.append((box,0))
    excluded=[];unresolved=[];minsep=np.inf;cells=0;maxdepth=0
    with ph.book.span('whole_support_majorant_audit'):
        while todo:
            box,depth=todo.pop();cells+=1;maxdepth=max(maxdepth,depth)
            x0,x1,y0,y1=box;a=complex((x0+x1)/2,(y0+y1)/2);hr=(x1-x0)/2;hi2=(y1-y0)/2
            major=ph.majorant_derivative(a,hr,hi2)
            if major is not None:
                f=ph.reduced(a);rn=float(np.linalg.norm(pack(f-y))*np.sqrt(2)/SIGMA);distance=max(0,rn-major['white_output_radius']);lower=.5*distance**2;gap=lower-mode['phi']
            else:gap=-np.inf
            if gap>=50:
                minsep=min(minsep,gap);excluded.append(dict(box=box,depth=depth,anchor=[a.real,a.imag],white_anchor_residual_norm=rn,white_response_radius_numerical=major['white_output_radius'],diagonal_margin_numerical=major['minimum_diagonal_margin'],response_radius=major['radius'],phi_gap_lower_numerical=gap,prior_mass_upper='joint prior of disjoint excluded union <=1; no per-cell unit mass addition'))
            elif cells+len(todo)>=maxcells or depth>=30:
                unresolved.append(dict(box=box,depth=depth,reason='majorant_or_budget_unclosed'))
            else:
                if hr>=hi2:
                    mid=(x0+x1)/2;todo.extend([([x0,mid,y0,y1],depth+1),([mid,x1,y0,y1],depth+1)])
                else:
                    mid=(y0+y1)/2;todo.extend([([x0,x1,y0,mid],depth+1),([x0,x1,mid,y1],depth+1)])
            if cells%10000==0:print(json.dumps(dict(event='GLOBAL_SUPPORT_AUDIT',case=name,cells=cells,pending=len(todo),unresolved=len(unresolved))),flush=True)
    result=dict(case=name,status='NUMERICAL_WHOLE_SUPPORT_COVER_COMPLETE' if not unresolved else 'UNRESOLVED',qualification='ordinary float evaluation of proven triangle-majorant algebra; coefficient and rounding enclosures absent, NOT STRICT',global_support_rectangle=[0,RHI,0,IHI],retained_rectangle=[lo,hi],retained_source_layers=['transition','pure'],host_atom='explicit prior PH in source integral; outside retained rectangle and included in disjoint union bound',evaluated_cells=cells,excluded_cells=len(excluded),maxdepth=maxdepth,unresolved_cells=unresolved,numerical_scaled_omitted_evidence_upper=float(np.exp(-minsep)) if not unresolved else None,minimum_excluded_phi_gap= float(minsep) if np.isfinite(minsep) else None,strict_omitted_mass_upper=None,wall_seconds=time.perf_counter()-start,CPU_seconds=time.process_time()-cpu)
    write(RES/f'{name}_whole_support_audit.json',result)
    with (RES/f'{name}_excluded_cells.jsonl').open('w') as f:
        for row in excluded:f.write(json.dumps(plain(row))+'\n')
    return result


def run():
    if (RES/'integration_receipt.json').exists():raise RuntimeError('integration job already executed; no automatic replay')
    b=Book('integration');status='FAILED';error=None;results={}
    try:
        ph=Physics(b);ph.load_setup()
        for name in ph.obs:
            mode=locate(ph,name);write(RES/f'{name}_fresh_mode.json',mode)
            assert mode['legal']
            print(json.dumps(dict(event='MODE',case=name,center=mode['center'].tolist(),phi=float(mode['phi']))),flush=True)
            ghref=[gh(ph,name,mode,n) for n in (21,31)]
            # Different source-coordinate integration with different refinement.
            sources=[source_integral(ph,name,mode,no,ni) for no,ni in [(12,20),(20,32)]]
            audit=global_audit(ph,name,mode)
            meanerr=max(np.max(np.abs(np.array(sources[-1]['mean_increment'])-np.array(ghref[-1]['mean_increment']))),np.max(np.abs(np.array(sources[-1]['mean_increment'])-np.array(sources[-2]['mean_increment']))),np.max(np.abs(np.array(ghref[-1]['mean_increment'])-np.array(ghref[-2]['mean_increment']))))
            sd=ghref[-1]['material_ROI_SD'];keys=['source_layer_mass','latent_phase_gt_half','material_real_increment_gt_002','material_nonhost_event','ROI256_hard_event'];eventerr=max(max(float(np.max(np.abs(np.asarray(sources[-1][k])-np.asarray(ghref[-1][k])))),float(np.max(np.abs(np.asarray(sources[-1][k])-np.asarray(sources[-2][k])))),float(np.max(np.abs(np.asarray(ghref[-1][k])-np.asarray(ghref[-2][k]))))) for k in keys)
            if audit['numerical_scaled_omitted_evidence_upper'] is not None:
                numerical_Z_lower=min(sources[-1]['scaled_evidence'],ghref[-1]['scaled_evidence'])*.999
                omitted=audit['numerical_scaled_omitted_evidence_upper']/numerical_Z_lower
            else:omitted=None
            qualified=omitted is not None and meanerr+RHI*omitted<=2e-6 and (meanerr+RHI*omitted)/sd<=.002 and eventerr+omitted<=.0004
            qualification=dict(status='NUMERICAL_REFERENCE_QUALIFIED' if qualified else 'REFERENCE_UNQUALIFIED',mean_max_absolute_cross_refinement=meanerr,mean_error_over_SD=meanerr/sd,event_layer_max_absolute_cross_coordinate=eventerr,numerical_omitted_mass_upper=omitted,strict_task_error_certificate=False,strict_reference_halfwidth=None,evidence_positive_numerical_lower_estimate=numerical_Z_lower if omitted is not None else None,scope='same fixed R6 S3 conditional model; numerical reference only',source_prior_and_atom_complete=True)
            write(RES/f'{name}_reference_qualification.json',qualification)
            baselines=[]
            for n in (9,11):
                baselines.append(gh(ph,name,mode,n,'reduced'))
                baselines.append(gh(ph,name,mode,n,'direct'))
            flat=gh(ph,name,mode,31,'reduced',flat=True)
            flatgap=float(np.linalg.norm(np.array(flat['mean_increment'])-np.array(ghref[-1]['mean_increment']))/sd)
            results[name]=dict(reference=ghref[-1],qualification=qualification,baselines=baselines,flat_prior_wrong_mean_shift_SD=flatgap,whole_support_audit=audit,source_coordinate_integrals=sources)
            write(RES/'integration_results_partial.json',results)
            print(json.dumps(dict(event='CASE_COMPLETE',case=name,reference_status=qualification['status'],meanerr=meanerr,eventerr=eventerr,used_full_RHS=b.full_rhs)),flush=True)
        write(RES/'integration_results.json',results);status='COMPLETE'
    except BaseException as e:error=repr(e);write(AUD/'integration_failure.json',dict(error=error));raise
    finally:b.close(status,error)

if __name__=='__main__':run()
