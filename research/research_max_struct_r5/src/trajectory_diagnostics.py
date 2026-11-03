"""R5 fixed-state trajectory observations; no truth, tuning, or imaging corrector.

All material columns are real and shared across sources. Complex field residuals
are realified after multiplication by sqrt(ws); data whitening is separate.
"""
from collections import Counter
from copy import copy
import time
import numpy as np
from scipy import linalg as la
from research.research_max_struct_r4.src.cie_onehop import (
    OneHopAbsorption, onehop_material_value_gradient, onehop_material_readout)
from research.research_max_struct_r2.src.cie_global import _baseline_readout
from research.research_max_struct_r4.src.cooperative_geometry import inherited_hvp


def cosine(a, b):
    a, b = np.asarray(a).ravel(), np.asarray(b).ravel()
    den = np.linalg.norm(a)*np.linalg.norm(b)
    return float(np.vdot(a,b).real/den) if den > 0 else None


def vector_summary(x, step, negphys):
    norm = float(np.linalg.norm(x))
    return dict(norm=norm, cosine_actual_step=cosine(x,step),
                cosine_negative_physical_gradient=cosine(x,negphys),
                dot_actual_step=float(np.dot(x,step)),
                dot_negative_physical_gradient=float(np.dot(x,negphys)))


def difference_counter(after, before):
    return {k:after.get(k,0)-before.get(k,0) for k in set(after)|set(before)
            if isinstance(after.get(k,0),(float,int)) and isinstance(before.get(k,0),(float,int))}


def phase_blocks(base, sign):
    blocks = copy(base)
    blocks.goff = sign*base.goff
    return blocks


def frozen_ls(blocks, m, Q, r, beta, ws, kind, cost):
    """Block 4x4 real normal equations, without a residual-space projector."""
    ids = blocks.indices
    q = Q[ids]
    sw = np.sqrt(ws)
    local_den = np.sum(abs(Q)**2*ws, axis=1)
    local_rhs = np.sum(np.real(Q.conj()*r)*ws, axis=1)
    d0 = np.divide(local_rhs,local_den,out=np.zeros_like(m),where=local_den>0)
    unexplained = r-Q*d0[:,None]
    d = 1+(beta-blocks.diag_g)*m
    a = m/d
    if kind == 'ordinary':
        p = np.zeros_like(blocks.goff)
        p[:,np.arange(4),np.arange(4)] = (beta/(1+beta*m))[ids]
    else:
        sign = {'self':0.,'plus':1.,'minus':-1.}[kind]
        p = sign*a[ids][:,:,None]*blocks.goff
        p[:,np.arange(4),np.arange(4)] += 1
        p *= (beta/d)[ids][:,None,:]
    # K[b,output cell,source,real material column]. Real normal equations use
    # Re K* K, so complex sources never become independent material parameters.
    K = np.einsum('boi,bis,s->bosi',p,q,sw)
    pr = np.einsum('boi,bis,s->bos',p,r[ids],sw)
    pe = np.einsum('boi,bis,s->bos',p,unexplained[ids],sw)
    gram = np.einsum('bosi,bosj->bij',K.conj(),K).real
    rhs = np.einsum('bosi,bos->bi',K.conj(),pr).real
    transport_rhs = np.einsum('bosi,bos->bi',K.conj(),pe).real
    eig = np.linalg.eigvalsh(gram)
    condition = np.divide(eig[:,-1],eig[:,0],out=np.full(len(ids),np.inf),where=eig[:,0]>0)
    solutions=[];fallback=[];ranks=[]
    for b,(g,z) in enumerate(zip(gram,rhs)):
        # Ordinary LAPACK precision only; no custom rank threshold or exact-kernel claim.
        if eig[b,0] <= np.finfo(float).eps*max(eig[b,-1],0)*4:
            x,_,rank,_ = la.lstsq(g,z,cond=None,lapack_driver='gelsd')
            fallback.append(b);ranks.append(int(rank));cost['small4_real_Gram_lstsq_fallbacks']+=1
        else:
            try:
                x=la.solve(g,z,assume_a='pos',check_finite=False)
            except la.LinAlgError:
                x,_,rank,_=la.lstsq(g,z,cond=None,lapack_driver='gelsd')
                fallback.append(b);ranks.append(int(rank));cost['small4_real_Gram_lstsq_fallbacks']+=1
        solutions.append(x)
    cost['small4_real_Gram_solves']+=len(ids)
    cost['small4_real_Gram_eigendecompositions']+=len(ids)
    cost['small4_real_Gram_builds']+=len(ids)
    cost['small4_metric_complex_RHS']+=2*len(ids)*Q.shape[1]
    ds=np.empty_like(m);ds[ids]=solutions
    tr=np.empty_like(m);tr[ids]=transport_rhs
    raw_transport=ds-d0
    raw_rhs=tr.copy()
    analytic_zero=kind in ('ordinary','self')
    if analytic_zero:
        # Nonzero scalar cell weights cancel in unrestricted shared-real LS.
        # Do not normalize the subtraction roundoff into a spurious direction.
        transport=np.zeros_like(ds);tr=np.zeros_like(tr)
    else:
        transport=raw_transport
    return dict(d=ds,d0=d0,transport=transport,transport_rhs=tr,
        transport_raw_roundoff=raw_transport if analytic_zero else np.zeros_like(ds),
        transport_rhs_raw_roundoff=raw_rhs if analytic_zero else np.zeros_like(tr),gram=gram,
        gram_eigenvalues=eig,gram_condition=condition,local_C_denominator=local_den,
        fallback_blocks=np.array(fallback,int),fallback_ranks=np.array(ranks,int))


def endpoint_geometry(model, state, Y, step, cost):
    wall,cpu=time.perf_counter(),time.process_time()
    Jr=model.jacobian(state['m'])
    Jc=state['Jcomplex']
    scale=float(np.vdot(Y,Y).real)
    real_gram=(Jr.T@Jr)/scale
    complex_gram=(Jc.conj().T@Jc)/scale
    er=la.eigvalsh(real_gram,check_finite=False)
    ec=la.eigvalsh(complex_gram,check_finite=False)
    cost['physical_real_Jacobian_Gram_builds']+=1
    cost['physical_complex_Jacobian_Gram_builds']+=1
    cost['physical_Jacobian_Gram_eigendecompositions']+=2
    cost['physical_Jacobian_SVD']=0
    vectors=dict(real_gram_eigenvalues=er,complex_gram_eigenvalues=ec)
    if step is None:
        gn=exact=prior_curvature=exact_data=residual_feedback=None
    else:
        hv=inherited_hvp(model,state,Y,step)
        vectors['HVP_actual_step']=hv
        gn=float(np.linalg.norm(Jr@step)**2/scale)
        exact=float(np.dot(step,hv))
        reg=model.reg
        delta=reg.D@state['m']
        smooth_q=np.sqrt(delta*delta+reg.epsilon**2)
        prior_curvature=float(np.sum((2*reg.smooth+reg.tv*reg.epsilon**2/smooth_q**3)/reg.edges*(reg.D@step)**2))
        exact_data=exact-prior_curvature
        residual_feedback=exact-gn-prior_curvature
    return vectors,dict(
        evaluation_state='checkpoint iteration500 m_before, before the last material readout; not stage endpoint m_after',
        gn_data_directional_curvature=gn,exact_Phi_directional_curvature=exact,
        prior_directional_curvature=prior_curvature,
        exact_data_directional_curvature=exact_data,
        residual_feedback_directional=residual_feedback,
        full_Phi_minus_data_GN_directional=exact-gn if step is not None else None,
        real_Gram_eigenvalue_min=float(er[0]),real_Gram_eigenvalue_max=float(er[-1]),
        complex_Gram_eigenvalue_min=float(ec[0]),complex_Gram_eigenvalue_max=float(ec[-1]),
        real_Gram_negative_eigenvalues=int(np.sum(er<0)),
        real_Gram_numerical_rank=int(np.sum(er>np.finfo(float).eps*len(er)*max(er[-1],0))),
        rank_semantics='ordinary floating precision report; no exact-nullspace claim',
        complex_J_shape=list(Jc.shape),real_J_shape=list(Jr.shape),
        packing='Jcomplex rows source then receiver; real J rows source then all receiver real rows then all receiver imaginary rows. This differs only by row permutation from model.pack; Gram and norms invariant.',
        data_whitening='1/||Y||; Omega_s is not used for physical Jacobian',
        HVP_scope='inherited exact full-Phi directional HVP, including exact prior Hessian action; no full physical Hessian',
        wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu)


def analyze_checkpoint(model,Y,checkpoint,beta,config,endpoint=False):
    wall,cpu=time.perf_counter(),time.process_time()
    before=dict(model.cost);cost=Counter()
    has_step='m_before' in checkpoint and 'm_after' in checkpoint
    m=np.asarray(checkpoint['m_before'] if has_step else checkpoint['m_state'],float)
    ma=np.asarray(checkpoint['m_after'],float) if has_step else m.copy()
    step=ma-m  # zero placeholder is removed from all fixed-state outputs below
    J=np.asarray(checkpoint['J_after_current']);Q=np.asarray(checkpoint['Q'])
    ws=1/np.sum(abs(model.E)**2,axis=0)
    Os=np.diag(ws);Ls=np.diag(np.sqrt(ws));rho=float(config.get('rho',1.))
    lo,hi=[np.broadcast_to(x,m.shape).copy() for x in config['material_box']]
    reg=lambda z:model.reg.evaluate(z)[:2]
    prior,gprior=reg(m)
    phi,gphysical=model.data_value_gradient(m,Y)
    state=model.state(m)
    residual_before=float(np.linalg.norm(state['Y']-Y)/np.linalg.norm(Y))
    blocks=OneHopAbsorption(model.G,model.shape,cost)
    b,d,a=blocks.parameters(m,beta)
    r=J-m[:,None]*Q
    u=blocks.scaled_rhs(r,beta,d)
    Lu=a[:,None]*blocks.off_action(u)
    crosscell=2*np.sum(np.real(u.conj()*Lu)*ws,axis=1)
    unorm=np.sqrt(np.sum(abs(u)**2*ws,axis=1))
    vnorm=np.sqrt(np.sum(abs(Lu)**2*ws,axis=1))
    phaseden=unorm*vnorm
    phasecos=np.divide(crosscell/2,phaseden,out=np.full_like(m,np.nan),where=phaseden>0)
    arrays=dict(m_before=m,m_after=ma,delta_m=step,r=r,u=u,Lu=Lu,
        cross_cell=crosscell,u_cell_amplitude=unorm,Lu_cell_amplitude=vnorm,
        phase_cosine_cell=phasecos,ws=ws,R_ordinary=beta*r/(1+beta*m)[:,None],
        R_self=u,R_plus=u+Lu,R_minus=u-Lu,g_physical=gphysical,
        g_physical_data=gphysical-gprior,g_prior=gprior)
    summary=dict(stage=int(checkpoint['stage']),iteration=int(checkpoint['iteration']),
        beta_material=float(beta),truth_access='NONE',actual_step_available=has_step,
        source_metric='fixed Omega_s=diag(1/sum_cells |Einc|^2), all source columns jointly share real material',
        physical_data_whitening='1/||Y|| in Phi; distinct from Omega_s and from noise sigma',
        defect_norm=float(np.linalg.norm(r*np.sqrt(ws))),
        cross_value=float(np.sum(crosscell)),
        cross_negative_cells=int(np.sum(crosscell<0)),
        u_norm=float(np.linalg.norm(u*np.sqrt(ws))),Lu_norm=float(np.linalg.norm(Lu*np.sqrt(ws))),
        global_phase_cosine=cosine(u*np.sqrt(ws),Lu*np.sqrt(ws)),
        prior_value=float(prior),frozen_LS={},gradients={})
    exact={};frozen={}
    for kind,sign in [('self',0),('plus',1),('minus',-1)]:
        cb=phase_blocks(blocks,sign)
        value,g=onehop_material_value_gradient(cb,J,Q,m,beta,Os,rho,reg,cost)
        R=u+sign*Lu
        w=R*ws
        adj=w+sign*blocks.off_action(a.conj()[:,None]*w,adjoint=True)
        ptr=(beta/d).conj()[:,None]*adj
        gf=-rho*np.sum(np.real(Q.conj()*ptr),axis=1)
        cost['exact_state_gradient_diagnostic_evaluations']+=1
        exact[kind]=g-gprior;frozen[kind]=gf
        arrays['g_exact_state_'+kind]=exact[kind]
        arrays['g_frozen_'+kind]=gf
        arrays['g_connection_'+kind]=exact[kind]-gf
        summary['gradients'][kind]=dict(state_value=float(value-prior),
            exact_state=vector_summary(exact[kind],step,-gphysical),
            frozen=vector_summary(gf,step,-gphysical),
            connection=vector_summary(exact[kind]-gf,step,-gphysical),
            negative_exact_cosine_step=cosine(-exact[kind],step),
            negative_exact_cosine_negative_physical=cosine(-exact[kind],-gphysical))
    denominator=1+beta*m
    ordinary_exact=-rho*beta**2*np.sum(np.real((Q+beta*J)*r.conj())*ws,axis=1)/denominator**3
    ordinary_frozen=-rho*beta**2*np.sum(np.real(Q*r.conj())*ws,axis=1)/denominator**2
    exact['ordinary']=ordinary_exact;frozen['ordinary']=ordinary_frozen
    arrays['g_exact_state_ordinary']=ordinary_exact
    arrays['g_frozen_ordinary']=ordinary_frozen
    arrays['g_connection_ordinary']=ordinary_exact-ordinary_frozen
    summary['gradients']['ordinary']=dict(
        state_value=.5*rho*float(np.linalg.norm(arrays['R_ordinary']*np.sqrt(ws))**2),
        exact_state=vector_summary(ordinary_exact,step,-gphysical),
        frozen=vector_summary(ordinary_frozen,step,-gphysical),
        connection=vector_summary(ordinary_exact-ordinary_frozen,step,-gphysical))
    summary['state_objectives']={kind:.5*rho*float(np.linalg.norm(arrays['R_'+kind]*np.sqrt(ws))**2) for kind in ('ordinary','self','plus','minus')}
    for kind in ('ordinary','self','plus','minus'):
        connection=exact[kind]-frozen[kind]
        negative=connection<0
        positive_step=step>0
        summary['gradients'][kind]['connection'].update(
            negative_fraction=float(np.mean(negative)),
            cosine_uniform_material_increase=cosine(-connection,np.ones_like(m)),
            negative_component_at_lower_active=int(np.sum(negative&(m<=lo))),
            negative_component_at_upper_active=int(np.sum(negative&(m>=hi))),
            negative_component_and_positive_actual_step=int(np.sum(negative&positive_step)),
            positive_actual_step_fraction_with_negative_component=float(np.mean(negative[positive_step])) if np.any(positive_step) else None,
            negative_component_fraction_with_positive_actual_step=float(np.mean(positive_step[negative])) if np.any(negative) else None,
            negative_component_actual_step_projection=float(np.dot(-np.minimum(connection,0),step)),
            b_real_min=float(np.min(b.real)))
    gcross=(exact['plus']-exact['minus'])/2
    gquad=(exact['plus']+exact['minus'])/2-exact['self']
    arrays.update(g_cross=gcross,g_quad=gquad)
    for kind,g in [('cross',gcross),('quad',gquad),('prior',gprior)]:
        summary['gradients'][kind]=vector_summary(g,step,-gphysical)
        summary['gradients'][kind]['negative_gradient_cosine_step']=cosine(-g,step)
        summary['gradients'][kind]['negative_gradient_cosine_negative_physical']=cosine(-g,-gphysical)
    for kind in ('ordinary','self','plus','minus'):
        ls=frozen_ls(blocks,m,Q,r,beta,ws,kind,cost)
        for key,value in ls.items(): arrays['LS_'+kind+'_'+key]=value
        finite=ls['gram_condition'][np.isfinite(ls['gram_condition'])]
        summary['frozen_LS'][kind]=dict(
            d=vector_summary(ls['d'],step,-gphysical),
            d0=vector_summary(ls['d0'],step,-gphysical),
            transport=vector_summary(ls['transport'],step,-gphysical),
            transport_rhs=vector_summary(ls['transport_rhs'],step,-gphysical),
            Gram_eigenvalue_min=float(ls['gram_eigenvalues'].min()),
            Gram_condition_max_finite=float(finite.max()) if finite.size else None,
            Gram_infinite_condition_blocks=int(np.sum(~np.isfinite(ls['gram_condition']))),
            fallback_blocks=ls['fallback_blocks'].tolist(),fallback_ranks=ls['fallback_ranks'].tolist(),
            kernel_claim='NONE; weak numerical columns and ordinary LAPACK fallback recorded')
        if kind in ('ordinary','self'):
            summary['frozen_LS'][kind]['transport_semantics']='ANALYTIC_ZERO: nonzero per-cell scalar metric cancels; numerical subtraction and projection roundoff retained under raw_roundoff arrays'
            for label in ('transport','transport_rhs'):
                desc=summary['frozen_LS'][kind][label]
                desc.update(norm=0.,cosine_actual_step=None,cosine_negative_physical_gradient=None,
                    dot_actual_step=None,dot_negative_physical_gradient=None)
                desc['raw_roundoff_norm']=float(np.linalg.norm(ls[label+'_raw_roundoff']))

    pg=m-np.clip(m-gphysical,lo,hi)
    summary['physical']=dict(Phi_before=float(phi),relative_residual_before=residual_before,
        stepnorm=float(np.linalg.norm(step)),negative_gradient_cosine_step=cosine(-gphysical,step),
        linear_prediction=float(np.dot(gphysical,step)),box_projected_gradient_inf=float(np.max(abs(pg))),
        active_lower_before=int(np.sum(m<=lo)),active_upper_before=int(np.sum(m>=hi)),
        active_lower_after=int(np.sum(ma<=lo)),active_upper_after=int(np.sum(ma>=hi)))
    arrays['box_projected_gradient']=pg
    if endpoint:
        arr,geometry=endpoint_geometry(model,state,Y,step if has_step else None,cost)
        arrays.update(arr);summary['endpoint_geometry']=geometry
    state_after=model.state(ma)
    phia=float(.5*np.linalg.norm(state_after['Y']-Y)**2/np.linalg.norm(Y)**2+reg(ma)[0])
    actual=phia-phi
    summary['physical'].update(Phi_after=phia,
        relative_residual_after=float(np.linalg.norm(state_after['Y']-Y)/np.linalg.norm(Y)),
        actual_deltaPhi=float(actual),prediction_error=float(actual-np.dot(gphysical,step)))
    if endpoint:
        cf={}
        for kind in ('ordinary','self','plus','minus'):
            c0=dict(cost);p0=dict(model.cost);tw,tc=time.perf_counter(),time.process_time()
            if kind=='ordinary':
                chosen,info=_baseline_readout(J@Ls,Q@Ls,m,beta,np.ones(model.ns),rho,
                    (lo,hi),reg,True,config,cost)
                cost['counterfactual_source_whitening_products']+=2
                cost['counterfactual_ordinary_LBFGS_iterations']+=info['material_inner_nit']
                cost['counterfactual_ordinary_readouts']+=1
            else:
                cb=phase_blocks(blocks,{'self':0,'plus':1,'minus':-1}[kind])
                chosen,info=onehop_material_readout(cb,J,Q,m,beta,Os,Ls,rho,(lo,hi),reg,config,cost)
                cost['counterfactual_onehop_readouts']+=1
            st=model.state(chosen)
            rr=float(np.linalg.norm(st['Y']-Y)/np.linalg.norm(Y))
            score=float(.5*rr**2+reg(chosen)[0])
            arrays['counterfactual_m_'+kind]=chosen
            arrays['counterfactual_delta_m_'+kind]=chosen-m
            cf[kind]=dict(Phi=score,deltaPhi=score-phi,relative_residual=rr,
                step=vector_summary(chosen-m,step,-gphysical),inner=info,
                original_forward_evaluation=True,
                physical_model_cost=difference_counter(model.cost,p0),
                auxiliary_operations=difference_counter(cost,c0),
                wall_seconds=time.perf_counter()-tw,cpu_seconds=time.process_time()-tc)
        summary['counterfactual_readouts']=cf
        summary['counterfactual_scope']='same fixed material,J,Q,beta,prior,box and public maximum of 60 steps with each original compatible readout and deterministic seed policy; ordinary uses mature scalar CIE readout and is not claimed to solve the same nonlinear subproblem as the onehop arms; true physical G unchanged; self is onehop formula with zero Goff, not historical blockself byte replay; no GT selection'
    physical_delta=difference_counter(model.cost,before)
    # Frozen physics receipts count LU/RHS, but do not explicitly count the
    # already executed dense Green products. Derive these from the unchanged
    # call paths; no extra matrix action is performed for this accounting.
    cost['derived_physical_Green_forward_columns']=int(physical_delta.get('LU',0))*model.ns+(model.ns if endpoint and has_step else 0)
    cost['derived_physical_Green_adjoint_columns']=model.ns if endpoint and has_step else 0
    cost['derived_physical_Green_matrix_apply_calls']=int(physical_delta.get('LU',0))+(2 if endpoint and has_step else 0)
    summary['cost']=dict(auxiliary_operations=dict(cost),physical_model_cost=physical_delta,
        wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,
        money=None,scope='diagnostic only; original replay charged separately; Green RHS, small Gram solves, physical LU/RHS distinct; no formula identity tests')
    if not has_step:
        summary['state_scope']='Historical stage endpoint m/J/Q only; actual previous material, readout step and PR trajectory NOT_AVAILABLE. No stage-to-stage difference is used as a single step.'
        for key in ('m_before','m_after','delta_m','box_projected_gradient'):
            if key in arrays: del arrays[key]
        arrays['m_state']=m
        summary['physical'].update(Phi_after=None,relative_residual_after=None,
            stepnorm=None,negative_gradient_cosine_step=None,linear_prediction=None,
            actual_deltaPhi=None,prediction_error=None,active_lower_after=None,active_upper_after=None)
        if endpoint:
            summary['endpoint_geometry']['HVP_scope']='NOT_AVAILABLE: historical actual step not saved; no directional HVP executed'
        def clear_step_links(value):
            if isinstance(value,dict):
                for key in value:
                    if key in ('cosine_actual_step','dot_actual_step','negative_exact_cosine_step','negative_gradient_cosine_step','negative_component_and_positive_actual_step','positive_actual_step_fraction_with_negative_component','negative_component_fraction_with_positive_actual_step','negative_component_actual_step_projection'):
                        value[key]=None
                    else: clear_step_links(value[key])
        clear_step_links(summary)
    return summary,arrays
