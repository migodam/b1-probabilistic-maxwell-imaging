"""Exact absorption of estimated wall cells, followed by genuine reopening."""
from __future__ import annotations
import numpy as np


def relative(a,b):return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),1e-14))


def absorb_wall(reference,z,currents,fields,fraction=.5,tolerance=1e-8,selected=None):
    from .currents import CurrentGreen
    model=reference.model;nw=model.nwall
    # Estimated material alone selects absorption, independently of truth/labels.
    strength=np.hypot(z[:nw]-1,z[nw:2*nw])
    number=max(1,int(np.ceil(fraction*nw)))
    selected=np.argsort(-strength,kind='stable')[:number] if selected is None else np.asarray(selected,int)
    if selected.ndim!=1 or not len(selected) or np.any((selected<0)|(selected>=nw)):raise ValueError('nonempty legal selected wall cells')
    frozen=np.zeros(model.nparam,bool);frozen[selected]=True;frozen[nw+selected]=True
    greens=[];js=[];efs=[];checks=[]
    for ii,(g0,jfull,field) in enumerate(zip(reference.greens,currents,fields)):
        mask=g0.wall&np.isin(g0.cells,selected)
        chiA=np.where(mask,g0.chi(z),0)
        # At most three persistent factors: release original factor before LA.
        # Original GS0 stays fixed and is used in the transformation check.
        g0.factor=None
        ga=CurrentGreen(model,ii,absorbed=chiA,reference='absorbed_wall')
        ga.absorbed_mask=mask.copy()
        jr=(g0.chi(z)-chiA)[:,None]*field
        ea=ga.inc+ga.apply_GD(jr)
        ja=jr+chiA[:,None]*ea
        original_data=g0.apply_GS(ja)
        effective_data=ga.dA+ga.apply_GS(jr)
        check=dict(frequency=g0.p.freq,field_relative=relative(field,ea),
                   current_relative=relative(jfull,ja),
                   original_GS0_vs_absorbed_data_relative=relative(original_data,effective_data),
                   removed_quadrature_points=int(mask.sum()),factor_bytes=ga.factor_bytes)
        if max(check[k] for k in ('field_relative','current_relative','original_GS0_vs_absorbed_data_relative'))>tolerance:
            raise RuntimeError('exact absorption transformation failed '+str(check))
        greens.append(ga);js.append(jr);efs.append(ea);checks.append(check)
    event=dict(kind='EXACT_ABSORB',selected_estimated_wall_cells=selected.tolist(),
               material_variables_removed=int(frozen.sum()),selection='caller-selected estimated wall cells',
               original_GS0_retained=True,background_data_retained=True,effective_illumination_retained=True,
               checks=checks,tolerance=tolerance)
    model.book.event(model.phase,'CSI exact absorption',event_kind=event['kind'],
                     **{key:value for key,value in event.items() if key!='kind'})
    return greens,js,efs,frozen,event


def reopen(reference,absorbed,currents,fields,z,tolerance=1e-8):
    js=[];efs=[];checks=[]
    for g0,ga,jr,field in zip(reference.greens,absorbed,currents,fields):
        # Exact physical current transformation, including absorbed-wall current.
        jfull=jr+ga.absorbed[:,None]*field
        before=ga.dA+ga.apply_GS(jr);after=g0.apply_GS(jfull)
        ga.factor=None
        original_field=g0.inc+g0.apply_GD(jfull)
        check=dict(frequency=g0.p.freq,data_relative=relative(before,after),
                   field_relative=relative(field,original_field))
        if max(check['data_relative'],check['field_relative'])>tolerance:
            raise RuntimeError('exact REOPEN transformation failed '+str(check))
        js.append(jfull);efs.append(original_field);checks.append(check)
    event=dict(kind='REOPEN',all_material_parameters_active=True,
               material_variables_restored=2*len(np.unique(np.concatenate([g.cells[g.absorbed_mask] for g in absorbed]))),
               original_fixed_reference_restored=True,checks=checks,tolerance=tolerance)
    reference.model.book.event(reference.model.phase,'CSI material REOPEN',event_kind=event['kind'],
                               **{key:value for key,value in event.items() if key!='kind'})
    return js,efs,event


def risk_probe(model,z,y,prior,reference=None):
    """Actual full-FEM projected objective gradient; every RHS is metered."""
    from .fwi import regularization
    dv,g,pred=(model.gradient(z,y) if reference is None else reference.gradient(z,y))
    pv,pg,parts=regularization(z,model.nwall,model.target_shape,model.wall_shape,prior)
    gradient=g+pg;lo,hi,scale=model.bounds_and_scale()
    projected=gradient.copy()
    projected[((z<=lo+1e-10)&(gradient>0))|((z>=hi-1e-10)&(gradient<0))]=0
    return dict(objective=dv+pv,data_objective=dv,prior=pv,
                scaled_projected_gradient=projected*scale,diagnostics=pred.diagnostics)


def risk_selection(reference,z,fields,history,probe,prior,options):
    """Frozen empirical surrogates, explicitly not posterior uncertainty/bias bounds."""
    model=reference.model;nw=model.nwall;_,_,scale=model.bounds_and_scale()
    recent=np.asarray(history[-(options.risk_stability_window+1):])
    stability=np.full(nw,np.inf)
    if len(recent)>=options.risk_stability_window+1:
        changes=np.max(np.abs(np.diff(recent,axis=0))/scale[None,:],axis=0)
        stability=np.maximum(changes[:nw],changes[nw:2*nw])
    # Conditional state least-squares diagonal plus the exact diagonal wall prior.
    neighbors=np.zeros(model.wall_shape)
    neighbors[:,:-1]+=1;neighbors[:,1:]+=1;neighbors[:-1,:]+=1;neighbors[1:,:]+=1
    re_curvature=(1+prior.wall_spatial_weight*neighbors.ravel())/prior.wall_eps_sigma**2
    loss_curvature=(1+prior.wall_spatial_weight*neighbors.ravel())/prior.wall_loss_sigma**2
    for g,field in zip(reference.greens,fields):
        density=options.state_weight*g.weights*np.sum(np.abs(field)**2,axis=1)
        np.add.at(re_curvature,g.cells[g.wall],density[g.wall])
        np.add.at(loss_curvature,g.cells[g.wall],density[g.wall]/g.p.freq**2)
    re_sd=1/np.sqrt(re_curvature);loss_sd=1/np.sqrt(loss_curvature)
    uncertainty=np.maximum(re_sd/scale[:nw],loss_sd/scale[nw:2*nw])
    pg=np.abs(probe['scaled_projected_gradient']);wall_pg=np.maximum(pg[:nw],pg[nw:2*nw])
    # Public GS0 action on a frozen-field wall-current tangent. This omits
    # self-consistent feedback and is only a data-impact screening surrogate.
    impact_squared=np.zeros(nw)
    for g,field in zip(reference.greens,fields):
        for cell in range(nw):
            mask=g.wall&(g.cells==cell)
            tangent=np.asarray(g.GS[:,mask]@field[mask]).T
            impact_squared[cell]+=float(np.sum(np.abs(tangent)**2))*(re_sd[cell]**2+(loss_sd[cell]/g.p.freq)**2)/model.sigma**2
    impact=np.sqrt(impact_squared)
    eligible=(stability<=options.risk_max_material_step)&(uncertainty<=options.risk_max_conditional_sd)&(wall_pg<=options.risk_max_projected_gradient)&(impact<=options.risk_max_data_impact)
    selected=np.flatnonzero(eligible)
    return selected,dict(kind='RISK_SCREEN',selected_estimated_wall_cells=selected.tolist(),
        stability_scaled=stability.tolist() if np.isfinite(stability).all() else None,
        conditional_sd_scaled=uncertainty.tolist(),scaled_projected_gradient=wall_pg.tolist(),
        conditional_sd_real=re_sd.tolist(),conditional_sd_loss=loss_sd.tolist(),
        frozen_field_data_impact_noise_units=impact.tolist(),
        conditional_uncertainty_scope='state-LS conditional curvature plus wall prior; not posterior uncertainty',
        impact_scope='public original GS0 frozen-field wall-current tangent; excludes self-consistent feedback and is not target-bias control',
        decision='PRELIMINARY_ELIGIBLE' if len(selected) else 'KEEP_ACTIVE')
