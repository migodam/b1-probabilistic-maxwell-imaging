"""Measurement-only bounded pair construction; OFFLINE scoring comes last.

These are feasible local compensation witnesses, never global nonlinear minima.
No measurement-dependent likelihood weights are called Bayesian posteriors.
"""
from pathlib import Path
from time import perf_counter
import json
import numpy as np
from scipy import ndimage
from scipy.optimize import minimize,LinearConstraint
from scipy.stats import norm
from .data import release_root,safe,case_metadata,offline_targets,train_deltas
from .physics import FrozenMaxwellAdapter
from .budget import Book,physical_totals
from .hull_audit import csv_write,topology_metrics
from .frozen_a26_metrics import image_metrics
from .capacity import serial


def pack(y):
    y=np.asarray(y);return np.concatenate((y.real.T,y.imag.T),axis=1).reshape(-1)


def real_material(z):return np.r_[np.asarray(z).real,np.asarray(z).imag]


def chart(image,background,protected,quantiles,contrast_radius=.30,boundary_radius=.35):
    """8 contrast + 4 continuous edge-displacement coordinates, per endpoint.

    Regions are fixed from the online estimate. Both parts of each material
    contrast may change independently. The proposed cavity interior is held
    at the background, so compensation cannot erase the nominal edit.
    """
    delta=(image-background).reshape(24,24);mag=np.abs(delta);maxmag=float(mag.max())
    q=np.asarray(quantiles)*maxmag
    labels=np.digitize(mag,q);occupied=mag>.15*maxmag
    allowed=~np.asarray(protected).reshape(24,24)
    regions=[(labels==i)&occupied&allowed for i in range(4)]
    contrast=[]
    for region in regions:contrast.extend([contrast_radius*region.ravel(),1j*contrast_radius*region.ravel()])
    # First-order cell boundary displacement, bounded to .35 pixel per axis.
    gy,gx=np.gradient(delta);yy,xx=np.mgrid[-1:1:24j,-1:1:24j]
    boundary=[gx,gy,gx*yy,gy*xx]
    boundary=[boundary_radius*v.ravel()*allowed.ravel() for v in boundary]
    return np.asarray(contrast+boundary,np.complex128).T


def bounded_pair(a,b,ua,ub,delta0,operator):
    da=operator(ua);db=operator(ub)
    design=np.c_[np.column_stack([pack(x) for x in np.moveaxis(da,-1,0)]),
                  np.column_stack([pack(x) for x in np.moveaxis(db,-1,0)])]
    p=ua.shape[1];design[:,:p]*=-1
    target=pack(delta0);scale=max(float(target@target),1.)
    hessian=design.T@design/scale;g=design.T@target/scale
    aa=np.vstack((ua.real,ua.imag));bb=np.vstack((ub.real,ub.imag))
    constraints=np.zeros((4*576,2*p));constraints[:1152,:p]=aa;constraints[1152:,p:]=bb
    low=np.r_[np.full(576,-.5),np.zeros(576)];high=np.r_[np.full(576,4.),np.full(576,.5)]
    lower=np.r_[low-real_material(a),low-real_material(b)]
    upper=np.r_[high-real_material(a),high-real_material(b)]
    def objective(z):return float(.5*z@hessian@z+g@z)
    result=minimize(objective,np.zeros(2*p),jac=lambda z:hessian@z+g,
        method='SLSQP',bounds=[(-1.,1.)]*(2*p),constraints=LinearConstraint(constraints,lower,upper),
        options=dict(maxiter=100,ftol=1e-12))
    z=result.x.copy();displacement=constraints@z
    # Repair rounding by contracting parameters, without pixel clipping or changing charts.
    factor=1.
    for _ in range(4):
        if np.all(displacement>=lower) and np.all(displacement<=upper):break
        ratios=np.r_[upper[displacement>0]/displacement[displacement>0],
                     lower[displacement<0]/displacement[displacement<0]]
        factor=min(1.,float(np.min(ratios)))*(1-1e-12);z*=factor;displacement=constraints@z
    if not (np.all(displacement>=lower) and np.all(displacement<=upper)):
        raise ValueError('COMPENSATION_PHYSICAL_BOX_FAILED')
    ca=a+ua@z[:p];cb=b+ub@z[p:]
    residual=target+design@z
    unrestricted=np.linalg.lstsq(design,-target,rcond=None)[0]
    unbounded_norm=float(np.linalg.norm(target+design@unrestricted))
    return dict(a=ca,b=cb,parameters=z,linear_distance=float(np.linalg.norm(residual)),
        unrestricted_distance=unbounded_norm,
        unrestricted_out_of_bounds=bool(np.any(np.abs(unrestricted)>1) or
            np.any(constraints@unrestricted<lower) or np.any(constraints@unrestricted>upper)),
        status=str(result.message),converged=bool(result.success),iterations=int(result.nit),
        parameter_contraction=factor,design=design,scope='FEASIBLE_LOCAL_LINEAR_WITNESS')


def run(root):
    root=Path(root).resolve();r=release_root(root);out=root/'results/a30/compensation_corrected';out.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((root/'a30_config.json').read_text());protocol=cfg['E1']['implementation']
    # TRAIN alone fixes region thresholds; no development truth enters pair construction.
    train,_=train_deltas(root);nz=np.abs(train)[np.abs(train)>.15]
    quantiles=np.quantile(nz,[.25,.5,.75]);quantiles=(quantiles/max(float(nz.max()),1e-12)).tolist()
    metadata=case_metadata(root);byparent={}
    for row in metadata:byparent.setdefault(row['parent_id'],[]).append(row)
    order=[]
    for j in range(4):
        for p in (250065+j,250069+j,250077+j,250073+j):
            order.extend(sorted(byparent[p],key=lambda x:0 if x['background_id']=='weak' else 1))
    prepared=[]
    for meta in order:
        with np.load(safe(r,'results/a29/pools/'+meta['case_id']+'.npz'),allow_pickle=False) as z:
            ids=z['ids'].astype(str)
            chosen=[i for i,s in enumerate(ids) if s=='geometry_piecewise:cavity:n3:r2:s0:phase_n3']
            if not chosen:chosen=[i for i,s in enumerate(ids) if s=='geometry_piecewise:cavity:n3:r1:s0:phase_n3']
            if not chosen:
                prepared.append(dict(meta=meta,status='NO_REGISTERED_CAVITY_EDIT'));continue
            index=chosen[0];a=z['chi'][0].copy();b=z['chi'][index].copy();bg=z['background'].copy()
            # Protect actual removed material, not globally changed phase-projection pixels.
            protected=(np.abs(b-bg)<1e-9)&(np.abs(a-bg)>.15*float(np.max(np.abs(a-bg))))
            if not np.any(protected):
                prepared.append(dict(meta=meta,status='NO_REGISTERED_CAVITY_INTERIOR'));continue
            # Cavity directions are zero throughout the edited interior for both endpoints.
            ua=chart(a,bg,protected,quantiles,protocol['contrast_radius'],protocol['boundary_pixels'])
            ub=chart(b,bg,protected,quantiles,protocol['contrast_radius'],protocol['boundary_pixels'])
            prepared.append(dict(meta=meta,status='READY',a=a,b=b,bg=bg,ua=ua,ub=ub,
                sigma=float(z['sigma']),candidate_id=str(ids[index]),protected=protected))
    (out/'PAIR_CHOICES.json').write_text(json.dumps(serial(dict(region_quantiles=quantiles,
        target_access=False,pairs=[dict(**x['meta'],status=x['status'],candidate_id=x.get('candidate_id')) for x in prepared])),indent=2)+'\n')
    import csv
    prior_path=root/'results/a30/compensation/COMPENSATION_PAIRS.csv'
    previous={}
    if prior_path.exists():
        with prior_path.open() as f:
            previous={row['case_id']:row for row in csv.DictReader(f) if row.get('valid')=='True'}
    book=Book(root,stage='e1');adapters={};gnorm={};snorm={};records=[];raw={};failed_certificates=0
    legal_complex_total=sum(x['status']=='READY' and x['meta']['family']!='smooth' for x in prepared)
    stopped=False
    for i,item in enumerate(prepared):
        meta=item['meta'];case=meta['case_id'];bgid=meta['background_id']
        if item['status']!='READY' or stopped:
            records.append(dict(**meta,status=item['status'] if not stopped else 'NOT_RUN_CERTIFICATE_VALUE_STOP',valid=False));continue
        tic=perf_counter()
        try:
            if bgid not in adapters:
                allowed=[safe(r,'inputs/a28/results/a28/exact_states/'+x['meta']['case_id']+'.npz') for x in prepared if x['meta']['background_id']==bgid]
                adapters[bgid]=FrozenMaxwellAdapter.from_background_npz(safe(r,'inputs/a25/background/'+bgid+'.npz'),book,allowed_cache_paths=allowed)
                gnorm[bgid]=float(np.linalg.norm(adapters[bgid].G,2));snorm[bgid]=float(np.linalg.norm(adapters[bgid].S,2))
            model=adapters[bgid];a,b=item['a'],item['b'];sigma=item['sigma'];white=2**.5/sigma
            anchor=model.cached_anchor(a,safe(r,'inputs/a28/results/a28/exact_states/'+case+'.npz'))
            transfer=model.receiver_transfer(anchor)
            def opF(u):
                return np.einsum('qn,nk,ns->qsk',transfer,u,anchor.field)*white
            def opL(u):
                return np.einsum('qn,nk,ns->qsk',model.Db,u,model.eb)*white
            la=model.L_response(a);lb=model.L_response(b)
            nominal_L=float(np.linalg.norm((lb-la)*white))
            linF=opF((b-a)[:,None])[:,:,0];linL=(lb-la)*white
            if case in previous:
                nominal_F=float(previous[case]['nominal_F'])
                nominal_remainder=float(previous[case]['nominal_F_tangent_remainder'])
                nominal_backward=float(previous[case]['core_backward_max']);new_nominal=0
            else:
                nominal=model.full_state(b,purpose=case+':nominal_cavity')
                nominal_F=float(np.linalg.norm((nominal.data-anchor.data)*white))
                nominal_remainder=float(np.linalg.norm((nominal.data-anchor.data)*white-linF))
                nominal_backward=nominal.residual;new_nominal=1
            ua,ub=item['ua'],item['ub']
            contrastF=bounded_pair(a,b,ua[:,:8],ub[:,:8],linF,opF)
            contrastL=bounded_pair(a,b,ua[:,:8],ub[:,:8],linL,opL)
            fullF=bounded_pair(a,b,ua,ub,linF,opF)
            fullL=bounded_pair(a,b,ua,ub,linL,opL)
            # Two compensated endpoints, plus nominal b: 3 fresh states per condition.
            # Both nominal/compensated a changes have to be paid; F endpoints total 2.
            fa=model.full_state(fullF['a'],purpose=case+':compensated_a')
            fb=model.full_state(fullF['b'],purpose=case+':compensated_b')
            distance_F=float(np.linalg.norm((fb.data-fa.data)*white))
            distance_L_on_F=float(np.linalg.norm((model.L_response(fullF['b'])-model.L_response(fullF['a']))*white))
            # A fixed first parent per family provides matched L-witness full-F checks.
            # The remaining L distances stay explicitly linear-only under the state cap.
            if meta['parent_id'] in protocol['L_witness_full_F_parent_ids']:
                lfa=model.full_state(fullL['a'],purpose=case+':L_compensated_a')
                lfb=model.full_state(fullL['b'],purpose=case+':L_compensated_b')
                distance_F_on_L=float(np.linalg.norm((lfb.data-lfa.data)*white));new_L=2
                L_backward=max(lfa.residual,lfb.residual)
            else:distance_F_on_L=None;new_L=0;L_backward=0.
            # One conservative, uniform fixed-chart bound. No extra solve/probe is used.
            # All points in the coordinate box satisfy |h|<=h_max around cached anchor.
            hmax_a=np.sum(np.abs(ua),axis=1)
            hmax_b=np.abs(b-a)+np.sum(np.abs(ub),axis=1)
            invnorm=1/float(np.linalg.svd(anchor.matrix,compute_uv=False)[-1])
            bounds=[];qs=[]
            for hmax in (hmax_a,hmax_b):
                q=invnorm*float(hmax.max())*gnorm[bgid];qs.append(q)
                bounds.append(white*snorm[bgid]*invnorm*q/(1-q)*
                    float(np.linalg.norm(hmax[:,None]*anchor.field)) if q<1 else None)
            uniform_b=sum(bounds) if all(x is not None for x in bounds) else None
            lower=max(0.,fullF['linear_distance']-uniform_b) if uniform_b is not None else None
            upper=fullF['linear_distance']+uniform_b if uniform_b is not None else None
            point_remainder=float(np.linalg.norm((fb.data-fa.data)*white-
                opF((fullF['b']-fullF['a'])[:,None])[:,:,0]))
            bound_valid=(uniform_b is not None and point_remainder<=uniform_b+1e-9)
            useful=bound_valid and (lower>=3 or upper<=1)
            if uniform_b is not None and not bound_valid:raise ValueError('FALSE_UNIFORM_REMAINDER_BOUND')
            if meta['family']!='smooth' and not useful:failed_certificates+=1
            rng=np.random.default_rng(protocol['independent_detection_seed']+i)
            # Fixed pair, independent test-noise draws. Not end-to-end candidate selection.
            draws=rng.standard_normal((1000,2))
            empirical_nominal=float(np.mean(np.r_[draws[:500,0]>nominal_F/2,draws[500:,0]<-nominal_F/2]))
            empirical_comp=float(np.mean(np.r_[draws[:500,1]>distance_F/2,draws[500:,1]<-distance_F/2]))
            rec=dict(**meta,status='COMPLETE',valid=True,candidate_id=item['candidate_id'],
                nominal_F=nominal_F,nominal_L=nominal_L,calibrated_L_pair_distance=nominal_L,
                nominal_slope_F=float(np.linalg.norm(linF)),nominal_F_tangent_remainder=nominal_remainder,
                contrast_only_F_linear=contrastF['linear_distance'],contrast_only_L_linear=contrastL['linear_distance'],
                geometry_contrast_F_linear=fullF['linear_distance'],geometry_contrast_L_linear=fullL['linear_distance'],
                compensated_F=distance_F,compensated_L=fullL['linear_distance'],compensated_F_at_L_witness=distance_F_on_L,
                compensated_L_at_F_witness=distance_L_on_F,unbounded_F=fullF['unrestricted_distance'],
                unbounded_L=fullL['unrestricted_distance'],unbounded_F_outside_legal=fullF['unrestricted_out_of_bounds'],
                uniform_q_a=qs[0],uniform_q_b=qs[1],uniform_b=uniform_b,
                certified_lower=lower,certified_upper=upper,endpoint_remainder=point_remainder,
                uniform_bound_valid=bound_valid,certificate_decision_useful=useful,
                material_pair_rms=float(np.linalg.norm(fullF['a']-fullF['b'])/24),
                material_nominal_rms=float(np.linalg.norm(a-b)/24),
                nominal_detection_error=float(norm.cdf(-nominal_F/2)),compensated_detection_error=float(norm.cdf(-distance_F/2)),
                independent_empirical_nominal_error=empirical_nominal,independent_empirical_comp_error=empirical_comp,
                convex_solver_converged=fullF['converged'] and fullL['converged'],
                solve_status_F=fullF['status'],solve_status_L=fullL['status'],
                core_backward_max=max(anchor.residual,nominal_backward,fa.residual,fb.residual,L_backward),
                protected_cavity_pixels=int(item['protected'].sum()),
                contrast_boundary_matrix_rank=int(np.linalg.matrix_rank(np.c_[item['ua'],item['ub']])),
                wall_seconds=perf_counter()-tic,actual_fresh_states=2+new_nominal+new_L,
                L_witness_full_F_verified=bool(new_L),nominal_F_reused_from_paid_run=not bool(new_nominal),
                scope='ONLINE_CONSTRUCTED_PAIRS_MECHANISM_DIAGNOSTIC_NOT_RECONSTRUCTION')
            records.append(rec)
            for key,image in [('nominal_a',a),('nominal_b',b),('compensated_a',fullF['a']),('compensated_b',fullF['b']),
                              ('L_compensated_a',fullL['a']),('L_compensated_b',fullL['b'])]:raw[case+'__'+key]=image
            raw[case+'__F_parameters']=fullF['parameters'];raw[case+'__L_parameters']=fullL['parameters']
            np.savez_compressed(out/'COMPENSATION_ENDPOINTS.npz',**raw)
            print(json.dumps(dict(case=case,nominal_F=nominal_F,compensated_F=distance_F,
                certificate_useful=bool(useful),fresh_states=physical_totals(root)['material_state_evaluations'])),flush=True)
            # Aggregate original gate is evaluated after this corrected fixed cohort.
            # There is no further physical expansion or certificate campaign.
        except Exception as error:
            failure=dict(**meta,status='FAILED',valid=False,error=type(error).__name__+': '+str(error))
            records.append(failure)
            with (root/'results/a30/FAILURE_LEDGER.jsonl').open('a') as f:f.write(json.dumps(failure)+'\n')
            stopped=True
        csv_write(out/'COMPENSATION_PAIRS.csv',records)
    csv_write(out/'COMPENSATION_PAIRS.csv',records)
    # Only now read development truth for material/edge/topology consequence.
    targets,tbg,_=offline_targets(root,[x['case_id'] for x in metadata]);lookup={x['case_id']:i for i,x in enumerate(metadata)}
    offline=[]
    for rec in records:
        if not rec['valid']:continue
        i=lookup[rec['case_id']];case=rec['case_id']
        for key in ('nominal_a','nominal_b','compensated_a','compensated_b','L_compensated_a','L_compensated_b'):
            image=raw[case+'__'+key]
            offline.append(dict(case_id=case,parent_id=rec['parent_id'],family=rec['family'],endpoint=key,
                **image_metrics(image,targets[i],tbg[i]),**topology_metrics(image,targets[i],tbg[i])))
    csv_write(out/'OFFLINE_ENDPOINT_MATERIAL.csv',offline)
    no_value=failed_certificates>legal_complex_total/2
    summary=dict(status='STOP_AFTER_CORRECTED_COHORT_CERTIFICATE_NO_PRACTICAL_VALUE' if no_value else ('FAILED' if stopped else 'COMPLETE'),
        conditions_complete=sum(x['valid'] for x in records),conditions_total=32,
        legal_complex_conditions=legal_complex_total,useful_complex_certificates=sum(x.get('certificate_decision_useful',False) for x in records if x['family']!='smooth'),
        physical_cost=physical_totals(root),full_F_global_minimum='NOT_ESTABLISHED',
        witnesses='FEASIBLE_LOCAL_ONLY',E2='NOT_RUN_USER_SCOPE',
        calibrated_L_note='The single anchor offset cancels exactly in every candidate-pair secant distance.',
        nuisance_note='Known background retained. Additional unknown calibration nuisance NOT_RUN; no nuisance difficulty attributed to the original known-background task.')
    (out/'E1_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path('.'));run(p.parse_args().root)
