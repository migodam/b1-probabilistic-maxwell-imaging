"""Fixed observation000 public two-branch probability handoff, zero new PDE.

No work runs at import. main reads completed observation-only inference and
point validation plus paid PUBLIC bank/prior materials. No truth file is read,
no Family.point/Book is used. Branches are prior events, not a claim of modes.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from .probability_family import Family

ROOT=Path(__file__).resolve().parents[2]
OBSERVATION=0
FAMILY_ID='W_FIXED_SUPPORT_JOINT_V1'


def _encode(v):
    if isinstance(v,Path):return str(v)
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    raise TypeError(type(v).__name__)


def _write(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(obj,default=_encode,indent=2,allow_nan=False)+'\n')
    temporary.replace(path)


def weighted_finite_summary(nodes,weights):
    """Conditional moments and finite-atom quantiles; no continuous quantile."""
    x=np.asarray(nodes,float);w=np.asarray(weights,float)
    if x.ndim!=1 or w.shape!=x.shape or not np.isfinite(x).all() or not np.isfinite(w).all() or np.any(w<0) or not np.sum(w)>0:
        raise ValueError('Finite1D nodes and positive nonnegative weights required')
    w=w/np.sum(w);mean=float(w@x);sd=float(np.sqrt(w@((x-mean)**2)))
    order=np.argsort(x,kind='stable');xs,ws=x[order],w[order];cdf=np.cumsum(ws);cdf[-1]=1.
    q=xs[np.searchsorted(cdf,[.05,.5,.95],side='left')]
    return dict(mean=mean,SD=sd,q05=float(q[0]),q50=float(q[1]),q95=float(q[2]),
        semantics='FINITE_QUADRATURE weighted atomic quantiles; not exact continuous quantiles')


def joint_factors(weights,grid_shape,Fapprox):
    """Preserve gamma/material dependence: E[gamma F], not E[gamma]E[F]."""
    shape=tuple(int(n) for n in grid_shape)
    if len(shape)!=5 or shape[0]!=2 or shape[3]!=shape[4]:raise ValueError('Expected(branch,u,c,a,phase)grid')
    m=int(np.prod(shape[:3]));k=shape[3]*shape[4]
    w=np.asarray(weights,float)
    if w.size!=m*k or not np.isfinite(w).all() or np.any(w<0) or abs(float(w.sum())-1)>1e-10:
        raise ValueError('Invalid normalized stored finite joint weights')
    F=np.asarray(Fapprox)
    if F.shape!=(m,5,4,12) or not np.iscomplexobj(F) or not np.isfinite(F).all():raise ValueError('Invalid material-response approximation')
    return w.reshape(m,k)


def dependent_prediction(joint_weights,Fapprox,amplitude_axis,phase_axis):
    a,p=np.meshgrid(amplitude_axis,phase_axis,indexing='ij');gamma=(a*np.exp(1j*p)).ravel()
    w=np.asarray(joint_weights,float);F=np.asarray(Fapprox)
    if w.ndim!=2 or w.shape[1]!=len(gamma) or w.shape[0]!=len(F):raise ValueError('Incompatible joint factors')
    material_masses=w.sum(axis=1)
    gamma_weighted_material_factors=(w@gamma.real)+1j*(w@gamma.imag)
    predictive=np.tensordot(gamma_weighted_material_factors,F,axes=(0,0))
    # Diagnostic only: equality is not assumed and the product is never used.
    factorized=(w.sum(axis=0)@gamma)*np.tensordot(material_masses,F,axes=(0,0))
    return dict(material_masses=material_masses,gamma_nodes=gamma,
        gamma_weighted_material_factors=gamma_weighted_material_factors,
        posterior_predictive_approx=predictive,
        factorized_predictive_difference_norm=float(np.linalg.norm(predictive-factorized)))


def _odds(p,other,qualified):
    if not qualified:return dict(value=None,log_value=None,status='UNQUALIFIED')
    if p<=0:return dict(value=0.,log_value=None,status='ZERO_ESTIMATED_EVENT_MASS_NEGATIVE_INFINITY_LOG_LIMIT')
    if other<=0:return dict(value=None,log_value=None,status='POSITIVE_INFINITY_ZERO_ESTIMATED_COMPLEMENT_LIMIT')
    logodds=float(np.log(p)-np.log(other))
    value=float(np.exp(logodds)) if logodds<np.log(np.finfo(float).max) else None
    return dict(value=value,log_value=logodds,status='FINITE_APPROXIMATION' if value is not None else 'OVERFLOW_LOG_VALUE_RETAINED')


def _path(path):return str(Path(path).relative_to(ROOT))


def main():
    inference=ROOT/'results/P1B/inference';validation=ROOT/'results/P1B/point_validation'
    for path in (inference/'COMPLETE.json',validation/'COMPLETE.json'):
        if not path.exists():raise RuntimeError('Handoff requires completed inference AND completed point validation')
    completion=json.loads((inference/'COMPLETE.json').read_text())
    point_completion=json.loads((validation/'COMPLETE.json').read_text())
    if point_completion.get('status')!='FULL_POINT_VALIDATION_COMPLETE':raise RuntimeError('Point-validation completion status incomplete')
    rows=json.loads((inference/'RESULTS.json').read_text())
    if not rows or rows[0].get('observation')!=OBSERVATION:raise RuntimeError('Fixed observation000 must be first; no retrospective example selection')
    row=rows[0];obsroot=(inference/'observation_000').resolve();attempt=(obsroot/row['attempt']).resolve()
    if attempt.parent!=obsroot or not attempt.is_dir():raise RuntimeError('Invalid fixed observation attempt')
    finite_path=attempt/'FINITE_REFERENCE.npz';grid_path=inference/'PUBLIC_MATERIAL_GRID.npz'
    with np.load(attempt/'OBSERVATION.npz',allow_pickle=False) as z:observed=z['y'].copy()
    full9_path=ROOT/'results/P1A/rules/full_9.npz'
    with np.load(finite_path,allow_pickle=False) as z:finite={k:z[k].copy() for k in z.files}
    with np.load(grid_path,allow_pickle=False) as z:material={k:z[k].copy() for k in z.files}
    with np.load(full9_path,allow_pickle=False) as z:paid={k:z[k].copy() for k in z.files}
    shape=tuple(int(n) for n in finite['grid_shape']);m=int(np.prod(shape[:3]))
    if tuple(int(v) for v in material['shape'])!=shape[:3] or m!=len(material['branches']):
        raise RuntimeError('Finite weights and public material rule shape disagree')
    response_source=row['finite_ref'].get('response_source') or row.get('response_source')
    if response_source is None:raise RuntimeError('Finite-response source metadata missing; do not relabel predictions')
    response_kind=response_source['response_kind']
    response_key=str(np.asarray(finite.get('material_response_key','F_full' if response_kind=='FULL_CACHED_NODE' else 'F_approx')).item())
    if response_key not in material:raise RuntimeError('Public material response key absent')
    Fmaterial=material[response_key]
    for saved,key in (('material_branches','branches'),('material_u','u'),('material_contrast','contrast'),('material_x','x')):
        if not np.array_equal(finite[saved],material[key]):raise RuntimeError('Finite factors and public material grid disagree')
    w=joint_factors(finite['posterior_weights'],shape,Fmaterial)
    amplitude_axis=np.asarray(finite['amplitude_axis'],float);phase_axis=np.asarray(finite['phase_axis'],float)
    if amplitude_axis.shape!=(shape[3],) or phase_axis.shape!=(shape[4],):raise RuntimeError('Wrong nuisance axes')
    factors=dependent_prediction(w,Fmaterial,amplitude_axis,phase_axis)
    aa,pp=np.meshgrid(amplitude_axis,phase_axis,indexing='ij')
    continuous=row['continuous_ref'];finite_info=row['finite_ref']
    upstream=ROOT/'results/P1A/reference/GATE.json'
    upstream_gate=json.loads(upstream.read_text()) if upstream.exists() else {}
    qualified=(continuous.get('qualified') is True and finite_info.get('qualified') is True
               and upstream_gate.get('probability_reference_qualified') is True)
    status='QUALIFIED_CONTROLLED_NUMERICAL_EXAMPLE' if qualified else 'UNQUALIFIED'
    probabilities=np.asarray(continuous.get('branch_mass',[np.nan,np.nan]),float)
    if qualified and (probabilities.shape!=(2,) or not np.isfinite(probabilities).all()
        or np.any(probabilities<0) or np.any(probabilities>1) or abs(float(probabilities.sum())-1)>1e-10):
        raise RuntimeError('Qualified continuous branch probabilities are invalid')
    family=Family.load()
    if family.config.get('family_id')!=FAMILY_ID:raise RuntimeError('Unexpected public family')
    point_file=validation/'observation_000/RESULT.json'
    point_results=json.loads(point_file.read_text())
    if point_results.get('observation')!=OBSERVATION:raise RuntimeError('Wrong point-validation observation')
    full9_metadata=ROOT/'results/P1A/rules/full_9.json'
    paid_metadata=json.loads(full9_metadata.read_text()) if full9_metadata.exists() else None
    events=[];representatives=[];branch_predictives=[];conditionals=[];paid_indices=[]
    for branch in (0,1):
        use=material['branches']==branch;normalizer=float(factors['material_masses'][use].sum())
        if normalizer<=0:raise RuntimeError('No positive finite conditional weight in displayed branch')
        mw=factors['material_masses'][use];gw=w[use].sum(axis=0)
        stats={name:weighted_finite_summary(nodes,weights) for name,nodes,weights in (
            ('x',material['x'][use],mw),('u',material['u'][use],mw),('contrast',material['contrast'][use],mw),
            ('amplitude',aa.ravel(),gw),('phase',pp.ravel(),gw))}
        mean_u=stats['u']['mean'];mean_c=stats['contrast']['mean']
        candidates=np.flatnonzero(np.asarray(paid['branches'])==branch)
        if not len(candidates):raise RuntimeError('Paid full9 lacks requested public branch')
        distances=((paid['u'][candidates]-mean_u)/.02)**2+((paid['contrast'][candidates]-mean_c)/.12)**2
        index=int(candidates[int(np.argmin(distances))]);paid_indices.append(index)
        b=int(paid['branches'][index]);u=float(paid['u'][index]);c=float(paid['contrast'][index])
        a=stats['amplitude']['mean'];p=stats['phase']['mean'];gamma=a*np.exp(1j*p)
        raw=np.array(paid['F'][index],copy=True)
        if raw.shape!=(5,4,12) or not np.iscomplexobj(raw) or not np.isfinite(raw).all():raise RuntimeError('Invalid paid full9 raw prediction')
        source=dict(rule_file=_path(full9_path),node_index=index,raw_array_key='F',
            source_type='DIRECT_PAID_FULL_MAXWELL_PUBLIC_PRIOR_NODE',new_RHS=0,original_node_evaluation_RHS=20)
        if paid_metadata is not None:
            node_path=paid_metadata['nodes'][index].get('path')
            source['bank_node_response_path']=str(Path('results/P1A/bank')/node_path) if node_path else None
        representative=dict(branch=b,u=u,contrast=c,amplitude=a,phase=p,
            target_center=family.coordinates(b,u,c)['center'],target_radius=.12,target_conductivity=.02,
            selection='existing paid full9 node nearest finite conditional(u,c)mean by prior-width normalized distance',
            normalized_distance_squared=float(np.min(distances)),gamma_definition='E[a|branch]*exp(i*E[phase|branch])',
            representative_not_posterior_mean_material=True,point_probability=None,full9_prediction_source=source)
        observed_distance=float(np.linalg.norm(gamma*raw-observed)/.01)
        representative.update(d_observed=observed_distance,Phi=.5*observed_distance**2,
            residual_scope='Noisy-observation fit, not the two-world mean-distance<=1 criterion')
        representatives.append(dict(metadata=representative,Fraw=raw,gammaF=gamma*raw))
        pred=np.tensordot(factors['gamma_weighted_material_factors'][use],Fmaterial[use],axes=(0,0))/normalizer
        branch_predictives.append(pred);conditionals.append(stats)
        pb=float(probabilities[branch]) if qualified else None
        if qualified and pb>0:
            conditional_logz=float(continuous['logZ']+np.log(pb)-np.log(.5));logz_status='CONTINUOUS_REFERENCE_QUADRATURE_APPROXIMATION'
        else:
            conditional_logz=None;logz_status='UNQUALIFIED' if not qualified else 'ZERO_ESTIMATED_EVENT_MASS_NEGATIVE_INFINITY_LIMIT'
        center=family.config['target']['branch_centers'][branch]
        events.append(dict(event_id=f'public_wall_branch_{branch}',object_type='PUBLIC_WALL_AND_TARGET_INTERVAL_EVENT',
            definition=dict(wall_parameter_file='inputs/W_FAMILY_WALLS.npz',wall_key=f'q{branch}',
                target_x_interval=[center[0]-.01,center[0]+.01],target_y=center[1],
                target_radius=.12,contrast_interval=[2.94,3.06],conductivity=.02,
                remaining_coordinates='all declared shared-nuisance prior support'),
            prior_mass=.5,posterior_mass=pb,
            posterior_mass_source='same-observation independent full5/full9 continuous-reference integration',
            prior_odds_vs_other=1.,posterior_odds_vs_other=_odds(pb,probabilities[1-branch] if qualified else None,qualified),
            conditional_logZ=conditional_logz,conditional_logZ_status=logz_status,
            conditional_logZ_formula='global_continuous_logZ+log(p_branch)-log(.5)',
            numerical_error_diagnostics=dict(continuous_order5_vs9_event_abs_change=continuous.get('event_order_error'),
                finite_vs_continuous_branch_abs_error=(abs(normalizer-float(probabilities[branch])) if probabilities.shape==(2,) and np.isfinite(probabilities[branch]) else None),
                not_a_certified_error_bound=True),
            finite_conditional_uncertainty=stats,uncertainty_status='FINITE_QUADRATURE' if qualified else 'UNQUALIFIED_FINITE_DIAGNOSTIC_ONLY',
            representative=representative,point_probability=None,status=status))
    x_axis=np.linspace(-2.,2.,400);y_axis=np.linspace(0.,3.,300)
    xx,yy=np.meshgrid(x_axis,y_axis);grid=np.stack((xx,yy),axis=-1)
    decoded=[]
    for rep in representatives:
        r=rep['metadata'];decoded.append(family.material(r['branch'],r['u'],r['contrast'],grid,f=1.))
    epsilon=np.asarray(decoded);raw=np.stack([r['Fraw'] for r in representatives]);gained=np.stack([r['gammaF'] for r in representatives])
    out=ROOT/'results/P1B/handoff/observation_000';out.mkdir(parents=True,exist_ok=True)
    arrays=out/'PREDICTIONS_AND_MATERIALS.npz'
    np.savez_compressed(arrays,grid_x=x_axis,grid_y=y_axis,frequency=1.,
        epsilon_real=epsilon.real,conductivity_sigma=epsilon.imag,epsilon_imag_loss=epsilon.imag,
        Fraw=raw,gammaF=gained,posterior_predictive_approx=factors['posterior_predictive_approx'],
        branch_conditional_predictive_approx=np.asarray(branch_predictives),
        paid_full9_indices=np.asarray(paid_indices),material_posterior_masses=factors['material_masses'],
        gamma_weighted_material_factors=factors['gamma_weighted_material_factors'])
    points={}
    stems={'direct joint point':'direct_joint_point','progressive block point':'progressive_block_point'}
    for name,point in point_results['methods'].items():
        points[name]=dict(parameters={k:point[k] for k in ('branch','u','contrast','amplitude','phase')},
            Phi=point['Phi'],full_validation_status=point['full_validation_status'],fallback_used=point['fallback_used'],
            physical_source=point['physical_source'],prediction_file=_path(validation/'observation_000'/(stems[name]+'_FINAL.npz')),
            posterior_probability=None,meaning='point estimate only; probability mass of a singleton is not assigned')
    samples=[]
    for schedule,records in row['SMC_runs'].items():
        for record in records:
            path=attempt/f'{schedule}_seed{record["seed"]}_TRAJECTORY.npz'
            if not path.exists():raise RuntimeError('Completed inference lacks a required trajectory artifact')
            samples.append(dict(schedule=schedule,seed=record['seed'],path=_path(path),
                final_particle_indexes_key='particle_indexes',final_weights_key='weights',
                all_repeated_states_key='trajectory_indexes',all_stage_weights_key='trajectory_weights',
                acceptance_key='mutation_accepted',resampling_key='resampling_indexes',grid_shape=list(shape),
                status=status,trajectory_scope='changing stage targets; only final particle_indexes/weights describe the final finite target'))
    example=dict(schema='W_PROBABILITY_HANDOFF_V1',observation=OBSERVATION,example_selection='fixed observation000 before data; no retrospective selection',
        status=status,family_id=FAMILY_ID,acquisition=family.config['acquisition'],target_identity=family.config['target'],
        input_complete_files=[_path(inference/'COMPLETE.json'),_path(validation/'COMPLETE.json')],
        reference_qualification=dict(observation0_continuous=continuous.get('qualified'),observation0_finite=finite_info.get('qualified'),
            upstream_P1A=upstream_gate.get('probability_reference_qualified'),
            all_P1B_finite_references=completion.get('all_finite_references_qualified')),
        event_partition=dict(events=events,exhaustive_within_controlled_two_wall_family=True,
            represented_mass=float(probabilities.sum()) if qualified else None,
            unshown_mass_within_partition=0. if qualified else None,
            topK_renormalization=False,selection='both declared branch events are shown'),
        prior_branch1_vs_branch0_odds=1.,posterior_branch1_vs_branch0_odds=_odds(
            probabilities[1] if qualified else None,probabilities[0] if qualified else None,qualified),
        global_logZ=float(continuous['logZ']) if qualified else None,
        points=points,point_plot_probability=None,
        prediction_arrays=dict(path=_path(arrays),Fraw='direct paid full9 nodes only',gammaF='representative gamma times paid Fraw',
            posterior_predictive_approx='sum_m Fmaterial_m * sum_gamma(gamma*w_mgamma), preserving nuisance-material dependence; finite integration approximation',
            material_response_source=response_source,
            branch_conditionals='same joint finite weights, normalized within branch; not reweighted to continuous branch masses',
            factorized_prediction_difference_norm=factors['factorized_predictive_difference_norm'],
            factorized_prediction_difference_norm_units='raw complex-response Euclidean norm; diagnostic only',
            loss_definition='Im(epsilon)=conductivity/f; f1 here, no dissipated-power estimate',
            material_grid=dict(nx=400,ny=300,x_domain=[-2,2],y_domain=[0,3],sampling='fixed endpoint-inclusive nodal grid')),
        compensation_association=dict(historical_target_center_translation=.1,historical_instrument_gain_fixed=1.,
            historical_witness_role='conditional construction history only; does not assign branch probability',
            current_prior_role='new explicit branch masses .5/.5 plus continuous target/nuisance prior',
            current_shared_gain='marginalized once for all frequencies/sources/receivers',
            representative_pair_distance_raw_gain1=float(np.linalg.norm(raw[1]-raw[0])/.01),
            representative_pair_distance_with_branch_representative_gain=float(np.linalg.norm(gained[1]-gained[0])/.01),
            representative_pair_status='Full paid production-grid predictions, not automatically a newly nominated independently grid-validated witness',
            representative_target_x_difference=representatives[1]['metadata']['target_center'][0]-representatives[0]['metadata']['target_center'][0]),
        downstream=dict(status=status,finite_joint_weights=dict(path=_path(finite_path),key='posterior_weights',grid_shape=list(shape),
                ravel_order=['branch','iu','ic','ia','iphi'],material_coordinate_keys=['material_branches','material_u','material_contrast','material_x'],
                nuisance_axis_keys=['amplitude_axis','phase_axis'],prior_factor_keys=['material_prior_weights','gamma_axis_weights'],
                prior_factorization='.5*(wu/2)*(wc/2)*(wa/2)*(wphase/2); same continuous controlled prior'),
            public_material_response=dict(path=_path(grid_path),key=response_key,response_kind=response_kind,
                source_metadata=response_source,interpretation='paid full-Maxwell node vectors' if response_kind=='FULL_CACHED_NODE' else 'explicit legacy response approximation'),
            correlation_factors=dict(path=_path(arrays),material_masses_key='material_posterior_masses',gamma_weighted_material_key='gamma_weighted_material_factors'),
            weighted_samples=samples),
        model_adequacy='controlled synthetic family only; no external-world adequacy claim',
        original_minimum_feature='UNRESOLVED_ORIGINAL_SCALAR',outside_family_mass=None,
        unexplored_original_wall_prior='UNQUANTIFIED',general_unknown_wall_online_claim=False,
        new_PDE=0,new_RHS=0,new_LU=0,truth_access=False)
    _write(out/'EXAMPLE.json',example)
    _write(out/'COMPLETE.json',dict(status=status,observation=OBSERVATION,example=_path(out/'EXAMPLE.json'),
        arrays=_path(arrays),new_PDE=0,new_RHS=0,new_LU=0,truth_access=False))
    return example


if __name__=='__main__':main()
