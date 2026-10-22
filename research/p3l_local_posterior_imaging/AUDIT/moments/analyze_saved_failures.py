"""Read existing inference outputs only; no target/physics/ledger imports."""
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text())


def triples(values):
    return dict(min=float(np.min(values)),median=float(np.median(values)),max=float(np.max(values)))


groups=[]
for stage in ['A','C']:
    for path in sorted((ROOT/stage/'references').glob('*/*/REFERENCE_STATS.json')):
        d=read(path)
        groups.append(dict(path=str(path.relative_to(ROOT)),status=d['status'],
            reference_qualified=d['reference_qualified'],gates=d['gates'],
            cross_mechanism_status=d['cross_mechanism_status'],
            full_grid_mean_rms_mcse=d['full_grid_mean_rms_mcse'],summaries=d['summaries']))

methods=[]
initializations={}
for stage,method in [('A','mala_retracted'),('C','mala_decoder')]:
    chains=[]
    for chain in range(4):
        folder=ROOT/stage/'jobs'/f'reference_DEV_030_{method}_c{chain}'
        closed=read(folder/'CLOSED.json');summary=read(folder/'TRANSITION_SUMMARY.json')
        checkpoint=read(folder/'CHECKPOINT.json');init=read(folder/'INITIALIZATION.json')
        initializations[(method,chain)]=init
        with np.load(folder/'CHAIN.npz',allow_pickle=False) as data:
            w=int(data['warmup']);allq=data['q'];allx=data['x']
            q=allq[w:];x=allx[w:];phi=data['Phi'][w:]
            accepted=data['accepted'][w:];lp=data['log_density'][w:]
            quarter=len(q)//4
            sameq=np.all(q==allq[w-1:-1],axis=1)
            samex=np.all(x==allx[w-1:-1],axis=1)
            counters=closed['counters']
            components={k:int(counters.get(k,0)) for k in
                        ['full_forward_RHS','receiver_transfer_adjoint_RHS','likelihood_gradient_adjoint_RHS']}
            assert sum(components.values())==counters['linear_solve_rhs']
            assert not np.any((~accepted)&(~samex))
            chains.append(dict(chain=chain,job=str(folder.relative_to(ROOT)),
                closed_status=closed['status'],warmup=w,production_draws=len(q),
                completed_job_solver_RHS=int(counters['linear_solve_rhs']),RHS_components=components,
                fresh_states=int(counters['full_forward_calls']),LU=int(counters['LU_factorizations']),
                job_controlled_wall_seconds=closed['controlled_wall_seconds'],
                production_acceptance=float(accepted.mean()),final_step_size=checkpoint['step_size'],
                Phi_first_quarter_mean=float(phi[:quarter].mean()),
                Phi_last_quarter_mean=float(phi[-quarter:].mean()),
                Phi_first_to_last_quarter_drift=float(phi[-quarter:].mean()-phi[:quarter].mean()),
                u_min=q[:,64:].min(0).tolist(),u_max=q[:,64:].max(0).tolist(),
                u_first_quarter_mean=q[:quarter,64:].mean(0).tolist(),
                u_last_quarter_mean=q[-quarter:,64:].mean(0).tolist(),
                u_first_to_last_quarter_drift=(q[-quarter:,64:].mean(0)-q[:quarter,64:].mean(0)).tolist(),
                latent_mean_first_to_last_quarter_L2=float(np.linalg.norm(q[-quarter:,:64].mean(0)-q[:quarter,:64].mean(0))),
                rejected_production_transitions=int((~accepted).sum()),
                rejected_nonrepeat_x=int(np.sum((~accepted)&(~samex))),
                rejected_nonrepeat_q=int(np.sum((~accepted)&(~sameq))),
                transformed_log_density_identity_max_abs_error=float(np.max(abs(lp+phi+.5*np.sum(x*x,axis=1)))),
                production_adaptation=summary['production_adaptation']))
    methods.append(dict(stage=stage,method=method,chains=chains,
        completed_job_solver_RHS=sum(c['completed_job_solver_RHS'] for c in chains),
        summed_job_controlled_wall_seconds=sum(c['job_controlled_wall_seconds'] for c in chains)))

initial_match=[]
for c in range(4):
    a=initializations[('mala_retracted',c)];b=initializations[('mala_decoder',c)]
    initial_match.append(dict(chain=c,weak_coordinates_identical=a['weak_coordinates']==b['weak_coordinates'],
        properties_identical=a['properties']==b['properties'],
        conditional_initial_Phi_difference=b['conditional_initial_Phi']-a['conditional_initial_Phi']))

with np.load(ROOT/'A/references/DEV_030/mala_retracted/REFERENCE_ARRAYS.npz',allow_pickle=False) as a:
    ma=a['material_mean'].copy()
with np.load(ROOT/'C/references/DEV_030/mala_decoder/REFERENCE_ARRAYS.npz',allow_pickle=False) as a:
    mc=a['material_mean'].copy()

pilot=read(ROOT/'B/CONDITIONAL_PILOT_RESULTS.json');rows=pilot['rows']
dimensions=[]
for d in [2,4,8,16,24]:
    rr=[r for r in rows if r['weak_dimension']==d]
    dimensions.append(dict(weak_dimension=d,strong_dimension=64-d,batches=len(rr),
        converged_modes=sum(r['mode_converged'] for r in rr),
        importance_ESS=triples([r['importance_ESS'] for r in rr]),
        maximum_weight=triples([r['max_weight'] for r in rr]),
        GN_mean_minus_mode_RMS=triples([r['GN_material_mean_minus_mode_RMS'] for r in rr]),
        SNIS_mean_minus_GN_RMS=triples([r['importance_mean_minus_GN_RMS'] for r in rr]),
        finite_logZ_IS_minus_GN=triples([r['importance_logZ']-r['GN_logZ'] for r in rr])))

volume_pairs=[]
for parent in pilot['protocol']['parents']:
    for d in [2,4,8,16,24]:
        z=next(r for r in rows if (r['parent'],r['weak_dimension'],r['context'])==(parent,d,'zero'))
        s=next(r for r in rows if (r['parent'],r['weak_dimension'],r['context'])==(parent,d,'shift'))
        profile_delta=s['profile_log_density']-z['profile_log_density']
        volume_delta=-.5*(s['logdet_precision']-z['logdet_precision'])
        gndelta=s['GN_logZ']-z['GN_logZ']
        assert abs(gndelta-profile_delta-volume_delta)<1e-9
        volume_pairs.append(dict(parent=parent,weak_dimension=d,profile_shift_minus_zero=profile_delta,
                                 GN_log_volume_shift_minus_zero=volume_delta,
                                 GN_logZ_shift_minus_zero=gndelta))

prefix={}
for k in [8,32]:
    prefix[str(k)]=dict(mean_difference_to_K128_RMS=triples([
        next(v for v in r['importance_prefix_convergence'] if v['K']==k)['mean_difference_to_K128_RMS'] for r in rows]),
        absolute_logZ_difference_to_K128=triples([abs(next(v for v in r['importance_prefix_convergence'] if v['K']==k)['logZ_estimate']-r['importance_logZ']) for r in rows]))

result=dict(status='ACTUAL_PAID_OUTPUT_FAILURE_AUDIT_NOT_A_POSTERIOR_REFERENCE',
    new_physics_calls=0,new_RHS=0,GT_accessed=False,ledger_modified=False,
    source_cost_note='Completed job CLOSED counters, with exact RHS-component reconciliation; no stage-cumulative totals summed. Authoritative ledger not independently reread.',
    saved_reference_groups=groups,common_parent='DEV_030',methods=methods,
    matched_initializations=initial_match,
    decoder_over_static_wall_ratio=methods[1]['summed_job_controlled_wall_seconds']/methods[0]['summed_job_controlled_wall_seconds'],
    decoder_minus_static_RHS=methods[1]['completed_job_solver_RHS']-methods[0]['completed_job_solver_RHS'],
    descriptive_unqualified_material_mean_difference_RMS=float(np.linalg.norm(mc-ma)/24),
    conditional=dict(batches=len(rows),K=128,all_modes_converged_count=sum(r['mode_converged'] for r in rows),
        strong_gradient_inf=triples([r['strong_gradient_inf'] for r in rows]),
        mode_gradient_by_context={ctx:triples([r['strong_gradient_inf'] for r in rows if r['context']==ctx]) for ctx in ['zero','shift']},
        ESS_below_2_count=sum(r['importance_ESS']<2 for r in rows),
        max_weight_above_point99_count=sum(r['max_weight']>.99 for r in rows),
        dimension_summaries=dimensions,prefix_comparisons=prefix,
        matched_volume_pairs=volume_pairs,
        observed_GN_log_volume_shift_minus_zero=triples([r['GN_log_volume_shift_minus_zero'] for r in volume_pairs]),
        maximum_positive_finite_logZ_gap_row=max(rows,key=lambda r:r['importance_logZ']-r['GN_logZ'])['parent']+' d'+str(max(rows,key=lambda r:r['importance_logZ']-r['GN_logZ'])['weak_dimension'])+' '+max(rows,key=lambda r:r['importance_logZ']-r['GN_logZ'])['context']),
    claim_limits=['Both common-parent MALA methods are unqualified; higher acceptance is not better mixing.',
                  'No calibrated posterior, Bayes action, ground-truth risk improvement or net acceleration established.',
                  'Gaussian decoder moments are exact only for the supplied Gaussian reference.',
                  'Finite SNIS normalizers/moments and GN profile volumes are not independent exact references.'])
(OUT/'ACTUAL_REFERENCE_AND_CONDITIONAL_FAILURES.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:result[k] for k in ['decoder_over_static_wall_ratio','decoder_minus_static_RHS','descriptive_unqualified_material_mean_difference_RMS']},indent=2))
print(json.dumps(result['conditional']['mode_gradient_by_context'],indent=2))
print(json.dumps(prefix,indent=2))
print(json.dumps(result['conditional']['observed_GN_log_volume_shift_minus_zero'],indent=2))
print('ESS<2',result['conditional']['ESS_below_2_count'],'maxweight>.99',result['conditional']['max_weight_above_point99_count'])
