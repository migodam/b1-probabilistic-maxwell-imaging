"""B1 E0: audit and repair summaries using zero new physical evaluations."""
from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
import argparse
import json
import time
import traceback
import numpy as np
from scipy.special import logsumexp
from .continuous import FrozenDensity, compare_summaries

CASES=['distinguishable','competitive','information_limited','gain_mismatch',
       'gain_marginalized','representation_mismatch']

def json_default(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)

def dump(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,default=json_default,allow_nan=False)+'\n')

def ledger_snapshot(root):
    rows=[json.loads(l) for l in (root/'results/PHYSICS_LEDGER.jsonl').read_text().splitlines() if l.strip()]
    contract=json.loads((root/'results/BUDGET_CONTRACT.json').read_text())
    return dict(reservations=sum(r.get('event')=='BEGIN' for r in rows),
                completed=sum(r.get('event')=='END' for r in rows),records=len(rows),
                high_water=contract['high_water_reservations'],cap=contract['cap'])

def read_observation(root,name):
    with np.load(root/f'results/cases/{name}/posterior.npz') as q:
        if 'gain_sd' not in q:raise ValueError('canonical observation lacks gain metadata')
        return {k:np.array(q[k]) for k in ['data','noise_real_sd','topology_prior','gain_sd']}

def topk(result,observation,config,error):
    p=np.array(result['topology_probabilities']);pi=np.array(observation['topology_prior'])
    indices=np.argsort(-p,kind='stable');entries=[]
    for h in indices:
        others=1-pi[h]
        entries.append(dict(H=int(h),event_definition=dict(template=config['topologies'][int(h)],
                       whole_continuous_material_region=config['material_prior'],center='fixed origin',
                       outer_radius=config['outer_radius']),representative_material=result['conditional_material'][h]['mean'],
                       representative_type='conditional mean epsilon, legal within fixed template',
                       representative_full_forward_verified=False,posterior_mass=float(p[h]),
                       mass_error_estimate=float(error),error_is_certified=False,prior_mass=float(pi[h]),
                       prior_odds_vs_rest=float(pi[h]/others),posterior_odds_vs_rest=None if p[h]>=1 else float(p[h]/(1-p[h])),
                       Bayes_factor_vs_solid=float(np.exp(result['log_bayes_factors_vs_solid'][h])),
                       conditional_uncertainty=result['conditional_material'][h]))
    return dict(items=entries[:2],remaining_entries=entries[2:],
                remaining_modeled_mass=float(p[indices[2:]].sum()),topk_renormalized=False,
                outside_adopted_prior_support_mass=0.,unmodeled_geometry_mass=None,
                unmodeled_geometry_status='Broader material model not defined; unquantified, not data-excluded',
                unexplored_modeled_mass=None,unexplored_modeled_status='No uniform mode-coverage certificate',
                exact_image_probability=0.)

def run(root,output=None):
    root=Path(root).resolve();stage=Path(output).resolve() if output else root/'research/b1_next/e0'
    out=stage/'results';out.mkdir(parents=True,exist_ok=True);(out/'cases').mkdir(exist_ok=True)
    if (out/'SUMMARY.json').exists():raise RuntimeError('E0 complete; preserve run, choose a new output for a separate iteration')
    config=json.loads((root/'configs/protocol.json').read_text());before=ledger_snapshot(root)
    authorization=json.loads((stage/'AUTHORIZATION.json').read_text())
    if authorization['new_stage_extra_authorization']!=0:raise ValueError('This E0 entry accepts only zero physical authorization')
    stats=dict(polynomial_response_evaluations=0,by_kind={});started=time.perf_counter();cpu=time.process_time()
    (out/'FAILURES.jsonl').touch(exist_ok=True)
    @contextmanager
    def phase(name):
        t=time.perf_counter();c=time.process_time();n=stats['polynomial_response_evaluations'];status='PASS'
        try:yield
        except Exception as e:
            status='FAIL';row=dict(phase=name,error=repr(e),traceback=traceback.format_exc(),new_physics=0)
            with (out/'FAILURES.jsonl').open('a') as f:f.write(json.dumps(row,ensure_ascii=False)+'\n')
            raise
        finally:
            row=dict(phase=name,status=status,wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c,
                     polynomial_response_evaluations=stats['polynomial_response_evaluations']-n,new_physics_evaluations=0)
            with (out/'JOB_LEDGER.jsonl').open('a') as f:f.write(json.dumps(row,ensure_ascii=False)+'\n')
    dump(out/'PROTOCOL.json',dict(stage='B1-E0-continuous-summary-v1',physical_cap=0,
          numerical_target='Frozen 7x7 tensor response interpolation, original continuous uniform material prior and saved likelihood metadata',
          alternative_response='Frozen 5x5 tensor response interpolation, same prior/data/likelihood',
          polynomial_coordinate_degrees=dict(poly7=6,poly5=4),
          tensor_orders=[25,41,65,97],primary_cdf_grid_points=129,adaptive_cdf_grid_points=33,
          independent_reference='partitioned adaptive 5/9 quadrature, both dimensions subdivided',
          tolerances=dict(cdf_tested_grid=.01,mean_and_quantile_prior_width_fraction=.001,
                          sd_relative=.01,sd_prior_width_floor=1e-6,coherence_prior_width_squared=1e-12),
          physical_solver_interface=False,prior=config['material_prior'],topology_events=config['topologies'],
          new_physics_authorized=0,strict_full_Mie_probability_certificate=False))
    with phase('specified_probability_controls'):
        from .e0_checks import run_checks
        controls=run_checks(root);dump(out/'CONTROLS.json',controls)
    polynomial=np.load(root/'results/refinement/FORWARD_POLYNOMIAL.npz');bounds=polynomial['bounds'];width=np.diff(bounds,axis=1)[:,0]
    observations={n:read_observation(root,n) for n in CASES};cases={}
    primary=json.loads((root/'results/PRIMARY_SUMMARY.json').read_text())
    refined=json.loads((root/'results/refinement/SUMMARY.json').read_text())
    for name in CASES:
        checkpoint=out/'cases'/f'{name}.json'
        if checkpoint.exists():
            cases[name]=json.loads(checkpoint.read_text());print('Reused completed E0 checkpoint '+name,flush=True);continue
        obs=observations[name];reports={};arrays={};comparisons={}
        for label in ['poly7','poly5']:
            with phase(name+'_'+label+'_setup'):
                density=FrozenDensity(polynomial['coeff7' if label=='poly7' else 'coeff5'],bounds,obs,stats,label)
            for n in [25,41,65,97]:
                with phase(name+'_'+label+'_order'+str(n)):
                    result,curves=density.summary(n)
                    reports[f'{label}_{n}']=result
                    for k,v in curves.items():arrays[f'{label}_{n}_{k}']=v
            for a,b in [(25,41),(41,65),(65,97)]:
                ac={k:arrays[f'{label}_{a}_{k}'] for k in ['normalized_grid','conditional_cdf','marginal_cdf']}
                bc={k:arrays[f'{label}_{b}_{k}'] for k in ac}
                comparisons[f'{label}_{a}_vs_{b}']=compare_summaries(reports[f'{label}_{a}'],reports[f'{label}_{b}'],ac,bc,width)
            print(f'E0 {name}: {label} tensor integrations complete',flush=True)
        with phase(name+'_independent_adaptive'):
            independent=FrozenDensity(polynomial['coeff7'],bounds,obs,stats,'poly7')
            ar,aa=independent.summary(adaptive=True);reports['adaptive_poly7']=ar
            for k,v in aa.items():arrays['adaptive_poly7_'+k]=v
        fine=reports['poly7_97'];fine_curves={k:arrays[f'poly7_97_{k}'] for k in aa}
        comparisons['poly7_97_vs_adaptive']=compare_summaries(fine,ar,fine_curves,aa,width)
        coarse_curves={k:arrays[f'poly5_97_{k}'] for k in aa}
        comparisons['response_poly5_vs_poly7']=compare_summaries(reports['poly5_97'],fine,coarse_curves,fine_curves,width)
        passed=bool(comparisons['poly7_65_vs_97']['passed'] and comparisons['poly7_97_vs_adaptive']['passed'] and ar['adaptive_converged']
                    and all(r['coherence']['passed'] and r['cdf_monotone'] and r['cdf_range_valid'] for r in fine['conditional_material'])
                    and fine['marginal_material']['coherence']['passed'])
        error=max(comparisons['poly7_65_vs_97']['max_event_mass_difference'],comparisons['poly7_97_vs_adaptive']['max_event_mass_difference'])
        row=dict(case=name,result=fine,integration_results=reports,comparisons=comparisons,
                 status='EMPIRICALLY_CONVERGED_CONTINUOUS_SUMMARY' if passed else 'CONTINUOUS_REFERENCE_UNRESOLVED',
                 certified_bound=False,measurement_target=dict(real_noise_sd=float(obs['noise_real_sd']),
                 complex_noise_sd=float(np.sqrt(2)*obs['noise_real_sd']),gain_sd=float(obs['gain_sd']),topology_prior=obs['topology_prior'].tolist(),packing='all Re then all Im'),
                 legacy_original=primary['cases'][name]['result'],legacy_refined=refined['cases'][name]['result'],
                 original_full_grid_reference=refined['cases'][name]['original_full_grid_reference'],
                 top_k=topk(fine,obs,config,error),curves_file=f'results/cases/{name}_curves.npz')
        np.savez_compressed(out/'cases'/f'{name}_curves.npz',**arrays)
        dump(checkpoint,row);cases[name]=row
        print(f'E0 {name}: {row["status"]}',flush=True)
    sensitivity={}
    for name,row in cases.items():
        entries=[];z=np.array(row['result']['log_evidence'])
        for label,pi in [('original',observations[name]['topology_prior'].tolist()),('solid_favored',[.8,.15,.05]),('cavity_favored',[.05,.15,.8])]:
            p=np.exp(np.log(pi)+z-logsumexp(np.log(pi)+z))
            entries.append(dict(protocol=f'E0-prior-sensitivity-{label}',topology_prior=pi,posterior_probabilities=p.tolist(),
                                likelihood_and_material_prior_unchanged=True,recomputed_response_evaluations=0))
        sensitivity[name]=entries
    dump(out/'PRIOR_AND_NUISANCE_SENSITIVITY.json',dict(prior_sensitivity=sensitivity,
          gain_comparison=dict(observations_identical=bool(np.array_equal(observations['gain_mismatch']['data'],observations['gain_marginalized']['data'])),
             fixed_gain_probabilities=cases['gain_mismatch']['result']['topology_probabilities'],
             marginalized_gain_probabilities=cases['gain_marginalized']['result']['topology_probabilities'],
             fixed_gain_sd=0.,marginalized_gain_sd=.15,candidate_dependent_logdet_retained=True,
             likelihood_protocols_distinct=True),new_physics=0))
    after=ledger_snapshot(root)
    if before!=after:raise RuntimeError('historical physics ledger changed; do not claim zero-physics isolation')
    resources=dict(new_physics_evaluations=0,new_physics_cap=0,historical_physics=after,
                   polynomial_counts=stats,wall_seconds_inference_and_controls=time.perf_counter()-started,
                   cpu_seconds_inference_and_controls=time.process_time()-cpu,
                   GPU_training_updates=0,locked_truth_reads=0,sha256_checks=0,
                   scope='actual E0 process body; Python import/startup, agent services and report rendering measured separately; no full-forward cold/warm claim')
    dump(out/'RESOURCE_SUMMARY.json',resources)
    summary=dict(stage='B1 E0',cases=cases,controls_pass=controls.get('all_asserted_checks_pass'),resources=resources,
                 authorization=authorization,decision='CONTINUE_PROBABILITY_BENCHMARK_REPAIR',
                 all_continuous_cases_converged=all(r['status']=='EMPIRICALLY_CONVERGED_CONTINUOUS_SUMMARY' for r in cases.values()),
                 E1_status='NOT_RUN_NO_NEW_PHYSICS_AUTHORIZATION',P1_4000_authorized=False,
                 mechanism_value='nuisance-aware nonlinear secant new explanatory value not tested in E0',
                 historical_missing_oracle_arrays=4,future_oracle_code_saves_arrays=True,
                 uniform_forward_CDF_error_certificate=False,public_release_authorized=False)
    dump(out/'SUMMARY.json',summary)
    print(json.dumps(dict(decision=summary['decision'],all_continuous_cases_converged=summary['all_continuous_cases_converged'],
                         controls_pass=summary['controls_pass'],resources=resources),ensure_ascii=False,default=json_default),flush=True)
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',default='.');p.add_argument('--output');a=p.parse_args();run(a.root,a.output)
