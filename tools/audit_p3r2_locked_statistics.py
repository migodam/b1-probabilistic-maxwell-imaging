#!/usr/bin/env python3
"""Independent post-scoring audit of frozen risk arrays, never truth or models.

Run only after BOTH fresh scoring COMPLETE markers. This tool does not create
a Book: the caller must run it in the main D analysis-only accounting wrapper.
Only scoring summaries, the point RAW_RISKS archive and frozen manifest are read.
"""
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_key]='1'
import argparse
import json
import re
from pathlib import Path
from statistics import NormalDist
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_OUT=ROOT/'research/p3l_r2_risk_posterior_exploration'
SEEDS=(11,22,33)
SECONDARY_PAIRS=(
    ('local_geometry_cnn','local_anchor_cnn'),
    ('local_geometry_cnn','local_capacity_cnn'),
    ('local_geometry_cnn','local_permuted_geometry_cnn'),
    ('weak_strong_decoder','random_decoder'),
    ('weak_strong_decoder','all_decoder_matched'),
    ('weak_strong_decoder','u2_decoder'),
)

def read_json(path):
    return json.loads(Path(path).read_text())

def paired(candidate,baseline,indices,relative=True):
    candidate=np.asarray(candidate,dtype=float);baseline=np.asarray(baseline,dtype=float)
    if candidate.shape!=(256,) or baseline.shape!=(256,):
        raise ValueError('Comparisons require exactly256 aligned independent parents')
    delta=candidate-baseline;means=delta[indices].mean(1)
    result=dict(n_parent=256,mean_candidate=float(candidate.mean()),mean_baseline=float(baseline.mean()),
                mean_difference=float(delta.mean()),direction='negative difference favors candidate')
    if relative:
        if baseline.mean()<=0 or np.any(baseline[indices].mean(1)<=0):
            raise ValueError('Relative risk reduction undefined for nonpositive baseline')
        result['relative_risk_reduction']=float(1-candidate.mean()/baseline.mean())
        reductions=-means/baseline[indices].mean(1)
    for confidence,label in ((.95,'95'),(.975,'97_5')):
        tail=(1-confidence)/2
        result['difference_CI_'+label]=np.quantile(means,[tail,1-tail]).tolist()
        if relative:result['relative_reduction_CI_'+label]=np.quantile(reductions,[tail,1-tail]).tolist()
    return result

def wilson(successes,n=256):
    z=NormalDist().inv_cdf(.975);p=successes/n;den=1+z*z/n
    center=(p+z*z/(2*n))/den
    half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [float(max(0,center-half)),float(min(1,center+half))]

def compare_stored(calculated,stored):
    keys=('mean_difference','relative_risk_reduction','difference_CI_95',
          'difference_CI_97_5','relative_reduction_CI_95','relative_reduction_CI_97_5')
    return {key:bool(key in stored and np.allclose(calculated[key],stored[key],rtol=1e-10,atol=1e-12)) for key in keys}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=DEFAULT_OUT)
    parser.add_argument('--methods',type=Path)
    args=parser.parse_args();out=args.root.resolve()
    # Gate before reading ANY risk arrays/results or frozen method contents.
    markers={}
    for name in ('evaluation','trajectory_evaluation'):
        path=out/'L'/name/'COMPLETE.json'
        if not path.is_file():raise RuntimeError('Both fresh scores must be COMPLETE before audit: '+str(path))
        marker=read_json(path)
        if marker.get('status')!='COMPLETE' or marker.get('count')!=256:
            raise RuntimeError('Fresh scoring not COMPLETE for all256: '+name)
        markers[name]=marker
    dest=out/'AUDIT/locked_statistics'
    if dest.exists():raise FileExistsError('No overwrite of existing independent audit: '+str(dest))
    methods=read_json(args.methods or out/'METHODS_FROZEN_FOR_LOCKED.json')
    if methods.get('status')!='FROZEN':raise ValueError('Frozen manifest required')
    parents=methods['parents']
    if len(parents)!=256 or len(set(parents))!=256:raise ValueError('Manifest must have256 unique parents')
    for name,marker in markers.items():
        if marker.get('parents')!=parents or marker.get('source_revision')!=methods['source_revision'] or float(marker.get('frozen_unix',-1))!=float(methods['frozen_unix']):
            raise ValueError('Scoring marker cohort/freeze mismatch: '+name)
    point=read_json(out/'L/evaluation/RESULTS.json')
    trajectory=read_json(out/'L/trajectory_evaluation/RESULTS.json')
    for name,result in (('point',point),('trajectory',trajectory)):
        if result.get('status')!='COMPLETE' or result.get('count')!=256 or result.get('parents')!=parents or result.get('no_exclusions') is not True:
            raise ValueError('Scoring result cohort/no-exclusion mismatch: '+name)
        if result.get('source_revision')!=methods['source_revision'] or float(result.get('frozen_unix',-1))!=float(methods['frozen_unix']):
            raise ValueError('Scoring result freeze mismatch: '+name)
    registered=methods['ablation_checkpoints'];registered_names=[entry['name'] for entry in registered]
    if len(registered_names)!=77 or len(set(registered_names))!=77:
        raise ValueError('This prespecified audit requires exactly77 unique frozen ablations')
    if point.get('frozen_ablation_checkpoints')!=registered:
        raise ValueError('Scored ablation registry differs from frozen manifest')
    with np.load(out/'L/evaluation/RAW_RISKS.npz',allow_pickle=False) as archive:
        raw_parents=archive['parent'].astype(str).tolist();names=archive['method'].astype(str).tolist()
        risks={key:archive[key].copy() for key in ('fullgrid_MSE','ROI_MSE','fullgrid_MAE','ROI_MAE')}
        density={name:{key:archive[name+'_'+key].copy() for key in ('nll','unit_prior_nll','coverage90','used_prior_fallback')}
                 for name in ('density','weak_density')}
        fallback_flags=archive['method_used_fallback'].copy() if 'method_used_fallback' in archive else None
    if raw_parents!=parents or len(set(raw_parents))!=256 or len(names)!=len(set(names)):
        raise ValueError('Raw risk parents or method identifiers are not uniquely aligned')
    expected_ablation_names={'ablation__'+name for name in registered_names}
    if {name for name in names if name.startswith('ablation__')}!=expected_ablation_names:
        raise ValueError('Raw risk arrays must score all and only77 registered ablations')
    for key,values in risks.items():
        if values.shape!=(256,len(names)) or not np.isfinite(values).all() or np.any(values<0):
            raise ValueError('Invalid complete per-parent risks: '+key)
    index={name:j for j,name in enumerate(names)}
    for name in ('fused','probit'):
        if name not in index:raise ValueError('Missing frozen primary method '+name)
    indices=np.random.default_rng(2026102292).integers(0,256,size=(2000,256))
    checks={}
    coherence={}
    for metric in ('MSE','MAE'):
        actual=risks['fullgrid_'+metric];expected=risks['ROI_'+metric]*(256/576)
        mismatch=~np.isclose(actual,expected,rtol=1e-9,atol=1e-12)
        coherence[metric]=dict(passed=not bool(mismatch.any()),max_absolute_difference=float(np.max(np.abs(actual-expected))),
            mismatched_cells=int(mismatch.sum()),affected_methods=[name for j,name in enumerate(names) if mismatch[:,j].any()],
            condition='Factor256/576 requires exact common host outside256-cell ROI; failures need method-specific review, not automatic metric repair.')
        checks['fullgrid_ROI_'+metric]=coherence[metric]['passed']
    baseline=risks['fullgrid_MSE'][:,index['probit']]
    primary_point=paired(risks['fullgrid_MSE'][:,index['fused']],baseline,indices)
    primary_point['stored_reproduction']=compare_stored(primary_point,point['primary_comparison'])
    checks['point_primary_reproduced']=all(primary_point['stored_reproduction'].values())
    # Trajectory summaries already expose per-parent SCORED RISKS; no online
    # predictions, trajectory raw archives, checkpoints, contexts or GT are read.
    primary_name='dynamic_curvature__arc_length_trapezoid'
    if trajectory.get('primary_arm')!='dynamic_curvature' or trajectory.get('primary_method')!='arc_length_trapezoid':
        raise ValueError('Prespecified trajectory primary was changed')
    trajectory_rows=trajectory['per_parent_risks']
    if set(trajectory_rows)!=set(parents):raise ValueError('Trajectory risk summary excludes/duplicates parents')
    arc=np.array([trajectory_rows[p][primary_name]['fullgrid_MSE'] for p in parents])
    trajectory_baseline=np.array([trajectory_rows[p]['probit']['fullgrid_MSE'] for p in parents])
    if not np.isfinite(arc).all() or np.any(arc<0):raise ValueError('Invalid trajectory primary risks')
    checks['same_primary_baseline']=bool(np.allclose(baseline,trajectory_baseline,rtol=1e-10,atol=1e-12))
    primary_arc=paired(arc,baseline,indices)
    primary_arc['stored_reproduction']=compare_stored(primary_arc,trajectory['co_primary_vs_strong_probit'])
    checks['trajectory_primary_reproduced']=all(primary_arc['stored_reproduction'].values())
    bootstrap_protocol=methods.get('bootstrap',{})
    checks['frozen_bootstrap_protocol']=bootstrap_protocol.get('draws')==2000 and bootstrap_protocol.get('seed')==2026102292
    checks['two_frozen_primary_comparisons']=len(methods.get('primary_comparisons',[]))==2
    families={}
    for candidate,control in SECONDARY_PAIRS:
        for family in (candidate,control):
            if family in families:continue
            columns=[]
            for seed in SEEDS:
                # n512 and seeds are prespecified. Never choose size, seed or
                # control according to locked risk values.
                name=f'{family}_n512_s{seed}'
                matches=[entry for entry in registered if entry['name']==name and entry['kind']=='residual']
                if len(matches)!=1:raise ValueError('Missing prespecified residual family member '+name)
                columns.append(index['ablation__'+name])
            families[family]=risks['fullgrid_MSE'][:,columns]
    secondary={}
    for candidate,control in SECONDARY_PAIRS:
        c,b=families[candidate],families[control]
        value=paired(c.mean(1),b.mean(1),indices)
        seed_differences=c.mean(0)-b.mean(0)
        seed_reductions=1-c.mean(0)/b.mean(0)
        value.update(seeds=list(SEEDS),individual_seed_mean_difference=seed_differences.tolist(),
            individual_seed_relative_reduction=seed_reductions.tolist(),
            individual_seed_relative_reduction_range=[float(seed_reductions.min()),float(seed_reductions.max())],
            individual_seed_candidate_MSE_range=[float(c.mean(0).min()),float(c.mean(0).max())],
            individual_seed_control_MSE_range=[float(b.mean(0).min()),float(b.mean(0).max())],
            count_improved_training_seeds=int((seed_differences<0).sum()),n_training_seeds=3,
            estimand='Mean of3fixed training-seed RISKS per parent, not risk of an averaged-prediction ensemble.',
            uncertainty='Bootstrap256parents only; not n768. Conditional on these3seeds; training-seed uncertainty is not captured.',
            multiplicity='Prespecified secondary family comparisons; nominal intervals unadjusted across secondary tests; no primary reselection.')
        secondary[candidate+'_vs_'+control]=value
    density_audit={}
    for name,data in density.items():
        coverage=np.asarray(data['coverage90'])
        if coverage.shape!=(256,10) or not np.isin(coverage,[0,1]).all():raise ValueError('Coverage requires256x10 binary outcomes')
        coverage=coverage.astype(float);nll=np.asarray(data['nll']);prior=np.asarray(data['unit_prior_nll'])
        if nll.shape!=(256,) or prior.shape!=(256,) or not np.isfinite(nll).all() or not np.isfinite(prior).all():
            raise ValueError('Density NLL needs256 finite parent scores')
        macro=coverage.mean(1);draw_macro=macro[indices].mean(1)
        components=[dict(component=j,covered=int(coverage[:,j].sum()),n_parent=256,
                         coverage=float(coverage[:,j].mean()),Wilson_CI_95=wilson(int(coverage[:,j].sum()))) for j in range(10)]
        density_audit[name]=dict(components=components,macro_coverage=float(macro.mean()),
            macro_coverage_CI_95=np.quantile(draw_macro,[.025,.975]).tolist(),
            macro_coverage_CI_97_5=np.quantile(draw_macro,[.0125,.9875]).tolist(),
            macro_resampling_unit='256paired parents, preserving all10 dependent coordinates; not2560 independent Bernoulli draws',
            nll_vs_own_unit_prior=paired(nll,prior,indices,relative=False),
            fallback_count=int(np.asarray(data['used_prior_fallback']).sum()),
            interpretation='Finite1024-sample marginal90% intervals; Wilson uncertainty treats parents independent, not coordinates independent. No simultaneous10-coordinate coverage claim.')
    warning=['All-constant MCMC chains can have MCSE0 and remain wrong/unmixed. MCSE0 is not accuracy, mass or posterior-fidelity evidence; no chain data are read by this audit.',
       'Informative8 and risk-weak8 densities target different fixed-coordinate marginals. Their raw NLL values are not directly rankable across frames; compare each with its own unit-prior reference.',
       'These are supervised material-risk and fixed-coordinate density scores. No Bayes-mean identity, full posterior calibration, novelty or GO decision follows.',
       'Two co-primary comparisons use two-sided97.5% bootstrap intervals (Bonferroni for planned family size2). Secondary intervals do not provide simultaneous coverage.',
       'The method family, three training seeds, n512 and all pairings are fixed in this audit source; no held-out model/seed/size reselection is performed.',
       'Uncertainty is conditional on fixed trained checkpoints and1024interval draws; training-seed and finite density-sampling variability are not resampled.']
    result=dict(status='PASS_STATISTICAL_COHERENCE' if all(checks.values()) else 'REVIEW_REQUIRED',
       checks=checks,parents=parents,n_parent=256,n_registered_ablation=77,n_methods=len(names),no_exclusions=True,
       source_revision=methods['source_revision'],frozen_unix=methods['frozen_unix'],
       input_scope=['L/evaluation/COMPLETE.json','L/trajectory_evaluation/COMPLETE.json','METHODS_FROZEN_FOR_LOCKED.json',
          'L/evaluation/RESULTS.json','L/evaluation/RAW_RISKS.npz','L/trajectory_evaluation/RESULTS.json'],
       physics_calls=0,new_RHS=0,training_updates=0,GT_read=False,contexts_read=False,checkpoints_read=False,
       primary_fixed_two=dict(fused_vs_probit=primary_point,curvature_arc_vs_probit=primary_arc),
       fullgrid_ROI_coherence=coherence,secondary_three_seed_family_mean_risks=secondary,
       density=density_audit,reported_point_fallback_counts=point.get('per_method_fallback_counts'),
       baseline_contains_fallbacks=point.get('primary_baseline_contains_fallbacks'),
       raw_method_fallback_counts={name:int(fallback_flags[:,j].sum()) for j,name in enumerate(names)} if fallback_flags is not None else None,
       bootstrap=dict(draws=2000,seed=2026102292,unit='paired independent parent',n=256),warnings=warning,
       final_scientific_decision='Reserved for parent; this tool reports statistical checks only.')
    dest.mkdir(parents=True,exist_ok=False)
    (dest/'RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    lines=['# Independent locked statistical audit','',result['status'],
       '', 'Both fresh scoring COMPLETE markers were required before reading scored risks. All256parents and77registered ablations are retained.',
       '', '## Two fixed primary comparisons','', '| Comparison | Relative risk reduction | Two-sided97.5% CI |','|---|---:|---|']
    for label,value in result['primary_fixed_two'].items():
        lines.append(f"| {label} | {value['relative_risk_reduction']:.7g} | {value['relative_reduction_CI_97_5']} |")
    lines+=['','## Secondary fixed-three-seed mean risks','','Intervals bootstrap256parents, conditional on3fixed seeds. These are mean risks, not prediction ensembles or n768. Secondary multiplicity is unadjusted.','',
       '| Comparison | Relative reduction | 95% CI | Improved seeds |','|---|---:|---|---:|']
    for label,value in secondary.items():
        lines.append(f"| {label} | {value['relative_risk_reduction']:.7g} | {value['relative_reduction_CI_95']} | {value['count_improved_training_seeds']}/3 |")
    lines+=['','## Metric coherence','']+[f'- {key}: {value}' for key,value in checks.items()]
    lines+=['','## Density coverage','', '| Frame | Macro90% coverage | Parent-bootstrap95% CI | Fallbacks |','|---|---:|---|---:|']
    for name,value in density_audit.items():
        lines.append(f"| {name} | {value['macro_coverage']:.7g} | {value['macro_coverage_CI_95']} | {value['fallback_count']} |")
    lines+=['','Per-component Wilson95% intervals and individual training-seed ranges/counts are in RESULTS.json.','']+['- '+v for v in warning]
    (dest/'REVIEW.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(status=result['status'],output=str(dest),checks=checks,parents=256,registered_ablations=77)))

if __name__=='__main__':main()
