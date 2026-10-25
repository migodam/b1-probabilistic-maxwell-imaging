#!/usr/bin/env python3
"""Offline saved-terminal pilot summary; no PDE, no new scientific GO gate.

Run through run_job to charge preparation CPU. Statistics use independent parent
means, two repeats first, and descriptive paired bootstrap(seed2026102510,n2000).
"""
from pathlib import Path
from collections import defaultdict,Counter
import argparse,csv,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
SEED=2026102510


def read(path):return json.loads(path.read_text()) if path.exists() else None


def num(x):
    try:
        x=float(x);return x if np.isfinite(x) else None
    except (TypeError,ValueError):return None


def scalar(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)


def parent_means(rows,category,metrics,expected_per_repeat=1):
    grouped=defaultdict(list)
    for row in rows:grouped[(row['parent'],row[category])].append(row)
    output=[]
    for (parent,method),rr in sorted(grouped.items()):
        row=dict(parent=parent,**{category:method},saved_rows=len(rr),expected_rows=2*expected_per_repeat,
                 status_denominators=dict(Counter(str(r.get('status','SAVED')) for r in rr)))
        for metric in metrics:
            values=[];complete=0;valid_rows=0
            for repeat in (0,1):
                block=[r for r in rr if r['repeat']==repeat]
                valid=[num(r.get(metric)) for r in block];valid=[x for x in valid if x is not None]
                valid_rows+=len(valid)
                if valid:values.append(float(np.mean(valid)))
                complete+=int(len(valid)==expected_per_repeat and len(block)==expected_per_repeat)
            row[metric]=float(np.mean(values)) if values else None
            row[metric+'_valid_rows']=valid_rows;row[metric+'_complete_repeats']=complete
        output.append(row)
    return output


def bootstrap(rows,category,metric,baseline):
    values=defaultdict(dict)
    for r in rows:
        if r.get(metric+'_complete_repeats')==2 and num(r.get(metric)) is not None:
            values[r[category]][r['parent']]=r[metric]
    output=[]
    for method in sorted(set(r[category] for r in rows)):
        ids=sorted(set(values[method])&set(values[baseline]))
        d=np.array([values[method][i]-values[baseline][i] for i in ids])
        row=dict(metric=metric,method=method,baseline=baseline,paired_parents=len(ids),seed=SEED,draws=2000,
                 estimate=float(d.mean()) if len(d) else None,low=None,high=None,descriptive_only=True,
                 status='NOT_RUN_INCOMPLETE_12_PARENT_PAIRS')
        if len(ids)==12:
            rng=np.random.default_rng(SEED);b=d[rng.integers(0,12,(2000,12))].mean(axis=1)
            row.update(low=float(np.quantile(b,.025)),high=float(np.quantile(b,.975)),status='DESCRIPTIVE_PAIRED_12_PARENT')
        output.append(row)
    return output


def method_summary(rows,category,metrics):
    output=[]
    for method in sorted(set(r[category] for r in rows)):
        rr=[r for r in rows if r[category]==method]
        row=dict(**{category:method},parent_denominator=12,saved_parent_rows=len(rr))
        for metric in metrics:
            valid=[r[metric] for r in rr if r.get(metric+'_complete_repeats')==2 and num(r.get(metric)) is not None]
            row[metric+'_mean']=float(np.mean(valid)) if valid else None;row[metric+'_complete_parents']=len(valid)
        output.append(row)
    return output


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--output-json',type=Path)
    parser.add_argument('--tables',type=Path)
    args=parser.parse_args();root=args.root.resolve()
    target=args.output_json or root/'results/SUMMARY_FOR_PARENT.json';csvtarget=target.with_suffix('.csv')
    tables=args.tables or root/'results/pilot_summary_tables'
    if target.exists() or csvtarget.exists() or tables.exists():raise FileExistsError('Use new summary paths; existing analysis outputs are preserved')
    target.parent.mkdir(parents=True,exist_ok=True);tables.mkdir(parents=True,exist_ok=False)
    report=dict(scope='OFFLINE_TERMINAL_SAVED_ARRAY_PILOT_ONLY',scientific_decision='PARENT_REVIEW_REQUIRED',
        no_new_gate=True,no_new_PDE=True,new_rhs=0,new_factorizations=0,
        finite_prior_only=True,continuous_prior_probability_certificate=False,
        statistics=dict(unit='12 independent parents; two repeats averaged first; MC seeds inside repeat',bootstrap_draws=2000,seed=SEED,descriptive_only=True),
        missing_inputs=[],E3=dict(status='NOT_RUN'),E4=dict(status='NOT_RUN'),E5=dict(status='NOT_RUN'),
        A0_ORDINARY_ADJOINT_POLISH=dict(status='NOT_RUN'))
    compact=[]
    def table(name,rows):
        if not rows:return
        keys=list(dict.fromkeys(k for r in rows for k in r))
        with (tables/(name+'.csv')).open('w',newline='') as f:
            writer=csv.DictWriter(f,keys);writer.writeheader()
            for r in rows:writer.writerow({k:json.dumps(v,default=scalar) if isinstance(v,(dict,list,np.ndarray)) else v for k,v in r.items()})
    def notrun(section,message):
        report[section]=dict(status='NOT_RUN',reason=message)
        report['missing_inputs'].append(section+': '+message)
        compact.append(dict(section=section,status='NOT_RUN',reason=message))
    bank=root/'results/e2e3/finite_bank';inf=root/'results/e2e3/inference'
    finalcost=read(bank/'online/JOB_COST_FINAL.json')
    setup=dict(rhs=11800,factorizations=2950,measured_CPU=None,summed_job_wall=None,status='NOT_RUN_NO_COMPLETE_FINAL_JOB_COST',
               source='results/e2e3/finite_bank/online/JOB_COST_FINAL.json',warning='Never use in-flight COST_RECEIPT/COMPLETE snapshots as total bank CPU')
    if finalcost and finalcost.get('complete') is True:
        setup.update(measured_CPU=num(finalcost.get('measured_process_cpu_seconds')),summed_job_wall=num(finalcost.get('summed_job_wall_seconds')),status='COMPLETE_WRAPPER_COST',
                     failed_jobs_included=finalcost.get('failed_jobs_included',[]))
    report['common_physical_bank_setup']=setup
    # Public coordinate inputs may be read first. Offline truth is not opened
    # until the complete saved inference mask/reference/record panel is checked.
    postpath=inf/'POSTERIORS.npz';gridpath=bank/'online/FINITE_GRID.npz';log=read(inf/'METHOD_RESULTS.json')
    complete=False
    if postpath.exists() and gridpath.exists() and log:
        with np.load(postpath,allow_pickle=False) as z:
            refs=z['weights'].copy();weights=z['method_weights'].copy();methods=z['methods'].astype(str).tolist()
            complete=bool(z['complete']) if 'complete' in z.files else False
            complete=complete and 'completed_mask' in z.files and z['completed_mask'].shape==(12,2,3) and bool(z['completed_mask'].all())
        complete=complete and refs.shape[:3]==(12,2,3) and weights.shape[:4]==(12,2,3,len(methods)) and len(log)==12*2*3*len(methods)
        complete=complete and np.isfinite(refs).all() and np.allclose(refs.sum(axis=-1),1.,rtol=0,atol=1e-10)
    if complete:
        with np.load(gridpath,allow_pickle=False) as z:grid={k:z[k].copy() for k in ['pose','support','epsinf','ell','support_id','wall_group_id']}
        pose=grid['pose'];support=grid['support'];wall=grid['wall_group_id'];sid=grid['support_id'];N=len(pose)
        if N!=2950:raise ValueError('Finite atom count changed')
        coords=np.column_stack((support*(grid['epsinf'][:,None]-1),support*grid['ell'][:,None],pose))
        ordinary=methods.index('INCREMENTAL_JOINT_SMOOTHER');hard=methods.index('HARD_FROZEN_WALL')
        ordinary_valid=np.isfinite(weights[:,:,:,ordinary]).all() and np.allclose(weights[:,:,:,ordinary].sum(axis=-1),1,rtol=0,atol=1e-10)
        ordinary_max=float(np.max(np.abs(weights[:,:,:,ordinary]-refs))) if ordinary_valid else None
        raw=[]
        for r in log:
            metrics=r['posterior_metrics'];tau=.05*metrics['reference_bayes_risk']+.005**2
            raw.append(dict(parent=r['parent'],repeat=r['repeat'],time=r['time'],method=r['method'],status='SAVED',
                target_CDF_sup=float(max(metrics['coordinate_cdf_sup'][-2:])),event_error=metrics['event_px_ge_0_error'],mode_error=metrics['mode_mass_max_error_ref_ge_005'],
                action_regret=metrics['action_regret'],task_tau=tau,
                joint_total_variation=metrics.get('joint_total_variation'),
                wall_group_total_variation=metrics.get('wall_group_total_variation'),
                joint_atom_max_probability_error=metrics.get('joint_atom_max_probability_error'),
                CDF_fraction=float(max(metrics['coordinate_cdf_sup'][-2:])<=.02),event_fraction=float(metrics['event_px_ge_0_error']<=.02),
                mode_fraction=float(metrics['mode_mass_max_error_ref_ge_005']<=.02),task_fraction=float(metrics['action_regret']<=tau),
                all_probability_fraction=float(max(metrics['coordinate_cdf_sup'][-2:])<=.02 and metrics['event_px_ge_0_error']<=.02 and metrics['mode_mass_max_error_ref_ge_005']<=.02),
                standalone_cpu_seconds=num(r.get('standalone_cpu_seconds')),action=r['action'],
                decision_reference_access=(r.get('criterion') or {}).get('current_reference_access_for_decision'),
                reference_threshold_requires_full_update=(r.get('criterion') or {}).get('reference_threshold_requires_full_update',False)))
        final=[r for r in raw if r['time']==3]
        metrics=['target_CDF_sup','event_error','mode_error','action_regret','joint_total_variation','wall_group_total_variation','joint_atom_max_probability_error','CDF_fraction','event_fraction','mode_fraction','task_fraction','all_probability_fraction']
        pg=parent_means(final,'method',metrics);ms=method_summary(pg,'method',metrics)
        table('E3_ALL_PREFIX_METHOD_RECORDS',raw);table('E3_FINAL_PARENT_MEANS',pg);table('E3_FINAL_METHOD_STATS',ms)
        boot=[b for metric in ['target_CDF_sup','action_regret'] for b in bootstrap(pg,'method',metric,'INCREMENTAL_JOINT_SMOOTHER')]
        table('E3_PAIRED_BOOTSTRAP',boot)
        by=defaultdict(list)
        for row in raw:by[(row['parent'],row['repeat'],row['method'])].append(row)
        timing=[]
        for (p,r,m),items in sorted(by.items()):
            valid=len(items)==3 and all(x['standalone_cpu_seconds'] is not None for x in items)
            cpu=sum(x['standalone_cpu_seconds'] for x in items) if valid else None
            timing.append(dict(parent=p,repeat=r,method=m,status='COMPLETE_CANONICAL' if valid else 'NOT_RUN_MISSING_CANONICAL_COST',
                standalone_cpu_seconds=cpu,common_setup_rhs=11800,common_setup_factorizations=2950,common_setup_CPU=setup['measured_CPU'],
                standalone_total_with_setup=cpu+setup['measured_CPU'] if cpu is not None and setup['measured_CPU'] is not None else None,
                warning='Standalone alternatives share physical preparation; not additive actual campaign charges or free online speedup'))
        tg=parent_means(timing,'method',['standalone_cpu_seconds','standalone_total_with_setup'])
        table('E3_METHOD_CANONICAL_TIMING',timing);table('E3_PARENT_TIMING',tg)
        # Prefix covariance is computed within each repeat, then parent averaged.
        prefix=[]
        for p in range(12):
            for r in range(2):
                for t in range(3):
                    w=refs[p,r,t];mu=w@coords;c=coords-mu;cov=(c*w[:,None]).T@c;mean=w@pose
                    prefix.append(dict(parent=p,repeat=r,time=t+1,occupancy=(w@support).tolist(),pose_mean=mean.tolist(),
                        target_risk=float(w@np.sum((pose-mean)**2,axis=1)),wall_target_cross_covariance=cov[:8,8:].tolist(),
                        wall_group_MAP=int(np.argmax(np.bincount(wall,weights=w))),support_MAP=int(np.argmax(np.bincount(sid,weights=w)))))
        table('E2_ALL_REFERENCE_PREFIXES',prefix)
        prefixstats=[]
        for t in range(3):
            pp=[]
            for p in range(12):
                rr=[r for r in prefix if r['parent']==p and r['time']==t+1]
                pp.append(dict(parent=p,time=t+1,target_risk=float(np.mean([r['target_risk'] for r in rr])),
                    occupancy=np.mean([r['occupancy'] for r in rr],axis=0).tolist(),
                    wall_target_cross_covariance=np.mean([r['wall_target_cross_covariance'] for r in rr],axis=0).tolist()))
            prefixstats.append(dict(time=t+1,parents=12,target_risk_mean=float(np.mean([r['target_risk'] for r in pp])),
                occupancy_mean=np.mean([r['occupancy'] for r in pp],axis=0).tolist(),
                wall_target_cross_covariance_mean=np.mean([r['wall_target_cross_covariance'] for r in pp],axis=0).tolist(),parent_rows=pp))
        # No selected success gallery: every one of the24 frozen noise cases is
        # retained, including early-correct/late-wrong and nonrecovery cases.
        causal=[];truthpath=bank/'offline/FINITE_PRIOR_TRUTHS_AND_NOISE.npz'
        if truthpath.exists():
            with np.load(truthpath,allow_pickle=False) as z:truth=z['truth_index'].copy()
            for p in range(12):
                gt=int(truth[p]);true_pose=pose[gt];true_wall=int(wall[gt]);true_support=int(sid[gt])
                for r in range(2):
                    old,late=refs[p,r,0],refs[p,r,2]
                    ew=int(np.argmax(np.bincount(wall,weights=old)));lw=int(np.argmax(np.bincount(wall,weights=late)))
                    es=int(np.argmax(np.bincount(sid,weights=old)));ls=int(np.argmax(np.bincount(sid,weights=late)))
                    old_mse=float(np.sum((old@pose-true_pose)**2));late_mse=float(np.sum((late@pose-true_pose)**2))
                    hw=weights[p,r,2,hard];ow=weights[p,r,2,ordinary]
                    hv=np.isfinite(hw).all() and np.isclose(hw.sum(),1.);ov=np.isfinite(ow).all() and np.isclose(ow.sum(),1.)
                    row=dict(case_id=2*p+r,parent=p,repeat=r,status='OFFLINE_SAVED_ALL24_CASE',true_wall_group=true_wall,true_support=true_support,
                        early_MAP_wall_group=ew,late_MAP_wall_group=lw,early_MAP_support=es,late_MAP_support=ls,
                        early_wall_wrong=ew!=true_wall,late_wall_correct=lw==true_wall,early_support_wrong=es!=true_support,late_support_correct=ls==true_support,
                        wall_wrong_to_correct=ew!=true_wall and lw==true_wall,support_wrong_to_correct=es!=true_support and ls==true_support,
                        old_pose_MSE=old_mse,late_pose_MSE=late_mse,late_minus_old_pose_MSE=late_mse-old_mse,late_pose_MSE_improved=late_mse<old_mse,
                        hard_late_pose_MSE=float(np.sum((hw@pose-true_pose)**2)) if hv else None,
                        ordinary_late_pose_MSE=float(np.sum((ow@pose-true_pose)**2)) if ov else None,
                        hard_action_regret_under_reference=float(np.sum((hw@pose-late@pose)**2)) if hv else None,
                        hard_minus_ordinary_pose_MSE=float(np.sum((hw@pose-true_pose)**2)-np.sum((ow@pose-true_pose)**2)) if hv and ov else None,
                        interpretation='offline mechanism-pattern flags; not a causal identification proof or online truth-driven selection')
                    causal.append(row)
            table('E3_ALL24_OFFLINE_MECHANISM_CASES',causal)
            cp=parent_means(causal,'status',['wall_wrong_to_correct','support_wrong_to_correct','late_minus_old_pose_MSE','hard_minus_ordinary_pose_MSE'])
            table('E3_MECHANISM_PARENT_MEANS',cp)
        else:report['missing_inputs'].append('E3 offline truths absent; all24 mechanism table NOT_RUN')
        report['E3']=dict(status='SAVED_COMPLETE_INFERENCE_SUMMARIZED',method_count=len(methods),methods=methods,
            expected_final_denominator_per_method=24,expected_prefix_denominator_per_method=72,record_count=len(log),
            method_final_stats=ms,parent_timing_stats=method_summary(tg,'method',['standalone_cpu_seconds','standalone_total_with_setup']),
            ordinary_max_abs_weight_difference_to_reference=ordinary_max,ordinary_weight_check_status='COMPLETE_ALL_PREFIXES' if ordinary_valid else 'NOT_EVALUABLE_INVALID_ORDINARY_WEIGHTS',
            fixed_mechanism_case_ids=list(range(24)),all24_offline_mechanism_cases=causal,prefix_reference_stats=prefixstats,
            bootstrap=boot,local_proxy_scope='ONLINE_TASK_WEIGHTED_PROXY is separately added; original reference-dependent TASK_WEIGHTED remains preserved',
            online_policy_addition_record=read(root/'E2_ONLINE_POLICY_ADDITION.json'),
            no_acquisition_or_scheduling_superiority_assertion=True)
        compact.extend(dict(section='E3',**r) for r in ms)
    else:notrun('E3','Complete saved mask/reference/method panel absent; no offline truths read and no invented causal flags')
    # E4: materialize the entire predeclared384-row panel, including undefined,
    # failed and never-run methods. Numerical summaries do not hide denominators.
    e4=read(root/'results/e4/comparison/RESULTS.json')
    if e4 and e4.get('status') not in ('RUNNING',None):
        settings=[(m,None,None) for m in ['MAP_POINT','PRODUCT_10_COORDINATE_MARGINALS','OCCUPANCY_MIXTURE_OFFLINE_PROJECTION','MOMENT_GAUSSIAN_OFFLINE_PROJECTION']]
        settings += [(m,n,s) for m in ['TEMPERED_SMC','PRIOR_IMPORTANCE'] for n in [256,1024] for s in [0,1,2]]
        lookup={(r['parent_id'],r['repeat_id'],r['method'],r.get('particles'),r.get('seed')):r for r in e4['rows']};rows=[]
        for p in range(12):
            for repeat in range(2):
                for m,n,seed in settings:
                    saved=lookup.get((p,repeat,m,n,seed),{})
                    row=dict(parent=p,repeat=repeat,method=m,particles=n,seed=seed,category=m+(f'_N{n}' if n else ''),
                        status=saved.get('status','NOT_RUN_MISSING_SAVED_ROW'),error=saved.get('error'))
                    for key in ['target_cdf_max','event_error','important_mode_mass_error','action_regret','process_cpu_seconds','all_declared_tolerances_met']:
                        row[key]=num(saved.get(key)) if saved.get('status')=='EVALUATED' else None
                    # Failure/undefined counts remain in the declared denominator.
                    row['evaluated_fraction']=float(saved.get('status')=='EVALUATED')
                    row['all_tolerances_declared_denominator_fraction']=float(saved.get('status')=='EVALUATED' and bool(saved.get('all_declared_tolerances_met')))
                    rows.append(row)
        table('E4_ALL384_INCLUDING_FAILURES',rows);parents=[]
        metrics=['target_cdf_max','event_error','important_mode_mass_error','action_regret','process_cpu_seconds','evaluated_fraction','all_tolerances_declared_denominator_fraction']
        for category in sorted(set(r['category'] for r in rows)):
            block=[r for r in rows if r['category']==category]
            parents.extend(parent_means(block,'category',metrics,3 if block[0]['particles'] else 1))
        table('E4_PARENT_SEEDS_WITHIN_REPEAT',parents)
        stats=method_summary(parents,'category',metrics)
        for s in stats:s['all_status_denominators']=dict(Counter(r['status'] for r in rows if r['category']==s['category']))
        boot=[b for m in ['target_cdf_max','action_regret'] for b in bootstrap(parents,'category',m,'MAP_POINT')]
        table('E4_PAIRED_BOOTSTRAP',boot)
        report['E4']=dict(status='TERMINAL_RAW_COMPARISONS_SUMMARIZED',input_status=e4['status'],declared_rows=384,
            all_status_denominators=dict(Counter(r['status'] for r in rows)),method_stats=stats,bootstrap=boot,
            scope='Explicit finite prior only; exact-moment offline projections are not deployable inference algorithms')
        compact.extend(dict(section='E4',**r) for r in stats)
    else:notrun('E4','Terminal E4 comparison panel absent')
    # E5 truth access occurs only after COMPLETE online policy decisions. The MC
    # preposterior expectation, realized posterior risk and realized GT loss are
    # deliberately different columns, never substituted for one another.
    e5root=root/'results/e5/finite_sensing';manifest=read(e5root/'MANIFEST.json');online=read(e5root/'ONLINE_RESULTS.json')
    if manifest and manifest.get('status')=='COMPLETE' and online:
        policies=['Random','Uniform','EIG','SingularValue','PrincipalAngleSum','GeneralizedRetention','TaskEVSI']
        lookup={(r['parent'],r['repeat'],r['policy']):r for r in online};rows=[];mc_pairs=[]
        for p in range(12):
            for repeat in range(2):
                case=e5root/f'parent_{p:02d}_repeat_{repeat}';selection=read(case/'ACTION_SELECTION.json')
                gtrows=read(e5root/f'offline/GT_parent_{p:02d}_repeat_{repeat}.json') or []
                gt={r['policy']:r for r in gtrows}
                draws=np.asarray(selection.get('predictive',{}).get('draw_posterior_risk',[]),float) if selection else np.empty((0,0))
                mcvalid=draws.shape==(4,32) and np.isfinite(draws).all()
                for policy in policies:
                    saved=lookup.get((p,repeat,policy));action=int(saved['action']) if saved else None
                    row=dict(parent=p,repeat=repeat,policy=policy,status='SAVED' if saved else 'NOT_RUN_MISSING_POLICY_ROW',action=action,
                        realized_posterior_bayes_risk=num(saved.get('target_bayes_risk')) if saved else None,
                        realized_GT_squared_error=num(gt.get(policy,{}).get('target_squared_error')),
                        selection_cpu_seconds_standalone=num(saved.get('selection_cpu_seconds_standalone')) if saved else None,
                        inference_cpu_seconds=num(saved.get('inference_cpu_seconds')) if saved else None,
                        preposterior_MC_expected_risk=float(draws[action].mean()) if mcvalid and action is not None else None,
                        preposterior_MC_SE=float(draws[action].std(ddof=1)/np.sqrt(32)) if mcvalid and action is not None else None,
                        MC_draws=32 if mcvalid else None,MC_status='SAVED32_CONDITIONAL_DESCRIPTIVE' if mcvalid else 'NOT_RUN_MISSING_DRAW_ARRAY',
                        common_setup_rhs=11800,common_setup_factorizations=2950,common_setup_CPU=setup['measured_CPU'])
                    rows.append(row)
                if mcvalid and selection:
                    for baseline in ['EIG','Uniform']:
                        a=selection['choices']['TaskEVSI'];b=selection['choices'][baseline];delta=draws[a]-draws[b]
                        mc_pairs.append(dict(parent=p,repeat=repeat,method='TaskEVSI',baseline=baseline,draws=32,
                            expected_risk_difference=float(delta.mean()),paired_MC_SE=float(delta.std(ddof=1)/np.sqrt(32)),
                            scope='common predictive draws; selection-biased conditional MC uncertainty, not policy superiority certificate'))
        metrics=['preposterior_MC_expected_risk','realized_posterior_bayes_risk','realized_GT_squared_error','selection_cpu_seconds_standalone','inference_cpu_seconds']
        pg=parent_means(rows,'policy',metrics);stats=method_summary(pg,'policy',metrics)
        for s in stats:
            s['declared_policy_noise_cases']=24
            s['all_status_denominators']=dict(Counter(r['status'] for r in rows if r['policy']==s['policy']))
        boots=[]
        # Positive TaskEVSI-baseline differences mean higher loss/risk; show both
        # baselines and every pair, without choosing the favorable comparison.
        for baseline in ['EIG','Uniform']:
            for metric in metrics[:3]:
                boots.extend(r for r in bootstrap(pg,'policy',metric,baseline) if r['method']=='TaskEVSI')
        table('E5_ALL168_POLICY_ROWS',rows);table('E5_PARENT_METRICS',pg);table('E5_TASK_EIG_UNIFORM_BOOTSTRAP',boots);table('E5_SAVED32_MC_PAIRED_UNCERTAINTY',mc_pairs)
        report['E5']=dict(status='COMPLETE_SAVED_POLICY_PANEL_SUMMARIZED',expected_rows=168,actual_saved_online_rows=len(online),
            all_status_denominators=dict(Counter(r['status'] for r in rows)),policy_stats=stats,
            TaskEVSI_vs_EIG_and_Uniform_bootstrap=boots,saved32_MC_uncertainty_rows=len(mc_pairs),
            quantities_separate=['MC preposterior expected risk','realized posterior Bayes risk','offline realized GT squared error'],
            scope='conditional finite-prior pilot only;32-draw conditional MCSE is not a calibrated selection/superiority test',acquisition_superiority_asserted=False)
        compact.extend(dict(section='E5',**r) for r in stats)
    else:notrun('E5','COMPLETE saved online sensing panel absent; no offline GT inputs read')
    # Optional ordinary optimizer diagnostic: four clusters only, not part of the
    # twelve-parent method bootstrap and not a matched-runtime experiment.
    from polish_artifacts import load_polish
    polish=load_polish(root)
    table('A0_POLISH_ALL_ACTUAL_TRIALS',polish['trials'])
    optional=dict(status=polish['status'],offline_status=polish.get('offline_status','NOT_RUN'),
        expected_cases=8,independent_parent_denominator=4,fixed_parent_ids=[0,3,6,9],
        bootstrap='NOT_RUN_FOUR_CLUSTER_EXPLORATORY_NO_FORMAL_CI',
        holdout_scope=polish['holdout_scope'],matched_speed_claim=False,
        original_three_primary_arms_unchanged=True,A0_cap_reset=False,setup_scope=polish['setup_scope'],
        budget_start_snapshot=polish.get('budget_start_snapshot'),budget_end_snapshot=polish.get('budget_end_snapshot'),
        snapshot_warning=polish.get('budget_snapshot_CPU_warning'),actual_trial_count=len(polish['trials']))
    if polish.get('offline_status')=='SAVED_TERMINAL_OFFLINE_PAIR_SUMMARY':
        pairs=[dict(r,diagnostic='ordinary_adjoint_polish') for r in polish['paired_rows']]
        metrics=[field+suffix for field in ['objective','training_noise_energy','holdout_noise_energy','pose_error_lambda0',
            'pose_squared_error_lambda0_squared','material_eps_contrast_L2','material_ell_L2'] for suffix in ['_before','_after','_change']]
        parents=parent_means(pairs,'diagnostic',metrics)
        paired_stats=dict(independent_parent_denominator=4,expected_noise_cases=8,
            all_status_denominators=polish['all_status_denominators'],
            starting_arm_counts=dict(Counter(r['initial_method'] for r in pairs)))
        for metric in metrics:
            values=[r[metric] for r in parents if r[metric+'_complete_repeats']==2 and r[metric] is not None]
            paired_stats[metric+'_mean']=float(np.mean(values)) if values else None
            paired_stats[metric+'_complete_parents']=len(values)
        table('A0_POLISH_PAIRED_BEFORE_AFTER',pairs);table('A0_POLISH_FOUR_PARENT_MEANS',parents)
        optional.update(paired_stats=paired_stats,all_eight_saved_pairs=pairs,
            interpretation='ordinary known-intrinsic optimization diagnostic; exposed-development reevaluation; no new GO/posterior/speed claim')
        compact.append(dict(section='A0_ORDINARY_ADJOINT_POLISH',status=polish['status'],**paired_stats))
    else:
        optional['reason']=polish.get('reason') or 'All eight online and terminal offline outputs not complete'
        compact.append(dict(section='A0_ORDINARY_ADJOINT_POLISH',status='NOT_RUN',reason=optional['reason']))
    report['A0_ORDINARY_ADJOINT_POLISH']=optional
    table('SUMMARY_COMPACT',compact)
    # Required sibling CSV contains compact method/policy rows; detailed all-case
    # denominators and fixed IDs are retained in the separate CSV tables.
    keys=list(dict.fromkeys(k for r in compact for k in r)) or ['section','status']
    with csvtarget.open('w',newline='') as f:
        writer=csv.DictWriter(f,keys);writer.writeheader()
        for r in compact:writer.writerow({k:json.dumps(v,default=scalar) if isinstance(v,(dict,list,np.ndarray)) else v for k,v in r.items()})
    report['tables_directory']=str(tables);report['summary_csv']=str(csvtarget)
    report['cost_note']='Preparation/summary CPU is charged by wrapper; zeroPDE is not zeroCPU. Common bank charged once, standalone alternatives not additive.'
    target.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False,default=scalar)+'\n')
    print(json.dumps(dict(path=str(target),E3=report['E3']['status'],E4=report['E4']['status'],E5=report['E5']['status'])),flush=True)


if __name__=='__main__':main()
