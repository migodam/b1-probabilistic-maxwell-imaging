"""Offline stored-array failure diagnostic; no forward/solver/GT-based tuning."""
from pathlib import Path
import json,time,importlib.util,itertools,collections
import numpy as np
ROOT=Path(__file__).resolve().parents[2];REPO=ROOT.parents[1]


def load(path):return json.loads(Path(path).read_text())
def save(path,obj):Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')


def main():
    wall,cpu=time.perf_counter(),time.process_time()
    score=load(ROOT/'RESULTS/confirmation_scoring/SCORES.json')
    sealed_score=load(ROOT/'RESULTS/confirmation_scoring/RECEIPT.json')
    if sealed_score['status']!='SCORED_ALL66_CONFIRMATION_TASKS':raise ValueError('All66 scored confirmation tasks required')
    cases=[f'r9_confirmation_{i:04d}' for i in range(6)]
    if sealed_score['cases']!=cases or set(r['case_id'] for r in score['rows'])!=set(cases):raise ValueError('Unexpected truth case set')
    meta={r['case_id']:r for r in load(ROOT/'sealed_truth/confirmation/score_metadata.json')['records']}
    pub={r['case_id']:r for r in load(ROOT/'data_public/confirmation/manifest.json')['records']}
    # Existing pure scalar prior definition; does not import a physical solver.
    source=REPO/'research/research_max_struct_r7/experiments/score_primary.py'
    spec=importlib.util.spec_from_file_location('r9_failure_fixed_prior',source)
    metrics=importlib.util.module_from_spec(spec);spec.loader.exec_module(metrics)
    model_rows=[r for r in score['rows'] if r['origin']!='FROZEN_TWO_BRANCH_PORTFOLIO']
    if len(model_rows)!=66:raise ValueError('Expected66 actual algorithm endpoints plus separately labelled portfolios')
    records=[];materials={};truth_paths=[];atol=1e-8
    for case in cases:
        truth_path=ROOT/'sealed_truth/confirmation'/(case+'.npz');truth_paths.append(str(truth_path))
        with np.load(truth_path,allow_pickle=False) as z:m=z['m'].ravel().copy();clean=z['Y_clean'].copy()
        public_path=REPO/pub[case]['public_files'][1]
        with np.load(public_path,allow_pickle=False) as z:Y=z['Y'].copy();sigma=float(z['sigma'])
        prior=metrics.fixed_prior(m);data=.5*float(np.vdot(clean-Y,clean-Y).real/np.vdot(Y,Y).real)
        feasible=data+prior
        rr=[]
        for r in [x for x in score['rows'] if x['case_id']==case]:
            with np.load(r['prediction_path'],allow_pickle=False) as z:a=z['m'].ravel().copy()
            materials[r['task_id']]=a
            gap=float(r['original_Phi']-feasible);better=gap>1e-12+1e-8*abs(feasible)
            stop=r.get('corrector_stop') or ''
            kind='REPORTED_PGTOL' if 'PROJECTED GRADIENT' in stop.upper() else 'BUDGET_CAP' if 'LIMIT' in stop.upper() or 'MAXIMUM' in stop.upper() else 'REPORTED_RELATIVE_F' if 'RELATIVE REDUCTION' in stop.upper() else 'OTHER'
            rr.append(dict(task_id=r['task_id'],method=r['method'],dpcs_rule=r.get('dpcs_rule'),origin=r['origin'],Phi=float(r['original_Phi']),minus_known_feasible_GT_Phi=gap,known_lower_feasible_GT_point=bool(better),Phi_over_GT_Phi=float(r['original_Phi']/feasible),material_nrmse=r['material_nrmse'],material_success=r['material_image_success'],stop_category=kind,corrector_stop=stop,corrector_success=r['corrector_success'],corrector_iterations=r['corrector_iterations'],active_lower_fraction=float(np.mean(abs(a-.02)<=atol)),active_upper_fraction=float(np.mean(abs(a-20.)<=atol)),prediction_path=r['prediction_path']))
        records.append(dict(case_id=case,amplitude=meta[case]['amplitude'],GT_feasible_Phi=feasible,GT_data_term=data,GT_prior=prior,GT_relative_residual=float(np.linalg.norm(clean-Y)/np.linalg.norm(Y)),GT_whitened_residual=float(np.linalg.norm(clean-Y)/sigma),GT_legal_box=bool(np.min(m)>=.02 and np.max(m)<=20.),same_discretization=True,stored_clean_field_is_feasible_prediction='GenerationandInverseboth24; inheritedsameG/E/S; nofieldrecomputation',rows=rr))
    star=[r for r in model_rows if r['case_id']==cases[4] and r['method'] in ['cie_single','cie_dual','global_dual','a_plus']]
    if len(star)!=8:raise ValueError('Expected4methodsx2noise-rules for star')
    pairs=[]
    for a,b in itertools.combinations(star,2):
        x,y=materials[a['task_id']],materials[b['task_id']]
        pairs.append(dict(task_a=a['task_id'],task_b=b['task_id'],same_method=a['method']==b['method'],relative_material_L2_symmetric=float(np.linalg.norm(x-y)/max(np.linalg.norm(x),np.linalg.norm(y))),max_absolute_material_difference=float(np.max(abs(x-y))),bitwise_equal=bool(np.array_equal(x,y))))
    star_summary=dict(case_id=cases[4],endpoints=8,definition='norm(ma-mb)/max(norm(ma),norm(mb)); storedfullfloatmaterials, notPhi-onlycomparison',all_pair_count=len(pairs),maximum_pair_relative_L2=max(p['relative_material_L2_symmetric'] for p in pairs),maximum_pair_absolute_difference=max(p['max_absolute_material_difference'] for p in pairs),all_bitwise_equal=all(p['bitwise_equal'] for p in pairs),fixed_vs_noise_same_method=[p for p in pairs if p['same_method']],all_pairs=pairs)
    summaries=[]
    for method,rule in sorted(set((r['method'],str(r.get('dpcs_rule'))) for r in score['rows'])):
        rr=[r for r in score['rows'] if r['method']==method and str(r.get('dpcs_rule'))==rule]
        summaries.append(dict(method=method,rule=rule,n=len(rr),successes=sum(r['material_image_success'] for r in rr),independent_algorithm_runs=rr[0]['origin']!='FROZEN_TWO_BRANCH_PORTFOLIO'))
    failures=[]
    for rec in [r for r in records if r['case_id'] in [cases[2],cases[4]]]:
        rr=[r for r in rec['rows'] if r['origin']!='FROZEN_TWO_BRANCH_PORTFOLIO']
        failures.append(dict(case_id=rec['case_id'],GT_feasible_Phi=rec['GT_feasible_Phi'],actual_endpoints=len(rr),all_material_fail=all(not r['material_success'] for r in rr),all_have_lower_known_feasible_point=all(r['known_lower_feasible_GT_point'] for r in rr),minimum_final_Phi=min(r['Phi'] for r in rr),maximum_final_Phi=max(r['Phi'] for r in rr),minimum_Phi_over_GT=min(r['Phi_over_GT_Phi'] for r in rr),stop_counts=dict(collections.Counter(r['stop_category'] for r in rr)),lower_active_fraction_range=[min(r['active_lower_fraction'] for r in rr),max(r['active_lower_fraction'] for r in rr)],upper_active_fraction_range=[min(r['active_upper_fraction'] for r in rr),max(r['active_upper_fraction'] for r in rr)]))
    result=dict(status='OFFLINE_POST_ALL66_SEAL_FAILURE_DIAGNOSTIC',actual_algorithm_endpoints=66,derived_portfolio_rows=12,GT_witness_source='StoredY_cleanandpublicYwithR7fixed_prior; zero newforward',prior_source=str(source),active_bound_atol=atol,active_bound_rtol=0.,known_better_tolerance='gap>1e-12+1e-8abs(GTPhi)',records=records,star_0004_endpoint_comparison=star_summary,common_failures=failures,algorithm_success_summary=summaries,qualification='GTisafeasiblewitnessnotaglobaloptimum. Lowerknownfeasiblepointshowsnotattainingthatvalue, notuniqueness/nonuniqueness/localminimum/globalbarrier. Reportedstopmaybelongcorrectorendpointnotselectedincumbent; no gradients recomputed. Noonlinechangesornewtesttuning.',new_physical_calls=0,new_LU=0,new_complex_RHS=0)
    save(ROOT/'audit/CONFIRMATION_FAILURE_DIAGNOSTIC.json',result)
    lines=['# 六例确认：离线失败诊断','',
        '全部66条真实算法端点封存/评分后，仅读取六个确认真值、已保存Y_clean及公开Y。用原R7 fixed_prior计算 `Phi_GT=.5||Y_clean−Y||²/||Y||²+R(m_GT)`，没有新增forward/adjoint/LU、反演或调参。12个两支路portfolio只是派生选择行，不计作额外独立运行。','',
        '| 对象 | 已知可行GT Phi | data项 | 原prior |','|---|---:|---:|---:|']
    for r in records:lines.append(f"| {r['case_id'][-4:]} | {r['GT_feasible_Phi']:.8g} | {r['GT_data_term']:.8g} | {r['GT_prior']:.8g} |")
    lines+=['','## 两个共同失败并非没有已知更低可行点','']
    for r in failures:
        lines.append(f"- {r['case_id'][-4:]}：11/11材料失败且11/11有严格更低的GT可行Phi。最终Phi范围{r['minimum_final_Phi']:.6g}–{r['maximum_final_Phi']:.6g}，最小仍为GT可行Phi的{r['minimum_Phi_over_GT']:.2f}倍；stop计数{r['stop_counts']}。下界活动比例范围{r['lower_active_fraction_range'][0]:.3f}–{r['lower_active_fraction_range'][1]:.3f}，上界{r['upper_active_fraction_range'][0]:.3f}–{r['upper_active_fraction_range'][1]:.3f}（绝对容差1e−8，无相对容差）。")
    lines+=['','0002 annulus中，fixed四臂和noise-single均CAP；noise-Dual/Global/A+、direct、profile及weak-handoff虽然报告PGTOL，也仍未达到更低的已知可行值。0004 star的八条fixed/noise CIE/A+及profile均报告PGTOL，direct与weak-handoff报告relative-F；它们仍是错误材料。此处沿用receipt的停止声明，不声称最终选中incumbent本身已有独立投影梯度认证。','',
        f"0004的八条CIE/A+保存材料也近乎相同，不能只凭Phi接近下结论：28个成对比较的最大对称相对L2差为{star_summary['maximum_pair_relative_L2']:.4g}，最大单cell绝对差{star_summary['maximum_pair_absolute_difference']:.4g}；逐位相同={star_summary['all_bitwise_equal']}。每方法fixed↔noise的独立材料差保留在JSON。近同端点不是严格共同吸引盆地、局部极小或唯一性证明。",'',
        '## 对预期主张的限定','',
        '- 噪声规则在本六例将single、Dual、Global和A+都从2/6提高到4/6。它支持共同DPCS合同修正的有限收益，不支持A+独有增量或“广泛解决高contrast”；0002/0004仍共同失败。',
        '- A+旧弯条等DEV正例必须保留，但新六例没有超过ordinary noise-aware CIE的4/6。不能把共同修复计为A+的新算法优势。',
        '- 已知WRI/LSIE profile continuation在DEV40ac提供互补正例，但新六例仅2/6，weak-handoff也是2/6；两支路portfolio4/6仍不超过单noise-aware ordinary CIE。DEV互补不等于新确认增量，组合还须支付两支路。',
        '- 可行GT目标低于错误端点，只说明这些端点尚未达到一个已知较低值；不能推出材料唯一/非唯一、已认证局部极小、全局障碍或信息不足。GT自身也不是噪声/prior目标的已证全局最优。',
        '- 以上没有建立全局凸化、强原创算法或普适SOTA。确认暴露后不得以这些真值继续选择rank、阈值、延拓日程或停止点。','',
        '完整六例逐算法Phi/stop/活动边界、0004材料差异见 `audit/CONFIRMATION_FAILURE_DIAGNOSTIC.json`。本诊断仅输出原目标可行见证与解释边界，不改任何算法、预测、评分或报告主表。']
    (ROOT/'delegated/confirmation_review/FAILURES_ZH.md').write_text('\n'.join(lines)+'\n')
    receipt=dict(status='COMPLETED_OFFLINE_FAILURE_DIAGNOSTIC',CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,truth_paths_opened=truth_paths,authorized_scope='SixconfirmationcasesAFTERALL66sealed/scored; OFFLINEONLY',new_physical_calls=0,new_LU=0,new_complex_RHS=0,inverse_runs=0,algorithm_changes=0,new_test_tuning=0,files=['audit/CONFIRMATION_FAILURE_DIAGNOSTIC.json','delegated/confirmation_review/FAILURES_ZH.md'],cost_scope='Storedarrayarithmeticandwriting; importsnotfullymetered; no newphysics')
    save(ROOT/'audit/CONFIRMATION_FAILURE_DIAGNOSTIC_COST.json',receipt)
    print(json.dumps(dict(common_failures=failures,star_relative_max=star_summary['maximum_pair_relative_L2'],same_method_differences=star_summary['fixed_vs_noise_same_method'],CPU_seconds=receipt['CPU_seconds'])))


if __name__=='__main__':main()
