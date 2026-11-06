"""Aggregate two authorized corrected fresh starts; no physical work here."""
from pathlib import Path
import json,time
import numpy as np
from fresnel_report import main as report_main,ROOT,OUT

def main():
 start=time.process_time();checks={}
 for arm in ['cie_dual','a_plus']:
  old=OUT/f'FoamTwinDielTM_8GHz_n60_{arm}_v3';new=OUT/f'FoamTwinDielTM_8GHz_n60_{arm}_dev_air0'
  c0=json.loads((old/'FROZEN_CONFIG.json').read_text());c1=json.loads((new/'FROZEN_CONFIG.json').read_text())
  q0=dict(c0['config']);q1=dict(c1['config']);assert q0.pop('material_box')==[.02,20.] and q1.pop('material_box')==[0.,20.]
  q0.pop('progress_path');q1.pop('progress_path');checks[arm+'_all_other_auxiliary_configuration_equal']=q0==q1
  for k in ['case','frequency_GHz','side','width_m','arm','stopping_contract','corrector','model','initialization','training_channels','withheld_channels','heldout_selection']:checks[arm+'_'+k]=c0[k]==c1[k]
  checks[arm+'_empty_target_calibration_equal']=json.loads((old/'INCIDENT_CALIBRATION.json').read_text())==json.loads((new/'INCIDENT_CALIBRATION.json').read_text())
  a=np.load(old/'derived.npz');b=np.load(new/'derived.npz')
  for k in ['centers','train_mask','withheld_mask']:checks[arm+'_'+k]=np.array_equal(a[k],b[k])
  checks[arm+'_fresh_initial_point1']=bool(np.all(b['initial_material']==.1))
  aux0=json.loads((old/'AUXILIARY_RECEIPT.json').read_text());aux1=json.loads((new/'AUXILIARY_RECEIPT.json').read_text());checks[arm+'_same_current_PR_mode']=aux0['sourcewise_cg']==aux1['sourcewise_cg']
 assert all(checks.values()),checks
 report_main()
 corrected=json.loads((OUT/'AIR0_CORRECTED_COMPARISON_SUMMARY.json').read_text())['runs'];assert len(corrected)==2
 cost=json.loads((OUT/'COST_AND_SELECTION_AUDIT.json').read_text());base=json.loads((OUT/'FINAL_QA_BASE23.json').read_text())
 report=ROOT/'FRESNEL_REPORT.md';text=report.read_text();text=text.replace('现另完成1条事后air0可行域诊断，共23条实际执行。','先另完成1条事后air0可行域诊断、再完成2条fresh-start修正；原23条保留，现共25条实际执行。',1)
 lines=['## 两条fresh-start空气下界修正','',
 '单次续跑已经表明主合同的正背景约束显著影响预测，因此父线程进一步授权且仅执行Dual-CIE/A+两个最小修正。两条均从原初值.1重新开始，沿用60²/8GHz、原source mask及空目标校准、DPCS rank5、mf6/10/16、material beta6/2/.5、固定current beta6、原prior、每阶段最多1000（总3000）辅助步及原1000 accepted/4000eval corrector。仅将辅助/readout/corrector的材料下界由.02改为0，上界仍20。Dual仍为原sourcewise PR，A+仍为原global PR；没有另加epsilon floor。', '',
 '这是已查看主结果后的DEV_MODEL_CORRECTION，不是新blind比较或A+胜利声明；6/8GHz留出均此前已暴露，6GHz只在固定8GHz终点作事后预测。原23条及其费用/QA封存保留，新增两条不覆盖原五臂表。', '',
 'A+零端点满足 `z=m/(1+(beta-Gii)m)`、`z(0)=0`、`z′(0)=1`。实际8GHz60² `Gii=.0864483905+.0771235448i`，三个beta均满足Re(beta−Gii)>0，故实区间[0,20]分母模≥1；端点、候选评价、原prior及fallback均合法。详细数学/代码审计见 `delegated/fresnel/AIR0_MOBIUS_DOMAIN_AUDIT.json`。','']
 for r in corrected:
  md=r['material_diagnostics'];bg=md['outer_roi_band_distribution'];whole=md['whole_roi_distribution'];c=r['corrector'];counts=r['auxiliary']['operation_counts'];old=json.loads((OUT/f"FoamTwinDielTM_8GHz_n60_{r['arm']}_v3/RECEIPT.json").read_text())
  lines.append(f"- {r['arm']}：8GHz train/heldout {r['train_relative_complex_residual']:.6f}/{r['withheld_all_relative_complex_residual']:.6f}（原floor.02为{old['train_relative_complex_residual']:.6f}/{old['withheld_all_relative_complex_residual']:.6f}）；事后6GHz heldout {r['posthoc_6GHz_metrics']['heldout_all']:.6f}；peak {md['max_contrast']:.6f}；统一旧floor-excess有效面积 {md['foreground_effective_area_m2']:.8f}m²，原始材料有效面积 {md['effective_area_raw_m2']:.8f}m²。全ROI零下界比例 {100*whole['fraction_at_declared_lower_bound']:.2f}%，低于旧.02比例 {100*whole['fraction_below_old_0p02']:.2f}%；外85mm圈mean/median/95%为 {bg['mean']:.6f}/{bg['quantiles']['0.5']:.6f}/{bg['quantiles']['0.95']:.6f}。实际辅助 {counts['PR_CG_steps']}步，训练物理目标选择stage={r['selected_auxiliary_stage']}；各辅助stage停止为{[v['status'] for v in r['auxiliary_stage_stops']]}；corrector {c['iterations']} accepted/{c['evaluations']}eval，success={c['success']}（{c['message']}）；CPU {r['cpu_seconds']:.6f}s，其中辅助及准备{r['auxiliary_cpu_seconds']:.6f}s。")
 lines+=['','分布/峰值/有效面积不等于GT精度，外圈只作描述而非支持约束。优化器CAP/收敛及离散state精度按真实receipt保留；不同选择stage或更长动力学不是恢复认证。没有继续扩展算法、扫参、训练或反演。','','![Fresh-start corrected comparison](FIGURES/fresnel/air0_corrected_comparison.png)','','']
 if '## 两条fresh-start空气下界修正' not in text:text=text.replace('## 计算和失败','\n'.join(lines)+'## 计算和失败',1)
 needle='23条合计9112.502822 CPU s；';text=text.replace(needle,f"原23条合计9112.502822 CPU s；新增两条fresh-start修正{cost['air0_fresh_start_DEV_MODEL_CORRECTION_CPU_seconds']:.6f} CPU s，25条总计{cost['successful_all_measured_tasks_CPU_seconds']:.6f} CPU s；",1)
 if '`experiments/fresnel_foam_air0_correction.py`' not in text:text+='\n两条fresh-start修正入口为 `experiments/fresnel_foam_air0_correction.py --arm cie_dual` / `--arm a_plus`（固定默认dev_air0、60²，拒绝覆盖）；`experiments/fresnel_foam_air0_finalize.py` 只汇总和核对已完成端点，不运行物理反演。\n'
 report.write_text(text)
 receipts=json.loads((OUT/'SUMMARY.json').read_text())+[json.loads((OUT/'FoamTwinDielTM_8GHz_n60_independent_real_nie/RECEIPT.json').read_text())]+json.loads((OUT/'MULTIFREQUENCY_SUMMARY.json').read_text())+[json.loads((OUT/'FoamTwinDielTM_air0_control_8GHz_dev/RECEIPT.json').read_text())]+corrected
 forbidden=[];npzs=list(OUT.rglob('*.npz'))
 for p in npzs:
  with np.load(p) as z:
   forbidden += [str(p.relative_to(ROOT))+':'+k for k in z.files if k in {'Y','Y_train','measured_Y','incident_measured','rawI','scattered_measured'}]
 figs=['material_endpoints.png','neural_warm_starts.png','heldout_predictions.png','independent_nie.png','frequency_followup.png','air0_control.png','air0_corrected_comparison.png'];missing=[f for f in figs if not (ROOT/'FIGURES/fresnel'/f).exists()]
 total=sum(r['cpu_seconds'] for r in receipts)
 qa=dict(status='COMPLETED_ARTIFACT_QA',expected_completed_inversions=25,actual_completed_inversions=len(receipts),original_23_preserved=base,all_receipt_statuses_completed=all(r['status'].startswith('COMPLETED') for r in receipts),successful_CPU_seconds=total,aggregate_CPU_matches_receipts=bool(np.isclose(total,cost['successful_all_measured_tasks_CPU_seconds'],rtol=0,atol=1e-8)),fresh_start_only_bound_change_checks=checks,forbidden_raw_measurement_keys=forbidden,missing_derived_figures=missing,public_npz_files=len(npzs),image_QA='Originalsixfiguresviewed; newcorrectedcomparisonawaitsvisualreview',GT_metrics=base['GT_metrics'],NIE_scope=base['NIE_scope'],cost_scope='Original23preserved+twofreshstartDEV_MODEL_CORRECTION; knownfailureCPU>=405.52s+unclosedfirstlayoutCPU; sharedhost/noisolatedtimingclaim;2GHzAuncached')
 assert len(receipts)==25 and qa['all_receipt_statuses_completed'] and qa['aggregate_CPU_matches_receipts'] and not forbidden and not missing
 (OUT/'FINAL_QA.json').write_text(json.dumps(qa,indent=2)+'\n')
 (ROOT/'delegated/fresnel/AIR0_CORRECTED_FINALIZATION_RECEIPT.json').write_text(json.dumps(dict(status='COMPLETED_METADATA_QA_AWAITING_FIGURE_REVIEW',additional_physical_runs=2,maximum_physics_processes=2,BLAS_threads_per_process=1,old_runs_repeated=False,configuration_checks=checks,finalization_CPU_seconds=time.process_time()-start),indent=2)+'\n')
 print(json.dumps(dict(count=25,total_CPU=total,additional_CPU=cost['air0_fresh_start_DEV_MODEL_CORRECTION_CPU_seconds'],runs=corrected)))

if __name__=='__main__':main()
