"""Finalize one paid air0 diagnostic; no numerical inversion or observation fits."""
from pathlib import Path
import json,time
import numpy as np
from fresnel_report import main as report_main,ROOT,OUT

def main():
 start=time.process_time()
 old=OUT/'FoamTwinDielTM_multifrequency_control_8';new=OUT/'FoamTwinDielTM_air0_control_8GHz_dev'
 cfg0=json.loads((old/'FROZEN_CONFIG.json').read_text());cfg1=json.loads((new/'FROZEN_CONFIG.json').read_text())
 keys=['training_frequencies_GHz','validation_frequencies_GHz','shared_real_material','side','width_m','objective','prior','anchor','anchor_selection','maximum_accepted_iterations','maximum_function_evaluations','corrector','mask','training_channels_per_frequency','heldout_channels_per_frequency','heldout_selection','frequency_sweep','illumination','physics']
 checks={k:cfg0[k]==cfg1[k] for k in keys};assert all(checks.values());assert cfg0['box']==[.02,20.] and cfg1['box']==[0.,20.]
 a=np.load(old/'derived.npz');b=np.load(new/'derived.npz');checks.update({k:np.array_equal(a[k],b[k]) for k in ['anchor_m','centers','train_mask','withheld_mask']});assert all(checks.values())
 assert np.array_equal(b['anchor_m'],np.load(ROOT/cfg1['anchor'])['m'])
 for f in [6,8]:
  checks[f'identical_empty_target_calibration_{f}GHz']=json.loads((old/f'INCIDENT_CALIBRATION_{f}GHz.json').read_text())==json.loads((new/f'INCIDENT_CALIBRATION_{f}GHz.json').read_text())
 assert all(checks.values())
 report_main()
 s=json.loads((OUT/'AIR0_CONTROL_SUMMARY.json').read_text());r=s['air0_control'];o=s['old_floor_control'];d=r['material_diagnostics'];bg=d['outer_roi_band_distribution'];cost=json.loads((OUT/'COST_AND_SELECTION_AUDIT.json').read_text());base=json.loads((OUT/'FINAL_QA_BASE22.json').read_text())
 report=ROOT/'FRESNEL_REPORT.md';text=report.read_text()
 text=text.replace('共22条实际执行。','原22条实际执行已保留；现另完成1条事后air0可行域诊断，共23条实际执行。',1)
 marker='## 计算和失败'
 section=f'''## 空气下界的单次事后诊断

父线程在查看完整材料图后，针对主合同 `m≥.02` 在5.33λ大ROI中排除已知空气这一物理问题，授权且仅执行一次 `DEV_MODEL_DIAGNOSTIC`。它从多频对照使用的同一个已付费8GHz Dual锚点出发，保留原模型、各频空目标校准、真实source mask、prior、上界20、L-BFGS-B及300 accepted/1000 eval，唯一改变是材料下界由.02变0。原.02对照直接复用；没有新增辅助/训练、扫参或观测。两频留出此前已暴露，故不能把这一事后模型诊断当作新封存测试或A+优势。

8GHz train/heldout由 {o['metrics']['8']['train']:.6f}/{o['metrics']['8']['heldout_all']:.6f} 变为 {r['metrics']['8']['train']:.6f}/{r['metrics']['8']['heldout_all']:.6f}；未用于反演的6GHz heldout由 {o['metrics']['6']['heldout_all']:.6f} 变为 {r['metrics']['6']['heldout_all']:.6f}。峰值由 {o['material_diagnostics']['max_contrast']:.4f} 变为 {d['max_contrast']:.4f}。统一按 `max(m-.02,0)` 定义的excess有效面积由 {o['material_diagnostics']['foreground_effective_area_m2']:.8f} 变为 {d['foreground_effective_area_m2']:.8f} m²；这只是分布统计，不是GT面积或材料精度。

固定外圈 `max(|x|,|y|)≥85mm` 仅用于事后描述，不参与优化，也不等于真实物体支持。air0外圈均值/中位数/95%分位为 {bg['mean']:.6f}/{bg['quantiles']['0.5']:.6f}/{bg['quantiles']['0.95']:.6f}；全ROI有 {100*d['whole_roi_distribution']['fraction_at_declared_lower_bound']:.2f}% 在零下界、{100*d['whole_roi_distribution']['fraction_below_old_0p02']:.2f}% 低于原.02。完整原对照及新端点的分布/原始有效面积/下界比例保存在 `RESULTS/fresnel/AIR0_CONTROL_SUMMARY.json`。

新运行实际 {r['corrector']['iterations']} accepted/{r['corrector']['evaluations']} eval，optimizer success={r['corrector']['success']}（{r['corrector']['message']}），新增 {r['cpu_seconds']:.6f} CPU s。它隔离了原正背景约束这一可行域因素；预测和图像变化仍不能证明绝对材料恢复，也不能排除入射/收发端模型误差或其他不可辨识性。未再做任何搜索或追加反演。

![Known-air lower-bound diagnostic](FIGURES/fresnel/air0_control.png)

'''
 if '## 空气下界的单次事后诊断' not in text:text=text.replace(marker,section+marker,1)
 old_cost='22条成功执行累计8717.167454 CPU s；'
 text=text.replace(old_cost,f"原22条成功执行累计8717.167454 CPU s已原样保留；新增单次air0事后诊断{r['cpu_seconds']:.6f} CPU s，23条合计{cost['successful_all_measured_tasks_CPU_seconds']:.6f} CPU s；",1)
 if '`experiments/fresnel_air0_control.py`' not in text:text+='\n追加事后诊断唯一入口是 `experiments/fresnel_air0_control.py`，固定唯一输出目录并拒绝覆盖；`experiments/fresnel_air0_finalize.py` 只汇总已付费端点、检查同锚点/配置及更新QA，不运行物理反演。原22条费用及QA封存为 `COST_AND_SELECTION_BASE22.json` / `FINAL_QA_BASE22.json`。\n'
 report.write_text(text)
 receipts=json.loads((OUT/'SUMMARY.json').read_text())+[json.loads((OUT/'FoamTwinDielTM_8GHz_n60_independent_real_nie/RECEIPT.json').read_text())]+json.loads((OUT/'MULTIFREQUENCY_SUMMARY.json').read_text())+[r]
 forbidden=[];npz=list(OUT.rglob('*.npz'))
 for p in npz:
  with np.load(p) as z:
   for k in z.files:
    if k in {'Y','Y_train','measured_Y','incident_measured','rawI','scattered_measured'}:forbidden.append(str(p.relative_to(ROOT))+':'+k)
 figs=['material_endpoints.png','neural_warm_starts.png','heldout_predictions.png','independent_nie.png','frequency_followup.png','air0_control.png'];missing=[f for f in figs if not (ROOT/'FIGURES/fresnel'/f).exists()]
 qa=dict(status='COMPLETED_ARTIFACT_QA',expected_completed_inversions=23,actual_completed_inversions=len(receipts),original_22_preserved=base,all_receipt_statuses_completed=all(v['status'].startswith('COMPLETED') for v in receipts),successful_CPU_seconds=sum(v['cpu_seconds'] for v in receipts),aggregate_CPU_matches_receipts=bool(np.isclose(sum(v['cpu_seconds'] for v in receipts),cost['successful_all_measured_tasks_CPU_seconds'],rtol=0,atol=1e-8)),air0_only_bound_change_checks=checks,forbidden_raw_measurement_keys=forbidden,missing_derived_figures=missing,public_npz_files=len(npz),image_QA='Originalfivefiguresviewed; air0figureawaitsvisualreview',GT_metrics='UNAVAILABLE; no measuredpixelNRMSE/IoU',NIE_scope=base['NIE_scope'],cost_scope='Main19+NIE1+MF2 preserved; plusonepostmeasurementair0DEVdiagnostic; knownfailureCPU>=405.52s+unclosedfirstlayoutCPU; sharedhost/noisolatedtimingclaim;2GHzAuncached')
 assert len(receipts)==23 and qa['all_receipt_statuses_completed'] and qa['aggregate_CPU_matches_receipts'] and not forbidden and not missing
 (OUT/'FINAL_QA.json').write_text(json.dumps(qa,indent=2)+'\n')
 (ROOT/'delegated/fresnel/AIR0_FINALIZATION_RECEIPT.json').write_text(json.dumps(dict(status='COMPLETED_METADATA_QA_AWAITING_FIGURE_REVIEW',completed_physical_runs_added=1,maximum_physics_processes=1,BLAS_threads=1,old_control_rerun=False,configuration_checks=checks,finalization_CPU_seconds=time.process_time()-start),indent=2)+'\n')
 print(json.dumps(dict(count=23,new_CPU=r['cpu_seconds'],total_CPU=cost['successful_all_measured_tasks_CPU_seconds'],metrics=r['metrics'],material=d,corrector=r['corrector'])))

if __name__=='__main__':main()
