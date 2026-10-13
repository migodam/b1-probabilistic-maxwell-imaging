"""Finalize machine-readable/document delivery from paid outputs only."""
from pathlib import Path
import json,csv,time,platform
import numpy as np,scipy
from src.b1.budget import Ledger
from src.b1.pilot import dump
root=Path.cwd();c=json.loads(Path('configs/protocol.json').read_text());first=json.loads(Path('results/PRIMARY_SUMMARY.json').read_text());ref=json.loads(Path('results/refinement/SUMMARY.json').read_text());cost=Ledger('results',384).summary();audit=json.loads(Path('results/FORWARD_AUDIT.json').read_text());cal=json.loads(Path('results/CALIBRATION_RESULTS.json').read_text());native=json.loads(Path('research/delegated/probability_review/CALIBRATION_INTERVAL_REANALYSIS_NATIVE.json').read_text());compare=json.loads(Path('research/delegated/probability_review/DERIVED_METHOD_COMPARISON_NATIVE.json').read_text())
rows=[]
for name in ['distinguishable','competitive','information_limited','gain_mismatch','gain_marginalized','representation_mismatch']:
 r=first['cases'][name];rr=ref['cases'][name];p=r['result']['topology_probabilities'];pp=rr['result']['topology_probabilities']
 rows.append(dict(case=name,p_solid_original=p[0],p_single_void_original=p[1],p_double_void_original=p[2],p_solid_refined=pp[0],p_single_void_refined=pp[1],p_double_void_refined=pp[2],original_reference=r['refinement_5_vs_7']['status'],refined_reference=rr['alternative_reference']['status'],predictive_tail_original=r['result']['posterior_predictive_tail_probability'],refinement_P_error=rr['alternative_reference']['integration_17_vs_25']['max_probability_difference'],refinement_logBF_error=rr['alternative_reference']['integration_17_vs_25']['max_log_BF_difference']))
with Path('results/CASE_PROBABILITIES.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
comparison=[]
for name,v in compare['summaries'].items():
 m=v['metrics'];row=dict(method=name,cold_nominal_full_states=v['cold_nominal_material_frequency_states'],actual_new_for_reanalysis=0,Brier=m['Brier_multiclass']['estimate'],log_score=m['log_score']['estimate'],max_prob_error=v.get('empirical_max_probability_error'),mean_prob_error=m['mean_max_probability_error_vs_gauss7']['estimate'],both_reference_gate_fraction=v.get('reference_gate_vs_gauss7',{}).get('fraction_within_both_reference_tolerances'))
 comparison.append(row)
with Path('results/METHOD_COMPARISON.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(comparison[0]));w.writeheader();w.writerows(comparison)
alljobs=[json.loads(x) for x in Path('results/JOB_LEDGER.jsonl').read_text().splitlines()]
visual=json.loads(Path('research/delegated/visualization/receipt.json').read_text());oldvisual=json.loads(Path('research/delegated/visualization/receipt_before_visual_review.json').read_text());review=json.loads(Path('research/delegated/probability_review/REVIEW_NATIVE.json').read_text())
cpu=dict(host_wrapper_jobs=len(alljobs),host_wrapper_wall_seconds=sum(r['wall_seconds'] for r in alljobs),host_wrapper_cpu_seconds=sum(r['cpu_seconds'] for r in alljobs),
         native_visual_stages=visual['visualization_stage_measured_totals'],earlier_renderer_internal_wall_seconds=oldvisual['wall_seconds'],earlier_renderer_internal_cpu_seconds=oldvisual['cpu_seconds'],earlier_renderer_scope='internal timer; interpreter startup unknown; preserved separately from native process totals',
         native_review_artifact='research/delegated/probability_review/REVIEW_NATIVE.json',native_derived_comparison_wall_seconds=compare['wall_seconds'],native_derived_comparison_cpu_seconds=compare['process_cpu_seconds'],
         no_double_sum='physics wall is nested in host jobs; review/visual totals have their own scope; do not sum nested internal receipt plus process wrapper',
         excluded_from_measured_jobs='interactive reading/writing, tool waits/visual review, scholarly network calls, model-service reasoning; these are not claimed as zero or fully measured research time',
         external_model_services='early named-profile workers were stopped at user request; costs unmeasured; all later agents native Codex',gpu_seconds=0,new_nn_training_updates=0)
dump('results/COMPUTE_COST_SUMMARY.json',cpu)
with Path('results/FULL_CALL_COSTS.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['kind','new_material_frequency_evaluations']);w.writerows(cost['by_kind'].items());w.writerow(['TOTAL',cost['new_physics_evaluations']])
status=dict(overall='GO_LIMITED_TOPK_PROTOTYPE',scope='conditional on explicit three-template continuous-material prior, concentric 3D vector-Mie only',
 gates={
  'forward_interface':dict(status='PASS',evidence='results/FORWARD_AUDIT.json; 19 checks from16 calls; independent28 checks included in384'),
  'probability_target_and_topk':dict(status='PASS',evidence='proper topology/material prior, continuous quadrature, determinant, L/q, third-event mass; native review'),
  'original_reference_competitive_information_limited':dict(status='PASS',evidence='5/7 deltaP <=.02 and deltaLogBF<=.1'),
  'original_reference_strong_and_mismatch':dict(status='FAIL',evidence='4 original REFERENCE_UNRESOLVED cases retained'),
  'alternative_interpolated_reference':dict(status='EMPIRICAL_PASS',certified=False,evidence='results/refinement/SUMMARY.json; no uniform prior-support remainder'),
  'gauss5_development_tolerance':dict(status='PASS',data='24parents,720noise; exposed development; referenceGauss7; nominalcold75',both_gates_passed=720),
  'fixed_importance_development_tolerance':dict(status='FAIL',data='proposal fixed from competitivey0; samepanel, cold27+48=75',both_gates_passed=124),
  'model_calibration':dict(status='HOLD',evidence='exploratory parent-cluster scores/coverage; some bin deviations; no rankSBC or external reality calibration'),
  'conditional_material_intervals':dict(status='HOLD',evidence='quantile bin approximation; epsilon real17/25 change .01515 in a strong case; coverage cluster uncertainty'),
  'general_geometry_support':dict(status='NOT_RUN',negative_control='registered support misses .68-radius cavity; predictive fail; outside support broader mass undefined'),
  'native_A17_feedback':dict(status='NOT_RUN',no_acceleration_claim=True,evidence='generic FD/SVD affine audit only; rank2 in2 material dims'),
  'strong_TAP_novelty':dict(status='NOT_ESTABLISHED',evidence='standard Bayes/quadrature/Mie/importance/interpolation; limited literature ledger, no new Maxwell remainder theorem'),
  'physics_budget':dict(status='PASS',charged=384,cap=384,remaining=0,physical_failures=0,history_2d=74,combined_analytic_calls=458),
  'implementation_tests':dict(status='PASS',unit_count=24,evidence='results/jobs/final-tests-003/stderr.log; both native-reviewed bugs corrected'),
  'static_figures':dict(status='REVIEWED',count=8,evidence='native visualization audit; no human/scientific acceptance asserted'),
  'live_html_browser':dict(status='NOT_RUN',reason='CUA local fileURL security rejection; no workaround'),
  'public_release':dict(status='NOT_RUN',reason='local delivery, no automatic publication requested')},
 next_phase='new budget/permission before more physics; first add one geometric freedom or a proper other branch, then re-establish reference. No largeNN.',
 stop=['new physics beyond384','claim current sampler meetsreference','native A17 acceleration/lowrank compression in2Dparameters','beam/residual/sample frequency as posterior mass','largeNN without a demonstrated bottleneck','global/realworld calibrated probability or TAPnovelty claims'])
dump('GATE_DECISION.json',status)
source=dict(name='B1',date='2026-10-13',input_zip="/Volumes/migodam's-external-brain/archive/user-home/downloads/A32_Research_Package.zip",zip_crc_integrity='PASS,10members',new_sha_checks=False,
 original_package='inputs/source_package/A32_Research_Package',source_package_status='THEORY_AND_SMALL_CHECKS_ONLY;74historical2Dcalls retained',
 reused=[dict(path="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/A13/REVISION_R2/src/mie_validation.py",use='unmodified pure function extracts vendor/a13_mie_functions.py',license='local original research; no new license assumed'),dict(path="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/A15/research/delegated/layered_reference/run_reference_checks.py",use='existing plain helicity treams coefficient map tools/treams_one.py in isolated ABI'),dict(path="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/A6/vendor/treams",use='4independent layered references, MIT dependency, not rebundled'),dict(path="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/A17/CLOSURE_R1/manuscript/CORE_THEORY.tex",use='read-only exact Galerie Schur/material theory check; no native solver modifications'),dict(path="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/A32/a32-material-risk-directed-feedback/A32_E0_E1_DECISION.md",use='distinct frozen E0-FAIL source, not the provided probability theorypackage'),dict(path="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/A30/a30-maxwell-material-manifold",use='negative/capacity evidence and budget idioms; no old arrays read')],
 environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),cpu='Apple Silicon,BLAS1',optional_reference='treams0.4.7 +vendorednumpy2.2.6 scipy1.15.3 in isolated subprocess'),
 artifacts='code, reports, configs, append-only cost/job ledgers, generated observations/posteriors, plots, dashboard, corrected/failed snapshots, native reviews and literature records',published=False,
 upstream_modification=False,locked_arrays_read=0,old_truth_arrays_read=0)
dump('SOURCE_MANIFEST.json',source)
# Generate compact results sheet, maintaining original and alternative labels.
lines=['# B1 实验结果','', '2026-10-13。原始、派生与经验近似结果分列，不是严格概率认证。','', '|案例|实心%|单空腔%|双空腔%|原生5/7参考|经验细化|','|---|---:|---:|---:|---|---|']
for r in rows:lines.append(f"|{r['case']}|{100*r['p_solid_refined']:.6g}|{100*r['p_single_void_refined']:.6g}|{100*r['p_double_void_refined']:.6g}|{r['original_reference']}|{r['refined_reference']}|")
lines+=['','表中概率来自第二轮，百分比条件于三个模板世界；原生P数据逐项在CASE_PROBABILITIES.csv。Top2不重归一化，第三event剩余保留；模板之外的更大物理世界未量化。','', '## 同数据/名义冷成本对照','', '|方法|名义冷全波状态|Brier|log score|平均最大P差|完整P/logBF门槛比例|','|---|---:|---:|---:|---:|---:|']
for r in comparison:
 def fmt(v):return '∞ / 非概率点决策' if v is None else f'{v:.6g}'
 lines.append(f"|{r['method']}|{r['cold_nominal_full_states']}|{fmt(r['Brier'])}|{fmt(r['log_score'])}|{fmt(r['mean_prob_error'])}|{fmt(r['both_reference_gate_fraction'])}|")
lines+=['','参考Gauss7，在24parent×30noise的已暴露开发panel上；固定importance的proposal来自一个y0，没有逐数据适应或重抽样。对照setup名义状态从已支付结果重算，实际本次新增全波0；Gauss5冷75，最初pipeline新72+复用3中心。不可把该表写成一般sampler优势或冷wall实测。','', '## 校准与材料不确定性','', f"Brier={cal['Brier_multiclass']:.6g}，log={cal['log_score']:.6g}，Top2真H覆盖={cal['top2_empirical_coverage']:.6g}，平均报告mass={cal['top2_mean_reported_mass']:.6g}。", '', '独立review的parent分层bootstrap、边际/真实H条件coverage与分箱偏差见CALIBRATION_INTERVAL_REANALYSIS_NATIVE.json。部分分箱区间不含0，不作校准PASS；rankSBC未运行。quantile仍是prior-mass-bin近似。','', '## 物理与成本','', '|任务|新材料–频率状态|','|---|---:|']
for k,n in cost['by_kind'].items():lines.append(f'|{k}|{n}|')
lines+=[f"|合计|{cost['new_physics_evaluations']}/384|",'',f"真实forward包围计时之和{cost['full_forward_wall_seconds']:.6f}s，包含treams子进程；它嵌套于job wall，不可重复相加。GPU0，NN0，旧locked/truth0。历史原包2D74+本阶段3D384=458；A31既有微型物理另列。",'',f"当前10个wrapper jobs计时sum wall={cpu['host_wrapper_wall_seconds']:.6f}s，sumCPU={cpu['host_wrapper_cpu_seconds']:.6f}s（报告生成等后续记录按JOB_LEDGER为准）。native视觉阶段wall19.845079s/CPU18.896490s；早期renderer内部wall5.3861s/CPU5.2974s另保留；native数值review3.152s/CPU2.802s按自身receipt。网络、模型推理、交互读写与工具等待未完整测量，不称免费或全部研究wall。",'', '## 失败/未运行','', '原生四案例reference不足、单例importance误差.02565、冻结panelimportance两容差仅124/720、gain/模板错配PPC异常都保留。两个review问题（coverage标签、gain NPZ元数据）已修复；第一次修正误用系统解释器的失败与错误CLI重放亦保留。不存在新物理失败，但实现失败不是零。','', '原生A17反馈、未知位置/边界、复杂一般形状、uniform概率认证、rankSBC、锁定测试、NN、公开发布与HTML现场浏览器验收均未运行。']
Path('RESULTS.md').write_text('\n'.join(lines)+'\n')
Path('DECISION.md').write_text('''# B1 最终研究裁决

**GO_LIMITED_TOPK_PROTOTYPE**：已实现从散射数据到3D vector-Mie模板材料事件概率、Top-K代表图、第三分支余量与不确定性。竞争主例50.64%单空腔、40.43%双空腔、8.93%实心。它是有明确proper prior的条件材料posterior，不是多图美化。

**HOLD_GENERAL_GEOMETRY_AND_CERTIFICATION**：几何固定、模板支持有限，位置没有被推断。原生全网格4个reference失败；第二轮经验插值核验通过，但无全域余项证书。24parent的校准只是开发检查，部分分箱存在偏差、CDF细化未认证，无posterior-rank SBC或仪器校准。

**选择积分，暂停当前importance扩张**：相同720数据上，名义cold75的Gauss5通过720/720完整参考门槛；来自单竞争y0的冻结importance27+48只通过124/720。这个结论限定在冻结proposal/已暴露低维panel，不是一般高维方法的NO-GO。

**NO_GO_STRONG_TAP_CLAIM_AT_THIS_STAGE**：没有新Maxwell-specific认证余项/可变量几何理论或同概率精度部署加速。A17只作generic affine签名诊断，没有native交换，rank2/材料2维不压缩。现有组合来自已知Mie/Bayes/importance/插值工具。大型NN、beam confidence当posterior、无密度修正的频率、以residual赢者替代evidence均停止。

**预算收口**：新物理384/384，完整历史不删除，原包二维74另外保留，合计458解析材料–频率评估。GPU/NN/locked0。后续新物理或扩展support需要用户另行授权；本阶段代码/结果/报告/失败/成本完成本地交付。下一阶段优先增加一个可积分的几何变量或proper其他分支，重建reference，再考虑likelihood-informed/tempered/RJ方案。

详细PASS/FAIL/HOLD/NOT_RUN与证据见[GATE_DECISION.json](GATE_DECISION.json)、[中文报告](RESEARCH_REPORT_ZH.md)、[独立native review](research/delegated/probability_review/REVIEW_NATIVE.md)。公开发布未执行。HTML现场渲染被本地file访问策略拒绝，记NOT_RUN；8张静态科学图已视觉核对。
''')
print(json.dumps(dict(status='documents_and_decision_written',new_physics=0,cases=len(rows),methods=len(comparison))))
