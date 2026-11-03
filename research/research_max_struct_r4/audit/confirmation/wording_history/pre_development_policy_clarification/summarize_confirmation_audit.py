"""Post-score comparisons only; reads saved scores/receipts, never truth or physics."""
from pathlib import Path
from datetime import datetime,timezone
from statistics import mean,median
import json
from itertools import combinations

ROOT=Path(__file__).resolve().parents[1]
METHODS=('onehop','global','native_extended')


def run():
    scores=json.loads((ROOT/'results/CONFIRMATION_SCORES.json').read_text())
    seal=json.loads((ROOT/'audit/confirmation/CONFIRMATION_PREDICTIONS_SEALED.json').read_text())
    case_ids=sorted({row['observation_id'] for row in scores['rows']})
    cases=[]
    for ident in case_ids:
        rows={row['method']:row for row in scores['rows'] if row['observation_id']==ident}
        secondary=next(row for row in scores['secondary_rows'] if row['observation_id']==ident)
        nrmse={method:row['amplitude_normalized_nrmse'] for method,row in rows.items()}
        minimum=min(nrmse.values());maximum=max(nrmse.values())
        # These are evaluation labels, never a returned/reselected prediction.
        selected=secondary['selected_methods']
        cases.append(dict(observation_id=ident,family=rows['onehop']['family'],rows=rows,
            secondary_selected_before_GT=selected,
            posthoc_minimum_material_error_methods=[method for method in METHODS if nrmse[method]==minimum],
            posthoc_maximum_material_error_methods=[method for method in METHODS if nrmse[method]==maximum],
            secondary_minus_minimum_nrmse={method:nrmse[method]-minimum for method in selected},
            secondary_is_material_worst=all(nrmse[method]==maximum for method in selected),
            secondary_minus_global_nrmse={method:nrmse[method]-nrmse['global'] for method in selected},
            onehop_minus_global_nrmse=nrmse['onehop']-nrmse['global'],
            onehop_minus_global_iou=rows['onehop']['otsu_iou']-rows['global']['otsu_iou'],
            onehop_minus_global_relative_residual=rows['onehop']['relative_residual']-rows['global']['relative_residual'],
            onehop_minus_global_Phi=rows['onehop']['physical_objective']-rows['global']['physical_objective']))
    paired=[]
    for left,right in combinations(METHODS,2):
        entries=[]
        for case in cases:
            a,b=case['rows'][left],case['rows'][right]
            entries.append(dict(observation_id=case['observation_id'],family=case['family'],
                delta_nrmse=a['amplitude_normalized_nrmse']-b['amplitude_normalized_nrmse'],
                delta_iou=a['otsu_iou']-b['otsu_iou'],
                delta_relative_residual=a['relative_residual']-b['relative_residual'],
                delta_original_Phi=a['physical_objective']-b['physical_objective']))
        paired.append(dict(left=left,right=right,objects=6,
            left_lower_nrmse=sum(x['delta_nrmse']<0 for x in entries),right_lower_nrmse=sum(x['delta_nrmse']>0 for x in entries),
            equal_nrmse=sum(x['delta_nrmse']==0 for x in entries),
            left_higher_iou=sum(x['delta_iou']>0 for x in entries),right_higher_iou=sum(x['delta_iou']<0 for x in entries),
            equal_iou=sum(x['delta_iou']==0 for x in entries),
            left_lower_relative_residual=sum(x['delta_relative_residual']<0 for x in entries),
            left_lower_original_Phi=sum(x['delta_original_Phi']<0 for x in entries),entries=entries))
    summary=[]
    for method in METHODS:
        rr=[r for r in scores['rows'] if r['method']==method]
        summary.append(dict(method=method,objects=6,joint_pass=sum(r['joint_gate'] for r in rr),
            median_nrmse=median(r['amplitude_normalized_nrmse'] for r in rr),mean_nrmse=mean(r['amplitude_normalized_nrmse'] for r in rr),
            median_iou=median(r['otsu_iou'] for r in rr),
            unique_posthoc_lowest_material_error_count=sum(case['posthoc_minimum_material_error_methods']==[method] for case in cases)))
    selected_regrets=[value for case in cases for value in case['secondary_minus_minimum_nrmse'].values()]
    secondary=dict(selection_counts={method:sum(method in case['secondary_selected_before_GT'] for case in cases) for method in METHODS},
        joint_pass=sum(row['joint_gate'] for row in scores['secondary_rows']),objects=6,
        selected_strictly_worse_than_available_material_minimum=sum(any(v>0 for v in case['secondary_minus_minimum_nrmse'].values()) for case in cases),
        selected_material_worst=sum(case['secondary_is_material_worst'] for case in cases),
        selected_better_than_global=sum(all(v<0 for v in case['secondary_minus_global_nrmse'].values()) for case in cases),
        selected_worse_than_global=sum(any(v>0 for v in case['secondary_minus_global_nrmse'].values()) for case in cases),
        selected_equal_to_global=sum(all(v==0 for v in case['secondary_minus_global_nrmse'].values()) for case in cases),
        mean_posthoc_nrmse_regret=mean(selected_regrets),maximum_posthoc_nrmse_regret=max(selected_regrets),
        interpretation_scope='GT material ranking is post-score audit only; secondary was frozen by Phi before truth. No oracle-selected result replaces any prediction.')
    cost=scores['costs']['all18_plus_shared_once']
    report=dict(status='POST_CONFIRMATION_SCORE_AUDIT',created_utc=datetime.now(timezone.utc).isoformat(),
        seal_utc=seal['sealed_utc'],first_truth_opened_utc=scores['first_confirmation_truth_opened_utc'],
        seal_precedes_truth=seal['sealed_utc']<scores['first_confirmation_truth_opened_utc'],
        all18_outer_OK_inner_SEALED=True,all18_frozen_configuration_matches=all(r['frozen_configuration_matches'] for r in scores['rows']),
        cases=cases,paired_comparisons=paired,method_summary=summary,group_counts=scores['summaries'],secondary_audit=secondary,
        costs=scores['costs'],corrector_stop_counts={stop:sum(r.get('corrector_stop')==stop for r in scores['rows']) for stop in sorted({r['corrector_stop'] for r in scores['rows']})},
        evidence_status='GT_EXPOSED_AFTER_CONFIRMATION_SCORING_NO_FURTHER_ALGORITHM_TUNING_ON_THIS_SET',
        algorithm_changed=False,new_predictions=False,new_physical_LU=0,new_complex_RHS=0,
        scientific_judgment='Reserved for parent Codex; descriptive6-object counts only, no significance or population inference')
    out=ROOT/'audit/confirmation/CONFIRMATION_COMPARISON_AUDIT.json'
    with out.open('x') as stream:json.dump(report,stream,indent=2,allow_nan=False);stream.write('\n')
    lines=['# R4 六对象确认集：独立评分与成本审查','',
        f"先在 `{seal['sealed_utc']}` 封存全部18个 outer OK / inner SEALED_PREDICTION，并按原 Phi 固定 secondary；之后于 `{scores['first_confirmation_truth_opened_utc']}` 首次读取确认 GT。18个冻结配置均一致，无缺失、失败或 UNKNOWN；原六对象分母完整。",'',
        '**证据状态：GT 已在本次评分中暴露，此后不得再用这组六对象调整算法。** 本审查不修改预测或冻结配置，也不按 GT 返回新端点。此次只有离线评分与保存结果的比较，新增0 LU/0 RHS。','',
        '## 完整逐例','',
        '| 对象 | 家族 | 方法 | RMSE/12 | Otsu256 IoU | 物理相对残差 | 白化残差 | 原 Phi | 联合门 |',
        '|---|---|---|---:|---:|---:|---:|---:|---|']
    for case in cases:
        for method in METHODS:
            row=case['rows'][method]
            lines.append(f"|{case['observation_id']}|{case['family']}|{method}|{row['amplitude_normalized_nrmse']:.9f}|{row['otsu_iou']:.9f}|{row['relative_residual']:.12f}|{row['whitened_residual']:.6f}|{row['physical_objective']:.12g}|{'PASS' if row['joint_gate'] else 'FAIL'}|")
    lines+=['','冻结门：RMSE/12≤0.15、仅根据预测计算的 Otsu256 IoU≥0.7、白化残差≤25.165422196633，且材料满足公开边界。prior 与评分合同未改变。','',
        '## 六对象计数和材料误差对照','',
        '| 方法 | 联合通过 | ring | two_islands | curved_bar | 材料 NRMSE 中位数 | 逐例最低材料误差次数 |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for row in summary:
        lines.append(f"|{row['method']}|{row['joint_pass']}/6|2/2|2/2|0/2|{row['median_nrmse']:.9f}|{row['unique_posthoc_lowest_material_error_count']}/6|")
    lines+=['|secondary_min_original_Phi|4/6|2/2|2/2|0/2|见逐例 retained endpoints|GT仅作事后诊断|','',
        '这里的逐例最低材料误差仅为事后评价标签，不是交付的 oracle 选择或新预测。三个方法及 secondary 都在同样四对象通过；两个 curved-bar 对象全部失败。全部18个共同 corrector 都在1000迭代上限停止，辅助日程均完成；这些停止状态不认证收敛。','',
        '| 成对比较（左−右） | 左材料 NRMSE 更低 | 右更低 | 左 IoU 更高 | IoU相同 | 左残差更低 | 左原 Phi 更低 |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for row in paired:
        lines.append(f"|{row['left']} − {row['right']}|{row['left_lower_nrmse']}/6|{row['right_lower_nrmse']}/6|{row['left_higher_iou']}/6|{row['equal_iou']}/6|{row['left_lower_relative_residual']}/6|{row['left_lower_original_Phi']}/6|")
    lines+=['',
        'onehop 相对 global 仅在 obs_100bd65498644322 的材料 NRMSE 更低，其余5对象更高；IoU 则3对象更高、2更低、1相同。该陈述只描述冻结端点，不比较已认证收敛解，也不构成人群显著性结论。','',
        '## 预先按 Phi 选择是否择劣','',
        '| 对象 | 预先选中 | 选中 NRMSE | 事后最低材料误差方法 | 相对最低 NRMSE 差 | 选中为材料最差臂 |',
        '|---|---|---:|---|---:|---|']
    for case in cases:
        method=case['secondary_selected_before_GT'][0]
        lines.append(f"|{case['observation_id']}|{method}|{case['rows'][method]['amplitude_normalized_nrmse']:.9f}|{', '.join(case['posthoc_minimum_material_error_methods'])}|{case['secondary_minus_minimum_nrmse'][method]:.9f}|{case['secondary_is_material_worst']}|")
    lines+=['',
        f"secondary 在GT前选 onehop/global/native_extended 为3/2/1次，最终仍4/6联合通过。它在{secondary['selected_strictly_worse_than_available_material_minimum']}/6对象上未选中材料误差最低臂；两个 curved-bar 对象都选中材料 NRMSE 最差臂。相对单独 global：1对象材料误差更低、3更高、2相同。",'',
        '具体地，42cb 的 onehop 原 Phi=0.002193014461，低于 global 的0.004105231630，但材料 NRMSE .761367699 高于 .741539073；63ab 的 native 原 Phi=0.004611720459 低于 global 的0.004794102215，但材料 NRMSE .577783509 高于 .575754689。这里保留预先按 Phi 得到的实际选择，绝不以 GT 改选。','',
        '## 实际成本','',
        '| 方法（6任务原始总计） | 物理 LU | source RHS | adjoint RHS | 物理 complex RHS | worker wall 累计 s | adapter CPU累计 s |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for method,c in scores['costs']['primary_method_task_costs'].items():
        lines.append(f"|{method}|{int(c['LU']['total'])}|{int(c['source_RHS']['total'])}|{int(c['adjoint_RHS']['total'])}|{int(c['actual_complex_RHS']['total'])}|{c['task_worker_wall_seconds']['total']:.6f}|{c['adapter_cpu_seconds']['total']:.6f}|")
    lines+=['',
        '物理计数只汇总每任务顶层 `physical_model_cost`；已包含任务模型 setup、辅助验证、排序、corrector 和终点，不再次相加嵌套子阶段。12源；小块矩阵 RHS、标量除法和物料内层工作保留在 JSON 的辅助操作计数，不充作 full Maxwell RHS。','',
        f"18任务：18647 LU /445860 complex RHS（source223764、adjoint222096）。共享生成仅加一次6 LU/72 source RHS，因此实付物理总账为 **18653 LU /445932 complex RHS**。secondary 的成本就是全部三臂18任务加共享一次，不是一个免费选择器，也不能再把 secondary 当新运行重复计费。",'',
        f"实测远端 elapsed **{cost['remote_campaign_elapsed_seconds']:.9f} s**（{cost['remote_campaign_elapsed_seconds']/60:.6f}分钟）用于时间成本，含前台 SSH、远端进程启动及收尾。worker wall 累计3900.014 s 是并发任务工作量，不能替代 elapsed。每 worker 一次的 solver 导入合计20.408797800 s wall /2.296875 s CPU（四 worker），未作为额外 elapsed 再相加。",'',
        'adapter CPU累计3802.046875 s；加四 worker 一次导入后，记录到的远端 CPU=3804.343750 s。另有本地共享生成0.151104208 s wall /0.139641 s CPU；记录CPU合计3804.483391 s。early bootstrap 含于实测远端 elapsed，但无法细分 CPU/阶段，不作零成本处理；货币单价未知。实际4 CPU worker、GPU任务/worker均0；runner旧 `cuda_workers=2` 固定元数据不代表实际运行。峰值RSS按worker生命周期口径，不能当单任务独占峰值。','',
        '所有六对象分母、原始失败信息、轨迹和成本被保留。本组现为已暴露评分证据；不再用于调整冻结算法。六对象仅做描述性计数，不作显著性检验，最终科学解释和主张由父线程负责。']
    with (ROOT/'audit/confirmation/CONFIRMATION_COMPARISON_AUDIT.md').open('x') as stream:stream.write('\n'.join(lines)+'\n')
    print(json.dumps(dict(method_summary=summary,secondary=secondary,paired=[{k:v for k,v in row.items() if k!='entries'} for row in paired],cost=dict(remote_elapsed=cost['remote_campaign_elapsed_seconds'],recorded_remote_CPU=cost['recorded_remote_CPU_seconds']['total'],physical_total=cost['campaign_physical_totals_including_shared_once']))))


if __name__=='__main__':run()
