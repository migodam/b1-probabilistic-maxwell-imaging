"""OFFLINE_DIAGNOSTIC from saved receipts/history and exposed truth/noise only.

No forward, gradient, optimizer, or model module is imported or invoked.
"""
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter
import json
import numpy as np
from offline_score_common import public_metrics,prior

ROOT=Path(__file__).resolve().parents[1]
CASES=('obs_42cb233e4f224bb3','obs_63abafa750fc483c')
METHODS=('onehop','global','native_extended')


def read(path):return json.loads(path.read_text())


def tail_drop(history,count):
    running=np.minimum.accumulate(history[:,0])
    offset=min(count,len(running)-1)
    before=float(running[-1-offset]);after=float(running[-1])
    return dict(evaluation_intervals=offset,start_evaluation_one_based=len(running)-offset,
        end_evaluation_one_based=len(running),best_original_Phi_before=before,
        best_original_Phi_after=after,absolute_decrease=before-after,
        relative_decrease=(before-after)/abs(before) if before else None,
        scope='Running best over objective evaluations, including rejected line-search trials; not optimizer iterations')


def run():
    remote=ROOT/'results/confirmation_remote'
    rows=[];gt_rows=[]
    for ident in CASES:
        public_path=ROOT/'data_public/confirmation'/(ident+'.npz')
        truth_path=ROOT/'sealed_truth/confirmation'/(ident+'.npz')
        public=public_metrics(public_path)
        # Already exposed truth; parent specifically authorizes saved clean/noise use.
        with np.load(truth_path,allow_pickle=False) as data:
            truth=data['m'].copy()
            clean=data['Y_clean'].copy() if 'Y_clean' in data else None
            noise=data['noise'].copy() if 'noise' in data else (public['Y']-clean if clean is not None else None)
        white=float(np.linalg.norm(noise)/float(public['sigma'])) if noise is not None else None
        relative=float(np.linalg.norm(noise)/np.linalg.norm(public['Y'])) if noise is not None else None
        reg=prior(truth,public,(24,24))
        lo,hi=map(float,public['material_bounds'])
        violation=float(max(0.,lo-truth.min(),truth.max()-hi))
        gt_rows.append(dict(observation_id=ident,classification='OFFLINE_DIAGNOSTIC_EXPOSED_GT',
            material_min=float(truth.min()),material_max=float(truth.max()),material_legal=violation<=1e-12,
            material_constraint_violation=violation,GT_whitened_residual=white,
            noise_radius99=float(public['noise_radius']),GT_within_frozen_noise_radius=bool(white<=float(public['noise_radius'])) if white is not None else None,
            GT_relative_residual=relative,GT_data_Phi=.5*relative**2 if relative is not None else None,
            GT_original_prior=reg,GT_original_Phi=.5*relative**2+reg if relative is not None else None,
            source=str(truth_path.relative_to(ROOT)),noise_source='stored noise' if noise is not None and clean is not None else 'stored noise or publicY-Y_clean',
            no_forward=True,missing_reason=None if noise is not None else 'Neither stored noise nor cleanY available; no new solve permitted'))
        for method in METHODS:
            receipt_path=remote/'results/confirmation'/f'{ident}__{method}.json'
            outer=read(receipt_path);rec=outer['result']
            trajectory_path=remote/'predictions'/ident/(method+'.npz')
            with np.load(trajectory_path,allow_pickle=False) as data:
                history=data['corrector_history'].copy()
                candidate_phi=data['candidate_original_Phi'].copy()
                stage_ids=data['history_stage'].copy()
                inner_messages=data['history_material_inner_message'].copy() if 'history_material_inner_message' in data else None
                iterations=data['history_global_iteration'].copy()
            final_phi=float(rec['physical_objective'])
            matches=np.flatnonzero(history[:,0]==final_phi)
            endpoint_grad=float(history[matches[-1],1]) if len(matches) else None
            selected_index=int(np.argmin(candidate_phi))
            start_phi=float(candidate_phi[selected_index])
            last_stage=len(rec['configuration']['stages'])-1
            stage_mask=(stage_ids==last_stage)&np.isfinite(iterations)
            inner_stops=dict(Counter(str(x) for x in inner_messages[stage_mask])) if inner_messages is not None else None
            rows.append(dict(observation_id=ident,method=method,classification='OFFLINE_DIAGNOSTIC',
                outer_status=outer['status'],inner_status=rec['status'],
                receipt=str(receipt_path.relative_to(ROOT)),trajectory=str(trajectory_path.relative_to(ROOT)),
                corrector_stop=rec['corrector']['stop'],corrector_iterations=rec['corrector']['iterations'],
                corrector_evaluations=rec['corrector']['evaluations'],
                selected_auxiliary_stage=selected_index-1,corrector_seed_original_Phi=start_phi,
                returned_original_Phi=final_phi,corrector_original_Phi_decrease=start_phi-final_phi,
                last_evaluated_original_Phi=float(history[-1,0]),last_evaluated_raw_gradient_inf=float(history[-1,1]),
                returned_endpoint_raw_gradient_inf=endpoint_grad,
                returned_endpoint_raw_gradient_source='stored corrector evaluation with exactly matching Phi; last exact match' if len(matches) else None,
                returned_endpoint_projected_gradient_inf=rec.get('projected_gradient_inf'),
                projected_gradient_unavailable_reason='Receipt and trajectory retain only raw gradient infinity norms; no gradient vector/signs to reconstruct bound-projected KKT gradient',
                raw_gradient_semantics='Unscaled original objective gradient infinity norm; not projected gradient or stationarity certificate',
                original_PGTOL=1e-8,scaled_gtol=.01,objective_scale=1e6,
                tail100=tail_drop(history,100),tail10=tail_drop(history,10),
                final_auxiliary_stage=last_stage,final_auxiliary_stage_mf=rec['configuration']['stages'][-1]['mf'],
                final_auxiliary_stage_updates=int(np.count_nonzero(stage_mask)),
                final_auxiliary_stage_status=rec['stage_history'][-1]['status'],
                prior_stage_physical_Phi=float(candidate_phi[-2]),final_stage_physical_Phi=float(candidate_phi[-1]),
                final_auxiliary_stage_physical_Phi_decrease=float(candidate_phi[-2]-candidate_phi[-1]),
                final_auxiliary_stage_inner_stop_counts=inner_stops,
                final_physical_relative_residual=rec['physical_relative_residual'],
                final_whitened_residual=rec['physical_whitened_residual'],
                noise_radius99=rec['noise_radius99'],
                budget_insufficiency='NOT_EXCLUDED',converged_basin_certificate=False))
    result=dict(status='OFFLINE_DIAGNOSTIC',created_utc=datetime.now(timezone.utc).isoformat(),
        scope='Only2 failed exposed confirmation objects; original18 predictions/configs/scores unchanged',
        gradient_column_source='pre-generation executed_source_snapshot/.../cooperative_paths.py::refine stores [original Phi,||raw g||inf] per objective evaluation',
        rows=rows,GT_saved_data_diagnostics=gt_rows,
        budget_conclusion='NOT_EXCLUDED: every corrector and final auxiliary stage hit iteration cap; no projected-gradient certificate is stored. Tail descent varies and cannot predict eventual recovery.',
        identifiability_scope='A legal stored GT whose saved clean prediction/noise fits the public data is not a uniqueness proof. Failed candidate residuals do not prove information nonidentifiability or exclude optimization budget limitations.',
        development_policy='Frozen confirmation algorithm and18 results remain immutable and cannot be reclaimed as independent validation. Mechanism-driven revisions may be separately declared POST_CONFIRMATION_DEVELOPMENT with future independent validation.',
        parent_declared_case_selection='42cb selected solely by public failure and smallest opaque ID; no GT-ranked choice',
        no_forward=True,no_gradient_recomputation=True,new_physical_LU=0,new_complex_RHS=0,
        new_iterations=0,new_tests=False,original_scores_changed=False,original_predictions_changed=False,
        scientific_judgment='Reserved for parent Codex')
    out=ROOT/'audit/confirmation/FAILURE_STOP_DIAGNOSTIC.json'
    with out.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    lines=['# 两例失败的只读停止诊断','',
        '**OFFLINE_DIAGNOSTIC**：仅读已暴露确认集两例失败对象的6份原始 receipt、保存轨迹及 parent 授权的 GT cleanY/noise/material。没有新 forward、梯度、迭代或测试；原18预测、原onehop及原评分数字均不变。','',
        '## 研究循环边界澄清','',
        '不得用已暴露数据调整已冻结确认算法或再声称独立确认；新机制修正可以另立 **POST_CONFIRMATION_DEVELOPMENT** 版本，并须接受新的独立验证。这一边界不阻止用户授权的失败分析和研究循环。父线程已声明新 partition-average 版本独立于冻结算法，42cb 仅按 public不合格＋最小opaqueID选择；本任务未运行或评分新版本。','',
        '## 停止、梯度与 corrector 下降','',
        '冻结执行源码确认 `corrector_history[:,1]` 是未缩放原目标的 **raw gradient infinity norm**，不是投影梯度。六份 receipt/NPZ 均未保存投影梯度或完整梯度向量，故 projected-gradient/KKT 资格记为 null；材料存在活动边界，不能把 raw norm 直接与 PGTOL 比较后宣布非平稳。原目标乘1e6、缩放后gtol=.01，等效原投影梯度目标1e-8。','',
        '| 对象 | 方法 | corrector停止 | seed Phi → 返回 Phi | 末次raw梯度∞ | 末100评估best-Phi下降 | 末10评估best-Phi下降 |',
        '|---|---|---|---:|---:|---:|---:|']
    for row in rows:
        lines.append(f"|{row['observation_id']}|{row['method']}|MAXITER {row['corrector_iterations']}|{row['corrector_seed_original_Phi']:.12g} → {row['returned_original_Phi']:.12g}|{row['last_evaluated_raw_gradient_inf']:.9g}|{row['tail100']['absolute_decrease']:.9g}|{row['tail10']['absolute_decrease']:.9g}|")
    lines+=['',
        '末100/10为目标函数评估区间，不是优化迭代；包括 line-search rejected trials，因此采用已见最小 Phi 的 running best。端点梯度只来自保存历史中与返回 Phi 精确相同的记录，没有重算。下降接近零只能描述此有限窗口，不能证明全局/材料收敛；持续下降也不保证增预算后会恢复GT。','',
        '## 最后辅助阶段的原物理 Phi','',
        '| 对象 | 方法 | 最后mf/更新 | stage停止 | 前阶段 Phi → 最后阶段 Phi | 有符号下降 | corrector选中stage |',
        '|---|---|---|---|---:|---:|---:|']
    for row in rows:
        lines.append(f"|{row['observation_id']}|{row['method']}|mf{row['final_auxiliary_stage_mf']}/{row['final_auxiliary_stage_updates']}|{row['final_auxiliary_stage_status']}|{row['prior_stage_physical_Phi']:.12g} → {row['final_stage_physical_Phi']:.12g}|{row['final_auxiliary_stage_physical_Phi_decrease']:.9g}|{row['selected_auxiliary_stage']}|")
    lines+=['',
        '辅助阶段 Phi 仅取原已执行的物理 candidate 验证；跨阶段同时改变容量/状态设置，不能当固定子问题的逐迭代收敛曲线。42cb onehop 的mf10物理 Phi 比mf6更高，原物理选择实际保留mf6；其余路径的最后stage下降保存在表中。内层停止计数见JSON。','',
        '## 合法GT自身的已存数据证据','',
        '| 对象 | GT材料合法 | GT白化残差 / 原99%半径 | GT原Phi（含同prior） | failed三臂白化残差范围 |',
        '|---|---|---:|---:|---:|']
    for gt in gt_rows:
        failed=[r['final_whitened_residual'] for r in rows if r['observation_id']==gt['observation_id']]
        white='null' if gt['GT_whitened_residual'] is None else f"{gt['GT_whitened_residual']:.9f} / {gt['noise_radius99']:.9f}"
        phi='null' if gt['GT_original_Phi'] is None else f"{gt['GT_original_Phi']:.12g}"
        lines.append(f"|{gt['observation_id']}|{gt['material_legal']}|{white}|{phi}|{min(failed):.6f}–{max(failed):.6f}|")
    lines+=['',
        'GT指标只用已存noise（若缺失则publicY−Y_clean）和GT材料的原prior，完全没有求解。已有合法GT及其保存的clean观测可解释数据这一事实，与当前failed候选远超noise同时成立。它不证明唯一性，也不能从失败候选推得信息不可辨识。','',
        '**预算不足仍未排除（NOT_EXCLUDED）**：六个 corrector 和最后辅助阶段均由迭代上限停止，投影梯度资格未保存，部分末窗口仍下降；这些记录不足以认证收敛盆地或排除更长预算/不同路径的影响。没有据此自动增加预算或改变冻结协议。新修正与原确认结果必须分版本、分分母，最终机制解释和科学判断留父线程。']
    with (ROOT/'audit/confirmation/FAILURE_STOP_DIAGNOSTIC.md').open('x') as stream:stream.write('\n'.join(lines)+'\n')
    print(json.dumps(dict(status=result['status'],GT=gt_rows,stops=[{key:r[key] for key in ('observation_id','method','last_evaluated_raw_gradient_inf','tail100','final_auxiliary_stage_physical_Phi_decrease')} for r in rows],budget_insufficiency='NOT_EXCLUDED',new_physical_LU=0)))


if __name__=='__main__':run()
