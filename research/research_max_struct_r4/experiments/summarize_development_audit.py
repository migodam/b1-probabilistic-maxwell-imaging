"""Receipt-only development score/cost audit; never import a physics module or read GT."""
from pathlib import Path
from datetime import datetime, timezone
from collections import Counter
import json
import time

ROOT = Path(__file__).resolve().parents[1]
NAMES = ('block_absorption', 'self_absorption', 'block_absorption_cached',
         'native_extended', 'onehop_absorption', 'global_common_corrector')


def read(path):
    return json.loads(path.read_text())


def summarize():
    start, cpu = time.perf_counter(), time.process_time()
    scores = read(ROOT/'results/DEVELOPMENT_SCORES.json')
    old = read(ROOT/'results/scoring_history/previous_nine/DEVELOPMENT_SCORES.json')
    seal = read(ROOT/'audit/DEVELOPMENT_PREDICTIONS_SEALED_V3.json')
    unchanged = old['rows'] == scores['rows'][:9]
    rows = []
    for name in NAMES:
        directory = ROOT/'results'/name
        rec = read(directory/'RECEIPT.json')
        corr = read(directory/'CORRECTOR.json')
        reused = name == 'global_common_corrector'
        if reused:
            aux = rec['historical_auxiliary_execution']
            operations = aux['operation_counts']
            updates = operations['PR_CG_steps']
            history = read(ROOT.parent/'research_max_struct_r2/results/blind/obs_b20d5a55677d4eca/global/AUXILIARY.json')['history']
            stage_stops = [x['status'] for x in history if x.get('event') == 'stage_end']
            inner_stops = dict(Counter(x.get('material_inner_message', 'UNRECORDED')
                                      for x in history if 'global_iteration' in x))
        else:
            aux = read(directory/'AUXILIARY.json')
            operations = rec.get('auxiliary_operation_counts', rec.get('auxiliary_cost'))
            history = aux['history']
            updates = sum('global_iteration' in x for x in history)
            stage_stops = [x['status'] for x in history if x.get('event') == 'stage_end']
            inner_stops = dict(Counter(x.get('material_inner_message', 'UNRECORDED')
                                      for x in history if 'global_iteration' in x))
        rows.append(dict(
            arm=name, receipt=f'results/{name}/RECEIPT.json',
            physical_LU=rec['physical_model_cost']['LU'],
            physical_complex_RHS=rec['physical_model_cost']['actual_complex_RHS'],
            wall_seconds=rec['wall_seconds'], cpu_seconds=rec['cpu_seconds'],
            auxiliary_wall_seconds=operations['wall_seconds'],
            auxiliary_reused=reused, auxiliary_updates=updates,
            material_inner_evaluations=operations['material_inner_evaluations'],
            auxiliary_stage_stops=stage_stops, material_inner_stops=inner_stops,
            corrector_stop=corr['stop'], corrector_iterations=corr['iterations'],
            corrector_evaluations=corr['evaluations'], corrector_wall_seconds=corr['wall_seconds'],
            small4x4_factorizations=operations.get('small4x4_factorizations', 0),
            small4x4_solve_complex_RHS=operations.get('small4x4_complex_RHS', 0),
            small4x4_matrix_apply_complex_RHS=operations.get('small4x4_matrix_apply_complex_RHS', 0),
            scalar_complex_divisions=operations.get('scalar_complex_divisions', 0),
            onehop_block_matrix_complex_RHS=operations.get('onehop_block_matrix_complex_RHS', 0),
            operation_counts=operations, failed_attempts=rec.get('failed_attempts', []),
            timer_scope=rec['cost_scope'], truth_access_by_runner=rec['truth_access']))
    receipt_delta = max(abs(r['physical_objective_receipt_delta']) for r in scores['rows'])
    residual_delta = max(abs(r['relative_residual_receipt_delta']) for r in scores['rows'])
    result = dict(
        status='RECEIPT_AND_SCORE_AUDIT_COMPLETE', created_utc=datetime.now(timezone.utc).isoformat(),
        classification='POST_SEAL_PUBLIC_DEVELOPMENT', case_id='obs_b20d5a55677d4eca',
        previous_nine_rows_exactly_preserved=unchanged,
        seal_precedes_truth_access=seal['checked_utc'] < scores['first_truth_opened_utc'],
        seal_count=seal['count'], maximum_objective_receipt_delta=receipt_delta,
        maximum_relative_residual_receipt_delta=residual_delta,
        frozen_gates=scores['inherited_contract'], comparison_rows=rows,
        common_physical_objective='Phi=.5*||Y_prediction-Y||_F^2/||Y||_F^2+R(m)',
        common_prior='R=1e-4 mean(d^2)+1e-4 mean(sqrt(d^2+.01^2)-.01), d=nearest-neighbor differences/19.5',
        common_corrector=dict(source='experiments/cooperative_paths.py::refine',
            method='L-BFGS-B', objective_scale=1e6, maxiter=1000, maxfun=4000,
            ftol=1e-12, scaled_gtol=.01, original_gtol=1e-8, maxls=30,
            material_box=[.02,20.], incumbent='best original Phi among evaluated iterates including seed'),
        selection='Original physical Phi only, no GT-ranked reselection; global control additionally retains published incumbent',
        historical_global_auxiliary=read(ROOT/'results/global_common_corrector/RECEIPT.json')['historical_auxiliary_execution'],
        scientific_interpretation='Reserved for parent Codex thread',
        confirmation_GT_access=False, new_physical_LU=0, new_complex_RHS=0,
        hash_checks=False, new_gates=False,
        wall_seconds=time.perf_counter()-start, cpu_seconds=time.process_time()-cpu,
        timer_scope='Read saved score/seals/receipts/histories and build summaries; serialization excluded')
    (ROOT/'audit/DEVELOPMENT_COST_AUDIT.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    lines = [
        '# R4 开发评分与成本独立审查', '',
        '范围：仅历史 R2 公开开发案例 `obs_b20d5a55677d4eca`。分类为 **POST_SEAL_PUBLIC_DEVELOPMENT**；GT 已公开，不是盲测。科学解释及最终结论由父线程负责。', '',
        '## 封存与评分', '',
        '先保存原九端点结果、评分代码和既有 seal 至 `results/scoring_history/previous_nine/`，再新增四个已封存预测。`audit/DEVELOPMENT_PREDICTIONS_SEALED_V3.json` 登记全部 13 个预测的路径、大小、mtime 和 receipt 状态，写入后才读取公开开发 GT。未读取任何 confirmation GT。', '',
        f"原九行评分逐字段保持一致：`{unchanged}`。新 seal 时间 `{seal['checked_utc']}`，本轮首次开发 GT 读取时间 `{scores['first_truth_opened_utc']}`。", '',
        '沿用 R2 冻结评分函数 `score_blind.py::{otsu,iou,prior}`：RMSE/12 ≤ 0.15；预测自身 256-bin Otsu IoU ≥ 0.7；原观测 99% CN 白化残差门；公开材料边界及原数值容差。没有新增门，也没有重新选择预测。', '',
        '| 端点 | RMSE/12 | Otsu IoU | 相对残差 | 白化残差 | 原联合门 |',
        '|---|---:|---:|---:|---:|---|']
    wanted={'A_block_absorption','self_absorption','block_absorption_cached','native_extended','onehop_absorption','global_common_corrector'}
    for row in scores['rows']:
        if row['arm'] in wanted:
            lines.append(f"| {row['arm']} | {row['amplitude_normalized_nrmse']:.9f} | {row['otsu_iou']:.9f} | {row['relative_residual']:.12f} | {row['whitened_residual']:.6f} | {'PASS' if row['joint_gate'] else 'FAIL'} |")
    lines += ['', f"原 99% 白化半径为 {scores['inherited_contract']['noise_radius99']:.12f}。从保存的预测复核物理目标与 receipt 的最大差为 {receipt_delta:.3g}，相对残差差为 {residual_delta:.3g}；无需 forward。", '',
        '## 同 prior、物理目标及停止规则', '',
        '六个比较端点共用冻结标量 TM 模型、同一公开 Y、host=0.1、材料范围 [0.02,20]、历史 smooth=TV=1e-4 prior。原物理目标为 Phi=0.5||Ypred−Y||²/||Y||²+R(m)。R(m)=1e-4 mean(d²)+1e-4 mean(sqrt(d²+0.01²)−0.01)，d 为相邻材料差/19.5。runner 将同一个 `model.reg.evaluate(m)[:2]` 传入辅助读出；未用 GT 选端点。辅助变量中的状态残差目标随方法不同，不能因此称辅助目标相同。', '',
        '证据：`research_fbg_r3_experiments/src/physics.py::{MaxwellModel.physical_score,data_value_gradient}`；`research_z_alg_r4_diversity/src/physics.py::Regularizer`；各端点保存的配置、runner 源码、AUXILIARY 和 CORRECTOR receipt。', '',
        'A 直解、self、cached 和 onehop 使用 mf3/6/10 各 500 次全局 PR 更新；源指标、DPCS 和 beta 日程保持原设置。物料读出最多 60 次、ftol=1e-12、gtol=1e-8。self 仅保留对角吸收；cached 以逆缓存/矩阵乘代替小块重复求解，浮点轨迹不同；onehop 保留一次跨格残差交换。native_extended 保留原 sourcewise CIE，额外执行 mf12 的 6500 次更新，合计 8000 次，其辅助计算预算不同。', '',
        'global_common_corrector 复用历史 mf3/6/10 全局 CIE 辅助端点；原历史终点也作为 incumbent 保留。各方法都用原 Phi 排序。所有六个端点调用同一 `cooperative_paths.refine`：L-BFGS-B，目标乘 1e6，maxiter=1000、maxfun=4000、ftol=1e-12、缩放后 gtol=0.01（原目标 1e-8）、maxls=30。全部实际在 1000 次迭代上限停止。此记录保留预算终止状态，不认证收敛或局部最优。', '',
        '| 端点 | 辅助更新 | 内层物料评估 | corrector评估 | corrector停止 |',
        '|---|---:|---:|---:|---|']
    for row in rows:
        lines.append(f"| {row['arm']} | {row['auxiliary_updates']}{'（历史复用）' if row['auxiliary_reused'] else ''} | {row['material_inner_evaluations']} | {row['corrector_evaluations']} | MAXITER 1000 |")
    lines += ['', '辅助各阶段均记录 MAXITER；内层读出的具体停止计数保存在 `DEVELOPMENT_COST_AUDIT.json`，不将达到外层日程称为收敛。', '',
        '## 实际成本', '',
        '| 端点 | 新物理 LU | 新物理 complex RHS | 进程 wall s | 进程 CPU s | 辅助 wall s |',
        '|---|---:|---:|---:|---:|---:|']
    for row in rows:
        suffix='（历史）' if row['auxiliary_reused'] else ''
        lines.append(f"| {row['arm']} | {row['physical_LU']} | {row['physical_complex_RHS']} | {row['wall_seconds']:.6f} | {row['cpu_seconds']:.6f} | {row['auxiliary_wall_seconds']:.6f}{suffix} |")
    lines += ['',
        'global 的新 wall 仅含本轮加载、排序、公共 corrector 与保存；另有历史辅助 5.181085 s、3 LU/36 RHS，不重复计入新成本。它还保留旧全球终点；完整历史全球运行本身为 689 LU/16452 RHS、wall 9.220951 s，见 global receipt，不能将其全部算作新的辅助成本。其余完整新运行包含 imports/公共数据读取、辅助执行、物理排序、corrector 与产物写入；不含未重跑的历史算子组装，最终 receipt 自身写入除外。', '',
        '| 端点 | 局部工作（与物理 LU/RHS 分开） |', '|---|---|']
    for row in rows:
        work=[]
        if row['small4x4_factorizations']:
            work.append(f"4x4 分解 {row['small4x4_factorizations']}；求解 RHS {row['small4x4_solve_complex_RHS']}")
        if row['small4x4_matrix_apply_complex_RHS']:
            work.append(f"逆矩阵应用 RHS {row['small4x4_matrix_apply_complex_RHS']}")
        if row['scalar_complex_divisions']:
            work.append(f"标量复除法 {row['scalar_complex_divisions']}")
        if row['onehop_block_matrix_complex_RHS']:
            work.append(f"onehop 块矩阵应用 RHS {row['onehop_block_matrix_complex_RHS']}")
        if row['arm']=='native_extended':
            work.append('source-mode-iterations 48408000；该代理量不是 FLOPs 或等时指标')
        if row['auxiliary_reused']:
            work.append('复用辅助端点；当前不重跑局部操作')
        lines.append(f"| {row['arm']} | {'；'.join(work)} |")
    lines += ['',
        'direct 与 cached 的物理设置、辅助日程和 corrector 上限相同，但求解顺序改变了轨迹、内层次数与 corrector 评估数。这里保留每次实际 wall 及工作量，未做重复计时，也未归因全部 wall 差。不同方法的局部工作、物料内层次数和 native 额外容量不同，因此不是等算力比较。', '',
        '六个运行 receipt 均记录 failed_attempts=[]。既有启动失败位于原审计记录；本次评分未启动求解器、未重跑检查，评分脚本本身新增 0 LU / 0 RHS，无 hash 检查，无确认集 GT 读取。费用/货币成本未提供，保留为未知。', '',
        '本审查输出是评分、协议一致性与成本事实，不替代父线程对新颖性、机制、成像收益或外推范围的判断。']
    (ROOT/'audit/DEVELOPMENT_AUDIT.md').write_text('\n'.join(lines)+'\n')
    path=ROOT/'results/DEVELOPMENT_SCORES.md'
    pointer='\nFull receipt-only prior/stopping/cost comparison: [DEVELOPMENT_AUDIT.md](../audit/DEVELOPMENT_AUDIT.md); machine-readable costs: [DEVELOPMENT_COST_AUDIT.json](../audit/DEVELOPMENT_COST_AUDIT.json).\n'
    if pointer.strip() not in path.read_text():
        path.write_text(path.read_text()+pointer)
    print(json.dumps(dict(status=result['status'], previous_nine_unchanged=unchanged,
        seal_precedes_truth_access=result['seal_precedes_truth_access'], arms=len(rows),
        new_physical_LU=0,new_complex_RHS=0)))


if __name__ == '__main__':
    summarize()
