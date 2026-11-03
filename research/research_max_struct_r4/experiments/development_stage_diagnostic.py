"""OFFLINE_DIAGNOSTIC: disclosed development GT and saved stages; no physics."""
from pathlib import Path
from datetime import datetime, timezone
import json
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
R2=ROOT.parent/'research_max_struct_r2'
sys.path.insert(0,str(REPO))
from research.research_max_struct_r2.experiments.score_blind import otsu, iou, prior

CASE='obs_b20d5a55677d4eca'
ARMS=('onehop_absorption','block_absorption','self_absorption','global_common_corrector')


def read(p):
    return json.loads(p.read_text())


def run():
    started,cpu=time.perf_counter(),time.process_time()
    truth_opened=datetime.now(timezone.utc).isoformat()
    with np.load(R2/'sealed_truth'/f'{CASE}.npz',allow_pickle=False) as d:
        truth=d['m'].copy()
        support=d['mask'].ravel().astype(bool)
    with np.load(R2/'data_public'/f'{CASE}.npz',allow_pickle=False) as d:
        host=float(d['host'])
    scores=read(ROOT/'results/DEVELOPMENT_SCORES.json')
    rows=[]

    def append(arm,label,updates,m,material_source,physical,physical_source):
        if m is None:
            nrmse=overlap=threshold=reg=None
        else:
            mask,threshold=otsu(m)
            nrmse=float(np.linalg.norm(m-truth)/np.sqrt(truth.size)/12.)
            overlap=iou(mask,support)
            reg=prior(m)
        relative=physical.get('physical_relative_residual',physical.get('relative_residual'))
        data=physical.get('physical_nll')
        if data is None and relative is not None:
            data=.5*float(relative)**2
        phi=physical.get('objective',physical.get('physical_objective'))
        rows.append(dict(arm=arm,label=label,auxiliary_updates=updates,
            classification='OFFLINE_DIAGNOSTIC',amplitude_normalized_nrmse=nrmse,
            otsu_iou=overlap,otsu_threshold=threshold,relative_physical_residual=relative,
            data_Phi=data,total_original_Phi=phi,material_prior=reg,
            total_Phi_minus_data_and_recomputed_prior=(phi-data-reg if all(v is not None for v in (phi,data,reg)) else None),
            material_source=material_source,physical_metric_source=physical_source,
            missing_physical_reason=None if relative is not None else 'No saved exact physical residual; forward prohibited'))

    for arm in ARMS:
        global_arm=arm=='global_common_corrector'
        folder=(R2/'results/blind'/CASE/'global') if global_arm else ROOT/'results'/arm
        aux=read(folder/'AUXILIARY.json')
        candidates={int(x['stage']):x for x in aux['candidates']}
        append(arm,'host',0,np.full(truth.shape,host),'public host scalar',candidates.get(-1,{}),
               str((folder/'AUXILIARY.json').relative_to(REPO))+'::candidates[stage=-1]')
        for k in range(3):
            path=folder/f'stage_{k}.npz'
            if path.exists():
                with np.load(path,allow_pickle=False) as d:
                    m=d['m'].copy()
            else:
                m=None
            append(arm,f'mf{(3,6,10)[k]}_stage_end',500*(k+1),m,str(path.relative_to(REPO)),
                   candidates.get(k,{}),str((folder/'AUXILIARY.json').relative_to(REPO))+f'::candidates[stage={k}]')
        final=next(x for x in scores['rows'] if x['arm']==('A_block_absorption' if arm=='block_absorption' else arm))
        path=ROOT/final['prediction']
        with np.load(path,allow_pickle=False) as d:
            m=d['m'].copy()
        append(arm,'common_corrector_endpoint',1500,m,str(path.relative_to(REPO)),final,
               'research/research_max_struct_r4/results/DEVELOPMENT_SCORES.json::'+final['arm'])
    historical=read(R2/'results/blind'/CASE/'global/RECEIPT.json')
    path=R2/'predictions'/CASE/'global.npz'
    with np.load(path,allow_pickle=False) as d:
        m=d['m'].copy()
    append('global_common_corrector','historical_R2_corrector_reference',1500,m,str(path.relative_to(REPO)),historical,
           str((R2/'results/blind'/CASE/'global/RECEIPT.json').relative_to(REPO)))
    change=[]
    for arm in ARMS:
        before=next(x for x in rows if x['arm']==arm and x['label']=='mf10_stage_end')
        after=next(x for x in rows if x['arm']==arm and x['label']=='common_corrector_endpoint')
        change.append(dict(arm=arm,
            pre_corrector_nrmse=before['amplitude_normalized_nrmse'],
            post_corrector_nrmse=after['amplitude_normalized_nrmse'],
            pre_corrector_iou=before['otsu_iou'],post_corrector_iou=after['otsu_iou'],
            pre_corrector_relative_residual=before['relative_physical_residual'],
            post_corrector_relative_residual=after['relative_physical_residual'],
            pre_corrector_total_Phi=before['total_original_Phi'],post_corrector_total_Phi=after['total_original_Phi']))
    result=dict(status='OFFLINE_DIAGNOSTIC',case_id=CASE,
        classification='GT_DISCLOSED_PUBLIC_DEVELOPMENT_NOT_BLIND',
        created_utc=datetime.now(timezone.utc).isoformat(),public_development_truth_opened_utc=truth_opened,
        rows=rows,pre_post_corrector=change,
        physical_metric_policy='Only saved exact physical candidate-validation receipts; no use of auxiliary J/Q as an exact Maxwell state',
        data_Phi_definition='.5*relative_physical_residual**2; saved physical_nll used when available',
        material_metrics='unchanged R2 RMSE/12 and prediction-only Otsu256 IoU; descriptive, no new gate',
        online_algorithm_changed=False,candidate_reselection=False,new_iterations=0,
        diagnostic_failed_attempts=[dict(error='SyntaxError: unmatched closing parenthesis',
            disposition='fixed before successful diagnostic execution',truth_access=False,new_physical_LU=0)],
        new_physical_LU=0,new_complex_RHS=0,new_tests=False,hash_checks=False,
        confirmation_GT_access=False,scientific_judgment='Reserved for parent Codex thread',
        wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu,
        timer_scope='Read stored development truth/material/receipts and calculate offline metrics; serialization excluded')
    out=ROOT/'audit/DEVELOPMENT_STAGE_DIAGNOSTIC.json'
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    lines=['# R4 开发阶段离线诊断','',
        '**OFFLINE_DIAGNOSTIC / GT_DISCLOSED_PUBLIC_DEVELOPMENT_NOT_BLIND**。仅历史公开案例 `obs_b20d5a55677d4eca`；读取公开开发 GT 后计算各保存阶段材料误差和 IoU。未更改冻结算法，未选新端点，未增加迭代、测试或 forward，未读 confirmation GT。','',
        '三条吸收路径使用各自保存的 `stage_0/1/2.npz`；成熟 global 复用 R2 历史相同阶段。物理残差和数据 Phi 来自对应 `AUXILIARY.json` 已执行的精确物理候选验证；保存的辅助 J/Q 不被当成材料的精确 Maxwell 解。数据 Phi=0.5×相对物理残差²；总 Phi 还含原 prior。缺失值保留 null。','',
        '| 路径 | 阶段 | RMSE/12 | Otsu256 IoU | 物理相对残差 | 数据 Phi | 总 Phi（含 prior） |',
        '|---|---|---:|---:|---:|---:|---:|']
    def fmt(x):
        return 'null' if x is None else f'{x:.9g}'
    for row in rows:
        lines.append('| '+ ' | '.join([row['arm'],row['label'],*(fmt(row[k]) for k in
            ('amplitude_normalized_nrmse','otsu_iou','relative_physical_residual','data_Phi','total_original_Phi'))])+' |')
    lines+=['','## 进入公共 corrector 之前与之后','',
        '| 路径 | RMSE/12 前→后 | IoU 前→后 | 物理残差 前→后 | 总 Phi 前→后 |',
        '|---|---:|---:|---:|---:|']
    for r in change:
        lines.append(f"| {r['arm']} | {r['pre_corrector_nrmse']:.9f} → {r['post_corrector_nrmse']:.9f} | {r['pre_corrector_iou']:.9f} → {r['post_corrector_iou']:.9f} | {r['pre_corrector_relative_residual']:.9f} → {r['post_corrector_relative_residual']:.9f} | {r['pre_corrector_total_Phi']:.9g} → {r['post_corrector_total_Phi']:.9g} |")
    before={r['arm']:r for r in rows if r['label']=='mf10_stage_end'}
    one=before['onehop_absorption'];block=before['block_absorption'];ownself=before['self_absorption'];glob=before['global_common_corrector']
    lines+=['',
        f"已保存的 mf10 辅助端点：onehop 的 RMSE/12={one['amplitude_normalized_nrmse']:.9f}、IoU={one['otsu_iou']:.9f}、残差={one['relative_physical_residual']:.9f}；exact block 为 {block['amplitude_normalized_nrmse']:.9f}、{block['otsu_iou']:.9f}、{block['relative_physical_residual']:.9f}；self 为 {ownself['amplitude_normalized_nrmse']:.9f}、{ownself['otsu_iou']:.9f}、{ownself['relative_physical_residual']:.9f}；global 为 {glob['amplitude_normalized_nrmse']:.9f}、{glob['otsu_iou']:.9f}、{glob['relative_physical_residual']:.9f}。因此路径之间在公共 corrector 之前已经存在可量化的材料/数据端点差异；上表另列 corrector 后的进一步变化。",'',
        '这是同一公开开发对象的事后阶段诊断。阶段之间同时推进 Fourier 容量和迭代，四条路径的辅助残差度量/物料读出不同；这里不把关联解释为单一因素的因果证明，不构成盲测、额外成功率或收敛认证。全部公共 corrector 原先都在 1000 迭代上限停止。最终科学解释由父线程负责。','',
        'JSON 保存逐行材料路径、物理指标来源与 total Phi−data Phi−prior 差值。此次新增物理成本 0 LU / 0 RHS、0 新迭代、0 新测试；无 hash 检查。离线脚本首次启动有一处括号语法错误，执行前即终止（未读 GT、0 物理成本），修正后完成诊断。']
    (ROOT/'audit/DEVELOPMENT_STAGE_DIAGNOSTIC.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(status='OFFLINE_DIAGNOSTIC',rows=len(rows),pre_post_corrector=change,new_physical_LU=0,new_complex_RHS=0)))


if __name__=='__main__':
    run()
