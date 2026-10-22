"""Render truthful current summaries only: no solver, training or test truth."""
from pathlib import Path
import sys,json,html,datetime
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r2_risk_posterior_exploration'
def read(path):
    p=OUT/path
    return json.loads(p.read_text()) if p.exists() else None
def main():
    ledger=read('RESOURCE_LEDGER.json');strong=read('A/training/RESULTS.json');density=read('A/density_training/SELECTED.json');weak=read('A/weak_density_training/SELECTED.json');C=read('C/analysis/C_exposed_eval_analysis/RESULTS.json')
    rows=[]
    scale2=read('A/training32768/RESULTS.json');phase=read('FRESH_CHECKPOINT.json');freeze=read('METHODS_FROZEN_FOR_LOCKED.json')
    residual=read('A/residual_training/SELECTED.json');fusion=read('D/fusion/FROZEN.json')
    for r in (strong or {}).get('results',[])+(scale2 or {}).get('results',[]):rows.append((r['train_size'],r['kind'],r['seed'],r['selected_update'],r['validation']['complex_roi_mse']*256/576))
    sections=['# P3L-R2：材料风险、非 Gaussian 推断与真实轨迹','状态：**RUNNING**。本页是进行中的实际实验，尚未完成新的锁定泛化，不是最终 GO 结论。','原 P1、观测与噪声不变。监督风险、概率与几何三条路线独立推进；不合格 MCMC 不阻断神经或风险策略。','## 已完成证据','已生成 16,384 个独立 train 与 2,048 个 validation；与历史封账严格分开。旧128 locked已暴露，新256 locked尚未生成。']
    sections[-1]=('已完成32,768个独立train与2,048个validation的联合模拟和数据规模学习；512train+128validation有独立物理上下文。旧128 locked已暴露。' if scale2 else sections[-1])
    if freeze:sections+=['77个checkpoint、原P1与256个fresh parent已在第一份新观测产生前本地冻结，source revision `'+freeze['source_revision']+'`。']
    if phase:sections+=['本次新测试当前阶段：`'+phase.get('task',phase.get('status','UNKNOWN'))+'`；完成前不报告新泛化成绩。']
    if residual and fusion:
        r=residual['result'];sections+=['512训练对象的局部CNN为%d参数；validation128上，物理起点MSE %.7f，选定材料风险动作MSE %.7f。这个差异经过validation选择，不能当作新测试成绩。'%(r['parameter_count'],r['anchor_validation']['complex_roi_mse']*256/576,fusion['selected']['validation_MSE'])]
    if density:
        r=density['result'];sections+=['固定10维联合模拟密度：validation NLL %.5f（无条件 %.5f），最佳更新 %d。它是归一化 proposal/坐标预测，不是已认证 FULL posterior。'%(r['validation_nll'],r['unit_prior_validation_nll'],r['selected_update'])]
    if weak:
        r=weak['result'];sections+=['材料风险弱 frame 密度：validation NLL %.5f（无条件 %.5f）。不同坐标的联合NLL不能当成相同材料概率目标。'%(r['validation_nll'],r['unit_prior_validation_nll'])]
    sections+=['## 直接材料风险：已完成的 validation 比较','表中为全网格 complex MSE；ROI MSE已乘256/576。未结束的训练不算完成。','| train parents | 模型 | seed | 最佳update | validation MSE |','|---:|---|---:|---:|---:|']
    for n,k,s,u,m in rows:sections.append(f'|{n}|{k}|{s}|{u}|{m:.7f}|')
    if C:
        sections+=['## 已暴露8对象的真实轨迹','这是机制开发证据，不是新盲测；轨迹为 candidate-only，不是 posterior samples。','|方法|最终MSE|相对起点差|实付RHS（未含共享anchor setup）|','|---|---:|---:|---:|']
        for k,r in C['by_arm'].items():sections.append(f"|{k}|{r['mean_final_MSE']:.7f}|{r['mean_final_minus_anchor_MSE']:+.7f}|{int(r['measured_RHS'])}|")
    sections+=['## 资源与边界','新上限20,000,000 complex RHS / 24h受控CPU / 24h GPU；本轮当前GPU0，不使用付费外部计算。']
    if ledger:
        sections+=['|方向|累计RHS|state attempts|已结算受控秒|','|---|---:|---:|---:|']
        for k,v in ledger['stages'].items():sections.append(f"|{k}|{v['RHS']}|{v['states']}|{v['controlled_wall_seconds']:.2f}|")
        sections.append('进行中的作业另有预留，不把尚未结算的wall记作零成本。')
    sections+=['[概率与数学合同](THEORY_AND_METHODS.md) · [冻结方法](METHODS_FROZEN_FOR_LOCKED.json) · [迭代记录](ITERATION_LOG.md) · [资源快照](RESOURCE_LEDGER.json) · [原P3-L封账报告](../p3l_local_posterior_imaging/)','下一步：完成256个全新parent上的固定在线推断、轨迹与独立评分；不再调参或重新选择模型。FULL posterior仍不能凭这些结果称为合格。']
    text='\n\n'.join(sections)+'\n';(OUT/'PROGRESS_ZH.md').write_text(text)
    # Plain source is intentionally visible without remote math/runtime assets.
    document='<!doctype html><html lang="zh"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>P3L-R2 研究进行中</title><style>body{margin:40px auto;max-width:1100px;padding:0 20px;font:17px/1.65 system-ui;color:#183146;background:#f7fafc}pre{white-space:pre-wrap;font:15px/1.75 ui-monospace,monospace;background:white;padding:26px;border-radius:12px}a{color:#066b91}</style><h1>P3L-R2 · 研究进行中</h1><p>监督材料风险 · 非 Gaussian posterior · 几何轨迹</p><p><a href="PROGRESS_ZH.md">匿名 Markdown 进度</a> · <a href="THEORY_AND_METHODS.md">数学合同</a> · <a href="../p3l_local_posterior_imaging/">旧P3-L结果</a></p><pre>'+html.escape(text)+'</pre></html>'
    (OUT/'index.html').write_text(document);(OUT/'report.html').write_text(document)
    (OUT/'START_HERE.md').write_text('# P3L-R2：独立材料风险与非 Gaussian 研究\n\n状态：RUNNING，未完成新的锁定泛化。\n\n[当前实际进度](report.html) · [Markdown](PROGRESS_ZH.md) · [理论合同](THEORY_AND_METHODS.md) · [资源快照](RESOURCE_LEDGER.json)\n\n本轮不覆盖旧数据，不把不合格参考当teacher，也不要求FULL参考通过后才研究监督风险和轨迹。原P1/noise/Maxwell固定。\n')
    print('progress rendered from saved summaries;0newphysics,0GT')
if __name__=='__main__':main()
