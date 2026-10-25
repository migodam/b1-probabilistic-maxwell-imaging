"""Build final Chinese reports from CLOSED R4 evidence; no experiment or fit.

The narrative and scientific decisions are reviewed by the parent Codex agent.
This tool only reads saved scores, copies compact public evidence and renders it.
"""
from pathlib import Path
import argparse
import collections
import csv
import json
import shutil
import time

ROOT = Path(__file__).resolve().parents[1]
S = ROOT / 'research/p3l_r4_conditional_transport'
PUBLIC = S / 'evidence_summary'


def read(p):
    return json.loads((S / p).read_text())


def write(p, x):
    p = S / p
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(x, ensure_ascii=False, indent=2) + '\n')


def md(p, x):
    (S / p).parent.mkdir(parents=True, exist_ok=True)
    (S / p).write_text(x.strip() + '\n')


def mean_families(score):
    d = collections.defaultdict(list)
    for model, row in score['methods'].items():
        for action in ('raw', 'clipped'):
            d[(row['family'], action)].append(row[action])
    out = {}
    for (family, action), rows in d.items():
        out.setdefault(family, {})[action] = {
            k: sum(r[k]['mean'] for r in rows) / len(rows)
            for k in ('full_mse', 'full_mae', 'roi_mse', 'roi_mae')
        }
        out[family]['seeds'] = len(rows)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--overwrite-final-documents', action='store_true',
                        help='Explicitly rebuild the original closed draft; post-pack publication and interpretation additions are not regenerated.')
    args = parser.parse_args()
    if (S / 'THEORY_CALIBRATION_ERROR_BOUND.md').exists() and not args.overwrite_final_documents:
        parser.error('Final publication contains additive post-pack supplements. Use tools/replay_r4_public_risk.py for numerical replay, or explicitly rebuild the original draft in a separate checkout with --overwrite-final-documents.')
    assert read('results/FRESH_ISOLATED_SCORES_COMPLETE.json')['status'] == 'COMPLETE'
    PUBLIC.mkdir(exist_ok=True)
    risk, density, risk_tables = {}, {}, {}
    for group in ('locked', 'confirmation'):
        r = read(f'results/fresh/{group}/risk45/scores/RESULTS.json')
        d = read(f'results/fresh/{group}/density_scores/RESULTS.json')
        assert all(not rows for rows in r['all_failures'].values())
        assert r['n_parent'] == (512 if group == 'locked' else 256)
        risk_tables[group] = mean_families(r)
        risk[group] = dict(n_parent=r['n_parent'], families=risk_tables[group],
                           comparisons=r['comparisons'], all_failures=r['all_failures'])
        density[group] = dict(n_parent=d['n_parent'], sample_count=d['sample_count'],
                             families=d['families'], comparisons=d['comparisons'],
                             posterior_qualification=d['posterior_qualification'])
        write(f'evidence_summary/RISK_{group.upper()}.json', risk[group])
        write(f'evidence_summary/DENSITY_{group.upper()}.json', density[group])
    ppc = read('results/fresh/full_probability/scores/RESULTS.json')
    importance = [r['physical_summary']['importance']['ESS'] for r in ppc['records']
                  if r['physical_summary'].get('importance', {}).get('status') == 'EVALUATED']
    ppc_small = dict(status=ppc['status'], Gamma256=ppc['Gamma256'], costs=ppc['costs'],
                     importance_K64=dict(parents=len(importance), minimum=min(importance),
                                        median=sorted(importance)[len(importance)//2], maximum=max(importance)),
                     no_residual_selection=ppc['no_residual_selection'],
                     posterior_qualification=ppc['posterior_qualification'])
    write('evidence_summary/FULL_PROBABILITY.json', ppc_small)
    physical_risk = {}
    for group in ('locked', 'confirmation'):
        v = read(f'results/fresh/{group}/raw_forward/RESULTS.json')
        by = collections.defaultdict(list)
        for row in v['rows']:
            by[row['action']].append(row)
        physical_risk[group] = {
            k: dict(n=len(rows), data_feasible=sum(r.get('Phi', float('inf')) <= 282.873792 for r in rows),
                    mean_Phi=sum(r['Phi'] for r in rows)/len(rows),
                    prior_realizability='NOT_PROVEN', role='material risk action')
            for k, rows in by.items()
        }
    write('evidence_summary/RISK_VS_PHYSICS.json', physical_risk)
    chain = read('results/IV/interval_chain_locked_v1/COMPARISON.json')
    write('evidence_summary/DECISION_CHAIN_COMPARISON.json', chain)
    rep = read('results/fresh/locked/representative_scores_ndarray_recovery/RESULTS.json')
    rep_small = dict(status=rep['status'], physical_cost=rep['physical_cost'], comparisons=[
        {k: v for k, v in r.items() if k != 'MSE_change_candidate_minus_reference'}
        for r in rep['comparisons']])
    write('evidence_summary/REPRESENTATIVES.json', rep_small)
    ledger = read('RESOURCE_LEDGER.json')
    totals = {k: sum(d[k] for d in ledger['stages'].values()) for k in
              ('RHS', 'states', 'LU', 'controlled_wall_seconds', 'physics_wall_seconds', 'GPU_seconds')}
    decision = {
        'GO_A': 'GO: replicated supervised material-risk improvement',
        'GO_B_relative': 'GO: explicit full66 proper-score improvement; static physical context contributes',
        'GO_B_calibrated_FULL': 'NOT_QUALIFIED: fixed-draw Maxwell Gamma test and importance concentration fail',
        'GO_C': 'GO_LIMITED: float64 same-kernel decision acceleration; not a continuous-PDE certificate or mixing gain',
        'GO_D_static': 'CONDITIONAL_GO: static Maxwell learning metrics improve selected independent scores',
        'GO_D_dynamic_nonholonomic': 'NO_GO_CURRENT_METHOD: no matched-cost material-statistical increment',
        'GO_E_strong_TAP': 'NOT_YET: no established new Maxwell-specific theorem/algorithm with full probabilistic and scaling evidence',
    }
    metrics = dict(format='p3lr4_final_evidence_v1', updated_ns=time.time_ns(),
                   target='P1-CORRELATED-HOST-MIXTURE-TM24-v1-20261018',
                   sigma=0.00119211798514519, new_test_parents=[512, 256],
                   exposed_after_evaluation=True, models_frozen_before_new_observations=True,
                   independent_statistical_unit='parent; average seed losses, not ensemble images',
                   risk=risk, density=density, full_probability=ppc_small,
                   risk_vs_physics=physical_risk, decision_chain=chain,
                   representative=rep_small, resources=dict(stages=ledger['stages'], totals=totals),
                   decisions=decision,
                   geometry_candidate_archive=dict(development_observations=8,
                       counts_including_anchor=[15,12,7,16,28,16,15,10], new_candidates=111,
                       minimum_pair_material_RMS=0.02, probability=None,
                       role='finite-noise candidate discovery; not posterior modes'),
                   new_3D=0, new_Mie=0, paid_cloud=False,
                   source_scores='closed original scores and full arrays retained; this file is a compact derivative')
    write('RESULTS_AND_METRICS.json', metrics)
    write('GATE_SUMMARY.json', decision)
    for name in ('FAMILY_SEED_MEAN_TABLE.csv', 'PRIMARY_FIXED_COMPARISONS.csv', 'ALL45_METHOD_ACTIONS.csv'):
        shutil.copy2(S/'results/fresh/independent_risk_statistical_audit_v1'/name, PUBLIC/name)
    md('START_HERE.md', '''
# Research P3L-R4：从材料风险改善到完整条件概率的实测边界

本轮真实训练、Maxwell 计算与两批新对象评价已完成。**监督成像 GO；完整后验仍未合格；固定核的数值接受决定加速 GO（限定范围）；动态非完整探索没有独立增量。**

512/256 个全新 parent 上，主要风险估计的全网格 complex MSE 为 **0.006806 / 0.007378**；冻结 R2 CNN 为 **0.007832 / 0.008323**，分别改善 **13.10% / 11.36%**。与本轮重新训练、数据更多的强 CNN 相比仍改善 **5.36% / 4.90%**。这些是监督风险动作，不是经过认证的 posterior mean。

[完整中文报告](RESEARCH_REPORT_ZH.md) · [网页](report.html) · [图库](gallery.html) · [科研裁决](GO_NO_GO.md) · [机器指标](RESULTS_AND_METRICS.json) · [复现说明](REPRODUCE.md) · [证据包入口](EVIDENCE_PACKAGE.md)

建议按以下顺序阅读：

1. [材料风险结果](MATERIAL_RISK_RESULTS.md)：为什么很差的全局生成均值仍能帮助局部 CNN，以及普通 Gaussian 控制为何重要。
2. [完整条件概率](CONDITIONAL_FLOW_RESULTS.md)和[校准审计](PROBABILITY_CALIBRATION_AUDIT.md)：真正的相对概率进步，以及尚未解决的严重 Maxwell 失配。
3. [物理计算复用](MAXWELL_REUSE_RESULTS.md)：相同 Full-MH 路径上的实际节省，以及缓存、近似和数值范围。
4. [几何实验](NONHOLONOMIC_GEOMETRY_RESULTS.md)：有限材料候选确实存在，四段路径却没有胜过直接提案。
5. [理论与方法](THEORY_AND_ALGORITHMS.md)、[完整失败记录](FAILURE_ANALYSIS.md)、[下一轮 Pro 交接](GPT_PRO_HANDOFF/START_HERE.md)。

历史 P1/P2/P3-L/R2 不改写。新对象在冻结方法和两批预测完成后才进入独立 GT 评分；现在两批标签均为 EXPOSED，不能再称未来盲测。原始失败、重复状态、拒绝、过时中期日志和费用仍完整保留。
''')
    md('RESEARCH_REPORT_ZH.md', '''
# P3L-R4 完整研究报告：条件材料输运、风险成像与 Maxwell 自适应计算

## 1. 这轮实际解决了什么

我们完成了独立的训练、物理与泛化实验，而不是把 R3 理论重新整理成计划。结果需要分开看。

材料成像有可重复进步：新方法在 512 和 256 个全新对象上，比历史 R2 CNN 分别降低约 13.10% 和 11.36% 的材料 MSE；比本轮更强、训练数据更多的 CNN 降低约 5.36% 和 4.90%。但普通 Gaussian 产生的全局材料统计与 Spline Flow 一样有效，不能把这一成像收益归因于 Flow 已经学会正确 posterior。

完整 66 维条件密度也有真实的相对进步。Spline 在独立联合样本上的 proper NLL 优于同容量 diagonal Gaussian；固定 Maxwell frame 的三角模型优于随机 frame。有界的 Maxwell 几何 CFM 使材料 energy score 比同谱随机控制改善约 21%。这些是相对分布质量证据，不能抹掉绝对概率失败：固定生成样本的 Maxwell 残差仍比原噪声模型预期大一至两个数量级，重要性权重仍然退化。

物理计算方面，在 8 个全新对象、相同提案和随机数的 8,192 步重放中，收紧数值势区间把 Full 的 71.40 秒降到 24.15 秒，逐步决定和状态路径零差异。它没有改善这个短步长核的混合，也不是连续 Maxwell 的严格认证。普通 GMRES 同样带来收益，A17 在这里没有独立优势。

动态非完整几何没有建立独立算法收益。确实保存了 8 份观测下 111 个新的任务尺度、有限噪声竞争候选；它们来自多种普通和几何方法的共同池。它们不是 111 个 posterior modes，也没有概率权重。四段路径在距离保持与成本上弱于直接方向提案。

## 2. 物理与统计问题始终相同

本轮使用原 P1：24×24 scalar-TM/VIE、256 个材料 ROI cells、8 次照明和 32 个复接收通道。前 64 个独立 Gaussian latent 经相关场矩阵生成材料，后两个 Gaussianized 属性经标准正态 CDF 回到原 uniform 属性。权威 decoder 来自冻结仓库；没有用摘要重新写一个相似模型。

原标准差为 0.00119211798514519。记 Phi=||F(G(x))−y||²/sigma²，posterior 相对于 dx66 为 exp(−Phi−||x||²/2)。复数据沿 source-major 顺序打包，每个 source 先 32 个实部再 32 个虚部；实白化系数为 sqrt(2)/sigma。所有照明共享同一个被动复材料。没有改 noise、gain、采集、材料 prior 或物理后端。

decoder 在两侧精确饱和，多对一映射允许退化材料层；所以 latent density 与材料像素 Lebesgue density 是不同对象。本轮完整密度始终在 66 维原测度上定义，不把 material Gram 行列式误加进 prior。

## 3. 数据与独立性

恢复并复用了历史 32,768 个 joint TRAIN 和 2,048 个 validation 对象。历史昂贵 physics context 只有 512/128 子集，不被记成 32,768 个物理锚点。本轮另外付费形成 1,024/256 的强 context 训练与验证集。

新生成 128 个 development parent，随后在方法、checkpoint、NFE、样本数、失败规则和评分合同冻结后，生成 512 个第一批与 256 个第二批新 parent。两批使用相同方法；第一批结果没有用来重选第二批算法。两批所有风险与密度在线输出、Full PPC 和优化诊断封存后，才开启独立 GT evaluator。

风险面板包含 45 个固定方法实例，密度面板包含 26 个固定模型实例；二者在全部 768 个对象上均无推断失败或 fallback。随机种子不是独立对象。主表先对每个 parent 的三个 seed loss 求平均，再以 parent 做 4,000 次 paired bootstrap；它不是三模型集成图的成绩。预注册主比较使用 98.333% 区间，其他比较属于描述/探索，未宣称 45 方法的全族错误率控制。

## 4. 材料风险：强基线之后仍有收益

|全网格 complex MSE|新512对象|新256对象|
|---|---:|---:|
|独立 Probit physics point|0.01091244|0.01133359|
|冻结 R2 CNN|0.00783202|0.00832320|
|同一 5,842 参数类别重新训练，旧512 context|0.00741971|0.00798060|
|更强 CNN，1024 context、真实物理谱|0.00719170|0.00775824|
|Spline 全局矩 + 局部风险 head，预冻结主要方法|**0.00680623**|**0.00737782**|
|普通 diagonal Gaussian 全局矩 + 同风险 head|0.00682170|0.00739262|

对 R2 的相对改善为 13.10% [11.37%,14.89%] 和 11.36% [8.84%,13.96%]；对强 CNN 为 5.36% [4.02%,6.91%] 和 4.90% [3.14%,6.72%]。Spline 与 Gaussian fusion 仅差约 0.23%/0.20%，区间跨零，不能声称 Spline 在风险上获胜。

主要模型不是“只有 6,418 参数的完整系统”：6,418 是局部风险 head，另有约 0.68M 参数的冻结全局 provider，K=32 材料样本和原 physics context 的取得费用。这个成本需要与单独 R2 CNN 区分。

固定物理谱相对同容量随机空间控制的风险增量约 1.20%/0.83%，小但两批方向一致。它是静态 context 的贡献，既不是动态 weak flow 的贡献，也没有增加观测信息。

## 5. 为什么劣质全局生成均值仍然有用

完整 Spline 模型的直接样本均值，第一批 MSE 约 0.02278，远弱于 CNN。把这个图直接替换重建是坏方法；把它的均值和方差作为四个全局 context 通道，交给带物理基础图的局部 CNN，却有价值。

这与标准条件期望正交性一致：局部 patch 未必包含整个观测的所有任务信息，另一份仅由 y 计算的粗糙全局统计可以补充它。保持相同 head、打乱这四个通道的空间对应关系后收益下降；普通 Gaussian 或直接 field-moment provider 的控制又表明，收益主要来自有用的全局上下文，而不是 posterior 已校准。

按基础图的相区做离线分层，第一批 R2→新方法的全图风险下降中，过渡相贡献 0.00102460，其他单元贡献 0.00000118；第二批为 0.00094420 与 0.00000117。约 99.9% 落在基础图的过渡相。这是分层恒等式与描述证据，不是 Maxwell 独特因果机制或 GT 可达性证明。

## 6. 相阈值条件矩：理论正确，当前算法不胜基线

真实复属性与相事件配对监督，未独立化 E(delta) 与事件概率。先运行 K8/K16 点阈值，再把标签改为精确 h 质量 bins，逐样本求和恰好等于原 decoder；因此去掉了求积偏差。最后只改变为 anchor-induced mass 初始化，保留标签、容量和优化器。

三种版本均未超过强普通 CNN。第一批精确 bin8/bin16 的 MSE 为 0.007742/0.007746；anchor-bin 更差，约 0.010263/0.010019。宿主 plateau 的 mass/logit 梯度存在已验证抑制，但不能由此解释所有纯相及过渡相。停止继续调整 epsilon 或扩大此 head，不把恒等式当成算法成功。

## 7. 完整条件概率：相对改善是实测事实

Spline 更新全部 66 维，显式逆变换及 Jacobian 自动核验；三角模型真正为 q(a|y)q(b|a,y)，不补独立的剩余 prior。显式密度模型以原 joint simulation proper NLL 训练；CFM 以同一联合样本的全66维 velocity 回归训练。

|独立概率指标，三个 seed 的 parent loss 均值|新512|新256|
|---|---:|---:|
|同容量 diagonal Gaussian NLL，nat/parent|83.07743|83.44856|
|Full-66 Spline NLL|**77.50922**|**77.27755**|
|三角随机 frame NLL|77.95014|77.86049|
|三角 Maxwell frame NLL|**76.59216**|**76.25860**|
|anchor CFM 随机同谱控制 material energy score|0.165310|0.174586|
|有界 Maxwell 几何 CFM energy score|**0.130661**|**0.138822**|

Spline 对同容量 Gaussian 改善 5.568/6.171 nat；固定物理 frame 对随机 frame 改善 1.358/1.602 nat。相同原测度下，joint-test NLL 差估计两个条件模型相对真实 posterior 的期望 forward-KL 差，未知的 conditional entropy 相消；不需要 MCMC teacher。它不提供绝对 KL、全局覆盖或可信区间资格。

CFM 物理权重的 energy score 改善约 20.96%/20.49%。权重在实际网络标准化坐标里有界于 I 和 4I，所有 66 个速度分量保留；相同谱随机旋转控制享有相同容量、数据和 base。这个正面结果不等于 density 已有可用精确 logq。

值得注意：三角 Maxwell frame 虽改善 latent NLL，却没有同样改善 material energy score；第一批随机 frame 0.09964，Maxwell frame 0.10127。不同 proper score 衡量的目标不同，不能挑一个指标把所有结果叙述成一致胜利。

## 8. 为什么仍然不能称完整 posterior 合格

每个 parent 使用预冻结 seed1 的 index0 样本，不能选择最小 residual。原 posterior 的 prior-predictive 检验应给 Phi~Gamma(256,1)，均值 256。实际 Spline 的两批均值为 8,371.95 和 8,662.42，第一批 510/512、第二批 255/256 的 PIT 超过 0.975；几何 CFM 分别为 22,139.70 和 23,557.96，仍严重失配。

对固定前16对象的 K64 Spline 重要性修正，ESS 中位数为 1、最大仅 2.273。MCMC 校正核的正反比正确，拒绝重复状态完整保存，但有限链的 ESS/MCSE 不合格，没有据此输出 FULL posterior 均值或 calibrated mass。

实际训练了三类、成对 physics/supervised 的混合双向 KL 方案。Full Maxwell 力可以降低生成样本 Phi，但 NLL 恶化超过预先登记的复制门，覆盖明显受损；没有继续扩大种子。失败说明当前损失与模型在 response-fit/coverage 间存在冲突，不说明 posterior mean 不如 MAP。

CFM 的数值密度另有独立问题：双精度全迹 DOP853 审核用尽 895 秒上限，两个 tolerance 的某对象 NLL 差 0.136 nat，未过 0.05 稳定门；NFE32 离散值与更精细结果也有 0.342 nat 差。其材料样本 score 仍可评价，但不能把连续 CNF 近似密度当有限步精确 map density 用于 exact-MH。未完成的八个求解标记 NOT_RUN，而不是擅自推断全家族失败。

## 9. 风险图与合法物理世界不是同一个输出

对全部 768 个对象的六种风险图，进行了 4,608 次实际 Full forward。新主要 seed1 clipped 图只有 242/512 和 125/256 满足原 noise-ball；冻结 R2 为 256/512 和 129/256。风险改善没有自动带来更多数据相容图。逐像素 clipping 也不能证明来自原 66 维 decoder。

另对固定第一批前64对象、五个臂进行独立合法材料代表优化，实付 80,640 RHS。所有 320 臂都在预算处停止，未达到普遍驻点；84 个回到原合法 anchor。输出满足数据门的原因包括显式 fallback，不能把它冒充优化全部成功。

原普通 SLSQP 的 callback 只观察大迭代试探点，并非保证线搜索接受点，导致最优已付费合法查询未被选中。原失败保留，零新物理 best-query 修复和新的 Full-J trust-region 诊断在首次 GT 开启前完成。后者证明：data-only TR 可以继续降低全部64对象 Phi，却把 MSE 从约 0.010514 升到 0.016346；保持原 Gaussian prior 的 TR 则几乎不改变风险，约 0.010505。现有 NN 收益不能仅归因于一个严重未优化的 regularized baseline，但也没有全局最优性结论。

材料代表投影保留部分风险收益，Spline raw 约 0.006748 投影后约 0.009154；代表比原图更差。这不是“理想合法 Bayes action 必然更差”的证明，因为投影未充分收敛、表示也受 decoder 限制。

## 10. Maxwell 复用与接受决定：什么真的节省了成本

精确 active-row T-matrix 用正确本构增量、缓存 inverse action 和全照明场。与 Full 数值一致，但有效变化通常涉及 180–256 个 cells，低 latent 维不等于低端口秩。含 setup，Full 7.427 秒、复用 8.374 秒，停止加速主张。

数值势区间的首轮链同样失败过：虽然 RHS 下降，wall 从 229.78 秒升到 355.79 秒。随后只收紧可用逆范数上界并减少不必要动作，未改 prior、kernel、接受随机数和统计目标。前瞻新对象的安静重放得到：

|同一8对象×1024步|完整计费 RHS|全流程秒|
|---|---:|---:|
|optimized Full|65,536|71.397|
|区间 Richardson|4,608|24.153|
|同缓存全局 GMRES，决定驱动|14,691|37.605|
|同缓存 GMRES，固定精度|20,759|46.832|

所有逐步接受/拒绝、潜在区间检查及状态路径零差异。包含每对象 576 列缓存 setup、rebase、残差/低阶线性动作与 fallback 权利；减少 RHS 不是单独可信指标，实际 wall 也下降约 66.17%。结论范围是此 beta=0.002 核的数值重放，不是 posterior mixing 或材料 MCSE 的加速证明。更大 beta 的边界实验显示 fallback 增多，接受和实际移动迅速下降。

## 11. Learned Jacobian 与非完整探索

六次局部 value/JVP surrogate 训练全部未过最大白化输出误差和 Phi 误差 0.05 的门。已有 teacher 成本、普通 field-update、Broyden、frozen-J 和 Full-J 均公平保留。学习动作误差的排序器仍需要候选 Full 场特征，且弱于便宜的 total-field+exact-decoder 更新。它不是完整 LazyDINO 的复现，也不能据此贬低该文献。

四段 Euler 的二阶项混入自身场加速度，不能直接当 Lie bracket。修订为 Heun 后，原任务尺度尝试 288 项中仅 9 项满足2%材料距离匹配，最终全部支持门只有2项；直接方向法有94项，frozen/dynamic 共同阈值描述为59/83。原1% teacher 条件未改写，新的2%门只是离线比较适用域。未匹配平均 GT 误差不能用于算法排名。

共同 archive 按不读取 GT 的 Phi/UID 排序、完整材料 RMS≥0.02 去重，八份观测分别有 15、12、7、16、28、16、15、10 个候选（含基础图），即 111 个新候选。它们满足各自 noise-ball、ΔPhi≤2、prior_delta≤2 和进度门，不代表两个世界的预测距离≤1，不是 exact fiber，也不是 posterior mode 数。

## 12. 成本、结论与下一步

独立资源账本和全部 CLOSED 收据是费用权威。最终资源与失败明细见 [RESOURCE_SUMMARY](RESOURCE_SUMMARY.md)；复用历史训练数据没有把它们的冷生成费写成零。全部原始材料、latent、current/field、预测、模型、chains、失败和配置保留在证据包中。

本轮没有把全部研究押在 MCMC 或某个 Flow 上，三条路线都完成了真实实验和有实质变化的修订。没有启动三维、新 Mie、外部付费计算或未经门禁的48/64网格及分布外扩展。

**科研裁决：继续风险成像与静态物理条件学习；保留限定场景的决定驱动数值求值；停止扩大当前相阈值 head、active-row 加速、动态四段路径和当前 learned-J surrogate。完整 posterior 与强 IEEE TAP 核心创新仍未建立。**

下一轮最重要的不是再换一个网络名称，而是把已有相对概率进步转成同噪声尺度的 response calibration，同时保住覆盖，并检验在更大二维系统上相同决定精度是否仍有端到端成本优势。所有新目标与预算需另行冻结；本轮不再追加计算寻找正例。
''')
    md('MATERIAL_RISK_RESULTS.md', '''
# 材料风险：已执行的估计器与独立评价

主表和区间见 [完整报告](RESEARCH_REPORT_ZH.md) 第4节及 [全部45方法原始表](evidence_summary/ALL45_METHOD_ACTIONS.csv)。主要动作是逐组件 clipping 后图像；raw 同时报告，不用 clipping 证明 decoder 合法。

## 实际模型与公平归因

R2 冻结 baseline 为5,842参数局部 CNN。新首轮直接 CNN 与它的函数类别相同，只需将 R2 residual_scale 吸收到最后一层权重和 bias；数值复播误差约1.21e−7。因此重新训练的收益记为强基线改善，不能叫新架构。

strong_spectral 使用1,024个付费 TRAIN physics contexts、256个 validation；decoder/random控制共享数据、容量和训练预算。静态谱相对随机空间包络的两批改善为1.20%/0.83%，不是动态解流形贡献。

最终主要 fusion head 有三层16通道3×3卷积、GELU及2通道head；输入8个原重建/物理通道，加4个全局材料均值实/虚与方差实/虚通道，共6,418参数。TRAIN-only normalizer，batch64、lr5e−4、weight decay1e−3、4,000 updates，三个固定seed。原基础图加残差输出。另需冻结provider约0.68M参数、K32及原物理反演费用。

provider 比较包括显式 Spline、同容量 diagonal Gaussian、直接全场矩网络，以及同图同值但固定按 parent 打乱空间对应的控制。主要模型按原R2VAL256选择，不按新DEV或新GT重选。两批均保留此选择。

## 独立主结果

Spline fusion 的 MSE 为0.00680623/0.00737782；R2为0.00783202/0.00832320；strong CNN为0.00719170/0.00775824。Gaussian fusion为0.00682170/0.00739262。Spline与Gaussian没有显著风险差，主要成像增益不能称为Flow独有贡献。

直接生成均值很差，但其全局信息能帮助局部模型。这个现象通过空间打乱控制与普通provider控制检验，属于有限信息表示与监督任务风险，不要求输出已校准的posterior。

## 风险与物理输出分离

新主要seed1 clipped图仅242/512、125/256满足noise-ball；它们是风险动作，尚未证明在原decoder image。合法代表在固定64对象上保留部分风险收益，却损失大部分原图精度；所有优化预算停止和fallback都保留。

不能用GT挑生成样本或候选当在线结果；本轮没有把Best-of-K作为主指标。统计单位为parent，三个seed是训练重复而非新增样本。所有原始每parent loss及canonical batch预测已由独立进程复核，loss最大差5.55e−17，bootstrap最大差2.15e−16，失败0。
''')
    md('PHASE_MOMENT_RESULTS.md', '''
# 相阈值条件矩：已完成三种机制比较

原饱和函数满足 h(s)=∫h'(t)1{s>t}dt。必须回归复属性与事件的联合矩，不能用E(delta|y)P(s>t|y)代替真实条件相关。

实际先运行K8/K16点阈值，再采用 h 质量等分的精确bin标签 L_j=delta·clip(Kh(s)−j,0,1)，于是每个原TRAIN材料都有 ΣL_j/K=delta h(s)。这移除了数值积分偏差，未改prior或decoder。最后将均匀mass初始化改为已付费anchor-induced joint masses，eps固定1e−4；其他条件不变。

第一批512新对象的 MSE：点阈值8/16为0.00768661/0.00764516；精确bin为0.00774152/0.00774645；anchor-bin为0.01026344/0.01001919。它们没有超过strong普通CNN0.00719170，更没有超过fusion0.00680623。完整两批表见evidence_summary。

宿主平台处的真实moment-loss logit梯度相对均匀初始化约1.81e−6/2.67e−6，说明某种初始化确有优化障碍；纯夹杂端点通常仍有两个显著mass，不能推断所有饱和区都出现同样epsilon病态。

结论：保留正确恒等式、精确标签与条件误差凸性界，停止扩大这个特定单调mass head。约99.9%的fusion增益位于基础图的过渡相，但这不证明必须使用相阈值网络，也不证明有限模型能逼近真实条件矩。未来若研究此结构，应解决共享低容量回归与优化问题，而非继续减少求积误差。
''')
    md('CONDITIONAL_FLOW_RESULTS.md', '''
# 完整66维条件概率：模型、正面结果与资格边界

## 实际实现

Full-66 rational-quadratic Spline（8层、约0.683M）、同容量diagonal Gaussian（约0.683M）、较小Gaussian、固定16/50三角分解（约0.549M）、raw-y CFM及匹配容量CFM均真实训练。保留三seed、4k/16k/32k数据规模、learning-rate与容量检查。模型在原Gaussianized latent上有完整支持；uniform属性由probit变换正确处理。

三角法为 q(a|y)q(b|a,y)。固定TRAIN Maxwell frame 与随机正交frame享有相同数据、训练容量和概率体积。另运行真实付费anchor-conditioned none/GN/random坐标对照，满秩GN矩阵的logdet完整保留；GN预条件在当前拟合下恶化NLL，未强制作为主法。

有界几何CFM保持所有66维，以网络坐标 J_r=J_x T 构造 W=I+3Vdiag(s²/(1+s²))Vᵀ。W仅依赖y和其在线context且正定，不硬投影；随机控制保留同一谱。三seed×geometry/random实际六次拟合，固定3,000updates。原始坐标的度量是T^{-T}WT^{-1}，不能误报原Euclidean界。

## 两批独立评分

Spline NLL77.509/77.278，matchedGaussian83.077/83.449；物理三角NLL76.592/76.259，随机三角77.950/77.860。Spline相对差5.568/6.171nat，物理三角相对差1.358/1.602nat，parent区间均支持改善。

几何CFM energy score0.130661/0.138822，对随机控制0.165310/0.174586改善20.96%/20.49%。它仍远弱于好的材料风险动作；第一批模型样本均值MSE约0.03568。显式Spline样本均值约0.02278，不能用最接近GT的样本替代均值或在线决策。

这些proper scores可提供真实相对conditional-distribution证据，不能被“没有MCMC reference”抹掉。它们又不能给出绝对π距离或概率认证。[标准评分解释](AUDIT/JOINT_SIMULATION_PROPER_SCORE_INTERPRETATION.md)逐项说明。

## 实际替代机制与失败

三类physics-trained mixed bidirectionalKL与成对supervised控制真实运行。physics目标改善平均response却显著损害joint NLL和覆盖，未过复制gate，因此停止扩大。所用source-adjoint8与forward8、validation Full查询及失败均记费；不是只写理论式。

CFM完整散度审核采用CPU双精度DOP853、全66维trace、两预定tolerance，895秒预算后保留部分完成结果。某对象NLL0.136nat的容差差超过0.05；离散32步与更精细积分也不一致。CFM数值logq不用于exact-MH，未完成求解为NOT_RUN。样本energy score仍是合法的模型样本评分。

真正FULL资格见[概率审计](PROBABILITY_CALIBRATION_AUDIT.md)：本轮未通过。相对评分进步与绝对失配同时成立。
''')
    md('PROBABILITY_CALIBRATION_AUDIT.md', '''
# 概率审计：没有用风险改善代替posterior正确性

完整目标是原P1 latent posterior，不是局部切片或profile。模型明确区分完整可归一化Spline/高斯/三角密度与CFM数值生成law。对风险图、合法代表、几何candidate，不输出伪posterior权重。

## 固定生成样本的必要物理检验

若X'是真实posterior draw，Y来自原prior predictive，则Phi(X';Y)应服从Gamma(256,1)。每parent固定index0，不用best residual筛选。实际两批分别为：

|模型|512平均Phi|256平均Phi|结论|
|---|---:|---:|---|
|raw Spline seed1|8,371.95|8,662.42|严重失败|
|anchor普通CFM|36,282.18|40,924.89|严重失败|
|anchor几何CFM|22,139.70|23,557.96|比随机改善，绝对仍失败|
|anchor随机控制CFM|47,633.59|49,493.64|严重失败|

原期望为256。物理PPC共7,104个实际状态、56,832RHS、无失败，且未根据GT或response挑样本。主Spline的PIT超过.975为510/512、255/256。生成law没有经过FULL资格。

边际coverage偶尔接近标称不抵消此失败；plateau原子与ties使用随机化rank单列，不能将全图含大量宿主单元的覆盖率当ROI或完整joint校准。全部latent/material/projection/phase分层指标见evidence_summary/DENSITY_*.json。

K64重要性采样仅固定前16parent，ESS中位1、最大2.273。它是权重集中度，不是独立链ESS，不构成posterior teacher。短校正链保留所有拒绝重复状态，形式MH比正确，但material函数ESS/MCSE仍UNQUALIFIED。

## 仍然成立的相对概率证据

在同一joint law、原latent参考测度和独立parent上，显式模型proper NLL差等于未知posterior的期望forward-KL差；未知conditional entropy相消。由此可判Spline比同容量Gaussian更接近，不能说已知绝对KL或完全覆盖。CFM只用样本energy score，不给不稳定数值密度加资格。

数据依赖rank在固定前64有实际truth prediction的对象上进行；没有缓存truth prediction的其余对象标NOT_RUN，不新增物理补齐后伪称全768 rank审计。Gamma固定draw本身覆盖全部768。

## 输出身份

风险估计器：监督material-risk action；生成样本：qθ的材料推送；约束代表：candidate-only合法世界；几何端点：可达性见证。四者分开保存。posterior probability/mode mass默认null，禁止residual softmax。
''')
    md('MAXWELL_REUSE_RESULTS.md', '''
# Maxwell有限反馈、数值势区间与实际费用

## 有限状态更新的正确作用

所有候选使用同一材料与8照明。本构增量为实际ΔX，不是错误复相位或任意current。精确关系 Δj=(I−R0ΔXG)^{-1}R0ΔXE0 与原Full端点核验。Jacobian变化分为传播反馈、内部总场注入和decoder导数三个有序望远镜项；单独冻结后两项不正确。高阶恒等式本身不是多解存在性或加速证明。

active-row exact T-matrix 支持集通常180–256/256，含setup复用8.374秒、Full7.427秒，净成本失败。原记录保存reduced RHS、缓存内存与fallback，不将低latent维数误写为低端口秩。

## 决定驱动误差区间的真实迭代

首轮完整链与Full相同决定却更慢（355.79比229.78秒），已封存。修订只收紧可验证的逆范数包络，使用min(Frobenius,sqrt(1范数×∞范数))，并优化残差动作，不改变统计目标。

独立锁定8对象×1024步结果：Full65,536RHS/71.397秒；Richardson4,608/24.153；同缓存globalGMRES14,691/37.605；固定精度GMRES20,759/46.832。总105,594RHS全部计费；同一innovation、U、current-cache权利与落盘合同。所有接受、拒绝、状态路径及势包络重放零差异。

这支持约66.17% wall节省，范围是小步长beta=.002、当前TM24和缓存策略。普通Krylov也是强正对照，不能称A17独占。更大beta的边界实验：.005/.01/.02的3072次提案中接受302/95/6，实际移动快速变小，fallback增多。相同坏混合核的加速不等于得到合格posterior更快。

包络在普通float64及已验证离散后端条件下成立，没有定向舍入或连续PDE稳定常数证书。因此称“核验后的数值决定区间”，不称连续Maxwell certified posterior。

现有缓存取得需每parent576列inverse actions；setup/rebase/低阶线性求解/残差/fallback全部保留。未观测到稀疏active-set净优势，不隐藏failed runs。详见evidence_summary/DECISION_CHAIN_COMPARISON.json及原closed收据。
''')
    md('LEARNED_JACOBIAN_RESULTS.md', '''
# Learned-Jacobian与物理动作：真实训练后的裁决

共享几何teacher取得459,712RHS，含10,063个独特配对；公共physics anchors另计。不是每seed重新付同一teacher，复用不等于冷启动免费。

六次局部前向/JVP surrogate拟合（value/jvp×3seeds）都运行1,000updates。原最大白化L2≤.05与max|ΔPhi|≤.05门全部失败；实际最大白化误差约328–348，Phi误差约9.6e4–1.1e5。快速网络推理不能抵消必要Full回退。此处仅是小型局部DINO-like控制，不是完整LazyDINO文献复现。

三个动作误差排序器各6,000updates，误差约.00845–.00895，弱于普通total-field更新+exact decoder约.00724；frozen-J约.07474、Broyden约.02286。排序器还需要候选Full场特征8RHS，不是免费预测。

普通物理复用在材料action fidelity上有价值，但未证明提高最终风险、混合或可靠likelihood决定。停止扩大完整J回归，保留exact decoder、可核验残差与决定驱动Full refresh作为更简单基线。

完整算法、每seed选择update、训练与算子测试门和费用来源见AUDIT/PHYSICS_MECHANISM_FACTS.json与physics_mechanism_tables。原预发布summary可能写“等待locked”；最终已执行结果以本报告与results/IV/interval_chain_locked_v1为准，不改写中期封存文本。
''')
    md('NONHOLONOMIC_GEOMETRY_RESULTS.md', '''
# 动态弱几何：有限材料可达性与非完整路径的区别

真实prior material metric与Jacobian按rank2/4/8等原配置计算，包含decoder核、谱旋转、缝跨越和prior代价。训练、GT离线评分与candidate archive分开。

四段Euler不是纯commutator：二阶含 [X,Y]−DX·X−DY·Y。发现后保留原错误实现证据，改用Heun并做固定小步检查，不把错误Euler结果当理论反证。

原任务尺度尝试及离线2%距离适用门：Euler4/288匹配、全部支持0；Heun9/288匹配、全部支持2；直接方向288/288匹配、全部支持94（原首轮95）；frozen144/144、支持59；dynamic144/144、支持83。未完成和失败保留在原计划分母。原teacher1%门不改成2%；新表只回答哪些比较真正处于同一材料距离。

因此未匹配Heun/direct平均GT误差不能用来排名。直接方向更便宜且更容易保持距离；当前四段非完整路径没有成本或材料统计增量。

八份开发观测的共同candidate pool按在线Phi/UID排序、完整χRMS≥.02严格去重，最终每份观测有15/12/7/16/28/16/15/10个（含anchor），共111个新的prior-legal、有限噪声竞争材料。48个仅略低于.02的pair仍被拒绝，不用epsilon放宽。其支持含noise-ball、ΔPhi≤2、prior_delta≤2、signed progress及backward residual条件。

这些是不同观测下的candidate集合，不能合并称119个同一世界解，也不能称posterior modes/质量。它们未要求两世界白化预测距离≤1，未构造exact-fiber证书。weights/probabilities均null。

另外实际运行四种冻结满秩pCN核、8对象×1024steps、262,144RHS。Q正定且Q<I，原Gaussian可逆；完整MH保留拒绝重复。有限材料函数ESS约2.5–4.1，均UNQUALIFIED；bracket-enriched未建立独立混合收益。静态几何在监督评分上的正面结果不能用来追认动态路径成功。
''')
    md('LOCKED_GENERALIZATION.md', '''
# 两批真正独立的新对象：冻结与独立评分

第一批512：seed20261025020，prefixR4L1；第二批256：seed20261025021，prefixR4L2。二者都不是R2已暴露的256+256。128新DEV单独用于方法开发。

FROZEN_LOCKED、FROZEN_CONFIRMATION、FROZEN_DENSITY_BOTH_COHORTS以及risk fresh_freeze_v3在新观测产生前绑定。模型、checkpoint、CFM NFE32、K128、canonical CPU batch64风险推断、失败规则和主要动作均冻结。第一批结果未用于重新选择第二批模型。

所有768 truth forwards、noisy observations、physics-only native96+probit64 contexts完整保存。在线进程仅接收观测和合法context；GT分别放inputs/offline_only。两批45风险与26密度实例均在线完成、物理PPC和pretruth优化诊断封存后，才有FIRST_HELDOUT_GT_SCORE_RELEASE与独立评分PID。

独立审计复核全部45 raw/clipped损失、fallback、源绑定、chronology与4,000parent-bootstrap：失败0，source mismatch0，损失差≤5.55e−17，区间差≤2.15e−16。三个训练seed的loss逐parent平均，不当作2304个独立对象，也不当ensemble图。

Spline主要方法由历史validation的.00712673对Gaussian.00713016的近似平局确定。新DEV曾更偏向Gaussian，但没有据此换主要方法。Spline/Gaussian最终没有显著风险差；不得在两批测试后挑赢家重命名。

主要科学结果两批分别报告；合并768仅作描述。其他45/26方法消融不宣称全族错误率控制。现在两批GT已EXPOSED，未来需要新的独立parent，不把它们再次当locked。

首64物理代表/Full-J优化诊断为预GT封存的固定子集，与全部768主风险比较分开，不从其GT结果重调网络或stopping。
''')
    md('GO_NO_GO.md', '''
# 最终科研裁决

|问题|裁决|实际支持与边界|
|---|---|---|
|GO-A 材料风险成像|**GO**|512/256新对象相对R2改善13.10%/11.36%，相对强CNN改善5.36%/4.90%；普通Gaussian fusion一样有效，限定监督风险动作|
|GO-B 相对完整条件概率|**GO（相对）**|同容量Spline proper NLL改善5.568/6.171nat；静态Maxwell三角优于随机frame；几何CFM样本energy score改善约21%|
|GO-B 已校准FULL posterior|**NOT QUALIFIED**|Gamma残差严重失败、IS中位ESS1、校正短链不合格、CFM数值logq不稳定；不能给完整可信区间或mode mass|
|GO-C Maxwell计算|**GO，限定固定核**|quiet新8对象重放同决定同路径，wall降66.17%，含setup；非连续PDE认证，非混合增量；普通GMRES也有效|
|GO-D 静态几何学习|**CONDITIONAL GO**|NLL、energy score及小幅风险消融有重复改善，但不同任务指标并不全面一致|
|GO-D 动态Lie bracket与当前learned-J|**NO-GO 当前实现/规模升级**|距离保持差、路径昂贵、无材料统计增量；surrogate不过误差门，动作排序弱于便宜普通更新|
|GO-E 强IEEE TAP核心创新|**尚未建立**|有实质算法与独立实验，但通用CNN融合、Flow、LIS、resolvent和retrospective bounds均有成熟先例；没有新的Maxwell专有核心定理与跨尺度有效性|

停止扩大：当前相阈值mass head、active-row exact更新加速、四段路径、当前完整J回归，以及以反向KL强迫response-fit的当前参数。保留完整负结果，不改目标/噪声救gate。

下一步最有辨别力的正面实验：在同一P1/新独立对象上，构造既保留proper-NLL覆盖又把固定样本Phi降到原噪声尺度的可逆full66条件输运，并与同容量普通模型比较；或者在更大二维系统中，以相同统计误差/决定正确性证明包括setup与refresh的端到端收益。它们需要新的冻结合同，不能把本轮未用额度当永久授权。

风险成像已经值得继续；它应定位为物理基础图上的amortized material-risk estimation，而不是已解决的Bayesian posterior。A17在denseTM24没有独立净优势，删除其强制依赖不影响其他路线。
''')
    md('FAILURE_ANALYSIS.md', '''
# 失败、修复与未运行模块

本页保留科学失败与工程失败的区别。原CLOSED/ERROR、source snapshots、partial数组及失败账单均不改。

|对象|实际失败/局限|采取的不同机制|最终状态|
|---|---|---|---|
|点阈值材料矩|求积偏差、有限head风险差|精确h质量bin，再anchor mass初始化|三类均不胜强CNN，停止扩大|
|完整条件生成|相对NLL改善但样本Full残差很大|完整物理力+成对监督对照、三角和几何CFM|覆盖/response冲突未解决，FULL未合格|
|CFM显式密度|ODE容差和离散map不一致|双精度全trace adaptiveDOP853|895秒后partial，未完成项NOT_RUN；不用于exact-MH|
|active-row精确反馈|变化集大，setup不划算|与Full及普通goal/Krylov公平对照|加速NO-GO|
|首轮区间链|RHS少但wall更慢|收紧逆范数包络、同核重放|新locked8数值加速通过；原失败不改写|
|Euler commutator|二阶自身场加速度污染|Heun及直接bracket，距离门复核|Heun仍不能保持任务距离/成本，停止扩大|
|localDINO-like/J学习|误差门全部失败|value/JVP、动作ranker、普通fieldupdate/Broyden|不宣称文献复现或净加速|
|普通SLSQP代表|callback只记录major trial，未覆盖最佳有效query|零新物理best-query；普通Full-J TR双目标|原预算停事实保留，data-only更拟合却risk变差|
|fusion checkpoint元数据|闭账6份moment_K被参数数目覆盖|只读compat从封存TRAIN/VAL门恢复K32/0|权重/原元数据不改，canonical推断完整|
|fusion批次复播|2e−6技术容差原FAIL|另列1e−5只读技术核验；fresh统一batch64|原FAIL保留，不靠容差改写主评分|
|首次代表离线评分|Python list矩阵二维索引TypeError|附加wrapper仅转等价ndarray、新输出目录|新RHS0，原失败收据及partial保留|
|PPC/后拟合release接口|producer状态名/绑定键错误，物理前停止|独立source-only绑定修复|无伪造物理结果，真实失败wall仍计费|

未执行：连续PDE严格interval certificate、可信FULL posterior teacher/credible intervals、整个LazyDINO文献复现、完整动态全秩状态相关RL策略、48/64二维扩展、分布外noise/采集变化、3D和新Mie。任何尚未运行的原R3建议不标PASS。

没有因概率reference不合格停止监督或几何路线，也没有在测试GT公开后继续调参。现实限制来自高信噪比下完整law的response calibration、过渡相的非线性平均、局部谱的不稳定性、参数/field事件耦合和当前网络/优化机制；不把这些概括成“Maxwell无多解”或“posterior mean不如MAP”。
''')
    md('GPT_PRO_HANDOFF/START_HERE.md', '''
# 下一轮 GPT Pro：基于新实验证据，而非算法清单

先读[完整报告](../RESEARCH_REPORT_ZH.md)、[裁决](../GO_NO_GO.md)、[评分数学解释](../AUDIT/JOINT_SIMULATION_PROPER_SCORE_INTERPRETATION.md)与[机器指标](../RESULTS_AND_METRICS.json)。本轮真实完成多模型训练、768新对象独立评价及物理重放。

## 十项明确答复

1. Spline有稳定显式full66 proper-NLL收益；CFM有样本energy-score收益，但其数值log-density未过稳定门，二者不能按一个虚构统一NLL排名。
2. 完整posterior失败具体表现为原噪声尺度的response失配、模型center和dispersion共同贡献、重要性权重退化及物理训练后的coverage恶化；不是只报ESS小。
3. 相阈值恒等式与精确bin正确，当前head不胜CNN；anchor mass初始化还出现平台梯度抑制。
4. 直接生成材料均值不胜CNN；作为全局context却有独立风险增益。Gaussian一样有效，不能归因于Flow概率正确。
5. 静态Maxwell frame/权重对一些proper scores与小幅风险有重复增益；动态weak/bracket没有匹配成本优势。
6. 当前J/动作学习弱于便宜field+exactdecoder更新，且有昂贵teacher/候选Full费用，不值得直接扩大。
7. active-row finite feedback没有净优势；残差/势区间驱动数值重放有66%wall收益，范围限定且普通Krylov也有效。
8. 非完整路径没有材料统计增量；直接法更容易达到同RMS。111有限噪声候选不等于modes或posterior mass。
9. 本轮没有新Maxwell定理已被证明。标准conditional-risk、proper-score、resolvent和retrospective判定必须正确引用，不能包装原创性。
10. 最值得继续的问题在下方，必须围绕现象提出新机制。

## 三个具体数学问题

**一：response calibration与forward-KL覆盖为何冲突？** 同容量Spline已显著改进joint-test NLL，固定draw Phi仍约8千；加入真实Maxwell反向KL降Phi却使NLL增加逾10nat。能否根据decoder饱和image、response共维数、条件厚度与有限容量，构造不坍缩的可逆概率输运或可量化的response-tube factorization？需要给出比再次加大likelihood权重更具体的条件。

**二：局部patch缺失的任务信息是什么？** 全局Spline均值很差，却帮助6,418参数局部head；Gaussian几乎相同，空间打乱失去收益，约99.9%增益落在基础图过渡相。能否把此现象写成可计算的transition-cell conditional-risk decomposition，区分真正Maxwell sensitivity信息与普通空间/属性context？标准条件期望正交性只是起点，不是新定理。

**三：何时只为概率决定付物理精度有利？** 相同Full核上已实现数值决定零差异和66%wall节省，大步长则fallback/拒绝上升。能否把preconditioner收缩、current injection、acceptance margin与material-function mixing一起优化，而非只优化同一坏混合核的求值？需包括缓存setup、有效统计移动、浮点包络和更大二维系统成本。

## 不应重复的路线

不再把Gaussian强补全当FULL posterior；不拿纯trajectory频率当概率；不重复Euler四段错误；不把有限Neumann恒等式当多解定理；不将teacher/anchor成本视为免费；不使用两批已暴露GT做新blind claim。

本轮GO-A成立，GO-B仅相对，GO-C仅限定数值核，GO-D动态部分未成立，GO-E尚缺关键机制。新理论可以推翻原假设，但必须解释实际数组与负结果。
''')
    # Keep the initial theory audit verbatim and append final status in a separate document.
    old = S/'THEORY_AND_ALGORITHMS.md'
    if not (S/'AUDIT/THEORY_AND_ALGORITHMS_EXECUTION_SNAPSHOT.md').exists():
        shutil.copy2(old, S/'AUDIT/THEORY_AND_ALGORITHMS_EXECUTION_SNAPSHOT.md')
    initial = (S/'AUDIT/THEORY_AND_ALGORITHMS_EXECUTION_SNAPSHOT.md').read_text()
    md('THEORY_AND_ALGORITHMS.md', initial.replace('状态：研究执行中。本文记录已经实现或验证的机制，不替代最终独立测试及科研裁决。',
       '状态：实际实验与两批独立评价完成。初始数学审计保留于AUDIT；下列机制须结合最终报告的正负结果与限定范围。') + '''

## 8. 独立joint simulation可验证的概率和风险结论

平方材料损失下 R(T)−R(mu)=E||T−mu||²；对相同latent参考测度，NLL(Q1)−NLL(Q2)=E_Y[KL(pi||Q1)−KL(pi||Q2)]。样本energy score亦为相应材料energy-distance的proper比较，允许材料分布退化。它们都是标准身份；独立parent测试可估计相对差，却不提供未知绝对后验误差。

这个区别在本轮具有实际意义：有相对proper-score改善，同时有严重Full-PPC失败；二者并不矛盾。CFM的数值logq不稳定则不能调用显式density的NLL结论。详见AUDIT/JOINT_SIMULATION_PROPER_SCORE_INTERPRETATION.md。

## 9. 最终已完成与未成立部分

完整66维模型、三角条件、bounded geometry CFM、精确相bin与anchor-bin、mixed bidirectional KL、source-adjoint物理训练、同Full决定区间、direct/bracket/Heun和独立测试均已真实运行。相阈值head、active-row加速、当前learned-J与四段路径未建立独立优势。有限噪声候选确实存在，但没有exact fiber或mode-mass证明。

风险融合的主要解释是局部材料patch与全局观测统计的互补，不是生成模型已经校准；Gaussian控制与空间打乱控制均已执行。静态Maxwell几何有特定proper-score增量，动态几何没有由此自动继承成功。
''')
    md('REPRODUCE.md', '''
# 复现与证据读取

## 读取完成结果

先读RESULTS_AND_METRICS.json和evidence_summary。完整Release包含模型、联合训练数组、原始posterior生成样本、物理端点/current/field、J/谱、chain拒绝重复、source snapshots和追加式账本。release分卷按EVIDENCE_PACKAGE说明恢复；不做SHA检查。记录文件size/mtime_ns用于溯源，不能把它们当内容密码学证明。

代码依赖Python3、NumPy、SciPy、PyTorch、matplotlib以及原A30后端；实际版本与CPU/MPS资源见Stage0源清单。原运行解释器为../../.venv_nn/bin/python；所有数值库单线程。已有本地MPS训练，未使用CUDA、云或新3D。

## 可执行入口

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONPATH=src
python -m p3r4.train_density --help
python -m p3r4.train_anchor_density --help
python -m p3r4.train_geometry_cfm --help
python -m p3r4.train_risk_fusion --help
python tools/prepare_r4_fresh_risk.py --help
python tools/infer_r4_fresh_density.py --help
python tools/run_r4_interval_chain_locked.py --help
```

精确as-run argv、环境、source snapshot和预算收据在各stage/jobs/<job>/START.json、COMMAND.json及CLOSED.json；模型RESULTS/TRAINING日志包含训练seed、容量、选定update及normalizer。不得清空原Book重新运行；复现实验应建立新的独立输出及授权，或只读取已完成数组。

原封存manifest含当时绝对路径。下载证据后先按packaging清单把原根映射至新根，保留原封存文件，使用新副本/明确路径映射执行replay；不改原snapshot身份。源文件的新运行与旧mtime身份核验应分开，不用touch掩盖版本差异。复现关闭账本的正确含义是读取证据和生成相同算法的新结果，不是把新机器伪装成原进程。

输入身份：原joint x/q是TRAIN监督标签；anchor_conditions/context_x/context_J才是观测驱动context，二者严禁互换。新heldout GT只在offline_only目录，在线过程不导入其标签。模型元数据K修复见AUDIT/fusion/providers/CLOSED_K_METADATA_CLARIFICATION.json，只读compat保留原checkpoint字节。

## 报告再生

tools/build_r4_final_reports.py从closed评分和账本产生紧凑指标及中文报告，不重新求解、训练或选择模型。tools/plot_r4_final_risk.py及最终各figures目录的code/INDEX提供画图来源；图片均来自实际数组，不用生成图代替实验。

原P1/A30/R2依赖与实际Git revision列入SOURCE_INVENTORY和证据包。发布不额外运行物理smoke campaign，读回网页或Release不当作科学gate。
''')
    md('RESOURCE_SUMMARY.md', '# R4 资源快照与费用语义\n\n' +
       '|阶段|complex RHS|state attempts|LU/factorization|受控作业wall秒|物理动作wall秒|MPS秒|\n|---|---:|---:|---:|---:|---:|---:|\n' +
       '\n'.join(f"|{k}|{v['RHS']:,}|{v['states']:,}|{v['LU']:,}|{v['controlled_wall_seconds']:.2f}|{v['physics_wall_seconds']:.2f}|{v['GPU_seconds']:.2f}|" for k,v in ledger['stages'].items()) +
       f"\n\n当前闭账总RHS **{totals['RHS']:,}**；受控作业wall之和{totals['controlled_wall_seconds']/3600:.3f}h；其中物理动作实测{totals['physics_wall_seconds']/3600:.3f}h；MPS {totals['GPU_seconds']/3600:.3f}h。并行作业wall求和，不用日历间隔代替；训练wall和MPS时间并列，不相加当额外耗时。最后打包/绘图/审计仍各有analysis收据，最终以RESOURCE_LEDGER.jsonl及封账快照为准。\n\n" +
       '硬上限58M RHS /94h受控wall/96h GPU；stage分别独立约束，不借历史余额。所有forward、receiver/source adjoint、LU、reduced小系统RHS、失败、拒绝、Full fallback、cache setup、数据与模型训练计入对应作业。离线读取与打包不产生新RHS，但实际wall另记。\n\n' +
       '历史32,768+2,048 joint pairs与旧R2模型是复用证据，不声称从零部署免费；新增contexts、768新observations/anchors和所有物理控制计入本轮。多seed共享同一teacher的拥有费只记一次。小系统RHS与full状态分列，不能仅看RHS节省宣称加速。\n\n' +
       '旧GPU行有序列化字段省略，但after.stages.GPU_seconds已真实累计；仅修复将来字段，不重复补记。clarification见AUDIT/GPU_LEDGER_SERIALIZATION_CLARIFICATION.md。')
    md('EVIDENCE_PACKAGE.md', '''
# 完整可复现证据包

本页由最终发布步骤追加真实Release地址与分卷清单；未出现真实地址前不称已发布。

紧凑结果与全部方法表在evidence_summary，模型/输入/每端点的原始大型数组通过现有B1仓库Release分卷提供。所有原始文件仍本地保存，不因发布删去失败或重复状态。

具体文件选择、大小、mtime_ns、复用数据身份、压缩方式和恢复命令见packaging。代码与报告通过原B1分支/main发布，复用现有Pages，不新建仓库。发布与匿名读回是交付核验，不是新增物理或统计资格。
''')
    md('PRIOR_ART_AND_CLAIM_LEDGER.md', '''
# 相关先例与本轮主张边界

这不是一次穷尽文献综述。此前实际ScholarQA检索与身份核验保留于AUDIT/literature_collect.json、literature_verify.json；429限制也保留，不能把检索空缺当原创性证明。下面使用原论文或作者版本。

|成熟工具/先例|原始来源|本轮不能据此单独声称的创新|
|---|---|---|
|Flow matching|[Lipman等，2022](https://arxiv.org/abs/2210.02747)|用CFM回归全维速度不是新Maxwell理论|
|Simulation-based conditional flow matching|[Dax等，2023](https://arxiv.org/abs/2305.17161)|用joint模拟标签训练条件posterior不是新概率合同|
|Simulator feedback posterior flow|[Holzschuh与Thuerey，2024](https://arxiv.org/abs/2410.22573)|加入物理/模拟器力不是本轮首次提出|
|电磁条件flow逆散射|[Zhang等，2022，IEEE TMTT](https://doi.org/10.1109/TMTT.2022.3205890)|电磁概率Flow本身已有直接先例|
|Amortized Gaussian/Flow variational inference|[2024原作者版本](https://arxiv.org/abs/2305.20004)|混合方向KL/监督与后验拟合需对照其已有机制|
|预条件normalizing flows|[2021原作者版本](https://arxiv.org/abs/2101.03709)|引入局部GN尺度不自动形成新算法贡献|
|DINO|[JMLR 2025](https://jmlr.org/papers/v26/24-0745.html)|本地小型value/JVP模型失败不等于复现并推翻DINO|
|LazyDINO|[JMLR 2026](https://jmlr.org/papers/v27/25-0858.html)|本轮不是完整LazyDINO复现，不能宣称成本胜过该文献|
|误差界驱动的retrospective接受决定|[2014作者机构全文](https://eprints.gla.ac.uk/96718/)|用概率阈值决定是否加精不是本轮新原则|
|Delayed acceptance效率分析|[原作者版本](https://arxiv.org/abs/1506.08155)、[adaptive DA](https://arxiv.org/abs/1509.00172)|普通近似或训练误差模型不能自动叫certified，setup也不能省略|

原P1/R2/A17-X/F3的具体本地来源与固定commit见SOURCE_INVENTORY.json。有限resolvent、Schur、Woodbury、条件期望正交性、proper-score身份和Lie bracket基础均作为工具，不以首次提出的形式叙述。

## 当前有证据的增量

1. 两批真正新parent上的material-risk改善，并能由强CNN、Gaussian/空间打乱控制分离归因。
2. 原66维测度下的相对proper-NLL与样本energy-score改善，同时保留Full-PPC失败，避免把风险与概率混同。
3. 同缓存、同完整核、同随机数上的Maxwell数值势区间费用收益，含setup和失败边界。
4. 对phase-bin参数化、物理训练coverage冲突、Euler路径污染和动态bracket距离失败给出可复核反例。

这些是实质实验现象。尚未证明它们足以形成新的Maxwell特有理论与算法，尤其缺少response-calibrated完整条件分布、与最强同成本方法的统计效率优势，以及较大二维系统的独立证据。最终强IEEE TAP结论仍为尚未建立，而不是把工程规模当新颖性。
''')
    with (S/'RESEARCH_REPORT_ZH.md').open('a') as f:
        f.write('\n## 13. 相关工作与可发表性\n\n'+
            '完整条件Flow、simulation-based learning、DINO、局部预条件与retrospective接受判定都有成熟先例。'+
            '本轮不把这些工具拼接后的名称当贡献；逐项原始来源和不支持的主张见[相关先例与主张账本](PRIOR_ART_AND_CLAIM_LEDGER.md)。'+
            '当前最有价值的是可被独立核验的新实验现象及失败机制，仍需要新的Maxwell专有机制与决定性比较。\n')
    md('ITERATION_LOG.md', '''
# R4 实质迭代与最终状态

## 风险路线

假设：完整属性–相事件矩会胜过直接回归。实际运行K8/16点阈值→精确h质量bin→anchor-induced bin初始化。依次改变求积误差与初始化机制，保持P1、标签/容量等共同条件。精确bin仍弱于强CNN，anchor-bin出现局部平台梯度抑制且风险更差，停止扩大该家族。

发现：首轮CNN与R2函数类别相同，收益先计为训练基线修复。随后扩大到1024context并保持decoder/random控制。另将全局model均值/方差作为四通道context，加入同容量Gaussian、直接field与空间打乱控制。新512/256独立确认主要风险改善13.10%/11.36%；Gaussian几乎相同，拒绝Flow特有风险归因。

## 完整概率路线

三seed Full66 Gaussian/Spline/三角/CFM以及数据规模、容量、lr真实比较。显式Spline和静态物理frame NLL有独立增量。满秩GN坐标在当前训练下失败，不强制使用。bounded geometry CFM改变训练误差权重而非硬投影，随机同谱控制、相同base与数据；两批energy score改善约21%。

进一步用真实Full source-adjoint做paired physics/supervised mixed-KL拟合：response改善但NLL/覆盖恶化，三族均未过复制门；停止seed扩张。CFM数值密度以双精度全trace自适应ODE独立核验，预算895秒后只得部分稳定性失败，其余NOT_RUN。不是重复Gaussianreference后结束其余路线。

## 物理计算路线

active-row精确复用因变化集合大而净wall失败，保留原结果。首轮完整区间链同Full决定却更慢；只修改逆范数包络与动作流程，在新DEV和预冻结8新对象quiet重放。新locked同核wall71.397→24.153秒、逐步决定/路径零差异；普通GMRES同时作为强控制。较大步长失效边界保留，未宣称posterior mixing改善。

六次localvalue/JVP和三次动作排序器真实训练，误差或相对普通fieldupdate控制不胜，停止扩大，不冒充完整LazyDINO复现。

## 动态几何路线

主线程发现四段Euler二阶自身加速度污染；保留旧实验，改Heun并做真实Full配对。小步核验不足以保证任务尺度，Heun288项仅9距离匹配、2全部支持；直接法94支持，停止四段路线升级。满秩冻结pCN四核实际运行，ESS仍UNQUALIFIED。共同archive严格RMS.02形成111新有限噪声候选，无后验质量声明。

## 优化诊断与工程修复

原SLSQP callback选点局限并非物理无解。零新物理best-query与新Full-J TR在GT开放前封存；λ0降Phi却增risk，λ1与原baseline相近。所有预算停和原失败保留，不事后改旧资格。

Closed fusion K元数据、PPC release状态与离线list矩阵评分故障使用独立只读compat或新输出恢复；原weights、partial、ERROR、预测和成本不覆盖。最终报告构建也保留首次schema assertion失败收据，仅修复对failure字典的读取，未改变科学数据。

本轮不再新训模型或付费物理。完整费用以RESOURCE_LEDGER和CLOSED收据为准，旧中期文本保留于AUDIT，最终主张由GO_NO_GO.md约束。
''')
    print(json.dumps(dict(status='FINAL_REPORTS_FROM_CLOSED_EVIDENCE',risk_parents=[512,256],resources=totals),ensure_ascii=False))


if __name__ == '__main__':
    main()
