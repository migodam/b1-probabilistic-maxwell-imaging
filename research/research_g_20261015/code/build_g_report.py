"""Write research delivery from completed arrays/receipts; no experiments."""
from pathlib import Path
import json
import shutil
import html

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = 'https://github.com/migodam/b1-probabilistic-maxwell-imaging'
WEB = 'https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_g_20261015/'
RAW = 'https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_g_20261015/RESEARCH_G_REPORT_ZH.md'
RELEASE = REPOSITORY + '/releases/tag/research-g-20261015'


def read(name):
    return json.loads((ROOT/name).read_text())


def write(name, text):
    (ROOT/name).write_text(text.strip() + '\n')


def table(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |',
                      '| ' + ' | '.join(['---']*len(headers)) + ' |'] +
                     ['| ' + ' | '.join(map(str, row)) + ' |' for row in rows])


def main():
    online = read('results/E1_LEGACY/RESULTS.json')
    offline = read('results/E1_OFFLINE/METRICS.json')
    discovery = read('results/E1_DISCOVERY_SUMMARY.json')
    witness = read('results/STRUCTURAL_WITNESS_AUDIT.json')
    resources = read('RESOURCE_SUMMARY.json')
    ledger = [json.loads(s) for s in (ROOT/'COST_LEDGER.jsonl').read_text().splitlines()]
    episodes = [e for m in online['methods'] for e in m['episodes']]
    qps = [r for r in ledger if r.get('event') == 'MATERIAL_QP']
    history = read('results/E0_HISTORY_SIGNAL.json')
    pooled = discovery['pooled_same_observation']['per_case']
    ptable = table(['同一观测', '材料代表（含base）', '结构代表（含base）', '新增结构'],
                   [[r['case_id'],r['material_base_plus_unique'],r['structure_base_plus_unique'],
                     r['structure_additional_unique']] for r in pooled])
    method_rows = []
    for name in ['multistart','static_scale','QD','bandit']:
        rows = [r for r in discovery['rows'] if r['round_id']=='first' and r['method']==name and r['prefix_t']==3]
        score = next(s for s in offline['summary'] if s['round_id']=='first' and s['method']==name and s['K']==4)
        method_rows.append([name,sum(r['attempted_states'] for r in rows),sum(r['complex_RHS'] for r in rows),
                            f"{sum(r['endpoint_valid'] for r in rows)}/9",
                            sum(r['new_material_events'] for r in rows),sum(r['new_structure_events'] for r in rows),
                            f"{score['mean_Best_complex_MAE_conditional_on_valid']:.8f}",
                            f"{score['all_available_case_historical_coverage']:.1%}"])
    mtable = table(['方法（首轮）','新states','复RHS','合法端点','材料新颖events','结构新颖events','Best-of-4复MAE','旧定义Coverage@4'],method_rows)
    valid = sum(e['valid'] for e in episodes)
    retained = sum(e['raw_structural_edit'] and e['retained_structural_edit'] for e in episodes)
    raw_edits = sum(e['raw_structural_edit'] for e in episodes)
    collapse = sum(e['correction_collapse_to_prefix'] for e in episodes)
    new_structural_events = sum(e['new_structural_world'] for e in episodes)

    first_answers = '''
1. **发现多少：**三个观测各有3、3、1个材料代表和1、2、1个结构代表（都含base）。新增结构只有一个；两轮重复同一world不能算两次。它是边界/支撑变化，不是新孔洞拓扑。
2. **单位成本：**bandit在首轮36states/864RHS中命中这一个结构，另外三方法为0；这是单对象描述性命中，没有稳定优越性证据。
3. **Best-of-4：**所有方法两轮平均复材料MAE均0.04224802，和base相同；没有最佳新方法。
4. **corrector保留：**唯一新结构从不相容raw提案修正为相容，投影retention约.138；大量其余提案失败或回到既有邻域，不能只数raw多样性。
5. **历史帮助：**径向32-parent分组验证中，增加历史没有稳定改善下一结构预测。失败/材料重复存在关联，不能推出更好的反事实动作。
6. **是否超越baseline：**没有建立持续、跨对象地超过continuous-scale static/multi-start；G实验冻结时未接入完整T包，T-guided未运行，因此完整五臂比较未完成。
7. **PPO/GRPO：**不值得现在训练；先积累共同corrector下的跨对象结构发现轨迹。
8. **进入非平凡3D：**尚不具备资格。保留真实二维witness与失败，不由一个近门槛边界变化升级为一般3D机制。
'''
    write('RESEARCH_G_REPORT_ZH.md', f'''# Research G：统一材料空间中的下一解搜索

2026-10-15 · candidate-only · 已执行E0与有界历史后端E1 · 不训练大模型

[中文网页]({WEB}) · [原始Markdown]({RAW}) · [实际材料图库](MULTISOLUTION_GALLERY/index.html) · [完整证据]({RELEASE})

**本轮科研判断为NO_GO_SEARCH：尚未建立值得扩大的持续新结构搜索优势。** 这不否定本轮找到的一个合法边界变化world，也不意味着完整Maxwell材料空间不存在更多解。主要缺口是跨对象重复发现与机制归因，而不是posterior未校准。

{first_answers}

## 1. 为什么不能把很多随机图当作多个材料世界

同一散射观测可能容许材料边界、损耗与介电分布改变。我们要保留的是不同的**物理材料**，且这些材料通过相同完整离散方程与观测规则。不同latent编码、轻微纹理或不同优化器状态都不能单独构成新世界；反过来，散射预测接近也不能把材料不同的候选删掉，因为这种歧义正是研究对象。

本轮完整阅读Pro包三份文档，保留其feasible-world、post-correction discovery、不可变success ledger与成本原则。直接用户预算384states/12288复RHS/10800秒取代Pro的1664-state建议，不叠加T的192state额度。G既不是posterior calibration项目，也没有必须使用NN的前提。

## 2. 与T的分工及本轮正式物理身份

E0时T尚无已提交physical handoff。因此，本轮使用用户明确允许的**完全相同、已经核验的历史冻结后端**：A30 commit `717d3d26a290248bd2b9aeb55f3d03694b407b84`，原adapter与Research A commit `f1b674eb4893e39126701ed48f3d3dcef3098c3c`中的GN函数，函数未改。不是另写T forward/corrector、nucleation或A17消融，也不是T新机制的试验。

物理是24×24标量TM/VIE、8源32接收器、complex128/float64、cell mass1/576，Reχ∈[-.5,4]、Imχ∈[0,.5]；逐source Re32后Im32，原观测和sqrt(2)/σ白化，gain固定1，相容阈值23.785449004057032。绝不混入B1径向Mie的gain/噪声约定。已保存Green/receiver矩阵与source-order descriptor是target身份的一部分。

历史八条件只有三个weak base合格，按原残差/材料合法性而非truth选中250069/73/77；新G各做一次付费full验证，预测重放、solver与c=χE诊断通过才继续。其余五条件保留E0失败，不用有利truth替换对象。三个对象已曝光，不是新的blind confirmation。

G冻结与数值收口时未导入T完整candidate/backend包，时间点记录见`results/T_HANDOFF_STATUS_AT_G_CLOSE.json`。发布整理时已看到T的公开提交`2102f6a`及candidate index，另记`results/T_HANDOFF_STATUS_AT_PUBLICATION.json`；这不追溯改变已经完成的实验。C=T-guided仍为NOT_RUN/本轮未接入，不能拿legacy结果冒称完成了T/G五臂共同corrector实验。没有向T chat发消息，也没有修改T工作树。

## 3. E0先解释旧成功与旧失败

旧径向1920个proposal的raw多样性很高，但multi-start有331/384条raw新颖提案修正后落回全部prefix合法材料的旧邻域；static为145/384。该邻域指标不证明相同动力学吸引盆或连通解分支。静态原本使用八个径向保护box，因此它的正面结果不能直接替代G的无孔数分类全field搜索。

旧二维有限编辑没有invalid→valid；保护格最大实际偏差0。唯一满足≥6格数字新洞定义的finite提案被保留却仍不相容；不能解释成“corrector把洞全填了”。C1训练proposal是受限四维径向材料，A30 Sparse64的35.69%是truth-fit容量，二者都没有提供已训练好的A30 full-field生成器。

654个合法duplicate曾改善留存代表full residual。我们把first-member改为pairwise-safe truth-free质量替换，在相同缓存池/成本上重放：369次替换、不丢原结构邻域覆盖，但没有新物理world。二维旧6格计数的3→4在新多尺度标准下仍3→3，不能挑有利门槛宣布新结构。完整[E0报告](E0_REPORT_ZH.md)、JSON与逐trajectory表保留。

## 4. 一套连续材料动作与四个搜索选择器

权威对象是576格复χ，而非固定孔洞专家。可用动作包含相关全域重启、连续位置/尺度/幅度的局部材料加减、多尺度分布改变。当前全域restart以合法base加独立Haar系数为起点；它不是全材料域均匀prior，也不是完整1152维支持。局部mask可以在不同位置创建或消除support，但有限带宽和corrector表示仍有限制。

static在连续尺度/有符号幅度上预定分层，不用0/1/2孔类别。QD用固定emitter混合及active archive的低使用/新颖parent；bandit用同动作库Beta(.5,.5)，只有**修正并full验证后相对全部prefix合法ledger的新结构**才获正reward。材料contrast新颖与结构reward分开。每三步第0步强制global，比例1/3，不能被筛掉；本轮无proxy screening，所有生成提案均full求值，screened-out false-negative audit不适用。

所有方法原样共享32实Haar GN两步、trust=.35、拒绝后×.3，最后最小full residual材料再做独立freshfull验证。既有fixed-step corrector可能继续移动已相容状态；它是历史regression corrector，不冒称T的proposal-retentive设计。

active archive为pairwise-separated材料/结构代表；success ledger保留全部合法端点含duplicates，frontier存失败而不安装。最低full residual选Top1；TopK保留合法base后做固定field/structure覆盖。少于K保持短集合，不填invalid。所有在线结果先冻结，独立进程才读已暴露truth。

## 5. 两轮同预算结果

每轮每方法三parent×三轨迹，每轨迹raw+两次corrector+独立验证，共36个新state。正常cache复用保留，没有为凑相同RHS额外求解；这次实际每方法都是864复RHS。共享三个base验证另计一次，不在每方法重复累加。

{mtable}

两轮只变first/quality archive，动作、随机种子、corrector、阈值保持配对。第二轮同一表，当前小pilot没有触发质量替换，Best-of-2/4/8与Coverage均未改善。E0修复在旧池中有用，并不保证新search trajectory受益。

bandit描述性结构效率为1/36 states（不含共享base）或1/864 RHS，其余三方法为0；只有一个对象命中，不能当作统计胜利。三步bandit轨迹只含极少正reward，没有足以识别学习价值的行为重叠/反事实实验。QD产生了额外材料contrast代表，但没有新结构事件。

## 6. 同一观测究竟保留多少世界

{ptable}

上表对同一观测跨四方法、跨两轮重新做固定顺序pairwise packing。两轮重复世界只算一次；不同观测的数字也不能相加成解释同一个观测的世界。18个合法端点ID加三个base共21个记录，不等于21个不同材料。

### 唯一结构witness的物理意义

`p250073`的bandit局部移除提案raw残差89.9456，两次原GN后25.7212→22.3180，独立full验证与重复预测一致，base残差20.8919也合法。solver/backward与共享本构缺陷约1.4e-15。它确实是从不相容提案修正得到的另一个合法field。

它相对base dM={witness['dM']:.8f}，低于材料分离.03；但dS={witness['dS']:.8f}、support改变{witness['primary_mask_XOR_cells']}格，通过预定.08/6格结构规则。support由42到52格，约是base support的24%增加，图像是**边界/支撑扩张**，两者数字孔数仍0。mask幅度门槛.20/.25/.30下改变13/10/8格；dS门槛加倍到.16则不再计结构不同。这是近工程门槛的有限witness，不作不连通分支、全局新拓扑或普遍重要world的证明。

![同一观测的base与新full-valid边界材料](figures/E1_structural_witness.png)

两轮精确重复这一world不是独立复现。它的Best-of-K未胜过base，不能归功于历史学习；但也不能仅因没有新孔洞类别，就把真实边界材料差异删掉。

## 7. correction retention与失败归因

72条新轨迹中raw full-valid为0，末端full-valid={valid}；{72-valid}条仍不相容，所有付费计算保留。{raw_edits}条raw满足对parent结构编辑规则，只有{retained}条在终点仍保留该结构差异（包含不相容终点）；{collapse}条合法终点从raw结构新颖回到prefix结构邻域。新结构events={new_structural_events}，跨paired重复去重后只有1。

144个材料QP中{sum(q['converged'] for q in qps)}个报告收敛，{sum(not q['converged'] for q in qps)}个达迭代上限/未收敛。未收敛flag保留，正式材料仍由独立full验证决定；不声称达到优化stationarity或全局最小。唯一结构witness的两步QP均报告收敛。不能把归因全部推给generator：有限corrector带宽、两步预算、初始化与实际物理约束仍混杂。

## 8. 离线成像与几何覆盖

全部方法、两轮的Best-of-2/4/8复MAE相同，均值.0422480213；三个parent分别.04351438/.03477679/.04845289。Top1也全部选base。所有Top4/Top8都短于请求K，不补图。

原RMS/1.2≤.1的Coverage@K为1/3，support IoU≥.8的几何coverage也是1/3。G新尺度距离≤.1给3/3，但那是不同归一化的**另一指标**，不能替换旧目标宣称全部成功。`p250077`真值有一个≥6格洞，所有合法候选仍无洞，hole recall为0；其它两个truth无洞时该recall为NA，不填1。软件phase-field的开孔演示没有改善这些真实Maxwell结果。

## 9. 历史信息与RL判断

E0按32个parent五折，结构B−A log-loss差{history['comparisons']['corrected_structural_novelty']['B_minus_A']['mean_log_loss_delta']:+.6f}、区间跨0；B对matched-history shuffle也无稳定优势。history-only D有预测关联，合法性/材料重复也能预测；这些可能反映动作分层和搜索耗竭阶段，不代表能选择更好提案。

没有logged propensity与跨动作counterfactual物理端点，不做off-policy价值判断。新二维轨迹只有三个对象，bandit的一次局部命中不能证明learning收益。暂缓PPO/GRPO、GRU、Diffusion和大型field generator；本轮没有提交其物理训练计划。

## 10. 理论与文献边界

Pro包的数学说明保留：生成器分布经corrector与validation后的推前才决定发现；raw多样性不能代替endpoint多样性。防御global探索的漏发现界依赖未知的正到达概率，不证明我们有限Haar/mask空间覆盖任意结构。连续参数可以改变support拓扑，但不等于某个训练decoder可达任意材料。

LITERATURE_MATRIX中已有topological sensitivity+GN、生成器后修正、unstructured QD/NSLC、surrogate illumination、PPO/GRPO先例。此处阅读了提供的矩阵，未把34篇源包阅读层级冒称本轮独立全文核验，也未重新做宽泛文献综述。普通archive、Beta bandit、材料约束GN与这些基本数学不是新算法定理；本轮未建立强IEEE TAP原创性。

## 11. 费用与下一步

新states {resources['new_full_state_attempts']}/384、复RHS {resources['new_complex_full_solver_RHS']}/12288，GPU0。受控数值wall {resources['controlled_wall_seconds']:.3f}/10800秒、CPU {resources['CPU_seconds']:.3f}秒，含E0、物理/优化/验证、失败评分尝试、统计与图库；没有将交互写代码、阅读或上传耗时伪装成数值执行。原A30/A/B/C1/P2/T账本不变。

先停止本轮付费运行，保留93states/5352RHS余额，不为挽救优势继续加提案。后续先只读接入T最终提交的physical handoff，核对target/packing/单位/共同corrector，再注册同预算五臂比较及跨对象独立复现。当前单witness只值得保存和机制复查，不够进入非平凡3D或RL训练。

**NO_GO_SEARCH**
''')
    shutil.copy2(ROOT/'SEARCH_HISTORY_SIGNAL_DRAFT.md',ROOT/'SEARCH_HISTORY_SIGNAL.md')
    write('QD_VS_BASELINES.md', f'''# QD、bandit与强基础搜索

{mtable}

表为首轮，第二轮相同。每方法9proposal/36state/864复RHS，公共3base/24RHS全局记一次。所有方法同合法支持/动作库/原32实Haar corrector/full validator；static分连续尺度与幅度，不分孔数。multistart为base-centered相关全域扰动，不冒称全材料域prior。

bandit描述性命中一个结构，QD没有；不是稳定超越static/multistart的证据。仅3个已曝光合格parent、每arm每parent3步，T-guided未运行，五臂不完整。三步beta更新不足以判断学习机制。

Best-of-2/4/8各方法全相同，Top1都是base；原Coverage@K=1/3，几何coverage=1/3。没有用G新归一化coverage=3/3替换旧目标。短K不填充。完整逐parent、prefix、实际RHS/wall见E1_DISCOVERY_SUMMARY和E1_OFFLINE CSV。

当前不能称QD/CEM/RL比强baseline更好，亦不能由小pilot证明它们永远无效。优先接入真正共同T handoff并扩大有效轨迹证据，未获新授权不启动下一campaign。
''')
    write('CORRECTOR_ATTRIBUTION.md', f'''# proposal → correction → full validation → archive

## 已执行与原E0的区别

原二维保护格缺陷0；31/32有限编辑未生出声明尺度数字洞，唯一新洞保留却不相容。原有限编辑invalid→valid=0。不能用“填洞”解释全部失败。

G重新生成72个真实field，raw full-valid0，endpoint full-valid{valid}，不相容{72-valid}；{raw_edits}个raw结构编辑、{retained}个终点保留（含invalid），{collapse}个合法终点折返prefix结构邻域。不要把保留而invalid的结构算成功，也不要把材料contrast reward算结构。

144QP：{sum(q['converged'] for q in qps)}收敛、{sum(not q['converged'] for q in qps)}未收敛flag。每个raw/两步/独立端点都保存χ、current、field与32×8复预测。QP未收敛不抹掉；端点合法性独立判断，未建立全局优化保证。

新结构world只有p250073 bandit一次：raw89.9456→25.7212→22.3180<23.78545，独立重算预测差0，solver和c=χE缺陷~1.4e-15；retention=.138013。support42→52，孔数0→0，是边界/支撑变化，不是新孔；dS=.081767，.16敏感性下不计数。两轮同world不算复现。Best-of-K仍由base取得。

新发现不能全归生成器，也不能说corrector毫无作用：这条提案是实际invalid→valid；其余大量提案由有限corrector折返或仍失败。区分物理不可行与32实chart/两步优化不足需要T的共同机制试验，本轮没有实现第二套corrector。

E0质量替换重放有369次替换，但G当前新轨迹中替换未触发；两轮结果相同。修复一个archive缺陷不自动带来新结构或MAE改善。
''')
    decision = '''# Research G 科学判断

本文件记录科学结论，不是额外审批流程。

**NO_GO_SEARCH：当前pilot未建立持续、跨对象的新结构搜索优势，不进入RL或非平凡3D。**

- E0：已执行，0新物理，历史归因和parent分组预测有证据。
- 统一field/动作/archive：已实现，材料与结构计数分开，正式candidate-only。
- 历史后端E1：已执行，3合格parent、四可用方法、first/quality两轮；不是T新机制的完整五臂比较。
- 真正新增：一个近预定门槛的full-valid边界材料world；仅p250073，paired重复只算一次，非新孔洞拓扑/已证连通分支。
- QD/bandit优越性：未建立；bandit单次命中只作描述性，Best-of-2/4/8和coverage均无改善。
- corrector：有真实invalid→valid witness，也有大量折返/不相容；不能由有限失败推材料空间唯一。
- 历史：有失败/重复的预测关联，未建立额外结构预测或反事实动作收益。PPO/GRPO暂缓，NOT_RUN。
- T-guided：NOT_RUN/本轮冻结时未接入；发布时T已有公开交接提交，后续需另作兼容接入。一般非同心3D、A17新消融、learned-field训练均NOT_RUN。
- 预算：291/384states、6936/12288复RHS；没有任何硬cap迫使提前截断，按冻结小pilot完成后主动停止，不追加寻找有利结果。

NO_GO不抹掉合法witness，也不以posterior校准或单一MAE作为项目成立条件。它只限制当前算法优越性与规模化主张。
'''
    write('GATE_DECISION.md',decision)
    write('E1_DECISION.md',decision)
    write('NEXT_STAGE_PLAN.md', '''# 下一阶段：先共同物理接口与结构发现复现

1. 只读导入T最终Git提交的完整physical handoff；版本化核对forward、单位/packing/噪声、合法base、corrector/validator和proposal callback。不得拿正在写的工作树输出自动开新实验。
2. 冻结同一corrector后补齐global multistart、连续尺度static、T-guided、QD与轻bandit五臂；对目前p250073边界witness作独立起点/跨对象复查，不能只复跑同seed。
3. 区分raw未产生结构、保留但不相容、修正折返、有限chart不可达与archive错误剔除；一次只变proposal或archive，不能把改corrector带来的收益归给search policy。
4. 保留原Best-of-K和几何coverage目标，短K如实报告。p250077真洞尚未覆盖，应成为困难暴露回归案例，不按truth指导在线提案。
5. 当前不提交PPO/GRPO训练方案。先有简单策略在多个对象上持续发现有意义新结构，并显示prefix历史可用于选择动作，再讨论学习成本摊销。

本轮E2大factorial、E3 RL、E4一般3D均NOT_RUN。Pro建议512/1944/4096等更大预算没有授权，本轮余额也不自动变为下一阶段预算。这里是未来方案，未执行。
''')
    shutil.copy2(ROOT/'NEXT_STAGE_PLAN.md',ROOT/'NEXT_STAGE_PROTOCOL.md')
    write('FAILURES.md', '''# 失败与限制

54/72个新endpoint不相容，仍保留raw、corrected、full预测及费用；这些是成像失败，不是未执行solver。QP24/144未达优化器收敛flag，完整QP状态在COST_LEDGER，另索引FAILURE_LEDGER；full候选不能凭优化器flag安装。

E0初次无阻尼Newton数值饱和，完整初次输出保留，仅做backtracking数值修复。独立offline reader第一次将phase-relative路径当repo-relative，失败回执保留并计费，修正reader后只重跑零物理评分，物理搜索和online数组未重跑/改写。

唯一结构witness dS略高于.08，门槛加倍不计结构，且不是新孔。p250077已有真实洞仍未覆盖。两轮paired相同world不是独立复现，三parent和三步bandit不支持泛化/策略优越性。

所有方法Top4/8短缺，MAE最优仍base。G归一化coverage较宽不能替代原RMS目标。只有合法base的三个weak条件被纳入G，新任务不撤销原其它五条件失败。

G实验冻结时未接入T-guided，不能声称完成五臂比较；发布时T已公开交接，记录独立保留。有限Haar/局部mask并非全材料域完备支持；无新3D/GPU/NN/RL/A17消融，没有完整posterior概率或全空间coverage。
''')
    write('RESOURCE_SUMMARY.md', f'''# Research G 独立成本

|费用|实际|硬上限|
|---|---:|---:|
|新full-state attempts|{resources['new_full_state_attempts']}|384|
|complex full-solver RHS|{resources['new_complex_full_solver_RHS']}|12288|
|受控数值wall秒|{resources['controlled_wall_seconds']:.6f}|10800|
|CPU秒|{resources['CPU_seconds']:.6f}|单列，不与wall互换|
|GPU|0|0|

三个独立base验证=3states/24primalRHS；72轨迹各raw+两步+独立验证=288states。共291次LU；2328primal RHS+4608receiver-transfer adjoint RHS=6936。没有隐藏JVP/VJP/newtruth solve。已保存background数组只读setup，wall仍计入，不称免费物理加速。

E0新state/RHS0。每个方法每轮36states/864RHS，公共setup/base只记一次；per-method随机顺序、原始wall及全部prefix见CSV。完整job wall包含准备/优化/序列化/独立full验证，nested ACTION不再相加。

总数值时间包含E0的独立回执、首次失败offline路径桥接、最终独立评分与75图、侧车导出及witness分析；子回执和rollup不重复计费。交互阅读/写代码/文档/上传时间不伪作数值CPU。当前物理solver异常0，数据不相容54端点、QP未收敛24条另列，所有尝试已付费。

剩93states/5352RHS，不用来救结果。T192states、旧A30/A/B/C1/P2账本不动，无预算转账。source/COST_LEDGER+FAILURE_LEDGER为追加记录，RESOURCE_SUMMARY为最新快照。没有SHA检查或发布前重跑实验。
''')
    write('REPRODUCIBILITY.md', '''# 复现与数据布局

本仓库保持历史A/B/P2，G位于research/research_g_20261015。现有Gaussian/.venv_nn Python/NumPy/SciPy/Matplotlib，BLAS threads1，nice10，无安装GPU/新框架。完整源包、raw/corrected/independently-validated χ/current/field/data、observations、所有失败及cost随同仓库Release ZIP；解压到仓库根保持相对路径。

code/e0_trajectory_audit.py、e0_archive_repair.py与e0_history_signal.py只读旧保存资产；cache replay不能用于产生新物理速度/反事实搜索主张。原资产revision见SOURCE_MANIFEST。历史A30源和输入需要同级Gaussian/A30/a30-maxwell-material-manifold；不得因找不到资产偷偷新跑solver。

code/run_g_e1.py是唯一G付费入口，复用原A30与A corrector，必须使用冻结协议；已完成RESULTS存在时拒绝重跑。不要重置COST_LEDGER或覆盖旧协议来获得新额度。

code/offline_g_score.py必须在独立进程、online冻结后运行，只有它读已暴露OFFLINE_TARGETS；不由online模块导入。online数组与manifest只读。首次reader路径错误保留在results/E1_OFFLINE_path_bridge_failure_v1，后来的评分没有重算物理。

code/export_g_episodes.py对真实保存数组生成GSearchEpisode侧车和合法候选索引，不改A/T native schema。软件35项合同测试通过只证明预算/archive/接口语义，不是成像成功；check_assets只查路径/array key，不是新物理验证。

从头复现实验须新独立输出/授权，所有states/RHS/failed trials计费；普通报告/数组阅读不触发求解。不能将本阶段未用余额自动借给新3D或RL。
''')
    write('START_HERE.md', f'''# Research G：Generative Maxwell Solution Search

[中文网页]({WEB}) · [完整报告](RESEARCH_G_REPORT_ZH.md) · [原始Markdown]({RAW}) · [全部证据]({RELEASE}) · [实际图库](MULTISOLUTION_GALLERY/index.html)

已执行零物理E0与历史同一A30后端的小型E1，291/384states、6936/12288复RHS、GPU0。三个观测的material代表3/3/1、structure代表1/2/1，新增结构只有一个近门槛边界world；两轮重复不加算。全部方法Best-of-4=.04224802，未改善base；p250077真洞仍未覆盖。

判断为NO_GO_SEARCH：不扩大当前算法优越性主张，不训练PPO/GRPO，不进入新3D。保留合法witness、全部失败及未来T接口入口。T-guided尚未在G运行；本轮不是完整五臂T/G比较。

先读RESEARCH_G_REPORT_ZH.md首屏八问，再看E0_REPORT_ZH.md、SEARCH_HISTORY_SIGNAL.md、QD_VS_BASELINES.md和CORRECTOR_ATTRIBUTION.md。GATE_DECISION仅科学判断；NEXT_STAGE_PLAN中更大factorial/RL/3D均NOT_RUN。

原始新候选引用在CANDIDATE_MATERIALS/index.json，72份侧车在GSearchEpisode/。21个valid endpoint ID包含跨轮重复，不能写成21个不同world。全部full预测和内部fields/current为真实已执行数组，candidate-only、概率质量null。
''')
    summary = dict(status='BOUNDED_G_RESEARCH_EXECUTED',decision='NO_GO_SEARCH',
                   same_observation_counts=pooled,method_first_round=method_rows,
                   resources=resources,structural_witness=witness,
                   T_guided='NOT_RUN_NOT_IMPORTED_DURING_FROZEN_G_EXPERIMENT',
                   Best_of_4_all_methods=.04224802130330251,RL='NOT_RUN',general_3D='NOT_RUN')
    write('PUBLIC_SUMMARY.json',json.dumps(summary,ensure_ascii=False,indent=2))
    write('E1_RESULTS.json',json.dumps(dict(status=online['status'],online_results='results/E1_LEGACY/RESULTS.json',
                                          offline_results='results/E1_OFFLINE/METRICS.json',
                                          summary='PUBLIC_SUMMARY.json',resources=resources),indent=2))
    write('SOURCE_MANIFEST.json',json.dumps(dict(
        package='/Volumes/migodam\'s-external-brain/archive/user-home/downloads/Research_G_Handoff_20261015.zip',
        research_A_revision='f1b674eb4893e39126701ed48f3d3dcef3098c3c',
        A30_revision='717d3d26a290248bd2b9aeb55f3d03694b407b84',
        C1_revision='83056881e89329212b16f5be47a4d0dbc4465943',
        G_preexperiment_freeze_commit='355990a',backend='Original A30 frozen saved Green/receiver arrays and source-order descriptor',
        corrector='Original ResearchA e2._gn_proposal,32realHaar',
        A17_revision=None,A17_reason='No independent Git revision invented; source snapshot provenance in source audit',
        source_audit='research/delegated/g_source_audit/SOURCE_AUDIT_ZH.md',
        no_file_digests_computed=True,historical_directories_modified=False,
        newer_T_status='results/T_HANDOFF_STATUS_AT_G_CLOSE.json'),ensure_ascii=False,indent=2))
    answers = ''.join('<li>'+html.escape(s.split('**',2)[-1].replace('**',''))+'</li>' for s in first_answers.splitlines() if s[:1].isdigit())
    trows = ''.join('<tr>'+''.join('<td>'+html.escape(str(x))+'</td>' for x in row)+'</tr>' for row in method_rows)
    write('index.html',f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Research G · 下一份相容材料</title>
<style>body{{font:16px/1.85 system-ui,-apple-system,"PingFang SC",sans-serif;margin:0;background:#f4f6f9;color:#172534}}main{{max-width:1120px;margin:auto;padding:38px 22px}}h1{{font-size:35px;line-height:1.4}}section{{padding:26px;background:white;border:1px solid #dce3ea;border-radius:12px;margin:22px 0}}nav{{display:flex;gap:18px;flex-wrap:wrap}}a{{color:#165c98}}img{{width:100%;height:auto}}table{{display:block;overflow:auto;border-collapse:collapse;font-size:14px}}th,td{{padding:10px;border:1px solid #dbe3eb;white-space:nowrap}}th{{background:#edf3f8}}.warning{{background:#fff3db}}</style><main>
<nav><a href="../">研究总览</a><a href="{REPOSITORY}/blob/main/research/research_g_20261015/RESEARCH_G_REPORT_ZH.md">完整中文报告</a><a href="{RAW}">原始Markdown</a><a href="{RELEASE}">全部证据ZIP</a><a href="MULTISOLUTION_GALLERY/">实际材料图库</a></nav>
<h1>生成更多图，还是发现更多材料世界？</h1><p>Research G · 2026-10-15 · topology-class-free search · candidate-only</p>
<section class="warning"><strong>当前判断：NO_GO_SEARCH，保留一个真实物理witness，暂不扩大。</strong><p>三观测的结构代表为1/2/1，新增只有一个边界world，未发现新孔洞。所有方法Best-of-4仍为0.04224802。仅三个已曝光合格parent，T-guided未运行，不能宣称复杂搜索或学习策略获胜。</p><p><b>291/384 states · 6936/12288 complex RHS · GPU0</b></p></section>
<section><h2>八个直接回答</h2><ol>{answers}</ol></section>
<section><h2>匹配完整物理成本的四臂结果</h2><table><thead><tr>{''.join('<th>'+x+'</th>' for x in ['方法','states','复RHS','合法端点','材料新颖','结构新颖','Best-of-4复MAE','Coverage'])}</tr></thead><tbody>{trows}</tbody></table><p>首轮：各方法3parent×3proposal。第二轮只改archive质量替换，结果相同。共享3base/24RHS全局记一次。单次bandit命中不是跨对象优越性；TopK不足不填invalid。</p></section>
<section><h2>一个实际full-Maxwell相容的边界变化</h2><img src="figures/E1_structural_witness.png" alt="实际base、新full-valid材料与差分"><p>同一观测，base残差20.892，新材料22.318，阈值23.785。新提案经两步修正从89.946降至相容；共享本构与solver缺陷约1.4e-15。dS=.081767、support差10格，孔数均0；门槛加倍则不再计为结构不同。两轮相同world只算一次。</p></section>
<section><h2>raw多样性为何不能替代修正后发现</h2><img src="figures/E0_stage_attribution.png" alt="历史raw、corrected、validated和archive的不同计数"><p>历史径向multi-start的331/384个raw新颖提案回到已知邻域。二维保护格原样保留但未相容，不能笼统解释为填孔。历史对下一结构预测无稳定额外增益。</p></section>
<section><h2>明确保留的失败</h2><p>54/72个新endpoint不相容；24/144个QP未达收敛flag。p250077真实洞仍未覆盖。三个对象是按合法base预先筛选的暴露development，不能与原八条件失败率混同。新尺度coverage不是旧目标达成。</p><p>代码复用原A30二维TM，非新3D。没有训练RL/NN，也没有修改T或借其预算。75张实际原图包含失败；所有current/field/predictions与账本随证据包保留。</p></section>
<section><h2>复现与下一步</h2><nav><a href="{REPOSITORY}/tree/main/research/research_g_20261015">阶段目录</a><a href="{REPOSITORY}/blob/main/research/research_g_20261015/RESOURCE_SUMMARY.md">完整费用</a><a href="{REPOSITORY}/blob/main/research/research_g_20261015/NEXT_STAGE_PLAN.md">下一阶段</a><a href="MULTISOLUTION_GALLERY/">全部实际材料图</a></nav><p>先接T最终Git版本化完整handoff，再用同corrector补齐五臂与跨对象复现。当前不进入PPO/GRPO或非平凡3D。发布状态不改变科学结论。</p></section></main></html>''')


if __name__ == '__main__':
    main()
