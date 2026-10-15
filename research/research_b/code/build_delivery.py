"""Render delivery documents from completed results. Never runs experiments."""
from pathlib import Path
import csv
import html
import json

ROOT = Path(__file__).resolve().parents[1]
REPO_URL = 'https://github.com/migodam/b1-probabilistic-maxwell-imaging'
WEB = 'https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_b/'
RELEASE = REPO_URL + '/releases/tag/research-b-20261015'


def load(name):
    return json.loads((ROOT / name).read_text())


def write(name, value):
    (ROOT / name).write_text(value.strip() + '\n')


def table(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |',
                      '| ' + ' | '.join(['---'] * len(headers)) + ' |'] +
                     ['| ' + ' | '.join(map(str, row)) + ' |' for row in rows])


def main():
    summary = list(csv.DictReader((ROOT / 'results/SUMMARY_METRICS.csv').open()))
    paired = list(csv.DictReader((ROOT / 'results/PAIRED_COMPARISONS.csv').open()))
    e3 = load('results/E3/RESULTS.json')
    cost = load('results/RESOURCE_SUMMARY.json')
    supplement = []
    specs = [
        ('research/delegated/cache_audit/CACHE_MAP.json', 'audit_read_script_wall_seconds', 'audit_read_script_cpu_seconds'),
        ('research/delegated/cache_audit/NUMERICAL_CONTRACT_RESULTS.json', 'wall_seconds', 'cpu_seconds'),
        ('research/delegated/ab_adapter/TEST_RECEIPT.json', 'total_wall_seconds', 'total_cpu_seconds'),
        ('research/delegated/ab_adapter/RESEARCH_A_ASSET_STATUS.json', 'audit_wall_seconds', 'audit_cpu_seconds'),
        ('research/delegated/ab_adapter/REAL_A_SNAPSHOT_RECEIPT.json', 'wall_seconds', 'cpu_seconds'),
        ('research/delegated/ab_adapter/FIRST_FOUR_A_IMPORTS.json', 'audit_wall_seconds', 'audit_cpu_seconds'),
        ('research/delegated/ab_adapter/B_MENU_EXTENSION_RECEIPT.json', 'wall_seconds', 'cpu_seconds'),
    ]
    for name, wall_key, cpu_key in specs:
        data = load(name)
        supplement.append({'receipt': name, 'wall_seconds': data[wall_key],
                           'cpu_seconds': data[cpu_key], 'new_Maxwell': 0,
                           'already_in_main_compute_ledger': False})
    for item in load('research/delegated/ab_adapter/PLOT_REVIEW.json')['figures']:
        supplement.append({'receipt': 'research/delegated/ab_adapter/PLOT_REVIEW.json:' + item['stage'],
                           'wall_seconds': item['wall_seconds'], 'cpu_seconds': item['cpu_seconds'],
                           'new_Maxwell': 0, 'already_in_main_compute_ledger': False})
    extra_wall = sum(r['wall_seconds'] for r in supplement)
    extra_cpu = sum(r['cpu_seconds'] for r in supplement)
    complete = {'main_ledger': cost, 'separate_zero_physics_receipts': supplement,
                'supplement_wall_seconds': extra_wall, 'supplement_cpu_seconds': extra_cpu,
                'known_total_numerical_wall_seconds': cost['all_controlled_wall_seconds'] + extra_wall,
                'known_total_numerical_cpu_seconds': cost['cpu_seconds'] + extra_cpu,
                'conservative_all_known_wall_below_E1_E2_cap': cost['all_controlled_wall_seconds'] + extra_wall < 7200,
                'scope': 'Measured controlled numerical work plus named separate receipts; not LLM time, coding, reading, document composition or upload latency',
                'double_counting': 'Individual TEST_RECEIPT_01..03 and plot stdout files are nested duplicates and are not added again',
                'new_Maxwell_attempts': 3284, 'GPU': 0}
    write('results/RESOURCE_COMPLETE.json', json.dumps(complete, ensure_ascii=False, indent=2))
    write('RESOURCE_SUMMARY.md', f'''# Research B 实际资源汇总

独立授权为本次Research B的4096个full-Mie material-frequency attempts、E1/E2受控wall7200秒。E0新Maxwell必须为0。历史P2/B1/C1/A余额只读、没有转账或重置；GPU0、一般非同心3D新campaign0、硬件测量0。

| 实际物理用途 | 新attempts |
| --- | --- |
| E0缓存排名/噪声模拟 | 0 |
| 新四频625节点bank | 2500 |
| 独立64材料×4频n14审计 | 256 |
| 前8材料×4频n18与独立递推 | 64 |
| 固定noise锚 | 4 |
| 新32parent×6频truth | 192 |
| 非线性共享补偿端点（缓存复用后） | 22 |
| 同一个小孔世界六频审计 | 6 |
| 四个真实A包未来field/新truth/双频corrector/参考audit | 240 |
| **合计** | **{cost['used_attempts']} / 4096** |
| **剩余，不再为消耗额度继续运行** | **{cost['remaining_attempts']}** |

物理调用失败{cost['failed_attempts']}、未完成{cost['unfinished_attempts']}；所有失败/审计/truth调用先记BEGIN、统一计费。kernel本身wall={cost['kernel_wall_seconds']:.6f}秒，不等于完整选测计算成本。

主要COMPUTE_LEDGER共{cost['jobs']}个job，E1/E2受控wall={cost['E1_E2_controlled_wall_seconds']:.6f}/7200秒；包括E0/统计/导出等全部已入账数值wall={cost['all_controlled_wall_seconds']:.6f}秒，CPU={cost['cpu_seconds']:.6f}秒。一个adaptive报告导出job失败，已完整算入；仅修正摘要键，积分没重跑。

另有单独回执的缓存/adapter/早期合同测试/绘图（未进主ledger）：wall={extra_wall:.6f}秒、CPU={extra_cpu:.6f}秒。去重后已知数值总wall={complete['known_total_numerical_wall_seconds']:.6f}秒，CPU={complete['known_total_numerical_cpu_seconds']:.6f}秒。即使把全部这些任务保守归入7200秒上限，仍低于它。没有将重复TEST_RECEIPT_01..03或plot stdout再加一次。

按后端metadata共有{cost['boundary_RHS']}个boundary/interface线性系统与{cost['factorizations']}次factorization记账；recursive实现早期遗漏计数通过追加ACCOUNTING_CORRECTION补1584，原END行保留。这是此Mie实现的内部系统计数，不能叫额外硬件测量，也不是通用PDE source-RHS成本。

环境为现有Gaussian/.venv_nn Apple Silicon CPU，BLAS threads1，受控job nice10，串行物理。没有安装基础设施、启动GPU或干扰并行A。记录只覆盖受控数值任务，不假称LLM、写代码、阅读、写报告或Git上传延迟已计为数值CPU。

物理台架一次共享所有菜单预测与模拟setup，再比较策略；逐方法的measurement blocks是真实模拟消费数，portfolio墙钟不能当各独立policy的优化运行速度。calibration/外部geometry各算一个block，没有真实切频/时间/能量报价。

原始追加账本：results/PHYSICS_LEDGER.jsonl、COMPUTE_LEDGER.jsonl；主摘要RESOURCE_SUMMARY.json，去重附录RESOURCE_COMPLETE.json，递推内部系统说明RHS_ACCOUNTING_APPENDIX.json。E3的FINAL_COST.json为当时快照，阶段最终值以上述最新主ledger为准。
''')
    e0 = load('results/E0_ACTION_RANKINGS.json')
    first = load('results/E1_RANKINGS.json')[0]
    methods = ['none', 'fixed_1.2', 'fixed_1.6', 'fixed_2', 'fixed_2.4', 'fixed_2.8', 'fixed_3.2', 'random', 'uniform', 'raw',
               'classical_Fisher', 'joint_projected_Fisher', 'event_EIG',
               'expected_material_risk', 'robust', 'conflict_diagnostic',
               'calibration_only', 'calibrate_then_risk', 'geometry_only',
               'risk_cost_stop', 'robust_abstain', 'two_step_bundle']

    def metric(stage, mode, policy):
        return next((r for r in summary if r['stage'] == stage and r['mode'] == mode
                     and r['policy'] == policy and r['horizon'] == 'final'), None)

    rows = []
    for p in methods:
        a = metric('E1', 'mode1_candidate_only', p)
        b = metric('E2', 'mode1_candidate_only', p)
        c = metric('E2', 'mode2_reference_selection', p)
        if b:
            rows.append([p, f"{float(a['selected_MAE']):.6f}" if a else '未运行',
                         f"{float(b['selected_MAE']):.6f}",
                         f"{float(c['selected_MAE']):.6f}" if c else '—',
                         f"{float(b['measurement_blocks']):.6f}"])
    policy_table = table(['策略', 'E1 固定625池选择 MAE', 'E2 固定625池选择 MAE',
                          'E2 有限参考辅助选择 MAE†', 'E2 新增blocks均值'], rows)
    pairs = []
    selected_pairs = []
    for r in paired:
        if (r['stage'] == 'E2' and r['horizon'] == 'final'
                and r['baseline'] == 'fixed_3.2' and r['metric'] == 'selected_MAE'
                and r['policy'] in ['robust', 'expected_material_risk', 'event_EIG']):
            pairs.append([r['mode'], r['policy'], f"{float(r['paired_delta']):+.6f}",
                          f"[{float(r['CI95_low']):+.6f}, {float(r['CI95_high']):+.6f}]"])
            selected_pairs.append(r)
    pair_table = table(['选择池/信息来源', '策略', 'MAE差：策略−固定3.2', 'parent配对95%区间'], pairs)
    e3_table = table(['真实A交接方法', 'Selected-MAE', 'Best-of-K MAE', 'Coverage@K ≤.05', '平均K'],
                     [[r['method'], f"{r['selected_MAE']:.6f}", f"{r['best_of_K_MAE']:.6f}",
                       f"{r['Coverage_at_K_005']:.0%}", f"{r['mean_K']:.2f}"] for r in e3['summary']])

    risk_rows = []
    for case in dict.fromkeys(r['case'] for r in e0):
        d = {r['action']: r for r in e0 if r['case'] == case}
        risk_rows.append([case] + [f"{d[a]['risk_reduction']:.6f}" for a in
                          ['repeat_1.6', 'frequency_2.4', 'reference_calibration']])
    risk_table = table(['旧暴露案例', 'repeat1.6 预测风险下降', '2.4 预测风险下降', '校准 预测风险下降'], risk_rows)
    ranking_rows = []
    for action, s in first['set_scores'].items():
        f = first['Fisher'][action]
        u = first['utilities'][action]
        ranking_rows.append([action, f"{s['raw']:.6f}", f"{s['robust']:.6f}", f"{s['conflict']:.6f}",
                             f"{f['joint_projected_Fisher']:.6f}", f"{u['event_EIG']:.6f}",
                             f"{u['risk_reduction']:.6f}"])
    ranking_table = table(['行动', 'raw', 'robust', 'conflict', 'joint Fisher', 'event EIG†', '预测材料风险下降†'], ranking_rows)

    science = '''
本轮结论：**继续研究候选补全与新观测后的corrector；停止在当前径向菜单上扩大复杂选测策略。**

1. **仍难分的世界：**小中心孔与实心可沿合法介电常数和gain补偿保持接近；固定625池本身也遗漏大量接近truth的材料。三步后多数有限兼容集合仍跨结构family。
2. **有用的新测量：**在本台架中，高频尤其3.2通常比repeat1.6、2.4或单独校准更有材料选择价值；这不是任意目标/硬件的频率推荐。校准能通过限制gain间接提供材料信息。
3. **新策略的优势：**32-parent比较没有证实robust、材料风险或EIG胜过强固定3.2。普通Fisher首步与固定3.2重合。四个真实A交接案例中robust比固定3.2好，但差异主要来自一个parent，不能升级为总体优势。
4. **MAE是否下降：**固定池三次3.2使Selected-MAE从0.162279降至0.090885；真实A四例加入连续先验多起点corrector后从0.125480降至0.036170。这两组不是同一人群；后者混合了新信息与额外重建计算，不能单独归因于选测。
5. **真正共享补偿贡献：**full-Mie端点支持部分跨频gain冲突，但很多对照冲突很小。不能把全部频谱收益都归给gain共享；原raw-conflict指标倾向校准且材料结果较差。
6. **遗漏candidate的失败：**删去中心孔family后29/30相关noise episodes未报警。低残差或唯一剩余候选都不能保证没有漏世界。
7. **是否继续测量：**不为消耗余额继续campaign；对仍难分的小孔保留歧义。当前优先改善候选支持，并保留简单固定策略作强对照。
8. **潜在TAP贡献：**可复算的跨频补偿见证、候选缺漏失败机制和实际A→B→corrector闭环是研究素材。复杂策略收益、一般非同心3D、真实仪器优势与强原创性尚未建立。
'''

    report = f'''# Research B：补偿感知主动 Maxwell 测量的实际实验

2026-10-15 · 独立Research B分支 · CPU · 原始A/P2资产只读

[在线入口]({WEB}) · [完整证据]({RELEASE}) · [数据与费用索引](START_HERE.md)

{science}

## 1. 从物理问题到可检验任务

一个散射观测可以由几种不同的材料世界解释。新增频率的价值不只在散射变化大，而在候选还能否使用**自身同一组合法材料/校准参数**解释全部历史。不同世界允许不同gain；同一世界不允许每频任意重置gain。我们先用已保存响应排序，再计算新频率响应，在独立模拟观测后更新候选，最后按真实材料而非熵评分。

范围是三维矢量Mie的**同心、径向、各向同性球体台架**，不是任意3D成像。R=1；四维连续系数prior为两个几何系数[-3,1]、Re ε∈[2.2,2.8]、Im ε∈[.04,.16]上的均匀测度。材料decoder自然含实心、中心孔、环形孔及组合；inactive几何对应真实材料原子，积分保留其prior质量。topology标签仅在事后分析，不输入候选生成/优化。

观测为16复通道，packing为Re16后Im16。每frequency采用固定solid锚点的12% RMS噪声，不能由truth/candidate大小调整。每个episode的一个实gain来自截断N(1,.15²)、支持[.7,1.3]，跨历史持续不变。不同episode独立gain。宽带ε常数是本台架假设，未经真实色散验证。

材料主指标是R=1球内体积加权的Re/Im相对介电常数**分量L1和**，s_ref=1；在两材料壳层边界并集上精确积分。它对极小孔体积不敏感，因此另报事件/结构选择，不能用MAE小代替正确结构。

## 2. 两种信息语义与算法

Mode 1输入full-validated有限材料集合，不赋posterior权重。每个材料在完整历史上profile一个有界共享gain，用最小joint χ²选材料；raw/robust使用full缓存响应。robust在当前best+9兼容材料中最多48个筛选点、材料L1≥.02的所有跨/同事件对上，比较合法gain区间内最难区分的未来白化距离。它是有限集合启发式，**不是严格连续minimax regret**。整个625池始终保存。

Mode 2另用原continuous prior、原likelihood和响应proxy的有限条件近似（2048样本、128预测draw、96初始decision材料）计算EIG/材料风险。正确更新为joint-history likelihood减old-history likelihood，gain只积分一次；不得独立乘两个已边缘化gain的likelihood。Fisher先堆叠旧新Jacobian，再消去同一nuisance。采用参考挑action再对固定池profile的流水线仍是**reference-assisted**，不能称纯candidate-only。

没有训练NN，没有ROM。A17的实际墙钟与A23/A24失败不支持默认加速；本轮直接小型full solver即可。候选数量、profile residual、design prior和概率后验分别标注，绝不转换成residual softmax posterior。

## 3. E0：真正缓存排序，新增Maxwell=0

11项源包理论/概率检查通过后，对七个已暴露P2案例运行了实际排名与新独立噪声模拟。原1.6已在posterior中，不能作为新证据；repeat使用新noise ID。原2.4未用于当前C参考，所以可加入。2.4极化/接收子集、analytic校准与no-measurement也列入inventory；缺缓存动作明确不可用，未偷偷补求解。

{risk_table}

这是有限连续proxy参考的预测风险下降，不是实测MAE。七例中raw、Fisher、joint Fisher、EIG、材料风险和robust均选2.4，conflict均选校准。Mode1七例的单次realized平均MAE却由0.124276变为0.132117；暴露小面板不证明主动测量获益，也不据此删难例。E0随机基线共享同一菜单draw、uniform一步等于预设2.4，均保留限制。

## 4. E1/E2：受控新物理与最多三步闭环

先冻结六频菜单{{1.2,1.6 repeat,2.0,2.4,2.8,3.2}}与方法家族，再生成32个新parent，每parent两份噪声。新增四频full bank、独立审计、truth与noise锚共3016次attempt。初始13/17参考事件差最大.000277；这只支持初始观测。

E1是一新增block；E2复用同32对象，是**暴露后续分析**。E2最多三block，uniform固定1.2→2.4→3.2；random在六频均匀抽取；固定六频全部报告。外部void-volume几何传感器SD=.01、校准各计一个新block，不免费，也不假装来自原散射。测量block只是相同dwell/并行通道的归一化成本，没有硬件秒/能量报价。

{policy_table}

†这些列是已执行算法的实际truth评分，但参考分布/效用数值并未得到完整收敛保证。`mode1_candidate_only`是输出选择池的字段名，EIG/材料风险行的action仍需要Mode2参考。

E2 risk-cost-stop平均2.984375 blocks，几乎未节省测量；二步calibration bundle有明确成本，但不能以两步结果和三步结果直接宣称等成本获胜。固定3.2是E1观察到的强描述性基线，E2前保留，不能冒称事先未知的独立最优频率。

按parent先平均两噪声，再2000次固定seed26101569 bootstrap，差为策略−baseline，负数较好：

{pair_table}

有限32parent和多重描述性比较不支持总体支配。共享prediction setup与portfolio计时不能换算成各策略独立运行速度。所有逐步/逐parent的结构、regret、coverage与成本在CSV/Release，未只保留最优案例。

## 5. 参考数值问题：不能称完整full posterior

固定四parent×两噪声、E1/E2和两政策标签对应16条唯一history。13/17比较中12条未稳，事件差最大.090781；2048样本对17阶事件差最大.186853、已选candidate期望L1风险差.016276。17阶本身也未全部收敛。

一次独立physical-coordinate自适应积分保留同prior Jacobian、inactive原子与shared gain：16条全部完成，只有2条达到原relative2e-6规则，14条在512 cells后仍UNRESOLVED。adaptive对17阶最大事件差.075959、对2048样本.186789。重复频率inverse-variance压缩与原目标身份误差约2e-14。没有放宽精度、21阶/全案例追加或重选获胜policy。

因此，本轮Mode2是独立有限proxy比较，**不是可信完整Maxwell posterior，也不是严格上界**。实际full-Mie truth和固定池选择的MAE仍可复算，负的复杂策略比较不依赖把这些参考概率当真。没有用未收敛的概率数值宣布概率校准成功。

## 6. 补偿机制与应该保留的歧义

12个固定small/large/annular结构、四新频的bounded ε/g joint优化端点经full-Mie验证。小孔r=.08的joint白化距离约.039–.142；large r=.25在2.8/3.2约3.225/3.176；ring旧频已很容易分，不能当作破解原歧义的成功。每个优化值只是best-found上界。大孔3.2同ε端点gain冲突平方约.589；许多其他组接近零。校准本身direct material derivative为零，却通过限制gain给出joint Fisher信息，说明补偿分析有物理意义，但单一冲突分数并非材料任务收益。

又用6次full调用核对**同一个**小孔ε/g世界在全菜单的近似相容性，避免拿不同频率的不同近对证明序贯不可能。对固定两个Gaussian世界及最多三次新增block，链式KL/Pinsker的等两世界设计先验分类错误率下界约.439742。它不是posterior质量或全材料空间定理；只是固定数值物理/噪声/菜单中的可检查见证。完整假设与字段见[SHARED_COMPENSATION_TEST.md](SHARED_COMPENSATION_TEST.md)。

## 7. 候选不完整与实际保留的世界

625池Best-of-K平均.062224，Coverage@K(MAE≤.05)=.40625；只重加权或选择，oracle下限不变。删central/annular/combined支持后重新排名，robust动作分别改变33/13/9个noise episodes。删annular后robust MAE约.28291。central漏解29/30、combined漏解15/16相关episodes未触发残差报警，sentinel不能作coverage证明。

三步后固定3.2/robust/risk的启发式兼容K范围为2–27/1–35/3–29；192条episode-policy中160条仍跨family。这不是192个独立对象。完整625材料、共享gain、profile χ²、history与mask已导出到`REMAINING_WORLDS.npz`，概率为none，未知支持质量为null；不只交一个点估计。

## 8. E3：真实Research A包与corrector，已执行

按lexicographic预先固定confirm_00..03，27个真实full-validated候选。原native schema、execution profile与B schema用version adapter衔接，保留材料数组/单位/solver/observation identity；gain-applied旧prediction先恢复未乘gain的F，防止重复乘。B自己扩展已授权菜单，不改A目录。新噪声与truth只供模拟/评分。

135候选未来响应、24truth、54保护几何ε-corrector端点、24三个连续prior起点的独立4D corrector端点、3参考端点核对，共240次attempt。新材料必须同时通过原1.6与所选新频full validation和原A兼容阈值才安装；原候选始终保留。12个prior起点均产生可用端点，其中2次optimizer未达其收敛标志，保留失败状态；通过full/data检查不等于找到优化全局最小。

{e3_table}

参考comparator†未作完整双频端点验证、不安装A、参考UNRESOLVED。`A_base_only_same_action_eps_corrector`复用base对应的已付费端点，是同action的候选容量对照；选action仍看完整A集合，不能当作纯base主动策略。

robust选3.2/2.8/3.2/2.8；固定3.2的四例均值较差主要由confirm_03贡献。只有四例，不作统计胜利。protected安装20个、prior-multistart安装12个新fullvalidated端点；二者分别与原池合并，不把不同更新结果混成同一方法。新增ε/geometry候选可以改善Best-of-K，固定池选择不可以。prior-multistart的Selected结构只对2/4；小材料MAE不代表全部结构正确。

尚缺“同容量corrector但不新增测量”的四例对照，故改善不能完全拆成信息/优化/初始化收益。本轮不为补故事继续花费剩余额度，也不把B独立corrector冒称Research A已接收并重跑。

## 9. 费用、失败与科研去向

新full-Mie尝试{cost['used_attempts']}/4096，失败物理调用0，GPU0；剩余{cost['remaining_attempts']}不消费。E1/E2受控wall {cost['E1_E2_controlled_wall_seconds']:.3f}/7200秒；全部已入账数值任务CPU {cost['cpu_seconds']:.3f}秒。最终值见持续追加账本与[RESOURCE_SUMMARY.md](RESOURCE_SUMMARY.md)，文档/上传时间不冒充物理执行时间。一次adaptive报告字段KeyError计为失败job，45秒已算入，导出修复不重算积分。

32-parent策略结论支持“新增测量有用”，不支持“复杂选测更优”；四个真实A案例支持候选补全值得做更公平的后续试验，尚未证明普遍优越。下一步优先有界候选发现/修正与保持歧义；当前径向任务不继续扩展复杂policy或ROM。

## 10. 文献与新颖性边界

源包已存在主动电磁/Bayesian design近邻，本次只定向核对五篇，没有重做宽泛综述。Fouda/Teixeira的[作者全文](https://arxiv.org/html/1401.1092)已有adaptive sensing与uniform比较；Sandhu等的[正式论文](https://doi.org/10.1109/ACCESS.2023.3288997)讨论EM imaging EIG设计；Bartuska等的[作者摘要](https://arxiv.org/abs/2112.06794)已有nuisance-aware design；Naghshvar/Javidi的[作者摘要](https://arxiv.org/abs/1203.4626)已有active sequential hypothesis testing与成本。摘要阅读不冒充全文复核。

标准Schur、Bayes risk、复合假设、KL链式法则和Pinsker不是本轮新定理。最终新颖性结论保持未建立；可再研究的是Maxwell-specific补偿和漏候选机制是否在更广问题中产生真实、付费值得的材料收益。历史A17/A23/A24/A30的负结果没有被本实验覆盖或推翻。
'''
    write('RESEARCH_REPORT_ZH.md', report)
    write('MATERIAL_RISK_AND_MAE.md', f'''# 材料风险、实际MAE与费用

{policy_table}

{pair_table}

Mode1字段指固定625池的profile选择；event EIG/material-risk行的选测仍依赖外部参考。Mode2有限连续proxy参考未通过新history的完整积分收敛，不能称upper bound。raw full-Mie truth评分与参考概率可靠性分别判断。

Material MAE是精确体积加权Re/Im ε分量L1和、s_ref1。Selected-MAE不是Best-of-K；固定池后者均值.0622242099、Coverage@K≤.05为.40625，选测不能改变它。结构/事件正确率另存CSV。材料loss是决策任务，本轮没有引入其他工业下游损失。

32parent是统计单位；先平均两noise，再固定2000次bootstrap，不把64episodes当64对象。E2同E1人群且已曝光；四例A不与之合并。block成本不是秒/能量，当前没有硬件成本最优性主张。

## 真实A候选的独立四例

{e3_table}

protected/prior-corrector只安装双频fullvalidated且数据兼容的新端点，才允许改变Best-of-K。原池始终保留。四例不足以证明policy优势，也缺同容量无新测量corrector对照。参考comparator未验证完整history，不安装材料。

逐状态原表为results/E1_POLICY_RESULTS.json、E2_POLICY_RESULTS.json、E3/RESULTS.json；逐parent与配对表见results/PARENT_METRICS.csv和PAIRED_COMPARISONS.csv。全部随Release提供。
''')
    write('ACTION_RANKING_RESULTS.md', f'''# 实际行动排名

## E0：缓存full端点与连续proxy参考

{risk_table}

七例的最高材料风险收益都为2.4；raw/Fisher/joint Fisher/EIG/robust相同，conflict都选calibration。预测风险来自有限proxy posterior，不是实测MAE。全部动作含分通道控制的数值见results/E0_ACTION_RANKINGS.json与E0_CRITERION_COMPARISON.json。所有0新增Maxwell。

## E1固定第一个parent、第一noise、第一block的完整数值

这是预先固定首例的展示，未按表现挑选；所有64episode×全部动作均在raw rankings。

{ranking_table}

分数单位/尺度不同，不能跨列比大小。每列内部越大越偏好；risk reduction/EIG属于有限参考†。本例raw、robust、Fisher、EIG、risk均选3.2，conflict选calibration。校准direct material sensitivity为0，但joint Fisher非0，因为约束gain。

E1中fixed3.2与Fisher/risk首步完全同选；E2后续会分开，但没有稳定材料MAE优势。零额外measurement是baseline，停止准则与成本假设见configs/E2_FROZEN.json。

## 真实A交接

confirm_00..03 robust动作依次3.2、2.8、3.2、2.8；reference-assisted全部3.2。每个导入包有完整rankings、新观测与updated_data_history。没有A posterior weights；未把任意等权视为真实posterior。

[主材料结果](MATERIAL_RISK_AND_MAE.md) · [完整库存](E0_ACTION_INVENTORY.json)
''')
    write('SHARED_COMPENSATION_TEST.md', '''# 共享补偿：物理见证与边界

每个世界可以有自己的ε和gain；该世界必须对整段history沿用同一个有界gain。首先堆叠旧新Jacobian再消去nuisance，不能分别消元后把信息当作joint信息。非线性远分支则优化合法补偿并验证full端点。

## 已运行的三类机制

**A：新增频率破坏部分合法补偿。** 固定central r=.25R，旧频优化proxy距离.89991；加入2.8后的joint full距离3.22530，3.2为3.17618。大孔3.2在同一ε端点的shared-vs-independent gain差的平方约.58907。但ε频谱本身也有贡献，不能将整个距离提升全部归给gain。

**B：散射变化而材料仍相近。** central r=.08R，分别新增1.2/2.0/2.8/3.2的joint full距离.03910/.06859/.14027/.14227。分频得到的优化材料不同，因此这些数值本身不证明三步也无效。另行同世界审计如下。

**C：间接校准价值。** known-reference校准不直接依赖未知材料，raw material差为0，但可通过限定shared gain改变joint Fisher与材料风险。例如E0 solid校准预测风险下降.019930、repeat .016637；2.4仍更好.032561。并非校准总优于频率，raw conflict单指标全选校准反而表现差。

共12个joint端点完整full-Mie验证。bounded优化是best-found，不是全局最小；独立gain/只新频等额外优化对照为proxy，未宣称全部fullvalidated。详见results/compensation/COMPENSATION_MECHANISMS.json/.csv及原优化trial。29次callback请求、缓存复用后22次真实新attempt，失败0。

## 固定同一小孔世界的三步上限见证

独立6次full审计拟合一个共同ε/g的r=.08R世界，对全六频的白化距离约[.02590,.05705,.07193,.09519,.13138,.13521]。calibration与独立geometry动作也列入同一有限菜单。

对这两个固定世界、已声明的相同Gaussian协方差以及只依赖历史的策略，KL链式法则给：

`KL(P_history || Q_history) ≤ (d_old² + T max_a d_a²)/2`。

取T≤3，Pinsker给等两世界**设计先验**最小分类错误率下界`.5*(1-min(1,sqrt(D²)/2)) = .4397423865`。不把design prior叫posterior。该见证是固定两世界与菜单的数值结论，有限Mie误差和所声明噪声/constant ε/gain模型均是条件；不是一般材料空间、现实仪器或所有可能测量的不可能定理。

完整原始材料、距离、full字段与预算：results/SMALL_CAVITY_HORIZON_WITNESS.json/.npz。体积MAE对微孔很不敏感，结构判别错误率应单列。
''')
    write('CANDIDATE_COVERAGE_STRESS.md', '''# 候选覆盖压力与保留歧义

本轮不把有限候选集合成功当成全材料空间覆盖。full625池的Best-of-K MAE=.0622242099，Coverage@K(MAE≤.05)=.40625；只有新增材料才可改变这两个值。

## 已执行删除与模型错配

同32parent、每parent两噪声，明确exposed followup；按删除后的池重新排名。full625、central_deleted125、annular_deleted125、combined_deleted225。固定2.4/3.2、raw、robust、sentinel与不测量都保留。

robust一block平均Selected-MAE：full约.103714，central_deleted约.12301，annular_deleted约.28291。对比完整池，推荐动作分别改变33/64、13/64、9/64个noise episodes。它们不是64独立对象。

central缺失29/30、combined缺失15/16相关episodes没有residual alarm。sentinel阈值来自χ²启发式，经历有限池选择与bounded gain，不具有声称的nominal coverage；无报警不保证世界齐全。固定3.2 sentinel不是新ranking胜利。

另有独立instrument drift负控制：future gain = old gain*(1+.03*(k−1.6))，原field与noise不变，拟合器故意仍假设shared gain。该表测试模型错配，不和in-model收益混合，不冒称材料色散实验。

未知结构family不在此有限径向prior支持内，真正的unseen nonconcentric family未运行。当前没有系统的prior-sensitivity新campaign，只有保留原prior和有限候选删除/扰动控制；不能据此称prior robust。所有原始删除、ranking、parent配对与drift记录在results/E2_COVERAGE_STRESS.json。

## 实际交付多个剩余世界

results/REMAINING_WORLDS.npz/.json保存整个625材料及固定3.2、robust、risk三策略的32parent×2noise末历史；每世界shared gain、χ²、材料参数、posthoc events/families和best+9兼容mask均保留。768个已执行trajectory step重放的χ²差为0；没有新噪声或新物理。

兼容K范围2–27、1–35、3–29；192条episode-policy中160条仍跨family。mask仅是有限目录的数据兼容启发式，不是posterior/confidence region。概率none，outside mass=null，整个未筛目录也保留，方便后续独立探索。

E3三个连续prior起点正确补入了一些接近truth的fullvalidated材料，使四例Best-of-K改善；仍不能保证发现全部结构。MAE改善与结构正确分开报告，不能把新ε变体个数称为新结构family。
''')
    write('GATE_DECISION.md', '''# Research B 科学判断

本文件仅记录证据与下一步选择，不是额外执行/发布审批。

| 科学问题 | 本轮实际判断 |
|---|---|
| E0真实缓存ranking | 已执行；新增Maxwell0；七例暴露、非独立确认 |
| E1新频率数值与材料收益 | 已执行；限定径向模型中新增测量降低Selected-MAE；非硬件实验 |
| 复杂compensation-aware策略优于固定/Fisher | 未获支持；当前径向菜单不扩大复杂policy |
| E2三步闭环与停止 | 已执行；同E1暴露parent；risk停止几乎不节约block |
| shared compensation机制 | 部分full端点支持；不等同全部材料收益；存在难分小孔 |
| 全材料世界/遗漏family保证 | 未建立；漏候选常无报警，保留歧义 |
| 新history完整reference posterior | UNRESOLVED；14/16 adaptive history未达到原准则；不作上界/校准成功主张 |
| 真实A导入与独立corrector | 已执行四例27原候选；有fullvalidated新增材料和Best-of-K改善；描述性且信息/优化收益未完全分离 |
| Research A已接收新数据并自行再生成 | 未执行；B只交付可导入的新history与独立corrector输出 |
| ROM inclusive加速 | 未运行；当前直接solver足够，不假定ROM收益 |
| 一般非同心3D/真实硬件/强TAP原创性 | 未验证/未建立 |

结论为：**GO：有限候选发现/补全与保留歧义；NO-GO：为当前径向菜单扩大复杂选测campaign。** 未用812次额度不自动转给下阶段。允许PRESERVE_AMBIGUITY / STOP_NOT_WORTH_COST，不为复杂策略追求胜出而继续修改规则。
''')
    write('START_HERE.md', f'''# Research B：从多个材料世界到有价值的新测量

[在线中文报告]({WEB}) · [完整报告](RESEARCH_REPORT_ZH.md) · [全部原始证据ZIP]({RELEASE})

**实际完成E0、E1、E2和四例真实A包交接。** 新Maxwell 3284/4096，GPU0。复杂选测在32parent上未稳定超过固定3.2；加入独立corrector才明显改善四例A的候选支持。新观测history参考仍UNRESOLVED，不宣称full posterior或一般3D成像。

建议按以下顺序阅读：

1. [中文研究报告](RESEARCH_REPORT_ZH.md)：首屏八个问题、真实结果和解释。
2. [科学判断](GATE_DECISION.md)、[下一阶段](NEXT_STAGE_PROPOSAL.md)、[失败](FAILURES.md)。
3. [实际ranking](ACTION_RANKING_RESULTS.md)、[材料MAE](MATERIAL_RISK_AND_MAE.md)、[补偿机制](SHARED_COMPENSATION_TEST.md)、[漏候选](CANDIDATE_COVERAGE_STRESS.md)。
4. [费用](RESOURCE_SUMMARY.md)、[复现](REPRODUCIBILITY.md)、[源资产](SOURCE_MANIFEST.json)、[协议差异](PROTOCOL_RECONCILIATION.md)。
5. [A→B数据交接](AcquisitionRecommendation/README.md)、`code/ab_adapter.py`及tests；`results/REMAINING_WORLDS.json`保留多个候选，权重none。

结果入口：`results/SUMMARY_METRICS.csv`、`PARENT_METRICS.csv`、`PAIRED_COMPARISONS.csv`；E3另有`results/E3/SUMMARY_METRICS.csv`及`PARENT_METRICS.csv`。大JSON/数组/forward caches随证据ZIP提供，仓库小索引链接不等于缺失实验。

历史P2、A、B1/C1配置/账本不改。源Research Pack原指令从属于本次用户任务；原schema保留。没有SHA检查、额外发布前实验或审批流程。
''')

    public = {'status': 'EXECUTED_BOUNDED_RESEARCH_NOT_GENERAL_3D',
              'E1_E2_final': [r for r in summary if r['stage'] in ['E1', 'E2'] and r['horizon'] == 'final'],
              'selected_paired_comparisons': selected_pairs,
              'E3': e3['summary'], 'E3_parents': e3['parent_summaries'],
              'new_history_reference': 'UNRESOLVED; finite proxy not upper bound',
              'cost': cost,
              'decision': 'candidate support/refinement first; no escalation of complex policy on this menu'}
    write('results/PUBLIC_SUMMARY.json', json.dumps(public, ensure_ascii=False, indent=2))

    def htable(headers, rs):
        return '<table><thead><tr>' + ''.join('<th>'+html.escape(x)+'</th>' for x in headers) + '</tr></thead><tbody>' + ''.join('<tr>' + ''.join('<td>'+html.escape(str(x))+'</td>' for x in r)+'</tr>' for r in rs) + '</tbody></table>'

    screen_answers = ''.join('<li>'+html.escape(line.strip().split('**', 2)[-1].replace('**',''))+'</li>'
                             for line in science.splitlines() if line[:1].isdigit())
    write('index.html', f'''<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Research B · 多世界材料与主动Maxwell测量</title>
<style>body{{margin:0;background:#f4f6f9;color:#172330;font:16px/1.85 system-ui,-apple-system,"PingFang SC",sans-serif}}main{{max-width:1160px;margin:auto;padding:38px 24px}}h1{{font-size:36px;line-height:1.35}}h2{{margin-top:40px}}nav{{display:flex;gap:18px;flex-wrap:wrap}}a{{color:#175b94}}section{{background:white;padding:28px;margin:22px 0;border:1px solid #dfe5eb;border-radius:12px}}.notice{{background:#fff4dd}}.metrics{{display:flex;gap:26px;flex-wrap:wrap}}.metrics b{{font-size:29px;display:block}}table{{border-collapse:collapse;width:100%;display:block;overflow:auto;font-size:14px}}th,td{{border:1px solid #dbe2e9;padding:8px 12px;text-align:left;white-space:nowrap}}th{{background:#edf3f8}}img{{width:100%;height:auto}}.small{{font-size:14px;color:#526271}}code{{background:#edf1f5;padding:2px 5px}}</style></head>
<body><main><nav><a href="../">研究总览</a><a href="{REPO_URL}/blob/main/research/research_b/RESEARCH_REPORT_ZH.md">完整中文报告</a><a href="https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_b/RESEARCH_REPORT_ZH.md">原始Markdown</a><a href="{RELEASE}">完整证据ZIP</a></nav>
<h1>多个材料世界，需要怎样的新测量？</h1><p>Research B · 补偿感知主动Maxwell测量 · 2026-10-15</p>
<section class="notice"><strong>结论：优先改善候选覆盖，保留难以区分的世界。</strong><p>新增测量在径向矢量Mie台架中有用，但32-parent比较没有证实复杂策略优于固定3.2。真实A四例加入独立corrector后明显改善材料MAE；不能把四例或未收敛参考概率当作普遍优势。</p><div class="metrics"><div><b>3284 / 4096</b>新增物理attempt</div><div><b>32 + 4</b>B parent + 真实A交接</div><div><b>0</b>GPU / 硬件测量</div></div></section>
<section><h2>八个直接回答</h2><ol>{screen_answers}</ol></section>
<section><h2>相同测量成本的真实材料结果</h2><p>精确体积加权Re/Im ε分量L1；E1一新增block，E2最多三block，同32parent的暴露后续。有限参考列†未完成新history积分收敛。</p>{htable(['策略','E1固定池MAE','E2固定池MAE','E2有限参考MAE†','E2 blocks'], rows)}<p class="small">Fisher/risk选择与固定3.2首步重合。固定3.2是必须保留的强描述性baseline，不代表一般频率最优。</p><img src="figures/E2_selected_MAE_modes.png" alt="E2材料MAE比较"><img src="figures/E2_measurement_cost_curve.png" alt="测量blocks与材料MAE曲线"></section>
<section><h2>真实A候选：选择和发现必须分开</h2>{htable(['方法','Selected-MAE','Best-of-K','Coverage≤.05','平均K'], [[r['method'], f"{r['selected_MAE']:.6f}", f"{r['best_of_K_MAE']:.6f}",f"{r['Coverage_at_K_005']:.0%}",f"{r['mean_K']:.2f}"] for r in e3['summary']])}<p>固定候选选测不改变Best-of-K；新增fullvalidated corrector端点才可改善。仅四parent，缺同容量无新测量corrector对照，材料误差下降不等于所有结构正确。参考comparator未完整history验证，不能当上界。</p></section>
<section><h2>四例材料选择与候选覆盖</h2><img src="results/E3/E3_MATERIAL_COMPARISON.png" alt="四个真实A案例的材料MAE和Best-of-K；描述性比较，参考未收敛"></section>
<section class="notice"><h2>明确的负结果与数值边界</h2><p>漏掉中心孔family时29/30相关noise episodes没有残差报警。三步后多数兼容集合仍跨结构family。独立adaptive审计14/16 history未达原积分准则；有限2048样本事件质量差可约0.187。因此不宣称完整posterior、全空间coverage、一般3D成功或复杂策略优势。</p><p>只读复用P2和A资产；失败、数值audit、truth响应统一计费；剩余812 attempts不为填满预算继续运行。</p></section>
<section><h2>报告、代码与数据</h2><nav><a href="{REPO_URL}/tree/main/research/research_b">阶段目录</a><a href="{REPO_URL}/blob/main/research/research_b/GATE_DECISION.md">科学判断</a><a href="{REPO_URL}/blob/main/research/research_b/RESOURCE_SUMMARY.md">实际费用</a><a href="{REPO_URL}/blob/main/research/research_b/FAILURES.md">失败与限制</a><a href="{RELEASE}">全部原始证据</a></nav><p>代码、冻结配置、科学合同测试、逐parent统计、实际推荐与更新history均保留。数组、完整rankings和缓存打包在同一B1仓库Release。发布状态不改变科学结论。</p></section>
</main></body></html>''')


if __name__ == '__main__':
    main()
