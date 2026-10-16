"""Render research documents from closed, saved diagnostics; zero physics."""
from pathlib import Path
import json,sys,subprocess
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from finite_feedback.common import OUT,dump

def write(name,s):(OUT/name).write_text(s.strip()+'\n')
def main():
    m=json.loads((OUT/'METRICS.json').read_text());p=json.loads((OUT/'POSTHOC_CONSISTENCY.json').read_text());e=json.loads((OUT/'E0_COMPLETE.json').read_text())
    dec={r['anchor']:r for r in m['decomposition']};fb={r['endpoint']:r for r in m['feedback']}
    legal={r['anchor']:r for r in m['paths'] if 'rows'in r};geo={r['anchor']:r for r in m['geometry']}
    replay=json.loads((OUT/'evidence/TRAJECTORY_REPLAY.json').read_text())
    table='| Anchor | GT材料RMS | 白化一阶响应 | 实际full差异 | 弱能量@.001/.01/.1 | 原观测残差 |\n|---|---:|---:|---:|---|---:|\n'
    for k in ['physics_R_base','N3']:
        r=dec[k];a='/'.join(f"{100*x['weak_including_rectangular_null_fraction']:.2f}%" for x in r['spectral_thresholds'])
        table+=f"| {k} | {r['GT_material_RMS']:.6f} | {r['linear_response']:.6f} | {r['actual_full_prediction_difference']:.6f} | {a} | {r['anchor_whitened_residual']:.6f} |\n"
    orders='| GT全程 | 一阶相对误差 | 二阶 | 三阶 | 六阶 | 二十阶 | rational |\n|---|---:|---:|---:|---:|---:|---:|\n'
    for k in ['physics_R_base','N3']:
        r=fb[k+'_oracle_t1'];a=r['relative_prediction_errors']
        orders+=f"| {k} | {a['order_1']:.6g} | {a['order_2']:.6g} | {a['order_3']:.6g} | {a['order_6']:.6g} | {a['order_20']:.6g} | {r['rational_relative_error']:.3g} |\n"
    path='| TRAIN anchor，prior步长.75 | 方法 | 材料RMS位移 | 白化prediction位移 | 观测残差 | 显著新结构 |\n|---|---|---:|---:|---:|---|\n'
    for k in legal:
        for method in ['frozen_weak','second_order_curvature','updated_J_corrector']:
            r=next(x for x in legal[k]['rows'] if x['method']==method and x['t']==.75)
            path+=f"| {k} | {method} | {r['material_RMS']:.8f} | {r['data_displacement']:.6f} | {r['whitened_residual']:.6f} | {'YES' if r['structure_qualified'] else 'NO'} |\n"
    curve='| TRAIN chart | Maxwell二阶项范数 | decoder曲率项范数 | decoder/Maxwell | t=1冻结chart响应 | t=1实际响应 |\n|---|---:|---:|---:|---:|---:|\n'
    for k in geo:
        g=geo[k];r=next(x for x in p['chart_vs_Maxwell'] if x['anchor']==k and x['t']==1)
        curve+=f"| {k} | {g['Maxwell_curvature_norm']:.6g} | {g['decoder_curvature_norm']:.6g} | {g['decoder_curvature_norm']/g['Maxwell_curvature_norm']:.1f} | {r['frozen_chart_response']:.6g} | {r['actual_full_response']:.6g} |\n"
    decision=dict(scientific_status='FINITE_RESPONSE_VERIFIED_NO_GT_FIBER_CONNECTION',new_2D_states=45,complex_RHS=1880,
        significant_new_structures=0,GT_connection_established=False,idealized_physics_reconstruction_condition_tested=False,
        candidate_only=True,posterior_probabilities=None,exact_fiber_existence_claim=False,
        next_priority='legal decoder and observable-component-qualified physics anchor before learned tangent/transport',
        new_training=0,new_Mie=0,new_3D=0)
    dump(OUT/'DIAGNOSTIC_DECISION.json',decision)
    write('START_HERE.md','''# Maxwell有限反馈与GT差量诊断 · 2026-10-16

**局部弱能量不等于可连通的材料解族。** 同一公开development对象上，N3–GT差量的98.09%–99.06%属于预设阈值下的弱材料方向；但N3和physics-only R基线均未成为原P0、原噪声下的合格anchor。没有建立通向GT的合法等数据路径。

真实full-Maxwell的有限反馈已核验；高阶显著提高预测精度。两个合法TRAIN chart只得到小幅连续材料变化、已采样端点noise-compatible，0个显著新结构。主要有限步障碍是材料decoder曲率，不是此处的Maxwell共振。

- [完整中文报告](DIAGNOSTIC_REPORT_ZH.md) · [网页报告](report.html) · [数值图与材料图](gallery.html)
- [GT强弱分解](STRONG_WEAK_GT_DECOMPOSITION.md) · [有限反馈](FINITE_FEEDBACK_ORDER_COMPARISON.md)
- [局部路径与曲率](LOCAL_PATH_AND_CURVATURE.md) · [anchor/先验审计](PRIOR_SUPPORT_AND_ANCHOR_AUDIT.md)
- [下一轮Pro应回答的问题](THEORY_QUESTIONS_FOR_PRO.md) · [机器可读数据](METRICS.json)
- [独立资源账本](RESOURCE_LEDGER.jsonl) · [复现及证据数组](REPRODUCE.md)

45/48新增二维states，1880/2048 complex RHS，受控wall 11.146秒；GPU训练/Mie/3D均0。逐endpoint立即保存，没有LXZ原始数组丢失的重演。已完成预设诊断，剩余3 states/168 RHS不追加寻找正结果。旧研究与预算未改。

TRAIN已知材料不是未知成像成功；GT方向只属ORACLE_GEOMETRY_DIAGNOSTIC。无posterior/mode概率声明。
''')
    write('DIAGNOSTIC_REPORT_ZH.md',f'''# Maxwell高阶反馈能否把局部弱方向延伸到Ground Truth？

**本轮回答：高阶反馈可以准确预测有限材料变化，但没有给出连接重建与GT的合法局部解族。** 这个区分是实验的主要结论。原N3仍是强Top-1比较，未训练新网络、未修改A17冻结稿，也没有重新执行T/G/LXZ。

## 从可检验的问题开始

材料变化小而散射变化更小，可能来自观测不敏感，也可能来自材料参数之间重复表示同一物质场。即使真正存在弱材料方向，有限步仍会改变内部场和材料chart；因此要分开比较材料差量、Jacobian预测、完整有限响应，以及先验是否合法。

本轮选择在看新结果之前冻结：原A17顺序首个共同development观测p250065_weak_a1_snr20；原A26 BASE_N2作为非NN的R闭包一步物理基线；同观测N3；原LXZ两份TRAIN smooth材料250033/weak和250034/ring。没有按GT误差或最后结构挑样本。

物理基线从背景adjoint初始化、经一次有界ROM材料更新得到，**不是已收敛full-GN**。同案例独立背景起点的迭代full-GN/TSOM缓存未找到。旧A17 full-GN从N3启动，只能叫PHYSICS_REFINED_N3，不能填补这一空缺。因此，用户提出的“物理重建已理想恢复敏感信息”前提仍未完整实现，不能用本轮负结果否定所有满足该理想条件的重建。

## 真实差量有多弱？

全部使用原24×24 scalar TM/VIE、多源共享材料、原固定gain=1和复噪声。材料度量为sqrt(1/576)[Re h;Im h]，数据按每source的32个Re接32个Im打包，并乘sqrt(2)/sigma白化。阈值为sigma_i/sigma_max≥.001/.01/.1，没有事后调阈值。

{table}

N3剩余材料误差确实大部分落在全复材料空间的弱方向，这与“学习器恢复了较强观测信息”相容，但不构成该机制的因果验证。physics R一步基线的弱比例仅71.91%–76.59%，其115.28量级的一阶响应并不弱；该anchor显然尚未恢复敏感信息。两者都不满足原P0固定增量相位和数据相容性。谱切点附近的gap约1.08–1.28，并无清楚的单一谱断层。

此处的弱比例含矩形Jacobian的未解析null分量，且不是合法P0 tangent能量比例。N3差量在必要固定相位子空间内的能量99.89%，仍不等于原rasterizer/chart可实现。TRAIN材料与自己的GT相同，hGT=0，相关比例为null而非100%恢复率。

## 高阶反馈确实改善了什么？

二维本构X=diag(chi)严格线性，有限增量必须用实际Delta chi。复用完整锚点resolvent后，第一项为R0 DeltaX E0；每一更高项再作用R0 DeltaX G。完整rational形式通过已付费endpoint LU求解，另收8 RHS，并与独立full-forward差值比较。

{orders}

N3→GT的非线性余项范数1.284，与一阶响应夹角cos=-.683，确有部分抵消，但只把响应范数16.402降至15.553，并未抵消成等散射。physics基线→GT的余项23.696主要与一阶响应近正交，完整响应反而由115.283增至117.714。不能把“有高阶项”当作必然存在非线性补偿解。

全部16个有限反馈点的rational最差相对误差{p['rational_worst_relative_error']:.3g}；独立旧N3 Jacobian与本轮相对差{p['N3_J_vs_independent_historical_relative_error']:.3g}；45个新endpoint全部保存，最大方程backward残差{p['max_saved_backward_residual']:.3g}。这些是本轮真实Maxwell诊断，未冒用旧L/Z/X合成代数成绩。

两条GT全程的Neumann范数**上界**分别1.129与1.032，不满足所用充分条件；20阶仍达到数值精度，完整resolvent稳定。上界失败不等于实际算子范数>1、级数必发散或resolvent不存在。沿GT直线更新Jacobian的Simpson积分，4段→8段相对误差：physics 7.37e-6→4.55e-7，N3 2.39e-8→1.50e-9，支持有限差量等于沿途Jacobian积分；这条直线不在已认证原先验内，也不是等数据路径。

## 合法弱方向能走多远？

两个TRAIN anchor分别残差22.313和22.970，低于冻结tau=23.785449，且可由原P0精确重放。沿原8个Gaussian生成参数探索，去掉共同权重缩放这个同材料gauge，只剩7个独立材料方向。两例响应/材料奇异值最小值12.849和297.076，均非零。本轮只检验有限噪声相容集合，不宣称精确等数据流形。

{path}

所有26个合法路径记录均在原noise球内、保持P0_EXACT，三个更新步都完成；但0个显著新结构。表中prior步长一致，材料距离和费用不一致，不能据此给三个方法排同成本胜负。updated-J每anchor另付2×32 adjoint RHS，而冻结路径每点额外做有限阶诊断；这些诊断开销与部署搜索开销分列。

## 有限步的主要弯曲来自哪里？

在chart坐标中，二阶响应同时包括“Maxwell重新反馈”和“材料decoder本身弯曲”。两者不能混成一个物理反馈增益。

{curve}

第一例弱chart方向在t=1的冻结线性响应只有.006429，实际full响应3.73365，幅值相差580.75倍；但把**实际材料增量**交给同一个冻结Maxwell Jacobian，误差只有0.446%。因此这里的巨大失真主要是有限chart步偏离切向，而不是Maxwell高阶散射把弱方向突然放大。另一例chart响应向量误差69.15%，实际增量的Maxwell线性化误差只有1.716%。

二阶normal补偿和更新Jacobian可减少部分数据漂移，但也减少/转动材料位移，尚无新结构与同成本成像收益。Maxwell L条件数约1.61–10.58，没有本轮近奇异/共振失稳证据；chart权重gauge和有限曲率更直接。

## 为什么没有证明通向GT？

首例历史噪声使GT本身残差24.090925，略高于预先冻结的95%界23.785449；这可以是合法噪声尾部，不能放宽sigma、换观测或称model mismatch。N3起点残差28.034654，physics基线120.432293，预定直线的所有采样点仍未进该noise球。更根本地，两个重建不是P0_EXACT，GT直线为OFF_PRIOR_ORACLE_PATH。合法TRAIN路径则从已知材料自己出发，不能成为从未知重建恢复GT的证据。

## 七个问题的裁决

1. **GT差量的弱比例？** physics一步基线71.91%–76.59%，N3 98.09%–99.06%；仅全材料、预定截断下的oracle诊断。
2. **能否有限变化？** 两个合法TRAIN chart上可以作小幅连续材料变化，全部已采样端点保持原noise球相容；最大冻结RMS位移.002636/.013231，0显著新结构，不能叫精确fiber。
3. **高阶feedback改变什么？** 显著修正有限响应预测；并未把已测GT差量变成等散射方向，亦未解决原先验合法性。
4. **一阶弱与有限等散射差多远？** 第一例chart线性响应被实际finite响应超过580倍，主要因decoder曲率。两例没有非零精确material null被验证。
5. **重建与GT连接证据？** 没有。ORACLE直线的Jacobian积分一致性不是合法等数据路径。
6. **障碍来自哪里？** 重建anchor不合格、原P0表示/相位、chart曲率/gauge；本轮Maxwell算子稳定，高阶预测可控。未完整实现理想physics-only恢复前提。
7. **下一步学什么？** 先建立合法decoder及敏感分量合格的physics anchor，再问短程finite transport或tangent field；当前不足以启动新Flow/神经生成器或声称应扩大posterior kernel。

## 费用、失败和可复现性

45/48 states、1880/2048 complex RHS，其中full-state forward360、receiver adjoint256、反馈系列/有理式1120、oracle沿路Jv128、chart曲率16；总forward RHS为1624，adjoint为256；LU45次。受控wall {m['cost']['controlled_wall_seconds']:.6f}秒，包含SVD与endpoint保存；不是纯solver latency。额外零物理分析与绘图另列，GPU/Mie/新3D均0。剩余额度不追加追逐正结果。

代码review发现方向定号使用了初始Dg而非前一个Dg；保存数组的独立重放表明6次sign均相同、12个predictor/corrector材料endpoint逐值重现，**实测结果未受影响**。弱/环形的真实切向转角分别约.65°/1.06°与15.38°/6.62°。曲率中心差的影响<8e-8。公开修正后的代码、实际执行revision d74f747和原源码快照均保存，不重跑来隐藏问题。原执行器在解算准入处检查wall，未有独立中断；实际11.15秒远低1200秒，当前代码补充deadline，关闭账本仍禁止自动再运行。

本轮没有新多解算法、posterior calibration或IEEE TAP贡献完成的主张。最有用的新诊断是：**把prior decoder曲率与Maxwell反馈曲率分开，否则会把chart离开弱切向误判成高阶Maxwell纤维机制。** 需在合格反演anchor上另行检验，不能外推到一般三维。
''')
    write('STRONG_WEAK_GT_DECOMPOSITION.md',f'''# GT差量的强/弱材料分解

仅ORACLE_GEOMETRY_DIAGNOSTIC。GT不进入任何在线方向、corrector或候选选择。

{table}

## 测度与阈值

材料坐标m=sqrt(v)[Re h;Im h]，v=1/576，因此范数为复材料RMS。J=W DF/sqrt(v)，W使用sqrt(2)/原复sigma。512行、1152列；数据packing每source的32Re/32Im。对J=U diag(s) Vᵀ，强投影V_s V_sᵀm，强能量为其平方范数/||m||²，弱能量为补集，包含矩形nullspace。相对cutoff .001/.01/.1预先冻结。不是参考后验measure，也不是“可自由生成”的方向数。

N3和physics基线的原P0 chart缺失且phase必要条件失败，故真正prior tangent占比NOT_ESTABLISHED。固定增量相位子空间只是必要条件：N3差量99.8896%、physics差量99.2166%的能量在该线性空间内，并不解决有限先验支持/峰值/rasterizer约束。TRAIN chart存在、rank7，但其hGT=0，能量比例null。

各切点的谱隙和可观测condition在METRICS.json逐项保存；不存在明显唯一cutoff。完整谱、V和GT系数在results/*_DECOMPOSITION.npz。Full J与旧N3独立缓存相对差6.75e-16。不能因全空间弱比例高就推定原P0里的真实等数据流形。
''')
    write('FINITE_FEEDBACK_ORDER_COMPARISON.md',f'''# 一阶、有限阶与完整反馈的真实Maxwell比较

本轮新二维证据，非旧L/Z/X合成核验。完整场均共享同一chi，无任意receiver-dark current候选。

L0=I-X0G，R0=L0^-1，E0=incident+G j0。A30的a(chi)=chi，DeltaX=diag(chi1-chi0)。

q1=R0 DeltaX E0；q(n+1)=R0 DeltaX G qn；Delta j=(I-R0 DeltaX G)^-1 q1=L1^-1 DeltaX E0。

每个q项支付8复RHS；rational复用endpoint的LU，再支付8。没有显式R0逆、576列隐藏求解。1/2/3/6阶作用于真实finite材料增量，不是对非线性decoder坐标的Taylor阶数。GT终点额外记录20阶。

{orders}

N3的高阶项与一阶部分抵消，幅度减少约5.18%，但最终响应15.55非零；physics基线的高阶项并未帮助抵消。当前实验没有出现抵消至等散射的GT方向。

## 收敛与稳定性

采用可计算充分上界||R0 DeltaX G||2≤||DeltaX G||2/sigma_min(L0)。GT全程两上界1.129/1.032，检查不通过，只能INCONCLUSIVE；不是实际norm或spectral radius>1的测量。20阶的source-response数值收敛不等于所有RHS一致收敛证明。Rational无需这个Neumann条件，实际endpoint条件数、backward残差均记录。未遇到近奇异失稳，不能声称已测试高非正规极端情况。

GT路径0/.25/.5/.75/1及额外midpoints通过独立full-JVP检验平均Jacobian，Simpson4/8差分别支持数值一致；无等数据、先验合法或路径存在性推论。

A17原三维polarizability是非线性的；其有限DeltaX应为a(chi1)-a(chi0)。这里a=chi恰好线性，不把二维结论冒充原三维测试。A17固定材料的resolvent已含全部多重散射；本轮研究的是材料改变后的finite反馈。
''')
    write('LOCAL_PATH_AND_CURVATURE.md',f'''# 合法局部路径与材料chart曲率

两个原P0 TRAIN smooth材料只用于数学校准，不是新反演或GT发现。原8个独立uniform生成参数经Gaussian CDF坐标化，不扩support；共同Gaussian权重缩放在峰值归一化中抵消，所以材料rank7。

{path}

冻结弱方向取每单位材料RMS响应最小的非gauge方向，再归一成prior坐标单位。初始rank7奇异值全非零；两个最小值12.849和297.076。有限noise球内小宽度与精确等数据fiber不同。

二阶恒等式针对F(g(u))：D²(F∘g)[v,v]=D²F[Dg v,Dg v]+DF D²g[v,v]。第一项用真实R0反馈，第二项为decoder曲率；decoder二阶用两步长有限差分，非解析exact curvature。

{curve}

更多tangent→finite材料偏离和向量误差见POSTHOC_CONSISTENCY.json。第一例t=1材料secant/tangent cosine仅.5343，原弱切向已显著转动；Maxwell对实际材料增量的一阶误差却仅.446%。有限弱方向失效主要是decoder弯曲，不是Maxwell突然近奇异。

方法均不读取GT。updated-J每步刷新后保持弱方向定号，再作一次有界法向corrector；非收敛水平集追踪器。第二阶曲线修正和updated-J的材料位移不同；不能在相同prior t下宣称相同材料距离，更不能忽略其额外Jacobian/RHS。26条记录皆noise-compatible但0显著新结构，updated-J累计材料路径长度.000791/.009167，端点净位移.000791/.008743。

独立轨迹重放验证了执行代码定号缺陷没有改变本次六次sign，12个endpoint maxabs重放误差0。真实tangent角变化约.65°/1.06°与15.38°/6.62°。rank与active max cell在采样点不变，不是连续路径全区间的严格rank/noise证明。路径中点未穷尽检查，不称certified path。

known TRAIN自己的GT就是起点，离开anchor只会增大“到GT”的距离；这不构成成像风险改进。N3与physics R基线没有合法P0局部path启动条件，不运行其E3。
''')
    write('PRIOR_SUPPORT_AND_ANCHOR_AUDIT.md',f'''# 先验支持与anchor来源

E0新物理0；原缓存和此次E1真实测量严格分开。

1. **Physics-only**：同观测BASE_N2，背景adjoint初值+一次R-ROM梯度box投影，无NN/GT初始化。E0只找到材料图，无独立完整F/J缓存。此次付费full-state/J后残差120.4323；先验phase失败。不是已充分收敛的物理重建。独立迭代full-GN/TSOM缓存MISSING。
2. **N3对照**：冻结A26已暴露validation，非fresh验证。E0已有完整状态/J，缓存残差28.0347、N3–GT预测距15.5526；新付费重建用于统一factor/receiver接口，明确重新收费。phase失败，不是原P0 anchor。
3. **Known P0**：原LXZ固定两份TRAIN smooth材料，原generator参数逐点重放；残差22.3130/22.9696，均过原95%界。只用于数学合法路径校准。没有将TRAIN已知材料作为development重建。
4. **PHYSICS_REFINED_N3**：A17/N3_APPLICATION_R1 full-GN/ROM均从N3起点，不能冒充独立physics-only。旧结果与MAE不覆盖本轮anchor资格。

原P0有四个family与固定增量相位/原峰值归一化，不等价于任意合法复pixel box。rho=0是原背景，不是绝对空气。smooth条件chart是原P0的条件子支持，不同于整个先验分布；没有新增P1。

首例truth noise残差24.090925>tau23.785449：这是历史公开噪声尾部，维持原sigma/gain/阈值。不通过挑新noise或放宽阈值获得合格GT。GT直线路径除最后endpoint外无原先验证书，统一OFF_PRIOR_ORACLE_PATH。

已读取“逆散射中的Dyson级数”原会话的finite-resolvent/Neumann/rational-Schur讨论，只将其视为理论上下文；私人会话正文不随仓库/Release发布。原始closure_final与A26/A30/LXZ源码证据见evidence/CONSTITUTIVE_CONTRACT.md和CACHE_AUDIT.md，未改冻结论文/历史账本。
''')
    write('THEORY_QUESTIONS_FOR_PRO.md','''# 下一次Pro理论研究应回答什么

先读本次实测数据，再讨论新理论。不要把经典resolvent或常秩结论当作原创贡献。

## 可由本轮约束的问题

1. **先验内有没有精确局部fiber？** 对合法商chart g去掉同材料gauge后，若DF Dg在某邻域满列秩，精确水平集局部为离散点；不能把全pixel复空间的巨大nullity搬过来。本轮两份7维材料tangent数值满秩，只支持有限噪声相容tube。需要声明噪声tube半径、材料度量和均匀余项，不能从一个小singular value推断有限族。
2. **二阶可延续的障碍在哪里？** 等数据曲线需A v=0，并存在a使A a+D²(F∘g)[v,v]=0。后二阶必须包括decoder曲率，且其在Range(A)外的部分不能由normal修正消除。若只近似Av≈0，应给出先验合法步长下的一阶漂移、未补偿曲率和高阶余项共同界。该经典必要条件不是新定理，本轮没有计算其uniform bound或宣称精确fiber存在。
3. **哪种条件能把“GT弱能量高”连接到可恢复性？** 必须先有独立physics-only、原先验内且观测合格的anchor，再度量GT差量是否位于同prior chart的finite可达区域。本次R基线不满足理想前提，N3 off-prior，二者不能验证该假设。需要比oracle能量比例更可观测的敏感分量恢复诊断。
4. **何时学习tangent field，何时学习finite transport？** 当Dg曲率远大于Maxwell反馈曲率，训练tangent vector本身可能无助有限步。先检验material-native retraction/商chart和合法decoder；若几何仍弯曲，才论证有监督或无监督finite transport的收益与物理成本。
5. **A17能增加什么？** 原A17是固定二次材料任务下current representation fidelity，不负责产生物理解。下一步若使用rational Schur加速，应首先比较完整finite响应精度和总成本，而非再次只改善frozen GN gap；当前小型dense solver已稳定廉价，未建立A17独立增量。
6. **分离分支与高阶反馈是什么关系？** 本次resolvent预测精确但没有新增结构，不能把高阶抵消视为跨分支机制。要先有同一冻结prior/观测的合格竞争材料区域，再问全局proposal或posterior kernel。否则先建立ambiguity benchmark或active sensing问题。

## 现阶段停止与保留

保留：finite反馈数值oracle、正确prior/质量度量、decoder-vs-Maxwell曲率分解、原N3 Top-1基线。

暂缓：新NN/Flow/Diffusion、非径向3D、多解成功叙事、通GT路径保证、A17新增成像收益。不得用剩余额度追加正结果。

决定性下一实验应是：从背景独立反演获得同原先验内合格anchor，验证其敏感分量，然后在同material距离/同full成本下比较普通与material-native曲率延续。须新授权；本轮不执行。
''')
    write('REPRODUCE.md','''# 复现与数组

实际物理执行版本：`d74f747`，源码快照 `evidence/source_executed/`；Python见SOURCE_MANIFEST.json。使用已有Gaussian/.venv_nn，不安装新的solver。

历史输入：A30 `inputs/a29_release/inputs/a25/{validation,train,background}`；A26 `results/a26/evaluation/261009/PREDICTIONS.npz`；A17/N3_APPLICATION_R1的同case旧Jacobian；原LXZ worktree两份TRAIN校准NPZ。

E0（只读缓存，无物理）：`PYTHONPATH=src python -m finite_feedback.e0`。

实际物理入口：`PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m finite_feedback.diagnose`。账本已封闭，当前入口自动拒绝重复执行，不得删除RUN_CLOSED/ledger重新计数。新物理复现须独立授权与原始历史保留。

仅分析/画图：`PYTHONPATH=src python tools/analyze_finite_feedback.py`、`python tools/report_finite_feedback.py`、`python tools/plot_finite_feedback.py`。只读已支付数组，无新的F/J求解。

结果数组在原B1公开Release的本轮证据包，含原weak/ring完整已支付background G/S/incident/L_b数组，复制不支付新物理。新授权复现时优先读取本轮inputs/background；仍复用固定版本A30 solver源码。`results/endpoints/*.npz` 45份，逐端点即时保存chi/current/field/prediction/LU/pivots/backward residual。`*_JACOBIAN.npz` 全Jacobian、transfer、质量实化；`*_DECOMPOSITION.npz` 谱与GT系数；`*_FEEDBACK.npz` 1/2/3/6/20/rational实际响应；`*_ORACLE_PATH.npz` 9点Jv与复合积分；`*_CHART.npz` 初始坐标/Dg/弱方向。26份legal路径记录对应13个endpoint/anchor，包含共同的frozen/predictor重复记录，不能叫26个解。

review重放脚本只用保存prediction、Jacobian和原decoder，重现12个updated路径端点，不重新解物理。当前源码修复了方向定号、curvature中心与未来wall deadline；实测数据仍归属原执行源码，未暗中重跑。

没有新probabilities，所有候选仅作diagnostic/candidate-only。已暴露development与TRAIN不称独立验证，原始私人Dyson聊天不公开。
''')
    sources={}
    for name,path0 in [('B1',ROOT),('A26',ROOT.parents[1]/'A26/a26-om-feedback-prior-completion'),('A30',ROOT.parents[1]/'A30/a30-maxwell-material-manifold')]:
        sources[name]=dict(path=str(path0),commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=path0,text=True).strip())
    dump(OUT/'SOURCE_MANIFEST.json',dict(executed_B1_commit=subprocess.check_output(['git','rev-parse','d74f747'],cwd=ROOT,text=True).strip(),sources=sources,
        source_snapshot='evidence/source_executed',background_arrays='inputs/background/{weak,ring}.npz copied from original paid A30 input',native_subagents='medium evidence/plots; high independent formula review; fresh concise context, no model/provider override or DeepSeek',post_execution_fixes='direction transport, decoder curvature center, future deadline; no physical rerun',
        data_scope='exposed development + known TRAIN calibration only',private_chat_published=False,hashes_computed=False,
        Python=sys.version.split()[0],numpy=__import__('numpy').__version__,scipy=__import__('scipy').__version__,
        code_freeze_before_paid_run=True,no_changes_old_ledgers=True))
    dump(OUT/'RESOURCE_SUMMARY.json',dict(**m['cost'],caps=dict(states=48,RHS=2048,physical_wall_seconds=1200),
        remaining_not_spent=dict(states=3,RHS=168),closed=True,GPU_training=0,new_Mie=0,new_vector_3D=0,
        all_forward_RHS=1624,all_adjoint_RHS=256,solver_rhs_breakdown=dict(full_state_forward=360,receiver_adjoint=256,finite_feedback=1120,oracle_Jv=128,chart_curvature=16),
        new_endpoints_durable=45,missing_paid_endpoint_arrays=0,posthoc_analysis_wall=p['wall_seconds'],
        wall_includes_persistence_SVD=True,old_balances_borrowed=0,
        cost_comparison='No matched wall-speed claim; method endpoint/J setup/audit costs separated in figures/endpoint_costs.json'))
    print(json.dumps(dict(reports_written=True,new_physics=0,decision=decision)))

if __name__=='__main__':main()
