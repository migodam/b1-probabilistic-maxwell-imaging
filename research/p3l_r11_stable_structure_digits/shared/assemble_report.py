from core import *
import csv,shutil,datetime

def md(path,text):
 path=P/path;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text.strip()+'\n')

def main():
 r=json.loads((P/'EXPERIMENT_RESULTS.json').read_text());c=json.loads((P/'COST_SUMMARY.json').read_text());rows=list(csv.DictReader((P/'DIGIT_PER_NUMBER_RESULTS.csv').open()));manifest=json.loads((P/'data/DATA_GENERATED.json').read_text());v=r['visual']
 table='|数据|方法|全网格复MSE|IoU|SSIM|noise-ball|九图总wall(s)|RHS|全阶LU|\n|---|---|---:|---:|---:|---:|---:|---:|---:|\n'
 for ds in ['N','P']:
  for m,x in r['comparisons'][ds].items():table+=f"|{ds}|{m}|{x['complex_MSE']:.8f}|{x['IoU']:.4f}|{x['SSIM_full_grid']:.4f}|{x['noise_ball_pass']}/9|{x['total_wall_seconds']:.4f}|{x['RHS']:,}|{x['LU']:,}|\n"
 vt='|数据|数字|盲读|视觉0–5|IoU|GT孔洞→重建孔洞|完整Φ|\n|---|---:|---|---:|---:|---|---:|\n'
 for z in rows:
  if z['method']=='FULL256':vt+=f"|{z['dataset']}|{z['digit']}|{z['blind_read']}|{z['visual_score_0_to_5']}|{float(z['IoU']):.4f}|{z['GT_holes']}→{z['reconstructed_holes']}|{float(z['Phi']):.3f}|\n"
 rep='|数字|P拟合occupancy MSE|阈值IoU|GT latent范数|有限拟合是否收敛|\n|---:|---:|---:|---:|---|\n'
 for z in manifest['rows']:
  if z['dataset']=='P':rep+=f"|{z['digit']}|{z['projection_MSE']:.6f}|{z['projection_IoU']:.4f}|{z['latent_norm']:.2f}|{z['fit_success']}|\n"
 aP=r['comparisons']['P'];aN=r['comparisons']['N'];pSave=1-aP['AGENT_A']['total_wall_seconds']/aP['FULL256']['total_wall_seconds'];nSave=1-aN['AGENT_A']['total_wall_seconds']/aN['FULL256']['total_wall_seconds']
 gates=dict(overall='NUMERICAL_CAMPAIGN_COMPLETE_WITH_RESOURCE_PROTOCOL_FAILURE',resource_protocol='FAIL_CAP_EXCEEDED',new_physics_stopped=True,agent_A_theory='EXACT_STANDARD_SUPPORT_ELIMINATION',agent_A_kernel='SCOPED_NUMERICAL_BENEFIT',agent_A_whole_imaging='HOLD_NOT_UNIVERSAL_SPEEDUP',agent_B='NO_GO_FOR_TESTED_FIXED_POLYNOMIAL_BASES',agent_C_theory='GLOBAL_CONDITIONAL_BOUND_POSITIVE',agent_C_acceleration='NO_GO_FOR_TESTED_BASIS_AND_GATES',imaging_quality_increment='NOT_ESTABLISHED',controlled_digit_results='72_REAL_INVERSE_PREDICTIONS_WITH_LIMITATIONS',full66_posterior='NOT_ESTIMATED_UNQUALIFIED',TSOM_actual_inverse='NOT_RUN',AS_FTRC_fusion='NOT_RUN_NO_POSITIVE_COMPLEMENTARITY_AND_RESOURCE_STOP',paper='HOLD_NO_STRONG_TAP_CLAIM',COST=c)
 write(P/'RESULTS_SUMMARY.json',gates)
 md('RESEARCH_REPORT_ZH.md',f'''
# Research P3L-R11-S3：稳定 Maxwell 结构与数字 1–9 真正逆成像

本轮完成三条可运行结构算法、两套各九个数字对象，以及 **72 次真实非线性 Maxwell 材料反演**。所有预测先封存，再由主线程实际打开匿名 PNG 盲读，随后解封GT评分。最可靠的正面结果是精确活动支撑的物理核收益；没有证明新的材料成像能力，固定总场与全形状耗散控制没有在本轮形成有效加速。

**必须同时保留资源协议失败。** 自设全阶数LU硬上限为12,000，实际合计{c['LU_all_orders']:,}，超出{c['LU_excess']:,}。根因是只在分支/单对象执行RHS检查，未建立跨作业LU原子预留。发现后停止全部新增物理求解，未改写上限，也未把all-order LU偷换为LU256。保留数值证据不等于整轮预算PASS。RHS合计{c['RHS_conservative_all_orders']:,}/400,000，含144个新增576阶RHS；历史宿主264 RHS576另列。详见[COST_SUMMARY.json](COST_SUMMARY.json)、[失败记录](FAILURE_ANALYSIS.md)。

## 1. 先看实际图像

这是逐格材料反演，算法没有输入数字标签、字体模板或GT latent，没有枚举九个模板。显示仅最近邻放大24×24物理网格。

### N：物理合法的原生数字，新的258维反演模型

![N九数字实际逆成像](figures/DIGITS_SIDE_BY_SIDE.png)

### P：原P1 decoder可实现的近似数字，原Gaussianized66反演

![P九数字实际逆成像](figures/P_DIGITS_SIDE_BY_SIDE.png)

四个后端输出几乎相同。这首先是同一物理与反演目标的一致性，不是四种方法同时改善了图像。数字1和7的大笔画较清楚；4的交叉结构可见但孔洞丢失。数字8最明显失败：N的两个孔均未保留，P仅保留一个孔。6和9的闭环也未可靠恢复。主线程盲读把9记作“9/3”、8记作“0/8”，未在解封后修改。详见[视觉审查](DIGITS_VISUAL_REVIEW_ZH.md)。

## 2. 数据、逆问题与盲法

固定原P1 TM24/VIE、256格圆形ROI、8源×32复接收、原self-term与Green、host=0.1+0.04i、sigma=0.00119211798514519及source-major Re32/Im32白化。先用原576阶完整电流方程生成每份truth响应，并与精确ROI合同复核；再加原圆复Gaussian噪声。所有方法共享同一y。

N直接使用字体的格内占比，物性固定生成于1.3+0.105i。反演是占比256+Gaussianized属性2的box约束新模型，固定TV/IRLS、ridge和属性惩罚。P使用原Gaussian64+quintic和原属性变换，离线拟合只构造合法truth；在线初始化x=0，不读取拟合latent。原Gaussianized点目标为Phi+0.5||x||²，不是native uniform坐标的点目标。

P的九个拟合均达到预注册图形误差门，但GT latent范数47.85–94.48，典型Gaussian64半径约8。**它们是可表达而极不典型的设计对象，不是九个独立prior-predictive P1样本。** 拟合达到1200步上限且未全部满足优化收敛；实际truth仍严格取decoder输出，不取目标字体。

{rep}

所有方法按各自N/P合同共用初始化、正则、阻尼GN/IRLS、线搜索及80-J/240-forward上限。N九图均到80-J上限，不能称已充分收敛的全局最优反演；P按预定objective stopping结束，也不证明全局最优。两套数据不能互相当成纯求解器优劣比较。truth物性恰为共同初始化的区间中点，是有利的设计条件，未知物性鲁棒性未建立。

## 3. A：精确支撑消元实际解决了什么

只删严格a_i=0的电流状态行，解|S|阶系统；全部源共用一次分解。**没有删去这些格点的材料导数列。** 非活动格点的导数包含直达和重散射，所以原生占比可出生新支持。保留全部256列后，正演与完整Jacobian仍是原离散Maxwell的精确计算。

数字truth的活动状态数为71–112，相对256降维56.25%–72.27%。同一静态数字核：forward节时84.92%，完整J核节时63.54%；这是36个配对节点的有限计时窗口，不是整套逆成像加速。R10不同稠密度形状的初始配对窗口forward节时38.80%、完整J节时19.52%；后续优化代码的独立窗口为27.83%、15.65%，不可跨窗口挑最优比值。

真正逆问题中支持会变密。N终点平均141.44个活动格点，P159.56。A在N九图整个作业wall为12.1738s，完整LU为9.5617s，观察到增加27.32%；P为2.2544s对2.8884s，观察到减少{100*pSave:.2f}%。窗口未按新独立重复作严格性能确认，且存在并行资源竞争，不能宣称普遍20%全系统优势。数值图像没有改善。P平台内的专用有限幅度support-birth算法 **NOT_RUN**；这里的支持变化来自普通原生projected GN完整导数，不冒称新反演原理。

## 4. B：固定总场克服了哪个障碍，又为什么失败

总场输入E_R固定，确实移除了R10 weighted输入sqrt(A)E_R在零支持激活时必变的障碍。但当前CA仍随形状变化，固定注入不推出共同低维可达空间。全部新形状仍支付CAV、OAV、投影与完整residual成本。

在原R10的24形状×8节点中，固定anchor总场rank32/64/96/128分别只有0/1/5/14个通过原双门；固定Green多项式rank32/64/96全部0/192。旧weighted64为0/192，每形状重建weighted64为184/192，余8个为实际rank61未执行小核。anchor128有119个实际白化误差≤0.01，但只有14个通过完整bound+potential门；仍有73个实际超限，所以既存在截断误差，也存在上界保守性，不能全归结为其中之一。

主数字反演B的N共1317次、P共383次forward trial全部实际回退。需要Jacobian时直接用精确一致F/J，因为响应误差门不自动控制J。没有免费省掉32个导数RHS。真正的current/SOM/TSOM联合材料反演 **NOT_RUN**，本轮不能以总场Galerkin失败否定TSOM。

## 5. C：全形状界存在，但没有自动形成加速

原固定C的Hermitian虚部最小本征值=-4.5739299047e-5；beta*=同值。原幅度矩形给b*=0.009210948658>0。因此，对全部0≤a≤1与原物性矩形，标准耗散强制性推导给出统一resolvent范数界。在实际实现中是float64条件数值界，没有区间向外舍入或连续PDE资格。

用完整总场residual和端口范数计算white≤0.01、potential≤0.1门。不通过有限增阶再精确回退。同192节点、同固定Green基中，全局与局部两臂都0/192通过。全局边界setup约0.00548s、逐形状局部约0.19036s，减少了谱分析，却没有减少最终精确LU。128阶实际误差也全部超0.01，因此该组失败不只是全局常数太松。数字C同样所有1700个forward trial回退，成本上升。

## 6. 逐对象质量与真实费用

{table}

RHS逐阶数另见JSON；保守合计RHS不是不同矩阵阶数的等价FLOP。每个失败的小LU和回退均入账。B/C原始physics_kernel_wall曾重复计入current projection；原始回执保留，表中corrected_kernel_wall由父级完成事件重算，整个作业wall、RHS、LU不变。

每个dataset的noise-ball均为8/9。N的6、P的7未通过。noise-ball不等于恢复细笔画/闭环，正确truth对随机噪声也不保证落入固定球。本轮许多模糊图像满足数据球，因此只能确认当前材料、观测与优化目标的有限结果，不能从模糊图像证明物理不可辨识或唯一性失败。

{vt}

GT与重建的原目标值另见[离线诊断CSV](audit/PRIOR_AND_CONVERGENCE_DIAGNOSTICS.csv)。P的高范数truth受到很大的原prior点惩罚，重建通常用更平滑的低范数材料解释数据；它不是decoder无法表达数字的证据。N仍存在80-J截断、TV平滑、全局物性—占比补偿与有限数据的共同影响，尚未把原因完全分离。

## 7. 是否合并A+B+C

B+C已经作为C-controlled实际实现并运行九图；没有费用收益。A+B与A+B+C **NOT_RUN**：固定基没有可重复的增量，且发现资源越限后禁止继续增加物理。没有将三条线路强行命名为成功AS-FTRC。

保留A作为可选精确物理核、C作为可计算的数值误差门；停止扩大本轮固定多项式总场基和当前全局门驱动的组合。理论上还能检验小规模受控降阶导数或真正TSOM变量投影，但它们在本轮是未执行工作，不能用未来可能性救本轮NO-GO。

## 8. 原创性与下一项最小问题

精确支撑消元、块导数、Krylov/Galerkin、TSOM投影、耗散范数界均有标准背景。本轮新证据是：保留非活动材料导数的具体支持可变实现、同合同九数字反演、输入障碍消除后仍不足的总场负例，以及全形状界有效却未能加速的分项裁决。它们值得作为方法学和实现边界保留，尚不足以支持强IEEE TAP独立算法主张。

[作者MR-FFT-TSOM全文](https://renchengsong.github.io/papers/TGRS2329032.pdf)已包含电流子空间、内部传播筛选与数据/状态残差；[Forestiere–Miano](https://arxiv.org/abs/1611.04501)已有固定形状材料无关模式和幅度有理响应。不能把“固定物理载体”或“形状材料分离”本身重新宣称首创。有限ScholarQA检索有429缺口，指定三篇身份已核验；没有全球查新保证。

下一轮唯一优先问题应是：在明确定义的更小预算内，真正的共享材料TSOM/变量投影是否能减少Full-J反演费用并保留闭环，而非继续扩大已经失败的固定Green多项式基。先修复全局原子费用预留，再开新物理实验；本轮不自动启动下一轮。

## 9. 可追溯交付

- [开始阅读](START_HERE.md)、[科学裁决](GO_NO_GO.md)、[统一数字合同](UNIFIED_DIGIT_BENCHMARK.md)
- [A详细报告](agent_a/REPORT_ZH.md)、[B详细报告](agent_b/REPORT_ZH.md)、[C详细报告](agent_c/REPORT_ZH.md)
- [逐数字CSV](DIGIT_PER_NUMBER_RESULTS.csv)、[方法CSV](METHOD_COMPARISON.csv)、[结果JSON](EXPERIMENT_RESULTS.json)
- [盲读与GT后视觉审查](DIGITS_VISUAL_REVIEW_ZH.md)、[可复现入口](REPRODUCE.md)、[源码与输入manifest](SOURCE_MANIFEST.json)
- [12条新推导](../maxwell_formula_atlas/11_ACTIVE_SUPPORT_TOTAL_FIELD_AND_GLOBAL_BOUNDS.md)，累计158条可追溯公式。

所有原R10/R9报告与原始账本只读。R11最终状态为数值实验已完成、资源协议FAIL、成像新增收益未建立；这些资格不得合并成单一PASS。
''')
 md('START_HERE.md','''
# R11-S3 从这里开始

已完成N/P各九个数字、四个物理后端共72次真实逆散射，不是模板分类。主线程已先盲读图像后解封GT。

**整轮资源协议FAIL_CAP_EXCEEDED：28,727全阶LU超过12,000。已停止新增物理，不回写预算。** 数值结果与失败仍完整保留。

1. [中文完整报告](RESEARCH_REPORT_ZH.md)，含两套九数字对照图。
2. [逐数字盲视觉审查](DIGITS_VISUAL_REVIEW_ZH.md)，特别检查8的两个孔。
3. [A精确支撑](agent_a/REPORT_ZH.md)：稀疏数字完整J核节时63.54%，不等于全成像收益。
4. [B固定总场](agent_b/REPORT_ZH.md)、[C统一损耗界](agent_c/REPORT_ZH.md)：可运行但未加速。
5. [方法CSV](METHOD_COMPARISON.csv)、[逐对象CSV](DIGIT_PER_NUMBER_RESULTS.csv)、[资源汇总](COST_SUMMARY.json)。
6. [复现说明](REPRODUCE.md)、[科学决策](GO_NO_GO.md)、[下一轮交接](GPT_PRO_HANDOFF/START_HERE.md)。

N为新258维占比模型；P为原66维decoder中高范数的设计真值。两者不能互称独立P1统计泛化。Full66 posterior与精确材料纤维均未认证。
''')
 for root,child in [('AGENT_A_EXACT_SUPPORT.md','agent_a/REPORT_ZH.md'),('AGENT_B_FIXED_TOTAL_FIELD.md','agent_b/REPORT_ZH.md'),('AGENT_C_GLOBAL_DISSIPATIVE_BOUND.md','agent_c/REPORT_ZH.md')]:md(root,f'# {root.removesuffix(".md")}\n\n完整可运行代码、结果与范围见[{child}]({child})。\n\n整轮资源状态为FAIL_CAP_EXCEEDED，不以分支数值结果冲抵预算失败。')
 md('UNIFIED_DIGIT_BENCHMARK.md',f'''
# 统一数字基准

冻结原P1物理与噪声。统一DejaVuSans-Bold字体112子像素，最大格内inkbox为9×13物理cells，8×8面积采样只构造占比，实际物理网格始终24×24。圆ROI不裁掉关键笔画。字体实际BLAKE2b、数据源和噪声seed在data/GENERATION_CONFIG.json与DATA_GENERATED.json。

N: 新占比256+Gaussianized属性2模型；固定TV/IRLS。P: 原Gaussianized66，原quintic和Gaussian点惩罚。共同材料幅度由原区间生成，truth取1.3+0.105i。该值等于统一属性初始化中点，这是有利设计条件，未知属性初值鲁棒性未检验。

在线只读data/online/*的y；不读truth、字体、数字标签、kernel_digit_shapes、拟合GTlatent。所有模型相同y、参数化内相同正则/初值/80J240F/停止。参考求解器FULL256为同一个非线性Maxwell。N达到迭代cap而非收敛证明。

物理truth由完整576阶方程生成；最终72个材料均独立ExactSolver检查。noise-ball为Phi≤282.8737921622988。候选图像不能因该门PASS就称拓扑恢复。固定二值阈值0.5，前景4连通/背景8连通，在原24网格统计。

{table}
''')
 md('DIGIT_REPRESENTABILITY_AUDIT.md',f'''
# P1 decoder数字表达能力

先拟合占比再用decoder输出作为物理truth，不把目标字体当作P真值。门在拟合前固定为ROI occupancyMSE≤0.025且阈值IoU≥0.75。全部九个达到该有限图形门，3与5误差明显较大；有限L-BFGS拟合并未全部满足优化收敛。

{rep}

范数48–94远高于Gaussian64典型半径约8。合法先验支持不等于高先验概率。P反演仍使用原Gaussianized点目标，低范数平滑解可比高范数truth有更低原目标；不能从P图像模糊推断decoder表达失败，更不能宣称九个典型P1随机parent。

![目标N数字](figures/DIGITS_GT_1_TO_9.png)
![真实P decoder数字](figures/P_DIGITS_GT_1_TO_9.png)
''')
 md('DIGITS_VISUAL_REVIEW_ZH.md',f'''
# 九数字实际视觉审查

主线程实际打开全部8张匿名后端九宫格PNG。先把A–I的盲读、0–5分和理由写入audit/BLIND_READ_N.json与BLIND_READ_P.json，随后才读BLIND_KEY_SEALED.json和GT。四种后端的图像无可见差别；盲读是主线程主观评价，不是独立人群识别实验。

{vt}

1和7较可靠。2/3/5的连接和弯折较模糊；4交叉可辨但内部孔丢失；6/9闭环未保留。8在N为0/2孔、P为1/2孔，均没有完整恢复。N9与P9盲读都记为9/3，8记为0/8，不能事后改成全部识别正确。

全部算法相同质量，不能凭A更快的小核将其视觉评分抬高。N与P的真值不同，二者SSIM或MSE不能作为纯算法因果比较。

![N对照](figures/DIGITS_SIDE_BY_SIDE.png)
![P对照](figures/P_DIGITS_SIDE_BY_SIDE.png)
![N真实材料误差](figures/DIGITS_RECONSTRUCTION_ERRORS.png)
![P真实材料误差](figures/P_DIGITS_RECONSTRUCTION_ERRORS.png)

显示是原24网格最近邻放大，无平滑显示或高分辨率重建伪装。mask、孔洞和连通数均在原网格计算。逐像素误差与拓扑不受视觉主观评分替代。
''')
 md('GO_NO_GO.md',f'''
# 分项科学裁决

|对象|裁决|依据|
|---|---|---|
|整轮资源协议|FAIL_CAP_EXCEEDED|28,727 LU超过12,000；发现后停止，未重释预算|
|A精确支撑|SCOPED_KERNEL_GO / WHOLE_IMAGING_HOLD|71–112阶数字truth核完整J节时63.54%；N全流程变慢27.32%，P观察节时{100*pSave:.2f}%且未做隔离重复确认|
|A专用latent支持出生|NOT_RUN|普通native projectedGN支持变化已运行；没有新有限latent birth算法证据|
|B固定总场|NO_GO_FOR_TESTED_BASES|anchor128仅14/192，数字forward trials全回退|
|C全形状界|CONDITIONAL_MATH_PASS|b*=0.009210948658>0，有限浮点条件界|
|C误差控制加速|NO_GO|全局/局部192节点均全回退；数字也全回退|
|新成像能力|NOT_ESTABLISHED|四后端相同材料输出；孔洞明显缺失|
|A+B/A+B+C|NOT_RUN|未观察互补增量且物理已STOP；不强行拼框架|
|B+C|RUN_NO_BENEFIT|C-controlled solver已集成运行|
|真正TSOM材料-电流优化|NOT_RUN|不得由总场Galerkin失败外推|
|Full66 posterior / exactfiber|UNQUALIFIED|未估计，不是本轮结果|
|强IEEE TAP独立贡献|HOLD|标准结构实现与有范围核收益，尚无强反演/泛化/同质量总成本证据|

下一项仅建议小预算共享材料TSOM/variable projection对普通Full-J，配原子全局费用预留。本轮不自动启动新物理。
''')
 md('FAILURE_ANALYSIS.md','''
# 失败、纠错与未执行记录

1. **FAIL_CAP_EXCEEDED。** 跨阶LU上限12,000，实际28,727。局部账本只保护RHS；多个对象与小阶增阶把LU累积漏过全局约束。72个预测已生成时由Sol-B审计发现，Main即时叫停三个分支。原始账本和初始合同保留，未将allorders改成LU256，未扩预算，也未继续物理计时重放。费用不合规是整轮限制，不是可以因数值好看而忽略的小问题。
2. **主线程时间标签错误。** 初始EXPERIMENT_CONTRACT手写20:07时间标签晚于真实文件mtime与首个inverse。审计后保留原字节于audit/INITIAL_CONTRACT_FUTURE_CLOCK_LABEL.json，修正为原文件实际mtime；阈值、方法和数据不变。INVERSE_CONFIG实际先于首个inverse。不能靠修正隐藏原错误。
3. **B/C嵌套计时重复。** CURRENT_SHAPE_PROJECTED已包含在EVALUATE_COMPLETE，原physics_kernel_wall重复计。原回执不覆盖；派生corrected_kernel_wall只求父级事件+最终精确验证。全wall、RHS、LU和数值输出不变。
4. **评分依赖缺失。** 初次评分因未安装skimage在读GT前导入失败；保留SCORING_DEPENDENCY_FAILURE.log，使用明确7×7 uniform-window SSIM公式和scipy滤波，未改二值阈值、主MSE或模型。无新增物理，无训练。
5. **Sol-C瞬时容量失败。** 首次Sol-C agent返回modelcapacity错误，重试同一个Sol最高强度后完成；没有Astra或DeepSeek替代冒名。
6. **固定总场失败。** 注入固定仅移除weighted支持障碍，CA不保持固定Krylov空间。anchor128仍73/192真实超0.01；Green128在C组全部真实超门。并非仅上界保守。
7. **全局耗散界未加速。** b*>0且便宜，但真residual不够小、端口/势门保守和精确J成本共同消除收益；不把理论界PASS改写成计算PASS。
8. **A核与整流程分离。** 稀疏truth上强正面核结果不能代替未知支持反演；N全流程较慢，P只单窗口观察节时。没有全系统普遍加速主张。
9. **图像与先验。** N80J未充分收敛；Ptruthlatent极不典型；unknownproperty初始化恰在truth中点。模糊与丢孔不能单因归结为物理不可辨识。9个设计对象不足大样本泛化。
10. **NOT_RUN。** 专用P1有限支持出生、真正TSOM电流材料联合反演、受控低阶Jacobian、A+B、A+B+C、大样本确认和Full66 posterior。不得以公式或设计替代结果。
''')
 md('PAPER_ORIGINALITY_AND_NARRATIVE.md','''
# 论文资格与可复用叙事

可保留的主线是“支持可变的精确离散Maxwell计算，完整非活动材料导数，以及跨形状降阶失效的可审计边界”。目前它仍是受控计算研究，不是新的概率成像突破。

候选标题：**支持可变非线性Maxwell反演的精确消元与残差受控降阶：数字成像的能力和边界**。

可以证实：严格零支持状态可以精确消去，非活动材料导数仍非零；数字占据稀疏时物理核明显节时；固定总场移除加权注入缺陷但没有自动跨形状加速；全形状耗散界存在而近似误差仍使其频繁回退。

不能证实：新解族、准确Full66 posterior、普遍20%整成像加速、比传统TSOM更强的反演、稳定高分辨率恢复数字闭环、原P1随机总体泛化或强IEEE TAP原创机制。资源LU越限更要求降级整轮实验资格。

TSOM已有固定电流与内部传播投影；作者公开全文提供直接依据：[Xu等2015](https://renchengsong.github.io/papers/TGRS2329032.pdf)。固定形状材料无关模态及幅度有理系数也有[Forestiere–Miano2016](https://arxiv.org/abs/1611.04501)先例。Schur、支撑消元、Galerkin与耗散强制性均不是本轮新定理。

ScholarQA两次关键词发现有429与弱相关记录，不能声称穷尽查新；三个已知DOI身份全部解析。资料在theory/，实质支持采用作者全文而非仅题名。

若继续，先用小预算真正共享材料TSOM/variable projection，对照同目标Full-J。以完整成本与闭环恢复裁决，不能再用“固定空间”这一已有概念本身当贡献。
''')
 md('GPT_PRO_HANDOFF/START_HERE.md','''
# R11给下一轮的交接

先读[中文报告](../RESEARCH_REPORT_ZH.md)、[GO/NO-GO](../GO_NO_GO.md)与[失败](../FAILURE_ANALYSIS.md)。本轮三位原生Sol/ultra实际并行，72次两模型九数字真inverse，不是模板输出；预测先封存，主线程实际盲读后GT评分。

资源协议FAIL：全阶LU28,727超过12,000，所有新物理已停止。RHS363,240/400,000；不能用RHS合规消除LU失败。

最有价值正面是A的严格支撑核与保留inactivecolumns；最重要负面是固定totalfield仍全场景不足，globalb*>0也没有实测加速。C并未提供新测量信息，四后端图像同质。

N为新258模型、P为原66decoder高范数设计真值，不是典型P1parent。N到80J上限；P低范数平滑解更受prior青睐。8双孔未恢复，不能称九图完整成功。

下一项唯一优先物理问题：真正共享材料TSOM/variable projection是否能比普通Full-J减少费用并保留拓扑？当前NOT_RUN。先修好全局费用原子预留，再另开预算，不继续扩大本轮失败Green多项式基。
''')
 md('GPT_PRO_HANDOFF/NEXT_THEORY_QUESTIONS.md','''
# 由真实结果提出的问题

1. 如何让状态支持消元保留零占比格点的可激活材料导数，同时避免大量inactive列消耗掉稀疏LU节省？输出任务/JVP形式是否比完整J更合适？
2. 固定总场输入已包含所有照明，但CA迭代不断产生新方向。能否用少量有理/经验采样而非高阶多项式基，在同样材料支持内满足严格端口误差？需要可反证的有限族条件。
3. 全局b*>0而拒绝率仍100%。何种低成本双边/dual加权残差能给任务相关界，减少不影响32端口的全状态误差惩罚？不能用经验误差改名证书。
4. 真正的TSOM有限罚变量投影不必每步精确Full-J；能否在原共享材料decoder内取得同质量总成本收益？本轮没有执行该实验。
5. 数字高范数truth受原prior点惩罚，如何区分材料支持、典型概率与设计成像目标，而不为了数字效果悄悄换prior？
6. 为什么拟合noise-ball仍丢闭环？需分开80J未收敛、TV/MAP偏好、占比—幅度补偿和物理分辨率；当前没有唯一性/不可辨识性证书。
''')
 # Portable evidence inputs copied losslessly, no Maxwell execution.
 shutil.copy2(R9/'results/stable_roi/SETUP.npz',P/'inputs/FIXED_HOST_SETUP.npz')
 bg=np.load(R4/'research/p3l_r4_conditional_transport/inputs/background/weak.npz');np.savez_compressed(P/'inputs/PRIOR_POINTS.npz',points=bg['points'])
 shutil.copy2(R4/'research/p3l_r4_conditional_transport/configs/TARGET.json',P/'inputs/PRIOR_CONFIG.json')
 shutil.copy2(R4/'src/p1_material/prior.py',P/'inputs/PRIOR_SOURCE.py')
 write(P/'audit/REPORT_ASSEMBLY.json',dict(utc=time.time(),new_physics_RHS=0,scientific_status=gates,root_source_before_R11='d6120d4bc0f9d8699b9e788da7abee345752b8f6'))
 print('REPORTS_ASSEMBLED',len(list(P.glob('*.md'))))
if __name__=='__main__':main()
