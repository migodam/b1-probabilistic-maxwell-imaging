from pathlib import Path
import csv,json,time,hashlib
import numpy as np
from accounting import P,write

def load(s):return json.loads((P/s).read_text())
def md(name,text):
    f=P/name;f.parent.mkdir(parents=True,exist_ok=True);f.write_text(text.strip()+'\n')
def table(rows,headers):
    return '|'+ '|'.join(headers)+'|\n|'+'|'.join(['---']*len(headers))+'|\n'+'\n'.join('|'+ '|'.join(str(v) for v in r)+'|' for r in rows)
def csvfile(name,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (P/name).open('w',newline='') as f:
        w=csv.DictWriter(f,keys);w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in r.items()})

def main():
    resource=load('RESOURCE_SUMMARY.json');primary=[];physical=[];coverage=[];extent=[]
    for parent in range(3):
        for method in ('A','B'):
            s=load(f'results/E1/P{parent}/{method}/SUMMARY.json');primary.append(s)
            for row in s['rows']:
                coverage.append(dict(stage='E1',**row))
                f=P/f'results/E1/P{parent}/{method}'/row['key']/'ENDPOINT.npz'
                if f.exists():
                    with np.load(f) as z:
                        m=z['m'];physical.append(dict(stage='E1',parent=parent,method=method,key=row['key'],material=m.tolist(),**{k:row.get(k) for k in ['qualified','legal','material_RMS','pair_relative','pair_whitened','Phi','state_relative','sigma_min','actual_RHS']}))
            if method=='A':
                good=[r for r in s['rows'] if r['qualified']]
                mats=[]
                for r in good:
                    with np.load(P/f'results/E1/P{parent}/A'/r['key']/'ENDPOINT.npz') as z:m=z['m'].copy()
                    if not any(np.linalg.norm(m-a)/np.sqrt(6)<1e-8 for a in mats):mats.append(m)
                extent.append(dict(parent=parent,sampled_qualified_nodes=len(good),numerically_distinct_materials_including_reference=len(mats),maximum_material_RMS=max(r['material_RMS'] for r in good),maximum_pair_relative=max(r['pair_relative'] for r in good)))
    costs=load('results/COST_CONTROL/SUMMARY.json');g=costs['grouped']
    totals={m:{key:sum(r[key] for r in g if r['method']==m) for key in ['RHS','total_seconds','solver_plus_validation_seconds','qualified','cases']} for m in ['DYNAMIC_FULLJ','FROZEN_J','CURRENT_FULL_KERNEL']}
    dyn,cur=totals['DYNAMIC_FULLJ'],totals['CURRENT_FULL_KERNEL']
    saving=1-cur['total_seconds']/dyn['total_seconds'];kernel=1-cur['solver_plus_validation_seconds']/dyn['solver_plus_validation_seconds']
    E2=load('results/E2/SUMMARY.json');local=load('results/E4_local/SUMMARY.json');legacy=load('audit/E4_legacy_replay/e4_two_cell_posterior.json')
    interval=load('audit/INTERVAL_CHART_CERTIFICATES.json');rank=load('audit/ENDPOINT_RANK_REVIEW.json');E3=load('results/E3/SUMMARY.json')
    closure=load('CLOSURE.json');res=dict(campaign='Z-ALG-E1',status='COMPLETED_AND_CLOSED',
        main_decision='NO_GO_FOR_CURRENT_GENERAL_ALGEBRAIC_ADVANTAGE',
        gates=dict(ASYMMETRIC_LOCAL_MATERIAL_FIBER='SCOPED_VALIDATION',CURRENT_EXTRA_COVERAGE='NOT_ESTABLISHED',CURRENT_TOTAL_COST_20_PERCENT='FAIL',
            REMOTE_DISCRETE_SECOND_BRANCH='NOT_FOUND_WITHIN_THIS_SEARCH',LOCAL_CONDITIONAL_MATERIAL_PROBABILITY='NUMERICAL_PASS_LOCAL_ONLY',GLOBAL_MISSING_MASS='NO_USEFUL_BOUND',
            CROSS_DECODER='COORDINATE_TRANSFORMATION_ONLY',FULL66_EXTENSION='DO_NOT_EXPAND',TAP_NOVELTY='NOT_ESTABLISHED'),
        public_reference_contexts=3,new_model_dimensions=dict(E1=[6,4,2],E2=[4,4]),
        extent=extent,primary_coverage=[{k:s[k] for k in ['parent','method','attempts','qualified','grid_qualified','total_actual_RHS']} for s in primary],
        optimized_cost=dict(totals=totals,current_total_saving=saving,current_solver_validation_saving=kernel,current_actual_RHS_saving=1-cur['RHS']/dyn['RHS'],
            gate=0.20,total_gate_pass=saving>=.20,shared_reference_RHS=27,whole_controlled_phase='RESOURCE_SUMMARY.jobs contains full RUN_END clocks; pooled method windows exclude unassigned common setup/report overhead'),
        E2=dict(starts=len(E2['starts']),qualified=sum(r['qualified'] for r in E2['starts']),found_roots=E2['root_count_found'],completed_retraced_cycles=sum(r['completed'] for r in E2['cycles']),
            maximum_return_RMS=max(r['return_RMS'] for r in E2['cycles']),coverage='No global uniqueness or all-branch certificate'),
        probability=dict(legacy_new32_event=legacy['results']['32']['original_box']['event_x1_gt_x2'],exact_pair_branch_probability=0.75,
            legacy_wrong_uniform_event=.5,fresh32_reused64_explicit=True,local_status=local['status'],local_convergence=local['convergence'],local_cross_chart=local['cross_chart'],
            local_n8_A=local['orders']['8']['A'],local_n8_B=local['orders']['8']['B'],local_n8_wrong_uniform_free=local['orders']['8']['wrong_uniform_free'],global_mass='UNKNOWN'),
        interval_certified_free_halfwidths=[r['attempts'][-1]['halfwidth'] for r in interval],
        posthoc_endpoint_rank=dict(points=rank['points'],min_response_sigma=min(r['response_sigma_min'] for r in rank['rows']),max_dependent_condition=max(r['dependent_condition'] for r in rank['rows']),qualification='Posthoc at saved points; not an online or between-point rank guarantee'),
        resource={k:v for k,v in resource.items() if k!='jobs'},qualifications=['Deterministic reference materials, not independent held-out prior parents','Small disk-pulse collocation model, no continuum PDE certificate','No Full66 posterior or Bayes-risk calibration','No new neural training'],
        initial_source_commit='f543d6cc',post_E0_fix_commit='9859af44',reproducibility='SOURCE_MANIFEST plus FINAL_SOURCE_MANIFEST; metadata corrections preserved')
    write(P/'RESULTS_SUMMARY.json',res)
    csvfile('PHYSICAL_EVIDENCE.csv',physical)
    coverage.extend(dict(stage='E2_multistart',**r) for r in E2['starts'])
    csvfile('BRANCH_AND_COVERAGE.csv',coverage)
    probs=[]
    for n in ('32','64'):
        for method in ('original_box','ordered_quotient'):
            z=legacy['results'][n][method];probs.append(dict(model='legacy_symmetric2',order=n,method=method,source='NEW_REPLAY' if n=='32' else 'IMMUTABLE_SUPPLIED_REFERENCE_NOT_RERUN',**{k:z[k] for k in ['scaled_evidence','mean_x1','mean_x2','event_x1_gt_x2','psi_x1_squared']}))
    for n,d in local['orders'].items():
        for method in ('A','B','wrong_uniform_free','wrong_equal_nodes'):probs.append(dict(model='new_local6_fiber_tile',order=n,method=method,source='NEW_LOCAL_ONLY',**d[method]))
    csvfile('PROBABILITY_EVIDENCE.csv',probs)
    costtable=table([[r['method'],r['RHS'],f"{r['total_seconds']:.6f}",f"{r['solver_plus_validation_seconds']:.6f}",f"{r['qualified']}/{r['cases']}"] for r in totals.values()] if False else
        [[m,v['RHS'],f"{v['total_seconds']:.6f}",f"{v['solver_plus_validation_seconds']:.6f}",f"{v['qualified']}/{v['cases']}"] for m,v in totals.items()],['方法','实际RHS','含轨迹保存的方法窗口/秒','求解＋验收＋方法setup/秒','合格'])
    mainrows=[]
    for i in range(3):
        a=next(s for s in primary if s['parent']==i and s['method']=='A');b=next(s for s in primary if s['parent']==i and s['method']=='B')
        mainrows.append([i,f"{a['qualified']}/{a['attempts']}",f"{b['qualified']}/{b['attempts']}",a['total_actual_RHS'],b['total_actual_RHS'],f"{extent[i]['maximum_material_RMS']:.4f}"])
    main_table=table(mainrows,['参考材料','直接材料A','电流B','A RHS','B RHS','最大材料RMS'])
    md('GO_NO_GO.md',f'''# Z-ALG-E1 分项裁决

**主裁决：NO-GO，针对“当前完整电流图册提供普通材料方法没有的独立覆盖或合格总成本优势”这一主张。**

非对称6单元的合法局部材料纤维与条件数值积分是本轮真实正面结果。它们可作为有条件的小模型基准继续使用；不把它们改写成 current 方法获胜。

|对象|裁决|理由|
|---|---|---|
|非对称正维材料图册|限定验证通过|6实材料、4实观测，真实Green；小区间收缩证书仅针对储存的binary64系数模型|
|current额外合法覆盖|未建立|三参考材料的固定网格与径向节点资格完全相同，配对材料RMS约1e-12|
|current总成本20%门|FAIL|统一完整轨迹后的方法窗口仅节省{saving*100:.3f}%；一例变慢|
|current求解＋验收子项|限定正面|节省{kernel*100:.3f}%，不可替代总费用门|
|非对称离散第二支|预算内未发现|33启动中31个合格并聚为1根；另两次失败保留；不是唯一性证明|
|局部coarea材料概率|数值局部PASS|材料图册与current隐式切向权重一致；全域遗漏质量无非平凡界|
|跨decoder|只完成换基|完整DFT/Haar6不改变材料与prior；截断改变支持|
|原P1 Full66|本轮不扩大|新支持和新采集；没有得到Full66新纤维或posterior证据|
|强IEEE TAP独立贡献|未建立|主要构件已有CSI/SOM、IFT和coarea先例；没有新的覆盖/总成本增量|

不因“没找到第二根”宣称不存在多解，不因“当前电流法不划算”否定全部材料几何，也不以重复次数代替独立parent。
''')
    md('RESEARCH_REPORT_ZH.md',f'''# Research Z-ALG-E1：非对称 Maxwell 材料图册、条件积分与成本反证

2026-10-29。完整实验已封账。实际支付 **{resource['actual_complex_RHS']:,} complex RHS**，另有闭式物理响应 **{resource['closed_form_source_RHS_equivalent']:,} source-RHS等价**；更保守合计 **{resource['conservative_combined']:,}/30,000**。原R6/R7只读，没有预算转移。

## 1. 核心答案

在没有镜面交换对称的六单元真实Hankel–Green模型中，确实构造了材料会实际改变、完整两照明复响应保持到严格浮点门内的局部二维材料图册。这个结果不只是潜变量gauge，也不是仅看一阶弱方向。直接材料IFT/Newton和受共享实材料约束的完整接收零空间电流法，都能构造这些节点。

**但本轮没有建立代数current图册独立胜过普通材料方法的总成本或覆盖收益。** 两法通过同样的147/147局部网格目标，径向节点验收也一致；相同自由坐标上的材料差最大约1.03e-12。统一完整轨迹保存并加入便宜的冻结Jacobian控制后，current方法窗口总耗时仅降低{saving*100:.3f}%，低于预注册20%门。求解＋验收子项的{kernel*100:.3f}%改善单独保留。

四单元满列秩子模型的33个启动中31个通过严格响应验收，合并为一个材料根；两次失败保留。两条各32节点的数据往返延拓回到原根，没有发现非对称远处第二支。有限搜索不是全域唯一性证明。

本轮裁决：**通用代数优势主张NO-GO；非对称局部材料图册与条件概率保留为限定正面基准；暂不进入原Full66。**

## 2. 附件输入、历史与协议差异

先完整读取38页Z-ALG报告及全部TeX章节、脚本、JSON、claim ledger和协议；PDF关键页另作视觉核对。输入ZIP共28条目，全部CRC检查并建立BLAKE2b清单。作者提供的SHA metadata原样保留，本轮不追加SHA-256检查。精简研究包没有R6/R7大原始数组；历史合同从固定B1源码/报告核对，不假装重新完成历史大实验。

原两单元材料(1,3)/(3,1)的交换正例、商(p,q)=(4,3)、dark-current反例、复三次monodromy和E4概率积分是**旧成果回归**，不是本轮新的非对称发现。E0重放实际20 RHS；之后仅报告序列化失败，从这些已支付矩阵/RHS/解逐字节缓存恢复，零新增物理。

默认旧E4同时重算32/64阶会有15,370次解析物理响应。本轮采用更保守预算：只新算32阶的两种积分，明确复用不可变64阶参考作为收敛回归，支付10实际RHS＋6,164解析等价。原报告的旧64阶时间不冒充本机新时间。新32阶原空间和商积分互相比较，同时各自对照交付64阶，原误差门没有放宽。

原P1仍是24×24、Full66、8×32复采集、原decoder与complex-noise约定。**以下新模型是新的P_Z材料支持、少量接收与real-component σ=0.001**，与P1统计合同不同。没有根据结果删通道、抬噪声或修改先验救正结果。

## 3. 冻结的非对称物理模型

新E1有6个半径0.04的圆盘常数pulse基，k=2π，exp(-iωt)，材料对比度为实数[0.5,4]，被动无损。内部Green包含真实圆盘积分自项和互项；前向为完整非线性VIE，而非Born。

两源为(-0.61,-0.72)、(0.63,-0.57)，单接收为(0.17,0.86)，标签固定，产生2个复数即4个实观测。三点不共线，不能由非恒等欧氏反射逐个固定；719个非恒等单元置换的完整G/Einc/S算子缺陷最小约0.559，未发现交换自同构。两个照明确实独立。

新原prior为独立实材料密度 p_Z(m)=∏[2m_i/(4²−0.5²)]，没有原子。材料RMS=||Δm||/√n。三份公开合法参考材料事前固定，用作**物理存在性构造**，不是从隐藏GT恢复的在线成像测试。它们也不是独立prior-predictive测试parents。噪声白化以每实分量0.001计，pair和Phi一直相对同一参考均值。

第一参考点完整材料Jacobian奇异值约(0.004576,0.003551,0.001906,0.001185)，数值满行秩4。全部6列在三档中心差分的最大相对误差9.26e-10；current方程Jacobian差分也通过。新非对称dark-current扰动的数据相对差约2.37e-16，但共享材料残差约0.03019，且恢复材料存在约0.0862虚部，不能是本问题的合法材料。

## 4. 两种真实图册怎样构造

A直接材料图册：以full material J作列主元QR，选4个受约束材料坐标、2个自由材料坐标，用带材料box的完整Newton和回溯求解固定完整响应。每个源forward RHS及接收伴随列都计费；全J仅额外1个复transpose-adjoint RHS，不人为弱化普通基线。

B完整current图册：每个源写j_s=j_s^0+Nc_s，S N=0；两个源共20个实current坐标，再加4个dependent实材料，联立24个实构成方程j_s−diag(m)(Einc_s+G j_s)=0。它和共享材料minors/实材料条件等价，但显式保留m避免场节点除零。G N预计算与所有代数Jacobian、实求解、回溯和保存均计入CPU/墙钟。不是给current系数另设平坦prior。

每个候选都用新的完整forward验收，不用拟合电流给自己作证。所有合格端点的数据、状态残差均不超过1e-10相对门，合法box和resolvent稳定性也检查。current构成残差最差约2.20e-14，最小两照明场平方和约0.01263，没有场节点。另以70位精度独立重算三个参考和三个事前指定径向点，最大相对差约2.01e-15；这不是连续PDE或Hankel区间认证。

## 5. 合法覆盖、有限位移与rank资格

{main_table}

上表是**采样目标节点资格**，不是互异根数，也不是posterior质量。网格和径向目标有重复；逐材料去重数及实际RMS范围在RESULTS_SUMMARY/CSV保留。A/B共748次目标尝试，728合格、738有完整独立验证数组，其余越界/停滞点同样保留。没有只挑漂亮端点。

对738个已保存端点，后来额外支付738个伴随RHS检查实际J；数值秩都为4，最小奇异值约0.00104123，固定dependent block最大条件数约7.403。这是**事后离散节点诊断**，不冒充原在线rank-stop，也不保证节点之间全部正则。没有确证的chart奇异事件，未宣称实现了跨奇异域或整个连通分量覆盖。

对三个中心另做Neumann逆包围＋带权收缩/自映射检验。得到自由坐标半宽0.00625、0.00625、0.0015625的小矩形：在储存的binary64 Green系数被视为精确常数的有理离散模型中，每个自由点都存在声明dependent box中的唯一根。这给出局部连续二维材料图册依据。真实超越Hankel系数没有严格区间包围，有限pulse离散也没有连续PDE误差证书。较大数值网格/径向范围不借小证书升级为全覆盖。

## 6. 为什么RHS优势没有变成合格总成本优势

原记录窗口中，A总计4,853 RHS，B747 RHS；受控job墙钟约3.846与2.179秒，表面节省约43.3%。但A付费物理中间态逐次落盘，B中间current代数态后来通过0-RHS精确重放完整补存；该窗口不能单独支撑公平总成本主张。

因此事前追加明确成本核验合同：三参考材料、各9个相同自由目标、3次重复，旋转方法顺序；全部真实试探与current代数状态使用同一完整buffer/逐端点批量无损保存策略；双方均复用精确中心，均支付新的完整端点验收。另加入普通冻结材料J基线。所有243个方法案例都合格。

{costtable}

上表是setup＋求解＋验收＋batch持久化的**逐方法窗口合计**；阶段整体RUN_END还包含共用参考和未归属循环/checkpoint开销，单列账本，不称上表为整个阶段wall。共同参考另实际支付27 RHS。给双方添加同样非负共同开销只会使13.6%总节时更小，不会救过20%门。RHS优势是真实的，但小型6×6 LU很便宜，current的24维实校正、代数状态构建和轨迹体量会抵消它；此处不能外推大型Maxwell加速倍数。

## 7. 远处分支：独立问题，结果为未发现

E2另取前4个非对称单元，保留原两源与接收，n=k=4。公开参考(1,3,2,1.4)。一个参考启动加32个冻结seed均匀启动，完整Full-J/边界回溯：31/33严格合格，合并为1个根；两次未达标原样保存。局部满秩和这次未发现不能证明全域唯一。

两条预定正弦数据路径各32节点，用普通continuation跟踪后返回；最大材料回返RMS为{res['E2']['maximum_return_RMS']:.3g}。它们是沿线段往返，不是非平凡monodromy。循环节点有新的完整forward/J和保存数组，但不是另付numpy独立重解；独立审计用保存矩阵乘法复核，不虚构新求解收据。有限覆盖/判别式前提未建立，复杂代数monodromy标NOT_RUN。

## 8. 概率：旧全域小模型与新局部核分开

旧两单元精确商值(p,q)=(4,3)的条件分支权重仍是1/4、3/4，事件概率0.75。有限噪声、完整四复通道下的新32阶回归事件为{res['probability']['legacy_new32_event']:.10f}，与原64阶一致；错误均匀branch仍为0.5。两者不是同一条件对象。

新六单元只报告第一个参考的自由坐标半宽0.05图册片，在同一P_Z上的coarea条件核。直接材料权重为 p_Z/|det J_dep|；current隐式恢复材料切向T后，权重为 p_Z sqrt(det(TᵀT))/∏s_j(J)。不是current均匀先验，不是候选likelihood softmax。使用4、6、8阶，物理材料阈值事件在两半矩形分别积分。

末两阶质量相对差约{local['convergence']['mass_relative']:.2g}，最大材料均值差约{local['convergence']['mean_max']:.2g}；两种图册质量差约{local['cross_chart']['mass_relative']:.2g}。这些是数值收敛与独立恢复权重的一致性，没有严格积分余项界。0.05的数值积分片大于0.00625的小收缩认证域；不能用后者保证前者整个图册或其它分支。

**全材料域的未覆盖条件质量没有非平凡上界；没有计算完整带噪posterior。** 等权节点、均匀自由坐标及打乱条件耦合的控制单独保存。不能把本局部核接到原Full66后就声称posterior恢复。

## 9. 跨decoder与成熟先例

完整原生坐标、满足实共轭条件的完整DFT、正交不平衡Haar6，在相同已保存材料上往返通过，原prior按同一换基推前。截断3系数有非零材料误差，是支持变化，不是gauge或新材料等价。native只留前三坐标后其余cell置零，728/728投影违反box；DFT/Haar的截断投影虽仍在box，却已改变材料与prior支持，未作新的forward。完整DFT位于满足共轭约束的6实维子空间，不能给ambient C6强设普通Lebesgue密度。没有把原P1 decoder强行逆到这个新的6-cell材料类。

CSI/SOM/TSOM早已有共同材料/电流lifting和确定/模糊电流分解；IFT、coarea与monodromy是标准数学。Everett1996/2026已有小型电磁/Helmholtz多项式路径寻找等价根的直接先例。ScholarQA两查询四操作429保留，5篇身份核验通过；SOM/TSOM与coarea有全文，Everett仅按出版社摘要级证据使用，CSI原文captcha未取得。详见audit/prior_art，未认证全球新颖性。

## 10. 费用、失败和未执行项

总actual RHS {resource['actual_complex_RHS']:,}，解析等价{resource['closed_form_source_RHS_equivalent']:,}，合计{resource['conservative_combined']:,}；实际LU{resource['LU_factorizations']:,}；最大RSS约{resource['peak_RSS']/1024**2:.1f} MiB。受控物理/代数job wall约{resource['physical_and_algebraic_job_wall_seconds']:.3f}秒、CPU约{resource['job_CPU_seconds']:.3f}秒；写作/检索/绘图/交付另有记录，不伪造为物理钟。每run合并费用小于10,000，零训练、零GPU/云。

两次实现错误为E0 np.bool_ JSON序列化和E2响应/状态sigma_min重名。旧失败日志、原源码和795已付RHS保留。修复通过缓存恢复E0、复用E2已付33启动及第一个循环节点，只完成剩余循环；没有为了元数据错误重跑或删除旧成本。六个高精度合并snapshot原来引用END序号，另有只追加的paid-source链接收据。原试探current过程由同一矩阵代数重放补存，374目标、2615代数调用，所有原端点m/z一致，新增RHS0。

未执行：原Full66新探索、4–9 cell全域根数认证、复杂有限覆盖monodromy、跨新采集/频率泛化、连续PDE认证、完整新6维带噪后验。未用预算不作为追加授权。

## 11. 七个最终问题

1. **找到什么结构？** 新六单元少接收系统中二维局部材料图册；小binary64系数模型有局部收缩证书，较大范围是保存并验证的采样覆盖。
2. **current独立价值？** 它先精确消掉接收方程，内迭代不需Full solve，实际RHS下降；但未知量从4增到24，未增加覆盖。
3. **是否超过普通方法？** 同目标覆盖相同；本轮两种普通对照中动态Full-J更快，公平总成本只省13.6%，未过20%。不宣布净加速突破。
4. **哪些只是gauge/换基/对称？** 旧两cell交换是人工对称阳性控制；旧QD是decoder gauge；完整DFT/Haar是坐标。新的原生实材料图册m改变，属于不同材料。
5. **概率与漏质量？** 旧两cell全域低维数值参考复现；新模型只在声明图册片正确coarea数值积分，全域漏质量未知。
6. **交给R8什么？** 非对称算子/材料/完整电流/精确账本、可复核小图册证书、相同覆盖的成本反证、以及局部积分与全域质量之间的明确缺口。
7. **现在进入Full66？** 不进入。先在小模型解决一个限定域的分支排除/覆盖证书或真正低于强普通基线的表示，不能用模型复杂度替代这一步。

## 12. 可复现证据

入口为START_HERE、REPRODUCE、PHYSICAL_EVIDENCE/BRANCH_AND_COVERAGE/PROBABILITY_EVIDENCE三CSV、RESULTS_SUMMARY、中央COST_LEDGER和完整无损证据包。独立审阅脚本不导入物理求解器；它从现有数组重算资格与费用。所有主要图来自本轮真实数组。原R6/R7结论保持不变。
''')
    md('THEORY_TO_EXPERIMENT_AUDIT.md',f'''# 理论到实验：哪些命题真正被检验

输入报告全文、源代码、协议与历史审计见inputs及audit/input_evidence、history_contract。原两cell对称、branchprior及monodromy仅回归，不算新发现。新实验编号按用户任务：E1非对称连续图册，E2独立离散分支，E3跨表示，E4概率；与附件原建议编号差异明确保留。

## Maxwell与材料约束

采用真实disk Green：自项 iπka H1(ka)/2−1；互项 iπka J1(ka)H0(kd)/2。完整状态L=I−diag(m)G，j=L⁻¹diag(m)Einc，y=Sj。不是Born。Full-J为(SL⁻¹)_i E_is；transpose解支付一个receiver RHS，source-major实化与sigma一致。

在接收方程Sj=u下，j=j0+Nc；shared-material minors为0只是非节点条件的一部分。实际算法保留实m并解j−mE=0，所以同时处理所有源、实值支持和节点。当前构造没有节点，节点一般仍须保留m。单独dark-current已经得到不可实现反例。

## 维数、局部证书与全域覆盖

n=6,k=4且合法内点、邻域满行秩时，标准IFT给二维局部纤维；这不是原创维数定理。interval_chart使用储存系数为精确常数的复矩形区间算术，R近似逆的Neumann余项包围L⁻¹，再包围J。T(d)=d−R_J(F(d,a)−u)将dependent box映到自身且带权导数范数<1时，每个自由a有唯一dependent根。三参考半宽为{res['interval_certified_free_halfwidths']}。不包含真实Hankel系数的严格区间误差、不认证整个材料域或continuum。

所有738端点的后付伴随J均满行秩，最大dependent条件数7.403；这只是节点诊断。原在线corrector未执行逐节点rank-stop，事后无越过阈值不追认为原算法在线证书。两个方法都没有认证跨过奇异chart、整个连通分量或全部原prior质量。

## 测度

材料图册m(a)=[d(a),a]的coarea密度p(m)/J_F * sqrt(det(DmᵀDm))，与p(m)/|det J_dep|相等；current恢复需要的是材料切向metric，不是高维c的平坦密度。E4只检验同一原prior的局部核。固定精确(p,q)时的0.75，与有限噪声的0.6751644709不同。均匀branch不能由群的存在推出。

## 原任务差异与未实现

原默认32/64双层E4的解析响应预算不适合保守30k合并账；本轮只新32并明确缓存64。最强普通基线具有最小接收伴随RHS；成本后补控制挑战而非挽救current结果。复杂monodromy因有限覆盖/判别式前提缺失NOT_RUN。完整Fourier/Wavelet成功只记换基。新模型不是P1，原Full66后验仍UNQUALIFIED。
''')
    md('FAILURE_ANALYSIS.md','''# 失败、修复、边界与未执行

1. E0原物理20RHS全部成功；写入np.bool_失败。原失败job/log/源码保留；np.generic.item修复只影响序列化。通过A/RHS/解精确缓存匹配恢复完整命名空间，0新RHS，不伪称新独立solve。
2. E2已完成33启动和第一个循环点时，响应与resolvent的sigma_min键重名。保留775RHS及INTERRUPTED.json，改名response_sigma_min；复用原33启动和首节点，只付870RHS完成其余循环。阈值、种子、方法不变。
3. E1二法各10个目标未合格，包含材料越界与边界/正确器停滞；不删除失败，不把它们算成低费用成功。33多启动仅31严格合格；其余2保留。
4. 小图册区间bound在较大box太松，按事前固定halving序列得到更小认证box。不能用这些小box认证0.05积分片或径向全路径。
5. 主cost窗口看似43.3%节时，但轨迹保存不匹配。公平buffered完整trace窗口总节时仅13.6%：保留旧数值并撤回其独立总成本主张，不改20%门。记录和代数变量膨胀会抵消RHS优势。
6. 六个70dps合并snapshot原action写END序号。原文件及ledger未改，HIGH_PRECISION_SNAPSHOT_LINK_ADDENDUM与追加事件明确关联两源已支付BEGIN；无费用丢失或额外物理。
7. 补存current中间c/j/E通过纯代数重放374目标，原端点一致，0新RHS；重放CPU另记，不把它当新发现或加速。
8. 文献两查询四操作429、CSI原文captcha、绘图system Python缺NumPy等行政失败保留。文献身份核验和可得全文层级分开；绘图换已有环境，不影响物理。

NOT_RUN：完整Full66/576研究、严格Hankel或连续PDE证书、全域根数/组件覆盖、复杂monodromy、新噪声/采集泛化、完整新模型带噪后验。没有通过降低容差、改noise、换先验或删接收救结果。
''')
    md('REPRODUCE.md','''# 复核与复现实验

本轮CLOSURE已封账，Book拒绝再次启动。不要删closure/ledger，也不要覆盖raw或input重置预算。

## 零新物理审查

在还原完整证据、安装相同NumPy/SciPy环境后，使用独立审阅脚本：

```bash
python audit/results_independent/audit.py --output /tmp/zalg-independent-review
```

具体脚本相对本研究目录；读取保存的A/RHS/current/J和中央账本，重算材料/数据/成本与概率资格。它不调用forward、不构造Book。图脚本只读已有数组。原输入报告和源码保持只读。

## 原执行顺序

E0解析/两cell→序列化缓存恢复→非对称probe→P0 A/B、P1 B/A、P2 A/B→小区间图册certificate→E3换基→E2多启动与往返数据continuation→元数据恢复→旧E4新32/缓存64→新局部coarea积分→70dps独立检查→公平完整trace成本控制→current代数trace重放→738端点伴随rank诊断→封账。所有失败也计费。

若要重新支付物理进行独立复现，须在新的输出目录、独立新授权与账本下保留相同EXPERIMENT_CONTRACT和源版本；不能把当前封账额度转给新项目。env限定OPENBLAS/OMP/MKL/VECLIB线程1，float64/complex128，单物理worker、一个live物理LU。库版本及SOURCE_MANIFEST保存。没有NN/GPU/cloud依赖。

原作者E4默认32/64不能直接重跑后只按10RHS记账；新适配只算32，明确复用64。闭式响应也按source-RHS等价单列。物理数组是无损压缩NPZ，可直接np.load(allow_pickle=False)。完整归档逐字节回读，阅读包与完整过程包分开标注。
''')
    md('GPT_PRO_HANDOFF/START_HERE.md',f'''# 给P3L-R8/下一轮Pro：先保留成本反证

本轮完整执行E0–E4，实际{resource['actual_complex_RHS']:,} RHS＋{resource['closed_form_source_RHS_equivalent']:,}闭式等价；不进入原Full66。先读RESEARCH_REPORT_ZH、GO_NO_GO和VERIFIED_RESULTS。

真实正面：非对称6-cell两照明，合法材料二维局部图册，两法同覆盖；小binary64系数模型有收缩证书；局部prior-coarea积分一致。

决定性负面：current24实校正变量对普通4变量没有新覆盖；公平同trace总成本仅省13.6% <20%；四cell31/33合格只回1根，无全域唯一性证明。不要把solver子项28.6%包装为总加速，也不要把局部数值概率包装为Full posterior。

完整原始数组和失败随证据包保存；新问题见OPEN_MATHEMATICAL_QUESTIONS。主要构件是CSI/SOM、IFT、coarea已知工具，不以代数名字宣称TAP创新。
''')
    md('GPT_PRO_HANDOFF/VERIFIED_RESULTS.md',f'''# 可复核结果

{main_table}

{costtable}

float64之外有70dps6点独立复核，最大相对差2.01e-15。738后付J诊断均rank4，maxdependent condition7.403；不认证节点之间。31/33多启动合格、1根、两32节点往返路径。E4旧32事件{res['probability']['legacy_new32_event']:.10f}与缓存64一致；新local4/6/8和两图册数值一致，全域漏质量未取得非平凡上界。

精确固定(p,q)两支权重1/4、3/4与有限噪声积分分开。完整DFT/Haar往返为换基。3确定性公开参考不是prior随机独立测试，不给总体置信区间。

复核文件：RESULTS_SUMMARY、三CSV、audit/results_independent/final_refresh、RESOURCE_SUMMARY、COST_LEDGER、配置与完整NPZ。初始冻结f543d6cc；元数据修复独立保留。CLOSURE禁止偷偷续算。
''')
    md('GPT_PRO_HANDOFF/OPEN_MATHEMATICAL_QUESTIONS.md','''# 三个最小且可证伪的后续问题

## 1. 能否认证一个非对称四cell材料域内的全部实根？

本轮31个合格启动仅回1根；没有域内排除或全根数证书。先固定本模型与盒，用区间分支排除或适当的稀疏多项式工具，区分“确实单根”与“存在漏分支”。只需回答这个限定域，不先上Full66；不能把多启动频率当概率。

## 2. 能否消掉电流lift的变量膨胀，真正减少总费用？

当前共享材料关系能省RHS，却要解24个实变量，普通材料只有4个。可研究受约束Schur消元或小n完整multiaffine有理响应缓存，但前者可能只是普通material-J，后者有2^n系数增长。两者均属成熟代数工具方向，尚未在本轮执行，不能写成已验证新算法。下一轮先比较冷启动/边际费用和同trace总时间，门仍20%。

## 3. 小局部coarea核如何变成有遗漏质量界的可计算材料条件核？

当前收缩认证box比数值积分片小，完整prior上可能还有其它图册/分支。需要可查的boundary atlas、overlap去重和pushforward证据下界，才能把未覆盖条件质量上界写成小数。仅把root绘成一片、或观察多次Gauss阶数一致，都不能排除未见分支。

在上述之一出现决定性优势前，不扩到原Full66，不追加Flow/Diffusion，也不把已有coarea/determinant或resolvent恒等式改名为新Maxwell理论。
''')
    md('START_HERE.md',f'''# Research Z-ALG-E1：非对称 Maxwell 材料图册

**已完成、已封账。主结论NO-GO：未证明当前current图册比强普通材料方法有独立覆盖或合格总成本优势。**

- [完整中文报告](RESEARCH_REPORT_ZH.md)：真实6-cell局部材料纤维、4-cell分支搜索与条件积分。
- [分项裁决](GO_NO_GO.md)：图册存在、覆盖、概率和成本分别判断。
- [理论到实验审计](THEORY_TO_EXPERIMENT_AUDIT.md)与[失败记录](FAILURE_ANALYSIS.md)。
- [GPT Pro交接](GPT_PRO_HANDOFF/START_HERE.md)与[复核入口](REPRODUCE.md)。

实际{resource['actual_complex_RHS']:,} complex RHS＋{resource['closed_form_source_RHS_equivalent']:,}闭式source等价＝{resource['conservative_combined']:,}/30,000保守预算。不是P1 Full66；三份公开确定性参考不是独立盲测parent。正确新局部coarea概率不等于全域posterior；全域漏质量未认证。

真实数组、完整失败与试探、中央账本、代码/配置和原输入均保留；外部文献全文只作本地阅读，不替代原出版社链接。
''')
    # Keep the original prephysics manifest intact while adding current-source identity.
    source_paths=list((P/'src').glob('*.py'))+[P/'build_site.py',P/'audit/results_independent/audit.py',P/'audit/solver_review/algebraic_qa.py']
    source={f.relative_to(P).as_posix():dict(bytes=f.stat().st_size,blake2b=hashlib.blake2b(f.read_bytes()).hexdigest()) for f in sorted(source_paths) if f.is_file()}
    write(P/'FINAL_SOURCE_MANIFEST.json',dict(initial_manifest='SOURCE_MANIFEST.json',source_files=source,closure='CLOSURE.json',core_changes='JSON conversion and duplicate-key reporting fixes; statistical/physical contract unchanged',generated_utc=time.time()))
    print(json.dumps(dict(status='REPORTS_WRITTEN',actual_RHS=resource['actual_complex_RHS'],combined=resource['conservative_combined'],optimized_total_saving=saving)))
if __name__=='__main__':main()
