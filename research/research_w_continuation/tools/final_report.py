"""Cached-only reports and ledger snapshots; never invokes Maxwell."""
from pathlib import Path
import json,csv,shutil,time
ROOT=Path(__file__).resolve().parents[1]
def read(p):return json.loads((ROOT/p).read_text())
def write(p,s):(ROOT/p).write_text(s.strip()+'\n')
q=read('P0_Q_QUALIFICATION.json');scores=read('results/p0b/SCORED_RESULTS.json');search=read('results/p0b/search-corrected/RESULTS.json')
checks=read('results/p0b/finer-endpoints/RESULTS.json');fitpath=ROOT/'results/p0b/fitted-task-endpoints/RESULTS.json'
fitted=read('results/p0b/fitted-task-endpoints/RESULTS.json') if fitpath.exists() else []
allrows=[];costs={};failures=[]
for stage,path in [('P0-Q','results/PHYSICS_LEDGER.jsonl'),('P0-B','results/p0b/PHYSICS_LEDGER.jsonl')]:
    rows=[json.loads(s) for s in (ROOT/path).read_text().splitlines() if s.strip()]
    reserved=[r for r in rows if r.get('event')=='reserve'];jobs=[r for r in rows if r.get('event')=='cost'];states=[r for r in rows if r.get('event')=='backend_state']
    costs[stage]=dict(reserved_complex_rhs=sum(r.get('rhs',0) for r in reserved),reserved_factorizations=sum(r.get('factorizations',0) for r in reserved),completed_backend_states=len(states),completed_factorizations=sum(r.get('event')=='backend_factor' for r in rows),completed_rhs=sum(r.get('rhs',0) for r in rows if r.get('event')=='backend_solve'),measured_job_process_cpu_seconds=sum(r.get('process_cpu_seconds') or 0 for r in jobs),measured_job_wall_seconds=sum(r.get('wall_seconds') or 0 for r in jobs),charged_cpu_seconds=sum(r.get('charged_cpu_seconds',0) for r in rows),peak_rss_bytes=max((r.get('peak_rss_bytes',0) for r in jobs),default=0),by_phase={phase:dict(rhs=sum(r.get('rhs',0) for r in reserved if r['phase']==phase),factorizations=sum(r.get('factorizations',0) for r in reserved if r['phase']==phase)) for phase in sorted(set(r['phase'] for r in reserved))})
    costs[stage]['recorded_process_cpu_seconds_total']=costs[stage]['measured_job_process_cpu_seconds']
    costs[stage]['measured_job_process_cpu_seconds']=sum(r.get('process_cpu_seconds') or 0 for r in jobs if r.get('command'))
    costs[stage]['measured_partial_region_cpu_seconds']=sum(r.get('process_cpu_seconds') or 0 for r in jobs if not r.get('command'))
    costs[stage]['CPU_scope_note']='Whole child-process jobs and partial worker numerical regions are separate; charged totals also include any explicitly declared preparation allocation.'
    for r in rows:
        allrows.append(dict(stage=stage,**r))
        if r.get('status')=='FAILED' or r.get('event') in ('budget_reject','backend_memory_reject','backend_residual_failure'):failures.append(dict(stage=stage,counts_already_in_physics_ledger=True,**r))
for r in search:
    if r['case']=='corrector_trial' and not r['accepted']:failures.append(dict(stage='P0-B',event='rejected_trial',counts_already_in_physics_ledger=True,**r))
failures.append(dict(stage='P0-B',event='implementation_repair',job='p0b-search-001',science_eligibility='EXCLUDED',reason='Corrector proceeded from rejected proposal; interrupted and repaired; all prior costs retained'))
failures.append(dict(stage='P0-B',event='protocol_implementation_gap',physics_actions=0,reason='Frozen minimum_feature=0.1 was not strictly enforced by the 48x32 pixel grid, Gaussian proposal smoothing and top-area projection; box/passivity/area/shared-frequency checks do not certify this condition.',candidate_positive_witness_eligibility=False,discovered='final review; no rerun or gate change'))
(ROOT/'PHYSICS_LEDGER.jsonl').write_text(''.join(json.dumps(r,ensure_ascii=False)+'\n' for r in allrows))
(ROOT/'FAILURE_LEDGER.jsonl').write_text(''.join(json.dumps(r,ensure_ascii=False)+'\n' for r in failures))
(ROOT/'COST_SUMMARY.json').write_text(json.dumps(dict(snapshot_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),stages=costs,scope='Merged derived snapshot; original stage ledgers append-only. Publication tail receipts may be later.',old_P0_costs_excluded=True,CPU_allowance_note='Q charged CPU includes 60s conservative uninstrumented preparation allocation; not measured process CPU.'),indent=2))
complete=[r for r in search if r['case']=='calibration_candidate'];missing=[r for r in search if r['case']=='context_not_run']
decision=dict(P0_Q=q['qualification'],P0_Q_gate='PASS',P0_B_calibration='NOT_FOUND_WITHIN_ALLOCATED_SEARCH_BUDGET',P0_B_joint='NOT_RUN_NO_CALIBRATION_WITNESS',P0_B_task_numerics='UNQUALIFIED_UNFITTED_TARGET_EXTRACTION' if not fitted else 'THREE_SELECTED_CIRCLE_FITTED_ENDPOINTS_CONSISTENT_NOT_CONTINUUM_BOUND',P1='HOLD_NOT_RUN',complete_contexts=len(complete)//4,planned_contexts=4,independent_reference_seeds=2,minimum_observed_dcal=min(r['d_cal'] for r in complete),unexecuted_contexts=missing,task_class_evaluated='fixed probe8 only; 1/24, no full supremum',qualified_witnesses=0,posterior_computed=False,physics_uniqueness_proved=False,task_aware_advantage='NOT_ESTABLISHED',source_commit='db7146cae2d74ba231577d42f225dbfc9a389268')
decision['P0_B_frozen_minimum_feature']=0.1
decision['P0_B_minimum_feature_enforced']=False
decision['P0_B_candidate_positive_witness_eligibility']=False
decision['P0_B_additional_limitation']='PROTOCOL_IMPLEMENTATION_GAP: box/passivity/area/frequency consistency checked; frozen minimum feature not certified'
(ROOT/'GATE_DECISION.json').write_text(json.dumps(decision,indent=2))
for name in ['P0_B_CALIBRATION_RESULTS.csv','P0_B_JOINT_RESULTS.csv']:shutil.copyfile(ROOT/'results/p0b'/name,ROOT/name)
write('P0_Q_GATE_DECISION.md',f'''# P0-Q 判决

**PASS / QUALIFIED_DISCRETE_FULLWAVE**。原两seed、两带、cal/task的8项完整数据最细差为0.006533–0.033224，全部低于原0.1。固定几何P2→P3、独立plane/line Mie、PML强度/厚度/外域、体积分/点接收、真实导数及伴随共同支持此有限资格。

这不是CERTIFIED_CONTINUUM_MAXWELL，也不是对全部未知材料、任意重网格或拓扑路径的统一保证。后续新候选需独立复核。没有放宽门槛、更换seed、噪声、频带或源强。旧P0结论保持原样。

逐项数据见[P0_Q_QUALIFICATION.json](P0_Q_QUALIFICATION.json)。本次Q快照：{costs['P0-Q']['reserved_complex_rhs']}复RHS、{costs['P0-Q']['reserved_factorizations']}次已预约因子、{costs['P0-Q']['charged_cpu_seconds']:.3f}计费CPU秒。发布尾项可能继续追加。
''')
table='|参考seed / 带宽|Full-J|随机匹配|双初始化FWI|Task-aware|\n|---|---:|---:|---:|---:|\n'
for seed,band in sorted(set((r['seed'],r['band']) for r in complete)):
    group={r['method']:r for r in complete if r['seed']==seed and r['band']==band}
    table+=f'|{seed} / {band}|'+'|'.join(f"{group[m]['d_cal']:.6f}" for m in ['FULL_J','RANDOM_MATCHED','MULTISTART_FWI','TASK_AWARE'])+'|\n'
write('P0_B_GATE_DECISION.md',f'''# P0-B 判决

**NOT_FOUND_WITHIN_ALLOCATED_SEARCH_BUDGET；joint为NOT_RUN_NO_CALIBRATION_WITNESS。** 已完成3个context×4方法；第二个seed的窄带因仅余480搜索RHS、不足完整800RHS四方法预留而NOT_RUN。没有把剩余endpoint/validation预算转为搜索，没有把旧失败费用退回。

{table}

最小cal距离19.912259，远未满足1。没有GO-MECHANISM或GO-JOINT；不能从搜索失败证明唯一性。固定目标8只提供1/24任务类诊断，而且目标与参考相同，质心差为0，不能冒充未知目标joint补偿。普通随机在首seed窄带更好，Task-aware在第二seed宽带略好，另一个context全回退BP；两seed不足以宣称独立算法优势。

最终审阅确认另有协议实施缺口：原minimum_feature=0.1并未由48×32像素、Gaussian平滑提案和top-area投影严格执行。dx=0.08333、dy=0.09375，单像素孤立结构可能小于原限制。box、passivity、面积和跨频材料律通过不等于全部合法材料条件通过。所有候选均不具备正向见证资格；这里没有临时改小最小特征，也不据该搜索排除原合法类中的多解。

首次搜索的正确器回退错误版本排除科学比较，760搜索RHS和133.673416 CPU秒保留。修正版仍只是短预算、一次共同Full-J正确器的局部探索，不是已收敛全局FWI、多启动覆盖或不可辨识性证明。

原像素网格cal端点的细化差0.000586–0.004319；未贴合目标圆盘的task差0.198–0.270，超过0.1。二者分开判断。目标界面补充复核见[端点验证](P0_B_PHYSICS_VALIDATION.md)，未达标结果不用于物理GO。

[机器判决](GATE_DECISION.json)；[逐候选结果](P0_B_CALIBRATION_RESULTS.csv)。P1、P0-C大型实验、A17 ROM、NN和三维均未启动。
''')
write('P0_B_WITNESS_SEARCH.md',f'''# 四方法见证搜索

逆算法读取测量、已知背景、噪声、完整48×32材料像素和固定数学任务，不读取wall_values、真实mask或孔洞数。初值为在线background Born/adjoint解及公开prior；所有源/频率共享同一contrast与conductivity。投影统一满足box、passivity和最多3面积，精确空气允许组件生灭。像素/平滑先验仍限制可表达的最小尺度，未建立任意拓扑搜索覆盖保证。

{table}

共同正确器用完整离散Full-J、同一ridge/trust幅度和最多4个标准半步trial；接受只依据cal全波失配和同prior。N3、truth最优步、未知墙体Green均没有作为输入。multistart在本pilot仅为共同BP与一个投影到box/面积的随机扰动初始化，不是假称充分多启动。

最小特征缺口：原0.1限制未被该表示/投影强制。本文的“合法”仅在明确的box、损耗、面积和共享频率律层面核验，不能读成全部冻结几何先验获认证。新一轮必须先实现该条件；本轮不事后改变材料类，也不重新搜索。

任务方向来自归一化材料坐标下K与J的正则广义比，只是local proposal。它没有稳定超过简单控制，也没有达到非线性数据相容。此轮不能判断局部危险方向能否沿真正的近等价纤维有限延拓。

原始[修正版轨迹](results/p0b/search-corrected/RESULTS.json)、[排除的错误版本](results/p0b/search/RESULTS.json)、[离线truth评分](results/p0b/SCORED_RESULTS.json)。后者在在线输出保存后运行，不参与选择。表中不同带宽属于同两parent，不是三独立对象；没有伪造bootstrap显著性或posterior概率。

费用表的method wall只覆盖增量proposal/corrector；共享setup、BP、F/J、任务方向与数学求解另包含在整个作业wall/CPU。RHS standalone归因为普通控制+80共享cal RHS、task-aware+160共享cal/task RHS。实际共享费用只在总账计一次；不能拿增量wall当端到端加速。
''')
finer='|seed / 带宽 /状态|细化白化差|\n|---|---:|\n'+''.join(f"|{r['seed']} / {r['band']} / {r['state']}|{r['mesh_white']:.6f}|\n" for r in checks)
fittext='尚未完成贴合目标界面补充检查；不得据原task结果宣布GO。'
if fitted:
    fittext='补充贴合目标界面的实际结果（同一像素世界、同一个probe8，P3，圆周256→512段）：\n\n|seed / 带宽|完整目标响应细化白化差|最细目标数据距离|最大RHS相对残差|\n|---|---:|---:|---:|\n'
    for r in fitted:
        if r['polygon_n']==512:fittext+=f"|{r['seed']} / {r['band']}|{r['fitted_refinement_white']:.9f}|{r['data_distance_to_reference']:.6f}|{r['maximum_rhs_relative_residual']:.3e}|\n"
    fittext+='\n三项均低于原0.1。为保持全部未知材料像素界面，圆周与像素网格交点会产生细瘦三角形；没有平滑物理几何、放宽1e-9求解残差或加jitter。网格最小角及aspect诊断公开，P3 DOF为35899/38203，复核峰值约1.065GB。这只是三端点的圆周细化一致性，不是全部12候选合格，也不是连续PDE上界。见[原始补充结果](results/p0b/fitted-task-endpoints/RESULTS.json)及[MESH_DIAGNOSTICS](results/p0b/fitted-task-endpoints/MESH_DIAGNOSTICS.json)。'
write('P0_B_PHYSICS_VALIDATION.md',f'''# 新候选数值复核

每个已完成context按在线cal失配最小、stable method ID打破并列选择一个端点，未用truth选择。原P3像素网格与全部三角形再细分后的P2做相同材料/源/噪声对比。

{finer}

cal一致不等于task一致。未贴合圆形目标的积分产生额外几何误差，不能把原参考世界的贴合网格资格扩展给所有新材料/任务。目标补充复核在既有validation预算，使用同一目标、相同连续总材料，只修正圆周界面处理。

{fittext}

[原复核结果](results/p0b/finer-endpoints/RESULTS.json)；实际NPZ复场/观测见Release的[数组清单](ARRAY_MANIFEST.json)。没有任何candidate达到cal<=1，没有最终GO端点微扰/全部24probe或joint复核；这些项目为NOT_RUN，不写成已成功。
''')
write('SCIENTIFIC_DECISION.md',f'''# 科学决策

- P0-Q：**PASS，QUALIFIED_DISCRETE_FULLWAVE**，取得可实际使用的二维高阶前向平台；不是新算法优势或连续认证。
- P0-B：**NOT_FOUND_WITHIN_ALLOCATED_SEARCH_BUDGET**。最好cal距离19.912259；3/4context完整比较，错误搜索消耗和缺失行保留。
- 算法增量：**NOT_ESTABLISHED**。Task-aware未稳定超过随机/Full-J，且所有结果仍未进入d<=1的数据相容区。
- joint：**NOT_RUN**。没有把同目标响应或cal差异替代未知目标补偿。
- 表示先验实施：**PROTOCOL_IMPLEMENTATION_GAP**。原0.1最小特征未被当前像素/投影严格执行，所有搜索端点不能作为正向物理见证；失败不排除原合法材料类的多解。
- P1：**HOLD**。本轮不支持启动posterior/progressive、A17加速、NN或三维。

当前最优先研究是让普通受约束FWI在原数据/先验下达到足够好的cal拟合，并把几何投影、优化收敛和表示误差区分开；随后再问有限非线性材料分支能否保持数据近等价。若这个基础点都未找到，调大task比值或构建复杂概率架构没有依据。

只有实际qualified joint witness并有扰动证据才建议P1-A；P1-B还必须有同数据同prior四格消融。若未来cal通过而宽带joint击破，则转向传播可辨识性/采集充分性；当前尚未取得这种证据，不能预先REDIRECT或证明唯一性。
''')
write('RESEARCH_W_CONTINUATION_REPORT_ZH.md',f'''# Research W 续接：先修准物理，再检验传播歧义

**本轮结论：P0-Q PASS；P0-B未发现合格见证；P1继续HOLD。** 这是一个数值基础取得进展、物理多解问题仍未解决的结果。原P0的失败结论不改写，本轮不把更准确的前向平台说成更好的成像算法。

## 1. 研究从哪里出发

旧P0的一维反例成立：单频墙反射相同而双程传播不同；有限带宽会打破某个指定的联合补偿。0.5%下cal=0.887、joint=3.038，是原已完成证据，本轮没有重新充当新结果。旧二维矩阵残差和伴随虽好，真实圆柱误差仍约7%，两墙网格差约1个白化单位，因而当时没有搜索。

本轮保持原两堵不规则有损墙、4源12接收器、两套5频点、sigma=0.01、有限尺寸和多次散射。每频率共享epsilon_infinity与conductivity，没有读真墙mask帮助逆算法，没有增加观测或调大噪声。

## 2. 旧7%主要是什么问题

直接提升非自Green积分GL4到16，nx96圆柱误差只从9.84519%变成9.84374%，说明核积分不是这个误差的主要来源。圆中心二值采样换成精确面积fraction有改善，但仍很差；墙体4/8/16/32子采样也会改变数据。双单元平均VIE没有改善，负结果保留。

独立FEM不使用旧VIE矩阵。它求散射场、用解析入射避免离散奇异点源，外域采用显式坐标拉伸PML。固定同一贴合圆周网格和积分，P2到P3把平面波相对误差从2.029%降到0.1367%；进一步细配置为0.03357%，实际线源为0.02662%、白化差0.001716。界面相容和场阶数都有直接作用，但没有精确百分比分摊，不能宣称全部误差只来自某一项。

![实际独立参考与受控改动](figures/01_cylinder_accuracy.png)

源/接收一致性也单独修复：体积分输出保留同一连续点接收物理，互易性约3.9e-15，真实材料差分最好1.49e-9、实伴随约3.15e-14。原FE点输出保留为对照，不拟合自由复增益。

## 3. 两堵原墙体真正过了哪些检查

两seed、wide/narrow、cal/同probe8共8个完整数据向量的最细网格差0.006533–0.033224，均<=0.1。PML强度、厚度/外域和输出提取复核也相容。两世界差分误差作为额外诊断，并没有替换单世界门槛。由此给出QUALIFIED_DISCRETE_FULLWAVE；没有严格连续PDE误差界，不能写CERTIFIED_CONTINUUM_MAXWELL。

![八项完整数据细化](figures/02_wall_convergence.png)
![原世界材料和实际总场](figures/04_known_world_and_field.png)
![已保存中心频率总场的幅值与相位](figures/11_saved_task_total_field.png)

图中是已知参考世界和真实数值场，不是逆成像成功图。资格只覆盖这些物理输入；新候选仍需检查。原墙source/receiver、无源真空、圆柱、平滑材料、被动性、损耗和跨频共享没有更换。

## 4. 搜索实际做了什么

全域48×32两实材料像素、同一box/面积/先验、同一online BP初值；Full-J、随机匹配、双初始化FWI、task-aware共享材料坐标与正确器权限。目标方向仅由当前材料的完整Maxwell切线构造，不使用reference材料或truth。所有在线输出先存盘，之后独立truth评分。

{table}

三个context完成12个候选，第四个context缺失。最好的19.912仍离1很远。首seed窄带随机比task-aware好；第二seed宽带task-aware略好，另一context全部回退BP。没有稳定同成本优势，也没有得到一组数据近等价的物理世界。

首次搜索有正确器回退实现错误，已停止、修复，旧版本排除科学比较；760搜索RHS和133.673416 CPU秒没有退还。修正版仍是短预算局部搜索，不是已收敛全局FWI。这是优化/探索不足，不能把它写成Maxwell唯一性。余480搜索RHS不足完整800RHScontext，未借端点/验证预算补搜索。

最终审阅还确认：原材料minimum_feature=0.1没有被48×32像素、提案平滑和top-area投影严格执行。像素边长0.08333/0.09375，可能产生较小孤立特征。box、损耗、面积和跨频一致不代表全部几何先验合法。此轮明确记PROTOCOL_IMPLEMENTATION_GAP，全部端点不具备正向见证资格；不改变原条件，也不把失败当成原合法类不存在多解的证据。下一轮应先修正该表示约束，再问普通反演与分支搜索的问题。

![实际cal距离和RHS归因](figures/06_p0b_calibration_and_rhs.png)
![参考墙、BP与在线选定端点的离线对照](figures/07_p0b_offline_material_comparison.png)
![观测误差与RHS费用散点](figures/09_p0b_calibration_cost_scatter.png)

## 5. 含目标和joint不能混淆

对所有最终候选实际计算固定probe8的完整强目标响应，只有1/24任务类，并另列目标诱导差。因为cal距离本来很大，即使同目标总响应差很大，也不构成传播歧义。相同目标的质心差为0，不是未知目标joint补偿。

进一步网格审核发现：最佳候选的cal差仅0.000586–0.004319，但未贴合目标圆盘的task差0.198–0.270，超过0.1。该项不合格结果明确保留，并在原验证额度补做贴合目标界面复核，详见[P0_B_PHYSICS_VALIDATION](P0_B_PHYSICS_VALIDATION.md)。没有任何cal端点达到1，故joint、全24任务、最终见证微扰均没有启动。

补充圆周256→512段、相同P3和像素材料的三个选定端点，完整目标响应细化差0.000644–0.000923，RHS相对残差最大2.80e-10。它支持把原task偏差归因于未贴合目标界面，但不把所有候选或任意新材料都自动升级为合格见证。

![同一端点贴合界面前后的数值一致性](figures/10_p0b_task_refinement_comparison.png)

![原未贴合任务输出及其限制](figures/08_p0b_fixed_probe8.png)

## 6. 成本、失败和证据边界

|新阶段|已预约复RHS|已预约因子|完成backend states|计费CPU秒|实际作业wall秒|
|---|---:|---:|---:|---:|---:|
|P0-Q|{costs['P0-Q']['reserved_complex_rhs']}|{costs['P0-Q']['reserved_factorizations']}|{costs['P0-Q']['completed_backend_states']}|{costs['P0-Q']['charged_cpu_seconds']:.3f}|{costs['P0-Q']['measured_job_wall_seconds']:.3f}|
|P0-B|{costs['P0-B']['reserved_complex_rhs']}|{costs['P0-B']['reserved_factorizations']}|{costs['P0-B']['completed_backend_states']}|{costs['P0-B']['charged_cpu_seconds']:.3f}|{costs['P0-B']['measured_job_wall_seconds']:.3f}|

这是报告生成快照，后续绘图/打包/网页读回仍追加。Q费用含60秒未仪表化准备的保守分配，不能把它叫测量CPU。预约包含在途/失败调用，不因SIGTERM退还；完成数单列。旧P0的959.069451计费秒完全分开。

所有LU、forward、receiver transpose、JVP/VJP和失败均先预约；Green核与材料积分/小矩阵CPU在setup/完整作业里，不冒充线性系统RHS。不同方法增量wall不含共享setup，不据此宣称完整加速。PML作业峰值2.160GB超过1.5GB工程目标，未发生系统OOM；搜索改为逐频释放，其修正版峰值约1.059GB。

保留输出序列化故障、内存预检拒绝、被拒绝trials及错误搜索中止；没有删除代价换新预算。[总账](PHYSICS_LEDGER.jsonl)、[失败](FAILURE_LEDGER.jsonl)、[费用快照](COST_SUMMARY.json)。本轮6项新合同测试通过及真实导数检查，只是实现证据，不替代Maxwell资格。没有NN、GPU campaign、posterior、三维或A17新ROM。

## 7. 十个问题逐项回答

1. **旧资格失败源于什么？** 主要是界面/系数与低阶场表示的空间误差组合，不是LU未解准；提高核积分几乎无效，匹配场阶数对照和界面细化给出直接证据。具体占比尚未量化。
2. **是否已有成本可接受的平台？** 对冻结两世界，是；高阶贴合FEM实际完成全部多频资格，CPU远低于本轮上限。不是任意高维搜索都已低延迟，也没有和成熟FEM/BEM竞速。
3. **能否支持d<=1判定？** 对已审核输入可做有数值裕度的离散判定；新的材料、目标/网格需要自己的审核，无统一连续误差证书。
4. **二维cal近等价且传播不同的合法世界是否找到？** 没有。最佳cal约19.91，有限搜索未进入相容区，不证明不存在。
5. **未知目标joint歧义是否持续？** 未检验。未把cal或同目标响应冒充joint。
6. **多次散射放大还是消歧？** 旧一维指定joint对会被带宽击破，旧离散四块显示Born遗漏强反馈；新的合格二维近等价对仍缺失，因此两种作用不能在此裁决。
7. **Task-aware同成本更好吗？** 未建立；观测到一优一劣一平，但均不合格，且task导数/共享setup有额外费用。
8. **目前支持哪类论文？** 目前只支持数值资格修复与受限搜索负结果，尚不足以支持概率传播或Progressive Imaging算法优越性论文。传播可辨识性是仍值得验证的物理问题。
9. **最可能的独立TAP理论贡献？** 在合法非规则材料与目标共同补偿下，连接局部观测弱方向、有限非线性任务歧义模数下界及宽带joint消歧条件。它是研究预测，尚未建立；四块/resolvent和同空间Schur恒等式不算新贡献。
10. **下一轮做什么、何时收缩？** 先证明普通受约束FWI能在原数据/先验下实际拟合到相容尺度，再比较同预算有限分支发现；不得让未收敛优化冒充物理限制。若cal分支被宽带joint击破，转测量充分性；有稳健joint分支才考虑低维P1-A，P1-B必须额外四格消融。没有收益就删除progressive/A17模块。

## 8. 后续只有方案，没有被自动执行

[P0-C相图](P0_C_IDENTIFIABILITY_PLAN.md)、[P1-A posterior参考](P1_A_POSTERIOR_PLAN.md)、[P1-B四格调度](P1_B_PROGRESSIVE_PLAN.md)、[P2-A费用比较](P2_A_A17_COMPUTATION_PLAN.md)、[P2-B拓扑扩展](P2_B_TOPOLOGY_SCALING_PLAN.md)、[P3矢量](P3_VECTOR_MAXWELL_PLAN.md)均为未来协议草案，不是结果。最终裁决见[SCIENTIFIC_DECISION](SCIENTIFIC_DECISION.md)。
''')
write('START_HERE.md','''# Research W continuation：高精度前向通过，物理歧义尚未得到见证

**P0-Q PASS / QUALIFIED_DISCRETE_FULLWAVE；P0-B NOT_FOUND_WITHIN_ALLOCATED_SEARCH_BUDGET；P1 HOLD。**

- 原两墙完整数据细化差0.0065–0.0332，低于原0.1；高阶贴合FEM的线源圆柱相对误差0.0266%。
- 三个context四方法完成，最好cal白化距离19.91，未达到1；第四个context缺失。没有joint、posterior、NN或三维。
- 搜索还存在明确实施缺口：没有严格执行原0.1最小特征条件；候选不能用于正向物理见证，不据失败推断唯一性。
- 不把搜索失败当唯一性，也不把数值平台通过当成成像优势。

[完整中文报告](RESEARCH_W_CONTINUATION_REPORT_ZH.md) · [科学裁决](SCIENTIFIC_DECISION.md) · [机器判决](GATE_DECISION.json) · [复现](REPRODUCE.md)

[圆柱图](figures/01_cylinder_accuracy.png) · [八项收敛](figures/02_wall_convergence.png) · [真实参考世界/场](figures/04_known_world_and_field.png)

代码在src/wq，原始指标与小数组在results，大型复场数组在同仓库Release research-w-continuation-20261023，路径见ARRAY_MANIFEST.json。旧P0保持原目录/commit，费用分开。
''')
shutil.copyfile(ROOT/'START_HERE.md',ROOT/'README.md')
print(json.dumps(dict(costs=costs,decision=decision),ensure_ascii=False))
