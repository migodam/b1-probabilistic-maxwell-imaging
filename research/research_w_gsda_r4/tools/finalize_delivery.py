"""Build Chinese evidence reports from sealed R4 artifacts; never calls physics."""
from pathlib import Path
import csv,json,math,statistics,time,sys,shutil

ROOT=Path(__file__).resolve().parents[1]

def read(path,default=None):
    p=ROOT/path
    return json.loads(p.read_text()) if p.exists() else default

def write(name,text):
    (ROOT/name).write_text(text.rstrip()+'\n')

def save(name,obj):
    write(name,json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False))

def number(x,d=4):
    return '—' if x is None else f'{x:.{d}f}'

def table(headers,rows):
    return '| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+'\n'.join('| '+' | '.join(str(x).replace('|','/') for x in row)+' |' for row in rows)

def csvsave(path,rows):
    if not rows:
        (ROOT/path).write_text('status\nNO_ROWS\n')
        return
    fields=sorted(set().union(*(x.keys() for x in rows)))
    with (ROOT/path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def main():
    records=read('evaluation/SCORES_ALL.json',{}).get('records',[])
    digits=[x for x in records if x.get('cohort')=='digit']
    cohorts=['B0_FULL_FWI','FULL_GN_TV_SUPPLEMENT','B4_GSDA']
    short={'B0_FULL_FWI':'B0 Full-FWI+TV','FULL_GN_TV_SUPPLEMENT':'B0+Full-GN+TV','B4_GSDA':'B4 原生CSI+风险筛选'}
    groups={m:sorted([x for x in digits if x.get('campaign')==m],key=lambda x:x['digit']) for m in cohorts}
    metric='target_complex_contrast_relative_L2_f1'
    def passed(x):
        v=x['metrics']
        return v['support_IoU']>=.7 and v[metric]<=.2 and v['contrast_centroid_error']<=.1
    def mean(g,key):return statistics.mean(x['metrics'][key] for x in g) if g else None
    summary=[]
    for name,g in groups.items():
        summary.append(dict(campaign=name,actual_digit_runs=len(g),passes=sum(passed(x) for x in g),mean_IoU=mean(g,'support_IoU'),mean_complex_relative_L2=mean(g,metric),centroid_RMSE=math.sqrt(statistics.mean(x['metrics']['contrast_centroid_error']**2 for x in g)) if g else None,mean_train_white_RMS=mean(g,'observation_white_RMS'),mean_heldout_white_RMS=statistics.mean(x['metrics']['prediction_heldout_clean_white_RMS'] for x in g if x['metrics'].get('prediction_heldout_clean_white_RMS') is not None) if any(x['metrics'].get('prediction_heldout_clean_white_RMS') is not None for x in g) else None))
    save('evaluation/FINAL_METHOD_SUMMARY.json',summary)
    summarytable=table(['方法','实际数字数','IoU均值','复对比度相对L2均值','质心RMSE / λ0','训练白化RMS','留出白化RMS','完整成像门通过'],[[short[x['campaign']],x['actual_digit_runs'],number(x['mean_IoU']),number(x['mean_complex_relative_L2']),number(x['centroid_RMSE']),number(x['mean_train_white_RMS']),number(x['mean_heldout_white_RMS']),str(x['passes'])+'/'+str(x['actual_digit_runs'])] for x in summary])
    detailed=[]
    for m,g in groups.items():
        for x in g:
            v=x['metrics'];c=x['cost'];top=v['predicted_topology']
            detailed.append(dict(digit=x['digit'],case=x['case_id'],method=short[m],campaign=m,actually_run=True,status=x['online_status'],IoU=v['support_IoU'],relative_complex_L2=v[metric],contrast_RMSE=v['target_real_contrast_RMSE'],loss_RMSE=v['target_loss_RMSE'],centroid_error=v['contrast_centroid_error'],true_holes=v['true_topology']['holes_bg4'],predicted_holes=top['holes_bg4'],true_holes_distinct=v['all_true_holes_distinct_enclosed'],predicted_components8=top['components8'],offsupport_contrast_mass_fraction=v['offsupport_positive_contrast_mass_fraction'],train_white_RMS=v['observation_white_RMS'],heldout_white_RMS=v.get('prediction_heldout_clean_white_RMS'),full_FE=c.get('full_fe'),complex_RHS=c.get('rhs'),CPU_seconds=c.get('wrapper_process_cpu_seconds'),wall_seconds=c.get('wall_seconds'),cost_scope=x.get('cost_basis','B4 per-case warm only, shared setup separately' if m=='B4_GSDA' else 'actual wrapper'),converged=x.get('converged'),projected_gradient_inf=x.get('selected_scaled_projected_gradient_inf'),all_image_gates=passed(x),endpoint=x['endpoint_dir']))
    csvsave('DIGITS_1_TO_9_METRICS_SUMMARY.csv',detailed)
    digit_tables=[]
    for m in cohorts:
        g=[x for x in detailed if x['campaign']==m]
        digit_tables.append('### '+short[m]+'\n\n'+table(['数字','状态','IoU','复材料相对L2','质心误差','孔洞 真/重建','训练RMS','留出RMS','FE / RHS','CPU秒'],[[x['digit'],'实际完成；'+('局部一阶达标' if x.get('converged') else '未认证收敛'),number(x['IoU'],3),number(x['relative_complex_L2'],3),number(x['centroid_error'],3),str(x['true_holes'])+'/'+str(x['predicted_holes']),number(x['train_white_RMS'],3),number(x['heldout_white_RMS'],3),str(x['full_FE'])+' / '+str(x['complex_RHS']),number(x['CPU_seconds'],2)] for x in g]))
    alltable='\n\n'.join(digit_tables)
    mesh=read('results/independent_validation/SUMMARY.json',{}).get('rows',[])
    meshtable=table(['场景','真实世界粗细网格白化差','原0.1门','估计终点网格差'],[[x['case'],number(x['truth_white_mesh_distance'],6),'PASS' if x['within_original_tolerance'] else 'FAIL',number(x.get('endpoint_white_mesh_distance'),6)] for x in mesh])
    controls=[x for x in records if x.get('cohort')!='digit' and (str(x.get('campaign','')).startswith('CONTROL_') or x.get('campaign')=='B0_PAIRED_CONTROLS')]
    controlrows=[]
    for x in controls:
        v=x['metrics'];c=x['cost'];mech=x.get('mechanism',{})
        controlrows.append(dict(case=x['case_id'],campaign=x['campaign'],method=x['method'],status=x['online_status'],objective=x.get('objective'),IoU=v['support_IoU'],relative_complex_L2=v[metric],centroid_error=v['contrast_centroid_error'],train_white_RMS=v['observation_white_RMS'],heldout_white_RMS=v.get('prediction_heldout_clean_white_RMS'),ABSORB=mech.get('ABSORB_count',0),REOPEN=mech.get('REOPEN_count',0),paid_target_screens=mech.get('paid_target_risk_screens',0),full_FE=c.get('full_fe'),complex_RHS=c.get('rhs'),CPU_seconds=c.get('wrapper_process_cpu_seconds'),wall_seconds=c.get('wall_seconds'),endpoint=x['endpoint_dir']))
    csvsave('METHOD_ABLATIONS.csv',controlrows)
    expected=[]
    actual={(x['campaign'],x['case']) for x in controlrows}
    for m in ['B1','B1_RESET','B2','B3','B4','B5']:
        for case in ['c010','c011','c012']:
            if ('CONTROL_'+m,case) not in actual:expected.append(dict(campaign='CONTROL_'+m,case=case,status='NOT_RUN_OR_UNSEALED',reason='No sealed scored output at report generation; never counted as success'))
    csvsave('METHOD_ABLATIONS_NOT_RUN.csv',expected)
    ctable=table(['场景','方法','全波Φ','IoU','复材料相对L2','吸收/重开','FE/RHS','CPU秒'],[[x['case'],x['campaign'],number(x['objective'],2),number(x['IoU'],3),number(x['relative_complex_L2'],3),str(x['ABSORB'])+'/'+str(x['REOPEN']),str(x['full_FE'])+'/'+str(x['complex_RHS']),number(x['CPU_seconds'],2)] for x in controlrows])
    cost=read('PHYSICS_COST.json',{})
    totals=cost.get('totals',{})
    csummary=read('results/B4_GSDA/BATCH_RECEIPT.json',{})
    gn=read('results/FULL_GN_TV_SUPPLEMENT/GN_COST_AND_STOPS.json',{})
    b0=read('results/B0_FULL_FWI/B0_COST_AND_NONCONVERGENCE.json',{})
    born=[]
    for case in ['c010','c011','c012']:
        o=read('results/absorption_controls/'+case+'/BORN_COMPARISON.json',{})
        born.append(dict(case=case,details=o))
    save('evaluation/FINAL_EVIDENCE_INDEX.json',dict(updated_unix=time.time(),digits=summary,controls_scored=len(controlrows),controls_pending_or_not_run=expected,mesh=mesh,cost_totals=totals,born=born,source_files=['evaluation/SCORES_ALL.json','PHYSICS_COST.json','results/FULL_GN_TV_SUPPLEMENT/GN_COST_AND_STOPS.json','results/B4_GSDA/BATCH_RECEIPT.json']))
    comparison=[]
    for x in detailed:
        comparison.append(dict(cohort='digit',**x))
    for x in controlrows:comparison.append(dict(cohort='control',**x))
    csvsave('BASELINE_COMPARISON.csv',comparison)
    experiments=[dict(experiment='E1 exact absorption',status='ACTUAL_SAME_DISCRETE_PASS',cases='c010,c011,c012 × three frequencies',evidence='results/absorption_controls/CONTROL_SUMMARY.json',limits='Same discrete operator; forced reopens do not prove automatic risk safety'),dict(experiment='E2 weak/strong mode truncation',status='ACTUAL_FAILURE_OF_48_MODE_PROJECTED_FIELD',cases='c010,c011,c012 at f0',evidence='results/spectral/SUMMARY.json',limits='Projected weak-Schur is exact only in retained48 space; main preconditioner retains all outside modes'),dict(experiment='E3 reversible controllers',status='ACTUAL_PARTIAL_OR_COMPLETE' if controls else 'NOT_RUN',cases=str(len(controlrows))+' sealed scored control endpoints',evidence='METHOD_ABLATIONS.csv',limits='Scheduled absorption/reopen and automatic risk-triggered actions must be reported separately'),dict(experiment='E4 independent unknown-wall free-material controls',status='ACTUAL_PARTIAL_OR_COMPLETE' if controls else 'NOT_RUN',cases='c010,c011,c012',evidence='BASELINE_COMPARISON.csv',limits='3x2 wall material class; air gap is boundary-connected, not arbitrary unknown enclosed cavity topology'),dict(experiment='E5 digits1..9 B0',status='ACTUAL_9_OF_9_RUNS_0_IMAGE_PASSES',cases='c001..c009',evidence='DIGITS_1_TO_9_B0_METRICS.csv',limits='0/9 B0 endpoints met frozen first-order tolerance'),dict(experiment='E5 digits1..9 B4',status='ACTUAL_9_OF_9_RUNS_0_IMAGE_PASSES',cases='c001..c009',evidence='DIGITS_1_TO_9_B4_METRICS.csv',limits='All KEEP_ACTIVE; automatic absorption and paid target-risk screen never executed'),dict(experiment='Full-GN+TV baseline refinement',status='ACTUAL_9_OF_9_DATA_EXPOSED_SUPPLEMENTS',cases='c001..c009',evidence='results/FULL_GN_TV_SUPPLEMENT/GN_COST_AND_STOPS.json',limits='B0 warm start, cumulative costs;1/9 local KKT, no blind/global optimum claim'),dict(experiment='9x13 material refinement',status='NOT_RUN',cases='',evidence='',limits='Time protected nine complete physical inversions and paired controller controls'),dict(experiment='multiple independent noise repeats',status='NOT_RUN',cases='',evidence='',limits='One noise realization per parent; no confidence interval/generalization-rate claim')]
    csvsave('PAPER_EXPERIMENTS.csv',experiments)
    nogo=dict(updated_unix=time.time(),RUNNABLE_UNKNOWN_WALL_FREE_MATERIAL_IMAGING='PASS_WITH_RESTRICTED_3X2_WALL_7X7_TARGET_SCOPE',DIGITS_ALL_NINE_ACTUALLY_RUN='PASS_B0_B4_AND_FULL_GN_SUPPLEMENT',DIGITS_ACCURATE_MATERIAL_AND_SHAPE='FAIL_0_OF_9_UNDER_FROZEN_IMAGE_GATES',NATIVE_FIXED_GS_CONTRAST_SOURCE='IMPLEMENTED_AND_PHYSICALLY_CHECKED',EXACT_GREEN_ASSIMILATION='PASS_SAME_DISCRETE',GD_SPECTRAL_PROCESSING='IMPLEMENTED_FULL_CURRENT_SPD_PRECONDITIONER;48_MODE_TRUNCATION_FAILS',AUTOMATIC_TARGET_RISK_ABSORPTION='NOT_EXERCISED_IN_NINE_DIGIT_RUNS',REOPEN='SEE_PAIRED_CONTROLLER_ACTION_LOGS;SCHEDULED_AND_AUTOMATIC_DISTINCT',SAME_QUALITY_PHYSICS_COMPUTATION_GAIN='NOT_ESTABLISHED',POSTERIOR_OR_GLOBAL_SAFETY_CERTIFICATE='NOT_ESTABLISHED',UNKNOWN_ARBITRARY_WALL_SUPPORT='NOT_ESTABLISHED',FULL_WAVE_PHYSICS='SAME_DISCRETE_VERIFIED;FURTHER_GRID_0_1_TOLERANCE_2_OF_5_TRUTHS',CONTINUUM_CERTIFICATE='NOT_ESTABLISHED',INDEPENDENT_ALGORITHM_NOVELTY='NOT_ESTABLISHED;STANDARD_CSI_SCHUR_AND_SPECTRAL_PRECONDITIONING',AUTO_PUBLISH='NOT_AUTHORIZED_NOT_DONE',NEXT_RESEARCH='STOP_THIS_ROUND_AND_HAND_RESULTS_TO_WEB_CHAT',counts=summary)
    save('GO_NO_GO.json',nogo)

    write('PHYSICS_VALIDATION.md',f'''# 物理正确性与数值限制

本轮主物理为已有 P3 scalar-TM FEM、非磁介质、exp(−iωt)，共享色散 epsilon=epsilon_inf+i ell/f。没有用 Born 作为主真值或反演前向。

## 实际核验

1. 同离散 j=chi(e_inc+GD0 j)、y=direct+GS0 j 与完整 FEM 使用同一真实积分点和权重。散射场右端 b_W+b_X 保留。完整材料微分包含 db+dV u_sc；B0每源一个加权伴随，Full-GN每频共享16接收伴随。
2. 一次必要的 weighted-adjoint 审核已完成；不把小代数误差当成本轮主科研结果。原生电流审计首次因配置元数据解析失败，费用和日志保留，修复解析器后通过。
3. c010/c011/c012、三频的部分墙吸收、全墙吸收及重新开放实际计算完成。所有背景响应、有效照明、接收返回路径与剩余多次散射均保留。数组在 results/absorption_controls/。
4. 强目标 Born 诱导响应相对误差约 c010=1.811、c011=3.494、c012=6.325。完整反馈不可省略。Born只是独立诊断。
5. 48个GD模态的投影系统内部精确Schur仍不能恢复完整场：三个场景完整测量白化差约20，场相对误差约0.89–0.93。再丢弱块交叉反馈通常更坏。实际B3/B4只用模态作全空间SPD预条件，不删除其余电流。

## 独立进一步网格

原网格h=.17,targetstep=.1,wall_refine=1；进一步网格h=.12,targetstep=.05,wall_refine=2。连续材料盒保持相同，比较完整192复观测的原始白化范数，不除通道数平方根。

{meshtable}

仅2/5真实场景达到旧0.1容差；数字8及其估计终点略超门，高对比c012达0.216。不得将本轮全部新世界写成已通过旧0.1资格，更不能称连续PDE认证。逆问题的同离散结论与这项误差限制分别保存。

所有响应RMS只作每实坐标描述，原始距离另存。噪声每实虚坐标sigma=.01，不曾增大噪声以挽救图像。没有后验覆盖证书、严格inf-sup全域下界或连续Maxwell误差上界。
''')
    write('DIGITS_1_TO_9_FAILURES.md',f'''# 九数字真实运行与失败

九个数字均有独立完整测量、独立未知墙与不同材料反演输出。B0九次、B4九次、Full-GN补充九次全部实际执行。数字ID不进入在线程序，未知变量为全部110个实材料参数。

{summarytable}

{alltable}

## 实际视觉复核

父Sol通过图像工具查看了真实的B0、B4、Full-GN统一色轴总图；详见 evaluation/VISUAL_REVIEW.json。不是只看IoU：B0/Full-GN有部分笔画证据，但1长竖缺失，2斜笔断裂，3/5横笔和侧边不连，4有外侧伪热点，6/9闭环不成立，7主笔偏移，8只留部分下环、两个孔洞未同时恢复。损耗大多接近零或弥散，未跟随真实笔画。B4为弱云团/近空白，不能因质心接近中心就判成像成功。

全部图显示实际7×7材料格，统一阈值chi>=1、统一色轴，没有逐数字调阈值或上采样冒充新增物理分辨率。真值本身也在相同7×7类中，因此本轮失败不能简单归为“真值完全不可表达”；但7×7限制了能研究的真实形状细节。未运行9×13，不能声称更高分辨率必定解决问题。

B0九例均未达冻结一阶停止条件；GN后仅c005局部达标。仍无法据此证明全局MAP、物理唯一或全局不可辨。低训练残差伴随错误材料图是真实局部结果，而不是精确非线性纤维证明。
''')
    controltext=f'''# 方法成对对照与归因

截至本次报告生成，已有{len(controlrows)}个封存并评分的控制终点；缺失行见 METHOD_ABLATIONS_NOT_RUN.csv。

{ctable}

B1是原生固定GS CSI；B1_RESET只在同一时刻重置到物理电流，隔离重置效应。B2加部分材料精确吸收和定时重新开放。B3再加模态Schur预条件。B4为同预条件加经验风险筛选。B5是同样的普通模态Schur预条件、不吸收。因此B1→B5检验预条件，B1→B1_RESET检验物理重置，B3→B5检验相同预条件下的吸收。每族独立冷setup，族内三个公开共同网格允许同样复用；不可只把setup计给某个baseline。

定时REOPEN只证明可逆实现；经验风险真正触发的REOPEN须见逐case动作日志。九数字B4无一次ABSORB/REOPEN，不能把其结果当成自动吸收收益。普通Schur与B3/B4预条件属于同一标准代数；若同样好，不归因于新GSDA机制。

全部控制仍是有限预算迭代，不称全局最优。较低FE但较差图像不构成同质量加速。完整成本在 PHYSICS_COST.json 与各批 BATCH_RECEIPT.json。
'''
    write('CONTROL_RESULTS.md',controltext)

    write('START_HERE.md',f'''# Research W-GSDA-R4：真实九数字全波反演，成像尚未通过

本轮已实现并实际运行未知墙与未知自由材料图的二维全波联合反演。数字1–9全部完成B0、B4与普通Full-GN补充，不是字形生成或标签识别。**运行交付成立，准确材料/形状成像未成立；当前不支持GSDA同质量加速。**

{summarytable}

先看 [九数字直接对照图](DIGITS_1_TO_9_CONTRAST.png)、[损耗图](DIGITS_1_TO_9_LOSS.png)、[完整实验报告](REPORT_ZH.md)、[网页端handoff](WEB_CHAT_HANDOFF.md)。全部图源自保存的材料数组，并经实际视觉检查。

实现正面证据：原生固定GS/GD contrast-source确实运行；精确吸收、背景照明及返回接收路径在同离散完整Maxwell下成立；可逆吸收/重开有实际控制运行。负面证据：九数字B4全部KEEP_ACTIVE，没有发生自动吸收；current数据拟合与材料全波预测存在明显差距；单纯GD模态截断失真。

墙体未知范围是公开3×2容器内12个材料参数；目标为7×7、98个自由实/损耗参数。未使用真墙、真实平均墙、字形模板或有限目标标签。不能推广为任意未知壁形/孔洞数。

复现入口见 [RUN_ALGORITHM.md](RUN_ALGORITHM.md)。原研究封账后拒绝新增物理调用；新部署需显式给独立本地时间额度。只重放已有数组和评分无需重新求解。没有自动commit、push或网页发布。

范围与失败：[GO_NO_GO.json](GO_NO_GO.json)、[PHYSICS_VALIDATION.md](PHYSICS_VALIDATION.md)、[DIGITS_1_TO_9_FAILURES.md](DIGITS_1_TO_9_FAILURES.md)、[CONTROL_RESULTS.md](CONTROL_RESULTS.md)、[PHYSICS_COST.json](PHYSICS_COST.json)。本轮结束后停止，不自动开启下一研究方向。
''')
    write('REPORT_ZH.md',f'''# Research W-GSDA-R4 实验报告

## 1. 结果先行

本轮完成真实二维 scalar-TM、未知墙与自由目标材料图反演：九数字在B0和B4各运行九次，另对九个B0终点做完整Jacobian GN+TV补充。全部选支由测量目标加冻结先验决定，输出后才读真值评分。**没有一个数字通过完整材料/形状工程门。** 此为有限表示、有限预算下的实际负结果，不是目标不存在或Maxwell唯一性的证明。

{summarytable}

## 2. 数据、未知量和开发冻结

全波源接收为四发16收、频率.85/1/1.15，192个复观测，实虚每坐标sigma=.01；全范数无通道归一化。独立留出频率.925不进入在线迭代或选支。墙材料epsilon_inf和ell跨频共享。墙容器3×2，目标7×7，两者所有实对比度与损耗共110实未知。材料box保证非磁、被动。墙未知限于宏格材料，不能宣称一般形状/拓扑。

独立非数字c000开发TV=.2与1，留出预测选择.2；空间3→5→7续接也实际尝试但比逐频自由7×7更差，保留失败。设置在数字反演/评分前冻结。每个数字独立噪声和非均匀墙，前景chi=3.5,ell=.06。c010/c011/c012额外覆盖非凸/空洞/强对比目标；c011墙缺格与外边界连通，不把它称未知封闭墙洞。

## 3. 真正运行的方法

- B0：完整Maxwell、每源一个加权伴随、两个公开均匀起点、频率continuation、box+相同TV/墙properGaussian先验、L-BFGS-B。不是R3无空间先验的弱baseline。
- Full-GN补充：从各B0测量选中的终点开始，不用真值；共享接收伴随构造110列材料J，解析TV曲率、box子问题、LM及完整全波line search。费用按B0+GN累加。数据已暴露，不能重标独立blind确认。只有c005达到局部一阶门，其余20轮上限；未证明全局MAP。
- B1：固定参考GS0、matrix-free GD0的物理积分点contrast-current CSI，材料跨源/频共享；非R3forcing chart。
- B2：实际LA求解精确材料吸收，包括bA、有效入射、直接墙响应与目标经墙返回接收器的路径；定时REOPEN保留恢复状态。
- B3：真实GD actions建模态，强8/弱24 Schur给全电流空间SPD预条件；Q外所有电流仍优化。并非把弱材料坐标锁死。
- B4：B3同预条件，按材料变化、条件曲率代理、完整梯度与GS影响筛选，必要时才付费局部目标风险。九数字全部筛选失败而KEEP_ACTIVE，因此没有执行自动吸收、REOPEN或付费目标风险J。
- B5：相同模态Schur预条件的普通不吸收对照。标准CSI/Schur/SVD机制不改名冒充原创。

## 4. 九数字：运行与成功严格分开

{alltable}

从统一色轴实际看图，B0和GN有部分笔画，但没有完整恢复6/8/9孔洞，其他数字亦有断裂、错位、伪热点。B4材料弱云团不能构成数字。原始物理格保持7×7，不用显示插值声称新增分辨率。详见 evaluation/VISUAL_REVIEW.json。

B0所有选中终点仍未满足1e−4投影梯度门；GN补充改善站定性但未解决形状。不能以未收敛baseline宣布成熟FWI已达上限。B4 current surrogate objective约80而相应材料Full-FEM目标约400–1400，说明电流拟合与共享物理材料一致性之间的优化缺口。最终仍用Full-FEM物理检查点选支，未把不可行电流解当成材料重建。

## 5. 精确吸收、全反馈与谱失效

c010/c011/c012三频实际完成部分墙吸收、全墙吸收及REOPEN。全波源项和所有反馈路径一致；该证据复用标准吸收代数，不是新成像算法收益。Born诱导响应相对误差约181%、349%、632%，真实强反馈不可忽略。

GD48模态投影的完整观测白化误差约20，内部场相对误差近.9，投影内精确Schur不能修复被丢失的外部空间。弱块交叉耦合再删掉往往更差。因此主B3/B4保留全电流空间，仅把Schur作预条件。谱rank下降不是物理计算加速证据。

## 6. 成对消融

{ctable}

完整解释与各控制缺失行在 CONTROL_RESULTS.md、METHOD_ABLATIONS.csv。B1_RESET避免把“重置为物理电流”带来的改善误认成吸收。B3/B5拥有完全相同谱和缓存权限。定时动作与自动风险动作分别计数，未触发的机制不记成功。

## 7. 独立网格及范围

{meshtable}

两真实场景通过旧.1白化容差，三个未通过。所有新图像结论主要为同离散探索，不能声称整个实验达到连续Maxwell认证。逆问题误差、forward离散差、统计误差和优化未收敛单列。

## 8. 成本与公平性

B0九数字实际6,176FE/49,408RHS，约2,932.19 wrapper CPU秒；GN额外1,206FE/13,704RHS，约430.22秒，累计7,382FE/63,112RHS。B4九数字含一次公共空气setup实际381FE/29,580RHS，368.23秒，另有6,480次小dense-Schur LU尝试。公共setup仅在相同公开网格/频率/源接收下复用；各case材料电流与选择状态重新初始化。

这些数字只证明实际工作量不同。B4材料和预测质量更差，不能报同质量20%加速。冷部署包含setup；warm分摊只适用于声明的相同几何输入。Full-GN不能把B0初始化成本消失。全研究还包括开发、独立网格、失败、控制、评分、输出I/O；其最终权威总计在 PHYSICS_COST.json。父子聚合不重复相加，kernel时间不叠加wrapper。

当前已记录全FE预约尝试={totals.get('full_fe_reserved_attempts','未汇总')}，complex RHS={totals.get('complex_rhs_reserved_attempts','未汇总')}；已完成wrapper+未分摊批CPU={number(totals.get('completed_wrapper_cpu_plus_unpartitioned_batch_seconds'),2)}秒。交互式推理/代码编辑CPU未伪造计入。峰值内存见账本，单重型FEM串行。

## 9. 科学裁决

工程交付：可运行、九数字全部实际反演、原生固定GS、电流本构、精确可逆吸收、谱预条件及实际负面数组成立。成像与效率：准确数字材料/形状恢复未成立，风险吸收未在主九数字触发，同质量净加速未成立。未知任意墙支撑、完整posterior和非线性安全证书均未建立。

本轮不开展下一理论研究、不引入NN/Astra/3D、不自动发布。将真实结果和具体障碍交给网页端复盘。标准Schur、CSI及谱预条件的正确实现本身不形成新论文贡献。
''')

    write('QUESTIONS_FOR_WEB_CHAT.md','''# 实验暴露的具体问题（本轮不继续理论研究）

1. 固定GS电流数据项已较小，但对应共享材料Full-FEM目标大得多：状态残差的尺度、材料闭式更新与全波反馈哪一项最阻碍本构一致性？如何保证改进的是物理材料而非不可行电流？
2. B0与Full-GN的低训练残差仍伴随错误笔画；同prior下真值目标函数和局部解目标函数的差异说明什么？如何区分正则偏差、局部极值与实际有限分辨率？
3. 48个GD模态已无法还原场，哪些目标相关方向是普通内部SVD未保留的？只用作SPD预条件为何仍未改善材料一致性？
4. 精确吸收物理上正确，但优化中重置、暂停材料与改变电流条件数的贡献怎样严格分开？B1_RESET/B2/B3/B5实际成对结果应怎样解释？
5. 主九数字风险规则始终KEEP_ACTIVE；是材料确实不够确定、代理过保守，还是筛选成本/尺度本身不合适？何种判据可在不偷读真值时允许有用吸收？
6. 若定时REOPEN恢复物理目标下降但不能恢复孔洞，说明可逆性只是必要条件，还是当前材料空间/信息合同缺少什么？
7. 为何部分笔画和质心可恢复而复材料/损耗与闭环孔洞不可恢复？7×7可表达的真值仍失败，哪种分辨率或正则变化值得作为独立新协议检验？
8. 普通Schur/CSI/Full-GN的哪些能力已经解释全部现象，尚有哪些目标风险调度的独立可证伪增量没有得到证据？
''')

    write('WEB_CHAT_HANDOFF.md',f'''# Research W-GSDA-R4 网页端深度复盘交接

## 本轮事实

九数字1–9均真实生成完整Maxwell测量并完成未知墙与自由材料图反演。B0/B4各9次，普通Full-GN补充9次；不是数字分类、固定标签或真墙Green。所有数字在线输入只有通用case测量与公共合同，输出封存后独立评分。

{summarytable}

## A. 算法到底实现了什么

固定GS0确已实现：参考FEM空气L0、真实积分点j、GS接收伴随、matrix-free GD。材料/局部场/J仍变动。精确吸收用LA及bA，照明、背景数据和接收返回路径完整。GD模态+精确小Schur实际用于SPD预条件，未删除外部电流；另有48模态截断失败实验。B2/B3的REOPEN为预定控制，B4主九数字未吸收/未重开/未支付目标风险J，不能写自动安全吸收已成功。CSI、Schur、SVD、TV、Full-GN属于标准方法。

## B. 正面证据

同离散完整吸收和可逆恢复在三个实际场景/三频通过。Born强反馈差异显著。完整批处理、独立输入/真值评分、逐case材料图/预测/成本可复现。B0/Full-GN部分笔画与质心信息可见；这些局部现象保留，未被低IoU一句话抹去。但没有完整数字通过材料/形状门。

## C. 负面证据与边界

B4速度较快但图像几乎为弱云团，未得到同质量成本优势。current surrogate拟合与Full-FEM材料预测不一致。48模态截断场失真，弱耦合不能随意删。B0原终点0/9达到一阶门，GN后仅c005达标，其余有明确迭代上限；不可用来宣布FWI已全局失败。粗细网格0.1门仅2/5真值通过，c012约.216；没有连续PDE认证。

## D. 九数字完整结果

{alltable}

每行FE/RHS为该方法部署账；B0+GN为累计账。B4逐case是warm，公共setup3FE/636RHS另列BATCH_RECEIPT，总九例381FE/29,580RHS。统一阈值chi>=1和共同色轴，孔洞采用背景4连通，前景同时报4/8连通。实际视觉检查见evaluation/VISUAL_REVIEW.json。

## E. 可直接检查的原始证据

- 主要参数/先验：EXPERIMENT_CONTRACT.json、configs/main.json、configs/current_main.json、METHOD_FREEZE.json、CSI_METHOD_FREEZE.json、GN_RUN_FREEZE.json。
- 独立测量：inputs/online/cNNN/observation.npz，只有y[frequency,source,receiver]；私有真值/heldout在inputs/private/cNNN/truth.npz。
- 结果：results/B0_FULL_FWI、results/B4_GSDA、results/FULL_GN_TV_SUPPLEMENT、results/CONTROL_*；各endpoint.npz保存z/预测/材料；target_material.npy/wall_material.npy为可直接读图的数组；ONLINE_SEALED.json记录选支身份；CSI另有best_currents.npz、trace与动作。
- 物理：results/absorption_controls（全/吸收/重开/目标Born arrays）；results/spectral/GD_MODES.npz和各case.npz；results/independent_validation/*npz与SUMMARY.json。
- 指标：DIGITS_1_TO_9_METRICS.csv、DIGITS_1_TO_9_METRICS_SUMMARY.csv、METHOD_ABLATIONS.csv、BASELINE_COMPARISON.csv、evaluation/SCORES_ALL.json。
- 费用：PHYSICS_LEDGER.jsonl、PHYSICS_COST.json、receipts/、各BATCH_RECEIPT.json。预约失败也计数；kernel不重复加入wrapper。六次控制启动曾误用解析后的系统Python而缺numpy，0PDE的失败原日志保留，随后只修正解释器路径。
- 图：DIGITS_1_TO_9_CONTRAST.png、DIGITS_1_TO_9_LOSS.png、figures/quality_cost/；图全部来自真实数组。
- 实现：src/gsda/、tools/、vendor/。部署命令在RUN_ALGORITHM.md；旧研究截止后拒绝物理求解，新部署需显式独立预算。重放已有结果只需score_all.py。

详细实验解释见REPORT_ZH.md；范围裁决见GO_NO_GO.json。需要新的理论解释的问题单列QUESTIONS_FOR_WEB_CHAT.md。**本轮到此停止，没有自动进入下一科研方向、commit、push或网页发布。**
''')
    write('COST_ACCOUNTING_SCOPE.md','''# 成本口径与已发现的计时缺口

FE/RHS以物理动作之前的追加式reservation为权威，失败尝试也计入。GD/GDH、正向、加权伴随、接收伴随与小dense-Schur分别记录。kernel CPU不再叠加到外层receipt。批总账只补分区未覆盖的CPU，不与所有分区再次相加。

**本轮已有receipt字段wrapper_process_cpu_seconds的精确范围是BudgetBook构造之后，或batch分区开始之后。** 它包含该范围内的网格、装配、LU、求解、日志、I/O、输出和检查，但不包括较早Python启动/依赖导入。科学作业的完整外部CPU不能从已经退出的旧进程恢复，不补猜测值。早期部分只读评分和交互编辑亦未计时。因此全研究所列已记录CPU是已覆盖范围总和，不是机器全部端到端CPU的严格总量。

该缺口在成对控制执行期间发现并保留。后续B0控制启动器如取得RUSAGE_CHILDREN user+system差，则单列whole-child记录，不与已有receipt重复加入。CSI控制外层有真实start/finish墙钟；其启动CPU仍标未测。六次错误解释器的preimport失败为0PDE，但墙钟实测、CPU未测，原日志与命令均保留。

小dense计数覆盖显式记录的模态Schur LU和condition-SVD。其它GN/Gram/投影诊断中的小矩阵求解计入外层CPU，但没有统一的逐调用factorization计数；不能把该字段冒充全部LAPACK调用。FE/RHS全预约账与这一小dense计数范围分别保存。

每方法冷部署必须包含公共参考setup；warm分摊以实际完成case数报告。B4公共setup可同样给普通B5，不属于吸收特有节省。Full-GN从B0终点开始，比较成本必须B0+GN累加。

由于成像质量未通过、旧CPU口径存在启动计时缺口且部分方法未收敛，本轮不提出严格端到端CPU加速结论，也不由较少FE/RHS宣称同质量优势。原始计时值不改写；PHYSICS_COST.json、每job receipt和外部计时文件共同供复核。
''')
    tv0=[x for x in records if x.get('campaign')=='B0_TV0_ABLATION']
    tvrows=[]
    for x in tv0:
        v=x['metrics'];tvrows.append(dict(campaign=x['campaign'],case=x['case_id'],digit=x.get('digit'),prior_TV_weight=0.,IoU=v['support_IoU'],complex_relative_L2=v[metric],centroid_error=v['contrast_centroid_error'],train_white_RMS=v['observation_white_RMS'],heldout_white_RMS=v.get('prediction_heldout_clean_white_RMS'),full_FE=x['cost'].get('full_fe'),RHS=x['cost'].get('rhs'),CPU_seconds=x['cost'].get('wrapper_process_cpu_seconds'),endpoint=x['endpoint_dir']))
    for x in groups.get('B0_FULL_FWI',[]):
        if x['case_id']!='c006':continue
        v=x['metrics'];tvrows.append(dict(campaign=x['campaign'],case=x['case_id'],digit=x.get('digit'),prior_TV_weight=.2,IoU=v['support_IoU'],complex_relative_L2=v[metric],centroid_error=v['contrast_centroid_error'],train_white_RMS=v['observation_white_RMS'],heldout_white_RMS=v.get('prediction_heldout_clean_white_RMS'),full_FE=x['cost'].get('full_fe'),RHS=x['cost'].get('rhs'),CPU_seconds=x['cost'].get('wrapper_process_cpu_seconds'),endpoint=x['endpoint_dir']))
    csvsave('TV0_VS_TV02_C006.csv',tvrows)
    tvtext='\n\n## 独立保留的同数据TV消融\n\n'
    if tv0:
        tvtext+=table(['数字','TV权重','IoU','复材料相对L2','质心误差','训练RMS','留出RMS'],[[x['digit'],x['prior_TV_weight'],number(x['IoU']),number(x['complex_relative_L2']),number(x['centroid_error']),number(x['train_white_RMS']),number(x['heldout_white_RMS'])] for x in tvrows])
        tvtext+='\n\nTV0与TV.2使用同一R4数据、两个公共起点、所有其它先验与完整调用上限。不同TV目标的原始Φ不可直接当成质量胜负；比较共同数据误差、实际图像及材料指标。该消融在数据已暴露后按预声明执行，不能升级为新盲测；不能据跨轮R3/R4差异直接归因。\n'
    else:tvtext+='TV0尚无封存评分输出，NOT_RUN；不得用跨轮R3数字8充当同数据消融。\n'
    for doc in ['REPORT_ZH.md','WEB_CHAT_HANDOFF.md']:
        with (ROOT/doc).open('a') as f:f.write(tvtext)
    paired=read('results/CONTROL_PAIRED_SUMMARY.json',{})
    if paired.get('status')=='ALL_18_COMPLETE_BOUNDED_ATTEMPTS':
        claim='''\n\n## 已完成的关键成对归因\n\n三个父场景中，B4与普通B5的最终参数z及复预测data逐元素完全相同，不仅是目标函数偶然相同；B4为经验诊断每例额外使用144个RHS、合计432个，FE未减少。自动风险控制在这些场景没有产生目标或材料增量。证据：results/CONTROL_PAIRED_SUMMARY.json。\n\nB2相对只重置的B1_RESET，三个全波目标分别降低36.020、29.083、15.373，但额外支付36FE/504RHS；这不是材料图更准确的充分证据。相同模态预条件下，B3相对B5的全波目标变化为+71.297、−21.386、+66.511，额外支付36FE/1080RHS；吸收并不稳定改善优化，更没有显示省物理计算。B2/B3合计真实12次吸收与12次预定REOPEN，自动触发次数为0。\n\n全局完整弱块消元未实现：已实现的是投影核心内精确Schur和全电流空间预条件；48模态截断另有实际失败证据。不能把这写成“弱传播自由度已被安全删除”。\n'''
        for doc in ['CONTROL_RESULTS.md','REPORT_ZH.md','WEB_CHAT_HANDOFF.md']:
            with (ROOT/doc).open('a') as f:f.write(claim)
        nogo['PAIRED_B4_VS_ORDINARY_B5']='IDENTICAL_ENDPOINTS_ALL3;EXTRA432_RHS_NO_FE_REDUCTION'
        nogo['GLOBAL_WEAK_BLOCK_ELIMINATION']='NOT_IMPLEMENTED;PROJECTED_SCHUR_AND_FULLSPACE_PRECONDITIONER_ONLY'
        nogo['REOPEN']='12_ACTUAL_SCHEDULED_REOPENS_B2_B3;AUTOMATIC_TARGET_RISK_REOPEN_NOT_TRIGGERED'
        nogo['CURRENT_ALGORITHM_INCREMENT']='NO_GO_CURRENT_IMPLEMENTATION;3_PAIRED_CONTROLS_PLUS9_DIGIT_PILOT_ONLY'
        nogo['COLD_COST']='NO_GO_CURRENT_IMPLEMENTATION_AT_REQUIRED_IMAGE_QUALITY'
        nogo['BASELINE_FAIRNESS']='TWO_PUBLIC_STARTS_TV_FREQUENCY_CONTINUATION_AND_FULL_GN_SUPPLEMENT;MOST_NONSTATIONARY;NOT_FWI_UPPER_BOUND'
        nogo['VISUAL_REVIEW']='ACTUAL_ALL_NINE_CONTRAST_AND_LOSS_PLUS_CONTROL_AND_TV_FIGURES'
        nogo['REPRODUCIBILITY']='SOURCE_INPUTS_SEALED_ARRAYS_AND_COMMANDS_SAVED;FRESH_INSTALL_AND_NEW_DEPLOYMENT_NOT_RUN'
        save('GO_NO_GO.json',nogo)
    prior_centroid=read('evaluation/PRIOR_ONLY_CENTROID.json',{})
    if prior_centroid:
        pc='\n\n## 质心指标的公开ROI负控制\n\n不查询任何物理、仅在公开目标ROI放均匀正材料，其质心恒为ROI中心。封存后用真实质心评分，九数字质心RMSE为'+number(prior_centroid.get('prior_only_centroid_RMSE'),6)+'λ0；B0为'+number(prior_centroid.get('B0_centroid_RMSE'),6)+'λ0。该零PDE描述性启发式没有字形、真墙或目标信息，却因本目标族大多位于ROI中部而取得较小质心误差。因此小质心误差本身不能证明数据驱动成像，尤其不能把B4弱云团的居中质心当作目标恢复。这里只评几何中心，不冒充proper prior的独立积分/Bayes action，也未伪造prior-only的PDE或预测分数。数据在evaluation/PRIOR_ONLY_CENTROID.json。\n'
        for doc in ['DIGITS_1_TO_9_FAILURES.md','REPORT_ZH.md','WEB_CHAT_HANDOFF.md']:
            with (ROOT/doc).open('a') as f:f.write(pc)
    for doc in ['START_HERE.md','REPORT_ZH.md','WEB_CHAT_HANDOFF.md']:
        with (ROOT/doc).open('a') as f:f.write('\n计时口径补充：已有wrapper CPU从BudgetBook/批分区开始，早于此的Python启动导入未覆盖；没有回填估计。完整说明见 [COST_ACCOUNTING_SCOPE.md](COST_ACCOUNTING_SCOPE.md)，不据此声称严格端到端加速。\n')
    for doc in ['DIGITS_1_TO_9_FAILURES.md','REPORT_ZH.md','WEB_CHAT_HANDOFF.md']:
        with (ROOT/doc).open('a') as f:f.write('\nB4选支补充：九例各实际执行120次交替迭代，但c008（数字3）最终由完整全波目标选择了第0次公开均匀起点；其余八例选择了第30–60次检查点。九份保存材料彼此不同，未复用一张图。不能把“都运行了”写成“每例都优于初始化”。\n')
    with (ROOT/'WEB_CHAT_HANDOFF.md').open('a') as f:f.write('''\n原始Jacobian补充：每个 results/FULL_GN_TV_SUPPLEMENT/cNNN/ 下的 FULL_J_FIRST.npz 与 FULL_J_LAST.npz 保存完整复J、实白化J、参数名/尺度/box、相应材料z、数据、梯度及GN+真实先验曲率；可供网页端后续分析，不需要重新支付PDE。它们仅对应各测量选中候选的局部信息，不能当真值Jacobian或全域posterior。\n\n实际接收器/发射器和PML的已固定历史默认值在RESOLVED_PUBLIC_PHYSICS.json，系事后解析补全出处，并未改变原合同。NPY跨入口格式见DATA_FORMATS.json；推荐用endpoint.npz的分开命名材料数组，避免把epsilon与contrast混读。\n''')
    for doc in ['REPORT_ZH.md','WEB_CHAT_HANDOFF.md']:
        with (ROOT/doc).open('a') as f:f.write('''\n统计单位补充：生成了13个不同物理世界，其中c000为单独开发、其余12个为9数字与3控制；每个只有一份噪声。49个非开发封存终点是同一批12个物理世界的不同方法/补充运行，不是49个独立parent。没有把方法、起点或噪声次数伪装成独立样本，也未据此提供泛化成功率/置信区间。仪器和收发位置已知，本轮没有拟合自由频率gain或注册nuisance。\n\n交付包：delivery/RESEARCH_W_GSDA_R4_20261030.zip保留全部本轮源码、原始数组、传播缓存、图表与账本；delivery/RESEARCH_W_GSDA_R4_WEB_CHAT.zip保留全部选中材料/预测/Jacobian/物理核验数组，但省去可在完整包取得的大型传播缓存，具体遗漏清单在WEB_CHAT_SCOPE.json。精简包不是所有缓存的完整副本；网页端理论复盘可先读WEB_CHAT_HANDOFF.md与QUESTIONS_FOR_WEB_CHAT.md。\n''')
    print(json.dumps(dict(scored_records=len(records),digit_runs={m:len(g) for m,g in groups.items()},controls=len(controlrows),reports_written=True),ensure_ascii=False))

if __name__=='__main__':main()
