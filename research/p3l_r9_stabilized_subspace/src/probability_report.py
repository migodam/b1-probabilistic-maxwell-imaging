import json,csv,platform,time
from pathlib import Path
P=Path("/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace")
R=P/'results/probability';A=P/'audit/probability'
final=json.loads((R/'final_numerical_comparison_v3.json').read_text());envelopes=json.loads((R/'task_envelope_results.json').read_text());one=json.loads((R/'original1d_results_v2.json').read_text())
receipts=[dict(file=f.name,**json.loads(f.read_text())) for f in R.glob('*_receipt.json')]
total=dict(full_complex_RHS=sum(x['counts'].get('full_complex_RHS',0) for x in receipts),full_LU=sum(x['counts'].get('full_LU',0) for x in receipts),Schur256=sum(x['counts'].get('Schur256',0) for x in receipts),reduced_nodes=sum(x['counts'].get('reduced_nodes',0) for x in receipts),reduced_complex_RHS256=sum(x['counts'].get('reduced_complex_RHS256',0) for x in receipts),job_wall_seconds_sum=sum(x['wall'] for x in receipts),CPU_seconds_sum=sum(x['CPU'] for x in receipts),peak_single_process_RSS_bytes=max(x['peak_RSS_bytes'] for x in receipts),reduced_CPU_cap_seconds=3600,full_RHS_cap=4000,remaining_full_RHS=4000-sum(x['counts'].get('full_complex_RHS',0) for x in receipts),OS=platform.platform(),machine=platform.machine(),threads_per_process=1,history='All R4/R6 assets read-only; no historical solve charged as a new solve; prior/Green construction cost is historical and not invented.',timing_scope='Job wall includes within-job I/O; sum is cumulative job wall and may exceed elapsed because reduced-only jobs briefly overlapped. Module/process imports excluded. All replay receipts included.')
(R/'aggregate_costs.json').write_text(json.dumps(total,indent=2)+'\n')
with (R/'probability_cost_comparison.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['case','physics','order','nodes','full_RHS','reduced_RHS256','integral_wall','task_status'])
 for name,r in final.items():
  for c in r['numerical_old_gate_comparisons']:
   n=c['order']**2;w.writerow([name,c['physics'],c['order'],n,c['cost_counts_delta'].get('full_complex_RHS',0),c['cost_counts_delta'].get('reduced_complex_RHS256',0),c['wall'],c['task_status']])
with (R/'probability_task_envelopes.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['model','cells','mean_halfwidth','event_halfwidth','wall','qualification'])
 for name,rows in envelopes.items():
  for c in rows:w.writerow([name,c['cells'],c['mean_halfwidth'],c['event_halfwidth'],c['wall_seconds'],'NUMERICAL_ONLY'])
 for c in one['envelope_trajectory']:w.writerow(['original_P1_offset_1d',c['cells'],c['mean_halfwidth'],c['event_halfwidth'],c['wall_cumulative'],'NUMERICAL_ONLY'])
failed=[dict(kind='REPORTING_REPAIR',attempt='original1d first run',issue='Floating underflow displayed event upper as zero. This is not zero finite-noise posterior mass.',action='Preserved first-run original1d_results.json/receipt; replayed with numerical upper likelihood floor1e-300 and atom upper bounds.',canonical='original1d_results_v2.json',additional_full_RHS=0,additional_reduced_work_receipt='original1d_underflow_reporting_repair_receipt.json'),dict(kind='QUALIFICATION_REPAIR',attempt='R6 first reference qualification',issue='Event uncertainty initially used only fine-source vs GH; did not include coarse-to-fine source discrepancy. PURE_LAYER final conservative numerical proxy2.551132942e-5, still under.0004.',action='Preserved v1/v2; v3 includes all relevant refinements and explicit non-guaranteed denominator proxy.',canonical='final_numerical_comparison_v3.json',additional_full_RHS=0),dict(kind='CACHE_LIMITATION',attempt='general nonuniform ROI transfer extension',issue='Initial setup persisted only Schur T,Oprime,ERprime,Fh; original Q,C,B,Eh,Op,ER no longer in memory. General nonuniform diagonal update cannot use scalar Schur representation.',action='No duplicate264-RHS setup paid; original nonuniform extension NOT_RUN.',additional_full_RHS=0),dict(kind='METHOD_NO_GO',attempt='R6 task envelope local levels8/16/32',issue='Mean halfwidth remains.00808/.01149 versus1e-5 gate; no strict coefficients/prior moments arithmetic enclosure.',action='Stopped at frozen reduced grid cap; preserved evidence/task envelope trajectories; no new direct or enlarged slice to rescue claim.',additional_full_RHS=0)]
with (A/'failed_attempts.jsonl').open('a') as f:
 for x in failed:f.write(json.dumps(x,ensure_ascii=False)+'\n')
lines=['# R9-S E4：概率积分与精确物理复用分支','',
'本分支在原 R6 的两个冻结 S3 子模型上通过数值物理一致性与旧数值任务门。可靠性等级为 **NUMERICAL_PHYSICS / NUMERICAL_REFERENCE / NUMERICAL_TASK_PASS**。所有 strict 区间字段均为 null；没有 CERTIFIED_CONDITIONAL_PASS。精确传递加速了同一普通求积，但任务中心化包络没有显示相对强普通求积的成本优势。',
'',f'完整新增费用：**{total["full_complex_RHS"]} complex RHS、{total["full_LU"]} 个576阶 LU、1个256阶 complex Schur**；另有 **{total["reduced_nodes"]:,}** 次256阶八源三角求解节点，即 **{total["reduced_complex_RHS256"]:,}** reduced complex RHS。包括两次报告修复/重放与所有参考、失败方法。总进程 CPU {total["CPU_seconds_sum"]:.3f}s，作业 wall 累计 {total["job_wall_seconds_sum"]:.3f}s；单进程峰值 RSS {total["peak_single_process_RSS_bytes"]/1024**2:.1f}MiB。剩余376个完整 RHS未消费。原历史 Green/接收算子的生成费用没有凭空估计，也没有把历史完整账当成新方法速度分母。',
'', '## 1. 冻结合同与物理资格','',
'R6 S3 活动 Gaussianized 索引为 [0,64,65]，独立标准 Gaussian 截断于[-4,4]，其余Gaussian场坐标为零。Q=saturation(1.3 z0)，D=(.8+u64)+i(.03+.15u65)，u=Φ(z)。材料为 χ=.1+.04i+V·1_ROI，V=QD。ROI256、8源、32复接收；源序打包每源Re32后Im32。σc=.00119211798514519；Φ=||F-y||²/σc²，绝不使用Z-ALG的实/虚σ=.001。',
'', '原截断 prior 的 host/pure 平台质量相等，均为 PH=(Φ(-1/1.3)-Φ(-4))/(Φ(4)-Φ(-4))。transition 的材料密度积分原 z0 prior，保留 Q^-2；pure 连续属性密度与transition在材料像重叠处相加，不将源层标签解释为不同材料模式。host是单独原子。',
'', '实际构造 Ah=I-GDh，Eh=Ah^-1Einc，B=Ah^-1GP，C=PᵀB，ER=PᵀEh，O=S(P+DhB)。计算完整 F(V)=Fh+O V(I-VC)^-1ER。非零host的setup计264个完整RHS=8+256；不是rank1。Schur后每节点计8个256阶三角RHS及32×256读出。',
'', '16个事前固定振幅的direct VIE核对：最大白化响应差4.90349e-12（门1e-8）；最大direct backward residual1.66755e-15（门1e-10）。host与basis残差约1.5e-15，Schur残差4.48e-15，非酉残差1.51e-14。普通浮点残差并非整个域的算术证书。旧GH21若干节点的phi差约1e-11，旧数据仍只按历史数值证据读取。',
'', '## 2. 全支持遗漏审查与reference资格','',
'每例先从统一起点(.9,.09)用真实reduced响应重新定位mode/GN尺度；不读取生成材料或把旧GH oracle免费交给方法。保留矩形D为mode各材料坐标±12个GN边缘SD。完整合法V支持包含于[0,RHI]×[0,IHI]，不是删除D之外的prior。D外八个初始矩形递归二分，所有单元均保留；d_i≤0的单元继续分割，未直接删除。当前两例没有未闭合单元。',
'', r'对单元中心a、外接圆半径r，d_i=|1-aT_ii|-r|T_ii|>0。vmax=|a|+r。对|E′|连续两次使用 x_i=(b_i+vmax∑_{j>i}|T_ij|x_j)/d_i 得到导数的三角majorant，继而 e=√2 r || |O′|x ||F/σc。由于 F′(v)=O′(I-vT)^-2E′，整个单元白化响应距离≤e。故 Φ(v)≥.5·max(0,||r_a||-e)²。这是精确代数下的整个单元结论；本程序使用普通浮点系数/运算，所以当前只是数值全域审查。',
'', '只在单元Φ下界比mode至少高50时排除。所有被排除单元的并集，其原prior总质量≤1，故缩放遗漏evidence≤max_cell exp[-(Φ_lower−Φ_mode)]；没有给每个cell分别赋质量1再错误相加。host原子与pure/transition的像全部在相同外矩形中；host也在source积分中显式保留。D内的独立source积分在(z0,u64,u65)原坐标执行，保留截断归一化、源层、相位事件与host原子；与正确Q^-2的GH21/31材料积分交叉。',
'', 'GH15/21本身只局部稳定，本报告没有仅凭它升级reference。GH21/31、原source GL12×20与20×32细化、全支持审查共同提供数值reference。数值不确定性代理取跨坐标与相邻细化的最大任务差；不是证明半宽。遗漏posterior代理的分母取所有source/GH缩放evidence最小值减其范围的两倍，保证本组数值为正；这个保守有限细化代理也不是已证明的evidence下界。',
'', '|案例|审查cells / 排除leaves|均值不确定性代理|代理/材料SD|源层/事件代理|遗漏posterior代理|', '|---|---:|---:|---:|---:|---:|']
for name,r in final.items():
 q=r['reference_qualification'];au=r['whole_support_audit'];lines.append(f'|{name}|{au["evaluated_cells"]}/{au["excluded_cells"]}|{q["mean_absolute_numeric_uncertainty_proxy"]:.6g}|{q["mean_uncertainty_over_SD"]:.6g}|{q["event_or_layer_numeric_uncertainty_proxy"]:.6g}|{q["numerical_omitted_posterior_upper_proxy"]:.6g}|')
lines += ['', '旧统计门仍为均值最大绝对误差≤1e-5且≤.01SD，概率/源层≤.002；reference代理采用再紧五倍的2e-6、.002SD、.0004。两例数值代理均通过。source粗→细evidence相对变化仍为约.58%/1.31%；最终source对GH31差为1.53e-6/3.78e-6。未把有限细化误差代理当严格界，也没有隐藏normalizer的较慢收敛。', '', '## 3. 强普通求积与物理复用二因子','', '|案例/规则|direct积分wall/s|Schur积分wall/s|Schur setup+积分/s|旧任务门|', '|---|---:|---:|---:|---|']
for name,r in final.items():
 for order,x in r['two_factor_physical_reuse'].items():lines.append(f'|{name} GH{order}|{x["direct_integral_wall"]:.6f}|{x["shared_integral_wall"]:.6f}|{x["shared_setup_plus_integral_wall"]:.6f}|NUMERICAL_TASK_PASS|')
lines += ['', 'direct与Schur执行完全相同的GN适应GH9/11节点、原prior fibre、任务与固定观测；包括reference代理余量后均通过旧任务门。这把物理复用收益与概率规则分开：单次记录积分wall比约5.7–9.8；计入一次setup的组件比约3.1–5.7。它们是单次运行的组件计时，**不是完整冷启动或统计速度结论**；共同fresh profile、全部验证/reference、读文件与imports不在该小表中，完整作业账另外保存。两种物理使用同一概率规则，因此此处的收益归于成熟精确消元/Schur复用。',
'', '完整reference用160k级reduced节点，尤其原source条件积分昂贵；没有把它当免费posterior oracle。最强普通GN适应求积GH9/11已经很便宜。全支持任务包络未在direct物理下重新全跑：它的reduced均值门尚未闭合，没有理由花最后376个RHS伪造完整二因子胜利。因此“概率方法×物理实现”的矩阵只完成普通正确quotient求积这一行，不能声称所有概率规则都有已执行direct对照。',
'', 'REGULAR移除material fibre prior会使均值移动.111977SD，违反旧均值门；PURE_OVERLAP正确源层质量为.56889634524/.43110365476。旧global source/quotient输出另做零新增物理的读取重现：即使prior归一化约1仍漏窄likelihood，REGULAR最细规则误差约.524/3.585SD。它是历史失败机制的复核，未包装为新发现。',
'', '## 4. 任务中心化包络与六类干预','',
'R6在D内用8²/16²/32²cell的R8P三角物理包络；prior正负任务矩采用GL3 cell节点与原fibre GL96，资格为NUMERICAL_PRIOR_MOMENTS。外域保留前述数值遗漏evidence代理。严格物理/系数/prior矩/舍入误差未封闭，所有strict字段null。',
'', '|案例/32²cells|均值数值半宽|256-ROI事件数值半宽|结果|', '|---|---:|---:|---|']
for name,rows in envelopes.items():
 c=rows[-1];lines.append(f'|{name}|{c["mean_halfwidth"]:.6g}|{c["event_halfwidth"]:.6g}|均值旧门未闭合；事件数值门通过|')
lines += ['', '该事件在均匀材料子模型中等价Re(V)≥.65，两例峰位远离阈值；容易事件不能作为一般不连续任务成本优势的证据。当前包络均值比1e-5门宽约三个数量级，不宣称达同精度更快。',
'', '|干预|实际执行与解释|', '|---|---|', '|忽略弱方向/体积|移除material fibre prior的flat-material错误对照，REGULAR偏移.111977SD；不是正确posterior。|', '|prior conditional积分|保留Q^-2、pure重叠源贡献与正确prior，适用于同材料gauge消元。|', '|真实likelihood条件积分|原(z0,u64,u65)source细化，独立验证正确quotient积分。|', '|便宜响应/局部近似|后述1D条件模型的匹配Laplace仅作局部近似；未用profile零残差赋posterior质量。|', '|提高物理精度|direct VIE与精确Schur在固定验证/同求积节点交叉；没有冒称任意精度或严格舍入已完成。|', '|探索另一合法区域|穷尽同一幅度支持的矩形审查；没有以候选数或优化重复收敛替代外域质量。|',
'', '六类干预不都是共享同一初始状态的完整最优调度试验；当前是可复核的组件对照与停止轨迹，不能声称已证明新策略最优。',
'', '## 5. 原P1一维固定正交条件模型','',
'另冻结原Gaussian坐标e0（全ROI均匀offset），其余63个Gaussian场坐标为0，两个Gaussianized属性为0，即native属性均为.5。D=1.3+.105i；z0服从**原始未截断标准Gaussian**。观测固定使用已公开R6 REGULAR controlled observation，不按slice拟合程度更换。此模型仍处于已验证均匀幅度族；不是新非均匀P1实验，更不是Full66材料边缘。',
'', f'两个无限Gaussian尾部准确推为host/pure原子，各质量{one["prior_host_atom"]:.14f}，transition质量{one["prior_transition"]:.14f}，没有删除tail或另立截断prior。原z0加权adaptive quad与t=Φ(z0)的uniform-prior坐标adaptive quad一致：mean={one["references"][0]["mean_task"]:.14f}，variance={one["references"][0]["material_variance"]:.12g}，mean差约5.6e-17。数值全支持任务包络首次在512cells达到均值1e-3门，半宽2.69220e-4；32cells已达到事件.002数值门。这展示了同一y下不同任务可在不同时间停止，posterior variance没有因多做计算而改变。',
'', '不过事件非常稀有。首次运行在浮点下溢后把数值上界写成0；原文件保留，追加重放用upper likelihood floor1e-300与atom upper bound，给出非零上界。512cells事件数值区间为[0,6.45969e-298]。这只是当前浮点矩阵与包络模型的数值资格；不把数值0当真正零概率，也不宣称极小上界的最后位数是严格精度。整个重放成本计入总账。',
'', '同模型强普通adaptive每个完整四矩/事件路线约.4s；任务包络到512cells约1.57s，未显示净成本优势。局部Laplace均值约.218911626，和reference差约5.20e-6，满足宽松1e-3数值均值门，但它没有完整尾域/余项资格，不能由此认证事件。',
'', '## 6. 保留的缺口与可复现入口','',
'- **严格证书缺口**：Schur、基列、浮点系数、Gaussian CDF/任务矩以及指数/舍入没有统一外向包围。数值全支持审查与数值任务通过不能解除这个缺口。',
'- **original nonuniform条件模型 NOT_RUN**：当前缓存未保存原Q/C/B/Eh/Op/ER，无法直接承接非均匀diag材料更新；没有重复消费264 RHS来扩展。',
'- **Full66_UNQUALIFIED**：一维/三维条件成功不提供Full66边缘或全域质量。',
'- **概率方法优势未建立**：物理复用有效，强普通求积已很快；task envelope在R6均值旧门未闭合，1D也未赢wall。',
'- **完整冷启动比较仍缺**：小表只是匹配组件；完整qualification/历史冷启动成本有单独账与明确未知项，不能以组件比冒充端到端算法速度。',
'', 'Canonical evidence: `results/probability/final_numerical_comparison_v3.json`, `original1d_results_v2.json`, `task_envelope_results.json`, `aggregate_costs.json`。输入与冻结：`config/probability_frozen.json`、`probability_reduced_addendum.json`、`probability_original_1d_frozen.json`；追加账：`audit/probability/ledger.jsonl`；失败/修复：`failed_attempts.jsonl`。详细excluded cells含anchor、白化半径、diagonal margin和phi gap。',
'', '复现应在新目录/ledger中执行，脚本检测已存在receipt并拒绝静默重放。历史R4/R6目录只读。本分支没有网络、训练、GPU或新增非对称图册。']
(P/'E4_PROBABILITY_REPORT_ZH.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(total))
