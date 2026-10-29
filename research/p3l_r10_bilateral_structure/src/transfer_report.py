"""Read-only postprocessing for the frozen R10 transfer campaign."""
import json,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).absolute().parents[1];O=P/'results/transfer';F=P/'figures'
r=json.loads((O/'RESULTS.json').read_text());receipt=json.loads((O/'RECEIPT.json').read_text());first=json.loads((O/'RECEIPT_attempt1.json').read_text())
groups=r['timing_groups'];conds=r['conditional_results'];colors=['#2379ac','#ed9d2d','#42a37b']
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
for j,g in enumerate(groups):
 m=np.array([a['M'] for a in g['M_comparison']]);ratio=np.array([a['Schur_p50_total']/a['direct_p50_total'] for a in g['M_comparison']]);ax[0].plot(m,ratio,marker='o',lw=1,label=g['parent'][-3:])
ax[0].axhline(1,color='k',lw=.8);ax[0].set_xscale('log');ax[0].set_yscale('log');ax[0].set(xlabel='Amplitude queries per fixed shape',ylabel='Schur / direct cold kernel total',title='Exact256 setup must amortize');ax[0].legend(ncol=2,fontsize=8)
x=np.arange(8);dd=np.array([g['LU_p50']*1e3 for g in groups]);qq=np.array([g['Schur_node_p50']*1e3 for g in groups]);ax[1].bar(x-.19,dd,.38,label='LU256');ax[1].bar(x+.19,qq,.38,label='Schur256 node');ax[1].set(xticks=x,xticklabels=[g['parent'][-3:] for g in groups],ylabel='Kernel median milliseconds',xlabel='Exposed DEV fixed shape',title='30 legitimate amplitudes per shape');ax[1].legend();fig.tight_layout();fig.savefig(F/'transfer_cost.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
for j,c in enumerate(conds):
 a=c['axis_references'];ref=c['independent_references'][-1]['log_evidence'];ax[0].plot([v['core_order'] for v in a],[abs(np.expm1(v['log_evidence']-ref)) for v in a],marker='o',label=c['parent'][-3:],color=colors[j])
ax[0].set_yscale('log');ax[0].set(xlabel='Original-coordinate core GL order',ylabel='Relative evidence gap vs rotated GH12',title='High-SNR quadrature resolution matters');ax[0].axhline(2e-4,color='k',ls='--',lw=.8);ax[0].legend()
c=conds[0];labels=[f"GL{v['core_order']}" for v in c['axis_references']]+[f"GH{v['order']}" for v in c['independent_references']];events=[v['event_probability'][1] for v in c['axis_references']+c['independent_references']];ax[1].bar(labels,events,color='#9b6eb6');ax[1].set(ylim=(0,1),ylabel='Numerical P(Im delta > .105 | fixed shape, y)',title='Parent000 event reference UNQUALIFIED');fig.tight_layout();fig.savefig(F/'transfer_conditional_convergence.png',dpi=180);plt.close(fig)
rows=[]
for c in conds:
 a=c['axis_references'][-1];p=c['profile_material_ROI_mean'];rows.append(dict(parent=c['parent'],status=c['status'],profile_phi=c['profile_phi'],profile_delta_real=c['mode'][0],profile_delta_imag=c['mode'][1],conditional_delta_real=a['amplitude_mean'][0],conditional_delta_imag=a['amplitude_mean'][1],conditional_material_ROI_real=a['material_ROI_mean'][0],conditional_material_ROI_imag=a['material_ROI_mean'][1],log_evidence_GL20=a['log_evidence'],log_evidence_GH12=c['independent_references'][-1]['log_evidence'],P_Re_delta_gt_1p3_GL20=a['event_probability'][0],P_Im_delta_gt_p105_GL20=a['event_probability'][1],evidence_refinement_gap=c['differences']['relative_evidence_axis'],event_refinement_gap=c['differences']['event_axis']))
with (O/'CONDITIONAL_SUMMARY.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
md='''# R10：固定形状后的复振幅传输与真实费用

这一步先问一个具体的物理问题：同一个物体的空间占比形状已经固定，只改变材料的实部、虚部属性时，是否值得反复做一次完整的256阶消元？答案取决于每个形状要查询多少次。一次复杂Schur分解需要预付费用，但之后每个属性点只需求解上三角系统。本轮在8个已经暴露的DEV形状上实际测量这种交换，未建立跨形状统一压缩、Full66后验或新的Maxwell物理定律。

## 精确公式与非厄米边界

原P1材料为 $\\chi=h+\\delta a(z)$，其中 $h=.1+.04i$，$a(z)\\in[0,1]^{256}$ 是原quintic占比，$\\delta=(.8+u_1)+i(.03+.15u_2)$，$u\\sim U([0,1]^2)$。固定主机背景的R9缓存给出 $C,E_R,O,F_h$。定义 $P=\\operatorname{diag}(\\sqrt a)$，$B=PCP$，$L=OP$，$R=PE_R$，则

$$F(\\delta)=F_h+\\delta L(I-\\delta B)^{-1}R.$$

$P$ 必须只使用非负占比的平方根；复属性幅度留在 $\\delta$ 中。即使材料被动，$B$ 仍为一般复非厄米矩阵，不能用Hermitian特征分解。使用复杂Schur $B=QTQ^*$ 后

$$F(\\delta)=F_h+\\delta (LQ)(I-\\delta T)^{-1}(Q^*R).$$

这是同一有限离散原模型的精确代数重排。解析属性导数为 $F'(\\delta)=L(I-\\delta B)^{-2}R$，独立实、虚方向分别为 $F'$ 和 $iF'$；优化阶段的第二次8列求解同样记账。

## 冻结测量与成本解释

配置为 `config/transfer_frozen.json`。形状与观测来自R9 `R9V_FRESH_000..007`；每个形状固定原Gaussian64，30次属性点均位于原始矩形内，两种方法使用完全相同的点、交替先后次序。每次8列求解均事前记账，BLAS1，未重新计算576阶背景。历史R9缓存已支付264 RHS576与1 LU576，这不属于本轮新增费用；新run确实读取缓存，其读取时间另列。

表中cold是“形状算子形成 + 当前方法的形状级预付 + M次节点”的同窗口内核费用，输入读取、保存、Python账本不在内核计时中。完整流程receipt另列，不能把缓存命中宣称为从原始采集起的零成本冷启动。M=1/4/16/64/256/1024是用实际30次节点时间分布推算的摊销表，并未为每个M重跑M次，也没有把历史576阶费用藏入新节点费用。

| DEV形状 | LU节点p50 ms | Schur节点p50 ms | Schur预付 ms | p50 break-even M |
|---|---:|---:|---:|---:|
'''
for g in groups:md+=f"| {g['parent']} | {g['LU_p50']*1e3:.3f} | {g['Schur_node_p50']*1e3:.3f} | {g['Schur_setup_wall']*1e3:.3f} | {g['break_even_M_p50']} |\n"
md+=f'''
![实际代价](figures/transfer_cost.png)

原始逐次范围、条件数估计、残差与完整M表见 `results/transfer/TIMINGS.csv` 和 `RESULTS.json`。240对节点的最大白化端口差为 **{r['max_white_port_error']:.4g}**；这只是有意义的费用/准确度比较中的同模型一致性检查，不作为单独机器精度科学结论。最大1范数条件数估计为 **{r['max_condition1_estimate']:.5g}**。质量未靠低阶近似换取。

## 账本、失败与结论边界

新增总数 **{receipt['RHS256']:,} RHS256 / {receipt['LU256']} LU256 / 0 RHS576**，包含首个失败尝试的 **{first['RHS256']:,} RHS256 / {first['LU256']} LU256**。首跑在旋转坐标积分中误用带符号det作为正Lebesgue体积Jacobian，严格JSON拒绝NaN，未生成合格结果；修正为 $|\\det|$ 后在原总预算内重放。失败源、日志与receipt保存在 `src/transfer_campaign_failed_v1.py`、`results/transfer/attempt1.log`、`RECEIPT_attempt1.json`；账本连续追加，没有清零重记。修正只恢复数学上的正积分测度，未改prior、门槛或物理数据。

两次控制流程wall合计 **{first['wall_seconds']+receipt['wall_seconds']:.3f}s**，CPU合计 **{first['CPU_seconds']+receipt['CPU_seconds']:.3f}s**；峰值RSS **{max(first['peak_RSS_bytes'],receipt['peak_RSS_bytes'])/1e6:.2f}MB**。时间包含数组IO、账本与积分工作，kernel表另有清楚窗口。

结果支持有限范围结论：固定形状重复属性查询时，256阶Schur三角节点比256阶LU快，但预付需要21–58次查询左右才能摊回。单点或少量候选并无冷成本优势；每变一个Gaussian64形状仍要重新形成、分解B。没有执行可选block-Krylov/有理降阶扫描，**经验最小阶数=NOT_RUN；认证最小阶数=未建立**。本轮未把256阶代数重排冒充低阶可控误差模型。
'''
(P/'SHAPE_AMPLITUDE_TRANSFER.md').write_text(md)
prob='''# R10：固定形状条件概率与随机结构的边界

固定一个空间形状后，两个材料属性仍然是随机变量。最优属性点只给出这个形状能有多好地拟合观测；积分还需要计算这一点周围有多少原prior质量。高信噪比使峰宽缩至约0.0015–0.0025；低维并不意味着粗网格自动可靠。本轮真正跑了三个暴露DEV父对象的原始Uniform2积分，保留不合格参考，而不把profile或看似稳定的均值当成后验认证。

## 原测度、噪声与材料量

固定原Gaussian64为暴露anchor值，使用其原观测与相同历史噪声，不读取GT。属性矩形为 $D=[.8,1.8]\\times[.03,.18]$，原密度为 $1/.15$；没有Gaussian化属性、没有新噪声、没有把prior密度替换为材料体积密度。圆对称复噪声的每个实分量方差为 $\\sigma^2/2$，$\\sigma=.00119211798514519$，所以

$$\\Phi(\\delta)=\\|F(\\delta)-y\\|_2^2/\\sigma^2,\\qquad Z_z(y)=\\int_D e^{-\\Phi(\\delta)}\\,d\\delta_Rd\\delta_I/.15.$$

固定形状下材料均值是 $E[\\chi\\mid z,y]=h+a(z)E[\\delta\\mid z,y]$，不是 $\\chi(E[z,\\delta])$ 的一般代替。报告ROI均值时只平均256个ROI单元，不误用576全域或未归一化求和。两个冻结事件为ROI均值实/虚部超过 $h_{R/I}+a_{ROI,mean}(1.3,.105)$，等价于 $\\delta_R>1.3$ 和 $\\delta_I>.105$；它们与形状拓扑事件不同，固定形状没有新的拓扑随机性。

$Z_z$ 可作为固定条件的似然积分；$\\min_\\delta\\Phi$ 是profile。本轮没有对原64维Gaussian形状prior再积分，因此 **不是Full66边缘概率、联合后验、全局证据或总体覆盖率**。形状本身由先前同一观测的anchor优化产生，本轮是暴露条件的计算研究，不是新锁定泛化试验。

## 实际积分与独立检查

先用原属性坐标有界最小二乘和解析Jacobian找profile点；Hessian使用解析gradient的差分计算，正定时用于自适应峰宽。以mode±8边际SD与prior矩形交集为core；原矩形剩余八个矩形同样数值积分，未删除原prior边界或尾部。依次core GL8/12/20、tail GL2/3/5，在原 $\\delta_R,\\delta_I$ 坐标保持正确Lebesgue权重。独立旋转Hessian坐标使用GH8/12，变换有正Jacobian，出prior节点不贡献。这是独立坐标的数值比对，不是区间包络或尾界证书。

冻结判据：相邻原坐标级与跨坐标相对证据差均<2e-4、属性均值差<2e-5、事件概率差<5e-4。三个对象结果如下；数字是候选数值输出，状态均未升级。

| 父对象 | profile Phi | GL20 log Z | GH12 log Z | GL12→20相对证据差 | GL20 P(Imδ>.105) | GH12同事件 | 状态 |
|---|---:|---:|---:|---:|---:|---:|---|
'''
for c in conds:
 a=c['axis_references'][-1];g=c['independent_references'][-1];prob+=f"| {c['parent']} | {c['profile_phi']:.4f} | {a['log_evidence']:.6f} | {g['log_evidence']:.6f} | {c['differences']['relative_evidence_axis']:.4%} | {a['event_probability'][1]:.6f} | {g['event_probability'][1]:.6f} | REFERENCE_UNQUALIFIED |\n"
prob+='''
![收敛与事件失败](figures/transfer_conditional_convergence.png)

三对象的GL20与GH12证据差只有约2.5e-5，均值差小于1e-8，GH8/12证据本身约机器精度一致，但GL12→20的证据变化仍有3.48%，故未通过预定refinement门槛。不能挑较好的一对参考宣称收敛。各tail规则给出的数值质量很小，但规则没有给出峰遗漏/舍入的严格上界，**完整矩形有节点覆盖不等于全支持认证**。

父对象000的虚属性事件阈值穿过窄峰。未按不连续事件边界切开的GL与GH指标积分分别给出约.905与.782，相邻GL级差约.404。这一失败具有直接解释：光滑似然积分的快速稳定，并不自动赋予不连续事件同样的误差控制。需要另外设计在事件边界切开的积分或适当可靠边界规则；本轮未延长预算、改门槛或追加采样强行通过。

完整属性与材料均值、方差、原始节点/权重/Phi、所有profile条件与收敛差保存在 `results/transfer/*conditional.json`、`*integration_nodes.csv`、`*conditional_arrays.npz` 和 `CONDITIONAL_SUMMARY.csv`。非合格数值均值不能作为已认证风险标签。没有模型训练、随机抽样覆盖认证、Full66均值/方差标签，也没有GT评分；原观测只用于条件计算。

## 与随机结构的连接

固定形状的有理传输消除了同一形状下反复LU的费用，却没有积分64维形状。若要让 $Z_z(y)$ 在跨形状概率中有意义，仍须保持原Gaussian64测度、正确计入属性体积因子、验证形状采样/积分与其权重；对每个形状只取最优属性会把profile与marginal混淆。当前正证据是“固定shape多节点查询的算子复用”，当前负证据是“高SNR两维积分及材料事件尚需单独可靠误差控制”。两者都不支持把现有近似Gaussian provider升格为Full66后验。
'''
(P/'PROBABILITY_AND_STOCHASTIC_STRUCTURE.md').write_text(prob)
print('REPORTS_AND_FIGURES_WRITTEN')

# Canonical final readout adds separately frozen revision without deleting V1 failure.
rev=json.loads((O/'RESULTS_revision2.json').read_text());rr=json.loads((O/'RECEIPT_revision2.json').read_text())
finalparents=rev['parents']
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
x=np.arange(3)
for j,c in enumerate(finalparents):
 ds=c['differences'];ax[0].bar(j-.2,ds['relative_evidence_refinement'],.2,color='#e79a32',label='GL4 to GL6' if j==0 else None);ax[0].bar(j,ds['relative_evidence_independent'],.2,color='#368ab1',label='Independent local check' if j==0 else None);ax[0].bar(j+.2,ds['relative_evidence_legacy_GH12'],.2,color='#59a475',label='Legacy GH12 smooth check' if j==0 else None)
ax[0].axhline(2e-4,color='k',ls='--',lw=.8);ax[0].set(xticks=x,xticklabels=['000','001','002'],ylabel='Relative evidence discrepancy',xlabel='Exposed fixed shape',title='REV2 numerical task tolerances pass');ax[0].set_yscale('log');ax[0].legend(fontsize=8)
legacy=conds[0];new=finalparents[0];labels=['V1 GL12','V1 GL20','V1 GH12','REV2 GL4','REV2 GL6','REV2 local'];vals=[legacy['axis_references'][1]['event_probability'][1],legacy['axis_references'][2]['event_probability'][1],legacy['independent_references'][-1]['event_probability'][1],new['axis_references'][0]['event_probability'][1],new['axis_references'][1]['event_probability'][1],new['independent_local_reference']['event_probability'][1]];ax[1].bar(labels,vals,color=['#d0b8d8']*3+['#875299']*3);ax[1].tick_params(axis='x',rotation=25,labelsize=8);ax[1].set(ylim=(0,1),ylabel='P(Im delta > .105 | fixed shape, y)',title='Event boundary split repairs the reference');fig.tight_layout();fig.savefig(F/'transfer_conditional_revision2.png',dpi=180);plt.close(fig)

path=P/'SHAPE_AMPLITUDE_TRANSFER.md';text=path.read_text()
text=text.replace(f"新增总数 **{receipt['RHS256']:,} RHS256 / {receipt['LU256']} LU256 / 0 RHS576**",f"最终新增总数 **{rr['RHS256']:,} RHS256 / {rr['LU256']} LU256 / 0 RHS576**（包括另外冻结的积分修订29,176 RHS256）")
text=text.replace('两次控制流程wall合计', '首失败、修正重放与积分修订三次控制流程wall合计').replace(f"{first['wall_seconds']+receipt['wall_seconds']:.3f}s",f"{first['wall_seconds']+receipt['wall_seconds']+rr['wall_seconds']:.3f}s").replace(f"{first['CPU_seconds']+receipt['CPU_seconds']:.3f}s",f"{first['CPU_seconds']+receipt['CPU_seconds']+rr['CPU_seconds']:.3f}s")
text+='''

## 零新增RHS的稳定性检查

将振幅写入 $K_\\delta=\\delta^{-1}I-B$，有 $F-F_h=LK_\\delta^{-1}R$。若 $\\beta=\\lambda_{min}(\\operatorname{Im}_H B)$ 且 $b_\\delta=\\operatorname{Im}\\delta/|\\delta|^2+\\beta>0$，那么 $-\\operatorname{Im}_H K_\\delta\\succeq b_\\delta I$，从而 $\\|K_\\delta^{-1}\\|_2\\le1/b_\\delta$，不要求矩阵正规。近似状态残差 $r=R-K_\\delta\\widehat u$ 给出 $\\|L(u-\\widehat u)\\|_F\\le\\|L\\|_2\\|r\\|_F/b_\\delta$。

8个缓存B的实际数值特征值显示 **beta为负**（约−4.55e-5到−3.50e-5），不能宣称Im_HB为半正定；合法复属性的被动项仍使全矩形数值下界b至少约.009211。该代数coercivity是有用的残差验收结构，详见 `results/transfer/COERCIVITY_ZERO_RHS.json`；它没有正确质量度量之外的物理解释，也没有浮点特征值的严格区间认证。本轮积分修订用完原重新分配预算的大部分，未追加低阶模型或回退扫描。
'''
path.write_text(text)

path=P/'PROBABILITY_AND_STOCHASTIC_STRUCTURE.md';text=path.read_text();text=text.replace('本轮真正跑了三个暴露DEV父对象的原始Uniform2积分，保留不合格参考，而不把profile或看似稳定的均值当成后验认证。','本轮真正跑了三个暴露DEV父对象的原始Uniform2积分。首版参考不合格；经过事前另行冻结的峰分区与事件边界修订，三者均达到指定任务的数值容差，首版失败完整保留。这仍不是全局后验认证。')
text=text.replace('三个对象结果如下；数字是候选数值输出，状态均未升级。','首版三个对象结果如下；这部分的UNQUALIFIED状态保留，后文给出独立冻结修订。')
text=text.replace('本轮未延长预算、改门槛或追加采样强行通过。','这一首版结果没有被覆盖；后续在根线程事前重新分配的总70,000 RHS256预算内，另行冻结针对积分域与事件边界的算法修订，数值容差保持不变。')
text=text.replace('当前负证据是“高SNR两维积分及材料事件尚需单独可靠误差控制”。','首版负证据是“高SNR两维积分及材料事件需要单独的可靠误差控制”；修订表明，事件边界切分在这个固定形状任务中能有效恢复数值稳定。')
text+='''

## 独立冻结的REV2：峰分区与事件边界修复

根线程在修订结果产生前授权总成本上限70,000 RHS256。第一次节点预检配置 `transfer_revision2_frozen.json` 预计31,480 RHS256，超过剩余29,392，未执行任何节点；保留该预检文件。实际执行配置为 `transfer_revision2_frozen_v2.json`，计划上限3,647节点/29,176 RHS256，已在新物理运行前冻结。

主积分继续保持原属性坐标与完整prior矩形；增加mode、mode±2SD、mode±8SD及**实事件1.3、虚事件.105边界**的分段，core逐块GL4/GL6、tail逐块GL1/GL2。每个矩形内事件恒定，因此不会用未切开的不连续指标乘一个为光滑函数设计的规则。独立检查保持原虚坐标，以协方差三角变换构造条件实坐标GH6，再在虚坐标各峰/事件段做GL6；它只覆盖虚坐标core，明确标作**LOCAL**检查，不能冒充全支持参考。另与旧旋转GH12比较光滑证据/均值；不再用旧未切开的GH指标作为事件真值。

| 父对象 | 分段GL6 log Z | 材料ROI均值Re | 材料ROI均值Im | P(Imδ>.105) | GL4→6证据差 | 事件refinement差 | 独立LOCAL事件差 | 修订状态 |
|---|---:|---:|---:|---:|---:|---:|---:|---|
'''
for c in finalparents:
 a=c['axis_references'][-1];d=c['differences'];text+=f"| {c['parent']} | {a['log_evidence']:.6f} | {a['material_ROI_mean'][0]:.8f} | {a['material_ROI_mean'][1]:.8f} | {a['event_probability'][1]:.8f} | {d['relative_evidence_refinement']:.3g} | {d['event_refinement']:.3g} | {d['event_independent']:.3g} | NUMERICAL_TASK_PASS |\n"
text+=f'''
![另行冻结的修订](figures/transfer_conditional_revision2.png)

三个对象的证据refinement差为3.07e-5–3.80e-5，跨独立LOCAL坐标证据差约6.25e-6，属性均值跨坐标差<1e-9；父对象000的事件概率稳定在约 **.7764753**，refinement差2.40e-5，独立LOCAL差9.39e-8。另两个对象这一事件在数值精度内接近1，不把浮点1解释为严格零补事件质量。

三者的准确状态是 **NUMERICAL_FIXED_SHAPE_TASK_REFERENCE_PASS**，只针对本次固定形状、原Uniform2条件的证据/均值/两个冻结事件达到事前数值门槛。完整矩形的尾块经过实际数值积分，细级数值tail质量约6.90e-42、7.63e-17、2.48e-16；没有严格峰遗漏、舍入或尾部上界，所以 **strict/full-support certificate仍为false**。数值通过与全局认证两项没有混用。

修订新增 **{rr['additional_RHS256']:,} RHS256 / 0 LU256**，wall **{rr['wall_seconds']:.3f}s**；总新增 **{rr['RHS256']:,} RHS256 / {rr['LU256']} LU256 / 0 RHS576**，包括最初失败9,552 RHS与全部重复验证，预算没有清零。`RESULTS_revision2.json`、`RECEIPT_revision2.json`、`*_conditional_revision2.json`、`*_integration_nodes_revision2.csv`、`*_conditional_revision2_arrays.npz` 是修订原始结果；旧文件保留其原状态。当前最小经验压缩阶数、严格尾界、Full66形状积分仍未建立。
'''
path.write_text(text)
print('FINAL_REVISION2_REPORTS_WRITTEN')
finalrows=[]
for c in finalparents:
 a=c['axis_references'][-1];d=c['differences'];finalrows.append(dict(parent=c['parent'],status=c['status'],log_evidence=a['log_evidence'],mean_delta_real=a['amplitude_mean'][0],mean_delta_imag=a['amplitude_mean'][1],material_ROI_mean_real=a['material_ROI_mean'][0],material_ROI_mean_imag=a['material_ROI_mean'][1],P_Re_delta_gt_1p3=a['event_probability'][0],P_Im_delta_gt_p105=a['event_probability'][1],relative_evidence_refinement=d['relative_evidence_refinement'],relative_evidence_independent_LOCAL=d['relative_evidence_independent'],event_refinement=d['event_refinement'],event_independent_LOCAL=d['event_independent'],fullsupport_certificate=False))
with (O/'CONDITIONAL_SUMMARY_revision2.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(finalrows[0]));w.writeheader();w.writerows(finalrows)
finalsummary=dict(status='COMPLETE_WITH_NUMERICAL_TASK_QUALIFICATION',canonical_resources=rr,all_attempts_wall_seconds=first['wall_seconds']+receipt['wall_seconds']+rr['wall_seconds'],all_attempts_CPU_seconds=first['CPU_seconds']+receipt['CPU_seconds']+rr['CPU_seconds'],initial_failure_RHS256=first['RHS256'],initial_failure_included_in_total=True,initial_failure_recount_additional=False,empirical_minimum_compressed_order='NOT_RUN',strict_certified_minimum_order=None,numerical_exact_transfer_order=256,LU_kernel_p50_ms_range=[min(g['LU_p50'] for g in groups)*1e3,max(g['LU_p50'] for g in groups)*1e3],Schur_node_p50_ms_range=[min(g['Schur_node_p50'] for g in groups)*1e3,max(g['Schur_node_p50'] for g in groups)*1e3],Schur_setup_ms_range=[min(g['Schur_setup_wall'] for g in groups)*1e3,max(g['Schur_setup_wall'] for g in groups)*1e3],p50_break_even_M_range=[min(g['break_even_M_p50'] for g in groups),max(g['break_even_M_p50'] for g in groups)],max_white_port_error=r['max_white_port_error'],V1_conditional_statuses=[c['status'] for c in conds],revision2_conditional_statuses=[c['status'] for c in finalparents],fullsupport_certificate=False,Full66_posterior='UNQUALIFIED_NOT_ESTIMATED',GT_access=False,exposure='all8DEV, no new locked generalization',reports=['SHAPE_AMPLITUDE_TRANSFER.md','PROBABILITY_AND_STOCHASTIC_STRUCTURE.md'],source_files=['src/transfer_campaign.py','src/transfer_campaign_failed_v1.py','src/transfer_revision2.py','src/transfer_report.py'],preserved_preflight_rejected_config='config/transfer_revision2_frozen.json',executed_revision2_config='config/transfer_revision2_frozen_v2.json')
(O/'FINAL_SUMMARY.json').write_text(json.dumps(finalsummary,indent=2)+'\n')
print('CANONICAL_FINAL_SUMMARY_WRITTEN')
path=P/'SHAPE_AMPLITUDE_TRANSFER.md';text=path.read_text();text+=f'''

## 复核与复现入口

只读重建两份报告、三个图与汇总表：从项目根目录使用原指定解释器 `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python` 运行 `research/p3l_r10_bilateral_structure/src/transfer_report.py`。它只读取结果，不触发Maxwell节点，也不改变冻结配置。

真实物理复现入口为 `transfer_campaign.py freeze/run` 与 `transfer_revision2.py freeze/run`。当前目录是已封账记录，不能在这里重跑覆盖。完整复现应使用另一个输出目录并先冻结配置、确认只读原R4 prior/background、原R9 host SETUP及R9V暴露anchor/observation路径；仍须单物理worker token与相同BLAS1。V1的带符号Jacobian错误源单独留作审计，不是推荐重现时再次付费执行的算法；想复现历史成本时应保留其失败9552 RHS，而非伪装成科学必需费用。REV2没有新增LU或Schur分解，直接复用前阶段持久化的固定形状Schur数组。报告脚本与原始CSV/NPZ足以独立重算表中数值统计。

成本拆分为：首失败9,552 RHS256/30 LU256；更正后主流程新付 **{receipt['RHS256']-first['RHS256']:,} RHS256/240 LU256**（receipt里40,608是至该阶段的累计值）；事件/峰积分修订新付29,176 RHS256/0 LU256；总69,784 RHS256/270 LU256。成功计时表是8形状×30属性点=240配对节点，失败阶段另有父对象000的30对计时，全部费用进入总账。每个LU前向8列、每个Schur前向8列；属性优化的解析导数额外8列；所有积分节点都是8列。
''';path.write_text(text)
print('REPRODUCTION_AND_COST_SPLIT_ADDED')
