"""Render the completed E0 arrays and write a self-contained research report."""
from pathlib import Path
import json,time
import numpy as np
from src.b1.inference import material_image
from src.b1.e0 import dump

def write(root):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root=Path(root);stage=root/'research/b1_next/e0';out=stage/'results';figures=stage/'figures';figures.mkdir(exist_ok=True)
    started=time.perf_counter();cpu=time.process_time()
    summary=json.loads((out/'SUMMARY.json').read_text());controls=json.loads((out/'CONTROLS.json').read_text())
    probes=json.loads((out/'CONTINUOUS_CDF_PROBES.json').read_text());config=json.loads((root/'configs/protocol.json').read_text())
    sensitivity=json.loads((out/'PRIOR_AND_NUISANCE_SENSITIVITY.json').read_text())
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
    names=['solid','central cavity','central + annular cavities'];images={}
    for i,(name,case) in enumerate(summary['cases'].items(),1):
        fig=plt.figure(figsize=(13,12),layout='constrained');gs=fig.add_gridspec(3,2,height_ratios=[1,1,.48])
        for row,item in enumerate(case['top_k']['items']):
            h=item['H'];mu=item['representative_material'];sd=item['conditional_uncertainty']['sd'];p=item['posterior_mass']
            axis,image,_=material_image(h,mu,config,151)
            intervals=config['topologies'][h]['void_intervals'];odds=item['posterior_odds_vs_rest']
            odds_label='rounded dominant event' if odds is None else f'{odds:.5g}'
            for col,values in enumerate([image.real,image.imag]):
                ax=fig.add_subplot(gs[row,col]);im=ax.imshow(values,origin='lower',extent=[axis[0],axis[-1],axis[0],axis[-1]],
                    cmap='viridis' if col==0 else 'magma',vmin=1 if col==0 else 0,vmax=2.8 if col==0 else .16)
                ax.set_title(f'H{h}: {names[h]} | '+('Re epsilon' if col==0 else 'Im epsilon')+
                  f'\nmass {100*p:.6g}%; integration delta <= {100*item["mass_error_estimate"]:.2e} pp'+
                  f'\nmean eps=({mu[0]:.6f}, {mu[1]:.6f})'+
                  f'\nSD (Re,Im)=({sd[0]:.5f}, {sd[1]:.5f}); model error unbounded',fontsize=9)
                label='x / R'
                if col==0:label+=f'\nVacuum radial intervals: {intervals}\nPrior odds vs rest=0.5; posterior odds={odds_label}'
                ax.set_xlabel(label,fontsize=8);ax.set_ylabel('y / R');fig.colorbar(im,ax=ax,shrink=.72)
        ax=fig.add_subplot(gs[2,:]);p=np.array(case['result']['topology_probabilities']);ax.bar(names,p,color=['#4c78a8','#f58518','#54a24b'])
        ax.set_ylim(0,1.12);ax.set_ylabel('whole-event mass')
        for j,v in enumerate(p):ax.text(j,v+.025,f'{v:.8g}',ha='center')
        ax.set_title(f'Top-2 modeled remainder={case["top_k"]["remaining_modeled_mass"]:.8g}; broader geometry mass: unquantified',fontsize=10)
        fig.suptitle(f'B1 E0 | {name}\n3D vector-Mie frozen polynomial target; conditional-mean material sections\n'
                     'Exact image mass=0; empirical numerical error only; representative not newly full-forward verified',fontsize=12)
        file=f'fig{i:02d}_{name}_topk.png';fig.savefig(figures/file,dpi=140);plt.close(fig);images[name]=file
    strong=summary['cases']['distinguishable'];h0=strong['result']['conditional_material'][0]
    curve=np.load(out/'cases/distinguishable_curves.npz');x=curve['poly7_97_physical_grid'][0];y=curve['poly7_97_conditional_cdf'][0,0]
    fig,ax=plt.subplots(figsize=(10,5),layout='constrained');ax.plot(x,y,label='same-target continuous CDF, tensor97',linewidth=2)
    ax.scatter(curve['adaptive_poly7_physical_grid'][0],curve['adaptive_poly7_conditional_cdf'][0,0],s=16,label='independent partitioned adaptive CDF')
    chosen=[r for r in probes['rows'] if r['case']=='distinguishable' and r['H']==0 and r['coordinate']==0 and r['saved_claimed_quantile']==.05]
    for index,row in enumerate(chosen):
        ax.scatter([row['cutoff']],[row['continuous_target_CDF']],s=60)
        position=(2.386,.34) if index==0 else (2.51,.16)
        ax.annotate(f'{row["source"]}: reported q05\nactual CDF={row["continuous_target_CDF"]:.5f}',
                    (row['cutoff'],row['continuous_target_CDF']),xytext=position,
                    fontsize=9,arrowprops={'arrowstyle':'->'})
    q=h0['quantiles_05_50_95'][0];ax.scatter(q,[.05,.5,.95],marker='x',s=55,label='direct-integral quantiles')
    ax.axhline(.05,color='gray',linestyle=':');ax.set(xlabel='Re epsilon',ylabel='conditional CDF',xlim=(2.38,2.64),ylim=(-.015,1.02),
           title='Strong-data H0: correcting CDF approximation without changing the prior or likelihood')
    ax.legend(loc='lower right',fontsize=8);fig.savefig(figures/'fig07_strong_material_CDF.png',dpi=150);plt.close(fig)
    native=strong['legacy_original']['conditional_material'][0];old=strong['legacy_refined']['conditional_material'][0]
    sd=[np.sqrt(native['covariance'][0][0]),np.sqrt(old['covariance'][0][0]),h0['sd'][0],strong['integration_results']['adaptive_poly7']['conditional_material'][0]['sd'][0]]
    fig,ax=plt.subplots(figsize=(10,5),layout='constrained');ax.bar(['legacy full-grid7\nnode moments','legacy polynomial25\nnode moments','continuous tensor97','independent adaptive'],sd)
    ax.set_ylabel('conditional SD of Re epsilon');ax.set_title('Topology mass saturation does not certify continuous material uncertainty')
    for j,v in enumerate(sd):ax.text(j,v+.0005,f'{v:.9f}',ha='center',fontsize=10)
    ax.set_ylim(0,.033);fig.savefig(figures/'fig08_material_SD.png',dpi=150);plt.close(fig)
    rows=[];references=[];methods=[]
    for name,case in summary['cases'].items():
        p=case['result']['topology_probabilities'];rows.append(f'|{name}|{100*p[0]:.8g}|{100*p[1]:.8g}|{100*p[2]:.8g}|{case["status"]}|')
        gl=case['comparisons']['poly7_65_vs_97'];ad=case['comparisons']['poly7_97_vs_adaptive'];model=case['comparisons']['response_poly5_vs_poly7']
        maxi=lambda r:max(max(q['cdf_tested_grid_max']) for q in r['rows'])
        references.append(f'|{name}|{maxi(gl):.3e}|{maxi(ad):.3e}|{maxi(model):.3e}|{ad["max_event_mass_difference"]:.3e}|')
        for method,r in case['integration_results'].items():
            methods.append(f'|{name}|{method}|{r["wall_seconds"]:.6f}|{r["conditional_material"][0]["sd"][0]:.9g}|{r["conditional_material"][0]["quantiles_05_50_95"][0][0]:.9g}|')
    stats=controls['checks']['saved_complete_case_moment_quantile_coherence'];bad=[r for r in stats['rows'] if r.get('violations')]
    failures='\n'.join(f'- {r["source"]} / {r["case"]} / H{r["H"]} / {r["coordinate"]}: {r["violations"]}' for r in bad)
    probe_native=next(r for r in chosen if r['source']=='legacy_original');probe_refined=next(r for r in chosen if r['source']=='legacy_refined')
    resources=json.loads((out/'RESOURCE_SUMMARY.json').read_text())
    def md(name,text):(stage/name).write_text(text+'\n')
    table='|案例|H0实心%|H1中心空腔%|H2中心＋环形空腔%|连续摘要状态|\n|---|---:|---:|---:|---|\n'+'\n'.join(rows)
    ref_table='|案例|65/97 CDF最大差|97/自适应CDF最大差|5×5/7×7 response的CDF差|独立事件质量最大差|\n|---|---:|---:|---:|---:|\n'+'\n'.join(references)
    report=f'''# B1 E0：让材料均值、方差和分位数属于同一个后验

2026-10-13。该阶段新增物理预算与实际调用均为0，工作位于独立分支`research/b1-e0-posterior-audit`。历史384次账本保留，P1的4,000次只是未授权建议。

## 从基础问题到发现

材料事件的概率是likelihood在整个事件及连续材料参数上的积分；一张连续材料图的精确概率质量为0。传统MAP与最低残差不能代替该积分。同样，稳定的拓扑质量不能保证事件内部的材料不确定性正确。

B1历史模型是三个固定同心模板、连续复介电常数和真实3D矢量Mie响应。主prior为Re ε均匀[2.2,2.8]、Im ε均匀[.04,.16]，拓扑prior各1/3。当前阶段是这个已支付benchmark的数值修复；未把三个模板升级为开放材料生成架构，也没有强制下一阶段按孔洞数分设生成专家。

强信号distinguishable/H0的原生节点方差为{native['covariance'][0][0]:.12g}，但原分箱q05/q95与均值要求至少约3.1851109e-4。约98.64倍的违反来自两套近似摘要混用，不能归因于evidence归一化错误。72个历史条件坐标记录中检测到{stats['incompatible_row_count']}项必要条件不相容；没有违反也不是准确性认证。

## 明确的posterior target与修复

本阶段使用已经付费的7×7 response网格对应的tensor-Legendre多项式；每坐标最高degree6。5×5对照最高degree4。它们仍是明确的近似forward，不是新增full-Mie求解。

对每个H，Z_H=∫L(F_poly(H,θ);y)dμ_H，p_H=π_H Z_H/Σπ_G Z_G。所有moments、协方差、CDF和quantiles来自同一个连续density。CDF在截断prior矩形上直接积分，quantile通过反解CDF求得，不再将节点质量铺到prior bins。Gauss权重是连续prior积分权重，不是离散候选prior。

实化采用全部16个复通道Re后接Im；σ是固定的每实分量SD，复RMS为√2σ。g~N(1,s_g²)的边缘协方差为σ²I+s_g²FFᵀ，保留candidate-dependent logdet。保存与重放显式携带gain_sd；新NPZ入口拒绝缺失该字段时静默变为known gain。

## 决定性的连续材料结果

|strong H0 / Re ε|均值|SD|q05|q95|
|---|---:|---:|---:|---:|
|历史原生节点/bin混合|{native['mean'][0]:.12f}|{np.sqrt(native['covariance'][0][0]):.12f}|{native['quantiles_05_50_95'][0][0]:.12f}|{native['quantiles_05_50_95'][0][2]:.12f}|
|历史经验细化节点/bin混合|{old['mean'][0]:.12f}|{np.sqrt(old['covariance'][0][0]):.12f}|{old['quantiles_05_50_95'][0][0]:.12f}|{old['quantiles_05_50_95'][0][2]:.12f}|
|同一连续density积分97|{h0['mean'][0]:.12f}|{h0['sd'][0]:.12f}|{q[0]:.12f}|{q[2]:.12f}|

在修复target下，原生所报q05实际CDF为{probe_native['continuous_target_CDF']:.12f}，细化所报q05实际CDF为{probe_refined['continuous_target_CDF']:.12f}，分别不是5%。这些值来自新的直接截断积分，不是绘图插值。细化moments已经接近收敛，但它的bin-CDF仍有可检测偏差。

![直接CDF核对](figures/fig07_strong_material_CDF.png)

![材料SD比较](figures/fig08_material_SD.png)

## 独立数值参考与误差边界

主积分按25/41/65/97阶细化，覆盖六案例、三个条件事件和材料边际。独立算法用二维分区自适应5/9点积分，在不同partition节点上重新计算同一density。主CDF网格129点，独立CDF33点并直接反解quantiles。CDF比较只针对已测网格，不能宣称全域sup误差证书。

预定容差为CDF差≤.01、均值/quantile差≤prior宽度×.001、SD差≤max(参考SD×.01,prior宽度×1e-6)，必要variance下界浮点容差为prior宽度²×1e-12。全部六例通过本阶段连续摘要的经验核对，未修改容差。response插值5×5/7×7的差异与likelihood积分误差分列：

{ref_table}

独立积分与相同polynomial target的一致支持数值积分正确性，不证明整个prior上full-Mie插值误差、连续Maxwell极限、漏模式质量或真实仪器模型正确。

## Top-K成像与概率语义

{table}

图的材料代表为事件内条件均值ε，在固定模板内合法，但未新增full-forward代表核验。图示为3D材料的2D中心切面。每项附prior/posterior odds、条件SD、数值mass误差及事件定义；Top2保留第三事件余量，不重新归一化。prior外质量0只描述这个封闭prior；更大几何世界未定义，质量保持null/未量化。极小质量的第二模板不能称为发现了第二个重要mode。

![竞争案例的两个材料解释](figures/{images['competitive']})

## prior、nuisance与必要负控制

信息不足数据下，改变拓扑prior显著改变质量；计算不因此矛盾。改变prior的结果拥有独立protocol标识，不能混作原target。相同gain数据分别用known gain与proper marginalized gain处理，likelihood协议不同；前者集中于错误H2，后者H1/H2竞争。补偿profile minimum仍不能替代integrated evidence。

附件零物理代数重放以及Cholesky独立核对通过；保留slogdet的除零/overflow/invalid警告。高ESS漏模式反例、prior-only posterior、同测度坐标变换、radius/volume不同prior、event partition、metadata回放和TopK余量均执行。toy和随机矩阵不作为新增Maxwell物理证据。历史材料coverage使用旧bin-CDF，本阶段没有重算720个data的coverage，不把旧数值沿用为修复后的coverage声明。

## 成本、失败与可复现性

主推断与controls实际body wall={resources['wall_seconds_inference_and_controls']:.6f}s，CPU={resources['cpu_seconds_inference_and_controls']:.6f}s，实际polynomial evaluations={resources['polynomial_counts']['polynomial_response_evaluations']:,}。这些不是物理调用。直接CDF probes和public CLI replay另计，绘图另计；网络、模型服务、交互时间不完整测量。没有新GPU/NN/旧locked数据读取，也没有SHA-256或额外smoke。

四个历史treams oracle数组仍缺失：保存了输入、误差和环境，但不能从旧原始oracle数组重算。后续refine源码保存数组不补齐历史。没有重跑来掩盖缺口。三次原生报告worker尝试因模型容量失败，父Codex接手，不切换provider/model，服务成本未测。全部记录见FAILURES及资源文件。

## 科研裁决与明确停止

继续保留并修复概率benchmark。E0证实了一个实用反例：拓扑质量饱和与连续材料摘要失真可并存；这属于数值正确性发现，不是新Maxwell可辨识定理。

E1 raw/full-Mie secant、ε及ε+gain补偿、多频、统一implicit radial field和算法机制比较全部NOT_RUN：没有新增物理授权。nuisance-aware nonlinear secant的新增解释力、物理proposal全成本优势、active sensing收益都尚未测试；本阶段不能选择它们为已验证赢家。4,000预算不自动开启，不能借用C1余额。

当前没有IEEE TAP强理论/算法创新证据。停止用高ESS、低residual、拓扑mass稳定或漂亮图替代连续posterior参考；停止把fixed-template benchmark说成统一开放几何成像。下一步在有真实授权时先验证结构信号与nuisance compensation；机制失败则转active sensing或仅保留benchmark，不能无限调噪声制造成功。
'''
    md('RESEARCH_REPORT_ZH.md',report)
    md('START_HERE.md','''# B1 E0 — 零新增物理的连续posterior修复

[中文研究报告](RESEARCH_REPORT_ZH.md)、[参考审计](REFERENCE_AUDIT.md)、[理论状态](THEORY_STATUS.md)、[方法对照](METHOD_COMPARISON.md)、[prior/nuisance](PRIOR_AND_NUISANCE_SENSITIVITY.md)、[资源](RESOURCE_SUMMARY.md)、[失败](FAILURES.md)、[裁决](GATE_DECISION.md)、[复现](README.md)。

真实结果：原生节点moments与prior-bin quantiles不相容；新摘要对同一个连续polynomial likelihood积分。全部六案例完成独立数值核对，仍是empirical reference。新增物理0，旧384/384不变；E1未授权、未运行，新阶段未公开发布。

数据入口：[SUMMARY.json](results/SUMMARY.json)、[controls](results/CONTROLS.json)、[直接CDF probes](results/CONTINUOUS_CDF_PROBES.json)、[实际CLI回放](results/PUBLIC_INTERFACE_REPLAY.json)。

![竞争材料解释](figures/fig02_competitive_topk.png)
''')
    md('THEORY_STATUS.md','''# 理论状态

已核对：proper prior上的evidence、Gaussian实/复噪声约定、gain解析边缘化与logdet、同测度坐标变换、事件分区与TopK余量；节点moments与bin-CDF混合的coherence反例。必要方差下界仅为必要条件，不是收敛充分条件。

附件finite material secant、receiver-dual、likelihood tilt、Gaussian coupling和A17 affine/remainder恒等式仅做零物理随机矩阵/代数核验。Maxwell小空腔系数、余项、nuisance后剩余信号、full posterior odds解释力、solver稳定常数、一般三维支持未测。

本阶段新贡献是明确数值审计及修复，未申报基础Bayes/积分恒等式为原创。只有明确的同一numerical target被检查；没有continuum或现实概率认证。见完整报告。
''')
    md('REFERENCE_AUDIT.md','# 连续参考审计\n\n'+ref_table+'\n\n使用129点主CDF、33点独立CDF及直接quantile反解；全部网格差为经验差。主积分25/41/65/97，独立二维分区5/9。5×5/7×7 response差单独记录。四项历史原生full-grid失败保留；四个历史treams原始oracle数组缺失保留。\n\n详细逐事件/坐标、moments、SD、quantile、evidence、CDF误差及adaptive收敛状态见results/SUMMARY.json和cases/*_curves.npz。')
    prior_lines=[]
    for name,entries in sensitivity['prior_sensitivity'].items():
        for r in entries:prior_lines.append(f'|{name}|{r["protocol"]}|{r["topology_prior"]}|{r["posterior_probabilities"]}|')
    md('PRIOR_AND_NUISANCE_SENSITIVITY.md','# Prior与nuisance敏感性\n\n各prior变动具有独立protocol；likelihood和材料条件prior不变，直接重新加权evidence，不新增response或物理。\n\n|案例|protocol|prior|posterior|\n|---|---|---|---|\n'+'\n'.join(prior_lines)+'\n\n相同gain观测的known-gain和g~N(1,.15²)结果为不同likelihood协议，保留logdet。完整值见results/PRIOR_AND_NUISANCE_SENSITIVITY.json。profile minima不是evidence；新raw/full-Mie补偿对照未运行。')
    md('METHOD_COMPARISON.md','# 同一density的积分方法对照\n\n|案例|方法|body wall秒|H0 Re SD|H0 Re q05|\n|---|---|---:|---:|---:|\n'+'\n'.join(methods)+'\n\n这些wall是实际本次summary过程，包含相应CDF与quantile计算，输入response coefficients已付费。不能作为full-solver cold/warm加速或新sampler优势；generic tempered/physics-guided full correction均未新增运行。')
    md('FAILURES.md','# 失败与未解决事项\n\n历史必要条件违反：\n\n'+failures+'\n\n原生四案例full-grid reference failure保留。历史bin-CDF与node covariance不同测度；原生/细化不是一个精确reference。四个旧treams oracle数组缺失。slogdet警告保留，Cholesky有限交叉核对通过。原生报告worker容量失败三次，父Codex完成，无模型/provider变化；成本未测。\n\nE1/非同心/大型NN/locked/新发布NOT_RUN。完整失败日志results/FAILURES.jsonl与controls warnings；任何未解决数值条件不得通过放宽门槛抹去。')
    md('GATE_DECISION.md','# 科研裁决\n\nCONTINUE_PROBABILITY_BENCHMARK_REPAIR。六案例的同一连续polynomial target摘要通过指定经验数值核对；不代表新Maxwell机制或一般成像成功。\n\n连续材料摘要修复：完成。历史coherence失败：保留。全域forward/CDF认证：未建立。新nuisance-aware secant、active sensing和全成本物理proposal收益：NOT_RUN。修复后population material coverage：NOT_RUN。\n\n预算历史384/384、新阶段0/0，P1建议4000未授权。维持当前概率benchmark，不启动重复campaign，不自动发布。后续机制若无新增解释力或成本收益，停止优越性主张；不能无限调模型。')
    md('README.md','''# 复现B1 E0

在本隔离worktree根目录安装原requirements。主入口：`python -m src.b1 e0`。已完成阶段禁止原地覆盖；独立复现需用新output目录并复制AUTHORIZATION.json（新增物理授权仍0）。不会导入/调用物理solver或写旧账本。

已付费观测成像：`python -m src.b1 infer --model continuous --observation results/cases/gain_marginalized/posterior.npz --output research/b1_next/e0/results/new_observation_output`。输出目录必须不存在。新refined/continuous入口使用同一连续density的moments与直接CDF；full-grid保留legacy节点诊断。NPZ必须有gain_sd，JSON可明确选择known/marginalized gain。

实际执行程序与原始数组见src/b1/continuous.py、e0.py、e0_checks.py及本目录results。报告程序为tools/e0_report.py，仅读取已完成结果。禁止在旧384目录再次运行audit/pilot/refine，不删除或降低账本。
''')
    render=dict(wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu,figures=8,new_physics=0,
                scope='measured report/plot body, excludes interpreter startup and initial imports; native service costs unmeasured')
    dump(out/'REPORT_RENDER_RECEIPT.json',render)
    with (out/'REPORT_RENDER_LEDGER.jsonl').open('a') as f:f.write(json.dumps(render)+'\n')
    resources['report_render']=render
    jobs=[json.loads(x) for x in (out/'JOB_LEDGER.jsonl').read_text().splitlines() if x.strip()]
    renders=[json.loads(x) for x in (out/'REPORT_RENDER_LEDGER.jsonl').read_text().splitlines() if x.strip()]
    resources['render_attempts']=renders
    resources['measured_job_wall_seconds_sum']=sum(r['wall_seconds'] for r in jobs)+sum(r['wall_seconds'] for r in renders)
    resources['measured_job_cpu_seconds_sum']=sum(r['cpu_seconds'] for r in jobs)+sum(r['cpu_seconds'] for r in renders)
    resources['measured_body_wall_total']=resources['wall_seconds_inference_and_controls']+resources['supplemental_CDF_and_public_replay']['wall_seconds']+sum(r['wall_seconds'] for r in renders)+resources.get('specified_slack_rediagnosis',{}).get('wall_seconds',0.)
    resources['scopes_not_double_counted']='Job phases are nested in main inference wall. Sum phases plus supplemental and render once; do not add whole inference wall again.'
    dump(out/'RESOURCE_SUMMARY.json',resources);summary['resources']=resources;dump(out/'SUMMARY.json',summary)
    md('RESOURCE_SUMMARY.md','# 完整资源说明\n\n新增物理0，历史384 BEGIN/384 END和high-water384保持。CPU-only、GPU/NN0、旧locked0、SHA0；未执行额外smoke或完整旧test suite。\n\n主body wall '+f'{resources["wall_seconds_inference_and_controls"]:.6f}'+'s，CPU '+f'{resources["cpu_seconds_inference_and_controls"]:.6f}'+'s；补充probes/CLI '+json.dumps(resources['supplemental_CDF_and_public_replay'],ensure_ascii=False)+'；报告绘图 '+json.dumps(render,ensure_ascii=False)+'。\n\npolynomial response evaluations与adaptive次数全部记在results/RESOURCE_SUMMARY.json，不算Maxwell调用。阶段JOB_LEDGER相加一次，不重复相加嵌套总wall。模型服务、人工交互、startup/import及旧历史测量按各自scope保留，不能把未测成本称0或称完整科研wall。')
    decision=dict(overall='CONTINUE_PROBABILITY_BENCHMARK_REPAIR',new_physics=0,historical_physics=resources['historical_physics'],
                  same_target_empirical_summary=summary['all_continuous_cases_converged'],specified_controls=summary['controls_pass'],
                  historical_summary_coherence_failures_preserved=stats['incompatible_row_count'],
                  uniform_certificate=False,full_Mie_new_explanatory_mechanism='NOT_RUN',material_coverage_recalculated=False,
                  E1='NOT_RUN_NO_AUTHORIZATION',P1_4000_authorized=False,public_release='NOT_RUN_NOT_AUTHORIZED')
    dump(stage/'GATE_DECISION.json',decision)
    print(json.dumps(dict(report=str(stage/'RESEARCH_REPORT_ZH.md'),figures=8,new_physics=0,render=render)),flush=True)

if __name__=='__main__':write(Path.cwd())
