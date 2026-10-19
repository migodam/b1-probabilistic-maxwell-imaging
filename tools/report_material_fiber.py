from pathlib import Path
import sys,json,html,platform,subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'tools'))
from material_fiber.core import OUT,write
from p1_material.prior import MaterialPrior
from analyze_material_fiber import endpoint
from publish_lxz import STYLE
from publish_research_f import rendered
REL='research/p2_material_fiber_20261019/'
WEB='https://migodam.github.io/b1-probabilistic-maxwell-imaging/'+REL
RELEASE='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p2-material-fiber-20261019'

def md(name,text):(OUT/name).write_text(text.strip()+'\n')
def figure(fig,name):
    for ext in ['png','pdf']:fig.savefig(OUT/f'figures/{name}.{ext}',dpi=165,bbox_inches='tight')
    plt.close(fig)
def table(header,rows):return '| '+' | '.join(header)+' |\n| '+' | '.join(['---']*len(header))+' |\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in rows)

def main():
    m=json.loads((OUT/'METRICS.json').read_text());s=json.loads((OUT/'SUPPLEMENT_RESULTS.json').read_text());r=json.loads((OUT/'RESOURCE_SUMMARY_FINAL.json').read_text());ac=json.loads((OUT/'ACCELERATION_FULL_COST.json').read_text());checks=json.loads((OUT/'SCIENTIFIC_NUMERICAL_CHECKS.json').read_text());parents=json.loads((OUT/'configs/PARENTS.json').read_text());anchors=json.loads((OUT/'A/ANCHORS.json').read_text());metrics=json.loads((OUT/'OFFLINE_PARENT_METRICS.json').read_text());cfg=json.loads((OUT/'configs/TARGET.json').read_text());P=MaterialPrior(cfg,json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']);selection=json.loads((OUT/'ONLINE_SELECTED_CANDIDATES.json').read_text());(OUT/'figures').mkdir(exist_ok=True)
    gallery=[]
    for parent in [p['parent'] for p in parents if p['split']=='heldout']:
        base=endpoint(anchors[parent]);truthkey=json.loads((OUT/f'offline_only/{parent}_truth.json').read_text())['key'];keys=[truthkey]+selection[parent][:5]
        fig,ax=plt.subplots(3,len(keys),figsize=(3*len(keys),8),squeeze=False,constrained_layout=True)
        for j,k in enumerate(keys):
            a=endpoint(k);chi=a['chi'];frac=P.decode(a['q'])[1]['fraction']
            images=[1+chi.real,chi.imag,frac]
            for i,(im,cmap,lo,hi) in enumerate(zip(images,['viridis','magma','cividis'],[1.1,.04,0],[2.9,.22,1])):
                obj=ax[i,j].imshow(im.reshape(24,24).T,origin='lower',extent=[-.5,.5,-.5,.5],cmap=cmap,vmin=lo,vmax=hi);ax[i,j].set_xticks([-.4,0,.4]);ax[i,j].set_yticks([-.4,0,.4]);ax[i,j].set_xlabel('x [m]')
                if j==len(keys)-1:fig.colorbar(obj,ax=ax[i,:],shrink=.75,label=['Re epsilon_r','Im epsilon_r','Inclusion fraction'][i])
                if j==0:ax[i,j].set_ylabel(['Re material','Loss','Phase fraction'][i]+' / y [m]')
            title='GT / OFFLINE only' if j==0 else ('Physics-only base' if j==1 else f'Candidate {j-1} / {k}')
            if j>=1:
                with np.load(OUT/f'inputs/online/{parent}.npz') as z:y=z['observation']
                phi=np.linalg.norm(a['prediction']-y)**2/cfg['noise']['sigma']**2
                title+=f'\nwhite residual={np.sqrt(2*phi):.3f}'
            ax[0,j].set_title(title,fontsize=10)
        valid=next(v['base_qualified'] for v in metrics if v['parent']==parent and v['method']=='ALL')
        fig.suptitle(parent+(' | compatible material alternatives; probabilities NOT estimated' if valid else ' | INVALID BASE: no compatible hypothesis found'),fontsize=13);figure(fig,'materials_'+parent);gallery.append(parent)
    # Strong optimizer history on the original exposed B.
    fig,ax=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
    for method in ['fixed_LM','adaptive_TR','material_TR','coarse_TR']:
        hh=[]
        for step in [8,16]:
            f=OUT/f'A/histories/P1_obs_B_{method}_checkpoint{step}.json'
            if f.exists():hh+=json.loads(f.read_text())
        ax[0].semilogy([v['white_residual'] if 'white_residual' in v else np.sqrt(2*v['Phi']) for v in hh],label=method)
        ax[1].semilogy([v.get('scaled_KKT',np.nan) for v in hh],label=method)
    ax[0].axhline(23.785449,color='k',ls='--');ax[0].set(title='Exposed old B: original physics and noise',xlabel='Measured linearization',ylabel='Whitened residual');ax[1].set(title='Feasibility does not certify stationarity',xlabel='Measured linearization',ylabel='Scaled projected-gradient diagnostic');ax[0].legend(fontsize=8);figure(fig,'optimizer_B')
    # Common qualified heldout panel for method and cost comparisons.
    br=json.loads((OUT/'B/SUMMARY.json').read_text())['blocks'];labels=[];cost=[];gain=[];div=[]
    for method in ['hard_prior','soft_frozen','soft_dynamic','soft_adaptive','full_J']:
        rr=[v for v in m['methods'] if v['method']=='B:'+method][0];blocks=[v for v in br if v.get('method')==method and v['parent'].startswith('P2_heldout') and 'RHS' in v];labels.append(method);cost.append(sum(v['RHS'] for v in blocks));gain.append((rr['mean_base_MAE']-rr['mean_BestOf4_MAE'])/rr['mean_base_MAE']*100);div.append(rr['mean_distinct_alternatives'])
    fig,ax=plt.subplots(1,3,figsize=(15,4),constrained_layout=True)
    for i in range(len(labels)):ax[0].scatter(cost[i],div[i],s=70);ax[0].annotate(labels[i],(cost[i],div[i]),xytext=(5,5),textcoords='offset points',fontsize=8)
    sup_blocks=json.loads((OUT/'B/BASELINE_COMPLETION.json').read_text())['results']
    for v in s['summary']:
        cc=sum(b.get('RHS',0) for b in sup_blocks if b.get('method')==v['method'] and b['parent'].startswith('P2_heldout'))
        ax[0].scatter(cc,v['heldout_distinct_alternatives']/4,s=90,facecolors='none',edgecolors='k',marker='D')
        ax[0].annotate(v['method'].replace('_sphere_only',' sphere (exposed)'),(cc,v['heldout_distinct_alternatives']/4),xytext=(5,5),textcoords='offset points',fontsize=7)
    ax[0].margins(x=.3,y=.2)
    ax[0].set(xlabel='Actual heldout exploration RHS (4 qualified parents)',ylabel='Distinct alternatives / parent',title='Hard split and sphere-only Full-J are strong controls')
    ax[1].bar(labels,gain,color=['#8e9b55','#538b9b','#2c6379','#394f83','#c6b69b']);ax[1].tick_params(axis='x',rotation=35);ax[1].set(ylabel='Best-of-4 MAE decrease [%]',title='Conditional on 4/8 qualified anchors')
    dl=[v for v in m['methods'] if v['method'].startswith('D:')];ax[2].bar([v['method'][2:] for v in dl],[v['mean_distinct_alternatives'] for v in dl]);ax[2].set(ylabel='Distinct alternatives / qualified parent',title='Fixed proposal-count pilot; actual RHS differ');figure(fig,'cost_accuracy_discovery')
    # Actual measured points; no invented continuous likelihood interpolation.
    parent=m['demo_parent'];fig,ax=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
    for method in ['hard_prior','soft_frozen','soft_dynamic','soft_adaptive']:
        bb=next((v for v in br if v.get('parent')==parent and v.get('method')==method and v.get('seed')==2026101910),None)
        if not bb:continue
        rr=[v for v in bb['rows'] if 'DeltaPhi' in v];ax[0].plot([0]+[v['material_RMS'] for v in rr],[0]+[v['DeltaPhi'] for v in rr],'o--',label=method);ax[1].plot([0]+[v['material_RMS'] for v in rr],'o--',label=method)
    ax[0].axhline(2,color='k',ls=':');ax[0].set(xlabel='Actual material RMS from base',ylabel='Measured DeltaPhi',title=parent+' / measured endpoints only');ax[1].set(xlabel='Saved predictor/corrector endpoint index',ylabel='Actual material RMS',title='No continuous-path certificate');ax[0].legend(fontsize=8);figure(fig,'measured_transport')
    fig,ax=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
    ax[0].bar([v['method'] for v in ac],[1000*v['amortized_seconds'] for v in ac]);ax[0].set(ylabel='Setup-amortized milliseconds / evaluation',title='16 fixed queries, three passes, equal cache rights');ax[0].tick_params(axis='x',rotation=25)
    ax[1].bar([v['method'] for v in ac],[v['fallbacks'] for v in ac]);ax[1].set(ylabel='Paid full fallback / 48 evaluations',title='Same-space identity is not an accuracy advantage');ax[1].tick_params(axis='x',rotation=25);figure(fig,'acceleration')
    demo=selection[parent][:5]
    with np.load(OUT/f'inputs/online/{parent}.npz') as z:y=z['observation']
    fig,ax=plt.subplots(1,2,figsize=(12,4),constrained_layout=True);ax[0].plot(y.real.T.ravel(),color='k',lw=1,label='Observed Re field')
    for k in demo:
        a=endpoint(k);ax[0].plot(a['prediction'].real.T.ravel(),lw=.7,alpha=.65,label=k);ax[1].hist(np.r_[(a['prediction']-y).real.ravel(),(a['prediction']-y).imag.ravel()]*np.sqrt(2)/cfg['noise']['sigma'],bins=24,histtype='step',label=k)
    ax[0].set(xlabel='Source-major receiver index',ylabel='Complex field real part',title='Exact Maxwell predictions');ax[1].set(xlabel='Whitened residual component',ylabel='Count',title='Fit diagnostics, not posterior probabilities');ax[0].legend(fontsize=7);ax[1].legend(fontsize=7);figure(fig,'forward_predictions')
    oldB=[]
    for method in ['fixed_LM','adaptive_TR','material_TR','coarse_TR']:
        f=OUT/f'A/runs/P1_obs_B_{method}_checkpoint16.json'
        if not f.exists():f=OUT/f'A/runs/P1_obs_B_{method}_checkpoint8.json'
        v=json.loads(f.read_text());oldB.append([method,f'{v["white_residual"]:.4f}',v['data_feasible'],v['stop_reason']])
    budgettable=table(['阶段','计费states / 上限','complex RHS / 上限','受控wall秒'],[[k,f'{r[k]["states"]}/{r[k]["caps"][0]}',f'{r[k]["RHS"]}/{r[k]["caps"][1]}',f'{r[k]["controlled_wall_seconds"]:.3f}'] for k in 'ABCDE'])
    methodtable=table(['方法','合格held-out数','不同替代/对象','Best-of-4 MAE','该panel基础MAE'],[[v['method'],v['qualified_anchor_parents'],f'{v["mean_distinct_alternatives"]:.2f}',f'{v["mean_BestOf4_MAE"]:.6f}',f'{v["mean_base_MAE"]:.6f}'] for v in m['methods'] if v['method'].startswith(('B:','D:'))])
    ac_table=table(['评价器','query+setup RHS','fallback /48','cold中位ms','warm中位ms','setup摊销ms'],[[v['method'],v['total_including_setup_RHS'],v['fallbacks'],f'{1e3*v["median_cold_seconds"]:.3f}',f'{1e3*v["median_warm_seconds"]:.3f}',f'{1e3*v["amortized_seconds"]:.3f}'] for v in ac])
    mainreport=f'''# Research P2：存在有限噪声材料自由度，但尚无独立的软谱生成算法优势

2026-10-19。最终裁决 **REDIRECT**。A–E已实际运行；没有将局部连续变化包装成多个重要posterior模式。所有后验概率为null。

## 先看这次究竟得到了什么

24个事前固定的新物体，15个得到合格physics-only基础材料；held-out仅4/8合格。首次冻结的材料archive在全部24对象保存69个任务尺度不同的替代解释，其中held-out为29个；这是按每份观测累加、最多12张代表的材料去重计数，不是69个拓扑模式或概率事件。

主held-out archive有1个达到原结构变化尺度的替代材料，0个拓扑事件改变。训练/验证的少数拓扑标签变化也不稳定：同时满足结构检查并在.25/.5/.75三个fraction阈值都改变的案例为0。故“可探索连续材料不确定性”有正证据，“发现重要不同拓扑世界”尚未建立。

全部8个held-out的基础复材料MAE为 **{m['heldout_mean_base_MAE']:.8f}**，在线farthest-first代表的Best-of-4为 **{m['heldout_mean_BestOf4_MAE']:.8f}**，下降约 **{100*(1-m['heldout_mean_BestOf4_MAE']/m['heldout_mean_base_MAE']):.2f}%**。2,000次parent bootstrap的绝对改善区间为[{m['heldout_BestOf4_improvement_CI'][0]:.8f}, {m['heldout_BestOf4_improvement_CI'][1]:.8f}]。它混合了多种方法及不同预算，只是联合archive的oracle诊断，不是某个算法的同成本胜出。

Best-of-K包含原基础图，oracle最小值按定义不会恶化；其bootstrap区间不能作为在线估计改进或新算法独立贡献的证据。在线默认Selected-MAE仍为基础图的.06513030。

动态软谱在合格子集的Best-of-4稍好，但比hard-prior使用更多RHS；动态减hard的配对MAE区间跨0。小Gaussian/Flow没有改善Best-of-4。A17更慢且fallback更多。因此本轮不支持“Soft Explorer / learned proposal / A17-X产生强算法增量”的论文主张。

![真实材料图库](figures/materials_{parent}.png)

## 1. 物理与先验，没有改成更容易的任务

原P1的64维相关Gaussian加两项uniform，C2被动材料fraction、24×24 scalar TM、8 sources/32 complex receivers、原host、gain=1及sigma=.00119211798514519全部保留。宿主区域不是空气孔。观测使用真正circular-complex Gaussian noise，实虚部方差sigma²/2；没有逐物体调SNR或noise。

24个新parent为12 train/4 validation/8 held-out，seed2026101901；旧A/B只是已暴露回归条件。图像显示epsilon_r=1+chi，实部、损耗与fraction全使用共同色标。这里不是三维vector-Maxwell或一般未知几何成像。

生成truth的24次forward计费。在线读取sanitize后的观测、合法材料及已付费物理因子；GT只在commit2024d77固定A–E之后，由独立offline脚本读取。后续sphere-only强基线补充发生在主评分之后，明确为EXPOSED审查，不能称新增盲测。

## 2. 先把物理反演修好

旧B的白化残差从309.787降到以下水平。固定.5截断LM继续运行仍不合格；直接TR、材料TR和coarse-to-fine都得到相容材料，但部分QP停止为SUBPROBLEM_FAILURE或预算停止，不能追认stationary MAP。

{table(['优化器','旧B最终白化残差','noise-ball资格','停止原因'],oldB)}

新物体上按开发条件选择coarse_TR，最终15/24相容，9/24未恢复。这是优化与decoder数值条件问题的实证，不是“其余物体不存在合法材料”。全部真值本来由原prior生成；noise-ball只是有限相容判据，也不是Gaussian posterior硬支持。

![物理收敛](figures/optimizer_B.png)

## 3. 局部探索：真正的收益和代价

以M=[Re Dg;Im Dg]/24构造真实材料RMS frame，prior-aware矩阵同时包含数据项和原Gaussian penalty。soft cutoff由任务距离s决定：lambda_c=1/s²；没有根据held-out GT挑谱阈值。筛选只影响提案，未改变prior measure。

全部方法保持合法decoder、实际RMS .02/.04、原noise与property bounds。主探测器的signed progress阻止corrector回原点。每个接收端点都立即保存；.08只在development做stress。dynamic/adaptive的正结果是有限端点，不是精确fiber存在证明。

{methodtable}

上表方法只在已获得合格基础解的4个held-out对象执行。其余4个计入全部8对象的主风险，不能只拿这4个容易案例宣布一般恢复成功。K不足时保持实际K，未复制材料凑数。默认Top-1基础图始终保留；oracle Best-of-K不用于在线筛选。

![成本与精度](figures/cost_accuracy_discovery.png)

普通随机、native及带随机方向进度约束的Full-J主对照较弱；这不证明Full-J本身不能探索。独立审查因此在B剩余额度内增加sphere-only Full-J与hard-prior初始化：只保持非零材料距离，取消人为角向限制。

补充对照在held-out找到 **10个Full-J、9个hard-prior** 去重替代解释；Full-J有1个结构变化，仍无重要拓扑模式证据。其全部8对象Best-of-4分别为{s['summary'][0]['all8_mean_BestOf4_MAE']:.8f}和{s['summary'][1]['all8_mean_BestOf4_MAE']:.8f}。这是已暴露审查，不作为新确认性胜出；但足以否定“只有软谱机制才能得到替代材料”的叙事。

## 4. 弯曲材料路径，有什么证据？

实际验证了RMS .02/.04的predictor/corrector端点，有方向更新与非线性数据修正。点之间未做稠密full-forward验证，无法承诺整个插值路径都相容，更不能称为精确等散射流形。图中的虚线仅连接已测点。

![已测有限端点](figures/measured_transport.png)

GT弱方向能量与Best-of-4收益在4个合格held-out上的Spearman约{m['weak_GT_energy_vs_MAE_gain_spearman']:.2f}，样本太少，且没有正面可达性证据。大量GT差量投到弱空间，不足以推出Explorer能走到GT。

## 5. 全局搜索与正确概率核

在4 train+4 held-out上，独立multistart、prior jump+共同corrector与tempered pCN都真实运行。全局方法没有稳定恢复额外结构分支；高温与prior链的大量计算没有转化为冷链任务尺度移动。短链mixing不足，故不输出posterior质量、credible interval、ESS校准或无漏模式保证。

pCN与反射property proposal在原prior下可逆，温度t目标为p0(q) exp(-t Phi)。MH接受使用完整likelihood，交换使用正确温度比。deterministic corrector输出另记candidate-only。算法形式正确不表示这次16-sweep链已经收敛。

## 6. 学习提案没有增量

8个train parents提供54个真实合格、任务尺度、likelihood-competitive更新；validation有16个。三个seed的小Gaussian与四层Flow均训练3,000 updates，参数约3.29万/10.72万。没有GT监督，没有posterior distillation，没有训练期Maxwell调用，GPU为0。

学习的只是64维Gaussian更新分布，两项uniform property使用共同反射提案。held-out非学习对照平均每个合格parent有5.5个去重替代，Gaussian和Flow各.5个；二者Best-of-4均未优于基础图。D是相同提案次数/物理预算上限的pilot，decoder失败使实际RHS不同；不能把该计数当严格相同实际RHS的统计排名。失败生成的CPU成本保留在D阶段总wall，早退项没有单独可靠计时。

不能据此否定全部learned Explorer，但当前数据量、表示和小模型pilot没有值得扩大训练的实证。

## 7. A17-X：恒等式正确，实际增量否定

固定16个真实材料提案、相同噪声及缓存权利，包含所有端口setup、small core、pool、fallback和冷warm成本。exact T-matrix使用真实256-cell ROI，不把66 latent当66 ports。相同[V,C]下Schur exchange和直接Galerkin的最大差约1.91e-15。

{ac_table}

exact_T表中的1056 RHS包含四个anchor及256列传播；不是query-only的0 RHS。A17的17/48 fallback高于goal-ROM的4/48，setup摊销时间约为optimized full的11.1倍。ordinary goal-ROM也未胜过full。warm basis随历史增长，rank与所有矩阵动作均保留；这些是本小型dense平台的测量，不外推到大型三维。

![评价器实际成本](figures/acceleration.png)

## 8. 数学贡献与先验工作的边界

材料尺度归一化与soft matrix filter可避免单看latent弱奇异值；线性重参数化并正确变换prior Hessian后，材料空间滤波相对差为5.27e-14。这是实现与几何一致性，不是新首创定理。LIS/DILI早已利用prior/likelihood谱结构；Riemannian MCMC已有位置依赖几何。

本轮可检验的新机制候选是：在材料任务尺度下分离decoder secant与Maxwell余项，决定何时付费刷新并保持有限进度。然而本次soft/dynamic并未在完整成本上明确优于hard split与强Full-J，因此尚不能将它升级为强IEEE TAP贡献。

## 9. 全部成本、失败与下一步

{budgettable}

最终 **{r['total']['states']}计费states / 4000、{r['total']['RHS']} complex RHS /132000、{r['total']['controlled_wall_seconds']:.3f}秒受控执行wall**。其中{r['total']['full_native_forward_calls']}次真实fresh native full-forward/LU；E另保守计费固定材料ROM与factor-reuse尝试。物理operation spans合计{r['total']['physics_operation_wall_seconds']:.3f}秒；它不含完整人机研究、编码、发布时间。NN训练CPU、GPU0。所有余量封账，不借旧余额。

C曾因property=-6.94e-18触发严格检查。原8states/320RHS及失败保持；只读回放确认机器精度误差后，以显式latent边界恢复继续，同一总账不重置。基础反演的subproblem failure、局部distance失败、重复点、短MH混合失败和无收益NN均保存。

**下一步选择REDIRECT：**保留便宜hard-prior/普通Full-J作为局部有限噪声基准，停止当前Flow与A17-X的优势主张。优先建立能被强Full-J与全局基线共同核验、具有稳定结构区别的ambiguity benchmark，或改进measurement设计；没有这样的正面实验，不训练更大的Explorer，也不包装强TAP论文。

[完整数据与代码]({RELEASE}) · [逐项裁决](GO_NO_GO_DECISION.md) · [理论](THEORY_REFINEMENTS.md) · [复现](REPRODUCE.md)
'''
    md('RESEARCH_REPORT_ZH.md',mainreport)
    md('START_HERE.md',f'''# Research P2：真实材料探索已完成，整体REDIRECT

2026-10-19。原P1、24个fresh parents、A–E与补充强基线均实际执行。不是新的概率posterior论文。

[中文报告]({WEB}report.html) · [材料图库]({WEB}gallery.html) · [证据Release]({RELEASE}) · [GitHub目录](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/{REL})

- 15/24合格物理解，held-out 4/8；旧B由309.787恢复至23.69左右。
- 预冻结archive在held-out有29个材料去重替代、1个结构变化、0个拓扑改变；概率为null。
- 全8 held-out联合Best-of-4 MAE .06513030→.06424926，仅1.35%；不是单算法同成本胜出。
- 便宜hard-prior是强对照；去掉角向限制后普通Full-J同样能发现替代材料。
- 当前Gaussian/Flow没有增量；A17摊销约为full的11.1倍，停止集成优势主张。
- 最终{r['total']['states']}/4000 states、{r['total']['RHS']}/132000 RHS、受控wall{r['total']['controlled_wall_seconds']:.3f}s，GPU0。旧账本未改。

阅读顺序：RESEARCH_REPORT_ZH.md → GO_NO_GO_DECISION.md → METHOD_COMPARISON.md → THEORY_REFINEMENTS.md → RESOURCE_SUMMARY_FINAL.json。

METRICS.json是补充对照前的预冻结主分析；SUPPLEMENT_RESULTS.json是评分后EXPOSED强基线审查，二者不混为新的盲测。MaterialHypothesisPackage/FINAL为最终交接包；不相容base留在failed_baseline，绝不作为合格假设。原始payload和所有失败保留。
''')
    md('METHOD_COMPARISON.md','# 同物理目标、完整费用比较\n\n'+methodtable+'\n\n主局部表为4个合格held-out，不能隐去另外4个基础恢复失败。soft_dynamic相对hard的Best-of-4 MAE差−.000333，95%配对区间约[−.002487,.002063]，未确立优势。hard 136 RHS，dynamic504 RHS，adaptive440 RHS（相同已执行held-out blocks）。\n\n## 补充强基线\n\nFull-J sphere-only发现held-out10个替代，hard-prior sphere-only9个。该审查发生在主评分后，明确为exposed；原角向progress约束可能限制Full-J，不能据主对照零成功证明其能力不存在。\n\n## A17与标准评价器\n\n'+ac_table+'\n\n所有setup纳入，full允许自有factor重用，普通ROM允许增量basis，固定材料不被corrector改变。')
    md('GO_NO_GO_DECISION.md','''# 最终研究裁决：REDIRECT

- GO：原P1确实存在任务尺度的有限噪声材料替代；可保留为局部不确定性benchmark。
- NOT_ESTABLISHED：精确fiber、连续全路径相容、重要拓扑posterior modes、概率校准。
- NO-GO当前优势主张：动态软谱相对便宜hard-prior/强Full-J的完整成本增量；小型Gaussian/Flow；本TM24平台A17-X集成。
- REDIRECT：下一步先获得稳定结构差异、可被强基线核验的ambiguity条件或主动测量问题。
- NO-GO当前强TAP算法论文：缺多独立对象上的稳定结构机制、同成本优势及可信事件概率参考。

不再追加本campaign预算。不以噪声放宽、先验扩张、GT初始化或大网络救结果。旧A17局部GN论文结论不受本平台否定结果替代。
''')
    md('RESOURCE_SUMMARY.md','# 最终资源封账\n\n'+budgettable+f'\n\n总计{r["total"]["states"]}计费states、{r["total"]["RHS"]} complex RHS、{r["total"]["controlled_wall_seconds"]:.3f}s受控wall、{r["total"]["physics_operation_wall_seconds"]:.3f}s物理operation wall，GPU0。\n\n2012次fresh full-state/LU与224次保守计费的ROM/factor-reuse尝试分开。small-core RHS、矩阵动作、pool构造、训练CPU和缓存读写在RESOURCE_LEDGER逐项记录。历史26states/976RHS与其他研究均未改。端口setup1056 RHS不能从exact-T query-only列中隐藏。\n\nRESOURCE_SUMMARY.json保留主分析时的2152states/44552RHS；RESOURCE_SUMMARY_FINAL.json包含后续84states/2176RHS强基线审查，不是重复授权。')
    md('ITERATION_LOG.md','''# 实际研究迭代

1. 假设：固定Gaussian .5截断阻碍恢复。四优化器按train/exposedB比较；coarse_TR取得4/5相容并被冻结，24fresh最终15/24相容。未通过者保留，不用GT替换。
2. 假设：材料度量与prior-aware谱有助有限探索。八臂真实.02/.04对照；dynamic/adaptive成功率较高。单一修订为刷新触发分离decoder与Maxwell secant；paired v2保留，held-out随后冻结执行。
3. 假设：非局部prior jump或tempering补充分支。global方向未稳定发现重要结构，短cold chain移动不足；继续保留负结果。
4. C数值错误：property=-6.94e-18。暂停、只读回放、保存原8state/320RHS；仅机器精度latent边界恢复。同stage ledger累计，未重置，失败parent的首次启动计算仍计费。
5. 学习假设：54个真实train更新可摊销方向提案。Gaussian/Flow三seed实际训练并物理验证，未优于非学习对照；不扩大网络。
6. 求值假设：A17有限exchange降低费用。同空间误差1.91e-15，但A17更慢/fallback更多，停止增量主张。
7. 独立审查：随机angular progress可能削弱普通Full-J。在B剩余额度补充sphere-only两强基线；新增84states/2176RHS，明确exposed。Full-J可发现替代，削弱软谱独占机制叙事。原主分析未改。

最后不是再调参数：整体REDIRECT，停止当前学习与A17优势路线。阶段余量封账。
''')
    md('FAILURES.md','''# 完整失败及解释

9/24基础解未进入noise-ball，包含4/8 held-out。其真实材料由同P1生成；不支持先验无解结论。部分合格端点也在subproblem失败处停下，不称MAP驻点。

材料distance/progress、Gaussian局部半径或bounds不满足时不付费forward，普通数值时间记录；实际forward失败/拒绝、line search均计入。B blocks保留每次decoder失败与BLOCK_PHYSICS_LIMIT。

C原始ERROR、CLOSED、BOUNDARY_ERROR_DIAGNOSIS及恢复receipt都保留。机器精度恢复阈值64 eps；不是任意物理材料clip。

高温采样没有给出可信posterior，ESS或低残差不证明mode覆盖。Flow/Gaussian训练目标是少量corrector轨迹，不能称完整条件posterior或进行概率校准。

主Full-J对照有angular-progress confound，补充sphere-only揭示其不能被当弱陪跑。补充已exposed，不能回溯升级主blind结果。

11个原archive拓扑标签变化中，只有1个同时通过结构检查；没有同时在三个阈值稳定的结构拓扑改变。主held-out无拓扑改变。不能以阈值闪烁制造多解。

Semantic Scholar focused collect返回429，保留访问缺口；后续只核对少量原始arXiv/出版方记录，不称穷尽新颖性搜索。原生subagent因已有线程槽位满未启动；核心全部由主Codex执行，无DeepSeek。
''')
    md('THEORY_REFINEMENTS.md',r'''# 理论：哪些成立，哪些不是新定理

## 材料尺度与prior测度

令M=[Re Dg;Im Dg]/24，M=Q S V^T，T=V S^{-1}，则MT=Q。K=J_white T。加入T^T H_prior T描述原Gaussian局部代价；它不是新的材料volume density。Dg的数值核、uniform property边界和可行切锥仍须处理。

可逆线性坐标变换q=C a时，M、J、H分别推送为MC、JC、C^T H C；对相同物理材料速度，Q W Q^T不依赖该参数标签。实测一处缓存相对差5.27e-14。非线性重参数化须另处理坐标加速度，这里没有宣称一般hessian tensor不变。

## soft滤波的连续性

对PSD H，W(H)=(I+s²H)^{-1}。精确恒等式W(B)-W(A)=-s² W(B)(B-A)W(A)，故norm差≤s² norm(B-A)。它解释soft filter避免hard截断在谱交叉处跳变；这是标准resolvent分析，不是Maxwell新定理，也不说明一定更便宜。

## 有限数据竞争的条件

实际材料增量h=g(q+d)-g(q)，完整白化变化为J_chi h+e。令a=J_chi h，则

ΔPhi = r^T a + ||a||²/2 + (r+a)^T e + ||e||²/2。

若有严格||e||≤epsilon，线性预测的Phi误差不超过||r+a||epsilon+epsilon²/2。这里只测secant与full端点，没有uniform bound，因此刷新是经验规则，不能标certified early acceptance或路径证书。decoder余项h-Dg d应先与Maxwell余项分离。

## 先验边界与精确fiber

原P1局部数值满秩不提供非零exact-fiber tangent；finite-noise不确定性无需exact null。合法端点不推出整个路径合法/相容，GT差量弱能量不推出可达性，corrector低residual不推出新结构或posterior质量。

## 与已有研究的关系

[Cui et al., 2014, LIS](https://arxiv.org/abs/1403.4680)已用prior与likelihood相对影响构造低维推断。这里增加材料任务metric，不把谱分解本身写作原创。

[Cui, Law, Marzouk, 2016, DILI](https://arxiv.org/abs/1411.3688)已研究operator-weighted、likelihood-informed MCMC；本轮candidate corrector不是DILI或合法posterior transport的替代。

[Girolami and Calderhead, manifold MCMC](https://doi.org/10.1111/j.1467-9868.2010.00765.x)已有位置几何采样。本轮没有证明新的不变测度结果。

A17 Schur/current-space恒等式保留。相同空间精度相同；新增价值只能来自选择、复用或误差控制。本轮该费用收益为负。
''')
    md('POSTERIOR_VALIDATION.md','''# 概率状态：没有P1可信完整reference

完整目标为p0(z64,u2) exp(-Phi)，Gaussian×uniform原prior，gain固定，sigma原值。没有变换成材料Lebesgue均匀prior。

C温度链的pCN+反射property proposal为prior可逆，使用完整likelihood MH；交换目标也正确。代数检查最大prior抵消误差2.49e-14。只有16 sweeps，mixing不足；不报告posterior mass、KL、CDF credible coverage或calibration。

Flow/Gaussian raw delta-z log probability通过显式Gaussian/affine-coupling密度定义；这只是提案密度。material retraction、固定距离和deterministic corrector后的频率不是posterior。概率全部null，不用residual softmax。

Package剩余已建模mass、未探索mass与模型外mass均无法量化，不把展示K项重新归一化到100%。已有Joint P2 Mie reference属于另一target，未用于这里的P1概率背书。
''')
    md('LOCAL_FIBER_RESULTS.md','''# 局部有限噪声探索

实测.02/.04端点存在；它们是合法P1材料、精确Maxwell相容并具有竞争likelihood的连续变化。主held-out29个去重替代、1个结构变化、0拓扑改变。数值不同不等于分离mode。

高阶feedback不是必需主模块。便宜hard-prior找到替代，dynamic比它多付RHS但MAE差区间跨0。完整插值路径未验证；不声称GT通过局部fiber可达。
''')
    md('GLOBAL_MODE_RESULTS.md','''# 跨结构与全局探索

两种nonlocal corrected search和tempered pCN在8个固定对象实际执行，所有起点、温度、拒绝、失败计费。全局路线没有产生稳定held-out拓扑分支；少数train竞争端点不构成重要模式reference。

没有在新prior、放大噪声或自由gain下寻找正例。当前P1仍是66维C2相关场；转新decoder仅为下一阶段独立协议，不混入本轮统计。
''')
    md('NEXT_RESEARCH_DIRECTIONS.md','''# 下一步：不是训练更大的Explorer

选择REDIRECT。保留hard-prior和sphere-only Full-J为低成本基线，继续把原P1作为连续有限噪声不确定性benchmark。

下一决定性实验应在新授权和独立协议下建立一对稳定结构不同、同prior且具有可核验支持的Maxwell ambiguity条件；或者研究新增测量如何排除当前不同材料解释。先让强Full-J/global方法也能复现该机制，再比较Explorer。

新prior若必要，必须定义实际物理相关尺度、被动材料、生成测度及acquisition，先证明支持与参考质量，不根据失败改变P1。原P1 exact-fiber、当前小Flow与A17加速路线不继续追加预算。

目前缺少强TAP论文的关键正例：多独立对象的稳定结构发现、同完整成本优于强基线、可解释Maxwell机制与必要的概率参考。Schur、LIS、普通多起点、Flow均不能单独充当创新。
''')
    md('REPRODUCE.md','''# 复现与审计

本轮原始账本已封闭，运行器默认拒绝重跑已完成阶段。下载证据Release到现有B1根目录；不要删CLOSED或ledger来重置本轮额度。

依赖：Python3.13、NumPy2.5.3、SciPy1.18.1、Torch2.14.0、Matplotlib3.11.1。物理使用A30 commit717d3d26a290248bd2b9aeb55f3d03694b407b84；可用B1_P2_A30_ROOT指定已存在checkout。背景数组和selected-anchor caches在Release。当前源码入口tools/run_material_fiber.py，阶段A/B/C/D/E由统一Book计费。

tools/analyze_material_fiber.py为隔离GT评价；它在封账后拒绝覆盖既有主分析。报告重建入口tools/report_material_fiber.py只读已有数组，不产生Maxwell。tools/check_material_fiber_evidence.py检查已保存物理方程和测度代数，不支付新RHS。完整从零物理复算需要独立的新执行目录、明确新预算和保留原账本，而不是把本证据包当免费solver bank。

主要代码：src/material_fiber/core.py、optimize.py、explore.py、global_search.py、learn.py、acceleration.py、baseline_audit.py。所有随机种子与已暴露补充协议在configs。模型checkpoint小于1MB/个；生成过程候选/失败日志与raw endpoints保留。

Release保留所有raw材料/current/field/prediction、训练weights、原prior与背景、选中anchor的LU/receiver transfer；未交接的其余LU属于可重建计算缓存，仍在本地保存，未作为免费部署成本。任何重新分解、RHS或额外频率都应按新预算计费。
''')
    write(OUT/'DEPENDENCIES.json',dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=__import__('scipy').__version__,matplotlib=matplotlib.__version__,torch=__import__('torch').__version__,A30_commit='717d3d26a290248bd2b9aeb55f3d03694b407b84',C1_reference_commit='83056881e89329212b16f5be47a4d0dbc4465943',new_GitHub_repo=False))
    nav=f'<nav><a href="../../">B1首页</a><a href="index.html">本轮入口</a><a href="report.html">中文研究报告</a><a href="gallery.html">材料图库</a><a href="START_HERE.md">GPT阅读入口</a><a href="{RELEASE}">原始证据</a></nav>'
    top='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>P2 · Maxwell材料探索</title><style>'+STYLE+'</style><main>'+nav
    body=f'''<h1>材料确实能动，但“多拓扑世界”尚未建立</h1><p>Research P2 · 2026-10-19 · A–E已真实运行 · 最终REDIRECT</p><section><h2>这次得到的真实结果</h2><div class="metrics"><span><b>15 / 24</b>合格物理基础解</span><span><b>29</b>held-out材料去重替代</span><span><b>1.35%</b>联合Best-of-4 MAE下降</span><span><b>0</b>held-out拓扑改变</span></div><p>29个替代来自4个合格held-out对象，1个通过结构变化检查；它们不是29个posterior modes。另4个基础恢复失败仍计入全部8对象统计。概率全部未估计。</p><p><strong>更强Full-J同样找到替代，削弱了软谱独占机制的主张；小Flow和A17没有独立收益。</strong></p><a href="report.html">完整结论、失败和对照 →</a></section><section><h2>同一观测下的实际材料解释</h2><img src="figures/materials_{parent}.png"><p>GT只离线评分；其余材料由在线探索生成并经完整Maxwell验证。相同色标，2D scalar TM。不是多张图片的概率生成。</p><a href="gallery.html">全部8个held-out及失败图 →</a></section><section><h2>强基线决定研究方向</h2><img src="figures/cost_accuracy_discovery.png"><p>hard-prior更便宜；dynamic相对hard的配对MAE区间跨0。sphere-only Full-J补充审查标记为exposed，原主分析保持。</p></section><section><h2>完整费用与明确停止</h2><p>{r['total']['states']} /4000计费states；{r['total']['RHS']} /132000 complex RHS；受控执行{r['total']['controlled_wall_seconds']:.3f}秒；GPU0。编码、研究与发布时间不包含在物理执行wall里。</p><p>停止当前Flow和A17优势主张；保留低成本材料不确定性基准，下一步需要稳定结构ambiguity或主动测量正例。</p></section>'''
    (OUT/'index.html').write_text(top+body+'</main></html>');(OUT/'report.html').write_text(top+'<article>'+rendered(mainreport)+'</article></main></html>')
    gal=top+'<h1>全部held-out：真值、基础解与候选</h1><p>相同材料色标；GT只用于离线评分。不合格base明确标INVALID，不展示为合格候选。候选不等于posterior模式。</p>'
    for par in gallery:gal+=f'<section><h2>{par}</h2><img src="figures/materials_{par}.png"><p><a href="MaterialHypothesisPackage/FINAL/{par}.json">机器可读候选及完整物理prediction引用</a></p></section>'
    for name in ['optimizer_B','cost_accuracy_discovery','measured_transport','forward_predictions','acceleration']:gal+=f'<section><h2>{name}</h2><img src="figures/{name}.png"><a href="figures/{name}.pdf">PDF</a></section>'
    (OUT/'gallery.html').write_text(gal+'</main></html>')
    notice=f'## 最新：P2材料探索完整Campaign\n\n[中文报告]({REL}RESEARCH_REPORT_ZH.md) · [网页与材料图]({WEB})\n\n24个新物体，15个合格anchor；有限噪声替代存在，强基线削弱软谱优势。整体REDIRECT，概率不估计；{r["total"]["states"]}states/{r["total"]["RHS"]}RHS封账。\n\n'
    for name in ['README.md','START_HERE.md']:
        path=ROOT/name;text=path.read_text()
        if REL not in text:path.write_text(notice+text)
    path=ROOT/'index.html';text=path.read_text()
    if REL not in text:path.write_text(text.replace('<div class="grid">','<div class="grid"><section class="card"><h2>最新：P2材料探索</h2><p>24物体完整campaign；有限噪声替代存在，软谱、学习及A17优势未建立。</p><a href="'+REL+'">报告与实际材料图 →</a></section>',1))
    path=ROOT/'research/index.html';text=path.read_text()
    if 'p2_material_fiber_20261019/' not in text:path.write_text(text.replace('<main>','<main><section><h2>最新：P2材料探索</h2><p>整体REDIRECT；完整强基线、负结果与资源账本。</p><a href="p2_material_fiber_20261019/">报告与图库 →</a></section>',1))
    print(json.dumps(dict(web_prepared=True,physics_calls=0,figures=len(list((OUT/'figures').glob('*.png'))),overall='REDIRECT')))

if __name__=='__main__':main()
