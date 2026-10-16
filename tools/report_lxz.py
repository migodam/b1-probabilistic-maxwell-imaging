"""Regenerate LXZ reports from saved runs, never performs a physical solve."""
from pathlib import Path
import json,sys,html,csv,time
import numpy as np
from collections import Counter
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from lxz.proxy import OUT,dump,plain
from b1_joint.material import geometry,event_ids,decode,RMAX

def text(name,body): (OUT/name).write_text(body.strip()+'\n',encoding='utf-8')

def main():
    start=time.perf_counter();results=OUT/'results'
    r1=json.loads((results/'round1/METRICS.json').read_text())['rows']
    g1=json.loads((results/'round1_geometry_repair/METRICS.json').read_text())['rows']
    g2=json.loads((results/'round2_geometry_repair/METRICS.json').read_text())['rows']
    rows=[dict(x,method_label=x['method'],source='round1') for x in r1 if x['method']!='geometry']
    rows += [dict(x,method_label='geometry',source='round1_geometry_repair') for x in g1]
    rows += [dict(x,method_label='geometry_nuisance',source='round2_geometry_repair') for x in g2]
    methods=['local','geometry','geometry_nuisance','global','hybrid','tempered','prior'];summary=[]
    for n in [1024,4096]:
        for m in methods:
            a=[x for x in rows if x['method_label']==m and x['checkpoint']==n]
            summary.append(dict(method=m,checkpoint=n,runs=len(a),event_MAE=np.mean([x['event_MAE'] for x in a]),
                CDF_grid_error_mean=np.mean([max(x['CDF_max_error']) for x in a]),recall_at_point005=np.mean([x['significant_event_recall'] for x in a]),
                full4096_arm_seconds_mean=np.mean([x['inclusive_seconds'] for x in a]),
                acceptance=np.mean([x.get('acceptance_checkpoint',0) for x in a]) if m!='tempered' else None))
    boot=[];rng=np.random.default_rng(26101630)
    for alternative in ['local','hybrid','geometry_nuisance']:
        dif=[]
        for case in ['solid','large_cavity','annular_cavity','combined']:
            a={r['seed']:r for r in rows if r['checkpoint']==4096 and r['observation_id']==case and r['method_label']=='global'}
            b={r['seed']:r for r in rows if r['checkpoint']==4096 and r['observation_id']==case and r['method_label']==alternative}
            dif.append([a[s]['event_MAE']-b[s]['event_MAE'] for s in a])
        dif=np.array(dif);rep=[]
        for _ in range(2000):
            ids=rng.integers(0,4,4);rep.append(np.mean([np.mean(dif[i,rng.integers(0,4,4)]) for i in ids]))
        boot.append(dict(comparison='global minus '+alternative,mean=dif.mean(),interval95=np.quantile(rep,[.025,.975]),
            scope='descriptive nested bootstrap across four exposed observations and four algorithm seeds; no new population validation'))
    pilot=json.loads((ROOT/'research/joint_p2/results/PILOT_REFERENCE.json').read_text())['cases']
    functionals=[]
    for r in rows:
        if r['checkpoint']!=4096:continue
        case=r['observation_id'];m=r['method'];p=results/r['source']/f'{case}_{m}_{r["seed"]}.npz';z=np.load(p)
        if m=='tempered': x=z['4096_latent'];w=z['4096_weights']
        else:x=z['latent'][len(z['latent'])//4:];w=np.full(len(x),1/len(x))
        v=np.c_[geometry(x),x[:,2:]];eid=event_ids(x);conditional=[]
        for ref in pilot[case]['C']['continuous']:
            ev=ref['event_id'];mask=np.ones(len(x),bool) if ev is None else eid==ev
            if not mask.any():conditional.append(dict(event=ev,status='NO_SAMPLES'));continue
            a=v[mask];ww=w[mask];ww=ww/ww.sum();mean=ww@a;var=ww@((a-mean)**2);qs=[];cdferr=[]
            for j in range(4):
                ids=np.argsort(a[:,j]);cdf=np.cumsum(ww[ids]);cdf[-1]=1
                qs.append(a[ids[np.searchsorted(cdf,[.05,.5,.95])],j])
                cs=np.array([ww[a[:,j]<=t].sum() for t in ref['CDF_grid'][j]])
                cdferr.append(float(np.max(abs(cs-ref['CDF'][j]))))
            qs=np.array(qs);lower=.05*(np.maximum(mean-qs[:,0],0)**2+np.maximum(qs[:,2]-mean,0)**2)
            conditional.append(dict(event=ev,sample_count=int(mask.sum()),mean=mean,SD=np.sqrt(var),quantiles=qs,
                CDF_grid_error=cdferr,coherence_pass=bool(np.all(var+1e-12*np.array([*RMAX,.6,.12])**2>=lower)),
                summary_measure='same empirical weighted chain/particle measure; reference continuous density kept separate'))
        functionals.append(dict(case=case,method=r['method_label'],seed=r['seed'],conditional=conditional))
    dump(OUT/'POSTERIOR_FUNCTIONALS.json',functionals)
    e2=json.loads((results/'E2_completion/CANDIDATES.json').read_text());local=[];packages=[]
    for g in e2['geometry']:
        case=g['anchor'];raw=e2['rows'];z=np.load(results/f'E2_completion/{case}.npz');unique=[];archive=[]
        for m in ['local','geometry','prior','continuation']:
            rr=[r for r in raw if r['anchor']==case and r['method']==m]
            good=[r for r in rr if r['compatible']];actual=[]
            for r in good:
                material=z[f'chi_{r["id"]}']
                if any(np.max(abs(material-a))<1e-10 for a in actual):continue
                actual.append(material)
                if not any(np.max(abs(material-a))<1e-10 for a in unique):
                    unique.append(material);archive.append(dict(record=r,array_key=f'chi_{r["id"]}',posterior_mass=None))
            local.append(dict(anchor=case,method=m,proposals=len(rr),compatible_records=len(good),deduplicated_compatible=len(actual),
                MH_accepted=sum(r.get('MH_accepted',False) for r in rr),structure_qualified=sum(r['compatible'] and r['structure_qualified'] for r in rr),
                maximum_compatible_material_MAE_to_anchor=max([r['material_MAE_to_anchor'] for r in good],default=0),
                maximum_compatible_material_RMS_to_anchor=max([r['material_RMS_to_anchor'] for r in good],default=0),
                compatible_and_log10_competitive=sum(r['compatible'] and abs(r['delta_ell'])<=np.log(10) for r in rr)))
        packages.append(dict(observation_id='PUBLIC_TRAIN_'+case,target_identity='A25-A30-P0-smooth-conditional8:'+case,
            scope='known public TRAIN local calibration, not reconstructed unknown object',hypotheses=archive,
            material_arrays=f'../results/E2_completion/{case}.npz',full_Maxwell_predictions_same_NPZ=True,posterior_mass=None,
            represented_mass=None,unshown_modeled_mass=None,unexplored_mass_status='unquantified',outside_model_support_mass=None,
            weights_semantics='candidate_only; no residual softmax',structurally_new_worlds=0,numerical_status='full discrete Maxwell endpoint validation; no continuum certificate'))
    for pkg in packages:dump(OUT/'MaterialHypothesisPackage'/f'{pkg["observation_id"]}.json',pkg)
    ledger=[json.loads(x) for x in (OUT/'COST_LEDGER.jsonl').read_text().splitlines()]
    states=sum(r.get('amount',0) for r in ledger if r.get('event')=='STATE');rhs=sum(r.get('amount',0) for r in ledger if r.get('event')=='RHS')
    proxy=sum(r.get('proxy_likelihood_attempts',0) for r in ledger)
    spans=[r for r in ledger if r.get('event')=='ACTION'];firsttime=min(r['utc'] for r in ledger if r.get('stage')=='E2')
    failure=next(r for r in ledger if r.get('event')=='E2_JOB_FAILURE');fail_wall_upper=failure['utc']-firsttime
    successful_wall=sum(r['controlled_wall_seconds'] for r in ledger if r.get('event')=='JOB_END')
    resources=dict(stage='LXZ-20261016',new_2D_states=states,states_cap=96,new_complex_RHS=rhs,RHS_cap=3072,
        full_forward_RHS=sum(r.get('amount',0) for r in ledger if r.get('event')=='RHS' and r.get('role')=='full_forward'),
        adjoint_receiver_RHS=sum(r.get('amount',0) for r in ledger if r.get('event')=='RHS' and 'adjoint' in r.get('role','')),
        factorizations=sum(r.get('counts',{}).get('full_LU_factorizations',0) for r in spans),
        new_Mie_endpoints=0,E3_Mie_cap=64,new_3D_nonradial=0,GPU_training=0,proxy_likelihood_attempts=proxy,
        auxiliary_display_proxy_responses=sum(r.get('auxiliary_proxy_response_attempts',0) for r in ledger),
        proxy_responses_total_including_display=proxy+sum(r.get('auxiliary_proxy_response_attempts',0) for r in ledger),
        proxy_invalid_geometry_attempts=131072,proxy_numerical_repair_attempts=131072,
        measured_successful_physics_job_wall_seconds=successful_wall,failed_physics_job_wall_seconds=None,
        failed_job_conservative_elapsed_upper_seconds=fail_wall_upper,physics_wall_upper_seconds=fail_wall_upper+successful_wall,
        wall_cap_seconds=3600,wall_note='failed-job upper includes main-thread interval until failure receipt; does not claim exact pure solver time',
        all_proxy_campaign_wall_seconds=sum(json.loads((results/n/'COMPLETE.json').read_text())['wall_seconds'] for n in ['round1','round2','round1_geometry_repair','round2_geometry_repair']),
        historical_B1_384_and_P2_3459_unchanged=True,historical_credit=0,raw_failure_gap='first35 paid states lost arrays due end-of-anchor saving; ledger/source retained, not reconstructed as evidence',
        paid_states_with_saved_material_arrays=61,failed_without_endpoint_arrays=35)
    dump(OUT/'RESOURCE_SUMMARY.json',resources)
    metrics=dict(E1=summary,paired_descriptive_bootstrap=boot,E2=local,physics=resources,
        separated_modes_established=False,E3_status='NOT_RUN',A17_status='NOT_RUN_NO_QUALIFIED_COST_BENEFIT_EVIDENCE',
        N3_dev_anchor_count=0,structurally_new_P0_worlds=0,N3_BestOf4_improvement=None,scope='exposed proxy mechanism and known TRAIN local calibration')
    dump(OUT/'METRICS.json',metrics)
    ctr=json.loads((results/'PROBABILITY_CONTROLS.json').read_text());refine=json.loads((results/'E2_derivative_refinement/COMPLETE.json').read_text())
    target_registry=dict(version='LXZ-target-registry-v1',targets=[
      dict(id='A25-A30-P0',prior='original equal four-family mixture; original uniform parameter law and image peak normalization',material='chi_b+1.2*(p/maxp)*(1+.2i)/sqrt1.04',
       physics='24x24 scalarTM/VIE,8sources32receivers,complex128',noise='original complex RMS at20dB; fixed sigma; sqrt2/sigma per real channel',gain='known shared1',
       E0_development_qualified_P0_anchors=0,not_same_as_jointP2=True),
      dict(id='A25-A30-P0-smooth-conditional8',parent='A25-A30-P0',conditional=True,
       description='original smooth-family eight original uniform Gaussian coordinates; all others conditioned on public TRAIN material parameters',
       observations=['PUBLIC_TRAIN_250033_weak','PUBLIC_TRAIN_250034_ring'],known_calibration_only=True,local_chart_rank=7,
       posterior_probability_claim=False,new_material_support=False),
      dict(id='P1',status='NOT_DEFINED_NOT_RUN',reason='no posthoc support expansion'),
      dict(id='joint-p2-radial-v1:C:tensor-Chebyshev-five-per-physical-axis-radial-v1',
       prior='dx uniform [-3,1]^2 x[2.2,2.8]x[.04,.16]',decoder='original radial compact bumps, strictvacuum u<0; inactivepriorvolume retained',
       measurement='kR1.6 eightangles two polarizations,32 real channels,allRe thenallIm',noise='fixed original real sigma .05044714750236199',
       gain='shared real truncatedNormal(1,.15^2), [.7,1.3], exactmarginal logdet+truncation',
       forward='paid 5^4 Chebyshev response, no newfullMie',reference='existing adaptivephysical +continuous13/17'),
      dict(id='full-vector-Mie-matching-n14:C:bounded-shared-gain',prior='same radialprior asproxy',forward='full3Dvector layered concentric Mie nmax14',
       new_evaluations=0,reference='historical full-minus-proxy event correction onlysolid/large_cavity',not_identical_to_proxy=True)])
    dump(OUT/'TARGET_REGISTRY.json',target_registry)
    table='| 方法 | event MAE | 材料CDF grid误差 | recall* | 4096次平均arm秒 |\n|---|---:|---:|---:|---:|\n'
    for r in summary:
        if r['checkpoint']==4096:table+=f"| {r['method']} | {r['event_MAE']:.6f} | {r['CDF_grid_error_mean']:.6f} | {r['recall_at_point005']:.4f} | {r['full4096_arm_seconds_mean']:.4f} |\n"
    table+='\n*固定参考事件质量≥.05、经验估计≥.005；这是本轮false-exclusion诊断阈值，不同于旧P2的.01 discovery。CDF仅保存grid最大差。SMC为64个加权粒子，MCMC为保留拒绝重复的链平均，ESS不可直接互换。\n'
    local_table='| 校准anchor | 方法 | compatible记录/提案 | 去重compatible | 新结构 | 最大材料MAE位移 |\n|---|---|---:|---:|---:|---:|\n'
    for r in local:local_table+=f"| {r['anchor']} | {r['method']} | {r['compatible_records']}/{r['proposals']} | {r['deduplicated_compatible']} | {r['structure_qualified']} | {r['maximum_compatible_material_MAE_to_anchor']:.6g} |\n"
    text('E0_REPORT_ZH.md','''# E0：零新增Maxwell审计

三个新包已完整读取41项文本/数学/code/schema。L有65,730-byte正式完整正文，纠正X对旧supplement的历史描述。L代码仍只synthetic algebra；Z是接口和算法规格；X有通用冻结prior-pCN/独立MH与selector correction，尚无P2 adapter。未重跑旧algebra或物理campaign。

原development32个N3中只有2个进noise球，0个满足原增量phase必要条件；旧T/G full端点也不在原P0。N3–真clean白化prediction distance已缓存，范围8.7782–64.1673；clean只进入独立离线诊断，不用于proposal。现有J只在已知背景，不可称N3/P0 anchor的Jacobian。原P0硬raster参数零导数可能只表示相同材料；smooth分支另可作条件连续校准。

公开TRAIN材料250033/weak、250034/ring逐元生成重放误差0，旧observations残差22.3130/22.9696。它们是已知TRAIN校准，不是从y反演出的development解，不能填补N3 anchor缺口。E0全部state/RHS/LU/JVP/VJP/adjoint为0，G/S缓存矩阵算术与metadata wall另列。

JointP2目标、noise、boundedsharedgain与C1四维density完全匹配；C1P1二维模型不可直接用。十二事件有显著质量，但没有分离峰/路径barrier证据，不能以bins个数启动E3。缺失证据详见下文和delegated审计。''')
    text('E0_EVIDENCE_AUDIT.md', (OUT/'E0_REPORT_ZH.md').read_text()+'''

## 缺失证据

- 原P0内合格development重建及其persisted exactJ：MISSING。
- N3与合法P0参考之间的local fiber/连通性：UNRESOLVED。
- 同target两个分离重要mode及跨区概率reference：NOT_ESTABLISHED，多个event不等同多个mode。
- A17在同材料fixed-likelihood proposal上的可用ROM、可信output区间和成本优势：MISSING。
- E2故障前35状态原始endpoint arrays：MISSING，费用保持，不补跑掩盖。

可复用入口、每个源版本与NPZ header见delegated/{package_inventory,p2_assets,e0_2d}。本轮所有旧报告、真值、reference、T/G/384/P2 ledger未修改。''')
    text('LZ_X_THEORY_INTEGRATION.md','''# L/Z/X统一对象、但不混淆结论

MaterialPrior定义真实参考测度、原decoder与材料约束；MaxwellTarget固定同一forward/noise/nuisance/likelihood并计费；LocalGlobalSampler仅改变合法proposal与接受kernel。P0、条件smoothP0、未来P1、radial proxy及fullMie均独立identity。

L提供prior-informed tangent与有限噪声局部几何，要求chart光滑、材料方向非零、normal rank条件，不能把max归一化冗余当解族。Z提供全局latent proposal；由y确定的physics context没有新增观测信息，只可能改善算法表示。X负责固定材料的合法likelihood及密度修正，corrector不进入posterior kernel。

本轮local使用x=a+wΦ(u)，uniform-x prior推回standardGaussian，冻结R对称且RRᵀ+BBᵀ=I，因此pCN接受只需likelihood比。Global在dx下用q=.9normalizedC1+.1prior，接受min(1,Lnew qold/(Lold qnew))。固定.25混合两合法kernel保持target；状态相关selector需要反向选择概率比。TemperedSMC依次πβ∝p0L^β，带权重、重采样、prior-reversible mutation。

条件smoothP0的normal repair只做candidate discovery；不把多对一修正后的材料频率叫posterior。A17低秩Schur/feedback identities是可选加速关系，本轮没有合格实际增量，不强制拼成GN+Flow+MH。''')
    text('LOCAL_FIBER_RESULTS.md','# 原P0条件局部校准\n\n'+local_table+'''

已知TRAIN条件中存在多个非零但很小的noise-compatible材料位移；依原多尺度结构规则新增world为0。两anchor的8参数material Dg秩均7，1个common Gaussian-weight scale是同一材料的gauge。绝不能按latent维数报告8个物理解自由度。

最弱材料归一化方向的白化response singular scale约12.8487和297.0758；第一例转换为prior坐标的方向范数1998.55，第二例65.52。故单位材料线性位移可能需要巨大参数步，并没有可延拓的大局部fiber保证。continuation受.5 prior-coordinate trust半径后重复到相同端点，已在material-field空间去重。

第一例相对full中心差分误差8.097%，减小步长10倍后为.08097%，与二阶截断一致；精确含peak归一化Dg与原FD整体差仅8.44e−9。这说明是有限步高曲率放大，不是凭小整体Jacobian误差便可保证weak方向有限延拓。最后5个状态只做此诊断，没有重跑搜索。

E2不是N3附近反演结果；N3的合法anchor缺口仍然存在。5步链不足以估计条件posterior质量。Best-of4/SelectedMAE相对N3：NOT_RUN，而不是0改善的独立试验证明。''')
    text('GLOBAL_MODE_RESULTS.md','''# 全局材料区域与分离mode

E1在四个历史暴露观测上完成原radial C proxy target。固定twelve events的质量与CDF可以可靠对照，C1经合法MH后明显改善有限预算event误差。但这些events可能处于连续ridge，inactive系数还产生材料原像flat方向。没有保存路径barrier/独立mode证据，不能说global proposals已经跨越已证分离mode。

E3：NOT_RUN，新增Mie端点0/64。原P0没有找到观察相容的development anchor，更没有两个重要结构区域。本轮不改prior、不放宽noise、不以TRAIN已知材料冒称未知物体多模态。停止扩大Flow/Diffusion，不用额外物理寻找正结果。''')
    text('POSTERIOR_VALIDATION.md',f'''# 概率正确性与有限链误差

有限状态详细平衡核验：selector修正不变性误差{ctr['selector_correct_invariance']:.3g}，不修正偏差{ctr['selector_uncorrected_bias']:.5g}；global正确qratio误差{ctr['global_MH_correct_invariance']:.3g}，漏qratio偏差{ctr['global_omitted_q_bias']:.5g}。高ESS漏一半mass反例保留。未重跑旧合成检查。

独立审阅确认当前kernel数学target未发现新偏差。几何FD索引修复前仍是冻结prior-preserving kernel，但不能解释成Jacobian；旧两臂131072次全保留，再修复重跑131072次。采样初始点为64个prior中最佳likelihood，SMC从全部64点开始，早期差异包含初始化差异。

{table}

事件、mean、SD、CDF、quantiles均从同一chain/weighted empirical measure求取，条件参考来自原同一continuous target积分，见POSTERIOR_FUNCTIONALS.json。reference数值误差另保留；不将empirical CDF叫continuous supnorm。constant-event ESS为null，不能删掉后声称全mode混合。peak-probability excess只是signed peak差，不等同校准或错误自信发生率。

teacher/bank是旧物理+后验额外监督，C1不是免费或同信息量优势。有限4seed链不是iid exactposterior，SMC finiteN也可能漏区域。不声明posterior calibration、fullMie正确性、continuum或模型外coverage。''')
    text('A17_INCREMENTAL_VALUE.md','''# A17增量：NOT_RUN / NOT_ESTABLISHED

此次有效局部几何基线使用原full material derivative/receiver adjoint，radial几何是普通proxy Jacobian及gain Schur投影。它们不属于A17增量证据。没有同proposal、同精度且计完整setup的A17 likelihood实现，也没有output uniform bound，故没有certified early rejection。未重新跑frozenGN gap，也未修改A17冻结稿。

原P0 smooth校准规模full-LU单state很便宜；优先保留ordinary full。只有未来存在更大可行chart且具体reduced-likelihood有实测总成本优势，才值得追加A17。A17原局部GN结论保持独立。''')
    text('MATCHED_COST_COMPARISON.md','# 同proxy预算比较\n\n'+table+f'''

每个case/seed/arm共4096次likelihood（含初始化，geometry含8次差分）。1,024 checkpoint也在相同总调用位置。完整本轮proxy次数{proxy}，含两次坏geometry与两次repair；不把失败隐藏为免费。

temperedSMC批量64比逐步MH便宜，不能只比较总sample数；其eventMAE约.0161更差。C1全局MH约.00436，hybrid约.00724；混合并未自动改善。warm arm秒数不含共享cold C1导入/模型加载，后者在每case PROPOSAL.json，整轮wall包含。当前所有MH臂都预生成q draws，local/prior存在额外非必要setup；这些实测时间属于本实现，不宣称优化后部署优势。训练与bank amortized费用追溯原P2，不能算从零免费。

E2共{states}states/{rhs}complexRHS，{resources['factorizations']}次LU；原已支付cache矩阵算术另列。成功job wall{successful_wall:.4f}s，首失败job纯wall未测，保守elapsed上界{fail_wall_upper:.3f}s含主线程间隔，总physics上界{resources['physics_wall_upper_seconds']:.3f}s≤3600。所有失败/重复均保留，剩余状态0。''')
    text('FAILURE_ANALYSIS.md','''# 失败及不能宣称的事项

1. Geometry FD paired-row索引错误影响两轮机制解释；旧131072 proxy保留，单独repair131072。目标不变性并不意味着实现了正确geometry。
2. E2首任务35states/312RHS后continuation的GaussianCDF坐标达到数值端点而失败。这些35状态由于当时end-of-anchor保存设计丢失endpoint原数组；不能以ledger代替原始证据，也不重构成发现。后续改每endpoint先落盘，缩为每arm5步，56states完成；最后5states用于导数诊断，总96，未提高上限。
3. 原本material-normalized weak方向映射成1998倍prior参数方向；盲目有限步与一阶预测不相容。trust控制保护prior后多数计划幅度变成同一端点，field去重后不能把图数叫新结构数。
4. developmentN3无原P0合格anchor，不能把低MAE和noise球条件互换。已知TRAIN校准不能填补未知物体的inverse验证缺口。
5. Global MH胜过hybrid，gain-aware geometry只恢复普通local水准；当前无必须加入local/global复杂混合或A17的证据。
6. 多事件/ridge不是分离mode证明，E3不运行；fullMie budget未花。没有新训练、noise放宽、prior扩展或locked truth读取。

辅助记录写入一次stdin UTF8失败，physics/proxy费用均0。原异常日志与FAILURE_RECORDS.json保留。所有未测指标明确NOT_RUN，不以漂亮图库替代缺失证据。''')
    text('THEORY_REFINEMENTS.md','''# 来自实验的理论修订

对合法chart χ=g(u)，A=D(F∘g)。一阶weak response只控制Av；有限步还依赖D²F[Dg v,Dg v]+DF D²g[v,v]。归一化gauge应先在material metric中去掉；若材料归一化v对应极大prior坐标范数，不能从小||Av||推断大finite-noise fiber宽度。

这次smooth校准Dg秩7/8：共同Gaussian权重缩放不改变max归一化材料。第一例material-normalized v的prior范数1998.55，full secant中心差分相对误差随步长缩10倍从.08097降到.0008097。即使Dg整体误差8.44e−9，有限步高曲率仍会支配最弱方向。这是L二阶条件的真实反证式资格实验；不是新的常秩/Woodbury定理，也没有证明fold或disconnected branch。

Gain-aware proposal metric在target不变时明显改善混合，但只达到普通local水平；global correctedMDN仍较好。因此“physics geometry必然更高效”被本台架否定。X的密度修正在有限预算可产生概率改进，Z的已有提案可复用，不意味着新Maxwell物理创新。''')
    decision=dict(decision='KEEP_CORRECTED_C1_PROXY_BENCHMARK_STOP_CURRENT_LOCAL_GLOBAL_SUPERIORITY_CLAIM',
      N3_local_nontrivial_uncertainty='UNRESOLVED_NO_LEGAL_DEVELOPMENT_ANCHOR',E1_proxy='COMPLETED',
      E2='COMPLETED_KNOWN_TRAIN_CALIBRATION_ONLY',E3='NOT_RUN_NO_SEPARATED_MODE_EVIDENCE',new_structural_worlds=0,
      prior_support_relaxed=False,new_GPU_training=0,A17_increment='NOT_ESTABLISHED',TAP_strong_algorithm_paper='NO_CURRENT_EVIDENCE',
      next='restore same-target prior-compatible high-quality inverse anchor or independently justified ambiguity benchmark; no more same search under closedbudget')
    dump(OUT/'GO_NO_GO_DECISION.json',decision)
    text('GO_NO_GO_DECISION.md','''# 科研裁决

保留可求密度C1 proposal+正确MH作为受限proxy benchmark，停止当前“Local–Global优于简单global”“A17有新增成像价值”主张。原P0非平凡N3 fiber与分离模式均UNRESOLVED，不能宣称不存在其他解。E2已知TRAIN局部校准完成，新增显著结构world=0。E3 NOT_RUN，Mie端点0。所有96state预算封账，不再添seed寻找正结果。

强IEEE TAP算法论文：当前不足。最关键的缺失正实验是：同一真实prior与complex2D观测下，高质量合格inverse anchor及两个有独立证据的重要材料区域，并在同完整物理成本下改善材料风险/结构覆盖。仅减少proxy事件误差或更多图像不构成该贡献。''')
    text('NEXT_RESEARCH_DIRECTIONS.md','''# 下一步

优先解决原P0内高质量可复现inverse anchor；对hard raster family需要分层/离散prior导航，而不能假装光滑流形。另一个可授权方向是独立构造具有真mode/路径reference的Maxwell ambiguity benchmark或event-aware主动测量。必须物理说明support/noise，不能因本轮负结果临时扩宽prior。

不建议现在训练Flow/Diffusion、扩current-ROM rank或启动非同心3D。C1现有density经MH已足够解释本轮有限预算收益，先证明更复杂模块有不可替代的材料价值。任何新campaign需要新的明确授权，本轮96/96封口。''')
    text('RESEARCH_REPORT_ZH.md',f'''# LXZ：局部–全局Maxwell多解研究，2026-10-16

**本轮真正新增的显著结构材料world：0；原N3附近的合法局部fiber：仍未证实。** 已完成两个真正可运行的实验：四个暴露径向观测的概率采样比较，以及两个公开TRAIN原P0 smooth材料的full-TM局部校准。保留N3强Top1，未将TRAIN已知真材料包装成未知物体重建。

**最明确的正结果是概率修正，而非多解成像突破：** 同4096次proxy likelihood预算，C1全局MH eventMAE .004364，普通local .012039，hybrid .007243。Gain-aware geometry .011720比原geometry .026158好，但未胜过简单global。E3未获分离mode证据，新增fullMie为0。

## 从材料成像问题出发

同一电磁观测可以容许材料不确定性，但图像相似、数据noise相容、原先验合理、posterior mode是不同要求。N3低MAE值得保留，却不是原P0合法anchor的证明。旧T/G任意复像素补偿进入noise球仍破坏固定增量相位。不能从这些材料开始宣称原P0解族。

本轮完整读L/Z/X实际41份文本材料。L正式全文在本包确实存在，但只提供合成代数实现；Z给全局提案设计；X提供标准正确采样法则。统一实现没有强行加Flow或A17。代码分为真实prior、固定target和local/global kernel，corrector单独candidate-only。

## E0及目标隔离

development N3仅2/32过noise球，0/32过原phase必要条件，故无合格P0重建anchor。公开TRAIN smooth250033/weak、250034/ring合法且残差22.3130/22.9696，只能做已知材料校准。JointP2 C proxy和fullMie、原P0和条件smoothP0、未定义P1分别登记，不混比较概率。

## E1：概率机制真实运行

{table}

四观测solid/large_cavity/annular_cavity/combined已暴露，四随机seed，预算含64点初始化和8点几何差分。Global用原四维C1normalizedfullMDN的.9密度与.1prior防御混合，加准确MH ratio。Local用Gaussian-CDF prior coordinates；normal/gain-aware metric冻结。TemperedSMC64粒子，所有温度和likelihood调用计费。初始bestprior单链与SMC完整population的差异明确披露。

有correctkernel不等于finitechain精确posterior。所有CDF/mean/SD/quantiles来自同一empirical measure，与旧独立continuous参考分开比较，保存拒绝重复。多事件仍可能在一条连续ridge上，不能据此宣布跨mode成功。原reference有限数值误差保留，没有uniform fullMie/continuum保证。

## E2：真实full-Maxwell局部材料试验

{local_table}

8个原Gaussian参数只有7个材料方向，1个max-normalization gauge。允许小幅noise-compatible连续材料变化，但原多尺度结构门槛新增world0；重复field去重。普通local及geometry的最大compatibleMAE位移约.00035–.00134，continuation另一例最大.00242。这些是对knownanchor的位移，不是N3误差改进或未知物体重建。

第一例weak方向映射到prior参数范数1998.55，高曲率使大步退化。只缩小FD步长10倍，full derivative相对误差从8.097%降到.08097%；原材料Dg整体FD误差仅8.44e−9。说明localweak方向并不自动生成有意义finite material family。停止把一阶nullspace维数当多解能力。

## 失败、完整费用与边界

几何差分index错导致两原臂解释无效，保留原131072次proxy；修复后单独重跑，合计proxy {proxy}。E2首35state因chart rounding失败并丢失原始endpoint数组，账本仍保留；下一任务缩短至5步并逐endpoint保存。最后5state只做导数核验，总{states}/96 fullstates、{rhs}/3072 complexRHS，GPU0、新Mie0、非同心3D0。失败纯wall未测，成功job加保守失败elapsed上界{resources['physics_wall_upper_seconds']:.2f}s，低于3600。未借旧预算或覆盖历史。

## 八个最终问题

1. N3附近非平凡材料不确定性：**未证实**，缺同prior数据合格N3anchor；knownTRAIN只见小局部宽度。
2. A17增量：**NOT_RUN/未建立**，ordinary full-J已经足够便宜，无可信likelihood区间的早筛。
3. 其他重要结构mode：radial有多个有质量事件，**分离mode未建立**；复杂P0无合格竞争见证，非唯一性证明。
4. Physics-context latent：已有C1 density经MH有效，但它只条件于measurement，不证明新增feedback/currentcontext价值。
5. Local–Global：本轮混合劣于global，**没有独立混合优势**。Gain-aware local改进只恢复普通local级别。
6. 瓶颈：complexP0主要是合法inverseanchor、非光滑/狭窄representation及finite-chart曲率；radial是有限mixing与density近似。不是缺更大NN。
7. 下一步：优先合格anchor/有可信reference的新ambiguitybenchmark或active sensing；**暂不训练Flow/Diffusion、不扩ROM、不启动3D**。
8. 强TAP论文：**当前不足**。缺同一复杂物理target上可信重要材料区域及同完整成本的材料成像/风险正收益，不能用proxy误差或图数替代。

候选JSON/NPZ、完整MCMCchains、SMCcheckpoints及失败账本见本目录。后验未知质量与模型外质量保持null/unquantified，不把TopK重归一化成图像概率。''')
    text('START_HERE.md','''# LXZ 研究入口

2026-10-16。**显著新结构world=0，N3合法local fiber未证实；现有C1全局proposal经正确MH是本次proxy基线中最佳。**

- [中文研究报告](RESEARCH_REPORT_ZH.md) · [网页版](report.html) · [图库](gallery.html)
- [原P0局部物理结果](LOCAL_FIBER_RESULTS.md) · [概率核验](POSTERIOR_VALIDATION.md)
- [目标注册](TARGET_REGISTRY.json) · [指标](METRICS.json) · [资源](RESOURCE_SUMMARY.json)
- [失败](FAILURE_ANALYSIS.md) · [科学裁决](GO_NO_GO_DECISION.md)
- [原始数组与候选](MaterialHypothesisPackage/) · [复现](REPRODUCE.md)

物理预算96/96states、896/3072 complexRHS，包含35状态失败费用。失败前原始endpoint数组缺失明确保留；其余61物理状态数组可读。原P0/conditionalP0/radialproxy/fullMie严格隔离。E3/A17/Flow未运行。既有T/G/Prior-Constrained/P2保持历史结果与账本。

本轮只做暴露数据机制与公开TRAIN已知材料局部校准，不能称新blindvalidation或N3材料重建改进。''')
    text('REPRODUCE.md','''# 复现与费用

代码在src/lxz，冻结协议在configs。旧P2的COMMON_RESPONSE_BANK与DEVELOPMENT_CASES、ADAPTIVE_PHYSICAL_REFERENCE、PILOT_REFERENCE在既有仓库。C1权重使用c1-joint-p2@83056881的p2/results/models/full_distill-26101321.pt，原Release joint-p2-20261014可取；不得改成P1二维模型。

零物理重新出报告：PYTHONPATH=src <numpy/scipy环境python> tools/report_lxz.py；图库：tools/plot_lxz.py。只读已有rawchains/results与旧参考，不做Maxwell。

proxy.py controls为本轮必要的finite-state详细平衡检查。proxy.py run/--round2会调用proxy，完成目录禁止静默重跑；repair参数只用于已记录的索引修正。它不运行fullMie。

local2d.py、derivative_refine.py运行真实Maxwell，需要原A30/a26 source和背景缓存。用户本轮96states已全部消费，RUN_CLOSED与ledger阻止任何新增物理；勿删除/换runID重置额度。复核结果优先读NPZ，需重跑应先取得新授权。

首次失败35state原始endpoint数组缺失不补造。历史执行源码926f847/87c7680及全部异常、修复revision保留在Git和日志，rawgap不影响费用计数，不把这些状态纳入成功候选比较。''')
    receipt=OUT/'PROXY_PACKAGE_RECEIPT.json'
    if receipt.exists():
        aux=json.loads(receipt.read_text())['auxiliary_proxy_responses'];p=OUT/'RESEARCH_REPORT_ZH.md'
        p.write_text(p.read_text()+f'\n## 可交接proxy概率图\n\nMaterialHypothesisPackage/PROXY_*.json和图库提供实际correctedchains的事件代表、质量、4seed描述范围、prior、条件不确定性与带权样本。展示mass不重归一化，剩余质量保留，mode未认证。代表prediction另计{aux}次proxy response，fullphysics0。\n')
    dump(OUT/'RUN_CLOSED.json',dict(closed=True,reason='96state hard cap reached; no new meaningful structuralworld; E3 modecondition notmet',new_physics_forbidden_under_this_ledger=True,states=states,RHS=rhs))
    dump(OUT/'REPORT_BUILD.json',dict(wall_seconds=time.perf_counter()-start,physics=0,proxy_likelihood=0,source='savedarraysonly'))

if __name__=='__main__':main()
