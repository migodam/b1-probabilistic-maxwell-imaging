"""Completed3D tables/report from already-authorized scored CSV only;0FW/GTarrays."""
import csv,json,time,statistics,html,datetime
from pathlib import Path
ENTRY=time.perf_counter();V=Path(__file__).resolve().parent;ROOT=V.parents[1];S=V/'confirmation_postseal_scores';O=V/'completed_analysis';O.mkdir(exist_ok=True)
def readcsv(p):return list(csv.DictReader(Path(p).open()))
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def table(p,rs):
 keys=list(dict.fromkeys(k for r in rs for k in r))
 with Path(p).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rs)
def n(x):return None if x in ('',None) else float(x)
def mean(a,k):return statistics.mean(n(x[k]) for x in a if n(x.get(k)) is not None)
def med(a,k):return statistics.median(n(x[k]) for x in a if n(x.get(k)) is not None)
def fmt(x,d=4):return '—' if x is None else f'{x:.{d}f}'
methods=['sourceCIE','sourceDual','globalDual','fullGreen','contractedLagged'];names=dict(zip(methods,['CIE','Dual','Global','full-r120','contracted-r120']))
rows=readcsv(S/'MATERIAL_AND_PHYSICS.csv');prefix=readcsv(S/'COST_PREFIX_MATERIAL_SCORES.csv');assert len(rows)==50 and len(prefix)==400
D=json.loads((ROOT/'configs/vector_confirmation_material_design.json').read_text());design={d['design_id']:d for d in D['cases']}
family_names={'nonspherical_coating':'非球形包覆','open_rectangular_cavity':'开口矩形腔','adjacent_disjoint_cuboids':'相邻分离双体','necked_bridge':'有颈桥'}
for r in rows:r['family']=design[r['case_id']]['family'];r['regime']=design[r['case_id']]['regime']
metrics=['chi_relative_L2','amplitude_normalized_RMSE','fixed_IoU','Phi3','sigma_CN99_ratio'];summary=[]
for cond in ('matched12','cross18to12'):
 for m in methods:
  a=[r for r in rows if r['condition']==cond and r['method']==m];summary.append(dict(condition=cond,method=m,objects=len(a),**{k+'_mean':mean(a,k) for k in metrics},chi_relative_L2_median=med(a,'chi_relative_L2'),CN99_inside=sum(r['sigma_CN99_qualified']=='True' for r in a),corrector320=sum(int(r['corrector_accepted'])==320 for r in a),PG_inf_min=min(n(r['PG_inf']) for r in a),PG_inf_max=max(n(r['PG_inf']) for r in a)))
table(O/'ENDPOINT_SUMMARY.csv',summary)
case_rows=[];phase=[]
for r in rows:
 case_rows.append({k:r[k] for k in ['task_id','case_id','public_id','condition','family','regime','method']+metrics+['nested_CPU_seconds','nested_wall_seconds','physical_LU','infrastructure_retry','nested_cost_scope']})
 for region in ('background','low','middle','high'):
  phase.append(dict(task_id=r['task_id'],public_id=r['public_id'],family=r['family'],regime=r['regime'],method=r['method'],region=region,count=int(r[region+'_count']),bias=n(r[region+'_bias']),RMSE=n(r[region+'_RMSE']),empty_region_is_missing_not_zero=int(r[region+'_count'])==0))
table(O/'CASE_METHOD_METRICS.csv',case_rows);table(O/'PHASE_BIAS_AND_RMSE.csv',phase)
matched=[r for r in rows if r['condition']=='matched12'];seven=[r for r in matched if r['case_id']!='c17'];cross=[r for r in rows if r['condition']=='cross18to12'];pairs=[]
for scope,data in [('matched8',matched),('matched7_infrastructure_sensitivity',seven),('cross2',cross)]:
 for ref in ('sourceDual','globalDual'):
  for m in methods:
   if m==ref:continue
   ds=[];di=[]
   for group in dict.fromkeys(r['public_id'] for r in data):
    a={r['method']:r for r in data if r['public_id']==group};ds.append(n(a[m]['chi_relative_L2'])-n(a[ref]['chi_relative_L2']));di.append(n(a[m]['fixed_IoU'])-n(a[ref]['fixed_IoU']))
   pairs.append(dict(scope=scope,method=m,reference=ref,n=len(ds),chi_relative_L2_delta_mean=statistics.mean(ds),chi_relative_L2_delta_median=statistics.median(ds),chi_smaller_count=sum(d<0 for d in ds),chi_larger_count=sum(d>0 for d in ds),fixed_IoU_delta_mean=statistics.mean(di),note='Descriptive paired differences only,no significance/bootstrap/newsuccessgate'))
table(O/'ENDPOINT_PAIRWISE.csv',pairs)
ps=[]
for scope,select in [('matched8_attempt_views',lambda r:'cross' not in r['public_id']),('matched7_infrastructure_sensitivity',lambda r:'cross' not in r['public_id'] and not r['public_id'].startswith('c17_')),('cross2',lambda r:'cross' in r['public_id'])]:
 for axis,budgets in [('CPU_seconds',[480,960]),('wall_seconds',[480,960]),('physical_LU',[80,160,240,320])]:
  for b in budgets:
   for m in methods:
    a=[r for r in prefix if select(r) and r['axis']==axis and n(r['budget'])==b and r['method']==m and r['status']=='PAID_PREFIX_SELECTED']
    ps.append(dict(scope=scope,axis=axis,budget=b,method=m,n=len(a),end_to_end_qualified_n=sum(r['end_to_end_equal_cost_qualified']=='True' for r in a),**{k+'_mean':mean(a,k) for k in ('chi_relative_L2','amplitude_normalized_RMSE','fixed_IoU')},selection='Only originalsealedpublicminPhi selectors;noGT re-selection'))
table(O/'COST_PREFIX_SUMMARY.csv',ps)
paired_prefix=[]
for axis,budgets in [('CPU_seconds',[480,960]),('wall_seconds',[480,960]),('physical_LU',[80,160,240,320])]:
 for b in budgets:
  for m in methods:
   if m=='globalDual':continue
   ds=[]
   for case in [f'c{i}' for i in range(11,19) if i!=17]:
    a={r['method']:r for r in prefix if r['public_id']==case+'_matched12' and r['axis']==axis and n(r['budget'])==b};ds.append(n(a[m]['chi_relative_L2'])-n(a['globalDual']['chi_relative_L2']))
   paired_prefix.append(dict(axis=axis,budget=b,method=m,reference='globalDual',n=7,chi_relative_L2_delta_mean=statistics.mean(ds),chi_relative_L2_delta_median=statistics.median(ds),chi_smaller_count=sum(d<0 for d in ds),scope='Sameobject/hostcostprefixes,excludeallc17byinfrastructureonly'))
table(O/'MATCHED7_SAME_COST_PAIRWISE.csv',paired_prefix)
model=json.loads((S/'CROSSGRID_MODEL_ERROR.json').read_text());dump(O/'SUMMARY.json',dict(status='COMPLETED_ORIGINAL_SCORES_ONLY',endpoint_summary=summary,endpoint_pairs=pairs,cost_prefix_summary=ps,crossgrid_model_error=model,phase_bins=dict(background='truth=0',low='0<truth<=1',middle='1<truth<=3',high='3<truth<=8'),no_new_GT_array_reads=True,new_physical_RHS=0,bootstrap_calls=0,new_3D_success_gate=False))
# Markdown factual report; scientific novelty/final interpretation stays with root.
md=['# R11 真实三维有限确认：完成报告','',
'本轮50个冻结反演全部完成：8个独立材料×5方法=40个matched12任务，另2个预选材料的18→12跨网格观测×5方法=10个任务。两组分开统计；交叉观测不增加独立材料数。GT在50个预测与400个公开cost-prefix逐项封存并经root独立审查后才离线评分。未新增物理调用、算法、对象、预算或3D成功门。','',
'8例matched的材料主轮廓通常可辨认；这与冻结IoU数值、分相幅值是否准确是不同问题。full-r120相对普通Dual在5/8例的χ相对L2较小，但八例平均仍略高于Dual；有限预算下总体接近，不应说所有方法“未成像”，也不足以据此把Green分母设为新的普适默认。局部改善与退化均保留。','',
'full-r120与contracted-r120仅改材料块，current始终ordinary CIE；二者使用同一公开背景源基U=orth([Einc,G Einc])、rank120。它们不是二维G19、不是BothBlocks，也不声称rank120覆盖全dyadic G。full指材料依赖分母，非全秩G。full/contracted与Global匹配global-PR current，三个阶段的β_current均为[6,4,3]；Dual采用sourcewise-PR且β_current同为[6,4,3]，CIE采用sourcewise-PR且β_current=[6,3,1]。所有方法β_material=[6,3,1]。因此full对Global是材料块匹配对照；对Dual的差异还包含PR聚合方式，不能只归因于材料分母。','',
'## 冻结合同与指标','',
'12³、1728个实各向同性χ、5184×60复电流、180×60观测；CM+RR dyadic物理不变，FP64/complex128，m0=.25、box0..8、R0、3%独立圆复高斯噪声。MF2/3/4各max500/min20，corrector320 accepted/800requests；50/50最终均达到320 accepted上限，PG∞仍高于1e−8，不标作收敛。','',
'χ-relative L2与amplitude-normalized RMSE同时列出，避免空背景分母美化误差。fixed IoU是原先固定阈值的连续重叠指标，不补设3D成败阈值。Phi3是source-normalized拟合目标，CN99是另一个按σ定义的全观测噪声球一致性指标；CN99内不等于材料正确。','',
'## 8个matched材料终点','',
'|方法|mean χ-relative L2|median χ-relative L2|mean amp-NRMSE|mean fixed IoU|CN99内|','|---|---:|---:|---:|---:|---:|']
for r in summary:
 if r['condition']=='matched12':md.append(f"|{names[r['method']]}|{fmt(r['chi_relative_L2_mean'],6)}|{fmt(r['chi_relative_L2_median'],6)}|{fmt(r['amplitude_normalized_RMSE_mean'],6)}|{fmt(r['fixed_IoU_mean'],6)}|{r['CN99_inside']}/8|")
md+=['','下面列全部8例，不作GT选择删例。','', '|对象|几何|bulk对比度|主机|','|---|---|---|---|']
for case in [f'c{i}' for i in range(11,19)]:
 d=design[case];r=next(r for r in rows if r['case_id']==case and r['condition']=='matched12');md.append(f"|{case} {d['regime']}|{family_names[d['family']]}|{d['bulk_contrasts']}|{r['host_dispatch_role']}|")
for metric,title in [('chi_relative_L2','χ-relative L2（越低越好）'),('fixed_IoU','fixed IoU（连续值，越高越好）')]:
 md+=['','### '+title,'','|对象|'+'|'.join(names[m] for m in methods)+'|','|---|'+'---:|'*5]
 for case in [f'c{i}' for i in range(11,19)]:
  a={r['method']:r for r in rows if r['case_id']==case and r['condition']=='matched12'};md.append('|'+case+'|'+'|'.join(fmt(n(a[m][metric])) for m in methods)+'|')
md+=['','full相对Dual的局部χ误差改善出现在c12/c13/c14/c16/c17；c11/c15/c18方向相反。c16的χ误差full更小，但IoU低于Dual；c18的full IoU较高，χ误差反而更大。不同指标不能互相替代。c13开腔、c15双体、c17桥在固定切片/体图仍有可辨认主结构；c12混合包覆的内高值区明显展宽、幅度偏低，五方法χ误差都约.56。这里是具体误差描述，不加新的二元标签。','',
'## 分相偏差','',
'直接复制原评分表的分箱：background χ=0，low 0<χ≤1，middle 1<χ≤3，high 3<χ≤8。边界/有效体素也可能落入这些区间；它们不自动等于设计的bulk相。c11的4与8都归high，不能声称已经单独测得χ=4相偏差。空分箱保留缺失，不补0。全部50×4分箱见[PHASE_BIAS_AND_RMSE.csv](completed_analysis/PHASE_BIAS_AND_RMSE.csv)。','',
'以下为matched低分箱的bias（估计−真值）；完整background/middle/high bias与RMSE在原表中。','', '|对象|'+'|'.join(names[m] for m in methods)+'|','|---|'+'---:|'*5]
for case in [f'c{i}' for i in range(11,19)]:
 a={r['method']:r for r in rows if r['case_id']==case and r['condition']=='matched12'};md.append('|'+case+'|'+'|'.join(fmt(n(a[m]['low_bias'])) for m in methods)+'|')
md+=['','## 同成本前缀与工程敏感性','',
'400个原先封存的选择均为已付状态、公开minPhi，未按GT重选。c17 Dual原attempt0在315 accepted后因nvidia-smi查询超时丢进程，额外320个5184-LU已计；独立attempt1从同初值按同320/800完整重跑。该臂重跑prefix明确是attempt-only，不能当端到端成本。主材料表仍保留8例；成本敏感性整例排除c17，使用7例共同对象集。这个工程排除不是由GT决定，也不替代原8例。','',
'以下是7例共同集的mean χ-relative L2。CPU/wall是每对象在其固定主机上的预算视图，不把Mac秒与Windows秒相加称算法速度；LU为同5184阶、all60物理factor预算，额外120阶材料分母工作另有记账。','',
'|预算|'+'|'.join(names[m] for m in methods)+'|','|---|'+'---:|'*5]
for axis,budgets in [('CPU_seconds',[480,960]),('wall_seconds',[480,960]),('physical_LU',[80,160,240,320])]:
 for b in budgets:
  a={r['method']:r for r in ps if r['scope']=='matched7_infrastructure_sensitivity' and r['axis']==axis and r['budget']==b};md.append('|'+axis+' '+str(b)+'|'+'|'.join(fmt(a[m]['chi_relative_L2_mean'],6) for m in methods)+'|')
md+=['','在这个共同7例视图，full在480/960CPU和wall预算下均值较Global更差；contracted在960CPU/wall附近接近Global；320-LU视图的差异仍小，Dual均值最低。它们是冻结前缀描述，不是新停止规则或显著性结论。全部8例attempt视图、7例共同集和2cross分开的精确数值见[COST_PREFIX_SUMMARY.csv](completed_analysis/COST_PREFIX_SUMMARY.csv)。','',
'## 2个跨网格观测单列','',
'|对象|clean18与clean12相对差|CN99|解释范围|','|---|---:|---|---|']
for r in model:md.append(f"|{r['public_id']}|{100*r['clean_grid12_grid18_relative']:.4f}%|五方法均在球外|固定12³反演对18³生成的模型/离散风险，不统称分母失败|")
md+=['','|方法|mean χ-relative L2|mean amp-NRMSE|mean fixed IoU|mean CN99 ratio|','|---|---:|---:|---:|---:|']
for r in summary:
 if r['condition']=='cross18to12':md.append(f"|{names[r['method']]}|{fmt(r['chi_relative_L2_mean'],6)}|{fmt(r['amplitude_normalized_RMSE_mean'],6)}|{fmt(r['fixed_IoU_mean'],6)}|{fmt(r['sigma_CN99_ratio_mean'],6)}|")
md+=['','跨格Global的两例平均χ误差最低，Green两臂没有改善这个平均；不能将这两例并入matched成功比例。c13的17.68% clean forward差远高于3%噪声，c18的3.53%也非零；观测模型差、有限预算与算法误差在这里并存。c18跨格图出现亮颈/空心化等材料分配偏差，IoU相近仍可能有显著χ误差。','',
'## 成本与失败保留','',
'统一总账[COSTS.json](COSTS.json)覆盖本worker负责的全部vector生产者：公开basis/导数/固定状态、DEV12、一次6LU库选择、10观测生成、50确认、原失败与重跑、controller/helper/transport、评分与图、这份机械汇总；不是仅fresh。root环境清单、root理论/发布与2D生产者明确排除。','',
'已知CPU下界按平台分开：Mac @MAC_CPU_COST@ s，Windows @WINDOWS_CPU_COST@ s。全vector 5184阶物理LU20507次；17496阶生成为2个all60逻辑批/120RHS/6392 GMRES迭代，未从RHS别名重复相加。小32/120阶分母LU独立列示，不混入5184物理LU。wall不可跨并行lane相加为总历时，CUDA event区间也不是活跃kernel时间。','',
'原GPU failed phase包含15个已封存臂+第16臂315步，替代这些嵌套arm计时；独立重跑phase额外计费。原315目录归archived attempt0，sealed重跑attempt1入canonical，下一臂mkdir sentinel受控0FW停止，后9新进程用30s有效空闲资源快照缓存。查询等待改善仅为工程资源协作，不算算法速度。','',
'查询、导入/退出、部分旧remote pack/extract等未测量部分在null gaps中保留。所有历史日志/非零exit都保留，受控controller drain/c18 waiting取消不误列科学失败。Mac20个临时spool均在完整原数组归档后删除，当前临时目录为空。[ATTEMPT_PATH_REGISTRY.json](ATTEMPT_PATH_REGISTRY.json)给出路径语义。','',
'## 本轮收束建议','',
'保留三维ordinary基线与两种rank120材料分母的局部差异证据；本轮平均值与同成本前缀不支持把Green材料分母直接替换成统一默认。跨网格观测按明确模型/离散风险单列。这个判断只覆盖规定的8个材料与冻结预算，不推广到未知总体、全G近似、BothBlocks或收敛极限。本轮到此结束，不追加算法/对象/物理试验。理论与原创性及综合论文结论由root审查。','',
'## 全30图入口','',
'每个体图为真3D体素外表面χ>.5、统一χ=0..8与相机；内部材料分相须结合固定中心正交切片和区域误差，体图遮挡不证明内核准确。所有8matched与2cross、五方法及GT均保留。']
figs=sorted((S/'figures').glob('*.png'));assert len(figs)==30
for p in figs:md.append(f"- [{p.stem}](confirmation_postseal_scores/figures/{p.name})")
md+=['','原始评分[主表CSV](confirmation_postseal_scores/MATERIAL_AND_PHYSICS.csv)、[400前缀CSV](confirmation_postseal_scores/COST_PREFIX_MATERIAL_SCORES.csv)、[失败/预算类](confirmation_postseal_scores/FAILURE_CLASSES.csv)、[public SELECTIONS](confirmation_host_union/SELECTIONS.json)、[完成复现说明](COMPLETED_REPRODUCE.md)。','']
(V/'REPORT.md').write_text('\n'.join(md))
# Self-contained static HTML uses only local evidence images/links.
def ht_table(headers,data):return '<div class="scroll"><table><thead><tr>'+''.join('<th>'+html.escape(str(x))+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(str(x))+'</td>' for x in r)+'</tr>' for r in data)+'</tbody></table></div>'
body=['<h1>R11真实三维有限确认：50任务完成</h1>','<p>8 matched材料与2 cross观测分开；全50预测/400前缀经独立公共审查后授权离线GT评分。所有320 accepted cap，不标作收敛，不补设3D成功门。</p>',
'<p>full-r120相对普通Dual在5/8例χ误差略小，但8例平均Dual仍最低。主轮廓可辨认、冻结IoU和材料幅值误差是不同层面；有限预算接近不等于“所有方法未成像”。Green只改材料，current仍ordinary CIE；full/contracted匹配Global的global-PR current，Dual是sourcewise-PR。不是二维G19或BothBlocks；full对Dual的差异不能只归因于材料分母。</p>',
'<nav><a href="REPORT.md">完整中文Markdown</a> · <a href="COMPLETED_REPRODUCE.md">完成复现说明</a> · <a href="COSTS.json">全部vector去重总账</a> · <a href="ATTEMPT_PATH_REGISTRY.json">attempt路径</a> · <a href="confirmation_host_union/SELECTIONS.json">50公开选择</a></nav>',
'<h2>Matched8连续终点</h2>',ht_table(['方法','mean χ-relative L2','mean amp-NRMSE','mean fixed IoU','CN99内'],[[names[r['method']],fmt(r['chi_relative_L2_mean'],6),fmt(r['amplitude_normalized_RMSE_mean'],6),fmt(r['fixed_IoU_mean'],6),str(r['CN99_inside'])+'/8'] for r in summary if r['condition']=='matched12']),
'<p>c13/c16等局部改善保留，c12混合包覆的幅值/层次模糊也保留。完整逐对象、分相bias与同成本前缀在Markdown和机器表中。</p>',
'<h2>Cross2：明确模型/离散风险</h2>',ht_table(['对象','clean18 vs clean12','CN99'],[[r['public_id'],f"{100*r['clean_grid12_grid18_relative']:.4f}%",'五方法均球外'] for r in model]),
ht_table(['方法','mean χ-relative L2','mean fixed IoU','mean CN99 ratio'],[[names[r['method']],fmt(r['chi_relative_L2_mean'],6),fmt(r['fixed_IoU_mean'],6),fmt(r['sigma_CN99_ratio_mean'],6)] for r in summary if r['condition']=='cross18to12']),
'<p>不与matched合并。不能把17.68%/3.53%的clean forward模型差统称为优化或分母失效；未膨胀3%噪声。</p>',
'<h2>成本与范围</h2><p>全部vector已知CPU下界：Mac @MAC_CPU_COST@ s；Windows @WINDOWS_CPU_COST@ s。平台分账，wall不可并行相加。5184物理LU20507；fine17496为2逻辑批/120RHS，6392 GMRES迭代；320额外失败LU已含一次。</p>',
'<p>原315 query-timeout失败、独立320重跑、0FW mkdir drain、30s资源缓存、c18往返调度与spool均保留。只有c17 Dual重跑prefix是attempt-only；同成本敏感性整例排c17，主8表不删例。全部具体null缺口及排除生产者见COSTS。</p>',
'<h2>机器表</h2><ul>'+''.join(f'<li><a href="completed_analysis/{p.name}">{p.name}</a></li>' for p in sorted(O.glob('*.csv')))+ '</ul>',
'<h2>全部30图</h2><p>统一物理色标χ0..8；体图χ>.5只显示外表面，内部须看正交切片。点击图可查看原尺寸；没有删去失败或较差案例。</p>']
for p in figs:body.append(f'<section><h3>{html.escape(p.stem)}</h3><a href="confirmation_postseal_scores/figures/{p.name}"><img loading="lazy" src="confirmation_postseal_scores/figures/{p.name}" alt="{html.escape(p.stem)}"></a></section>')
body+=['<footer>本轮到此结束，不新增科学调用。综合科学/原创性判断由root完成；本HTML为本地归档，最终发布由root统一处理。</footer>']
(V/'REPORT.html').write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>R11真实三维完成报告</title><style>body{max-width:1350px;margin:36px auto;padding:0 24px;font:16px/1.65 system-ui;color:#17202a;background:#fafbfc}h1,h2{line-height:1.3}h2{margin-top:40px}nav,a{color:#14559a}nav{padding:16px;background:#eef4fb}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;background:white}td,th{padding:9px 12px;border:1px solid #d9e0e7;text-align:left}img{width:100%;height:auto;background:white}section{margin:35px 0;padding:12px;background:white;border:1px solid #dde3e8}footer{padding:25px;border-top:1px solid #ccd3da}</style><body>'+''.join(body)+'</body></html>')
repro='''# R11 vector完成复现与历史边界

本文件描述已经完成的生产者与只读复算入口。本轮科学已结束；不在当前目录重新启动任何Maxwell反演、benchmark、corrector扩展或新设计。建议在归档副本复算，保留原始预测与public seals。

## 真实入口与冻结输入

工作树：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r11_worktree`。实际Mac Python：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python`；Windows：`D:/python/python.exe`，历史remote repo `D:/AI/MAX_STRUCT_R10_VECTOR_20261107`。

真实依赖为现有R10 `vector_grid/vector_cuda_model/vector_physical/vector_cie/vector_cuda_dense/vector_resources`及R9资源job helper；R10不改写。物理是12³ dyadic CM+RR，3N interleaved电流与每cell一个实χ。新代码只在R11 `src/vector_denominator.py`、`src/vector_cie.py`材料回调、薄backend/storage/resource包装；forward未改变。

正式合同：`configs/vector_confirmation_contract.json`、`vector_confirmation_material_design.json`、`vector_confirmation_queue.json`。源基在Y生成前由公开Einc与G Einc构造，relative SVD阈值1e−10、实际rank120；不是全G谱截断。源基与传播遗漏诊断、唯一必要full/lagged/z新梯度接线检查均在本目录独立receipt保留，不需重复。

## 已完成科学次序与源码证据

1. ordinary DEV6（default库）→一次公开3状态×default/cusolver=6LU完整value/gradient/原方程残差→new DEV6（选cusolver）→fresh50。probe嵌套在new DEV6，库差不是纯算法速度。DEV12先封存后仅旧曝光GT评分。
2. 新观测生成：8个12³ dense-LU/all60与预选c13/c18的18³ GMRES/all60；3% iid CN无样本重归一化。10个public observations只含Y与公开noise metadata；private generation仅写GT，反演从未读GT。
3. host最初50GPU计划，后25/25，再基于已付冷首臂38.6min改35GPU/15Mac；Mac15完成后按实际10个spool臂中位21.11min和条件授权将完全未启c18 matched+cross10转回Mac。最终25/25，每对象五方法与两grid同host。所有旧manifest/config版本与时刻/预测保留；最终实际成本替代早期ETA。
4. Mac c11五臂旧RAM版自然结束；仅空闲controller受控drain，c14/c15与随后c18两组用单临时memmap spool保存历史J/Q视图，数学/FP64/cache-hit/预算不变。20个spool均在原数组/前缀封存后释放并删临时文件，原输出完整。
5. GPU original25完成15臂后，c17 Dual在315 accepted/317完成requests之后、下一FW之前nvidia-smi10s查询超时。原attempt0额外320×5184LU/19200primal/19020adj RHS全计；80/160/240快照保留，315态与拟牛顿历史未持久化。root授权一次同初值同预算独立重跑。
6. retry1 c17 Dual自然320封存；下一c17Global路径预建明确sentinel，既有exclusive mkdir在其初始FW前受控结束。finally全部LU/primal/adj计数等于唯一已完成arm，证明下一臂精确0FW。raw FileExistsError不是算法失败。
7. 后9原未启任务在新进程：仅查询超时杀自建nvidia-smi helper并有限重试；成功external=[]快照缓存30s，过期/外部出现先暂停并重查，最大新增外部探测延迟30s。资源等待改善不是科学算法提速。新增seal+归档后stop-request接口，本轮未使用该停止请求。
8. 全50预测/400prefix公共联合封存后，root独立逐值审查并于UTC2026-11-07T07:53:45.607075写`audit/VECTOR50_ROOT_POSTSEAL_AUTHORIZATION.json`，之后才运行唯一fresh scorer。原public seals的GTfalse保持历史记录，不回写。

源码/输入快照与成本证据在`remaining_transport/SOURCE.zip`及其各overlay、`c18_append_transport/SOURCE.zip`、`gpu_recovery_attempt1_transport/SOURCE.zip`、`gpu_recovery_attempt2_transport/SOURCE.zip`、各PROCESS/FROZEN_CONFIG/SOURCE_EVIDENCE。source快照/size/mtime与registry保留，不新增hash检查。

## 只读复算与离线产物

原评分数据已在`confirmation_postseal_scores/`：50主表、400冻结前缀、已付轨迹、固定accepted80/160/240/320、CN99与物理残差、全30图。fresh scorer真实入口为`experiments/vector_postseal_score.py`；外部wait4入口`experiments/vector_fresh_postseal_external.py`，必须有后来root独立授权与all50公共barrier。当前已完成，不在原输出目录再次跑；脚本exclusive输出保护会拒绝覆盖。

只从已评分CSV复算机器汇总/本报告与去重账本（0新FW/0GT数组访问），在归档副本的worktree中：

```sh
VECTOR_PYTHON="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python"
"$VECTOR_PYTHON" research/research_max_struct_r11/delegated/vector/compile_completed_analysis.py
"$VECTOR_PYTHON" research/research_max_struct_r11/delegated/vector/finalize_completed_costs.py
```

`compile_completed_analysis.py`只读已授权评分CSV/公开合同，输出completed_analysis与REPORT；不按GT重选状态，不新增成功门。`finalize_completed_costs.py`调用原metadata costcompiler并纳入本次报告生产者，所有重复元数据编译run的CPU历史保留。总账scope含本worker全部vector DEV/setup/diagnostic/generation/fresh/失败/transport/helper/score/plot/summary；root环境清单与2D/理论/发布单列排除。原failed wholephase替代15个arm而非叠加；Macwait4含子孙/RSShelper；含null gaps的平台已知CPU下界不称精确总费，Mac/Windows秒不相加称算法速度。5184 full-LU与32/120 small-LU、17496 iterative工作量分开。

科学历史runner的实际接口是`experiments/vector_run.py --config ... --queue ... --offset ... --max-tasks ... --output ...`。原50、Mac分组、recovery1/2的真实命令已完整保存在各controller PROCESS_RECEIPTS/ACTIVE、transport manifests与FROZEN_CONFIG；历史transport脚本固定exclusive输出，不能在已完成archive里直接重跑。若未来明确要求独立科学重现，须复制成隔离输出并保持相同公式、初值、数据、rank120、预算/FP64；本轮没有授权追加这种重跑。

## 路径与统计限制

`ATTEMPT_PATH_REGISTRY.json`将原c17 Dual历史canonical路径解析到archived attempt0，将recovery1 sealed路径映射到当前canonical；历史receipt路径原样保留。recovery2是后9个原未启动任务的执行代际，不是九次科学重跑。共50完成预测+1个基础设施失败attempt，0FW sentinel/取消waiting controller不算科学attempt。

8matched材料是规定设计有限集，不是iid总体；cross只重复其中2个材料，不扩N。主8表不删c17；同成本敏感性整例排c17，因唯一retry prefix非端到端。全部连续指标、分箱bias与图保留；CN99资格、IoU连续值、主轮廓可辨认与幅值准确分别解释。没有新科学调用、bootstrap或新3D成功阈值；本轮分析与交付后结束。
'''
(V/'COMPLETED_REPRODUCE.md').write_text(repro)
review=dict(status='VISUAL_REVIEW_RECORDED',orthogonal_images_reviewed=[p.name for p in figs if p.name.endswith('_orthogonal.png')],volume_images_reviewed=['c12_matched12_volume.png','c17_matched12_volume.png','c18_cross18to12_volume.png'],observed='All10orthogonal plates useGT+5methods andfixed0..8;reviewedvolumes true3D,externalfacescanoccludeinternalcore;mainstructuresvisiblewithcase-specificblur/phaseerror',all30_files_linked=True,no_human_acceptance_or_newbinarysuccess_claim=True)
dump(O/'VISUAL_REVIEW.json',review)
record=dict(UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-ENTRY,scope='Own CSV-only completedsummaries/REPORTmd+HTML/COMPLETED_REPRODUCE generation;noGTarrays/newphysics/bootstrap',new_forward_RHS=0,GT_array_reads=0)
with (O/'ANALYSIS_RUNS.jsonl').open('a') as f:f.write(json.dumps(record)+'\n')
dump(O/'ANALYSIS_COST.json',record);print(json.dumps(dict(status='COMPLETED_TABLES_REPORT_AND_ALL30_LINKS',CPU_seconds=record['CPU_seconds'],output=str(V))),flush=True)
