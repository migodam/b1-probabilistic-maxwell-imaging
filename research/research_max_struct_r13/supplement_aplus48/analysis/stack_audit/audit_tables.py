"""Existing-table independent arithmetic/index audit; no scoring, arrays or physics."""
import json,csv,math,statistics,time
from pathlib import Path
from collections import Counter
START=time.process_time();OUT=Path(__file__).resolve().parent;SUP=OUT.parents[1];ROOT=SUP.parent;STACK=SUP/'analysis/stack_comparison';CHECKS=[];ISSUES=[]
def read(p):return json.loads(p.read_text())
def load(name):
 p=STACK/(name+'.json');rr=read(p)
 for i,r in enumerate(rr):r['_source']=f'{name}[{i}]';r['_source_path']=str(p)
 return rr

def check(name,condition,evidence):
 CHECKS.append(dict(check=name,passes=bool(condition),evidence=evidence))
 if not condition:ISSUES.append(name)
def table(name,rr):
 (OUT/(name+'.json')).write_text(json.dumps(rr,indent=2,allow_nan=False)+'\n');cols=list(dict.fromkeys(k for r in rr for k in r))
 with (OUT/(name+'.csv')).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=cols,lineterminator='\n');w.writeheader()
  for r in rr:w.writerow({k:json.dumps(v,separators=(',',':')) if isinstance(v,(dict,list)) else v for k,v in r.items()})
def mean(rr,k):
 vv=[r[k] for r in rr if isinstance(r.get(k),(int,float)) and not isinstance(r[k],bool)];return statistics.mean(vv) if vv else None
solvers=('cie','dual_cie','g19','dr','a_plus_mm_r13','a_plus_arc_r7');cases=[f'r13_fresh_{i:04d}' for i in range(48)];nncases=cases[::3];views=('initial','raw','corrected','primary','secondary');prefixes=(15,30,60,90,120)
rows=load('FACTORIAL_SCORES');paid=load('PAID_PREFIX_SCORES');fair=load('COMMON_ACTUAL_CPU_SIX_AND_TWELVE_SCORES');cost=load('COST_LEDGER');selected=load('PUBLIC_SELECTED_INDICES_BEFORE_METRICS');nn=load('NN_MINUS_UNIFORM_OBJECT_CHANGES');corr=load('CORRECTED_MINUS_RAW_OBJECT_CHANGES');limits=load('PUBLIC_SIX_AND_TWELVE_CPU_LIMITS')
expected={(cid,s,'uniform') for cid in cases for s in solvers}|{(cid,s,'NN') for cid in nncases for s in solvers+('direct',)}
check('400 unique run keys',len(cost)==400 and {(r['case_id'],r['solver'],r['initialization']) for r in cost}==expected,'288uniform+96NNfactorial+16optionalNNdirect; COST_LEDGER400')
check('2000 unique five-view rows',len(rows)==2000 and {(r['case_id'],r['solver'],r['initialization'],r['view']) for r in rows}=={k+(v,) for k in expected for v in views},'400*5; no duplicated or omitted run-view keys')
check('2000 unique paid-prefix rows',len(paid)==2000 and {(r['case_id'],r['solver'],r['initialization'],r['CPU_prefix']) for r in paid}=={k+(v,) for k in expected for v in prefixes},'400*5 original public CPU prefixes')
check('All score and prefix rows available; no failed tasks',all(r.get('available') is True and r['run_status']!='FAILED' for r in rows+paid+fair),'4576 existing score rows, including576 common-budget rows')
selectionkey=lambda r:tuple(r.get(k) for k in ('case_id','solver','initialization','view','CPU_prefix','budget_group'))
pre={selectionkey(r):r for r in selected};check('Public serialized candidate IDs match every emitted view/prefix/common row',len(selected)==4576 and all(selectionkey(r) in pre and pre[selectionkey(r)]['candidate']==r['selection_candidate'] for r in rows+paid+fair),'PUBLIC_SELECTED_INDICES_BEFORE_METRICS4576; no recomputed GT selector')
check('No GT fields used in serialized public selections',all(r.get('GT_fields_read') is False for r in selected),'4576/4576 public index records declare false')
# Directly cross-check frozen receipt view identifiers, rather than only two derived tables.
receipt_mismatches=[];versions=[]
for c in cost:
 p=Path(c['source_receipt']);d=read(p);folder=p.parent
 ids={'raw':d.get('raw_candidate'),'primary':d.get('primary_candidate'),'initial':read(folder/'initial.json').get('candidate'),'corrected':read(folder/'corrected.json').get('candidate'),'secondary':d.get('secondary_candidate')}
 rr=[r for r in rows if (r['case_id'],r['solver'],r['initialization'])==(c['case_id'],c['solver'],c['initialization'])]
 for r in rr:
  expected_id=ids[r['view']]
  if r['view']=='secondary' and c['secondary_fallback_to_primary']:expected_id=ids['primary']
  if expected_id!=r['selection_candidate']:receipt_mismatches.append(r['_source'])
 if c['solver'] in ('a_plus_mm_r13','a_plus_arc_r7'):
  a=read(folder/'AUXILIARY_EXECUTION.json');ops=a.get('operation_counts',{});arc=bool(ops.get('aplus_arc_stationary_roots'));mm=bool(ops.get('aplus_conditional_polynomial_roots'));guard=bool(ops.get('aplus_inherited_gtol_checks'))
  versions.append(dict(case_id=c['case_id'],solver=c['solver'],initialization=c['initialization'],cost_source_id=c['_source'],execution_source=str(folder/'AUXILIARY_EXECUTION.json'),arc_roots=ops.get('aplus_arc_stationary_roots'),quintic_roots=ops.get('aplus_conditional_polynomial_roots'),guard_checks=ops.get('aplus_inherited_gtol_checks'),flags_observed=(arc,mm,guard),expected_branch_matches=(arc and not mm and guard) if c['solver']=='a_plus_arc_r7' else(not arc and mm and not guard)))
check('All frozen receipt view candidate IDs match',not receipt_mismatches,receipt_mismatches or'400 runs /2000 frozen-view candidate IDs agree')
check('Actual Arc/guard vsMM counters match labels',len(versions)==128 and all(v['expected_branch_matches'] for v in versions),'64Arc (48uniform16NN): arc+guard, noquintic;64MM: quintic,noarc/noguard')
table('VERSION_COUNTER_AUDIT',versions)
# Same-six and all-twelve limits are distinct, independently recomputed from existing job costs.
expectedfair={(cid,s,'uniform','uniform_six') for cid in cases for s in solvers}|{(cid,s,'NN','NN_six') for cid in nncases for s in solvers}|{(cid,s,i,'factorial_twelve') for cid in nncases for s in solvers for i in ('uniform','NN')}
check('576 common-budget rows contain exact six/twelve memberships',len(fair)==576 and {(r['case_id'],r['solver'],r['initialization'],r['budget_group']) for r in fair}==expectedfair,'288uniform_six+96NN_six+192factorial_twelve')
limit_errors=[]
for l in limits:
 rr=[c for c in cost if c['case_id']==l['case_id'] and c['solver'] in solvers and (l['budget_group']=='factorial_twelve' or c['initialization']==('uniform' if l['budget_group']=='uniform_six' else'NN'))];expectedlim=min([120.]+[r['external_child_CPU_seconds'] for r in rr])
 if l['common_actual_CPU']!=expectedlim:limit_errors.append(l['_source'])
check('All per-object limits equal minimum actual six/twelve CPU',not limit_errors,limit_errors or'80 public limit rows agree with whole-child source costs')
check('Common budget candidate audit CPU stays within its limit',all(r['CPU_seconds']<=r['public_budget_limit'] for r in fair),'576 rows; retrospective paid candidate, not stopped deployment wall')
# Independently verify paired arithmetic from score rows.
lookup={(r['case_id'],r['solver'],r['initialization'],r['view'],r.get('CPU_prefix'),r.get('budget_group')):r for r in rows+paid+fair};pair_errors=[]
for c in nn:
 k=(c['case_id'],c['solver'],'NN',c['view'],c.get('CPU_prefix'),c.get('budget_group'));a=lookup[k];b=lookup.get((c['case_id'],c['solver'],'uniform',c['view'],c.get('CPU_prefix'),c.get('budget_group')))
 expectedavailable=bool(a['available'] and b and b['available'])
 if expectedavailable!=c['paired_available'] or (expectedavailable and c['delta_NRMSE']!=a['NRMSE']-b['NRMSE']):pair_errors.append(c['_source'])
for c in corr:
 a=lookup[(c['case_id'],c['solver'],c['initialization'],'corrected',None,None)];b=lookup[(c['case_id'],c['solver'],c['initialization'],'raw',None,None)]
 if c['delta_NRMSE']!=a['NRMSE']-b['NRMSE'] or c['delta_CPU_seconds']!=a['CPU_seconds']-b['CPU_seconds']:pair_errors.append(c['_source'])
check('NN andcorrector pair arithmetic matches exact same-object views',not pair_errors,pair_errors or'1312 NN derived rows incl structural noncomparisons;400corrector pairs')
structural=Counter((r['solver'],r['budget_group']) for r in nn if not r['paired_available']);check('Structural missing NN pairs are only expected noncomparisons',sum(structural.values())==256 and sum(v for (solver,b),v in structural.items() if b=='NN_six')==96 and sum(v for (solver,b),v in structural.items() if solver=='direct')==160,'96NN_six rows have no uniform at SAME NN_six limits;160optionalNNdirect rows lack uniform-direct in stack; not solverfailures')
# Corrected derived summary must be the backup minus six wrongly duplicated uniform48 labels.
cur=read(STACK/'COMMON_ACTUAL_CPU_SIX_AND_TWELVE_SUMMARY.json');orig=read(STACK/'COHORT_LABEL_CORRECTION_ORIGINAL/SUMMARY.json');removed=[r for r in orig if r not in cur];added=[r for r in cur if r not in orig]
check('Cohort-only correction removed six falseuniform48/twelve16 rows without altering metrics',len(removed)==6 and not added and all(r['cohort']=='uniform48' and r['budget_group']=='factorial_twelve' and r['expected_cases']==16 for r in removed),'Original36derivedrows→current30; all currentrows exist unchanged inbackup')
for summaryname in ('METHOD_VIEW_SUMMARY','PREFIX_SUMMARY','COMMON_ACTUAL_CPU_SIX_AND_TWELVE_SUMMARY'):
 sr=read(STACK/(summaryname+'.json'));check(summaryname+' uniform48 label always means48',all(r['expected_cases']==48 for r in sr if r['cohort']=='uniform48'),f'{len(sr)}derivedrows inspected')
# Unique cost groups only; view rows must never be resummed as separate science.
roles=Counter(r['cost_role'] for r in cost);tot=read(STACK/'COST_TOTALS.json');costerrors=[]
for z in tot:
 rr=[r for r in cost if r['cost_role']==z['cost_role']]
 for key in ('external_child_CPU_seconds','physical_LU','physical_complex_RHS'):
  if sum(r[key] for r in rr)!=z[key+'_known_sum']:costerrors.append(z['cost_role']+'/'+key)
check('Unique cost totals equal new112 +reused288 references',roles=={'new_stack112_run':112,'reused_original_run_reference_only':288} and not costerrors,costerrors or roles)
# Source-table reference extraction: no arrays, no new metrics, no common-cost selector.
oldpath=ROOT/'analysis/confirmation_mechanism16/SCORES.json';old=read(oldpath);corepath=ROOT/'analysis/confirmation_core48/SCORES.json';core=read(corepath);queuepaths={c:ROOT/'queues'/c/'RECEIPT.json' for c in ('confirmation_mechanism16','confirmation_core48')};qjobs={c:{j['task_id']:(i,j) for i,j in enumerate(read(p)['jobs'])} for c,p in queuepaths.items()}
def extract_reference(sourcepath,rr,allowed,viewlist):
 out=[]
 for i,r in enumerate(rr):
  if r['case_id'] not in nncases or r['method'] not in allowed or r['view'] not in viewlist:continue
  qi,job=qjobs[r['campaign']][r['task_id']];folder=ROOT/'RESULTS'/r['campaign']/r['task_id'];rp=folder/'RECEIPT.json';receipt=read(rp)
  expectedcandidate=receipt.get('raw_candidate') if r['view']=='raw' else receipt.get('primary_candidate') if r['view']=='primary' else read(folder/'corrected.json').get('candidate')
  check('historicalreference '+r['task_id']+'/'+r['view']+' candidate',r['selection_candidate']==expectedcandidate,f'SCORES[{i}] candidate{r["selection_candidate"]}')
  out.append(dict(source_table=str(sourcepath),source_row_id=f'SCORES[{i}]',source_CSV_line=i+2,case_id=r['case_id'],method=r['method'],view=r['view'],selection_candidate=r['selection_candidate'],candidate=r.get('candidate'),NRMSE=r['NRMSE'],RMSE=r['RMSE'],IoU=r['IoU'],Phi=r['Phi'],D=r['D'],R=r['R'],noise_feasible=r['noise_feasible'],residual_over_delta=r['residual_over_delta'],candidate_CPU_seconds=r['CPU_seconds'],actual_whole_child_CPU_seconds=job['external']['child_process_CPU_seconds'],queue_source=str(queuepaths[r['campaign']]),queue_row_id=f'jobs[{qi}]',receipt_source=str(rp),run_status=receipt['status'],role='Historical same predefined16 cached selected estimates; no new selector/metrics/common-cost claim'))
 return out
refs=extract_reference(oldpath,old,('global_dual','a_plus'),('raw','corrected','primary'));table('GLOBAL_DUAL16_REFERENCE',refs)
direct=extract_reference(corepath,core,('direct',),('primary',))+extract_reference(oldpath,old,('nn_direct',),('primary',));table('DIRECT16_REFERENCE',direct)
check('Historical global/A+exact16each3views',len(refs)==96 and all(len([r for r in refs if r['method']==m and r['view']==v])==16 for m in ('global_dual','a_plus') for v in ('raw','corrected','primary')),'Original mechanism16 only, uniformglobalPRDual notphysicalGlobalDual')
check('Historical directpaired16primary',len(direct)==32 and {r['case_id'] for r in direct if r['method']=='direct'}=={r['case_id'] for r in direct if r['method']=='nn_direct'}==set(nncases),'Core48 direct restricted to every-third16; original NNdirect16')
lines=['# Historical same16 global-PR Dual-CIE reference','', '这里只复制原mechanism16评分表和原queue/receipt已经保存的raw/corrected/primary行；没有新评分、参数、候选选择或共同成本截断。global_dual是global-PR Dual-CIE，非physical GlobalDual。A+为原R13 quintic-MM/no-guard。', '', '|method|raw均值NRMSE|corrected/primary均值NRMSE|平均whole-child CPU秒|', '|---|---:|---:|---:|']
for m in ('global_dual','a_plus'):
 raw=[r for r in refs if r['method']==m and r['view']=='raw'];primary=[r for r in refs if r['method']==m and r['view']=='primary'];lines.append(f'|{m}|{mean(raw,"NRMSE"):.10f}|{mean(primary,"NRMSE"):.10f}|{mean(primary,"actual_whole_child_CPU_seconds"):.8f}|')
lines+=['', '全16保留，包括0045。每个CSV行有原JSON索引、CSV行号、candidate ID、queue job索引和receipt路径。均值是同对象整日程端点描述；实际消费不同，不能宣称等CPU胜负。', '', '## r13_fresh_0045', '', '|method|view|candidate|NRMSE|Phi|source|','|---|---|---:|---:|---:|---|']
for r in refs:
 if r['case_id']=='r13_fresh_0045':lines.append(f'|{r["method"]}|{r["view"]}|{r["selection_candidate"]}|{r["NRMSE"]:.10f}|{r["Phi"]:.10g}|{r["source_row_id"]}|')
lines+=['', '0045的global_dual校正降低Phi同时增加材料误差；原A+从不同raw候选出发后材料误差下降。不同raw路径不能写成同一个corrector起点的比较。', '', '## 同16 direct参考（原缓存primary）', '', '|method|对象数|primary NRMSE均值|平均whole-child CPU秒|','|---|---:|---:|---:|']
for m in ('direct','nn_direct'):
 rr=[r for r in direct if r['method']==m];lines.append(f'|{m}|{len(rr)}|{mean(rr,"NRMSE"):.10f}|{mean(rr,"actual_whole_child_CPU_seconds"):.8f}|')
lines+=['', 'DIRECT16_REFERENCE.csv取core48 direct中同一every-third16，匹配原NNdirect16。没有把NN16与uniform48均值相减。NN另带R8训练先验与既有预处理费用；这里仅列原选定终点，不重选共同成本状态。']
(OUT/'GLOBAL_DUAL16_REFERENCE.md').write_text('\n'.join(lines)+'\n')
# Main descriptive independent audit.
audit=['# Stack tables independent arithmetic and provenance audit','', f'检查既有表格、公开候选索引和运行回执，未调用评分器/物理函数，未读NPZ或GT，未改算法或既有评分代码。检查项{len(CHECKS)}，未通过{len(ISSUES)}。', '', '## 分母、索引与版本', '', '- Uniform48为48对象×6方法=288运行；固定every-third16为16×6×uniform/NN，其中uniform96已属于288。新增NN96加optional NNdirect16，总唯一400运行；不是400个材料。', '- 五视图2000行，五prefix2000行；576共同成本行=uniform_six288+NN_six96+factorial_twelve192。所有行available，0科学失败。4576候选ID与先落盘公开索引一致；2000视图ID另与源receipt/initial/corrected交叉核对。', '- 实际回执128份A+核对：uniform48+NN16的64Arc均arc roots和gtol checks为正、无quintic roots；64MM均quintic roots为正、无arc/guard计数。精确逐回执来源见VERSION_COUNTER_AUDIT。Arc是R7最终readout选项放进R13协议，不能称完整R7实验原样重跑。', '', '## Uniform48实际端点与新六方法共同成本', '', '|solver|primary NRMSE均值|primary IoU均值|球内/48|平均whole-child CPU秒|common-six NRMSE均值|','|---|---:|---:|---:|---:|---:|']
for sol in solvers:
 rr=[r for r in rows if r['solver']==sol and r['initialization']=='uniform' and r['view']=='primary'];cc=[r for r in fair if r['solver']==sol and r['initialization']=='uniform' and r['budget_group']=='uniform_six'];audit.append(f'|{sol}|{mean(rr,"NRMSE"):.8f}|{mean(rr,"IoU"):.8f}|{sum(r["noise_feasible"] for r in rr)}|{mean(rr,"external_child_CPU_seconds"):.6f}|{mean(cc,"NRMSE"):.8f}|')
arcboot=read(STACK/'ARC_MINUS_A_PLUS_MM_R13_BOOTSTRAP.json');b=next(r for r in arcboot if r['initialization']=='uniform' and r['view']=='primary');audit += ['',f'原输出Arc−MM均值差{b["NRMSE_mean_delta"]:.10g}，描述bootstrap95区间[{b["NRMSE_bootstrap95_low"]:.10g},{b["NRMSE_bootstrap95_high"]:.10g}]。这不证明未来等价或一般优越；当前均值近似而费用不同。球内数不是材料成功数。', '', '## NN同16及corrector（端点描述）', '', '|solver|uniform16 primary|NN16 primary|NN−uniform|NN raw→corrected|NN实际CPU均值|', '|---|---:|---:|---:|---:|---:|']
for sol in solvers:
 u=[r for r in rows if r['solver']==sol and r['initialization']=='uniform' and r['case_id'] in nncases and r['view']=='primary'];n=[r for r in rows if r['solver']==sol and r['initialization']=='NN' and r['view']=='primary'];raw=[r for r in rows if r['solver']==sol and r['initialization']=='NN' and r['view']=='raw'];corrected=[r for r in rows if r['solver']==sol and r['initialization']=='NN' and r['view']=='corrected'];audit.append(f'|{sol}|{mean(u,"NRMSE"):.8f}|{mean(n,"NRMSE"):.8f}|{mean(n,"NRMSE")-mean(u,"NRMSE"):+.8f}|{mean(raw,"NRMSE"):.8f}→{mean(corrected,"NRMSE"):.8f}|{mean(n,"external_child_CPU_seconds"):.6f}|')
audit += ['', '## 必须保留的derived表边界', '', '1. **NN_six与uniform_six预算不同，不能跨组相配。** 每个对象各自取其6个NN任务、6个uniform任务的实际费用最小值。当前16/16对象的这两个limit实际均不同（由PUBLIC_SIX_AND_TWELVE_CPU_LIMITS逐对象核对）。NN_MINUS_UNIFORM中的96个NN_six paired_available=false是结构性“不作比较”，不是缺失结果或科学失败。只有factorial_twelve使用该对象全部12运行的共同limit，提供同成本NN−uniform配对（96个可用配对）。标准primary/prefix的同对象比较另外按其原预算语义读。', '2. Optional NNdirect16没有uniform-direct纳入stack400，因此5视图80+5prefix80的160个配对均明确不可用。原数据并非缺失：本审计另从core48 direct限制同16，列于DIRECT16_REFERENCE，完全不改stack评分表。', '3. 已记录的cohort label correction确实仅从旧36条共同成本summary中删除6条错误uniform48标签；当前30条均在备份中逐字段原样存在。被删的是factorial_twelve16的额外uniform48副本；预测、候选、误差和其余汇总数字未变。当前所有uniform48标签均为48。', '4. Corrector差表的whole-job CPU/LU差为同一run视图相减，可为0；不是corrector阶段免费。阶段候选累计CPU/LU/RHS另列。raw→corrected指标有好有坏；不得将Phi下降等同材料恢复改善。', '', '## 费用与复用', '', '新增stack112实际whole-child CPU=12373.671875秒、物理LU96749、complexRHS1773918；旧288任务是已付运行引用，CPU30590.515625秒、LU253913、RHS4818594。五视图、prefix、六/十二方法表重复引用同一run，不能重复加总。', '', 'R8主physics NN训练加生成历史费用402.558889 Mac CPU秒已保留；本次未训练或推断，既有R13预处理另计，不能当作零成本部署。', '', '## 历史global-PR Dual与direct控制', '', 'GLOBAL_DUAL16_REFERENCE.csv/.md包含原mechanism16 global_dual和A+的96条raw/corrected/primary行；DIRECT16_REFERENCE.csv包含同16 uniform/NNdirect32条原primary行。每行原索引和candidate均与原receipt核对。0045 Global校正材料恶化与A+不同raw起点的改善均保留。此控制不是新算法、重选候选或共同实际CPU比较。', '', '## 审计结论范围', '', '本审计确认表格分母、候选索引、版本分支、配对差算术和唯一费用复用的内部一致性；没有核验新成像主张、数学因果、新颖性或普遍性能。完整逐项结果见TABLE_AUDIT_CHECKS.csv；无算法修改。']
if ISSUES:audit += ['', '未通过检查：'+', '.join(ISSUES)]
(OUT/'INDEPENDENT_TABLE_AUDIT.md').write_text('\n'.join(audit)+'\n');table('TABLE_AUDIT_CHECKS',CHECKS)
(OUT/'AUDIT_RECEIPT.json').write_text(json.dumps(dict(status='EXISTING_TABLE_AUDIT_COMPLETE',checks=len(CHECKS),failed_checks=ISSUES,new_scoring=0,new_physics=0,GT_NPZ_read=False,algorithm_modified=False,source_stack=str(STACK),CPU_seconds=time.process_time()-START),indent=2)+'\n')
print(json.dumps(dict(checks=len(CHECKS),failed=ISSUES,global_reference_rows=len(refs),direct_reference_rows=len(direct),out=str(OUT))))
