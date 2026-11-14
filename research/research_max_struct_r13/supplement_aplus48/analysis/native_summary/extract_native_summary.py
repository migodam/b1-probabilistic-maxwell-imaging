"""Read precomputed JSON/CSV diagnostics only; no arrays, models or score functions."""
import csv,json,math,statistics,time
from pathlib import Path
from collections import Counter,defaultdict
START=time.process_time();OUT=Path(__file__).resolve().parent;SRC=OUT.parent/'native_states';INDEX=OUT.parent/'comparison48/RUN_INDEX.json'
SOURCES={};TABLES={}
def load(name):
 p=SRC/(name+'.json');data=json.loads(p.read_text());SOURCES[name]=str(p)
 for i,r in enumerate(data):r['_source_id']=f'{name}[{i}]';r['_CSV_line']=i+2
 TABLES[name]=data;return data
coverage=load('TASK_COVERAGE');targets=load('RAW_TARGET_CHECKS');trans=load('READOUT_TRANSITIONS');force=load('CURRENT_FORCE_TRANSPORT');fields=load('FIELD_DEFECTS');dense=load('DENSE_NATIVE_HISTORY');audits=load('PAID_PHYSICAL_AUDITS');snapshots=load('SNAPSHOTS');issues=load('ISSUES')
receipt=json.loads((SRC/'RECEIPT.json').read_text());index=json.loads(INDEX.read_text());cases=sorted({r['case_id'] for r in coverage})
def key(r):return r['case_id'],r['method']
def num(x):return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)
def vals(rr,k):return [r[k] for r in rr if num(r.get(k))]
def stats(rr,k,prefix=None):
 v=vals(rr,k);p=prefix or k
 return {p+'_n':len(v),p+'_min':min(v) if v else None,p+'_median':statistics.median(v) if v else None,p+'_max':max(v) if v else None}
def sign(x):return 'missing' if not num(x) else 'negative' if x<0 else'positive' if x>0 else'zero'
def counts(rr,k,p):
 c=Counter(sign(r.get(k)) for r in rr);return {p+'_'+s:c[s] for s in ('negative','zero','positive','missing')}
def refs(rr):
 ids=[r['_source_id'] for r in rr];return ids

def compressed(rr):
 # JSON element ranges and CSV line ranges; no claim about missing intermediate rows.
 group=defaultdict(list)
 for r in rr:
  name,i=r['_source_id'].split('[');group[name].append(int(i[:-1]))
 output=[]
 for name,ii in group.items():
  ii=sorted(ii);ranges=[]
  for x in ii:
   if ranges and x==ranges[-1][1]+1:ranges[-1][1]=x
   else:ranges.append([x,x])
  output.extend(dict(table=name,json_index_range=[a,b],CSV_line_range=[a+2,b+2]) for a,b in ranges)
 return output

def table(name,rr):
 (OUT/(name+'.json')).write_text(json.dumps(rr,indent=2,allow_nan=False)+'\n');columns=list(dict.fromkeys(k for r in rr for k in r))
 with (OUT/(name+'.csv')).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=columns,lineterminator='\n');w.writeheader()
  for r in rr:w.writerow({k:json.dumps(v,separators=(',',':')) if isinstance(v,(dict,list)) else v for k,v in r.items()})

groups={name:defaultdict(list) for name in TABLES}
for name,rr in TABLES.items():
 for r in rr:groups[name][key(r)].append(r)

rawrows=[];overview=[];denseout=[];actualdetails=[];forcedetails=[];defectout=[];selection=[];extremes=[]
idx={(r['case_id'],r['method']):(i,r) for i,r in enumerate(index['runs'])}
for cov in coverage:
 k=key(cov);rr=groups['READOUT_TRANSITIONS'][k];ff=groups['CURRENT_FORCE_TRANSPORT'][k];tt=groups['RAW_TARGET_CHECKS'][k];dd=groups['DENSE_NATIVE_HISTORY'][k];aa=groups['PAID_PHYSICAL_AUDITS'][k]
 common=dict(case_id=k[0],method=k[1],task_id=cov['task_id']);rawrows.append(dict(common,coverage_status=cov['status'],requested_target_iterations=cov['requested_target_iterations'],reached_target_iterations=cov['reached_target_iterations'],unreached_target_iterations=cov['unreached_target_iterations'],target_rows=len(tt),target_status_counts=dict(Counter(x['status'] for x in tt)),exact_material_matches=sum(x.get('material_array_equal') is True for x in tt),target_request_matches_raw=all(x.get('request_matches_original_raw_iteration') is True for x in tt) if tt else None,material_max_abs_difference_max=max(vals(tt,'material_max_abs_difference')) if vals(tt,'material_max_abs_difference') else None,full_auxiliary_state_equality_established=False,source_coverage_id=cov['_source_id'],source_target_ids=refs(tt),source_rows=[{z:x.get(z) for z in ('_source_id','_CSV_line','requested_iteration','original_raw_iteration','original_raw_candidate','status','observer_declares_reached','target_snapshot_count','material_array_equal','material_max_abs_difference','material_RMS_difference','snapshot_id','qualification')} for x in tt],missing_reason='No raw iteration requested for this task; no fidelity claim' if not tt else None,replay_termination=cov.get('replay_termination')))
 o=dict(common,coverage_status=cov['status'],saved_snapshots=cov['saved_snapshots'],requested=cov['requested_target_iterations'],unreached=cov['unreached_target_iterations'],readout_rows=len(rr),readout_status='NOT_APPLICABLE_DR_NO_CIE_READOUT' if k[1]=='dr' else'AVAILABLE' if rr else'MISSING',current_force_rows=len(ff),field_rows=len(groups['FIELD_DEFECTS'][k]),native_history_rows=len(dd),native_paid_audit_rows=len(aa),source_coverage_id=cov['_source_id'],readout_source_ids=refs(rr),force_source_ids=refs(ff),fidelity_source_ids=refs(tt),scope='8 exposed selected objects; retained snapshots not uniform sample of iterations; DR has no CIE readout/PR rows')
 o.update(counts(rr,'actual_delta_Phi','actual_physical_DeltaPhi'));o.update(counts(rr,'physical_first_order_step','physical_linear_prediction'));o.update(counts(rr,'counterfactual_ordinary_loss_change','counterfactual_ordinary_Deltaloss'));o.update(counts(rr,'counterfactual_onehop_loss_change','counterfactual_onehop_Deltaloss'))
 for metric in ('actual_delta_Phi','physical_first_order_step','material_step_norm','cosine_step_negative_physical_gradient','cosine_step_negative_ordinary_gradient','cosine_step_negative_onehop_gradient'):o.update(stats(rr,metric))
 matrix=Counter(sign(r.get('postcurrent_prereadout_gradient_dot_native_direction'))+'→'+sign(r.get('postreadout_gradient_dot_native_direction')) for r in ff);o['PR_force_direction_sign_matrix']=dict(matrix)
 for metric in ('normalized_prereadout_gradient_dot_native_direction','normalized_postreadout_gradient_dot_native_direction','cosine_postreadout_prereadout_current_gradients','gradient_transport_relative_to_prereadout'):o.update(stats(ff,metric))
 overview.append(o)
 dense_row=dict(common,rows=len(dd),streams=dict(Counter(r.get('stream') for r in dd)),stage_counts=dict(Counter(str(r.get('stage')) for r in dd)),source_ranges=compressed(dd),CPU_values_known=len(vals(dd,'CPU_seconds')),CPU_null_rows=sum(r.get('CPU_seconds') is None for r in dd),direction_restart_true=sum(r.get('direction_restarted') is True for r in dd),direction_restart_false=sum(r.get('direction_restarted') is False for r in dd),PR_plus_positive=sum(num(r.get('PR_plus')) and r['PR_plus']>0 for r in dd),scope='Native scalar auxiliary history, not original physical Phi; per-update wall is not process CPU; DR sparse audit records not every-step objective')
 for metric in ('current_step_objective_change','readout_material_change','PR_plus','data_term','state_term_current','state_term_material','regularizer','shadow_gap'):
  dense_row.update(stats(dd,metric))
 dense_row.update(counts(dd,'current_step_objective_change','native_current_objective_change'));denseout.append(dense_row)
 for r in rr:
  actualdetails.append(dict(common,source_row_id=r['_source_id'],source_CSV_line=r['_CSV_line'],**{z:r.get(z) for z in ('snapshot_id','stage','iteration','global_iteration','snapshot_CPU_seconds','retained_slots','status','before_evaluation_id','after_evaluation_id','material_step_norm','actual_delta_Phi','physical_first_order_step','cosine_step_negative_physical_gradient','cosine_step_negative_physical_pg','counterfactual_ordinary_loss_change','counterfactual_onehop_loss_change','ordinary_first_order_step','onehop_first_order_step','cosine_step_negative_ordinary_gradient','cosine_step_negative_onehop_gradient','cosine_ordinary_onehop_gradients')},actual_physical_change_sign=sign(r.get('actual_delta_Phi')),physical_first_order_sign=sign(r.get('physical_first_order_step')),material_gradient_scope='Ordinary and onehop are fixed-current counterfactuals, including on G19 trajectories; not native G19 derivative'))
 for r in ff:
  forcedetails.append(dict(common,source_row_id=r['_source_id'],source_CSV_line=r['_CSV_line'],**{z:r.get(z) for z in ('snapshot_id','stage','iteration','global_iteration','snapshot_CPU_seconds','retained_slots','status','beta_current','postcurrent_prereadout_gradient_norm','postreadout_current_gradient_norm','native_precurrent_gradient_norm','native_current_direction_norm','current_force_transport_norm','cosine_transport_native_direction','cosine_postreadout_prereadout_current_gradients','postcurrent_prereadout_gradient_dot_native_direction','postreadout_gradient_dot_native_direction','normalized_prereadout_gradient_dot_native_direction','normalized_postreadout_gradient_dot_native_direction','gradient_transport_relative_to_prereadout')},prereadout_direction_sign=sign(r.get('postcurrent_prereadout_gradient_dot_native_direction')),postreadout_direction_sign=sign(r.get('postreadout_gradient_dot_native_direction')),scope='Same saved native direction; recomputed forces fixed at postcurrent J/Q/A and beta_current; saved native gradient is pre-current; not next PR iterate'))
 for block in sorted({r['field_block'] for r in groups['FIELD_DEFECTS'][k]}):
  fr=[r for r in groups['FIELD_DEFECTS'][k] if r['field_block']==block];z=dict(common,field_block=block,rows=len(fr),source_row_ids=refs(fr),scope='Observed saved field block vs same paid Jphys(m); directional norm ratio not resolvent operator norm')
  for metric in ('u_norm','v_norm','constitutive_plus_propagation_forcing_norm','auxiliary_data_residual_over_delta','current_to_physical_norm','port_current_to_physical_over_delta','physical_residual_over_delta','directional_amplification'):z.update(stats(fr,metric))
  for metric in ('port_current_to_physical_over_delta','physical_residual_over_delta','directional_amplification'):
   candidates=[r for r in fr if num(r.get(metric))];high=max(candidates,key=lambda r:r[metric]) if candidates else{};z[metric+'_max_source_id']=high.get('_source_id');z[metric+'_max_snapshot_id']=high.get('snapshot_id')
  defectout.append(z)
 for polarity in ('min','max'):
  candidates=[r for r in rr if num(r.get('actual_delta_Phi'))]
  if candidates:
   r=(min if polarity=='min' else max)(candidates,key=lambda x:x['actual_delta_Phi']);extremes.append(dict(common,example_rule=polarity+' actual_delta_Phi among retained snapshots; not selected by GT',source_row_id=r['_source_id'],source_CSV_line=r['_CSV_line'],**{z:r.get(z) for z in ('snapshot_id','global_iteration','actual_delta_Phi','physical_first_order_step','counterfactual_ordinary_loss_change','counterfactual_onehop_loss_change','material_step_norm','before_evaluation_id','after_evaluation_id')}))
 oldpair=idx.get(k)
 if oldpair:
  i,old=oldpair;v=old['view_metadata'];z=dict(common,original_run_index_source=str(INDEX),source_row_id=f'RUN_INDEX/runs[{i}]',original_run_status=old['run_status'],original_view_candidate_ids=old['view_candidate_ids'],original_secondary_fallback=old['secondary_fallback_to_primary'],original_external_CPU_seconds=old['external_cost_source'].get('external_child_CPU_seconds'),replay_paid_audit_rows=len(aa),replay_audit_source_ranges=compressed(aa),replay_audit_phase_counts=dict(Counter(r.get('phase') for r in aa)),raw_target_source_ids=refs(tt),scope='Original frozen public view metadata; replay audits separate. No new choice or GT metric. Corrector endpoint may differ in timed replay; raw material match alone does not validate replay corrected endpoint.')
  for view in ('initial','raw','corrected','primary','secondary'):
   for metric in ('candidate','phase','CPU_seconds','Phi','D','R','residual_over_delta','noise_feasible','legal','eligible'):z[view+'_'+metric]=v.get(view,{}).get(metric)
  for a,b in [('raw','corrected'),('raw','primary'),('corrected','primary')]:
   for metric in ('Phi','D','R','CPU_seconds'):
    x,y=v.get(a,{}).get(metric),v.get(b,{}).get(metric);z[a+'_to_'+b+'_delta_'+metric]=y-x if num(x) and num(y) else None
  selection.append(z)
for name,rr in [('CASE_METHOD_SUMMARY',overview),('RAW_FIDELITY',rawrows),('ACTUAL_READOUT_DETAILS',actualdetails),('CURRENT_FORCE_DETAILS',forcedetails),('FIELD_BLOCK_SUMMARY',defectout),('DENSE_HISTORY_SUMMARY',denseout),('SELECTION_CORRECTOR_PUBLIC',selection),('READOUT_EXTREME_EXAMPLES',extremes)]:table(name,rr)
# Only the exact row subsets used in summaries are reproduced; no new scoring.
table('SOURCE_SCHEMA',[dict(table=name,source_json=SOURCES[name],source_csv=str(SRC/(name+'.csv')),rows=len(rr),columns=sorted({k for r in rr for k in r if not k.startswith('_')}),row_id_rule=f'{name}[zero-based JSON index]; CSV line=index+2') for name,rr in TABLES.items()])
actual=Counter(sign(r.get('actual_delta_Phi')) for r in trans);cf_down_up=sum(r.get('counterfactual_onehop_loss_change',0)<0 and r.get('actual_delta_Phi',0)>0 for r in trans);cfOrd_down_up=sum(r.get('counterfactual_ordinary_loss_change',0)<0 and r.get('actual_delta_Phi',0)>0 for r in trans)
lines=['# 原生状态证据摘录（仅既有表格汇总）','',f'对象为事后指定的8个既有对象、40个case/method任务；原生分析保留561个snapshot、475个CIE系readout transition、475个current-force行、1208个field-block行。汇总未读取GT/NPZ、未新增forward/adjoint/optimizer/评分。原native analyzer此前已付844次LU、18924 complex RHS；不将这些当成本次新开销。','','## 一、回放对齐不是全状态等价','',f'39个已请求raw target中，36个是材料数组逐单元完全相等，3个UNREACHED。另有1任务未请求raw target：r13_fresh_0005__dr，不能把它计为通过或失败。40任务AVAILABLE仅代表轨迹文件可读。', '', '未到达：r13_fresh_0032__g19（1883）、r13_fresh_0032__a_plus（2060）、r13_fresh_0038__g19（1863）。完整原候选ID、iteration、snapshot和源行见RAW_FIDELITY。即使材料逐位相等，也没有核验辅助J/Q、PR历史或corrector端点相等。','','## 二、实际readout的物理变化与反事实梯度','',f'475个保留readout上，实际原物理DeltaPhi严格负/零/正为{actual["negative"]}/{actual["zero"]}/{actual["positive"]}。其中反事实onehop loss下降但物理Phi上升有{cf_down_up}行；ordinary反事实loss下降但物理Phi上升有{cfOrd_down_up}行。这些是保留点逐行描述，不是全部native迭代率，也不是独立样本统计。', '', '`actual_delta_Phi=Phi(m_after)-Phi(m_before)`来自既有同原物理评估；`physical_first_order_step=gPhi_before·(m_after−m_before)`只是局部一阶量。`counterfactual_ordinary/onehop_loss_change`与其梯度在同一个真实postcurrent J/Q处另算，G19轨迹上的ordinary/onehop梯度也不是G19原生材料导数。不能把其cosine当作“实际使用了这一梯度”。', '', 'ACTUAL_READOUT_DETAILS逐行保留方向量、物理评估ID与counterfactual值；READOUT_EXTREME_EXAMPLES给每case/method最大上升及最大下降（仅按公开DeltaPhi，非GT）。','','## 三、原生PR方向与readout后的current force','', '`CURRENT_FORCE_DETAILS`的方向是已保存native current方向。before/after force都在同一postcurrent J/Q/A上用beta_current计算，tau=gA_after−gA_before。saved native gradient来自pre-current位置，三者不能混称同一时刻。负gA·direction表示该已存方向对该重算force的一阶下降；并不表示已经执行下一步PR，也不证明下次line search结果。', '', 'CASE_METHOD_SUMMARY给严格符号矩阵及normalized dot、gradient rotation、tau/||gA_before||的范围。严格符号计数不设容差；before dot常极接近零，其浮点正负不能直接解释为稳定的下降/上升切换，必须连同具体量级读取。DR没有CIE readout/PR transport行，保持结构性缺失。','','## 四、u/v及端口差','', '`u=J−mQ`，`v=Q−Einc−GJ`，forcing为`u+m v`。current discrepancy为`||J−Jphys(m)||`；port discrepancy为`||S J−Yphys(m)||/delta`；auxiliary residual为`||S J−Y||/delta`，physical residual为`||Yphys(m)−Y||/delta`。`directional_amplification=||J−Jphys||/||u+m v||`是此状态方向上的观察比值，不是resolvent算子范数或普遍上界。', '', 'FIELD_BLOCK_SUMMARY按case/method/field block分列min/median/max与最大值源行。CIE fixed_current_before_readout和fixed_current_after_readout使用同J/Q，不同m；DR分别保留shadow/reflected-A/pre-iterator语义，不把不同field block池化。','','## 五、native auxiliary history、原选择与corrector','', 'DENSE_HISTORY_SUMMARY总结104754条原native scalar/audit记录：current_step_objective_change是原辅助current步目标变化，并不是DeltaPhi；state_term_current与state_term_material使用的beta可不同。CIE逐迭代只有wall，缺失CPU保持null；不能用wall画成实际CPU。DR sparse shadow gap只对已付审计记录成立。', '', 'SELECTION_CORRECTOR_PUBLIC取自此前comparison48/RUN_INDEX.json中原冻结initial/raw/corrected/primary/secondary的候选ID和公开D/R/Phi/ball/CPU，不附新GT分数。raw→corrected、raw→primary等差值只在同原运行内部比较。Replay的paid audit相位数另列，不把计时回放的corrector端点冒充原corrector。所有40任务保留，不按Phi改善筛任务。','','## 来源、定位与范围','', 'JSON源行ID为TABLE[零起索引]，对应CSV物理行号index+2（首行为header）。表内source_row_ids或压缩连续range可回到精确输入行；SNAPSHOT字段与physical_evaluation_id进一步定位。`SOURCE_SCHEMA`列完整实际schema。', '', '定义源码：supplement_aplus48/src/analyze_native_states.py:171–193（field defects）、219–264（current force）、546–580（readout/counterfactual）。详细合同见src/README_NATIVE_STATES.md。以上没有显著性检验、成像成功判定、新算法优越性或因果解释；科学判断留给父线程。','']
# Per-case compact observed table; counts refer sparse saved snapshots only.
lines += ['|case|method|readout下降/不变/上升|raw target状态|source coverage|','|---|---|---|---|---|']
for o,raw in zip(overview,rawrows):
 readout_counts=f'{o["actual_physical_DeltaPhi_negative"]}/{o["actual_physical_DeltaPhi_zero"]}/{o["actual_physical_DeltaPhi_positive"]}' if o['readout_rows'] else'不适用DR'
 lines.append(f'|{o["case_id"]}|{o["method"]}|{readout_counts}|{json.dumps(raw["target_status_counts"],ensure_ascii=False)}|{o["source_coverage_id"]}|')
lines += ['', '## 数值定位示例（按实际DeltaPhi极值，非GT挑选）', '', '|method|规则|case / snapshot|DeltaPhi|物理一阶项|ordinary反事实Delta loss|onehop反事实Delta loss|源行|', '|---|---|---|---:|---:|---:|---:|---|']
for method in ('cie','dual_cie','g19','a_plus'):
 rr=[r for r in trans if r['method']==method and num(r.get('actual_delta_Phi'))]
 for label,fn in [('最小',min),('最大',max)]:
  if not rr:continue
  r=fn(rr,key=lambda x:x['actual_delta_Phi'])
  lines.append('|'+method+'|'+label+'|'+r['case_id']+' / '+str(r['snapshot_id'])+'|'+f"{r['actual_delta_Phi']:.8g}"+'|'+f"{r['physical_first_order_step']:.8g}"+'|'+f"{r['counterfactual_ordinary_loss_change']:.8g}"+'|'+f"{r['counterfactual_onehop_loss_change']:.8g}"+'|'+r['_source_id']+'|')
lines += ['', '这些极值仅定位已保存轨迹片段。原生辅助目标、反事实目标和原Phi各列独立；不能把最大/最小行当作代表性分布。']
(OUT/'EVIDENCE_SUMMARY.md').write_text('\n'.join(lines)+'\n')
(OUT/'SUMMARY_RECEIPT.json').write_text(json.dumps(dict(status='EXISTING_TABLES_DESCRIPTIVE_EXTRACTION_COMPLETE',source_root=str(SRC),input_receipt=str(SRC/'RECEIPT.json'),source_run_index=str(INDEX),tasks=40,objects=8,source_tables={k:len(v) for k,v in TABLES.items()},new_forward=0,new_adjoint=0,new_optimizer=0,new_scoring=0,GT_NPZ_read=False,NPZ_read=False,source_native_analysis_paid_LU_reference=receipt['original_Physical_LU'],CPU_seconds=time.process_time()-START,scope='stdlib JSON/CSV aggregation only; source physical evaluation costs not newly incurred'),indent=2)+'\n')
print(json.dumps(dict(tasks=len(coverage),raw_target_status=dict(Counter(r['status'] for r in targets)),actual_readout_signs=dict(actual),output=str(OUT))))
