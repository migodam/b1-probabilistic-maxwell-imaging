"""Authorized exposed48-only offline recovery dynamics; never imports a solver.

All incumbents use original legal/eligible public Phi. Truth is an annotation,
not a selector. Large complete paid-candidate table is Release-only.
"""
from __future__ import annotations
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import argparse,csv,gzip,importlib.util,json,math,sys,time
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
import score_dualmetric as public

PAIRS=('00','01','10','11')
REPRESENTATIVES=(2,17,37,38,42,44)
FEATURES=('material_current_potential_work_sum','material_current_potential_work_positive_sum','material_current_potential_work_positive_fraction','delta_F_positive_sum','delta_F_positive_fraction','restarts','native_recorded_steps')

def read(p):return json.loads(Path(p).read_text())
def dump(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def table(out,name,rows):
 dump(out/(name+'.json'),rows)
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with (out/(name+'.csv')).open('x',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
  for r in rows:w.writerow({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()})
def num(v):
 try:
  x=float(v);return x if math.isfinite(x) else None
 except (TypeError,ValueError):return None
def csvrows(p):
 with Path(p).open() as f:
  rows=list(csv.DictReader(f))
 for r in rows:
  for k,v in r.items():
   if v=='':r[k]=None
   elif v in ('True','False'):r[k]=v=='True'
   elif k not in ('pair','branch','case_id','task_id') and num(v) is not None:r[k]=num(v)
 return rows

def native_work(path):
 rows=read(path);out=[];previous=None
 for r in rows:
  if 'global_iteration' not in r or 'data_before_current' not in r:continue
  prior=previous.get('regularizer') if previous else 0. if r['global_iteration']==1 else None
  vals=[r.get('data_term'),r.get('state_term_current'),r.get('regularizer'),r.get('data_after_current'),r.get('state_after_current'),prior]
  work=sum(vals[:3])-sum(vals[3:]) if all(num(x) is not None for x in vals) else None
  out.append(dict(CPU_seconds=r.get('CPU_seconds'),global_iteration=r['global_iteration'],stage=r.get('stage'),work=work))
  previous=r
 return out

def correlation(x,y,np,spearmanr):
 good=[(a,b) for a,b in zip(x,y) if num(a) is not None and num(b) is not None]
 if len(good)<3:return None
 a,b=np.array(good).T
 if np.ptp(a)==0 or np.ptp(b)==0:return None
 return float(spearmanr(a,b).statistic)

def figures(R,O,F,representatives,joined,np,plt):
 colors=dict(zip(PAIRS,['#0072B2','#D55E00','#009E73','#CC79A7']));figpaths=[]
 def save(fig,name):
  source='Source: existing paid intervention snapshots in dynamics_all192_v1; no new forward solves.' if name.startswith('SAME_STATE') else 'Sources: sealed R14 paid HISTORY/pools + dynamics_all192_v1; authorized old exposed48 truth annotations.'
  fig.text(.5,-.005,source,ha='center',fontsize=8)
  for ext in ('png','pdf'):fig.savefig(F/(name+'.'+ext),dpi=165)
  plt.close(fig);figpaths.append('FIGURES/'+name+'.png')
 for page in range(2):
  cases=REPRESENTATIVES[page*3:(page+1)*3];fig,axes=plt.subplots(3,4,figsize=(16,9),squeeze=False)
  for row,c in enumerate(cases):
   for pair in PAIRS:
    local,native=representatives[f'r13_fresh_{c:04d}__pair{pair}'];xs=[r['CPU_seconds'] for r in local]
    for col,key in enumerate(('Phi','NRMSE','residual_over_delta')):axes[row,col].step(xs,[r[key] for r in local],where='post',c=colors[pair],label=pair)
    xx=[];yy=[]
    for r in native:
     if r['work'] is not None and r['CPU_seconds'] is not None:xx.append(r['CPU_seconds']);yy.append(r['work'])
    axes[row,3].plot(xx,yy,c=colors[pair],label=pair)
   for col in range(4):axes[row,col].set_xlabel('Whole-process CPU seconds');axes[row,col].set_xlim(0,120);axes[row,col].grid(alpha=.2)
   axes[row,0].set_ylabel(f'{c:04d}\nPhysical Φ');axes[row,0].set_yscale('log');axes[row,1].set_ylabel('Incumbent NRMSE');axes[row,2].set_ylabel('Residual / δ');axes[row,2].set_yscale('log');axes[row,2].axhline(1,c='gray',ls=':');axes[row,3].set_ylabel('Per-update signed material work');axes[row,3].set_yscale('symlog',linthresh=1e-8);axes[row,3].axhline(0,c='gray',lw=.5)
  axes[0,0].legend(title='Metric pair',ncol=4);fig.suptitle('Fixed representatives: public minΦ incumbent; truth annotation only\nNative signed work stops at proposal end; later physical corrector has no native work',fontsize=11);fig.tight_layout(rect=(0,0,1,.94));save(fig,f'REPRESENTATIVE_DYNAMICS_{page+1}')
 # Compact case-cluster association figure, with all four pair rows preserved.
 fig,axes=plt.subplots(1,3,figsize=(14,4))
 for pair in PAIRS:
  g=[r for r in joined if r['pair']==pair]
  for ax,x,label in zip(axes,('material_current_potential_work_sum','delta_F_positive_sum','restarts'),('Signed material work (sum)','Positive ΔF (sum)','Native PR restarts')):
   ax.scatter([r[x] for r in g],[r['primary_NRMSE'] for r in g],c=colors[pair],label=pair,s=17,alpha=.65);ax.set_xlabel(label);ax.set_ylabel('Primary NRMSE');ax.grid(alpha=.2)
 axes[0].set_xscale('symlog',linthresh=1e-5);axes[1].set_xscale('symlog',linthresh=1e-12);axes[1].set_xticks([0,1e-10,1e-6,1e-2]);axes[0].legend(title='Pair');fig.suptitle('48 exposed case clusters × four dependent metric trajectories; posthoc only');fig.tight_layout();save(fig,'POSTHOC_TRAJECTORY_ASSOCIATIONS')
 snaps=csvrows(R/'audit/dynamics_all192_v1/INTERVENTION_SNAPSHOTS.csv');branches=csvrows(R/'audit/dynamics_all192_v1/INTERVENTION_BRANCHES.csv')
 ids=sorted({r['intervention'] for r in snaps})
 for page in range((len(ids)+3)//4):
  group=ids[page*4:(page+1)*4];fig,axes=plt.subplots(len(group),3,figsize=(14,len(group)*2.5),squeeze=False)
  for ri,ident in enumerate(group):
   for pair in PAIRS:
    for reset in (False,True):
     rows=sorted([r for r in snaps if r['intervention']==ident and r['pair']==pair and r['reset_PR']==reset],key=lambda r:r['relative_update']);ls='--' if reset else '-';label=pair+(' reset' if reset else ' retain')
     for ci,key in enumerate(('delta_physical_Phi_from_branch_start','material_step_angle_to_negative_physical_Phi_gradient_degrees','auxiliary_constitutive_relative_to_J')):
      rr=[r for r in rows if num(r.get(key)) is not None];axes[ri,ci].plot([r['relative_update'] for r in rr],[r[key] for r in rr],color=colors[pair],ls=ls,marker='.',label=label)
   for ax in axes[ri]:ax.set_xlabel('Fixed completed updates');ax.set_xticks([0,8,32]);ax.grid(alpha=.2)
   axes[ri,0].set_ylabel(ident.replace('r13_fresh_','')+'\nΔ physical Φ');axes[ri,0].axhline(0,c='gray',lw=.5);axes[ri,1].set_ylabel('Step angle to −∇Φ\n(degrees)');axes[ri,1].set_xticks([1,8,32]);axes[ri,1].axhline(90,c='gray',ls=':');axes[ri,1].set_ylim(0,180);axes[ri,2].set_ylabel('Constitutive defect / ‖J‖')
  handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,ncol=8,fontsize=8,loc='upper center',bbox_to_anchor=(.5,.93));fig.suptitle('Paid same-full-state interventions: four pairs × retain/reset PR\nAngles use previous-material paid physical gradient; initial/zero-step angles remain missing');fig.tight_layout(rect=(0,0,1,.87));save(fig,f'SAME_STATE_PHYSICAL_ANGLES_{page+1}')
 return figpaths,snaps,branches,ids

def refresh_figures(R,O,start,wall):
 _,tasks,_,account=public.authorize(R,R/'audit/PREDICTION_SEAL_AUTHORIZED.json')
 import numpy as np
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight'})
 inc=read(R/'release/recovery_dynamics_v1/INCUMBENT_CHANGES.json');groups=defaultdict(list)
 for r in inc:
  if int(r['case_id'][-4:]) in REPRESENTATIVES:groups[r['task_id']].append(r)
 representatives={}
 for t in tasks:
  if int(t['case_id'][-4:]) in REPRESENTATIVES:representatives[t['task_id']]=(groups[t['task_id']],native_work(account[t['task_id']]['native_history_file']['path']))
 paths,_,_,_=figures(R,O,O/'FIGURES',representatives,read(O/'TRAJECTORY_PRIMARY_ASSOCIATIONS.json'),np,plt)
 dump(O/'FIGURE_REFRESH_RECEIPT.json',dict(scope='Existing data only: per-update signed work replaces cumulative work that masked signs; full incumbent tables moved to ignored Release side',CPU_seconds=time.process_time()-start,wall_seconds=time.perf_counter()-wall,new_physical_LU=0,new_complex_RHS=0,optimizer_calls=0,figures=paths))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--history-root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--refresh-figures',action='store_true');ap.add_argument('--severe-cases',default='0007,0016,0029,0037,0044');a=ap.parse_args()
 start=time.process_time();wall=time.perf_counter();R=a.root.resolve();H=a.history_root.resolve();O=a.out.resolve()
 if a.refresh_figures:return refresh_figures(R,O,start,wall)
 if O.exists():raise FileExistsError(O)
 seal,tasks,protocol,account=public.authorize(R,R/'audit/PREDICTION_SEAL_AUTHORIZED.json')
 if {t['case_id'] for t in tasks}!={f'r13_fresh_{i:04d}' for i in range(48)}:raise ValueError('Only exposed48 authorized')
 O.mkdir(parents=True);F=O/'FIGURES';F.mkdir();release=R/'release/recovery_dynamics_v1';release.mkdir(parents=True,exist_ok=False)
 # Record the public rule and original sources before loading any truth payload.
 dump(O/'PUBLIC_PLAN_BEFORE_TRUTH.json',dict(authorized_seal=str(R/'audit/PREDICTION_SEAL_AUTHORIZED.json'),tasks=[dict(task_id=t['task_id'],case_id=t['case_id'],pair=t['pair'],history=account[t['task_id']]['history_path'],pool=account[t['task_id']]['candidate_pool_path']) for t in tasks],incumbent='At candidate CPU t<=120, min(Phi,candidate) among original legal=True eligible=True candidates already paid; GT never enters selection',truth_scope='48 old exposed r13_fresh cases only; independent24 forbidden'))
 import numpy as np
 from scipy.stats import spearmanr
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight'})
 spec=importlib.util.spec_from_file_location('reused_old_truth_metrics',H/'src/score_results.py');helper=importlib.util.module_from_spec(spec);sys.path.insert(0,str(H/'src'));spec.loader.exec_module(helper)
 truths={cid:helper.load_truth(helper.truth_spec(H,dict(case_id=cid))) for cid in sorted({t['case_id'] for t in tasks})}
 primary={r['task_id']:r for r in read(R/'analysis/exposed48_v2/ALL_VIEWS.json') if r['view']=='primary'}
 trajectories=csvrows(R/'audit/dynamics_all192_v1/TRAJECTORIES.csv');stages=csvrows(R/'audit/dynamics_all192_v1/NATIVE_STAGE_SUMMARY.csv')
 incumbent_rows=[];summaries=[];issues=[];count=0;invalid=0;representatives={};lasts={}
 columns=['campaign','task_id','case_id','pair','candidate','CPU_seconds','phase','eligible','legal','within120','Phi','D','R','relative_residual','residual_over_delta','NRMSE','IoU','material_legal','incumbent_candidate','incumbent_Phi','incumbent_NRMSE','incumbent_IoU','incumbent_residual_over_delta']
 big=release/'ALL_PAID_CANDIDATES_WITH_OFFLINE_TRUTH.csv.gz'
 with gzip.open(big,'xt',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=columns);writer.writeheader()
  for ti,t in enumerate(tasks):
   ident=t['task_id'];item=account[ident];history=read(item['history_path']);best=None;bestmetrics=None;best_candidate=None;updates=0;uphill=0
   with np.load(item['candidate_pool_path'],allow_pickle=False) as z:pool=z['m']
   if len(pool)!=len(history):raise ValueError('Pool mismatch '+ident)
   local=[]
   # Candidate index is native paid order; sort CPU to make explicitly causal prefix.
   ordered=sorted(history,key=lambda r:(r.get('CPU_seconds',float('inf')),r['candidate']))
   for row in ordered:
    i=row['candidate'];metrics=helper.material_metrics(pool[i],truths[t['case_id']],.02,20.);count+=1
    if not metrics.get('available'):invalid+=1;issues.append(dict(task_id=ident,candidate=i,issue=metrics.get('availability_reason')))
    if public.eligible(row,120.) and (best is None or (row['Phi'],i)<(best['Phi'],best_candidate)):
     if bestmetrics is not None and metrics.get('NRMSE') is not None and bestmetrics.get('NRMSE') is not None:uphill+=metrics['NRMSE']>bestmetrics['NRMSE']
     best=row;best_candidate=i;bestmetrics=metrics;updates+=1
     inc=dict(campaign=t['campaign'],task_id=ident,case_id=t['case_id'],pair=str(t['pair']),candidate=i,CPU_seconds=row['CPU_seconds'],phase=row['phase'],Phi=row['Phi'],D=row.get('D'),R=row.get('R'),residual_over_delta=row.get('residual_over_delta'),NRMSE=metrics.get('NRMSE'),IoU=metrics.get('IoU'))
     incumbent_rows.append(inc);local.append(inc)
    rr={k:row.get(k) for k in columns};rr.update(campaign=t['campaign'],task_id=ident,case_id=t['case_id'],pair=str(t['pair']),within120=row.get('CPU_seconds',float('inf'))<=120,NRMSE=metrics.get('NRMSE'),IoU=metrics.get('IoU'),material_legal=metrics.get('material_legal'),incumbent_candidate=best_candidate,incumbent_Phi=best.get('Phi') if best else None,incumbent_NRMSE=bestmetrics.get('NRMSE') if bestmetrics else None,incumbent_IoU=bestmetrics.get('IoU') if bestmetrics else None,incumbent_residual_over_delta=best.get('residual_over_delta') if best else None);writer.writerow(rr)
   p=primary[ident];mismatch=best_candidate!=p['selection_candidate'];lasts[ident]=local[-1] if local else {}
   if mismatch:issues.append(dict(task_id=ident,issue='PUBLIC_PRIMARY_CANDIDATE_MISMATCH',recomputed=best_candidate,saved=p['selection_candidate']))
   summaries.append(dict(task_id=ident,case_id=t['case_id'],pair=str(t['pair']),paid_candidates=len(history),incumbent_updates=updates,incumbent_NRMSE_increases=uphill,final_candidate=best_candidate,primary_candidate=p['selection_candidate'],primary_selector_consistent=not mismatch,primary_NRMSE=p.get('NRMSE'),primary_IoU=p.get('IoU'),primary_Phi=p.get('Phi'),primary_residual_over_delta=p.get('residual_over_delta')))
   if int(t['case_id'][-4:]) in REPRESENTATIVES:representatives[ident]=(local,native_work(item['native_history_file']['path']))
   if (ti+1)%24==0:print(f'Offline pools {ti+1}/192; paid candidates {count}',flush=True)
 table(release,'INCUMBENT_CHANGES',incumbent_rows);table(O,'TASK_CANDIDATE_SUMMARY',summaries)
 # All trajectory summaries joined to saved primary; no per-step independent samples.
 joined=[]
 for r in trajectories:
  p=primary[r['task_id']];j=dict(r,primary_NRMSE=p.get('NRMSE'),primary_IoU=p.get('IoU'),primary_Phi=p.get('Phi'),primary_residual_over_delta=p.get('residual_over_delta'))
  for s in stages:
   if s['task_id']==r['task_id']:j[f"stage_{int(s['stage'])}_recorded_steps"]=s['native_recorded_steps']
  joined.append(j)
 table(O,'TRAJECTORY_PRIMARY_ASSOCIATIONS',joined)
 stagejoined=[dict(s,primary_NRMSE=primary[s['task_id']].get('NRMSE'),primary_IoU=primary[s['task_id']].get('IoU')) for s in stages];table(O,'STAGE_PRIMARY_ASSOCIATIONS',stagejoined)
 case_rows=[]
 featurelist=list(FEATURES)+sorted({k for r in joined for k in r if k.startswith('stage_') and k.endswith('_recorded_steps')})
 for cid in sorted(truths):
  group=[r for r in joined if r['case_id']==cid];case_rows.append(dict(case_id=cid,**{k:float(np.mean([r[k] for r in group if num(r.get(k)) is not None])) if any(num(r.get(k)) is not None for r in group) else None for k in featurelist+['primary_NRMSE','primary_IoU']}))
 table(O,'CASE_CLUSTER_MEANS',case_rows);assoc=[]
 for unit,rows in [('case_mean_4pairs',case_rows)]+[('pair'+p,[r for r in joined if r['pair']==p]) for p in PAIRS]:
  for x in featurelist:
   for y in ('primary_NRMSE','primary_IoU'):
    rho=correlation([r.get(x) for r in rows],[r.get(y) for r in rows],np,spearmanr)
    loo=[correlation([r.get(x) for r in rows if r['case_id']!=cid],[r.get(y) for r in rows if r['case_id']!=cid],np,spearmanr) for cid in truths];loo=[v for v in loo if v is not None]
    assoc.append(dict(unit=unit,case_clusters=48,feature=x,outcome=y,Spearman=rho,leave_one_case_min=min(loo) if loo else None,leave_one_case_max=max(loo) if loo else None,scope='Exposed posthoc descriptive association; no p-value/predictor/causal recovery claim; 519199 dependent steps are not samples'))
 table(O,'CASE_CLUSTER_DESCRIPTIVE_CORRELATIONS',assoc)
 figpaths,snaps,branches,ids=figures(R,O,F,representatives,joined,np,plt)
 table(O,'SAME_STATE_PHYSICAL_ANGLES',snaps)
 refs=read(H/'analysis/objective_truth/TRUTH_OBJECTIVE_REFERENCES.json');refmap={r['case_id']:r for r in refs};refrows=[]
 for ident,p in primary.items():
  ref=refmap[p['case_id']];refrows.append(dict(case_id=p['case_id'],pair=p['pair'],primary_Phi=p.get('Phi'),primary_NRMSE=p.get('NRMSE'),primary_IoU=p.get('IoU'),GT_reference_Phi=ref['GT_Phi'],primary_Phi_minus_GT_reference=p['Phi']-ref['GT_Phi'],GT_reference_residual_over_delta=ref['GT_residual_over_delta'],matched_discrete_model=ref['matched_discrete_model'],scope='Known feasible better reference only when difference positive; NOT proof of strict local minimum, basin, or stationarity',source=str(H/'analysis/objective_truth/TRUTH_OBJECTIVE_REFERENCES.json')))
 table(O,'ALL_PRIMARY_MINUS_KNOWN_GT_REFERENCE',refrows)
 severe={f'r13_fresh_{int(s):04d}' for s in a.severe_cases.split(',')};table(O,'FIVE_FIXED_SEVERE_CASES_MINUS_REFERENCE',[r for r in refrows if r['case_id'] in severe])
 endpoint_path=R/'audit/endpoint_stationarity_v1/ENDPOINTS.json'
 if endpoint_path.is_file():
  endpoints={r['task_id']:r for r in read(endpoint_path)};erows=[]
  for r in refrows:
   if r['case_id'] not in severe:continue
   e=endpoints.get(r['case_id']+'__pair'+r['pair'])
   if e is None:continue
   rr=dict(r,existing_endpoint_source=str(endpoint_path))
   for k in ('PG_Linf','PG_L2','gradient_L2','directional_derivative','GN_data_curvature','prior_curvature'):rr['existing_endpoint_'+k]=e.get(k)
   erows.append(rr)
  table(O,'FIVE_FIXED_SEVERE_CASES_WITH_EXISTING_ENDPOINTS',erows)
 table(O,'ISSUES',issues)
 angle_valid=sum(num(r.get('material_step_angle_to_negative_physical_Phi_gradient_degrees')) is not None for r in snaps)
 receipt=dict(status='OFFLINE_COMPLETE',finished_utc=datetime.now(timezone.utc).isoformat(),scope='Authorized old exposed48 only; no independent24 truth',tasks=len(tasks),paid_candidates=count,invalid_material_candidates=invalid,incumbent_changes=len(incumbent_rows),primary_selector_mismatches=sum(not r['primary_selector_consistent'] for r in summaries),original_recorded_native_steps=sum(int(r['native_recorded_steps']) for r in trajectories),case_clusters=48,native_stage_summary_rows=len(stages),interventions=len(ids),branches=len(branches),snapshots=len(snaps),available_physical_step_angles=angle_valid,missing_physical_step_angles=len(snaps)-angle_valid,large_candidate_table=str(big),large_table_bytes=big.stat().st_size,CPU_seconds=time.process_time()-start,wall_seconds=time.perf_counter()-wall,new_physical_LU=0,new_complex_RHS=0,optimizer_calls=0,figures=figpaths,issues=len(issues),fixed_severe_cases=sorted(severe),source_code=str(Path(__file__).resolve()))
 dump(O/'RECEIPT.json',receipt)
 lines=['# 已曝光48案例的离线动态附注','',f"192条轨迹、{count:,}个已付费候选；原生逐步记录519,199条，统计单位为48个case。全部候选NRMSE/IoU只附注，不参与prefix选择。公开eligible/legal且CPU≤120的min(Φ,candidate)维护incumbent；与已封存primary不一致数为{receipt['primary_selector_mismatches']}。",'',f"完整候选大表仅在Release侧：`{big}`；Git忽略。完整INCUMBENT_CHANGES也仅在Release侧（多数候选均改进公开Φ）；Git内仅192行任务汇总和关联小表。",'','关联表保留有符号材料功、正ΔF、PR restart和各stage步数。CASE_CLUSTER_MEANS先在每case内平均四pair；相关另按每pair的48case计算，并给逐case删除范围。没有以519,199步作独立样本，没有显著性检验，也没有训练预测器。所有关系是已曝光后的描述。','', '辅助功使用各native stage的当前度量；stage更换后不构成单一跨stage势。公开物理Φ、误差和残差/δ与辅助功分图显示。未返回的超时更新仍由原receipt计费，不能补成已记录步骤。','',f"同state数据：{len(ids)}个完整checkpoint、{len(branches)}枝、{len(snaps)}个0/1/8/32快照；{angle_valid}个已有物理梯度角，{len(snaps)-angle_valid}个缺失（原因原样保留）。实线保留PR，虚线仅初始清空PR；固定更新数而非等CPU。角度是材料step与前一材料的负物理Φ梯度，不是raw defect向量角，也不是完整32步位移角。",'','GT参考来自旧generation-paid逆模型缓存，不新跑forward。ALL_PRIMARY_MINUS_KNOWN_GT_REFERENCE给全192行，五个固定严重失败case另列；正Φ差只说明存在已知更好参考，不证明strict local minimum、stationarity或basin。crossgrid参考的离散化限制原样保留。','', '记录中的direction_restarted均为0；对应相关为null，不代表PR+从未截断或stage换基从未重置memory。阶段0/1/2记录步数均500，仅末stage变化。', '固定代表0002/0017/0037及互补0038/0042、近失0044，全部四pair，没有事后按最佳误差挑图。物理残差/δ=1仅沿用原noise ball，不是新增成功阈值。','',f"本次离线CPU {receipt['CPU_seconds']:.6f}s，wall {receipt['wall_seconds']:.6f}s；新LU/RHS/optimizer均0。原运行与干预成本见源receipt，不能把离线制图成本当recovery成本。本附注不作最终科学或新颖性判断。"]
 (O/'REPORT_ZH.md').write_text('\n'.join(lines)+'\n')
 html='<html><meta charset="utf-8"><title>Exposed48 offline recovery dynamics</title><body><h1>Exposed48 offline recovery dynamics</h1><pre>'+ '\n'.join(lines)+'</pre>'+''.join(f'<img style="width:100%;max-width:1600px" src="{p}">' for p in figpaths)+'</body></html>'
 (O/'report.html').write_text(html);print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
