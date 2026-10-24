#!/usr/bin/env python3
"""Separately authorized isolated GT append for already-sealed fixedL1choices.

Only first64 diagnostic references/addenda; never duplicate primary45/26 score
registries or select/refit/rescale any output after truth.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,csv,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,write
from p3r4.evaluation import checkpoint_identity,material_losses,paired_bootstrap
from p1_material.prior import MaterialPrior

PARENTS=tuple(f'R4L1_{i:05d}' for i in range(64));KEYS=('spline_raw','spline_clipped','r2_raw','r2_clipped')
TR=OUT/'results/fresh/locked/ordinary_fullJ_TR_L1_64_addendum_v1';BEST=OUT/'results/fresh/locked/best_query_L1_64_addendum_v1';OLD=OUT/'results/fresh/locked/representatives'
DEST=OUT/'results/offline_only/L1_optimization_addenda_GT_v1';RELEASE=OUT/'configs/FIRST_HELDOUT_GT_SCORE_RELEASE.json'
def read(p):return json.loads(Path(p).read_text())
def check(v):
 if checkpoint_identity(v['path'])!=v:raise ValueError('sealedsourcechanged:'+v['path'])
def main():
 if DEST.exists():raise FileExistsError('norescoring/reselection')
 release=read(RELEASE)
 if release.get('released') is not True or release.get('all_online_both_cohorts_closed_before_truth') is not True:raise ValueError('explicitglobalGTrelease required')
 for v in release['closed_gates']:
  check(v)
  if v['mtime_ns']>RELEASE.stat().st_mtime_ns:raise ValueError('releaseprecedesprerequisite')
 release_paths={v['path'] for v in release['closed_gates']};source_pids=[];seals=[]
 for base in (TR,BEST):
  p=base/'CHOICES_SEALED_NO_GT.json';seal=read(p)
  if str(p.resolve()) not in release_paths or seal['GT_read'] is not False or seal['pid']==os.getpid():raise ValueError('releaseddifferentPIDchoicesrequired')
  check(seal['receipt']);closed=read(seal['receipt']['path']);start=read(Path(seal['receipt']['path']).parent/'START.json')
  if closed['status']!='COMPLETE' or start['pid']==os.getpid() or start['stage']!='IV':raise ValueError('closedmatchingStageIVsource')
  for v in seal['artifacts']:check(v)
  source_pids.append(start['pid']);seals.append(checkpoint_identity(p))
 oldjob=OUT/'IV/jobs/fresh_locked64_representatives_v1';old_start=read(oldjob/'START.json');old_closed=read(oldjob/'CLOSED.json')
 if old_closed['status']!='COMPLETE' or old_start['pid']==os.getpid():raise ValueError('closedindependentoriginalsource')
 reqpath=OUT/'results/fresh/locked/representatives_prepared/FROZEN_REQUEST.json';req=read(reqpath)
 if req['GT_read'] is not False or len(req['qualified'])!=64 or req['created_time_ns']>int(old_start['start_unix']*1e9):raise ValueError('originalonlinechronology')
 if read(OLD/'FROZEN_REQUEST.json')!=req:raise ValueError('originalrequestidentity')
 for v in req['source_identities']:check(v)
 for item in req['qualified']:check(item['mu_identity'])
 old_results=read(OLD/'RESULTS.json')
 if old_results['GT_read'] is not False or old_results['status']!='COMPLETE':raise ValueError('originalnoGTgate')
 for case in old_results['cases']:check(case['endpoint_identity'])
 actions=['original_anchor','original_ordinary',*[f'fixed_{k}' for k in KEYS],*[f'original_rep_{k}' for k in KEYS],
  'TR_lambda0','TR_lambda1','best_query_ordinary',*[f'best_query_{k}' for k in KEYS]]
 refs=['original_anchor','original_ordinary',*[f'fixed_{k}' for k in KEYS],*[f'original_rep_{k}' for k in KEYS]]
 pairs=[(a,r) for a in ('TR_lambda0','TR_lambda1') for r in refs]+[('best_query_ordinary',r) for r in ('original_anchor','original_ordinary','TR_lambda0','TR_lambda1')]
 pairs += [(f'best_query_{k}',r) for k in KEYS for r in ('original_anchor','original_ordinary',f'fixed_{k}',f'original_rep_{k}','TR_lambda0','TR_lambda1')]
 b=Book('IV','L1_optimization_addenda_isolated_GT_v1',max_wall=900,analysis_only=True);status='COMPLETE';error=None
 try:
  DEST.mkdir(parents=True)
  # All chronology/PID/seal validation and plan publication precede truth.
  write(DEST/'CHECKED_GATE_AND_FIXED_ANALYSIS_PLAN_BEFORE_GT.json',dict(status='ALL_CHECKS_PASS_BEFORE_FIRST_TRUTH_READ',scoring_pid=os.getpid(),source_pids=source_pids+[old_start['pid']],
   release_identity=checkpoint_identity(RELEASE),choices_seals=seals,parents=PARENTS,actions=actions,comparisons=pairs,
   statistical_unit='64pairedindependentparents,not128lambda armsor1088actionrecords',metrics=['full_mse','full_mae','roi_mse','roi_mae'],
   paired_parent_bootstrap=dict(draws=4000,seed=20261025211,CI_level=.95,scope='exploratorydescriptiveCI,notconfirmatoryprimary/multiplicityclaim'),
   no_reselection=True,no_rescale=True,no_new_models_or_physics=True,GT_read_at_plan_write=False,created_unix_ns=time.time_ns()))
  cfg=read(OUT/'configs/TARGET.json');prior=MaterialPrior(cfg,read(ROOT/'research/p1_material_target_20261018/inputs/geometry.json')['points']);roi=np.flatnonzero(prior.mask)
  items={r['parent']:r for r in req['qualified']};losses={a:{k:[] for k in ('full_mse','full_mae','roi_mse','roi_mae')} for a in actions};records=[]
  for parent in PARENTS:
   item=items[parent];arrays={};sources={}
   with np.load(item['anchor_cache'],allow_pickle=False) as z:arrays['original_anchor']=z['chi'].copy()
   sources['original_anchor']=item['anchor_cache']
   with np.load(item['frozen_mu'],allow_pickle=False) as z:
    for k in KEYS:arrays['fixed_'+k]=z[k].copy();sources['fixed_'+k]=item['frozen_mu']+':'+k
   for arm in ('ordinary_likelihood_shared',*KEYS):
    label='original_ordinary' if arm=='ordinary_likelihood_shared' else 'original_rep_'+arm;p=OLD/parent/arm/'FINAL.npz'
    with np.load(p,allow_pickle=False) as z:arrays[label]=z['chi'].copy()
    sources[label]=str(p)
    label='best_query_ordinary' if arm=='ordinary_likelihood_shared' else 'best_query_'+arm;p=BEST/parent/arm/'BEST_QUERY_CANDIDATE.npz'
    with np.load(p,allow_pickle=False) as z:arrays[label]=z['chi'].copy()
    sources[label]=str(p)
   for lam in (0,1):
    p=TR/parent/f'lambda{lam}/FINAL.npz'
    with np.load(p,allow_pickle=False) as z:arrays[f'TR_lambda{lam}']=z['chi'].copy()
    sources[f'TR_lambda{lam}']=str(p)
   if any(a not in arrays or arrays[a].shape!=(576,) or not np.isfinite(arrays[a]).all() for a in actions):raise ValueError('all64sourcesmustbefinite,nodropping')
   # FIRST heldout GT file read in this independently gated scoring process.
   with np.load(OUT/'inputs/offline_only'/f'{parent}.npz',allow_pickle=False) as z:truth=z['chi'].copy()
   if truth.shape!=(576,) or not np.isfinite(truth).all():raise ValueError('declaredGTmissing,nodropping')
   rows=[]
   for action in actions:
    values=material_losses(arrays[action][None],truth[None],roi)
    r=dict(parent=parent,action=action,source=sources[action],**{k:float(values[k][0]) for k in losses[action]},scope='fixed64exploratoryaddendum,choicessealedbeforeGT')
    rows.append(r);records.append(r)
    for key in losses[action]:losses[action][key].append(r[key])
   write(DEST/f'{parent}.json',dict(parent=parent,rows=rows))
  with (DEST/'PER_PARENT.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
  np.savez(DEST/'PER_PARENT_LOSSES.npz',parent=np.array(PARENTS),actions=np.array(actions),**{k:np.array([losses[a][k] for a in actions]).T for k in ('full_mse','full_mae','roi_mse','roi_mae')})
  comparisons=[]
  for candidate,reference in pairs:
   comparisons.append(dict(candidate=candidate,reference=reference,metric='full_mse',scope='exploratorydescriptive,notprimaryconfirmatory',
    **paired_bootstrap(losses[candidate]['full_mse'],losses[reference]['full_mse'],draws=4000,seed=20261025211,level=.95)))
  means=[dict(action=a,**{k:float(np.mean(v)) for k,v in losses[a].items()}) for a in actions]
  contexts=[read(OUT/'results/contexts/locked'/f'{p}.json') for p in PARENTS];tr_result=read(TR/'RESULTS.json')
  cost=dict(new_scoring_RHS=0,new_scoring_states=0,new_scoring_LU=0,new_scoring_training=0,
   original320_query_pool_RHS=old_results['actual_cost']['RHS'],original320_cost_owner=str(oldjob/'CLOSED.json'),
   additional_TR_RHS=tr_result['actual_cost']['RHS'],additional_TR_cost_owner=str(OUT/'IV/jobs/ordinary_fullJ_TR_L1_64_addendum_v1/CLOSED.json'),
   best_query_incremental_physics_RHS=0,best_query_search_cost_owner='original320paidquerypool,notfreecold',
   shared_context_actual_RHS_first64=sum(c['new_RHS'] for c in contexts),cached_initial_FullJ_cold_reference_RHS_perparent=40,
   setup_scope='originalcontextnative96/probit64acquiredonceperparent,sharedbyallactions;actualcostnotrebilled',
   full_preGT_risk_physics_queue_RHS=117504,rawdiagnostics_scope='L1+L2rawFullalreadyownedseparately;notrequirednewqueriesforscoring')
  result=dict(status='ISOLATED_GT_APPEND_COMPLETE_FIXED_CHOICES',means=means,comparisons=comparisons,cost_owners=cost,statistical_unit='64parents,not128arms',
   scoring_pid=os.getpid(),source_pids=source_pids+[old_start['pid']],GT_read=True,choices_changed=False,new_model_or_algorithm_selection=False,
   primary45risk26density_master_scores_duplicated=False,post_method_freeze_pre_GT_choice_exploratory_addendum=True,confirmatory_primary_claim=False)
  write(DEST/'RESULTS.json',result)
  text=['# 第一64个L1：封存优化补充的隔离GT附录','',
   '父释放在全部risk/density两cohort、原physical/PPC和两项choices封存之后。此不同PID先核对chronology/全部封存身份，写固定比较计划，才第一次读取本进程heldout真值；未重新选择、缩放、迭代、拟合或查询Full。统计单位为相同64个父对象，两个lambda臂不是128个独立样本。这里仅固定64子集补充比较，主45risk/26density评分由主队列拥有。','',
   '|固定动作顺序|full MSE均值|full MAE均值|ROI MSE均值|','|---|---:|---:|---:|']
  for r in means:text.append(f"|{r['action']}|{r['full_mse']:.9g}|{r['full_mae']:.9g}|{r['roi_mse']:.9g}|")
  text += ['', 'RESULTS.json给48个固定对照的按父对象配对bootstrap；95%区间仅探索性描述，不是主方法confirmatory或多重检验控制声明。全部64保留，无根据GT过滤对象或选择最优lambda。原callback-only代表元、固定Spline_s1/R2 raw/clipped和原锚点均原样作为参照。','',
   f"成本拥有：原320查询池{cost['original320_query_pool_RHS']}RHS，新Full-J TR {cost['additional_TR_RHS']}RHS；best-query新增Full0，原source成本不消失。公共64context实际拥有{cost['shared_context_actual_RHS_first64']}RHS，cached初始Full+J冷参考40RHS/parent，均不重计或当作免费冷启动。原raw诊断/PPC费用另保留原receipt。本评分新增Maxwell/LU/RHS/训练为0。",'',
   'callback major-iteration事件不等于已通过线搜索的下降步；原预算停滞存在选择语义与未缩放方向问题。Full-J TR降低Phi不等于降低GT风险，lambda1小GN decrement也不是全局最优、MMSE或后验均值证书。best-query角色是已经验证的查询候选，不是优化收敛。两项均为方法冻结后/GT前登记的探索性强控制，不改写primary结果；最终科学解释由主线程保留。']
  (OUT/'AUDIT/L1_OPTIMIZATION_ADDENDA_GT_APPEND_REPORT_ZH.md').write_text('\n'.join(text)+'\n');print(json.dumps(dict(status=result['status'],means=means,cost=cost)),flush=True)
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:b.close(status,error)
if __name__=='__main__':main()
