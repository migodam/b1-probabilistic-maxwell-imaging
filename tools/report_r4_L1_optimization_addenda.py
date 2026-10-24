"""Closed non-GT objective/KKT facts for fixedL1addenda; no scoring/physics."""
from pathlib import Path
import sys,json,csv,collections
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,write,pack
from p3r4.evaluation import checkpoint_identity

def main():
 dest=OUT/'AUDIT/L1_optimization_addenda_v2';tr=OUT/'results/fresh/locked/ordinary_fullJ_TR_L1_64_addendum_v1';best=OUT/'results/fresh/locked/best_query_L1_64_addendum_v1'
 if dest.exists():raise FileExistsError('preservefirstfacts')
 for base in (tr,best):
  seal=json.loads((base/'CHOICES_SEALED_NO_GT.json').read_text())
  if seal['GT_read'] is not False or json.loads(Path(seal['receipt']['path']).read_text())['status']!='COMPLETE':raise ValueError('closednoscorechoicesrequired')
  for v in seal['artifacts']:
   if checkpoint_identity(v['path'])!=v:raise ValueError('choiceartifactchanged')
 b=Book('IV','L1_optimization_addendum_nonGT_report_v2',max_wall=600,analysis_only=True);status='COMPLETE';error=None
 try:
  dest.mkdir();data=json.loads((tr/'RESULTS.json').read_text());selection=json.loads((best/'RESULTS.json').read_text());sigma=json.loads((OUT/'configs/TARGET.json').read_text())['noise']['sigma'];rows=[];ratios=[]
  for r in data['cases']:
   folder=tr/r['parent']/f'lambda{r["lambda_"]}';cost=r['actual_cost']
   if cost['RHS']>256 or cost['states']>16 or cost['RHS']!=8*r['completed_trials']+32*r['accepted_updates']:raise AssertionError('sameDEVcostguard')
   with np.load(folder/'FINAL.npz',allow_pickle=False) as z:F=z['prediction'].copy();x=z['x'].copy()
   with np.load(OUT/'results/contexts/locked'/f'{r["parent"]}.npz',allow_pickle=False) as z:y=z['observation'].copy()
   with np.load(folder/'derivatives'/f'{r["endpoint"]}.npz',allow_pickle=False) as z:J=z['Jcoordinate'].copy();savedg=z['gradient_total'].copy()
   rw=pack(F-y)*np.sqrt(2)/sigma;g=J.T@rw+r['lambda_']*x
   if not np.allclose(g,savedg,rtol=1e-12,atol=1e-10):raise ValueError('FullJgradientreplay')
   for line in (folder/'TRIAL_RATIOS.jsonl').read_text().splitlines() if (folder/'TRIAL_RATIOS.jsonl').exists() else []:
    q=json.loads(line);rho=(q['current_objective']-q['trial_objective'])/q['predicted_reduction']
    if not np.isclose(rho,q['ratio'],rtol=1e-12,atol=1e-12) or q['accepted']!=(rho>.1 and q['actual_reduction']>0):raise ValueError('ratioacceptance replay')
    ratios.append(dict(q,parent=r['parent']))
   rows.append(dict(parent=r['parent'],lambda_=r['lambda_'],RHS=cost['RHS'],states=cost['states'],LU=cost['LU'],accepted=r['accepted_updates'],rejected=r['rejected_trials'],
    anchor_Phi=r['anchor_Phi'],final_Phi=r['Phi'],Phi_change=r['Phi']-r['anchor_Phi'],objective_change=r['objective']-r['anchor_objective'],
    initial_GN_decrement=r['initial_kkt']['GN_Newton_decrement'],final_GN_decrement=r['kkt']['GN_Newton_decrement'],
    initial_gradient_inf=r['initial_kkt']['gradient_total_inf'],final_gradient_inf=r['kkt']['gradient_total_inf'],
    initial_scaled_gradient_inf=r['initial_kkt']['scaled_GN_gradient_inf'],final_scaled_gradient_inf=r['kkt']['scaled_GN_gradient_inf'],stop_reason=r['stop_reason']))
  with (dest/'PER_PARENT_KKT_OBJECTIVE_COST.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
  groups=[]
  for lam in (0,1):
   selected=[r for r in rows if r['lambda_']==lam];d=dict(lambda_=lam,parents=64,RHS=sum(r['RHS'] for r in selected),states=sum(r['states'] for r in selected),accepted=sum(r['accepted'] for r in selected),rejected=sum(r['rejected'] for r in selected),stops=dict(collections.Counter(r['stop_reason'] for r in selected)))
   for key in ('Phi_change','objective_change','initial_GN_decrement','final_GN_decrement','initial_gradient_inf','final_gradient_inf','initial_scaled_gradient_inf','final_scaled_gradient_inf'):
    a=np.array([r[key] for r in selected]);d[key]=dict(minimum=float(a.min()),median=float(np.median(a)),maximum=float(a.max()),mean=float(a.mean()))
   groups.append(d)
  owner=json.loads((OUT/'IV/jobs/ordinary_fullJ_TR_L1_64_addendum_v1/CLOSED.json').read_text());old=json.loads((OUT/'IV/jobs/fresh_locked64_representatives_v1/CLOSED.json').read_text())
  contexts=[json.loads((OUT/'results/contexts/locked'/f'R4L1_{i:05d}.json').read_text()) for i in range(64)]
  facts=dict(status='CLOSED_ADDENDUM_NON_GT_REPLAY_PASS',groups=groups,actual_new_TR_cost=data['actual_cost'],TR_controlled_seconds=owner['controlled_wall_seconds'],TR_profile=owner['profile'],
   new_Full8_RHS=8*sum(r['accepted']+r['rejected'] for r in rows),new_accepted_J32_RHS=32*sum(r['accepted'] for r in rows),
   best_query_original_paid_cost_owner=str(OUT/'IV/jobs/fresh_locked64_representatives_v1/CLOSED.json'),old_representative_search_RHS=80640,
   common_context_actual_RHS_first64=sum(v['new_RHS'] for v in contexts),common_context_owner='alreadypaidlockednative96/probit64;initialFullJcold40perparentreference,notfreeandrepeatednotrebilled',
   bestquery=dict(ordinary_Phi_lower=selection['ordinary_Phi_lower'],mu_objective_lower=selection['mu_objective_lower'],new_RHS=0),
   exploratory_post_method_freeze_pre_GT=True,primary_outputs_unchanged=True,GT_read=False,not_global_stationarity_or_optimization_certificate=True)
  write(dest/'NON_GT_FACTS.json',facts)
  text=['# 第一64个L1：预GT探索性优化补充诊断','',
   '两项choices均先封存并封账；算法字节与DEV as-run相同，所有模型、45risk/26density主输出、原320case/callback flags不变。此处只重放Full-J梯度/ratio、KKT与目标进展，不作GT风险评分。','',
   '|lambda|新RHS|Full状态|接受更新|拒绝试探|Phi变化中位|目标变化中位|GN decrement原/末中位|','|---|---:|---:|---:|---:|---:|---:|---|']
  for r in groups:text.append(f"|{r['lambda_']}|{r['RHS']}|{r['states']}|{r['accepted']}|{r['rejected']}|{r['Phi_change']['median']:.9g}|{r['objective_change']['median']:.9g}|{r['initial_GN_decrement']['median']:.6g}/{r['final_GN_decrement']['median']:.6g}|")
  text += ['', 'lambda1的decrement中位很小，说明多数初始点在GN曲率尺度上已近MAP；逐父对象的最大值/梯度仍必须保留，不能由中位宣告全部全局收敛。lambda0的Proper TR确能下降数据Phi，原unscaled SLSQP+callback-only返回停滞不能当强优化基线。所有TR仍受相同16状态/256RHS上限，停止不是全局最优证书。', '',
   f"实际新TR {data['actual_cost']['RHS']}/32768RHS，其中Full8={facts['new_Full8_RHS']}、acceptedJ32={facts['new_accepted_J32_RHS']}。受控{owner['controlled_wall_seconds']:.3f}秒，CPU{owner['profile']['CPU_seconds']:.3f}秒。cached初始J来自原96+64上下文，64对象实际context拥有{facts['common_context_actual_RHS_first64']}RHS，不作为免费冷启动。best-query复用原80640RHS查询池，新增Full0。",'',
   f"best-query ordinary 64对象中{selection['ordinary_Phi_lower']}个Phi比原返回更低，mu256case中{selection['mu_objective_lower']}个找到更低可行mu目标。它叫BEST_VALIDATED_QUERY_CANDIDATE，不叫优化器接受步或收敛解。原outputs和成本保留；没有mu缩放、额外迭代或参数搜索。",'',
   '旧callback的Python触发条件是major iteration增加，不是下降/可行证书。Full-J TR的精确梯度、GN decrement与trust ratio另行逐对象保存；GN不是完整非线性Hessian，不据此宣称全局最优。任何后续GT附录只能读取封存choices，64个父对象为统计单位，不能把两lambda臂当128个独立对象。补充发生于方法冻结之后、GT之前，仅探索性，不追认confirmatory primary结论。']
  (dest/'DIAGNOSIS_REPORT_ZH.md').write_text('\n'.join(text)+'\n');print('L1nonGTreplayfactscomplete.')
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:b.close(status,error)
if __name__=='__main__':main()
