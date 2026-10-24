"""Postclosed DEV FullJTR numerical/cost facts; no physics, mu or GT."""
from pathlib import Path
import sys,json,csv,collections
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,write,pack

def main():
 root=OUT/'results/II/ordinary_fullJ_TR_DEV16_v1';dest=OUT/'AUDIT/ordinary_fullJ_TR_DEV16';source=json.loads((root/'RESULTS.json').read_text());receipt=OUT/'II/jobs/ordinary_fullJ_TR_DEV16_v1/CLOSED.json'
 if dest.exists() or source['status']!='COMPLETE' or json.loads(receipt.read_text())['status']!='COMPLETE':raise ValueError('freshreportandclosedphysicalsourceonly')
 b=Book('II','ordinary_fullJ_TR_DEV16_report_v1',max_wall=600,analysis_only=True);status='COMPLETE';error=None
 try:
  dest.mkdir();sigma=json.loads((OUT/'configs/TARGET.json').read_text())['noise']['sigma'];rows=[];checks=[];allratios=[]
  for case in source['cases']:
   p=root/case['parent']/f'lambda{case["lambda_"]}';ratios=[json.loads(s) for s in (p/'TRIAL_RATIOS.jsonl').read_text().splitlines()] if (p/'TRIAL_RATIOS.jsonl').exists() else []
   with np.load(p/'FINAL.npz',allow_pickle=False) as z:F=z['prediction'].copy();x=z['x'].copy()
   with np.load(OUT/'results/contexts/dev'/f'{case["parent"]}.npz',allow_pickle=False) as z:y=z['observation'].copy()
   with np.load(p/'derivatives'/f'{case["endpoint"]}.npz',allow_pickle=False) as z:J=z['Jcoordinate'].copy();savedg=z['gradient_total'].copy()
   r=pack(F-y)*np.sqrt(2)/sigma;g=J.T@r+case['lambda_']*x;phi=float(.5*r@r);objective=phi+.5*case['lambda_']*float(x@x)
   err=float(np.linalg.norm(g-savedg));phierr=abs(phi-case['Phi']);objerr=abs(objective-case['objective']);expected=8*case['completed_trials']+32*case['accepted_updates']
   if err>1e-10*max(1.,np.linalg.norm(g)) or phierr>1e-9 or objerr>1e-9 or expected!=case['actual_cost']['RHS']:raise AssertionError('postclosegradient/unit/costfailure')
   for trial in ratios:
    predicted=trial['predicted_reduction'];actual=trial['current_objective']-trial['trial_objective'];rho=actual/predicted
    if not np.isclose(rho,trial['ratio'],rtol=1e-12,atol=1e-12) or trial['accepted']!=(np.isfinite(rho) and rho>.1 and actual>0):raise AssertionError('ratioacceptancecontract')
    if trial['accepted'] and actual<=0:raise AssertionError('acceptedobjectiveincrease')
    checks.append(dict(parent=case['parent'],lambda_=case['lambda_'],iteration=trial['iteration'],ratio=rho,accepted=trial['accepted'],step_norm=trial['step_norm'],
     radius=trial['radius_before'],predicted_reduction=predicted,actual_reduction=actual,subproblem_KKT_residual=trial['subproblem']['KKT_residual']))
   row=dict(parent=case['parent'],lambda_=case['lambda_'],anchor_Phi=case['anchor_Phi'],final_Phi=case['Phi'],Phi_change=case['Phi']-case['anchor_Phi'],
    objective_change=case['objective']-case['anchor_objective'],initial_gradient_inf=case['initial_kkt']['gradient_total_inf'],final_gradient_inf=case['kkt']['gradient_total_inf'],
    initial_scaled_GN_gradient_inf=case['initial_kkt']['scaled_GN_gradient_inf'],final_scaled_GN_gradient_inf=case['kkt']['scaled_GN_gradient_inf'],
    initial_GN_decrement=case['initial_kkt']['GN_Newton_decrement'],final_GN_decrement=case['kkt']['GN_Newton_decrement'],GN_rank=case['kkt']['GN_pseudoinverse_rank'],
    final_discarded_gradient_norm=case['kkt']['gradient_in_discarded_subspace_norm'],trials=case['completed_trials'],accepted=case['accepted_updates'],rejected=case['rejected_trials'],
    RHS=case['actual_cost']['RHS'],states=case['actual_cost']['states'],LU=case['actual_cost']['LU'],stop_reason=case['stop_reason'],gradient_replay_error=err,
    Phi_replay_error=phierr,objective_replay_error=objerr,source=str(p/'RESULTS.json'))
   rows.append(row);allratios+=ratios
  with (dest/'PER_PARENT.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
  with (dest/'ALL_TRUST_RATIOS.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(checks[0]));w.writeheader();w.writerows(checks)
  groups=[]
  for lam in (0,1):
   selected=[r for r in rows if r['lambda_']==lam];group=dict(lambda_=lam,parents=16,RHS=sum(r['RHS'] for r in selected),states=sum(r['states'] for r in selected),accepted=sum(r['accepted'] for r in selected),rejected=sum(r['rejected'] for r in selected),
    objective_decreased_parents=sum(r['objective_change']<0 for r in selected),Phi_decreased_parents=sum(r['Phi_change']<0 for r in selected),stops=dict(collections.Counter(r['stop_reason'] for r in selected)))
   for field in ('Phi_change','objective_change','initial_gradient_inf','final_gradient_inf','initial_scaled_GN_gradient_inf','final_scaled_GN_gradient_inf','initial_GN_decrement','final_GN_decrement'):
    v=np.array([r[field] for r in selected]);group[field]=dict(minimum=float(v.min()),median=float(np.median(v)),maximum=float(v.max()),mean=float(v.mean()))
   groups.append(group)
  closed=json.loads(receipt.read_text());result=dict(status='FULLJTR_DEV16_CLOSED_NUMERICAL_REPLAY_PASS',groups=groups,actual_cost=source['actual_cost'],physical_receipt=str(receipt),
   controlled_wall=closed['controlled_wall_seconds'],profile=closed['profile'],fresh_J32_RHS=32*sum(r['accepted'] for r in rows),Full8_RHS=8*sum(r['trials'] for r in rows),
   max_replayed_gradient_error=max(r['gradient_replay_error'] for r in rows),max_subproblem_KKT_residual=max(r['subproblem_KKT_residual'] for r in checks),
   no_extra_physics=True,GT_read=False,mu_used=False,old_320_cases_or_bestquery_changed=False,L1_control_launched=False,
   stationarity_scope='exactcurrentFullJgradientandGNdecrement; nodeclaredbinaryglobalstationaritypass;GNnotfullnonlinearHessian')
  write(dest/'RESULTS.json',result)
  text=['# 新Full-J信赖域普通控制：DEV16封账事实','',
   '预注册同一原full66 probit解码与数据白化，两臂lambda=0的Phi和lambda=1的Phi+.5||x||²。每臂每父对象16新状态/256RHS，所有accepted更新都新付J32，每个Full候选8RHS，初始J复用旧付费context。没有mu、GT、模型/原代表元改动、超参搜索或L1运行。','',
   '|目标lambda|RHS|Full状态|接受更新|拒绝Full试探|Phi变化中位|目标变化中位|初始GN decrement中位|末端GN decrement中位|','|---|---:|---:|---:|---:|---:|---:|---:|---:|']
  for g in groups:text.append(f"|{g['lambda_']}|{g['RHS']}|{g['states']}|{g['accepted']}|{g['rejected']}|{g['Phi_change']['median']:.9g}|{g['objective_change']['median']:.9g}|{g['initial_GN_decrement']['median']:.9g}|{g['final_GN_decrement']['median']:.9g}|")
  text += ['',f"合计{source['actual_cost']['RHS']}/8192RHS，原Full8部分{result['Full8_RHS']}、accepted新J32部分{result['fresh_J32_RHS']}，全部计费。受控{result['controlled_wall']:.3f}秒，峰值RSS{closed['profile']['RSS_peak_bytes']}字节。原cached锚点Full+J冷取得40RHS、96+64上下文优化拥有另披露，不当作免费冷启动。",'',
   '保存每次原Full current/field/prediction、accepted的完整Jx/梯度、实际与预测下降、信赖比、半径和TR子问题KKT。只读重放验证单位、梯度、接受条件和费用全部一致。GN decrement是gᵀ(JᵀJ+lambdaI)^daggerg平方根，曲率为Gauss-Newton而非真实非线性Hessian；截断rank与被舍弃梯度分量也披露，不追认全局最优或后验资格。','',
   '是否“MAP已经驻定”应按逐父对象原/末端精确梯度与decrement量级判断，不用目标微小下降代替数学证书。旧unscaled SLSQP第一步与callback-only返回问题另见ordinary_slsqp_diagnosis；best-query DEV附加候选及原320个L1case都保持不变。此新控制只能作为预GT诊断补充，任何L1强控制需父线程先审阅并另登记，程序不自动运行。']
  (dest/'REPORT_ZH.md').write_text('\n'.join(text)+'\n');print(json.dumps(result),flush=True)
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:b.close(status,error)
if __name__=='__main__':main()
