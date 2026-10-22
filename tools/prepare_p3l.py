from pathlib import Path
import sys,json,shutil,platform,subprocess
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,CAPS,write

def main():
    if (OUT/'configs/AUTHORIZATION.json').exists():raise RuntimeError('already frozen; do not replace authorization')
    cfg=json.loads((ROOT/'research/p1_material_target_20261018/configs/ONLINE_TARGET.json').read_text())
    write(OUT/'configs/TARGET.json',cfg)
    write(OUT/'configs/AUTHORIZATION.json',dict(campaign_id='P3L-20261022',explicit_user_authorization=True,authority='User Research P3-L full implementation and experimental campaign, section10/13; independently authorized new stage. Attachment protocol false flags remain historical theory-only.',client_date='2026-10-22',stage_RHS=CAPS,total_complex_RHS=12000000,total_controlled_CPU_seconds=172800,GPU_training_seconds=86400,new_3D=0,new_Mie=0,paid_cloud=False,old_balances_used=False,old_campaigns_read_only=True))
    protocol=dict(campaign_id='P3L-20261022',target_id=cfg['target_id'],master_seed=2026102201,
        target_families=['FULL','LOCAL-C','SLICE','PROFILE'],default_sampling_target='FULL',
        original_reference='dz64 du2; independent standard Gaussian64 and Uniform[0,1]^2',
        likelihood='Phi=norm(F-y)^2/sigma^2=.5 norm(white(F-y))^2',
        material_loss='full-grid mean abs(complex error)^2; MAE=mean abs(complex error)',
        GT_processes='truth generation exits before inference starts; independent offline scoring after method freeze; no imports from offline arrays in online modules',
        baseline=dict(methods=['full_TRF','coarse_TRF'],pilot_parents=['DEV_000','DEV_001','DEV_002','DEV_003'],max_full_evaluations=96,KKT_scaled_tolerance=.001,polish_method='bounded L-BFGS-B with source adjoint',polish_evaluations=32,selection='objective and stationarity, then paid RHS; no GT'),
        development=dict(parents=32,roles={'train':24,'validation':8},deep_references=16,selection='y-only anchor features: final Phi, saturation fraction, property-bound distance, log spectral condition; deterministic farthest-first after standardization, first parent DEV_000; never GT'),
        reference=dict(primary='fixed-preconditioned MALA on exact probit property coordinates',coordinates='x=(z,eta),u=normal_CDF(eta); transformed density=exp(-Phi-.5 norm(x)^2), exact uniform-prior transport with Jacobian, not a Gaussian replacement in q',chains=4,warmup=1024,production=4096,initial_dispersion=[.5,1,1.5,2],warmup_metric_refresh=[256,512,768],production_adaptation=False,primary_chain_checkpoint=256,crosschecks=['canonical blocked pCN+property reflection','elliptical slice+property block','independent tempered SMC if primary fails'],Rhat_max=1.01,material_function_min_ESS=1000,posterior_mean_RMS_MCSE_max=.001,qualification='all diagnostics plus independent-method agreement; never a global mixing proof',per_deep_RHS_planning=330000,reference_cap_may_mark_unqualified=True),
        local_region=dict(id='LOCAL-C-R12-M016',Gaussian_distance_from_y_anchor=12.,material_RMS_from_y_anchor=.16,properties='entire original box',hard_noise_ball=False,escape_mass='unknown without qualified independent FULL reference'),
        weak_dimensions=[2,4,8,16,24],gaussian_frame='orthogonal64, properties explicit; no Gaussian prior precision assigned to uniform coordinates',
        analytic_moments='exact for declared Gaussian field conditional with fixed explicit properties; not exact Maxwell posterior',
        independent_pool=dict(train=512,validation=128,locked=128,truths_shared_per_parent=True,locked_GT='not opened until model/rules frozen and outcomes persisted',locked_expansion_condition='development material-risk or credible UQ signal; otherwise report not run'),
        learning=dict(models=['MLP diagonal Gaussian','MLP full Gaussian','MLP MDN','MLP coupling flow','CNN+MLP Gaussian'],seeds=[2026102270,2026102271,2026102272],max_updates=3000,per_run_training_seconds=1200,teacher_distillation='qualified independent posterior only; otherwise NOT_RUN',joint_NLL='full weak marginal conditional on declared context, not LOCAL-C teacher'),
        geometry=dict(projector_differences=[.01,.005,.0025],weak_rank=8,reference_sampling='full-support kernels; no Frobenius assumption'),
        primary_success=dict(heldout_material_MSE_relative_improvement=.05,paired_parent_bootstrap=2000,CI=.95,reference_quality_required=True,probability_degradation='must report, not assume calibration'),
        decision_order='probability reference -> conditional material imaging -> geometry/cost -> normalized learned density -> bandit/RL if reference/benefit supports it',
        fallback='retain all parents, original prior-mean and best deterministic image; unqualified posterior not treated as teacher or certified UQ',
        historical=dict(B1_baseline=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),P2='read-only regression, never posterior teacher',A30='717d3d26a290248bd2b9aeb55f3d03694b407b84'),
        zero_new_physics_in_preparation=True,no_hash_checks=True)
    write(OUT/'configs/PROTOCOL.json',protocol)
    parents=[]
    for group,n,code in [('development',32,0),('train',512,1),('validation',128,2),('locked',128,3)]:
        for i in range(n):
            parent=f'DEV_{i:03d}' if group=='development' else f'SBI_{group}_{i:04d}'
            parents.append(dict(parent=parent,group=group,role=('train' if i<24 else 'validation') if group=='development' else group,index=i,seed=[2026102201,code,i],new_truth=True))
    write(OUT/'configs/PARENTS.json',parents)
    source=ROOT.parent/'b1-p1-material-target/research/p1_material_target_20261018/inputs/background/weak.npz'
    dest=OUT/'inputs/background/weak.npz';dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,dest)
    write(OUT/'SOURCE_INVENTORY.json',dict(B1_revision=protocol['historical']['B1_baseline'],frozen_A30_revision=protocol['historical']['A30'],background_paid_source=str(source),background_reused=True,all_new_RHS=0,platform=platform.platform(),prior_unchanged=True,attachment_theory_date='2026-10-19',execution_authorization_date='2026-10-22'))
    (OUT/'EXPERIMENT_PROTOCOL.md').write_text('''# P3-L 执行协议（2026-10-22冻结）

正式配置为configs/PROTOCOL.json。原P1 prior、noise、acquisition全部固定；本轮是用户独立授权12M complex RHS、48小时受控CPU、24 GPU-hours的新阶段，理论包的false执行标记不构成本轮授权。

先运行32个新development truth与强化物理反演；从仅依赖y的公开anchor特征中确定16个深reference。四条过分散链，1024 warmup、4096 production，所有拒绝保留。MALA对u=Phi(eta)使用精确uniform-prior transport，canonical q上的pCN/ESS作独立原理核验。未达到Rhat<1.01、材料函数ESS>=1000、mean-RMS MCSE<=.001且跨机制一致的reference标UNQUALIFIED，不能成为teacher。诊断通过也不证明全局混合。

所有成像比较以材料MSE/MAE为主，保留prior mean、充分优化TRF、MAP与全部失败。LOCAL-C采用固定Gaussian距离12及材料RMS .16，无noise-ball硬截断；FULL、LOCAL-C、SLICE、PROFILE严格分开。强条件Gaussian材料矩只对声明Gaussian参考准确。属性box完整保留。

最多512 joint-train、128 validation；若开发风险或可信UQ支持扩大，方法冻结后才生成并评价128 locked新parent。旧locked数据不访问。网络只用合法joint NLL或合格teacher；未合格teacher不蒸馏。几何/RL是计算机制对照，不改变目标。任何更新另存修订记录，不能改原P1来求正结果。

新增求解先由本轮独立ledger预约；失败、拒绝、LU、forward8 RHS、source-adjoint gradient8 RHS、receiver-J32 RHS均计费。每个新endpoint立即保存q、材料、prediction和backward residual；truth/anchor额外保存current/field/LU，链保留所有重复状态及RNG checkpoint。大规模采样的全部LU不持久保存，避免将可重建缓存当原始数据。48h为受控执行CPU wall，不含人工/agent思考及发布等待。

仅在可信概率与独立统计允许时主张>=5% MSE改善且parent paired95%区间排除0。没有改善也形成完整负结果；不把Best-of-K、候选数、短链或residual softmax用作posterior成绩。停止明确失败方法，预算是上限而非耗尽目标。
''')
    print(json.dumps(dict(output=str(OUT),parents_preregistered=len(parents),new_physics=0)))

if __name__=='__main__':main()
