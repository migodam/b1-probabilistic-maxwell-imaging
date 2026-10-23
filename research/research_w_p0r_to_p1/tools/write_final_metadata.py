"""Post-close, parent-reviewed summary from the existing receipts. No physics.

Raw endpoint/stage records are read-only. This creates current display gates
and human accounting, keeping historical R/old-P0 decisions in their own scope.
"""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT/name).read_text())


def write(name, value):
    (ROOT/name).write_text(json.dumps(value, indent=2, ensure_ascii=False,
                                  allow_nan=False)+'\n')


def main():
    closed = read('CLOSED.json')
    if closed['status'] != 'CAMPAIGN_CLOSED' or closed['scientific_decision'] != 'REDIRECT':
        raise RuntimeError('Requires the existing parent-reviewed REDIRECT close')
    resource = read('RESOURCE_LEDGER.json')
    prefix = 'results/B2/validation/'
    joint = read(prefix+'joint_seed2_narrow_fullj/RESULTS.json')
    cavity = read(prefix+'cavity_seed2_narrow_bestcal/RESULTS.json')
    argmax = {key: read(prefix+name+'/argmax_validation/TASK_ARGMAX_GRID_CHECK.json')
              for key, name in [('joint','joint_seed2_narrow_fullj'),
                                ('cavity','cavity_seed2_narrow_bestcal')]}
    if not all(r['task_argmax_numerical_pass'] for r in argmax.values()):
        raise RuntimeError('Independent task qualification failed; parent review is required')
    cal_pass = bool(cavity['numerical_pass'] and cavity['d_cal'] <= 1
                    and cavity['d_T_sup_24'] >= 3)
    cavity_probe8 = next(r['distance'] for r in read(prefix+'cavity_seed2_narrow_bestcal/ALL_TASKS.json') if r['probe']==8)
    with (ROOT/'results/C/phase_diagram/PHASE_TABLE.csv').open() as stream:
        phase = list(csv.DictReader(stream))
    wide_fixed = float(next(r['d_mean_all'] for r in phase
                            if r['case']=='joint_contrast3' and r['band']=='wide'
                            and float(r['common_wall_loss_scale'])==1.))
    wide_reopt = read('results/C/wide_reoptimized/RESULTS.json')['reoptimized_pair_wide_distance']
    b2 = dict(status='B2_ORACLE_ONLY_CONDITIONAL_FIXED_SUPPORT',
        joint_physical_numerical_pass=bool(joint['numerical_pass'] and joint['joint_physics_threshold']),
        online_discovery=False, calibration_pair_pass=cal_pass,
        d_joint=joint['d_joint_shifted'], d_cal=joint['d_cal'],
        d_T_sup_24=joint['d_T_sup_24'], target_shift=.1, mass_relative_change=0,
        task_argmax_independently_qualified=argmax['joint'],
        cavity=dict(d_cal=cavity['d_cal'], fixed_probe8_distance=cavity_probe8,
                    d_T_sup_24=cavity['d_T_sup_24'],
                    numerical_pass=cavity['numerical_pass'],
                    calibration_task_threshold_pass=cal_pass,
                    task_argmax_independently_qualified=argmax['cavity'],
                    selection='pre-frozen completed coarse calibration only; not task/GT selected',
                    controlled_P1_prior_changed=False),
        original_minimum_feature='UNRESOLVED_SCALAR_SEMANTICS_NOT_RETROACTIVE_PASS',
        original_broad_prior_GO=False, no_general_unknown_wall_online_claim=True,
        P1_prior='configs/W_CONTROLLED_PRIOR.json',
        new_probability_meaning='new normalized prior frozen after discovery, before fresh observations; no original broad-family coverage/mass claim',
        continued_independent_branches=['C_FIXED_WORLD_SENSITIVITY','P1_EXPLICIT_NEW_CONTROLLED_FAMILY'])
    write('P0_B2_GATE_DECISION.json',b2)
    summary = dict(status='CLOSED_REDIRECT', scientific_decision='REDIRECT',
        closed_utc=closed['utc'], original_prior_qualification=False,
        original_prior_reason='minimum_feature=0.1 lacks a complete mathematical definition',
        evidence_identities=['ORACLE_DIAGNOSTIC','ORACLE_CONSTRUCTED_PHYSICAL_WITNESS',
                             'ONLINE_INFERENCE_FAILED_WITHIN_BUDGET',
                             'NEW_FROZEN_CONTROLLED_PRIOR_INDEPENDENT_EVALUATION'],
        R=dict(gate=read('P0_R_GATE_DECISION.json'),
               snapshot_scope='original R seal; its not-yet-run P1 note is historical, not current'),
        B2=b2,
        C=dict(status='COMPLETED_FIXED_PAIR_AND_LOCAL_REOPTIMIZATION',
               narrow_joint=joint['d_joint_shifted'], wide_joint=wide_fixed,
               wide_local_reoptimized=wide_reopt,
               global_uniqueness_claim=False,
               primary_source='results/C/phase_diagram/PHASE_TABLE.csv'),
        P1A=read('P1_A_GATE_DECISION.json'), P1B=read('P1_B_GATE_DECISION.json'),
        resources='RESOURCE_LEDGER.json', outside_controlled_family_mass=None,
        no_new_physics_authorization_from_unused_budget=True)
    write('GATE_SUMMARY.json',summary)
    write('METRICS.json',dict(status='FINAL_REVIEWED_SCOPE',
          exact_source_metrics=dict(B2='P0_B2_GATE_DECISION.json',
            C='results/C/phase_diagram/PHASE_TABLE.csv',
            P1A='results/P1A/reference/RESULTS.json',
            P1B='results/P1B/offline_scores/SUMMARY.json',
            paired_CI='results/P1B/offline_scores/PAIRED_CI.json',
            Green='tables/green_control/RESULTS.json',
            material_distances='tables/public_world_distance/RESULTS.json'),
          main_MSE=summary['P1B']['main_MSE'], resource_total=resource['total'],
          scientific_decision='REDIRECT', original_broad_prior_GO=False,
          probabilities_apply_only_to='W_FIXED_SUPPORT_JOINT_V1',
          outside_family_mass=None))
    lines=['# 最终资源与封账\n',
      '唯一计费依据为追加式PHYSICS_LEDGER.jsonl，逐job RECEIPT只用于核对，不重复相加。失败、拒绝、额外精度、adjoint与准备工作均保留；所有阶段上限和scope没有转移或重置。\n',
      '|阶段|complex RHS / 上限|full FE LU / 上限|受控job CPU秒 / 上限|受控job wall秒|',
      '|---|---:|---:|---:|---:|']
    for stage,item in resource['stages'].items():
        t,c=item['used'],item['cap']
        lines.append(f"|{stage}|{t['rhs']:,} / {c['rhs']:,}|{t['factorizations']:,} / {c['factorizations']:,}|{t['cpu_seconds']:.3f} / {c['cpu_seconds']:,}|{t['controlled_job_wall_seconds']:.3f}|")
    t=resource['total']; cap=sum(x['cap']['rhs'] for x in resource['stages'].values())
    lines += ['',f"合计 **{t['rhs']:,} / {cap:,} complex RHS**，**{t['factorizations']:,} 次full FE材料–频率LU**；受控CPU **{t['cpu_seconds']:.3f}秒（{t['cpu_seconds']/3600:.3f}小时）**，job wall合计 **{t['controlled_job_wall_seconds']:.3f}秒**。剩余 **{cap-t['rhs']:,} RHS**封存，不是新授权。GPU=0，新三维=0。",'',
      '每阶段同时核对CPU和受控job wall，均未超其冻结上限。各阶段依次串行，wall不是人为把并行核时当单进程时间。wrapper CPU包含导入、输出和持久化，是保守的完整受控执行费用；solver内层时钟嵌套，不能再加一遍。', '',
      f"最大完成wrapper峰值RSS **{resource['max_completed_wrapper_peak_RSS_bytes']:,} bytes**；独立采样峰值 **{resource['max_sampled_wrapper_RSS_bytes']:,} bytes**，与冻结6GB过程上限分列。", '',
      '## 求解与复用的费用语义\n',
      'full FE LU是材料–频率因子分解，不是独立候选世界数。Green控制另有80次reduced dense-core LU（1228/16/32维），CPU已计入P1B；没有用warm Full-LU=0隐藏小矩阵求解。其1760 RHS/21次full FE LU、setup、核验、残差和instrumentation全部属于P1B总额。详见P1_B_GREEN_COMPUTATION_CONTROL.md。', '',
      '同一bank的历史精确节点复用另列来源：P1B实际650节点中10个坐标复用已支付P1A端点，200 RHS为历史setup，未重复记为新增；所有新noise parent、连续点控制760 RHS和独立验证都实际付费。P1A/P1B共同teacher、reference、bank不能当免费部署资产。', '',
      '## 失败与后处理\n']
    for r in resource['failed_jobs']:
        lines.append(f"- `{r['stage']}/{r['label']}`：CPU {r['CPU']:.6f}秒，wall {r['wall']:.6f}秒；失败费用和原目录保留。")
    for r in resource.get('pre_exec_failures',[]):
        lines.append(f"- `{r['stage']}/{r['label']}`：启动前失败，0 RHS/0 LU；{r['CPU']:.6f}秒为已保守计入CPU总额的launcher elapsed归因，并非实测子进程CPU。wall未独立测量，不写成0。原receipt及ledger均保留。")
    lines += ['', '纯数组独立复算、图库、结构审计和打包在AUDIT/tables中分别保留profile，不新增Maxwell RHS/LU，不冒充从零部署免费。人工代码审阅、文字编写及工具显示没有完整CPU采样，未声明为0。发布与网络用时在发布receipt另列，不挪入物理余额。', '',
      '本地保留所有现存原始状态。Release包含完整响应、材料、粒子、配置、源码和第一joint见证cal/joint×两世界×五频的明确20份完整复场子集；其余full fields保留本地逐文件库存，不声称全部上传。', '',
      '封账凭据：[CLOSED.json](CLOSED.json)；逐阶段机器表：[RESOURCE_LEDGER.json](RESOURCE_LEDGER.json)；追加式原账本：[PHYSICS_LEDGER.jsonl](PHYSICS_LEDGER.jsonl)。任何新物理查询均需独立授权。']
    (ROOT/'RESOURCE_SUMMARY.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(status='POST_CLOSE_METADATA_WRITTEN',new_physics=0,
                         total=resource['total']),ensure_ascii=False))


if __name__ == '__main__':
    main()
