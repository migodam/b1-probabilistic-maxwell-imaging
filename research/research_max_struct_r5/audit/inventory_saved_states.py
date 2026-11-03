"""Read-only R5 state inventory. Never imports solvers or opens truth files."""
from pathlib import Path
import json
import numpy as np

REPO = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
R2 = REPO/'research/research_max_struct_r2'
R4 = REPO/'research/research_max_struct_r4'
DEV = 'obs_b20d5a55677d4eca'
CASES = ['obs_42cb233e4f224bb3','obs_63abafa750fc483c','obs_100bd65498644322']
opened = []

def safe(path):
    path = path.resolve()
    assert not any(('truth' in part.lower() or 'private' in part.lower()) for part in path.parts)
    opened.append(str(path))
    return path

def npz_info(path):
    with np.load(safe(path), allow_pickle=False) as f:
        return {'path':str(path.resolve()), 'arrays':{k:{'shape':list(f[k].shape),'dtype':str(f[k].dtype)} for k in f.files}}

def read_json(path):
    return json.loads(safe(path).read_text())

def schema(value):
    if isinstance(value,dict): return {k:schema(v) for k,v in value.items()}
    if isinstance(value,list): return {'type':'list','length':len(value),'first_schema':schema(value[0]) if value else None}
    return type(value).__name__

def json_info(path):
    d=read_json(path)
    info={'path':str(path.resolve()),'schema':schema(d)}
    inner=d.get('result',d)
    for k in ['status','configuration','receiver_rank','source_metrics','current_solver','sourcewise_cg','stage_history','candidate_scores','selected_candidate','seed_loaded_from','selection','history_rows','auxiliary_status','corrector']:
        if k in inner: info[k]=inner[k]
    if 'history' in d:
        h=d['history']
        if h and isinstance(h[0],dict):
            info['history']={'rows':len(h),'iteration_rows':sum('iteration' in x for x in h),'keys':sorted({k for x in h for k in x}),
                'stage_ends':[x for x in h if x.get('event')=='stage_end'],
                'PR_state_saved':[k for k in ['PR_plus','line_alpha','direction_restarted'] if any(k in x for x in h)],
                'vector_state_saved':False}
    return info

runs=[]
for arm,folder in [('global',R2/'results/blind'/DEV/'global'),('onehop',R4/'results/onehop_absorption'),('phase_pi',R4/'results/onehop_phase_pi'),('self',R4/'results/self_absorption')]:
    finalfolder=R4/'results/global_common_corrector' if arm=='global' else folder
    run={'case':DEV,'arm':arm,'kind':'development','stage_folder':str(folder.resolve()),'final_folder':str(finalfolder.resolve()),
         'npz':[npz_info(p) for p in sorted(folder.glob('*.npz'))],
         'json':[json_info(p) for p in [folder/'AUXILIARY.json',finalfolder/'RECEIPT.json',finalfolder/'CORRECTOR.json'] if p.exists()],
         'missing':['per-iteration m/J/Q/A','old_grad','old_direction','material-inner L-BFGS memory','corrector material iterates'],
         'stage_auxiliary_state_saved':True,'iteration_vector_state_saved':False}
    if finalfolder!=folder: run['npz'] += [npz_info(p) for p in sorted(finalfolder.glob('*.npz'))]
    cfg=finalfolder/'FROZEN_CONFIG.json'
    if cfg.exists():run['config']={'path':str(cfg.resolve()),'value':read_json(cfg)}
    runs.append(run)
for case in CASES:
    for arm in ['global','onehop']:
        p=R4/'results/confirmation_remote/predictions'/case/(arm+'.npz')
        j=R4/'results/confirmation_remote/results/confirmation'/(case+'__'+arm+'.json')
        runs.append({'case':case,'arm':arm,'kind':'confirmation','npz':[npz_info(p)],'json':[json_info(j)],
            'stage_auxiliary_state_saved':False,'iteration_vector_state_saved':False,
            'missing':['stage auxiliary J/Q/A','per-iteration m/J/Q/A','old_grad','old_direction','material-inner L-BFGS memory','corrector material iterates']})
    for arm in ['phase_pi','self']:
        runs.append({'case':case,'arm':arm,'kind':'confirmation','status':'NOT_RUN_NO_MATCHING_ARTIFACT',
            'npz':[],'json':[],'stage_auxiliary_state_saved':False,'iteration_vector_state_saved':False})

public=[npz_info(R2/'data_public'/(DEV+'.npz')),npz_info(R4/'data_public/confirmation/operators.npz')]
public += [npz_info(R4/'data_public/confirmation'/(case+'.npz')) for case in CASES]
result={'status':'READ_ONLY_INVENTORY_COMPLETE','model_assignment':'Sol medium','physical_solves':0,'tests_run':0,
    'truth_access':'NONE','source_root':str(REPO),'runs':runs,'public_payloads':public,
    'read_paths':sorted(set(opened)),
    'reconstruction':{'stage_saved':'Use stored J and Q directly; do not replace auxiliary J with a physical solve from m.',
      'from_coefficients':'S=Us diag(sig) Vh; r=count(sig>0.001*sig[0]); Vd=Vh[:r].conj().T; JD=Vd@((Us[:,:r].conj().T@Y)/sig[:r,None]); construct U with original fourier_apcs(shape,mf,Vd); J=JD+U@A; Q=E0+G@J. This requires SVD/basis arithmetic, not a Maxwell material solve; not executed in this inventory.',
      'confirmation_candidate_m':'candidate_m rows 0,1,2,3 are initial host, stage0,stage1,stage2 in validate append order; these do not determine off-shell auxiliary J.',
      'final_j':'final lowercase j/Q are exact Maxwell fields at selected final material; they are not joint auxiliary stage states.'},
    'replay_requirements':{'old_stage_diagnostics':'No replay needed for m/J/Q at saved stage ends; use metadata from matching AUXILIARY history stage_end.',
      'old_boundary_restart':'Next stage can start from saved m and A using original frequency zero-padding; PR is reset at stage boundaries. Exact numerical reproduction remains environment-dependent.',
      'old_interior_trajectory':'Replay only auxiliary stages to requested point; no corrector needed. Save pre-current m/A plus PR old_grad/old_direction if exact interior restart is required.',
      'confirmation_stage_states':'Replay auxiliary only from public host m0 and original frozen config because only stage m survived. Stop at needed stage; do not rerun corrector to obtain auxiliary state.',
      'callback_gap':'Existing iteration_callback exposes post-readout m/J/A; stage_callback exposes m/J/Q/A at stage ends. A wrapper can retain previous callback m (initial host before first call) for m_before and reconstruct Q=E0+G@J at predetermined checkpoints only. No solver edit needed for these diagnostics. old_grad/old_direction remain unavailable for exact interior PR restart.',
      'absent_controls':'phase_pi/self on confirmation cases have no saved run; cannot be called replay. A new control execution requires explicit separate scope and cost.'}}
(OUT/'STATE_INVENTORY.json').write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n')

lines=['# R5 Stage I：只读状态库存','',
 '模型配置：Sol medium。仅读取公开输入、预测、receipt/history 和源码；物理求解 0、测试 0、GT 读取 NONE。原结果未修改。数组仅作键、shape、dtype 库存，没有求解或数值验证。',
 '',f'工作树：`{REPO}`','',
 '## 结论','',
 '- 旧开发 b20d 的 global/onehop/phase_pi/self 各有三个完整阶段端点 m/J/Q/A，均 mf3/6/10、500/500/500 步；PR 标量历史存在，但 PR 向量状态和逐步材料轨迹未保存。',
 '- 三个确认对象的 global/onehop 只有四行 candidate_m、1503 行标量历史、corrector 标量历史和最终精确物理 m/j/Q。confirmation_solver 没有输出 aux.stage_states，无法从阶段材料重建当时的辅助 J。',
 '- 三个确认对象没有 phase_pi/self 历史结果；不能把旧开发控制移植为已完成确认。',
 '- 最终 prediction 的小写 j 是材料对应精确物理解，stage 的大写 J 是优化中的辅助电流。二者不可互换。',
 '', '## 可重建 J 与状态语义','',result['reconstruction']['from_coefficients'],'',
 'fourier_apcs 的 U 使用投影后 Fourier 原系数坐标：F 的角频率枚举是每轴 `[0,…,mf−1] ∪ [size−mf,…,size−1]`，`indexing="ij"`，再用原 SVD 截断舍入级 APCS 零动作，不能随意正交化后复用 A。阶段 A shape 为 (36,12)、(144,12)、(400,12)。本库存不运行这些重建。',
 '', '阶段 NPZ 的 m 是当前更新后的材料，J/Q 是当次 current 更新的辅助状态；m 更新后不会重新求精确物理解。阶段 metadata（β、mf、rank、stage stop）在 AUXILIARY 的 stage_end，不在 NPZ 内。',
 '', '## PR 与最小重放边界','',
 '两个主源码在每阶段初始化 `old_grad=old_direction=None`；阶段间按频率零填充旧 Fourier A。PR+ 每步用当前梯度和上一梯度/方向，history 仅保留 PR_plus、line_alpha、direction_restarted。已有 iteration_callback 是 post-readout m/J/A；wrapper 保留上一次 callback m（首步用 host）即可获得 pre-readout m。Q=E0+G@J 只在预定 checkpoint 计算，避免每步额外 G 动作。stage_callback 已直接提供 m/J/Q/A。以上无需修改数值源码；精确内部 PR 连续重启仍缺 old_grad/old_direction。',
 '']
for k,v in result['replay_requirements'].items():lines.append(f'- **{k}**：{v}')
lines += ['', '## 逐例逐臂与真实路径','']
for r in runs:
    lines += [f'### {r["case"]} / {r["arm"]}','']
    if not r['npz']:
        lines += ['NOT_RUN / 无匹配已保存控制。',''];continue
    for n in r['npz']:
        lines.append(f'- `{n["path"]}`')
        arrays=n['arrays']; main={k:v for k,v in arrays.items() if not k.startswith('history_')}
        lines.append('  - '+ '; '.join(f'`{k}` {v["shape"]} {v["dtype"]}' for k,v in main.items()))
        h={k:v for k,v in arrays.items() if k.startswith('history_')}
        if h: lines.append('  - 标量列 '+', '.join(f'`{k}`' for k in h)+f'；shape {next(iter(h.values()))["shape"]}。完整 dtype 见 JSON。')
    for j in r['json']:
        lines.append(f'- JSON：`{j["path"]}`')
        if 'history' in j:lines.append(f'  - history {j["history"]["rows"]} 行 / {j["history"]["iteration_rows"]} 更新行；stage_end：'+json.dumps(j['history']['stage_ends'],ensure_ascii=False))
        if 'stage_history' in j:lines.append('  - result.stage_history：'+json.dumps(j['stage_history'],ensure_ascii=False))
    lines += ['- 缺失：'+', '.join(r['missing'])+'。','']
lines += ['## JSON 结构、公开载荷和来源','',
 '`STATE_INVENTORY.json` 保存每个已检查 NPZ 的所有键/shape/dtype、相关 JSON 完整递归结构、阶段 metadata、配置与绝对路径。确认 JSON 外层是 event/slot/pid/task_id/status/result/task；科学 receipt 在 result，配置在 result.configuration，阶段事件在 result.stage_history，candidate_scores 与 candidate_m 顺序一致。corrector_history 是每次目标梯度评估的 `[原 Phi, ||grad Phi||∞]` 对，不是物理残差或材料向量；行数可大于 optimizer iteration 数。',
 '', '公开输入：','']
for p in public: lines.append(f'- `{p["path"]}`：'+', '.join(f'{k}{v["shape"]}' for k,v in p['arrays'].items()))
lines += ['', '来源代码：R4 src/cie_onehop.py、src/confirmation_solver.py；R2 src/cie_global.py、experiments/run_blind.py；成熟基线 src/cie_tsom_baseline.py 的 fourier_apcs；R4 onehop/phase/self 的 executed_source_snapshot 与 FROZEN_CONFIG。global_common_corrector 复用 R2 阶段候选并保留历史 incumbent，故其最终选择还需查 RECEIPT.selected_candidate。',
 '', '只读脚本仅依赖 pathlib/json/numpy，不导入项目求解模块；只枚举明确 allowlist 的结果和公开载荷，不递归扫描 sealed_truth/private。执行环境按任务指定为 Gaussian/.venv_nn/bin/python。该库存不构成数值重放、科学确认或新结果发布。','']
(OUT/'STATE_INVENTORY.md').write_text('\n'.join(lines))
print(json.dumps({'status':result['status'],'runs':len(runs),'npz_files':sum(len(r['npz']) for r in runs),'outputs':[str(OUT/'STATE_INVENTORY.md'),str(OUT/'STATE_INVENTORY.json')]}))
