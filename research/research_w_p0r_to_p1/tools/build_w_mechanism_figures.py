#!/usr/bin/env python3
"""Plot ONLY five explicitly named, already completed B2 array/JSON artifacts.

No wr imports, new physical endpoints, observations, truth, queue, or accounting.
Run with the project's existing NumPy/Matplotlib Python; paths derive from this file.
"""
from pathlib import Path
import csv
import json
import math
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np

def _stage_root():
    # Supports the repository layout and the adjacent stage/history Release layout.
    name='research_w_p0r_to_p1'
    for parent in Path(__file__).resolve().parents:
        for candidate in (parent,parent/'research'/name):
            if candidate.name==name and (candidate/'configs/AUTHORIZATION.json').is_file():
                return candidate
    raise RuntimeError('Research W evidence root was not found')

ROOT = _stage_root()
FIG = ROOT / 'figures/mechanism'
TAB = ROOT / 'tables/mechanism'
INPUTS = {
    'compensation_arrays': 'results/B2/compensation_mechanism/FINITE_SECANTS.npz',
    'compensation_result': 'results/B2/compensation_mechanism/RESULTS.json',
    'binary_result': 'results/B2/binary_probability_diagnostic/RESULTS.json',
    'finite_arrays': 'results/B2/task_preserving_2026102202_narrow/finite_response_diagnostic/SECANT_ARRAYS.npz',
    'finite_result': 'results/B2/task_preserving_2026102202_narrow/finite_response_diagnostic/RESULTS.json',
}
COLORS = {'wall': '#b75b27', 'target': '#346caa', 'joint': '#318362',
          'linear': '#346caa', 'full': '#b75b27', 'remainder': '#8860a1'}


def close(a, b, name):
    if not np.allclose(a, b, rtol=2e-11, atol=2e-12):
        raise ValueError(f'Source consistency failure: {name}: {a!r} vs {b!r}')


def realify(a):
    a = np.asarray(a)
    if a.shape != (5, 4, 12) or not np.iscomplexobj(a) or not np.all(np.isfinite(a)):
        raise ValueError('Expected finite complex saved response with shape (5,4,12)')
    return np.concatenate((a.real.ravel(), a.imag.ravel()))


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def table(name, rows):
    path = TAB / f'{name}.csv'
    with path.open('w', newline='', encoding='utf-8-sig') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return str(path.relative_to(ROOT))


def setup():
    names = {f.name for f in font_manager.fontManager.ttflist}
    preferred = ['PingFang SC', 'Heiti SC', 'Arial Unicode MS']
    chosen = next((f for f in preferred if f in names), None)
    # These figures intentionally require an actual available CJK font.
    # Fail clearly instead of silently rendering missing-glyph boxes.
    if chosen is None:
        raise RuntimeError('No Chinese font found; install/select an existing CJK font before plotting')
    matplotlib.rcParams.update({'font.family': [chosen, 'DejaVu Sans'], 'font.size': 11,
        'axes.unicode_minus': False, 'svg.fonttype': 'path', 'axes.spines.top': False,
        'axes.spines.right': False, 'figure.facecolor': '#faf9f5',
        'axes.facecolor': '#faf9f5', 'savefig.facecolor': '#faf9f5'})
    warnings.filterwarnings('error', message='Glyph .* missing from font')
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    return chosen


def finish(fig, name, footnote, sources, csvs, details):
    fig.text(.055, .035, footnote, ha='left', va='bottom', fontsize=9, linespacing=1.65, color='#42464a')
    fig.savefig(FIG / f'{name}.png', dpi=180)
    fig.savefig(FIG / f'{name}.svg')
    plt.close(fig)
    provenance = {'figure': name, 'inputs': {k: INPUTS[k] for k in sources},
        'data_tables': csvs, 'new_rhs': 0, 'new_lu': 0,
        'scope': 'Completed B2 offline diagnostics only; source physics was already paid in B2.',
        'p1_data_truth_plan_access': False, 'details': details}
    write_json(FIG / f'{name}.sources.json', provenance)
    return provenance


def compensation(arrays, result):
    w = realify(arrays['delta_wall'])
    t = realify(arrays['delta_target_induced_increment'])
    j = realify(arrays['delta_joint'])
    close(j, w + t, 'vector identity')
    close(float(arrays['noise_quadrature_std']), .01, 'noise')
    norms = [np.linalg.norm(v) for v in (w, t, j)]
    for v, k in zip(norms, ('d_wall', 'd_target_induced_increment', 'd_joint')):
        close(v, result[k], k)
    e1 = w / norms[0]
    orthogonal = t - np.dot(t, e1) * e1
    e2 = orthogonal / np.linalg.norm(orthogonal)
    basis = np.stack((e1, e2), axis=1)
    close(basis.T @ basis, np.eye(2), 'orthonormal basis')
    projected = np.stack((w, t, j)) @ basis
    residual = max(float(np.linalg.norm(v - basis @ p)) for v, p in zip((w, t, j), projected))
    if residual > 1e-11:
        raise ValueError('The 2D plane fails to preserve saved vectors')
    cosine = float(np.dot(w, t) / (norms[0] * norms[1]))
    close(cosine, result['realified_response_cosine'], 'cosine')
    squared_sum = float(norms[0]**2 + norms[1]**2)
    cross = float(2 * np.dot(w, t))
    close(squared_sum + cross, norms[2]**2, 'quadratic interference')
    close(cross, result['quadratic_interference_term'], 'saved interference')
    rows = [{'vector': n, 'plane_e1': float(p[0]), 'plane_e2': float(p[1]),
        'full_480_real_dim_white_norm': float(d), 'unit': 'dimensionless'}
        for n, p, d in zip(('wall', 'target_induced_increment', 'joint'), projected, norms)]
    csv_path = table('01_compensation_vectors', rows)
    fig, ax = plt.subplots(figsize=(12, 7.1))
    fig.subplots_adjust(left=.085, right=.61, top=.82, bottom=.25)
    fig.suptitle('墙与目标诱导响应相消：两个大分量形成较小 joint 变化', x=.055, ha='left', y=.96, fontsize=19)
    fig.text(.055, .88, '真实保存向量张成的正交平面；480 维实化白化范数在此平面内完整保留', fontsize=11)
    W, T, J = projected
    for start, end, key, label in [(np.zeros(2), W, 'wall', '墙变化'),
            (W, W + T, 'target', '目标诱导增量'), (np.zeros(2), J, 'joint', 'joint 变化')]:
        ax.annotate('', xy=end, xytext=start,
            arrowprops={'arrowstyle': '-|>', 'color': COLORS[key], 'lw': 3, 'mutation_scale': 18})
        ax.plot([], [], color=COLORS[key], lw=3, label=label)
    ax.scatter([0, W[0], J[0]], [0, W[1], J[1]], s=28, color='#34424b', zorder=5)
    ax.text(W[0], W[1]-.09, '墙变化终点', ha='right', va='top', fontsize=10)
    ax.text(J[0]+.08, J[1], 'joint 终点', va='center', fontsize=10)
    ax.text(-.06, -.08, '0', ha='right')
    ax.set(xlim=(-.5, 4.05), ylim=(-.32, 1.1), xlabel='正交轴 e1：沿墙变化（无量纲）',
        ylabel='正交轴 e2（无量纲）')
    ax.set_aspect('equal', adjustable='box')
    ax.grid(alpha=.18)
    ax.legend(loc='upper center', bbox_to_anchor=(.5, 1.48), ncol=3, frameon=False)
    text = (f'墙变化范数          {norms[0]:.5f}\n目标诱导增量范数  {norms[1]:.5f}\njoint 变化范数      {norms[2]:.5f}\n\n'
        f'两分量余弦       {cosine:.5f}\n\n平方范数恒等式\n'
        f'{squared_sum:.5f} − {abs(cross):.5f}\n= {norms[2]**2:.5f}')
    fig.text(.66, .71, text, va='top', fontsize=13, linespacing=1.75)
    details = {'norms': dict(zip(('wall', 'target_induced_increment', 'joint'), map(float, norms))),
        'realified_cosine': cosine, 'squared_norm_sum': squared_sum, 'interference': cross,
        'projection_reconstruction_error': residual, 'projection_basis': 'e1=w/||w||; e2=normalized(t-<t,e1>e1)',
        'whitening': 'Arrays are already divided by sigma_real=.01; no second whitening.',
        'identity': 'DeltaJoint=DeltaWall+[(Joint1-Cal1)-(Joint0-Cal0)]'}
    return finish(fig, '01_compensation_plane',
        '单位：每实分量噪声标准差 σ=0.01 白化后的无量纲响应；5 频率 × 4 源 × 12 接收器，固定 γ=1。\n'
        '离线已付费端点分解；目标诱导项也随墙改变。此图不代表在线墙场减除、后验概率或新物理验证。',
        ['compensation_arrays', 'compensation_result'], [csv_path], details)


def binary(result, compensation_result):
    d = float(result['mean_distance_white'])
    close(d, compensation_result['d_joint'], 'shared fixed pair distance')
    error = .5 * math.erfc(d / (2 * math.sqrt(2)))
    close(error, result['optimal_classification_error'], 'Gaussian classifier')
    sep = float(result['target_position_separation'])
    prior = sep**2 / 4
    hard = sep**2 * error
    close(prior, result['zero_observation_prior_mean_squared_risk'], 'prior risk')
    close(hard, result['hard_endpoint_MAP_squared_risk'], 'MAP risk')
    orders = result['posterior_mean_squared_risk_by_GH_order']
    recomputed = {}
    for order, reported in orders.items():
        z, weights = np.polynomial.hermite.hermgauss(int(order))
        # Under fixed class 0: log posterior odds = d*Z-d²/2, Z~N(0,1).
        # Symmetry makes E[p(1-p)|class0] equal its two-class mixture expectation.
        logit = d * np.sqrt(2) * z - d*d/2
        posterior_p = np.exp(-np.logaddexp(0, -logit))
        risk = float(sep**2 * np.dot(weights, posterior_p*(1-posterior_p))/np.sqrt(np.pi))
        close(risk, reported, f'GH {order} MMSE risk')
        recomputed[order] = risk
    mmse = float(orders[str(max(map(int, orders)))])
    close((prior-mmse)/prior, result['posterior_mean_relative_gain_vs_prior'], 'MMSE gain')
    rows = [{'estimator': k, 'mean_squared_risk': v, 'unit': 'native_position_coordinate_squared',
        'model': 'two_fixed_means_equal_prior_known_common_noise_gain1'}
        for k, v in [('prior_mean_no_observation', prior), ('hard_endpoint_MAP', hard), ('posterior_mean_MMSE', mmse)]]
    csv_path = table('02_binary_risks', rows)
    gh_csv = table('02_binary_quadrature', [{'GH_order': int(k), 'mean_squared_risk': v} for k,v in orders.items()])
    fig, ax = plt.subplots(figsize=(12, 7.1))
    fig.subplots_adjust(left=.12, right=.94, top=.78, bottom=.29)
    fig.suptitle('二固定均值控制：平均决策可降低平方损失', x=.055, ha='left', y=.96, fontsize=19)
    fig.text(.055, .87, f'显式等先验 0.5 / 0.5；位置间距 {sep:g}；均值白化距离 {d:.5f}；最优分类错误率 {100*error:.2f}%', fontsize=12)
    vals = np.array([prior, hard, mmse])
    bars = ax.bar(np.arange(3), vals * 1000, width=.58, color=['#899096', COLORS['wall'], COLORS['joint']])
    for bar, value in zip(bars, vals):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+.09, f'{value:.7f}', ha='center', fontsize=13)
    ax.set_xticks(np.arange(3), ['无观测：先验均值', '观测后：硬端点 MAP', '观测后：后验均值 MMSE'])
    ax.set_ylabel('位置平方误差风险 / 0.001（原生位置坐标²）')
    ax.set_ylim(0, max(vals)*1000*1.25)
    ax.grid(axis='y', alpha=.18)
    ax.set_axisbelow(True)
    fig.text(.12, .22, f'MMSE 相对先验风险下降 {100*(prior-mmse)/prior:.2f}%；相对硬端点 MAP 下降 {100*(hard-mmse)/hard:.2f}%。', fontsize=12)
    return finish(fig, '02_binary_decision_risk',
        '单位：原生位置坐标的平方；条顶标出未缩放风险。已知共同高斯噪声 σ=0.01 / 实分量，固定 γ=1。\n'
        '仅二固定完整 Maxwell 均值的数值积分控制（GH 41/81/161）；不是 P1 连续先验风险或 held-out 成绩。\n'
        '较低平方损失可来自决策平均本身，不能归因于渐进调度或新 Maxwell 算法；没有新增物理端点。',
        ['binary_result', 'compensation_result'], [csv_path, gh_csv],
        {'white_mean_distance': d, 'optimal_classification_error': error, 'risks': dict(zip(
        ('prior_mean', 'hard_MAP', 'MMSE'), map(float, vals))), 'quadrature_by_order': orders,
        'independent_recomputed_GH_risks': recomputed,
        'zero_reported_last_order_change_is_not_a_universal_error_bound': True})


def finite(arrays, result):
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 7.4), sharey=True)
    fig.subplots_adjust(left=.08, right=.97, top=.76, bottom=.31, wspace=.15)
    fig.suptitle('有限步的线性预测与完整响应：校准门槛未通过', x=.055, ha='left', y=.96, fontsize=19)
    fig.text(.055, .87, '固定支撑上的材料 decoder 为仿射；其 secant 余项为零，但 Maxwell 响应仍可非线性', fontsize=12)
    for ax, kind, title, threshold in zip(axes, ('calibration', 'task8'), ('校准响应（要求 ≤ 1）', '任务 8 响应（要求 ≥ 3）'), (1., 3.)):
        for i, endpoint in enumerate(result['endpoints']):
            name = endpoint['endpoint']
            vectors = [realify(arrays[f'{name}_{kind}_{suffix}']) for suffix in ('linear', 'secant', 'remainder')]
            close(vectors[1], vectors[0]+vectors[2], f'{name}/{kind} secant identity')
            norms = [float(np.linalg.norm(v)) for v in vectors]
            cosine = float(np.dot(vectors[0], vectors[1])/(norms[0]*norms[1]))
            close(cosine, endpoint[kind]['linear_secant_realified_cosine'], 'linear/full cosine')
            for n, key in zip(norms, ('linear_prediction_norm','full_secant_norm','Maxwell_finite_remainder_norm')):
                close(n, endpoint[kind][key], f'{name}/{kind}/{key}')
            rows.append({'endpoint': name, 'response': kind, 'linear_white_norm': norms[0],
                'full_secant_white_norm': norms[1], 'remainder_white_norm': norms[2],
                'linear_full_cosine': cosine, 'remainder_fraction_full': norms[2]/norms[1],
                'material_RMS': endpoint['material_RMS'], 'box_max_abs_coordinate': endpoint['box_max_abs_coordinate'],
                'threshold': threshold, 'threshold_relation': '<=' if kind=='calibration' else '>=',
                'full_response_threshold_pass': bool(norms[1]<=threshold if kind=='calibration' else norms[1]>=threshold),
                'norm_unit': 'dimensionless_whitened', 'material_RMS_unit': 'source_material_chart_metric'})
            for j, (n, key, label) in enumerate(zip(norms, ('linear','full','remainder'), ('线性预测','完整有限差','Maxwell 余项'))):
                bar = ax.bar(i+(j-1)*.24, n, width=.215, color=COLORS[key], label=label if i==0 else None)
                ax.text(bar[0].get_x()+bar[0].get_width()/2, n+.08, f'{n:.3f}', ha='center', fontsize=10)
        ax.axhline(threshold, ls='--', lw=1.5, color='#343b40')
        ax.set_xticks([0,1], ['提议步 PROPOSAL', '最终步 FINAL'])
        ax.set_title(title, fontsize=14, pad=13)
        ax.set_ylim(0, 5.3)
        ax.grid(axis='y', alpha=.17)
        ax.set_axisbelow(True)
    axes[0].set_ylabel('白化响应范数（无量纲）')
    axes[0].legend(loc='upper left', ncol=3, bbox_to_anchor=(-.02,1.25), frameon=False, fontsize=11)
    final = next(r for r in rows if r['endpoint']=='FINAL' and r['response']=='calibration')
    fig.text(.08, .23, f'最终校准 {final["full_secant_white_norm"]:.5f} > 1：FAIL；任务距离 ≥ 3 不足以成立 calibration witness。', fontsize=12, color='#9c3629')
    csv_path = table('03_finite_step_norms', rows)
    return finish(fig, '03_finite_step_remainders',
        '单位：σ=0.01 / 实分量白化后的响应范数；全部来自已付费 B2 端点和基线 Jacobian。\n'
        '向量恒等式 full = linear + remainder；三个范数不能相加。材料仿射不意味着物理响应线性。\n'
        '离线搜索后的单个有限步诊断；不能推为普遍障碍、在线方法表现或新增物理验证。',
        ['finite_arrays', 'finite_result'], [csv_path],
        {'scientific_gate': result['scientific_gate'], 'material_decoder': result['material_decoder'],
         'decoder_secant_remainder': result['decoder_secant_remainder'], 'rows': rows})


def main():
    # Explicit allowlist: nothing else in results, configs, queues or truth is read.
    loaded = {}
    metadata = {}
    for name, rel in INPUTS.items():
        path = ROOT / rel
        metadata[name] = {'path': rel, 'bytes': path.stat().st_size}
        if path.suffix == '.npz':
            with np.load(path, allow_pickle=False) as z:
                loaded[name] = {k:z[k].copy() for k in z.files}
        else:
            loaded[name] = json.loads(path.read_text())
    for name in ('compensation_result','binary_result','finite_result'):
        if loaded[name]['new_rhs'] != 0:
            raise ValueError(f'{name} no longer matches zero-new-physics saved diagnostic')
    font = setup()
    sources = [compensation(loaded['compensation_arrays'], loaded['compensation_result']),
        binary(loaded['binary_result'], loaded['compensation_result']),
        finite(loaded['finite_arrays'], loaded['finite_result'])]
    write_json(FIG/'BUILD_RESULTS.json', {'status':'BUILT_SOURCE_ARRAY_IDENTITIES_CHECKED',
        'font':font, 'source_files':metadata, 'figures':sources, 'new_rhs':0, 'new_lu':0,
        'read_scope':'Only the five explicitly allowlisted completed B2 artifacts',
        'uncertainties':['No new physical validation or P1 posterior qualification.',
            'Binary risks concern two fixed equal-prior means, not a continuous prior or held-out observations.',
            'Finite-step failure is one posthoc endpoint diagnostic, not a universal impossibility theorem.',
            'Native spatial-coordinate units are retained; no unsupported conversion to metres.',
            'Original minimum_feature scalar semantics remain unresolved; these plots do not resolve them.']})
    (FIG/'REPORT.md').write_text('''# 已完成 B2 机制图\n\n本次仅从五个明确列出的已有数组/JSON 生成三组 PNG、SVG、逐图来源和 CSV，新增 RHS/LU 均为 0；没有访问 P1 观测、GT 或计划，也没有导入研究物理、费用或队列模块。\n\n1. **01_compensation_plane**：用真实 480 维实化白化向量构造正交二维平面，完整保留墙变化、目标诱导增量及 joint 的范数。墙与目标诱导增量余弦约 −0.97585；joint 范数约 0.79872。第二个箭头从墙终点出发，以准确表示向量和。该项是离线端点分解，不能当成在线墙场减除。\n2. **02_binary_decision_risk**：等先验的两个固定均值控制，位置间距 0.1；先验均值风险 0.0025、硬端点 MAP 约 0.00344815、MMSE 约 0.00215260。是已知噪声下的二点统计计算，不是 P1 连续先验或 held-out 成绩，也不能证明渐进调度有益。\n3. **03_finite_step_remainders**：同时显示提议步与最终步、校准与任务响应的 linear/full/remainder；最终校准距离约 1.92998，未通过 ≤1，任务距离约 3.47192 ≥3 并不能挽救该 witness。decoder 仿射余项为零不意味着 Maxwell 响应线性。\n\n每张图有同名 `.sources.json`，包含确切来源路径、推导合同及对应 `tables/mechanism` CSV；汇总见 `BUILD_RESULTS.json`。NPZ 中响应已白化，脚本没有重复除以 σ；位置风险保留原生坐标平方，不擅自标成米平方。SVG 使用本机字体路径，便于独立展示。\n\n校验包括数组形状/有限性、真实向量和、正交投影重建、范数/余弦与保存 JSON 一致、平方干涉恒等式，以及二固定高斯模型的分类错误和 MAP/先验风险公式。GH 41/81/161 的报告差为零不构成普遍误差界。源物理早已在 B2 付费，本次增量为零不代表整个研究成本为零。\n\n原始 minimum_feature 单标量语义仍未解决；这些图没有解决该问题。三组图仅覆盖已有离线机制事实与有限二点控制，不做新后验资格、一般未知墙或普遍算法性能主张。\n''', encoding='utf-8')
    print(json.dumps({'status':'BUILT', 'font':font, 'figures':len(sources),
        'output':str(FIG), 'tables':str(TAB), 'new_rhs':0, 'new_lu':0}, ensure_ascii=False))


if __name__ == '__main__':
    main()
