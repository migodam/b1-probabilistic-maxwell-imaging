"""Summarize authorized, completed old endpoint CSVs; never open predictions or GT."""
from __future__ import annotations

import argparse
import csv
import json
import os
import time
from pathlib import Path

for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_key] = '1'

START_WALL = time.perf_counter()
START_CPU = time.process_time()


def read_csv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))


def dump_csv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def ratio(a, b):
    return a / b if b > 0 else None


def f(row, key):
    return float(row[key]) if row.get(key) not in ('', None) else None


def fmt(value):
    return 'null' if value is None else f'{value:.6g}'


def table(rows, columns):
    return ['| ' + ' | '.join(label for label, _ in columns) + ' |',
            '| ' + ' | '.join('---' for _ in columns) + ' |'] + [
                '| ' + ' | '.join(fmt(row[key]) if isinstance(row[key], (float, int)) else str(row[key])
                                  for _, key in columns) + ' |' for row in rows]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    source = root / 'audit/old_endpoint_geometry_v1'
    out = (args.out or root / 'analysis/old_geometry_summary').resolve()
    opened = []

    def load_json(path):
        opened.append(str(path))
        return json.loads(path.read_text())

    auth = load_json(source / 'AUTHORIZATION_AT_ACTIONS.json')
    receipt = load_json(source / 'RECEIPT.json')
    if not auth.get('authorized') or receipt.get('completed_endpoints') != 48 or receipt.get('failed_endpoints') != 0:
        raise ValueError('Expected parent-authorized 48 completed old endpoints')
    datasets = {}
    for name in ('ENDPOINTS', 'CONDITIONAL_BANDS', 'BACKGROUND_MATERIAL_BANDS', 'PROPAGATED_CROSS_ENERGIES'):
        path = source / (name + '.csv')
        opened.append(str(path))
        datasets[name] = read_csv(path)
    endpoints, bands, background, cross = (datasets[n] for n in
        ('ENDPOINTS', 'CONDITIONAL_BANDS', 'BACKGROUND_MATERIAL_BANDS', 'PROPAGATED_CROSS_ENERGIES'))
    expected = {(f'r14_new_{i:04d}', method) for i in range(24) for method in ('g19', 'dr')}
    if len(endpoints) != 48 or {(r['case_id'], r['method']) for r in endpoints} != expected:
        raise ValueError('Endpoint inventory differs from fixed 24 x 2')
    if not all(r['status'] == 'COMPLETE' for r in endpoints):
        raise ValueError('Incomplete source endpoint')
    if (len(bands), len(background), len(cross)) != (192, 288, 384):
        raise ValueError('Unexpected source table inventory')
    idx = lambda rows, *keys: {tuple(r[k] for k in keys): r for r in rows}
    bi = idx(bands, 'case_id', 'method', 'band')
    bgi = idx(background, 'case_id', 'method', 'operator', 'band')
    ci = idx(cross, 'case_id', 'method', 'quantity')

    ranks = []
    groups = {}
    for i in range(24):
        item = load_json(root / f'audit/rank_all24/real/EXECUTION_{i:02d}.json')
        if item.get('status') != 'COMPLETE' or not item.get('real_material'):
            raise ValueError('Missing real geometry metadata')
        group = item['operator_group']
        if group not in groups:
            groups[group] = load_json(root / f'audit/rank_all24/real/{group}/METADATA.json')
        meta = groups[group]
        ranks.append(dict(case_id=f'r14_new_{i:04d}', operator_group=group, scalar_domain=meta['scalar_domain'],
                          material_dof=meta['N'], tau=meta['tau'], numerical_rtol=meta['numerical_rank_rtol'],
                          **meta['ranks'], S_sigma_max=meta['S_sigma_max'], G_sigma_max=meta['G_sigma_max'],
                          cross_S_sector_G_Gram_Frobenius=meta['cross_S_sector_G_Gram_Frobenius'],
                          global_projector_commutator_Frobenius=meta['global_projector_commutator_Frobenius']))
    ri = {r['case_id']: r for r in ranks}
    compact = []
    for row in sorted(endpoints, key=lambda r: (r['case_id'], 0 if r['method'] == 'g19' else 1)):
        case, method = row['case_id'], row['method']
        total = f(row, 'material_squared_error')
        result = dict(case_id=case, method=method, operator_group=ri[case]['operator_group'],
                      material_error_energy_denominator=total, material_L2_error=f(row, 'material_L2_error'),
                      physical_Phi=f(row, 'physical_Phi'), physical_relative_residual=f(row, 'physical_relative_residual'))
        for band in ('U11', 'U10', 'U01', 'U00'):
            br = bi[(case, method, band)]
            result[band + '_rank'] = int(br['material_rank'])
            result[band + '_material_error_energy'] = f(br, 'material_squared_error')
            result[band + '_material_error_fraction'] = f(br, 'material_error_fraction')
        result['material_fraction_sum'] = sum(result[b + '_material_error_fraction'] for b in ('U11', 'U10', 'U01', 'U00'))
        for operator in ('S', 'G'):
            for band in ('strong', 'weak_observable', 'numerical_null'):
                br = bgi[(case, method, operator + '_real', band)]
                result[f'background_{operator}_{band}_rank'] = int(br['rank'])
                result[f'background_{operator}_{band}_error_energy'] = f(br, 'material_squared_error')
                result[f'background_{operator}_{band}_error_fraction'] = f(br, 'fraction')
            result[f'background_{operator}_weak_including_null_error_energy'] = sum(
                result[f'background_{operator}_{b}_error_energy'] for b in ('weak_observable', 'numerical_null'))
            result[f'background_{operator}_weak_including_null_error_fraction'] = ratio(
                result[f'background_{operator}_weak_including_null_error_energy'], total)
        result['chord_total_port_energy_denominator'] = f(row, 'true_port_error_energy')
        result['chord_total_whitened_port_energy_denominator'] = f(row, 'true_whitened_port_error_energy')
        for key in ('exact_real_Snull_background_S_action_energy', 'exact_real_Snull_chord_port_energy',
                    'exact_real_Snull_chord_whitened_port_energy', 'exact_real_Snull_local_port_energy',
                    'exact_real_Snull_local_whitened_port_energy', 'exact_real_Snull_local_gain_per_material_energy',
                    'physical_Phi_gradient_dot_exact_Snull'):
            result[key] = f(row, key)
        result['Snull_isolated_chord_over_total_chord_port_energy'] = ratio(result['exact_real_Snull_chord_port_energy'],
                                                                           result['chord_total_port_energy_denominator'])
        local_total = f(ci[(case, method, 'local_S_port')], 'total_local_action_energy')
        result['local_total_port_action_energy_denominator'] = local_total
        result['Snull_isolated_local_over_total_local_port_energy'] = ratio(result['exact_real_Snull_local_port_energy'], local_total)
        for quantity in ('S_port', 'local_S_port'):
            cr = ci[(case, method, quantity)]
            result[quantity + '_signed_cross_energy'] = f(cr, 'signed_cross_energy')
            result[quantity + '_cross_over_total'] = f(cr, 'cross_over_total')
        compact.append(result)

    enriched_bands = []
    eidx = idx(compact, 'case_id', 'method')
    for row in bands:
        row = dict(row)
        ep = eidx[(row['case_id'], row['method'])]
        chord_denom = ep['chord_total_port_energy_denominator']
        local_denom = ep['local_total_port_action_energy_denominator']
        row.update(material_error_energy_denominator=ep['material_error_energy_denominator'],
                   chord_total_port_energy_denominator=chord_denom,
                   local_total_port_action_energy_denominator=local_denom,
                   chord_isolated_port_energy_over_total=ratio(f(row, 'chord_S_port_energy'), chord_denom),
                   local_isolated_port_energy_over_total=ratio(f(row, 'local_S_port_action_energy'), local_denom),
                   propagated_ratio_scope='Isolated energy / complete transported energy; not additive shares')
        enriched_bands.append(row)
    enriched_cross = []
    for row in cross:
        row = dict(row)
        denominator = f(row, 'total_local_action_energy') if row['quantity'].startswith('local_') else f(row, 'total_actual_chord_energy')
        row.update(complete_energy_denominator=denominator,
                   signed_cross_definition='complete_energy - sum_isolated_component_energies = 2 sum_{b<c} Re inner(v_b,v_c)')
        enriched_cross.append(row)

    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter

    summary = []
    metrics = [f'{b}_material_error_fraction' for b in ('U11', 'U10', 'U01', 'U00')] + [
        'background_S_weak_including_null_error_fraction', 'background_S_numerical_null_error_fraction',
        'background_G_weak_including_null_error_fraction', 'Snull_isolated_chord_over_total_chord_port_energy',
        'Snull_isolated_local_over_total_local_port_energy', 'S_port_cross_over_total', 'local_S_port_cross_over_total']
    for method in ('g19', 'dr'):
        rr = [r for r in compact if r['method'] == method]
        for metric in metrics:
            values = [r[metric] for r in rr if r[metric] is not None]
            if metric.endswith('_material_error_fraction'):
                numerator = sum(r[metric.replace('_fraction', '_energy')] for r in rr)
                denominator = sum(r['material_error_energy_denominator'] for r in rr)
            elif metric.startswith('background_'):
                numerator = sum(r[metric.replace('_fraction', '_energy')] for r in rr)
                denominator = sum(r['material_error_energy_denominator'] for r in rr)
            elif 'chord_over_total' in metric:
                numerator = sum(r['exact_real_Snull_chord_port_energy'] for r in rr)
                denominator = sum(r['chord_total_port_energy_denominator'] for r in rr)
            elif 'local_over_total' in metric:
                numerator = sum(r['exact_real_Snull_local_port_energy'] for r in rr)
                denominator = sum(r['local_total_port_action_energy_denominator'] for r in rr)
            else:
                prefix = 'local_S_port' if metric.startswith('local_') else 'S_port'
                numerator = sum(r[prefix + '_signed_cross_energy'] for r in rr)
                denominator = sum(r['local_total_port_action_energy_denominator' if prefix.startswith('local') else 'chord_total_port_energy_denominator'] for r in rr)
            summary.append(dict(method=method, metric=metric, endpoints=len(rr), unit='dimensionless energy ratio',
                                pooled_numerator=numerator, pooled_denominator=denominator, pooled_ratio=ratio(numerator, denominator),
                                case_median=float(np.median(values)), case_min=min(values), case_max=max(values),
                                positive_count=sum(v > 0 for v in values), negative_count=sum(v < 0 for v in values)))
    focus = [r for r in compact if r['case_id'] in {f'r14_new_{i:04d}' for i in (9, 20, 21, 23)}]
    focus_keys = {(r['case_id'], r['method']) for r in focus}
    focus_bands = [r for r in enriched_bands if (r['case_id'], r['method']) in focus_keys]
    out.mkdir(parents=True, exist_ok=True)
    for filename, rows in [('ENDPOINT_COMPACT.csv', compact), ('METHOD_SUMMARY.csv', summary), ('FOCUS_CASES.csv', focus),
                           ('CONDITIONAL_BANDS_WITH_DENOMINATORS.csv', enriched_bands), ('FOCUS_BANDS.csv', focus_bands),
                           ('BACKGROUND_MATERIAL_BANDS.csv', background), ('SIGNED_CROSS_ENERGIES.csv', enriched_cross),
                           ('REAL_RANK_CASES.csv', ranks)]:
        dump_csv(out / filename, rows)

    colors = {'U11': '#4477AA', 'U10': '#EEAA33', 'U01': '#66CCAA', 'U00': '#BB6677'}
    fig, axes = plt.subplots(2, 1, figsize=(13.5, 6.5), sharex=True, sharey=True)
    for ax, method in zip(axes, ('g19', 'dr')):
        rr = [r for r in compact if r['method'] == method]
        bottom = np.zeros(24)
        for band, color in colors.items():
            yy = np.array([r[band + '_material_error_fraction'] for r in rr])
            ax.bar(np.arange(24), yy, bottom=bottom, color=color, label=band, width=.82)
            bottom += yy
        ax.set_ylabel(method.upper() + '\nmaterial error energy')
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.set_ylim(0, 1.03)
        ax.grid(axis='y', alpha=.2)
    axes[0].legend(ncol=4, loc='upper center', bbox_to_anchor=(.5, 1.24), frameon=False)
    axes[1].set_xticks(np.arange(24), [f'{i:04d}' for i in range(24)], rotation=60)
    axes[1].set_xlabel('r14_new_ case suffix; all 24 old paid primary endpoints')
    fig.suptitle('Real S-conditioned G sectors: orthogonal material-error energy / complete material-error energy', y=.995, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, .96))
    for extension in ('png', 'pdf'):
        fig.savefig(out / ('REAL_FOUR_SECTOR_MATERIAL_ERROR.' + extension), dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
    for method, marker, color in [('g19', 'o', '#4477AA'), ('dr', '^', '#CC6677')]:
        rr = [r for r in compact if r['method'] == method]
        xx = [100*r['background_S_numerical_null_error_fraction'] for r in rr]
        for ax, metric in zip(axes, ('Snull_isolated_chord_over_total_chord_port_energy',
                                    'Snull_isolated_local_over_total_local_port_energy')):
            yy = [r[metric] for r in rr]
            ax.scatter(xx, yy, marker=marker, c=color, label=method.upper(), alpha=.8, s=35)
            for x, y, row in zip(xx, yy, rr):
                if (row['case_id'], method) in focus_keys:
                    ax.annotate(row['case_id'][-4:], (x, y), xytext=(3, 4 if method == 'g19' else -10),
                                textcoords='offset points', fontsize=7, color=color)
            ax.set_yscale('log')
            ax.grid(alpha=.2)
            ax.set_xlabel('Background real S numerical-null\nmaterial-error energy (%)')
    axes[0].set_ylabel('Isolated S-null chord port energy\n/ complete chord port energy')
    axes[1].set_ylabel('Isolated S-null local port-action energy\n/ complete local port-action energy')
    axes[0].legend(frameon=False)
    fig.suptitle('Transported response of background S-null material error; ratios are not additive shares', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, .95))
    for extension in ('png', 'pdf'):
        fig.savefig(out / ('SNULL_TRANSPORTED_PORT_VISIBILITY.' + extension), dpi=180)
    plt.close(fig)

    lines = ['# 旧端点几何：机械离线汇总', '',
        '范围固定为24个旧案例 × G19/DR两种旧primary端点，共48/48 COMPLETE。只读取原授权CSV、rank_all24实材料JSON及原费用回执；未打开任何预测NPZ、GT数组或public物理数组，未增加physics/优化或更改算法。所有结果为已曝光端点的事后描述，科学解释由父线程完成。', '',
        '## 定义、分母与适用范围', '',
        '- h = m_hat − m_coarse_reference；材料误差总能量为 ||h||²。表中份额是实材料正交投影 ||P_b h||²/||h||²，四区可加。这里不是NRMSE，也不是逐坐标支持率。',
        '- real算子以 [Re A; Im A] 实化；先按S的完整右基以tau=.1分强/弱，再在各S区内对G*G压缩矩阵分条件G强/弱。U11/U10为S强内的G强/弱，U01/U00为S弱内的G强/弱；不是全局S强/G强空间的四交集、共同SVD或四个可交换投影。',
        '- 背景S/G strong: sigma/sigma_max >= .1；weak_observable: 1e-12 < sigma/sigma_max < .1；numerical_null: <=1e-12。weak_including_null是后两者之和；背景G weak不同于条件U00。null名称始终按声明数值容差。',
        '- L_hat=I−diag(m_hat)G。有限chord分量为 deltaJ_b=L_hat^-1 diag(P_b h)Q_reference；local分量为 dJ_b=L_hat^-1 diag(P_b h)Q_hat。端口能量是 ||S deltaJ_b||_F² 或 ||S dJ_b||_F²；whitened列再除以原公开标量sigma²。source维按Frobenius求和，未额外除以source数。',
        '- 完整chord分母为 ||S(J_hat−J_reference)||_F²；完整local分母为 ||sum_b S dJ_b||_F²。分别列出原始能量及分母。孤立分量传播后不再正交：signed cross = 完整能量 − 各孤立能量之和 = 2 sum_{b<c} Re<v_b,v_c>。因此传播比率不强迫sum=100%，单分量比率可>1，signed cross可正可负。',
        '- exact_real_Snull是背景实S数值null内的h投影，与上述四区分割独立。S h_null约零不等于材料Jacobian作用 S L_hat^-1 diag(h_null)Q_hat为零，也不等于finite-chord端口响应为零；此处不作材料不可辨识性或恢复保证。',
        '- coarse_reference仅为原反演网格上的Maxwell参考。跨网格观测保留原acquisition mismatch；physical Phi/残差也不构成noise或GT证书。', '',
        '## real rank元数据', '']
    rankrows = []
    for group, meta in groups.items():
        rankrows.append(dict(group=group, cases=sum(r['operator_group'] == group for r in ranks), **meta['ranks'],
                             commutator=meta['global_projector_commutator_Frobenius']))
    lines += table(rankrows, [('组', 'group'), ('案例数', 'cases'), ('U11', 'U11'), ('U10', 'U10'), ('U01', 'U01'),
                            ('U00', 'U00'), ('S rank', 'S_numerical_rank'), ('S null', 'S_numerical_null'),
                            ('G strong', 'global_G_strong'), ('G null', 'G_numerical_null'), ('global commutator F', 'commutator')])
    sm = {(r['method'], r['metric']): r for r in summary}
    lines += ['', '## 全24案例的材料能量汇总', '',
        f'机械描述：G19/DR的U00材料误差份额逐例中位数为 {sm[("g19", "U00_material_error_fraction")]["case_median"]:.8%} / {sm[("dr", "U00_material_error_fraction")]["case_median"]:.8%}；pooled份额为 {sm[("g19", "U00_material_error_fraction")]["pooled_ratio"]:.8%} / {sm[("dr", "U00_material_error_fraction")]["pooled_ratio"]:.8%}。两种算法的24个端点均有负的finite-chord S-port signed cross。指定0021、0023的总材料误差能量远高于0009、0020，故本表不以单一pooled数替代逐例数；不从这些描述推出恢复机制或预测器。', '',
        'pooled_ratio明确是跨24案例的分子能量总和/分母能量总和；case_median是24个案例比率的中位数。大误差案例会影响pooled，因此并列保留，未把48个端点或每个网格点当独立统计检验样本。', '']
    lines += table(summary, [('方法', 'method'), ('指标', 'metric'), ('pooled ratio', 'pooled_ratio'),
                             ('case median', 'case_median'), ('min', 'case_min'), ('max', 'case_max')])
    lines += ['', '## 指定逐例：0009 / 0020 / 0021 / 0023', '',
              '下列份额以0–1表示；完整48行、原始分子与分母在 ENDPOINT_COMPACT.csv。', '']
    lines += table(focus, [('案例', 'case_id'), ('方法', 'method'), ('||h||²', 'material_error_energy_denominator'),
                          ('U11', 'U11_material_error_fraction'), ('U10', 'U10_material_error_fraction'),
                          ('U01', 'U01_material_error_fraction'), ('U00', 'U00_material_error_fraction'),
                          ('S weak+null', 'background_S_weak_including_null_error_fraction'),
                          ('S null', 'background_S_numerical_null_error_fraction'),
                          ('G weak+null', 'background_G_weak_including_null_error_fraction')])
    lines += ['', '同一指定端点的S-null传播与signed cross（原始端口能量，不white）：', '']
    lines += table(focus, [('案例', 'case_id'), ('方法', 'method'), ('||S h_null||²', 'exact_real_Snull_background_S_action_energy'),
                          ('null chord energy', 'exact_real_Snull_chord_port_energy'), ('total chord denominator', 'chord_total_port_energy_denominator'),
                          ('null local energy', 'exact_real_Snull_local_port_energy'), ('total local denominator', 'local_total_port_action_energy_denominator'),
                          ('chord cross', 'S_port_signed_cross_energy'), ('local cross', 'local_S_port_signed_cross_energy')])
    lines += ['', '各band的材料误差、finite-chord/local端口能量（全部192行在 CONDITIONAL_BANDS_WITH_DENOMINATORS.csv；此处固定32行）：', '']
    lines += table([{**r, 'short': r['case_id'][-4:], **{k: f(r, k) for k in
        ('material_error_fraction', 'chord_S_port_energy', 'local_S_port_action_energy')}} for r in focus_bands],
        [('案例', 'short'), ('方法', 'method'), ('band', 'band'), ('材料误差份额', 'material_error_fraction'),
         ('chord port energy', 'chord_S_port_energy'), ('local port energy', 'local_S_port_action_energy')])
    lines += ['', '## 文件与费用', '',
        '- REAL_FOUR_SECTOR_MATERIAL_ERROR.png/pdf：全24×2四区材料误差堆叠，合法可加的材料投影能量。',
        '- SNULL_TRANSPORTED_PORT_VISIBILITY.png/pdf：全24×2背景S-null份额与finite-chord/local端口能量比；对数纵轴，非信息百分比。',
        '- SIGNED_CROSS_ENERGIES.csv保留384条 current/G/port/whitened、chord/local原始能量和signed cross；BACKGROUND_MATERIAL_BANDS.csv保留288条背景谱分区。',
        '- METHOD_SUMMARY.csv保留每项pooled分子/分母与case median/min/max；REAL_RANK_CASES.csv保留24案例分组和rank定义。',
        f'- 原离线端点几何计算费用（不重复计费）：CPU {receipt["CPU_seconds"]:.9f}s，原whole-process receipt {receipt["whole_process_CPU_seconds_at_receipt"]:.9f}s，{receipt["physical_cost"]["LU"]} LU / {receipt["physical_cost"]["actual_complex_RHS"]} complex RHS。原物理counter子项与总RHS重叠，不能相加。',
        '- 本次新增仅CSV/JSON汇总和静态绘图，0 LU / 0 RHS / 0 SVD / 0优化；实际新增whole-process CPU与wall见 ANALYSIS_RECEIPT.json，外层进程结束后补记OS child CPU。无新科学通过阈值或恢复判断。', '']
    (out / 'COMPACT.md').write_text('\n'.join(lines))
    dump = dict(schema='R15_OLD_GEOMETRY_SUMMARY_V1', status='COMPLETE', endpoint_count=48, cases=24,
                methods=['g19', 'dr'], input_files_opened=opened,
                binary_arrays_opened=0, predictions_opened=0, truth_arrays_opened=0,
                physical_LU=0, physical_complex_RHS=0, SVD=0, optimization_calls=0,
                source_endpoint_receipt_CPU_seconds=receipt['CPU_seconds'],
                source_endpoint_physical_cost=receipt['physical_cost'],
                CPU_since_script_main_imports_seconds=time.process_time()-START_CPU,
                whole_process_CPU_at_receipt_seconds=time.process_time(),
                wall_since_script_main_imports_seconds=time.perf_counter()-START_WALL,
                source=__file__, output_files=sorted(p.name for p in out.iterdir()),
                conclusions='Mechanical aggregation only; no scientific validity/novelty or recovery decision')
    (out / 'ANALYSIS_RECEIPT.json').write_text(json.dumps(dump, indent=2) + '\n')
    print(json.dumps(dict(status='COMPLETE', out=str(out), endpoints=48, bands=192, cross=384)))


if __name__ == '__main__':
    main()
