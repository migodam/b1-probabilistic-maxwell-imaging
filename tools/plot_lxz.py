#!/usr/bin/env python3
"""Plot paid LXZ arrays only; no solver, proposal, training, or inference imports."""
import csv
import html
import json
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'research/research_lxz_20261016'
FIG = OUT / 'figures'


def read(path):
    return json.loads(path.read_text())


def save(fig, stem):
    for ext in ('png', 'pdf'):
        fig.savefig(FIG / f'{stem}.{ext}', dpi=165)
    plt.close(fig)


def e1():
    # Never use the superseded round1/round2 geometry rows.
    rows = [r for r in read(OUT / 'results/round1/METRICS.json')['rows']
            if r['method'] != 'geometry']
    for directory, label in [('round1_geometry_repair', 'geometry'),
                             ('round2_geometry_repair', 'geometry+nuisance')]:
        for r in read(OUT / f'results/{directory}/METRICS.json')['rows']:
            rows.append(dict(r, method=label))
    methods = ['prior', 'local', 'global', 'hybrid', 'tempered',
               'geometry', 'geometry+nuisance']
    records = []
    fig, axs = plt.subplots(1, 2, figsize=(13.8, 4.7), layout='constrained')
    for checkpoint, offset, color in [(1024, -.17, '#799eb7'), (4096, .17, '#194f73')]:
        for k, metric in enumerate(['event_MAE', 'CDF_max_error']):
            means, sds = [], []
            for method in methods:
                rr = [r for r in rows if r['method'] == method and r['checkpoint'] == checkpoint]
                values = np.asarray([r[metric] for r in rr], float)
                if values.ndim == 2:
                    values = values.mean(axis=1)
                means.append(float(values.mean()))
                sds.append(float(values.std(ddof=1)))
                records.append(dict(method=method, checkpoint=checkpoint,
                                    metric=metric, mean=means[-1], SD=sds[-1], n=len(values)))
            axs[k].bar(np.arange(len(methods)) + offset, means, width=.32,
                       yerr=sds, color=color, capsize=2, label=f'{checkpoint:,} proxy calls')
            axs[k].set_xticks(np.arange(len(methods)), methods, rotation=25, ha='right')
            axs[k].grid(axis='y', alpha=.18)
    axs[0].set_ylabel('Mean absolute event-mass error')
    axs[1].set_ylabel('Mean of 4 material-coordinate grid-CDF errors')
    axs[0].legend(frameon=False)
    fig.suptitle('Exposed P2 proxy target: 4 observations × 4 seeds per method\nBars = mean; whiskers = between-run SD, not confidence intervals')
    fig.supxlabel('Geometry uses corrected paired finite differences only; nuisance changes the frozen proposal metric, not the posterior.\nCDF is a saved-grid diagnostic, not a continuum supremum; this is not a full-Maxwell posterior claim.', fontsize=9)
    save(fig, 'e1_proxy_event_cdf')
    (FIG / 'E1_PLOT_SUMMARY.json').write_text(json.dumps(records, indent=2) + '\n')


def candidates():
    source = read(OUT / 'results/E2_completion/CANDIDATES.json')
    methods = ['local', 'geometry', 'continuation', 'prior']
    selections, all_materials = [], []
    for g in source['geometry']:
        anchor = g['anchor']
        with np.load(OUT / f'results/E2_completion/{anchor}.npz', allow_pickle=False) as z:
            known = z['anchor'].reshape(24, 24).copy()
            all_materials.append(known)
            picks = []
            for method in methods:
                eligible = [r for r in source['rows'] if r['anchor'] == anchor
                            and r['method'] == method and r['compatible']]
                # Fixed deterministic online criterion: maximize material distance
                # within the declared residual/prior-feasible set; tie by residual.
                selected = min(eligible, key=lambda r: (-r['material_RMS_to_anchor'], r['residual'])) if eligible else None
                image = z[f"chi_{selected['id']}"].reshape(24, 24).copy() if selected else None
                if image is not None:
                    all_materials.append(image)
                picks.append((method, selected, image))
            selections.append((g, known, picks))
    limits = [(min((1+x.real).min() for x in all_materials), max((1+x.real).max() for x in all_materials)),
              (min(x.imag.min() for x in all_materials), max(x.imag.max() for x in all_materials))]
    portable = {'selection': 'compatible=True; largest material_RMS_to_anchor; residual tie-break; no truth-error selection',
                'scope': '2D scalar known PUBLIC TRAIN material calibration; not N3 reconstruction or posterior inference',
                'posterior_mass': None, 'shared_color_limits': limits, 'anchors': []}
    for g, known, picks in selections:
        fig, axs = plt.subplots(2, 5, figsize=(15, 6.3), layout='constrained')
        labels = [('Known TRAIN anchor', None, known)] + picks
        for j, (name, r, image) in enumerate(labels):
            for i in range(2):
                ax = axs[i, j]
                if image is None:
                    ax.text(.5, .5, 'No compatible candidate', transform=ax.transAxes, ha='center')
                    if i == 0:
                        ax.set_title(name, fontsize=10)
                    ax.set_axis_off()
                    continue
                values = 1+image.real if i == 0 else image.imag
                ax.imshow(values, origin='lower', interpolation='nearest', cmap='viridis', vmin=limits[i][0], vmax=limits[i][1])
                ax.set_xticks([0, 12, 23]); ax.set_yticks([0, 12, 23])
                if j == 0:
                    ax.set_ylabel(('Re epsilon' if i == 0 else 'Im epsilon') + '\ny cell')
                if i == 1:
                    ax.set_xlabel('x cell')
                if i == 0:
                    subtitle = f"rW={g['base_residual']:.3f}" if r is None else f"rW={r['residual']:.3f}; material RMS={r['material_RMS_to_anchor']:.4g}"
                    ax.set_title(f'{name}\n{subtitle}', fontsize=10)
        for i in range(2):
            fig.colorbar(axs[i, 0].images[0], ax=axs[i, :], shrink=.75, label='relative permittivity')
        fig.suptitle(f"{g['anchor']}: full-forward verified local calibration candidates\nKnown original-P0 TRAIN material; eight-coordinate conditional chart, not general 2D material support")
        fig.supxlabel('Shared Re/Im scales across both anchors and all methods. Images selected by residual feasibility and material distance.\nCompatibility is not posterior mass; continuation is candidate-only; visible contrast is not a separated-topology certificate.', fontsize=9)
        save(fig, f"e2_materials_{g['anchor']}")
        portable['anchors'].append({'anchor': g['anchor'], 'anchor_epsilon_real': (1+known.real).tolist(),
            'anchor_epsilon_imag': known.imag.tolist(), 'candidates': [dict(method=m, record=r,
              epsilon_real=(1+x.real).tolist() if x is not None else None,
              epsilon_imag=x.imag.tolist() if x is not None else None) for m, r, x in picks]})
    (FIG / 'MATERIAL_GALLERY_DATA.json').write_text(json.dumps(portable, indent=2) + '\n')


def qualification():
    with (OUT / 'delegated/e0_2d/N3_ONLINE_CACHE_AUDIT.csv').open() as f:
        rows = list(csv.DictReader(f))
    ids = [r['case_id'].replace('_a1_snr20', '').replace('p2500', '') for r in rows]
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.7), layout='constrained')
    y = np.arange(len(rows))
    axs[0].barh(y, [float(r['whitened_residual']) for r in rows], color='#b66b5d')
    axs[0].axvline(float(rows[0]['tau']), color='black', ls='--', label='95% noise-ball radius')
    axs[0].set_yticks(y, ids); axs[0].set_xlabel('Whitened residual norm'); axs[0].legend(frameon=False)
    axs[1].barh(y, [float(r['phase_orth_relative']) for r in rows], color='#b66b5d')
    axs[1].set_yticks(y, ids); axs[1].set_xlabel('Relative phase-orthogonal material departure')
    fig.suptitle('E0: historical N3 is not a qualified original-P0 material anchor')
    fig.supxlabel('All audited N3 rows fail the noise ball and necessary phase rule. Prior membership is unknown, not proved by a box check.\nE2 therefore uses explicitly known TRAIN material calibration anchors; it does not repair or certify these N3 reconstructions.', fontsize=9)
    save(fig, 'e0_n3_qualification')


def derivative():
    d = read(OUT / 'results/E2_derivative_refinement/COMPLETE.json')
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.2), layout='constrained')
    axs[0].semilogy(np.arange(1, 9), d['material_singular_values'], 'o-', color='#194f73')
    axs[0].set_xlabel('Material-chart singular value index'); axs[0].set_ylabel('Singular value')
    axs[1].loglog([a['parameter_step_norm'] for a in d['audits']],
                   [a['relative_error'] for a in d['audits']], 'o-', color='#194f73')
    axs[1].set_xlabel('Prior-coordinate secant step norm'); axs[1].set_ylabel('Full secant vs analytic Jv relative error')
    fig.suptitle('Weak TRAIN anchor: gauge removed, nonlinear secant checked')
    fig.supxlabel('Material rank 7 of 8; one chart gauge. Refinement qualifies the derivative only and does not rerun method trajectories.\nSaved full-forward derivative check; no additional solve in plotting.', fontsize=9)
    save(fig, 'e2_derivative_refinement')


def main():
    start = time.perf_counter()
    FIG.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'figure.facecolor': 'white', 'savefig.facecolor': 'white'})
    e1(); candidates(); qualification(); derivative()
    panels = [
        ('代理后验的事件质量与材料 CDF', 'e1_proxy_event_cdf', '旧几何差分行已排除。四个已暴露观测、四个种子；代理后验的比较不等于全 Maxwell 校准。'),
        ('N3 的原先验与噪声资格', 'e0_n3_qualification', '残差可行和材料先验支持是两项不同要求；此处没有把后续已知 TRAIN 材料校准描述为 N3 重建。'),
        ('弱背景 TRAIN 材料', 'e2_materials_250033_weak', '只按在线残差兼容及与 anchor 的材料距离选代表；未按真值误差择优。候选概率未估计。'),
        ('环形背景 TRAIN 材料', 'e2_materials_250034_ring', '与上一张图共用实部和损耗色标。结构差异仍需独立资格判断，不能用低残差替代后验事件质量。'),
        ('导数与全非线性 secant', 'e2_derivative_refinement', '去除 chart gauge 后检查局部材料方向；细化检查没有追溯替换原方法轨迹。'),
    ]
    cards = ''.join(f'<section><h2>{html.escape(t)}</h2><p>{html.escape(s)}</p><a href="figures/{n}.pdf"><img src="figures/{n}.png" alt="{html.escape(t)}"></a></section>' for t,n,s in panels)
    page = '<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>LXZ 结果图库</title><style>body{font:17px/1.65 system-ui,sans-serif;max-width:1400px;margin:40px auto;padding:0 24px;color:#192f3f;background:#f5f7f8}section{background:white;padding:20px;margin:30px 0;border-radius:12px}img{width:100%;height:auto}h1,h2{line-height:1.3}a{color:#175e83}</style><h1>LXZ：概率代理与材料校准图库</h1><p>本页只读取已支付数组；绘图新增物理调用、训练和后验采样均为 0。二维材料图为已知公开 TRAIN anchor 周围的局部校准，不是重建 N3，也不赋予候选后验概率。</p><p><a href="START_HERE.md">研究入口</a> · <a href="delegated/p2_assets/AUDIT.md">资产审计</a> · <a href="delegated/probability_review/REVIEW.md">概率实现审阅</a></p>' + cards + '</html>'
    (OUT / 'GALLERY.html').write_text(page)
    receipt = dict(scope='saved-array visualization only', new_physics=0, new_proxy_evaluations=0,
                   new_training=0, new_posterior_sampling=0, wall_seconds=time.perf_counter()-start,
                   geometry_sources=['round1_geometry_repair', 'round2_geometry_repair'],
                   excluded_sources=['round1 geometry rows', 'round2 geometry rows'])
    receipt_path = FIG / 'PLOTTING_RECEIPT.json'
    previous = read(receipt_path) if receipt_path.exists() else None
    runs = previous['runs'] if previous and 'runs' in previous else ([previous] if previous else [])
    runs.append(receipt)
    receipt_path.write_text(json.dumps(dict(runs=runs, cumulative_wall_seconds=sum(r['wall_seconds'] for r in runs)), indent=2)+'\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
