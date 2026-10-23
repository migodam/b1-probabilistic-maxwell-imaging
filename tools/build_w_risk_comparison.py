#!/usr/bin/env python3
"""Plot two frozen offline summary files only: no GT, inference, PDE, or Book.

Fresh output directories are mandatory. Existing evidence is never overwritten.
Optional --run-id creates a new named subdirectory under final_comparison.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import platform
import resource
import sys
import time

def _stage_root():
    # Supports the repository layout and the adjacent stage/history Release layout.
    name='research_w_p0r_to_p1'
    for parent in Path(__file__).resolve().parents:
        for candidate in (parent,parent/'research'/name):
            if candidate.name==name and (candidate/'configs/AUTHORIZATION.json').is_file():
                return candidate
    raise RuntimeError('Research W evidence root was not found')

ROOT = _stage_root()
SUMMARY = 'results/P1B/offline_scores/SUMMARY.json'
PAIRED = 'results/P1B/offline_scores/PAIRED_CI.json'
ARMS = [
    ('zero_observation_prior_position_mean', '零观测先验均值', '#687584'),
    ('full_validated_continuous_point_control', 'Full 验证连续点控制', '#b86329'),
    ('continuous_quadrature_Bayes_mean', '连续求积后验均值', '#39846f'),
    ('direct_SMC_mean', '直接 SMC 后验均值', '#2c6a9b'),
    ('progressive_SMC_mean', '渐进 SMC 后验均值', '#795e9f'),
    ('direct_point', '直接网格点估计', '#b18b4a'),
    ('progressive_point', '渐进网格点估计', '#9b6b57'),
]
PAIRS = [
    ('multi_hypothesis', 'direct_point', 'direct_SMC_mean', '直接网格点 → 直接 SMC 均值'),
    ('multi_hypothesis', 'progressive_point', 'progressive_SMC_mean', '渐进网格点 → 渐进 SMC 均值'),
    ('strong_point_control', 'full_validated_continuous_point_control', 'direct_SMC_mean',
     'Full 验证连续点 → 直接 SMC 均值'),
    ('observation_information', 'zero_observation_prior_position_mean', 'direct_SMC_mean',
     '零观测先验均值 → 直接 SMC 均值'),
    ('progressive_schedule', 'direct_SMC_mean', 'progressive_SMC_mean',
     '直接 SMC 均值 → 渐进 SMC 均值'),
]
SCOPE = '受控新 fixed-support prior；原 minimum_feature 未解决；本图不作科学 GO 裁决。'


def number(value, where):
    if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value):
        raise ValueError(f'{where}: expected a finite number')
    return float(value)


def load_input(root, relative):
    path = root / relative
    before = path.stat()
    def reject(value):
        raise ValueError(f'Nonfinite JSON constant: {value}')
    value = json.loads(path.read_text(encoding='utf-8'), parse_constant=reject)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError('Frozen input changed while reading: ' + relative)
    return value, dict(path=relative, bytes=after.st_size, mtime_ns=str(after.st_mtime_ns))


def write_json(path, value):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')


def write_csv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open('x', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows({key: json.dumps(value, ensure_ascii=False) if isinstance(value, (dict, list)) else value
                         for key, value in row.items()} for row in rows)


def assemble(summary, paired):
    if not isinstance(summary, dict) or not isinstance(paired, dict):
        raise ValueError('SUMMARY and PAIRED_CI must both be objects')
    population = summary.get('actual_population', {})
    if (summary.get('parent_count') != 16 or population.get('observations') != 32
            or population.get('parents') != 16 or population.get('noise_replicas') != 2):
        raise ValueError('Expected the saved 16-parent, 32-observation, two-noise population')
    arm_rows = []
    for arm, label, _ in ARMS:
        record = summary['arm_summaries'][arm]
        mean = number(record['parent_means']['position_squared_error'], arm + '.MSE')
        coverage = record['complete_parent_coverage']['position_squared_error']
        if coverage != 16:
            raise ValueError(f'{arm}: missing complete parents; no silent selection')
        arm_rows.append(dict(arm=arm, label=label, parent_mean_MSE=mean,
            chart_MSE_times_1000=mean * 1000, complete_parent_count=coverage,
            qualified_subset_only=record.get('qualified_subset_only'),
            source=SUMMARY, source_pointer=f'#/arm_summaries/{arm}/parent_means/position_squared_error'))
    all_ci = []
    for family, rows in paired.items():
        if not isinstance(rows, list):
            raise ValueError(f'PAIRED_CI.{family} must be a list')
        for index, row in enumerate(rows):
            if not isinstance(row, dict):
                raise ValueError(f'PAIRED_CI.{family}[{index}] must be an object')
            interval = row.get('paired_percentile_CI95')
            if not isinstance(interval, list) or len(interval) != 2:
                raise ValueError(f'{family}[{index}]: missing saved two-sided interval')
            gain = number(row.get('gain_a_minus_b'), family + '.gain')
            lo, hi = (number(x, family + '.CI') for x in interval)
            if lo > hi:
                raise ValueError('Reversed paired interval')
            if (row.get('parent_count') != 16 or row.get('replicates') != 2000
                    or row.get('seed') != 2026102317 or row.get('missing_parent_ids') != []
                    or row.get('all_parent_ids') != list(range(16))):
                raise ValueError('Paired CI population or resampling contract differs')
            all_ci.append(dict(family=family, arm_a=row['arm_a'], arm_b=row['arm_b'],
                loss=row['loss'], gain_a_minus_b=gain, CI95_lower=lo, CI95_upper=hi,
                contains_zero=bool(lo <= 0 <= hi), parent_count=16, bootstrap_replicates=2000,
                bootstrap_seed=2026102317, status=row.get('status'), unit=row.get('unit'),
                source=PAIRED, source_pointer=f'#/{family}/{index}', plotted_primary=False))
    selected = []
    for family, arm_a, arm_b, label in PAIRS:
        matches = [row for row in all_ci if (row['family'], row['arm_a'], row['arm_b'], row['loss'])
                   == (family, arm_a, arm_b, 'position_squared_error')]
        if len(matches) != 1:
            raise ValueError(f'Missing/duplicate primary contrast: {family}/{arm_a}/{arm_b}')
        row = matches[0]
        row.update(plotted_primary=True, label=label)
        selected.append(row)
    if summary.get('contrasts') != paired:
        raise ValueError('SUMMARY.contrasts and PAIRED_CI disagree')
    return arm_rows, selected, all_ci


def plotting(table_dir):
    # Any Matplotlib font-cache output stays within the permitted output tree.
    os.environ['MPLCONFIGDIR'] = str(table_dir / '.mplconfig')
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import font_manager
    names = {font.name for font in font_manager.fontManager.ttflist}
    font = next((name for name in ('PingFang SC', 'Heiti SC', 'Arial Unicode MS') if name in names), None)
    if font is None:
        raise RuntimeError('No existing Chinese font is available; no fallback with missing glyphs')
    matplotlib.rcParams.update({'font.family': [font, 'DejaVu Sans'], 'font.size': 12,
        'axes.unicode_minus': False, 'svg.fonttype': 'path', 'figure.facecolor': '#fbfaf7',
        'axes.facecolor': '#fbfaf7', 'savefig.facecolor': '#fbfaf7',
        'axes.spines.top': False, 'axes.spines.right': False})
    import matplotlib.pyplot as plt
    return plt, font, matplotlib.__version__


def save(plt, figure, directory, stem):
    files = []
    for extension in ('png', 'svg'):
        path = directory / f'{stem}.{extension}'
        if path.exists():
            raise FileExistsError('Refusing figure overwrite: ' + str(path))
        figure.savefig(path, dpi=180, facecolor='#fbfaf7')
        files.append(path)
    plt.close(figure)
    return files


def draw_mse(plt, rows, directory):
    figure, ax = plt.subplots(figsize=(11.8, 7.5))
    figure.subplots_adjust(left=.30, right=.90, top=.81, bottom=.24)
    values = [row['chart_MSE_times_1000'] for row in rows]
    ax.barh(range(len(rows)), values, height=.65, color=[item[2] for item in ARMS])
    ax.set_yticks(range(len(rows)), [row['label'] for row in rows], fontsize=12)
    ax.invert_yaxis()
    ax.set_xlim(0, max(values) * 1.17)
    ax.set_xlabel('parent 平均平方位置误差 MSE（×10⁻³，越小越好）', labelpad=12)
    ax.xaxis.grid(True, color='#d9dfe3', alpha=.65)
    ax.set_axisbelow(True)
    for index, value in enumerate(values):
        ax.text(value + .045, index, f'{value:.4f}', va='center', fontsize=12)
    ax.axvline(values[0], color='#687584', linestyle=':', linewidth=1.3)
    figure.text(.06, .94, '七臂 MSE：完整对照', fontsize=22, weight='bold')
    figure.text(.06, .885, '16 个 parent；每 parent 先平均 2 次噪声；SMC 为 4 seed 混合后验均值行动。', fontsize=11.5, color='#46566a')
    figure.text(.06, .145, '虚线为零观测先验均值的本批风险；七臂全部展示，未按结果方向筛选。', fontsize=10.5)
    figure.text(.06, .10, '仅展示冻结 SUMMARY 的风险，不比较成本；bank / setup 均已付费，不代表端到端节省。', fontsize=10.5)
    figure.text(.06, .055, SCOPE, fontsize=10, color='#46566a')
    return save(plt, figure, directory, '01_seven_arm_MSE')


def draw_ci(plt, rows, directory):
    figure, ax = plt.subplots(figsize=(15.8, 7.2))
    figure.subplots_adjust(left=.31, right=.71, top=.79, bottom=.25)
    span = max(abs(row[key] * 1000) for row in rows for key in ('CI95_lower', 'CI95_upper')) * 1.17
    for index, row in enumerate(rows):
        lo, hi, gain = (row[key] * 1000 for key in ('CI95_lower', 'CI95_upper', 'gain_a_minus_b'))
        color = '#795e9f' if row['family'] == 'progressive_schedule' else '#2c6a9b'
        ax.hlines(index, lo, hi, color=color, linewidth=2.4)
        ax.vlines((lo, hi), index - .075, index + .075, color=color, linewidth=1.3)
        ax.scatter(gain, index, s=48, color=color, zorder=3)
        ax.text(1.04, index, f'{gain:+.4f}  [{lo:+.4f}, {hi:+.4f}]',
                transform=ax.get_yaxis_transform(), va='center', fontsize=10.5)
    ax.axvline(0, color='#596775', linestyle='--', linewidth=1.1)
    ax.set_xlim(-span, span)  # preserve negative and positive portions, no truncation at zero
    ax.set_ylim(len(rows) - .45, -.55)
    ax.set_yticks(range(len(rows)), [row['label'] for row in rows], fontsize=11.5)
    ax.xaxis.grid(True, color='#d9dfe3', alpha=.65)
    ax.set_axisbelow(True)
    ax.set_xlabel('配对 MSE 增益：左方法损失 − 右方法损失（×10⁻³）', labelpad=12)
    figure.text(.06, .94, '主 MSE 配对区间：保留全部正负结果', fontsize=22, weight='bold')
    figure.text(.06, .885, '正值有利右方法；16 个 parent 配对重采样 2000 次，seed = 2026102317。', fontsize=12, color='#46566a')
    figure.text(.735, .82, '增益  [95% 区间]，均 ×10⁻³', fontsize=10.5, color='#46566a')
    figure.text(.06, .145, '区间直接读取 PAIRED_CI；未重新 bootstrap、未裁掉负区间或跨零区间。', fontsize=10.5)
    figure.text(.06, .10, '配对单位是 parent，2 次噪声已先平均；4 个 SMC seed 不另算独立 parent。', fontsize=10.5)
    figure.text(.06, .055, SCOPE, fontsize=10, color='#46566a')
    return save(plt, figure, directory, '02_primary_MSE_paired_CI')


def main(argv=None):
    cpu0, wall0 = time.process_time(), time.perf_counter()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--run-id', help='optional fresh subdirectory name; never overwrite old comparisons')
    args = parser.parse_args(argv)
    root = args.root.resolve()
    figure_dir, table_dir = root / 'figures/final_comparison', root / 'tables/final_comparison'
    if args.run_id:
        if not args.run_id.replace('-', '').replace('_', '').isalnum():
            parser.error('--run-id must contain letters/digits/hyphens/underscores')
        figure_dir, table_dir = figure_dir / args.run_id, table_dir / args.run_id
    if figure_dir.exists() or table_dir.exists():
        parser.error('Output already exists; refusing overwrite. Use a fresh --run-id if a new build is authorized.')
    figure_dir.mkdir(parents=True, exist_ok=False)
    table_dir.mkdir(parents=True, exist_ok=False)
    log_path = table_dir / 'POSTPROCESS_LOG.json'
    log = dict(started_utc=datetime.now(timezone.utc).isoformat(),
        task='frozen offline risk comparison figures', read_allowlist=[SUMMARY, PAIRED],
        new_RHS=0, new_factorizations=0, new_PDE=0, Book_imported=False,
        production_main_imported=False, GT_read=False, inference_run=False,
        frozen_results_modified=False, ledger_modified=False,
        scientific_decision='NOT_ASSIGNED', cost_replay_read=False,
        end_to_end_cost_savings_claim=False)
    try:
        summary, source_summary = load_input(root, SUMMARY)
        paired, source_paired = load_input(root, PAIRED)
        rows, selected, all_ci = assemble(summary, paired)
        write_csv(table_dir / 'seven_arm_MSE.csv', rows)
        write_csv(table_dir / 'primary_MSE_paired_CI.csv', selected)
        write_csv(table_dir / 'all_saved_paired_CI.csv', all_ci)
        plt, font, version = plotting(table_dir)
        mse_files = draw_mse(plt, rows, figure_dir)
        ci_files = draw_ci(plt, selected, figure_dir)
        sources = dict(read_inputs=[source_summary, source_paired],
            array_or_GT_input=False, calculations='MSE and CI unchanged; axis values only scaled by1000',
            figure_sources=[
                dict(files=[str(p.relative_to(root)) for p in mse_files],
                     table=str((table_dir / 'seven_arm_MSE.csv').relative_to(root)), rows=rows),
                dict(files=[str(p.relative_to(root)) for p in ci_files],
                     table=str((table_dir / 'primary_MSE_paired_CI.csv').relative_to(root)), rows=selected)],
            all_contrasts_table=str((table_dir / 'all_saved_paired_CI.csv').relative_to(root)),
            omitted_by_sign=False, plotted_MSE_arms=7, plotted_MSE_contrasts=len(selected),
            all_saved_contrasts_exported=len(all_ci), scope=SCOPE)
        for item in sources['read_inputs']:
            current = (root / item['path']).stat()
            if (current.st_size, str(current.st_mtime_ns)) != (item['bytes'], item['mtime_ns']):
                raise RuntimeError('Frozen input changed during rendering: ' + item['path'])
        write_json(table_dir / 'FIGURE_SOURCES.json', sources)
        log.update(status='COMPLETE_0RHS_POSTPROCESS_ONLY', font=font,
            matplotlib_version=version, input_sources=sources['read_inputs'],
            generated_images=[dict(path=str(p.relative_to(root)), bytes=p.stat().st_size)
                              for p in mse_files + ci_files],
            tables=[str(p.relative_to(root)) for p in table_dir.glob('*.csv')])
    except BaseException as exc:
        log.update(status='FAILED_POSTPROCESS_ONLY', error_type=type(exc).__name__, error=str(exc))
        raise
    finally:
        usage = resource.getrusage(resource.RUSAGE_SELF)
        log.update(cpu_seconds=time.process_time() - cpu0,
            wall_seconds=time.perf_counter() - wall0,
            peak_RSS_bytes=int(usage.ru_maxrss * (1 if sys.platform == 'darwin' else 1024)),
            RSS_source='resource.getrusage(RUSAGE_SELF).ru_maxrss; Darwin bytes, Linux KiB converted',
            platform=platform.platform(), finished_utc=datetime.now(timezone.utc).isoformat(),
            output_directories_fresh=True)
        write_json(log_path, log)
    print(json.dumps(dict(status=log['status'], figures=str(figure_dir), tables=str(table_dir),
                         cpu_seconds=log['cpu_seconds'], wall_seconds=log['wall_seconds'],
                         peak_RSS_bytes=log['peak_RSS_bytes'], new_RHS=0), ensure_ascii=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
