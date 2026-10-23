#!/usr/bin/env python3
"""Render saved Research W records without importing or running any physics.

Inputs: RESULTS.json and complete records in ITERATIONS.json, plus FINAL.npz
only beside an existing RESULTS.json. Never reads trial field/J caches or logs.
--final changes the presentation label; it does not certify any scientific gate.
"""
from __future__ import annotations
import argparse
import csv
import html
import io
import json
import math
import os
from pathlib import Path
from datetime import datetime, timezone


def atomic(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(text, encoding='utf-8')
    os.replace(temporary, path)


def reject_constant(value):
    raise ValueError(f'nonfinite JSON number {value}')


def read_json(path, root, sources, errors):
    try:
        raw = path.read_text(encoding='utf-8')
        obj = json.loads(raw, parse_constant=reject_constant)
        def finite_tree(value):
            if isinstance(value, float) and not math.isfinite(value):
                raise ValueError('nonfinite numeric value in JSON tree')
            if isinstance(value, dict):
                for item in value.values(): finite_tree(item)
            if isinstance(value, list):
                for item in value: finite_tree(item)
        finite_tree(obj)
        sources.append(dict(path=str(path.relative_to(root)), bytes=len(raw.encode()),
                            modified_utc=datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat()))
        return obj
    except (OSError, ValueError) as exc:
        errors.append(dict(path=str(path.relative_to(root)), kind=type(exc).__name__, message=str(exc)))
        return None


def number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def scalar_rows(rows, source):
    return [dict(source=source, **{k: v for k, v in row.items()
                                 if not isinstance(v, (list, dict))}) for row in rows]


def write_csv(path, rows):
    if not rows:
        atomic(path, 'status\nNOT_RUN\n')
        return
    keys = list(dict.fromkeys(k for row in rows for k in row))
    f = io.StringIO(newline='')
    writer = csv.DictWriter(f, keys)
    writer.writeheader()
    writer.writerows([{k: json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v
                      for k, v in row.items()} for row in rows])
    atomic(path, f.getvalue())


def table(rows, columns=None):
    if not rows:
        return '<p class="empty">NOT_RUN — no completed saved records.</p>'
    keys = columns or list(dict.fromkeys(k for r in rows for k in r))
    def cell(v):
        if isinstance(v, float):
            return f'{v:.7g}'
        return str(v) if v is not None else '—'
    return '<div class="scroll"><table><thead><tr>' + ''.join(
        f'<th>{html.escape(k)}</th>' for k in keys) + '</tr></thead><tbody>' + ''.join(
            '<tr>' + ''.join(f'<td>{html.escape(cell(r.get(k)))}</td>' for k in keys) + '</tr>'
            for r in rows) + '</tbody></table></div>'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--final', action='store_true', help='Use final snapshot label; no stage/gate changes.')
    parser.add_argument('--online-prior-lambda', type=float, default=.002,
                        help='Frozen online prior weight (.002), used with recorded x and objective; override only with matching run metadata.')
    args = parser.parse_args()
    root = args.root.resolve()
    generated = root / 'figures' / 'saved_results'
    tables_dir = root / 'tables' / 'saved_results'
    generated.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    errors, sources, stage_rows, stage_data = [], [], [], {}
    runs, histories, iteration_rows = [], {}, []

    for stage in ('R0', 'R1', 'R2'):
        path = root / 'results' / stage / 'RESULTS.json'
        if not path.exists():
            stage_rows.append(dict(stage=stage, saved_status='NOT_RUN', records=0))
            stage_data[stage] = []
            write_csv(tables_dir / f'{stage}.csv', [])
            continue
        obj = read_json(path, root, sources, errors)
        if obj is None:
            stage_rows.append(dict(stage=stage, saved_status='PARSE_ERROR', records=0))
            stage_data[stage] = []
            continue
        rows = obj.get('rows') if isinstance(obj, dict) else obj
        if not isinstance(rows, list) or not all(isinstance(r, dict) for r in rows):
            errors.append(dict(path=str(path.relative_to(root)), kind='SchemaError', message='expected list of dictionaries, or dict.rows'))
            rows = []
        stage_data[stage] = rows
        stage_rows.append(dict(stage=stage, saved_status=obj.get('status', 'SAVED_RECORDS') if isinstance(obj, dict) else 'SAVED_RECORDS', records=len(rows)))
        write_csv(tables_dir / f'{stage}.csv', scalar_rows(rows, str(path.relative_to(root))))

    r3 = root / 'results' / 'R3'
    folders = sorted(set(p.parent for pattern in ('**/RESULTS.json', '**/ITERATIONS.json')
                         for p in r3.glob(pattern))) if r3.exists() else []
    for folder in folders:
        run = str(folder.relative_to(r3))
        identity = 'ORACLE_DIAGNOSTIC' if run.startswith('oracle_') else 'ONLINE_INFERENCE' if run.startswith('online_') else 'UNDECLARED'
        final_path = folder / 'RESULTS.json'
        summary = None
        if final_path.exists():
            summary = read_json(final_path, root, sources, errors)
            if summary is not None and not isinstance(summary, dict):
                errors.append(dict(path=str(final_path.relative_to(root)), kind='SchemaError', message='R3 final result must be a dictionary'))
                summary = None
        if summary is not None:
            identity = summary.get('identity', identity)
        result = dict(run=run, identity=identity, saved_status='SAVED_RESULT' if summary is not None else 'PARSE_ERROR' if final_path.exists() else 'IN_PROGRESS')
        if summary is not None:
            result.update({k: v for k, v in summary.items() if not isinstance(v, (list, dict))})
        runs.append(result)
        iterations = folder / 'ITERATIONS.json'
        if not iterations.exists():
            continue
        events = read_json(iterations, root, sources, errors)
        if events is None:
            continue
        if not isinstance(events, list) or not all(isinstance(e, dict) for e in events):
            errors.append(dict(path=str(iterations.relative_to(root)), kind='SchemaError', message='iteration stream must be a list of dictionaries'))
            continue
        history = []
        accepted = {e.get('residual_call_id', e.get('call_id')) for e in events if e.get('event') == 'trial_accepted'}
        for e in events:
            if e.get('event') != 'evaluation_complete':
                continue  # A started/in-flight/failed evaluation never becomes a plotted point.
            if not number(e.get('fun')) or not number(e.get('call_id')):
                errors.append(dict(path=str(iterations.relative_to(root)), kind='SchemaError', message='evaluation_complete lacks finite fun/call_id'))
                continue
            point = dict(run=run, identity=identity, call_id=e['call_id'], objective=e['fun'],
                         iteration=e.get('iteration'), phase=e.get('phase', ''),
                         accepted_trial=e['call_id'] in accepted,
                         objective_norm=math.sqrt(max(0., 2 * e['fun'])))
            direct = next((e[k] for k in ('data_residual_norm', 'residual_norm', 'd_mean', 'd_cal') if number(e.get(k))), None)
            if direct is not None:
                point.update(data_residual_norm=direct, residual_source='saved explicit field')
            else:
                lam = e.get('prior_lambda')
                basis = 'saved prior_lambda'
                if not number(lam) and identity == 'ORACLE_DIAGNOSTIC':
                    lam, basis = 0., 'oracle-fit implementation prior=0; objective reconstruction'
                if not number(lam) and identity == 'ONLINE_INFERENCE' and args.online_prior_lambda is not None:
                    lam, basis = args.online_prior_lambda, 'frozen online prior=.002, zero center; CLI override recorded' if args.online_prior_lambda == .002 else 'explicit CLI prior override, zero center'
                x = e.get('x')
                if number(lam) and lam >= 0 and isinstance(x, list) and all(number(t) for t in x):
                    squared = 2 * e['fun'] - lam * math.fsum(t*t for t in x)
                    if squared >= -1e-10 * max(1., abs(e['fun'])):
                        point.update(data_residual_norm=math.sqrt(max(0., squared)), residual_source=basis)
                    else:
                        errors.append(dict(path=str(iterations.relative_to(root)), kind='MetricError', message=f'negative reconstructed residual square at call {e["call_id"]}'))
            history.append(point)
        histories[run] = history
        iteration_rows.extend(history)
        result['completed_evaluations'] = len(history)
        result['unfinished_evaluations'] = sum(e.get('event') == 'evaluation_start' for e in events) - len(history) - sum(e.get('event') == 'evaluation_failed' for e in events)

    stage_rows.append(dict(stage='R3', saved_status='SAVED_RESULT' if any(r['saved_status'] == 'SAVED_RESULT' for r in runs) else 'IN_PROGRESS' if runs else 'NOT_RUN', records=sum(r['saved_status'] == 'SAVED_RESULT' for r in runs)))
    later_rows = []
    for stage in ('B2', 'C', 'P1A', 'P1B'):
        paths = sorted((root / 'results' / stage).glob('**/RESULTS.json'))
        valid = []
        for path in paths:
            data = read_json(path, root, sources, errors)
            if data is not None:
                valid.append(dict(source=str(path.relative_to(root)), status=data.get('status', 'SAVED_RECORDS') if isinstance(data, dict) else 'SAVED_RECORDS'))
        later_rows.append(dict(stage=stage, status='SAVED_RECORDS' if valid else 'NOT_RUN', saved_result_files=len(valid), sampler_metrics='reported only if explicitly saved; no samples inferred'))
    projection = [dict(seed=r.get('seed'), band=r.get('band'), grid='×'.join(map(str, r.get('grid', []))),
                       subsamples=r.get('subsamples'), d_repr=r.get('d_repr'), occupied_area=r.get('occupied_area'),
                       identity=r.get('identity'), minimum_feature=r.get('minimum_feature'), status=r.get('status', 'SAVED_RECORD'))
                  for r in stage_data.get('R1', [])]
    for name, data in [('stage_status', stage_rows), ('R3_runs', runs), ('completed_evaluations', iteration_rows), ('projection_distance', projection), ('later_stages', later_rows)]:
        write_csv(tables_dir / f'{name}.csv', data)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})
    figures = []
    def finish(fig, name, title, caption):
        fig.tight_layout()
        for suffix in ('png', 'svg'):
            fig.savefig(generated / f'{name}.{suffix}', dpi=170, facecolor='white')
        plt.close(fig)
        figures.append(dict(path=str((generated/f'{name}.png').relative_to(root)), svg=str((generated/f'{name}.svg').relative_to(root)), title=title, caption=caption))
    def placeholder(ax, text):
        ax.text(.5, .5, text, transform=ax.transAxes, ha='center', va='center', wrap=True)
        ax.set_xticks([]); ax.set_yticks([])

    fig, ax = plt.subplots(figsize=(9, 4.5))
    found = False
    for run, hist in histories.items():
        points = [p for p in hist if 'data_residual_norm' in p]
        if points:
            ax.plot([p['call_id'] for p in points], [p['data_residual_norm'] for p in points], 'o-', ms=4, label=run)
            found = True
    if found:
        ax.axhline(1., color='#b45309', ls=':', lw=1, label='d=1 reference line; no gate judgment')
        ax.legend(fontsize=8)
        values = [p['data_residual_norm'] for hist in histories.values() for p in hist if 'data_residual_norm' in p]
        if min(values) > 0:
            ax.set_yscale('log');ax.set_ylim(min(1., min(values))*.7, max(1.,max(values))*1.4)
        else:
            ax.set_yscale('symlog', linthresh=.1);ax.set_ylim(bottom=0)
    else:
        placeholder(ax, 'NOT_RUN: no completed evaluations with recoverable data residual.\nObjective values are shown separately.')
    ax.set_xlabel('Completed evaluator call ID (not RHS or accepted iteration)')
    ax.set_ylabel('Whitened data residual norm')
    ax.set_title('Saved completed evaluations: data residual')
    finish(fig, 'residual_vs_evaluation', '残差与已完成评价次数', '只画 evaluation_complete。横轴为 evaluator call ID，包括 accepted refresh；不是独立材料数、接受步数或物理 RHS。重建来源见 CSV residual_source。')

    fig, ax = plt.subplots(figsize=(9, 4.5))
    for run, hist in histories.items():
        if hist:
            ax.plot([p['call_id'] for p in hist], [p['objective'] for p in hist], 'o-', ms=4, label=run)
    if any(histories.values()):
        vals = [p['objective'] for hist in histories.values() for p in hist]
        if min(vals) > 0:
            ax.set_yscale('log');ax.set_ylim(min(vals)*.7,max(vals)*1.4)
        else:
            ax.set_yscale('symlog', linthresh=.1);ax.set_ylim(bottom=0)
        ax.legend(fontsize=8)
    else:
        placeholder(ax, 'NOT_RUN: no completed objective evaluations')
    ax.set_xlabel('Completed evaluator call ID');ax.set_ylabel('Saved objective fun (prior may be included)')
    finish(fig, 'objective_vs_evaluation', '目标函数轨迹', '直接使用已保存 fun；含正则先验时它与纯数据残差不同。')

    fig, ax = plt.subplots(figsize=(9, 4.5))
    points = [r for r in runs if r['saved_status'] == 'SAVED_RESULT' and number(r.get('d_mean'))]
    if points:
        colors = ['#c2410c' if r['identity'] == 'ORACLE_DIAGNOSTIC' else '#0369a1' for r in points]
        ax.bar(range(len(points)), [r['d_mean'] for r in points], color=colors)
        ax.set_xticks(range(len(points)), [r['run'] for r in points], rotation=25, ha='right', fontsize=8)
        ax.set_ylabel('Final saved d_mean')
    else:
        placeholder(ax, 'NOT_RUN: no completed R3 RESULTS.json endpoints.\nIn-progress traces do not count as final outcomes.')
    ax.set_title('Oracle diagnostics / measurement-only online endpoints')
    finish(fig, 'oracle_vs_online', 'Oracle 与在线端点', '橙色 oracle 使用真值派生初值或支持，仅为离线表达/优化诊断；蓝色在线结果单独标记。无一类结果时不作虚构对照或优势判断。')

    fig, ax = plt.subplots(figsize=(9, 4.5))
    groups = {}
    for r in projection:
        if number(r.get('d_repr')):
            groups.setdefault(f"seed{r['seed']} / {r['band']}", []).append(r)
    if groups:
        for label, rs in groups.items():
            ax.plot(range(len(rs)), [r['d_repr'] for r in rs], 'o-', label=label)
        first = next(iter(groups.values()))
        ax.set_xticks(range(len(first)), [f"{r['grid']} / ss{r['subsamples']}" for r in first])
        ax.legend(fontsize=8);ax.set_ylabel('d_repr: whitened forward data distance')
    else:
        placeholder(ax, 'NOT_RUN: no saved R1 projection-distance values')
    ax.set_title('Offline oracle projection diagnostics; no legal-witness certification')
    finish(fig, 'projection_distance', '真值投影后的数据距离', 'd_repr 来自 R1 原始结果，是正演数据距离而非材料图像误差。网格与亚像素平均不同；像素投影未认证 minimum_feature。')

    fig, ax = plt.subplots(figsize=(9, 4.5))
    cases = {}
    for r in stage_data.get('R2', []):
        if number(r.get('step')) and number(r.get('relative_error')):
            cases.setdefault(f"{r.get('case')} / {r.get('scheme')}", []).append(r)
    if cases:
        for label, rs in cases.items():
            ax.loglog([r['step'] for r in rs], [r['relative_error'] for r in rs], 'o-', label=label)
        ax.legend(fontsize=8)
    else:
        placeholder(ax, 'NOT_RUN: no saved R2 derivative errors')
    ax.set_xlabel('Finite difference step');ax.set_ylabel('Saved relative derivative error')
    finish(fig, 'derivative_errors', '固定离散导数检查', '中心差分与边界单边差分分开显示；这不是投影、拓扑或形状可微性认证。')

    # Optional completed endpoint material plots; no generator or field cache is read.
    for r in runs:
        folder = r3 / r['run']
        final = folder / 'FINAL.npz'
        if r['saved_status'] != 'SAVED_RESULT' or not final.exists():
            continue
        try:
            import numpy as np
            with np.load(final, allow_pickle=False) as arrays:
                if 'full_parameters' not in arrays:
                    continue
                parameters = np.asarray(arrays['full_parameters'])
                final_x = np.asarray(arrays['x']) if 'x' in arrays else None
            if final_x is not None and number(r.get('objective')) and number(r.get('d_mean')):
                lam = 0. if r['identity'] == 'ORACLE_DIAGNOSTIC' else args.online_prior_lambda if r['identity'] == 'ONLINE_INFERENCE' else None
                if lam is not None:
                    squared = 2*r['objective'] - lam*float(final_x@final_x)
                    reconstructed = math.sqrt(max(0., squared))
                    discrepancy = abs(reconstructed-r['d_mean'])
                    r['objective_residual_reconstruction_difference'] = discrepancy
                    r['objective_residual_reconstruction_status'] = 'MATCH' if discrepancy <= 1e-7*max(1.,r['d_mean']) else 'MISMATCH'
                    if r['objective_residual_reconstruction_status'] == 'MISMATCH':
                        errors.append(dict(path=str(final.relative_to(root)), kind='MetricError', message=f'final objective/prior reconstructed residual differs from saved d_mean by {discrepancy}'))
            if parameters.shape != (3072,) or not np.isfinite(parameters).all():
                raise ValueError('material rendering supports recorded full 48x32 two-field vector only; other shape requires explicit metadata')
            fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)
            for ax, a, label, upper in zip(axes, (parameters[:1536], parameters[1536:]), ('dielectric contrast', 'conductivity'), (8., .15)):
                im = ax.imshow(a.reshape(32,48), origin='lower', extent=(-2,2,0,3), vmin=0, vmax=upper, aspect='equal', cmap='viridis')
                ax.set_title(label);ax.set_xlabel('x');fig.colorbar(im, ax=ax, shrink=.8)
            axes[0].set_ylabel('y')
            fig.suptitle(f"{r['identity']} | {r['run']}", fontsize=9)
            safe = r['run'].replace('/', '_')
            finish(fig, f'material_{safe}', f"保存材料：{r['run']}", '读取 completed endpoint FINAL.npz/full_parameters 原值。ORACLE_DIAGNOSTIC 材料严格离线展示；颜色统一使用冻结幅度范围，不据图认证合法性。')
            sources.append(dict(path=str(final.relative_to(root)), role='saved material only'))
        except (OSError, ValueError, KeyError) as exc:
            errors.append(dict(path=str(final.relative_to(root)), kind=type(exc).__name__, message=str(exc)))

    write_csv(tables_dir/'R3_runs.csv', runs)  # include completed objective/residual checks
    status = 'FINAL_SNAPSHOT' if args.final else 'WORK_IN_PROGRESS'
    now = datetime.now(timezone.utc).isoformat()
    payload = dict(report_status=status, generated_utc=now, physics_actions_by_builder=0,
                   stage_status=stage_rows, R3_runs=runs, projection_distance=projection,
                   later_stages=later_rows, input_sources=sources, parse_errors=errors,
                   frozen_online_prior_lambda=args.online_prior_lambda, figures=figures,
                   limitations=['Saved-record presentation only; no scientific gate/novelty judgment.',
                                'Syntactically valid iteration snapshots may include unfinished starts; only evaluation_complete is plotted.',
                                'Residual reconstructed as sqrt(2 fun-lambda ||x||^2): recorded objective/x, zero prior center, frozen oracle lambda0 / online lambda.002. CLI override is explicit.',
                                'No samples, posterior weights, ESS, chains or coverage inferred from optimizer endpoints.'])
    atomic(root/'RESULTS_REPORT.json', json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False))
    atomic(root/'REPORT_PARSE_ERRORS.json', json.dumps(dict(status='ERROR' if errors else 'OK', errors=errors), ensure_ascii=False, indent=2))
    figure_html = ''.join(f'<section><h2>{html.escape(f["title"])}</h2><a href="{html.escape(f["svg"])}"><img src="{html.escape(f["path"])}" alt="{html.escape(f["title"])}"></a><p class="note">{html.escape(f["caption"])}</p></section>' for f in figures)
    document = f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Research W · saved results</title>
<style>body{{font:16px/1.65 system-ui,sans-serif;background:#f3f5f7;color:#17212e;margin:0}}main{{max-width:1120px;margin:auto;padding:36px 24px}}h1{{font-size:30px;margin:8px 0}}h2{{font-size:21px}}section{{background:white;border:1px solid #dbe2e9;border-radius:12px;padding:20px;margin:20px 0}}.badge{{display:inline-block;background:#fff0cf;color:#7c4900;padding:4px 10px;border-radius:6px;font-weight:700}}.note,.meta{{color:#506173;font-size:14px}}.scroll{{overflow:auto}}table{{border-collapse:collapse;width:100%;font-size:13px}}td,th{{text-align:left;border-bottom:1px solid #e2e8f0;padding:8px;vertical-align:top;white-space:nowrap}}th{{background:#f6f8fa}}img{{display:block;width:100%;height:auto}}a{{color:#0369a1}}.empty{{color:#6b7280}}.error{{color:#b91c1c}}</style>
<main><span class="badge">{status}</span><h1>Research W · 保存结果快照</h1><p>从 0 读起：本页先展示真值回放与离线表达诊断，再展示在线反演的已保存轨迹。所有图表由现有数据生成；报告构建新增物理动作 0。</p><p class="meta">生成时间 UTC：{now} · <a href="RESULTS_REPORT.json">机器可读报告</a> · <a href="REPORT_PARSE_ERRORS.json">解析错误记录</a></p>
<section><h2>阶段记录状态</h2>{table(stage_rows)}<p class="note">状态来自文件中的记录，不是本页作出的科学门判决。WORK_IN_PROGRESS 与 FINAL_SNAPSHOT 只区分报告展示模式。不存在的后续阶段保持 NOT_RUN。</p></section>
<section><h2>Oracle / 在线结果清单</h2>{table(runs, ['run','identity','saved_status','d_mean','projected_gradient','iterations','success','stop','stopping','completed_evaluations','unfinished_evaluations','objective_residual_reconstruction_status'])}<p class="note">Oracle 使用真值信息，仅用于离线诊断；在线结果独立列出。projected_gradient 单独显示，数据残差阈值与 KKT 驻点不等同。正在进行的轨迹不被补成最终端点。残差重建公式为 sqrt(2 fun − λ||x||²)，先验中心为零；当前冻结 oracle λ=0、online λ={args.online_prior_lambda}，已完成端点与保存 d_mean 交叉核对。<a href="tables/saved_results/R3_runs.csv">CSV</a></p></section>
<section><h2>原始真值回放</h2>{table([r for r in stage_data.get('R0',[]) if 'd_mean' in r], ['seed','band','label','d_mean','d_observed_truth','volume_probe','same_contract_pass','identity'])}<p class="note">d_observed_truth 是保存噪声向量的白化长度，d_mean 是原均值数据的回放距离；两者不混合。</p></section>
{figure_html}<section><h2>投影距离明细 · 严格离线</h2>{table(projection)}<p><a href="tables/saved_results/projection_distance.csv">CSV</a> · <a href="tables/saved_results/completed_evaluations.csv">完整已完成评价表</a></p></section>
<section><h2>后续见证 / 采样阶段</h2>{table(later_rows)}<p class="note">不存在真实结果时，没有 posterior、链、ESS、模态权重或覆盖率图。本页不由优化候选数量构造概率。</p></section>
<section><h2>数据解析与来源</h2><p class="{'error' if errors else 'note'}">{'存在解析错误；相关文件未用于对应图表。' if errors else '当前读取文件均通过 JSON 解析；这不替代科学审阅。'}</p>{table(errors, ['path','kind','message']) if errors else ''}{table(sources, ['path','bytes','modified_utc','role'])}</section></main></html>'''
    for name in ('report.html', 'index.html'):
        atomic(root/name, document)
    print(json.dumps(dict(status='ERROR' if errors else status, completed_evaluations=len(iteration_rows),
                          final_R3_results=sum(r['saved_status']=='SAVED_RESULT' for r in runs),
                          plots=len(figures), parse_errors=len(errors), physics_actions=0,
                          index=str(root/'index.html')), ensure_ascii=False))
    return 2 if errors else 0


if __name__ == '__main__':
    raise SystemExit(main())
