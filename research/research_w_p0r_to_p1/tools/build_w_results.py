#!/usr/bin/env python3
"""Render saved Research W records without importing or running any physics.

Inputs: RESULTS.json and complete records in ITERATIONS.json/FWI_HISTORY.json,
R_SEALED.json for resource status, plus FINAL.npz
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
    sealed_path=root/'R_SEALED.json'
    sealed=read_json(sealed_path,root,sources,errors) if sealed_path.exists() else None
    control_map={}
    controls_path=root/'results/R3/controls/RESULTS.json'
    if controls_path.exists():
        records=read_json(controls_path,root,sources,errors)
        if isinstance(records,list) and all(isinstance(r,dict) for r in records):
            control_map={f"controls/{r['kind']}/{r['method']}":r for r in records}
        else:
            errors.append(dict(path=str(controls_path.relative_to(root)),kind='SchemaError',message='controls RESULTS must be a list of dictionaries'))

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
    folders = sorted(set(p.parent for pattern in ('**/RESULTS.json', '**/ITERATIONS.json','**/FWI_HISTORY.json')
                         for p in r3.glob(pattern))) if r3.exists() else []
    for folder in folders:
        run = str(folder.relative_to(r3))
        if run=='controls':continue
        role='CONTROL' if run.startswith('controls/') else 'FREQUENCY_SUBSET_STAGE' if folder.name.startswith('stage_') else 'ALL_BAND_ENDPOINT'
        identity = 'CONTROLLED_OPTIMIZER_TEST_NOT_UNKNOWN_WALL' if role=='CONTROL' else 'ORACLE_DIAGNOSTIC' if run.startswith('oracle_') else 'ONLINE_INFERENCE' if run.startswith('online_') else 'UNDECLARED'
        final_path = folder / 'RESULTS.json'
        summary = control_map.get(run)
        if final_path.exists():
            summary = read_json(final_path, root, sources, errors)
            if summary is not None and not isinstance(summary, dict):
                errors.append(dict(path=str(final_path.relative_to(root)), kind='SchemaError', message='R3 final result must be a dictionary'))
                summary = None
        if summary is not None:
            identity = summary.get('identity', identity)
        result = dict(run=run, identity=identity, role=role, saved_status='SAVED_RESULT' if summary is not None else 'PARSE_ERROR' if final_path.exists() else 'UNCOMPLETED_R_SEALED' if sealed else 'IN_PROGRESS')
        if summary is not None:
            result.update({k: v for k, v in summary.items() if not isinstance(v, (list, dict))})
            for field in ('d_mean','all_wide_d_mean','d_mean_subset'):
                if number(summary.get(field)):
                    result.update(reported_data_distance=summary[field], data_distance_field=field);break
            if isinstance(summary.get('frequencies'),list):result['frequency_scope']=summary['frequencies']
        runs.append(result)
        iterations = folder / 'ITERATIONS.json'
        if not iterations.exists():
            iterations=folder/'FWI_HISTORY.json'
        if not iterations.exists():continue
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
                lam = e.get('prior_lambda', summary.get('prior_lambda') if summary else None)
                basis = 'saved prior_lambda'
                if not number(lam) and role=='CONTROL':lam,basis=0.,'known control optimizer prior=0'
                if not number(lam) and identity == 'ORACLE_DIAGNOSTIC':
                    lam, basis = 0., 'oracle-fit implementation prior=0; objective reconstruction'
                if not number(lam) and identity.startswith('ONLINE_') and args.online_prior_lambda is not None:
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
    later_rows = [];b2_rows=[];validation_rows=[];task_rows=[];array_checks=[]
    for stage in ('B2', 'C', 'P1A', 'P1B'):
        paths = sorted((root / 'results' / stage).glob('**/RESULTS.json'))
        valid = []
        for path in paths:
            data = read_json(path, root, sources, errors)
            if data is not None:
                valid.append(dict(source=str(path.relative_to(root)), status=data.get('status', 'SAVED_RECORDS') if isinstance(data, dict) else 'SAVED_RECORDS'))
                if stage=='B2':
                    rows=data if isinstance(data,list) else data.get('rows',[data]) if isinstance(data,dict) else []
                    if isinstance(rows,list) and all(isinstance(r,dict) for r in rows):
                        for r in rows:
                            rec=dict(source=str(path.relative_to(root)),**{k:v for k,v in r.items() if not isinstance(v,(list,dict))})
                            cost=r.get('actual_incremental_cost',{})
                            if isinstance(cost,dict):
                                for key in ('new_charged_complex_rhs','new_charged_factorizations','cpu_seconds','wall_seconds'):
                                    if key in cost:rec['cost_'+key]=cost[key]
                            b2_rows.append(rec)
        later_rows.append(dict(stage=stage, status='SAVED_RECORDS' if valid else 'RUNNING' if stage=='C' and (root/'results/C/phase_diagram').exists() else 'NOT_RUN', saved_result_files=len(valid), sampler_metrics='reported only if explicitly saved; no samples inferred'))
    validation_root=root/'results/B2/validation'
    for path in sorted(validation_root.glob('*/GRID_COMPARISON.json')):
        records=read_json(path,root,sources,errors)
        if isinstance(records,list) and all(isinstance(r,dict) for r in records):
            validation_rows.extend(dict(source=str(path.relative_to(root)),validation_id=path.parent.name,**r) for r in records)
        else:errors.append(dict(path=str(path.relative_to(root)),kind='SchemaError',message='GRID_COMPARISON must be a list of dictionaries'))
    for path in sorted(validation_root.glob('*/ALL_TASKS.json')):
        records=read_json(path,root,sources,errors)
        if isinstance(records,list) and all(isinstance(r,dict) for r in records):
            task_rows.extend(dict(source=str(path.relative_to(root)),validation_id=path.parent.name,**r) for r in records)
        else:errors.append(dict(path=str(path.relative_to(root)),kind='SchemaError',message='ALL_TASKS must be a list of completed task dictionaries'))
    projection = [dict(seed=r.get('seed'), band=r.get('band'), grid='×'.join(map(str, r.get('grid', []))),
                       subsamples=r.get('subsamples'), d_repr=r.get('d_repr'), occupied_area=r.get('occupied_area'),
                       identity=r.get('identity'), minimum_feature=r.get('minimum_feature'), status=r.get('status', 'SAVED_RECORD'))
                  for r in stage_data.get('R1', [])]
    for name, data in [('stage_status', stage_rows), ('R3_runs', runs), ('completed_evaluations', iteration_rows), ('projection_distance', projection), ('later_stages', later_rows),('B2_candidates',b2_rows),('B2_grid_validation',validation_rows),('B2_completed_tasks',task_rows)]:
        write_csv(tables_dir / f'{name}.csv', data)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})
    figures = []
    def finish(fig, name, title, caption):
        fig.tight_layout()
        for ax in fig.axes:
            for collection in ax.collections:collection.set_rasterized(True)
        for suffix in ('png', 'svg'):
            fig.savefig(generated / f'{name}.{suffix}', dpi=170, facecolor='white')
        plt.close(fig)
        figures.append(dict(path=str((generated/f'{name}.png').relative_to(root)), svg=str((generated/f'{name}.svg').relative_to(root)), title=title, caption=caption))
    def placeholder(ax, text):
        ax.text(.5, .5, text, transform=ax.transAxes, ha='center', va='center', wrap=True)
        ax.set_xticks([]); ax.set_yticks([])
    def concise_label(run):
        r=next((r for r in runs if r['run']==run),{})
        if r.get('role')=='CONTROL':return run.replace('controls/','Control ').replace('/',' ')
        if r.get('role')=='FREQUENCY_SUBSET_STAGE':return f"Continuation stage {r.get('stage')}"
        group='Oracle' if r.get('identity')=='ORACLE_DIAGNOSTIC' else 'Online'
        return f"{group} {str(r.get('seed',''))[-4:]} {r.get('band','')} {r.get('method','')}"

    fig, ax = plt.subplots(figsize=(9, 4.5))
    found = False
    for run, hist in histories.items():
        if run.startswith('controls/') or '/stage_' in run:continue
        points = [p for p in hist if 'data_residual_norm' in p]
        if points:
            ax.plot([p['call_id'] for p in points], [p['data_residual_norm'] for p in points], 'o-', ms=4, label=concise_label(run))
            found = True
    if found:
        ax.axhline(1., color='#b45309', ls=':', lw=1, label='d=1 reference line; no gate judgment')
        ax.legend(fontsize=8)
        values = [p['data_residual_norm'] for run,hist in histories.items() if not run.startswith('controls/') and '/stage_' not in run for p in hist if 'data_residual_norm' in p]
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
        if run.startswith('controls/') or '/stage_' in run:continue
        if hist:
            ax.plot([p['call_id'] for p in hist], [p['objective'] for p in hist], 'o-', ms=4, label=concise_label(run))
    if any(hist for run,hist in histories.items() if not run.startswith('controls/') and '/stage_' not in run):
        vals = [p['objective'] for run,hist in histories.items() if not run.startswith('controls/') and '/stage_' not in run for p in hist]
        if min(vals) > 0:
            ax.set_yscale('log');ax.set_ylim(min(vals)*.7,max(vals)*1.4)
        else:
            ax.set_yscale('symlog', linthresh=.1);ax.set_ylim(bottom=0)
        ax.legend(fontsize=8)
    else:
        placeholder(ax, 'NOT_RUN: no completed objective evaluations')
    ax.set_xlabel('Completed evaluator call ID');ax.set_ylabel('Saved objective fun (prior may be included)')
    finish(fig, 'objective_vs_evaluation', '目标函数轨迹', '直接使用已保存 fun；含正则先验时它与纯数据残差不同。')

    fig, ax = plt.subplots(figsize=(10, 6.5))
    points = [r for r in runs if r['saved_status'] == 'SAVED_RESULT' and r['role']=='ALL_BAND_ENDPOINT' and number(r.get('reported_data_distance'))]
    if points:
        colors = ['#c2410c' if r['identity'] == 'ORACLE_DIAGNOSTIC' else '#0369a1' for r in points]
        bars=ax.barh(range(len(points)), [r['reported_data_distance'] for r in points], color=colors)
        ax.set_yticks(range(len(points)), [concise_label(r['run']) for r in points],fontsize=9)
        ax.invert_yaxis()
        for bar,r in zip(bars,points):
            if r.get('method')=='LBFGSB' and r.get('iterations')==0:bar.set_hatch('//')
            ax.text(r['reported_data_distance']+.3,bar.get_y()+bar.get_height()/2,f"{r['reported_data_distance']:.4g}",va='center',fontsize=9)
        ax.set_xlim(right=max(r['reported_data_distance'] for r in points)*1.12)
        ax.set_xlabel('Final full-band saved data distance; hatched = historical interface failure')
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

    subsets=[r for r in runs if r['role']=='FREQUENCY_SUBSET_STAGE']
    if subsets:
        fig,axes=plt.subplots(1,len(subsets),figsize=(4*len(subsets),3.8),squeeze=False)
        for ax,r in zip(axes[0],subsets):
            hist=histories.get(r['run'],[]);points=[p for p in hist if 'data_residual_norm' in p]
            if points:ax.plot([p['call_id'] for p in points],[p['data_residual_norm'] for p in points],'o-',ms=3)
            else:placeholder(ax,'No completed history records')
            ax.set_title(f"stage {r.get('stage')} / {r.get('frequency_count')} frequencies / lambda={r.get('prior_lambda')}",fontsize=9)
            ax.set_xlabel('Evaluator call ID');ax.set_ylabel('Subset data residual')
        finish(fig,'frequency_continuation','分频率续接：每阶段单独展示','每阶段残差对应不同频率子集，不把从1频到5频的距离变化当成同一目标函数的改善或恶化。最后5频距离来自all_wide_d_mean。')

    controls=[r for r in runs if r['role']=='CONTROL']
    if controls:
        fig,ax=plt.subplots(figsize=(9,4.5))
        for r in controls:
            points=[p for p in histories.get(r['run'],[]) if 'data_residual_norm' in p]
            if points:ax.plot([p['call_id'] for p in points],[p['data_residual_norm'] for p in points],'-',label=r['run'])
        ax.set_yscale('symlog',linthresh=1e-7);ax.set_ylim(bottom=0);ax.legend(fontsize=7)
        ax.set_xlabel('Completed evaluator call ID');ax.set_ylabel('Whitened control-data residual')
        finish(fig,'controlled_optimizer_tests','已知几何控制：与未知墙体在线反演分开','vacuum/block的完成结果与disk已保存的部分轨迹分别列出。R_SEALED中的未完成controls不被补成成功，局部驻点也不保证全局拟合。')

    if b2_rows:
        fig,ax=plt.subplots(figsize=(8,4.8))
        for method in sorted({r.get('method','') for r in b2_rows}):
            rows=[r for r in b2_rows if r.get('method')==method and number(r.get('d_cal')) and number(r.get('d_T'))]
            if rows:
                ax.scatter([r['d_cal'] for r in rows],[r['d_T'] for r in rows],label=method)
                for r in rows:ax.annotate(str(r.get('coefficient_progress')),(r['d_cal'],r['d_T']),fontsize=8)
        ax.set_xlabel('Saved calibration distance d_cal');ax.set_ylabel('Saved task distance d_T');ax.legend(fontsize=8)
        ax.set_title('B2 saved oracle-constructed candidates; independent validation pending')
        finish(fig,'B2_cal_vs_task','B2实际候选：校准与任务距离','仅展示B2已保存候选；标注为系数延拓进度。ORACLE_CONSTRUCTED_PHYSICAL_WITNESS不表示在线成功，候选status和未执行项原样保留；不自动给出见证资格。')

    # Exact saved native joint material pair. No decoder, truth generator or PDE import.
    import numpy as np
    joint_folder=root/'results/B2/native_joint_direct_2026102202_narrow'
    joint_final=joint_folder/'FULL_J_FINAL.npz'
    joint_record=next((r for r in b2_rows if r.get('method')=='FULL_J' and r.get('source')==str((joint_folder/'RESULTS.json').relative_to(root))),None)
    if joint_record and joint_final.exists():
        try:
            with np.load(joint_final,allow_pickle=False) as arrays:
                pts=np.asarray(arrays['material_grid']);reference=np.asarray(arrays['reference_material']);candidate=np.asarray(arrays['material'])
            if pts.ndim!=2 or pts.shape[1]!=2 or reference.shape!=(len(pts),) or candidate.shape!=(len(pts),):raise ValueError('invalid saved native material grid/field shape')
            if not all(np.isfinite(a).all() for a in (pts,reference,candidate)):raise ValueError('nonfinite saved material')
            xs,ys=np.unique(pts[:,0]),np.unique(pts[:,1]);ix=np.searchsorted(xs,pts[:,0]);iy=np.searchsorted(ys,pts[:,1])
            if len(xs)*len(ys)!=len(pts):raise ValueError('material_grid is not a complete tensor grid')
            fig,axes=plt.subplots(2,2,figsize=(10,7),sharex=True,sharey=True)
            for row,(world,values) in enumerate((('Reference world',reference),('FULL_J retained world',candidate))):
                image=np.full((len(ys),len(xs)),np.nan+0j);image[iy,ix]=values
                if not np.isfinite(image).all():raise ValueError('material grid has missing/duplicate coordinates')
                for col,(field,lower,upper,name) in enumerate(((image.real,1.,9.,'Real permittivity'),(image.imag,0.,.15,'Conductivity (saved f=1)'))):
                    ax=axes[row,col];im=ax.pcolormesh(xs,ys,field,shading='nearest',vmin=lower,vmax=upper,cmap='viridis',rasterized=True)
                    ax.set_aspect('equal');ax.set_title(f'{world}: {name}',fontsize=10);ax.set_xlabel('x');ax.set_ylabel('y')
                    center=(-.9 if row==0 else -.8,2.15)
                    ax.annotate('original target' if row==0 else 'target translated +0.1',xy=center,xytext=(-1.7,2.65),fontsize=8,color='white',arrowprops=dict(arrowstyle='->',color='white'))
                    fig.colorbar(im,ax=ax,shrink=.8)
            fig.suptitle('ORACLE_EXISTENCE_NOT_ONLINE | target shift 0.1, same radius/contrast/mass',fontsize=11)
            finish(fig,'B2_joint_world_material_gallery','真实joint双世界材料图 · 严格Oracle','读取FULL_J_FINAL.npz中的material_grid/reference_material/material原值，保存频率f=1。两世界Re ε统一[1,9]，σ统一[0,.15]。目标仅平移0.1，质量不变；conditional fixed support，原minimum_feature语义UNRESOLVED，不是在线重建或最终GO。')
            sources.append(dict(path=str(joint_final.relative_to(root)),role='saved native joint world material pair only'))
            array_checks.append(dict(source=str(joint_final.relative_to(root)),identity='ORACLE_EXISTENCE_NOT_ONLINE',saved_frequency=1.,grid_shape=[len(ys),len(xs)],sample_count=len(pts),reference_real_range=[float(reference.real.min()),float(reference.real.max())],candidate_real_range=[float(candidate.real.min()),float(candidate.real.max())],reference_conductivity_range=[float(reference.imag.min()),float(reference.imag.max())],candidate_conductivity_range=[float(candidate.imag.min()),float(candidate.imag.max())],original_minimum_feature='UNRESOLVED',conditional_fixed_support=True))
        except (OSError,ValueError,KeyError) as exc:errors.append(dict(path=str(joint_final.relative_to(root)),kind=type(exc).__name__,message=str(exc)))

    # Independent coarse/fine receiver predictions: fixed saved worlds, no refit.
    for vid in sorted({r['validation_id'] for r in validation_rows}):
        rows=[r for r in validation_rows if r['validation_id']==vid]
        fig,ax=plt.subplots(figsize=(8,4.5));positions=np.arange(len(rows));width=.34
        ax.bar(positions-width/2,[r['coarse_distance'] for r in rows],width,label='coarse saved pair')
        ax.bar(positions+width/2,[r['fine_distance'] for r in rows],width,label='fine saved pair')
        ax.set_xticks(positions,[r['mode'] for r in rows]);ax.set_ylabel('Whitened distance between worlds');ax.legend()
        ax.set_title(f'{vid}: calibration and translated-target joint data',fontsize=10)
        finish(fig,f'B2_grid_distances_{vid}','校准与joint：独立两网格距离','cal为空气查询校准条件；joint8加入原目标/平移目标，是另一观测条件。数值来自GRID_COMPARISON，不混同d_cal与d_joint；不据一项小于1宣布全部物理门通过。')
        pair_path=validation_root/vid/'joint8_PAIR.npz'
        if not pair_path.exists():continue
        try:
            with np.load(pair_path,allow_pickle=False) as arrays:
                co0,co1,fi0,fi1=[np.asarray(arrays[k]) for k in ('coarse0','coarse1','fine0','fine1')]
            if not(co0.shape==co1.shape==fi0.shape==fi1.shape):raise ValueError('coarse/fine joint pair shapes differ')
            if not all(np.isfinite(a).all() for a in (co0,co1,fi0,fi1)):raise ValueError('nonfinite joint pair')
            dc=((co1-co0)/.01).ravel();df=((fi1-fi0)/.01).ravel()
            actual=dict(coarse_distance=float(np.linalg.norm(dc)),fine_distance=float(np.linalg.norm(df)),world0_grid_distance=float(np.linalg.norm(co0-fi0)/.01),world1_grid_distance=float(np.linalg.norm(co1-fi1)/.01))
            claimed=next((r for r in rows if r['mode']=='joint8'),None)
            if claimed:
                for key,value in actual.items():
                    if abs(value-claimed[key])>1e-9*max(1.,abs(value)):raise ValueError(f'saved joint pair disagrees with GRID_COMPARISON: {key}')
            fig,axes=plt.subplots(3,1,figsize=(10,7),sharex=True)
            for ax,part,title in zip(axes[:2],('real','imag'),('Real part of world1 - world0','Imaginary part of world1 - world0')):
                ax.plot(getattr(dc,part),label='coarse',lw=1.2);ax.plot(getattr(df,part),'--',label='fine',lw=1.1)
                ax.set_ylabel('difference / .01');ax.set_title(title,fontsize=10);ax.legend(fontsize=8)
            axes[2].plot(abs(df-dc),color='#a16207');axes[2].set_ylabel('|fine - coarse| / .01');axes[2].set_xlabel('Frequency-major, source-major receiver component')
            fig.suptitle(f'Independent saved joint predictions | d_coarse={actual["coarse_distance"]:.6f}, d_fine={actual["fine_distance"]:.6f}',fontsize=11)
            finish(fig,f'B2_joint_predictions_{vid}','独立coarse/fine实际joint预测差','直接读取joint8_PAIR.npz四组复数接收预测，展示两世界差及其网格变化；从数组复算的距离与GRID_COMPARISON一致。有限两网格一致性不是连续误差包络。')
            sources.append(dict(path=str(pair_path.relative_to(root)),role='independent saved coarse/fine complex predictions'))
            array_checks.append(dict(source=str(pair_path.relative_to(root)),shape=list(co0.shape),distance_replay=actual,GRID_COMPARISON_replay='MATCH',pair_difference_grid_change_white=float(np.linalg.norm(df-dc)),new_physics=0))
        except(OSError,ValueError,KeyError) as exc:errors.append(dict(path=str(pair_path.relative_to(root)),kind=type(exc).__name__,message=str(exc)))

    # Optional completed endpoint material plots; no generator or field cache is read.
    for r in runs:
        folder = r3 / r['run']
        final = folder / 'FINAL.npz'
        if r['saved_status'] != 'SAVED_RESULT' or r['role']!='ALL_BAND_ENDPOINT' or not final.exists():
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

    from w_results_extensions import build_extensions
    extensions=build_extensions(root,sources,errors,finish,table,read_json,atomic,write_csv)
    write_csv(tables_dir/'R3_runs.csv', runs)  # include completed objective/residual checks
    status = 'FINAL_SNAPSHOT' if args.final else 'WORK_IN_PROGRESS'
    now = datetime.now(timezone.utc).isoformat()
    payload = dict(report_status=status, generated_utc=now, physics_actions_by_builder=0,
                   stage_status=stage_rows, R3_runs=runs, projection_distance=projection,
                   later_stages=later_rows,B2_candidates=b2_rows,B2_grid_validation=validation_rows,B2_completed_tasks=task_rows,B2_saved_array_checks=array_checks,R_sealed=sealed,input_sources=sources, parse_errors=errors,
                   frozen_online_prior_lambda=args.online_prior_lambda,C_B2_extensions=extensions['data'],figures=figures,
                   limitations=['Saved-record presentation only; no scientific gate/novelty judgment.',
                                'Syntactically valid iteration snapshots may include unfinished starts; only evaluation_complete is plotted.',
                                'Residual reconstructed as sqrt(2 fun-lambda ||x||^2): recorded objective/x, zero prior center, frozen oracle lambda0 / online lambda.002. CLI override is explicit.',
                                'No samples, posterior weights, ESS, chains or coverage inferred from optimizer endpoints.'])
    atomic(root/'RESULTS_REPORT.json', json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False))
    atomic(root/'REPORT_PARSE_ERRORS.json', json.dumps(dict(status='ERROR' if errors else 'OK', errors=errors), ensure_ascii=False, indent=2))
    atomic(root/'tables/saved_results/B2_saved_array_checks.json',json.dumps(dict(status='SAVED_ARRAY_READBACK_ONLY',physics_actions=0,checks=array_checks),ensure_ascii=False,indent=2,allow_nan=False))
    # This is a data checkpoint, deliberately separate from main-owned gates.
    checkpoint=['# B2实际结果检查点\n','状态：**WORK_IN_PROGRESS / ORACLE_EXISTENCE_NOT_ONLINE**。本次构建只读已完成保存记录，新增Maxwell/LU/RHS为0；没有优化、预算操作或P1执行。\n',
                '## 1. 信息条件与几何限制\n',
                'native分支从已知参考世界构造108参数固定支持chart，属于Oracle物理存在性探索，不是在线重建。参考与FULL_J_FINAL实际材料见[双世界图](figures/saved_results/B2_joint_world_material_gallery.png)。目标沿x平移0.1，圆盘半径/介电对比/质量不变，target_mass_relative_change=0。原minimum_feature标量语义保持**UNRESOLVED**；conditional fixed support不认证全部原几何材料类，也不构成最终GO。\n',
                '## 2. 已保存native候选\n',
                '|来源|方法|d_cal|d_T|direct-joint distance|status|\n|---|---|---:|---:|---:|---|']
    for r in b2_rows:
        if '/native_' in r['source']:
            checkpoint.append('|'+ '|'.join(str(r.get(k,'—')) for k in ('source','method','d_cal','d_T','final_distance','status'))+'|')
    checkpoint+=['\n## 3. 独立coarse/fine完成项\n','|验证ID / 条件|coarse两世界距离|fine两世界距离|world0网格差|world1网格差|\n|---|---:|---:|---:|---:|']
    for r in validation_rows:
        checkpoint.append(f"|{r['validation_id']} / {r['mode']}|{r['coarse_distance']:.12g}|{r['fine_distance']:.12g}|{r['world0_grid_distance']:.12g}|{r['world1_grid_distance']:.12g}|")
    checkpoint+=['\n同一保留材料对在无目标校准条件与平移目标joint条件下具有不同数据距离：不能把joint小于1写成cal小于1。两网格差是经验一致性检查，不是连续PDE误差上界。实际复数预测差见[coarse/fine图](figures/saved_results/B2_joint_predictions_joint_seed2_narrow_fullj.png)，来自joint8_PAIR四数组；没有重新求解或重新拟合。\n','## 4. 任务类进度只按已完成记录公示\n']
    for vid in sorted({r['validation_id'] for r in task_rows}):
        ts=[r for r in task_rows if r['validation_id']==vid];count=len({r.get('probe') for r in ts});distances=[r['distance'] for r in ts if number(r.get('distance'))]
        checkpoint.append(f"* {vid}：当前保存{count}/24个probe结果；已完成子集最大距离{max(distances):.12g}。"+('尚不是完整24类supremum；其余项不补造。' if count<24 else '24项保存齐全；资格解释仍归主线程。'))
    checkpoint+=['\n原始来源：[native cal](results/B2/native_cal_2026102202_narrow/RESULTS.json)、[native joint](results/B2/native_joint_direct_2026102202_narrow/RESULTS.json)、[GRID_COMPARISON](results/B2/validation/joint_seed2_narrow_fullj/GRID_COMPARISON.json)、[ALL_TASKS](results/B2/validation/joint_seed2_narrow_fullj/ALL_TASKS.json)。其他seed/分支不由本seed2检查点推定完成；不据Oracle候选数量产生概率。\n',f'生成UTC：{now}。页面继续WORK_IN_PROGRESS；不修改科学门、物理源码或账本。\n']
    atomic(root/'B2_CHECKPOINT.md','\n'.join(checkpoint))
    figure_html = ''.join(f'<section><h2>{html.escape(f["title"])}</h2><a href="{html.escape(f["svg"])}"><img src="{html.escape(f["path"])}" alt="{html.escape(f["title"])}"></a><p class="note">{html.escape(f["caption"])}</p></section>' for f in figures)
    document = f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Research W · saved results</title>
<style>body{{font:16px/1.65 system-ui,sans-serif;background:#f3f5f7;color:#17212e;margin:0}}main{{max-width:1120px;margin:auto;padding:36px 24px}}h1{{font-size:30px;margin:8px 0}}h2{{font-size:21px}}section{{background:white;border:1px solid #dbe2e9;border-radius:12px;padding:20px;margin:20px 0}}.badge{{display:inline-block;background:#fff0cf;color:#7c4900;padding:4px 10px;border-radius:6px;font-weight:700}}.note,.meta{{color:#506173;font-size:14px}}.scroll{{overflow:auto}}table{{border-collapse:collapse;width:100%;font-size:13px}}td,th{{text-align:left;border-bottom:1px solid #e2e8f0;padding:8px;vertical-align:top;white-space:nowrap}}th{{background:#f6f8fa}}img{{display:block;width:100%;height:auto}}a{{color:#0369a1}}.empty{{color:#6b7280}}.error{{color:#b91c1c}}</style>
<main><span class="badge">{status}</span><h1>Research W · 保存结果快照</h1><p>从 0 读起：本页先展示真值回放与离线表达诊断，再展示在线反演的已保存轨迹。所有图表由现有数据生成；报告构建新增物理动作 0。</p><p class="meta">生成时间 UTC：{now} · <a href="RESULTS_REPORT.json">机器可读报告</a> · <a href="REPORT_PARSE_ERRORS.json">解析错误记录</a></p>
<section><h2>R阶段封账与阅读入口</h2>{table([dict(status=sealed.get('status'),**sealed.get('totals',{}),uncompleted_controls=sealed.get('uncompleted_controls'))]) if sealed else '<p>尚无R_SEALED.json。</p>'}<p><a href="P0_R_REPRESENTABILITY.md">可表达性</a> · <a href="P0_R_CODE_AUDIT.md">代码审计</a> · <a href="P0_R_PROJECTION_ANALYSIS.md">投影分析</a> · <a href="P0_R_FWI_CONVERGENCE.md">FWI收敛记录</a> · <a href="P0_R_ONLINE_VS_ORACLE.md">在线与Oracle</a></p><p class="note">R封账只描述该阶段预算；本页继续WORK_IN_PROGRESS，不作campaign最终裁决。</p></section>
<section><h2>阶段记录状态</h2>{table(stage_rows)}<p class="note">状态来自文件中的记录，不是本页作出的科学门判决。WORK_IN_PROGRESS 与 FINAL_SNAPSHOT 只区分报告展示模式。不存在的后续阶段保持 NOT_RUN。</p></section>
<section><h2>Oracle / 在线结果清单</h2>{table(runs, ['run','role','identity','saved_status','reported_data_distance','data_distance_field','projected_gradient','iterations','success','stop','stopping','completed_evaluations','unfinished_evaluations','objective_residual_reconstruction_status'])}<p class="note">Oracle使用真值信息，仅为离线诊断；controls与未知墙体在线结果分开。projected_gradient不等同数据残差。续接stage按其RESULTS.prior_lambda与实际频率子集恢复残差，最后全带端点使用all_wide_d_mean；不是所有stage都采用λ=.002。<a href="tables/saved_results/R3_runs.csv">CSV</a></p></section>
<section><h2>原始真值回放</h2>{table([r for r in stage_data.get('R0',[]) if 'd_mean' in r], ['seed','band','label','d_mean','d_observed_truth','volume_probe','same_contract_pass','identity'])}<p class="note">d_observed_truth 是保存噪声向量的白化长度，d_mean 是原均值数据的回放距离；两者不混合。</p></section>
{figure_html}<section><h2>投影距离明细 · 严格离线</h2>{table(projection)}<p><a href="tables/saved_results/projection_distance.csv">CSV</a> · <a href="tables/saved_results/completed_evaluations.csv">完整已完成评价表</a></p></section>
<section><h2>后续见证 / 采样阶段</h2>{table(later_rows)}<p class="note">不存在真实结果时，没有 posterior、链、ESS、模态权重或覆盖率图。本页不由优化候选数量构造概率。</p></section>
<section><h2>B2已保存候选明细</h2>{table(b2_rows,['source','seed','band','method','coefficient_progress','d_cal','d_T','final_distance','d_joint_shifted','target_shift','target_mass_relative_change','material_RMS','incremental_rhs','rhs','cost_new_charged_complex_rhs','status','identity','stop','stopping'])}<p><a href="tables/saved_results/B2_candidates.csv">B2 CSV</a>。支持list、rows及单独dict结果；数值保留字段原义。保留原始status，不从小量局部响应或候选个数推断概率与充分覆盖。</p></section>
<section><h2>B2独立验证检查点 · 持续更新</h2><p><a href="B2_CHECKPOINT.md">已完成事实说明</a>。ORACLE_EXISTENCE_NOT_ONLINE；原minimum_feature语义UNRESOLVED，conditional fixed support。目标平移0.1、半径/对比/质量不变；不将joint与无目标cal混为一项，不宣布最终GO或P1执行。</p>{table(validation_rows,['validation_id','mode','coarse_distance','fine_distance','world0_grid_distance','world1_grid_distance'])}<p class="note">ALL_TASKS仅按已保存probe列出，缺项不填；不足24项时不是完整任务类supremum。其他seed/分支进度不由此推定。</p>{table(task_rows,['validation_id','probe','center','radius','contrast','distance'])}<p><a href="tables/saved_results/B2_grid_validation.csv">网格数据CSV</a> · <a href="tables/saved_results/B2_completed_tasks.csv">已完成任务CSV</a></p></section>
{extensions['html']}<section><h2>数据解析与来源</h2><p class="{'error' if errors else 'note'}">{'存在解析错误；相关文件未用于对应图表。' if errors else '当前读取文件均通过 JSON 解析；这不替代科学审阅。'}</p>{table(errors, ['path','kind','message']) if errors else ''}{table(sources, ['path','bytes','modified_utc','role'])}</section></main></html>'''
    for name in ('report.html', 'index.html'):
        atomic(root/name, document)
    print(json.dumps(dict(status='ERROR' if errors else status, completed_evaluations=len(iteration_rows),
                          final_R3_results=sum(r['saved_status']=='SAVED_RESULT' for r in runs),
                          plots=len(figures), parse_errors=len(errors), physics_actions=0,
                          index=str(root/'index.html')), ensure_ascii=False))
    return 2 if errors else 0


if __name__ == '__main__':
    raise SystemExit(main())
