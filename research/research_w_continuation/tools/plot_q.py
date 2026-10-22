#!/usr/bin/env python3
"""Plot saved WQ qualification evidence; never import or run a field solver.

Run under the parent's job wrapper.  Reads only completed JSON rows and lazy
NPZ members.  Missing, partial and failed inputs remain visible in notes; each
figure is independent, so one plot failure does not erase other evidence.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import gc
import json
from pathlib import Path
import re
import traceback
import zipfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
import numpy as np


BLUE = '#2364A5'
ORANGE = '#D47A20'
PURPLE = '#7B4CA0'
GREEN = '#21866B'
GRAY = '#657283'


def finite(value):
    return isinstance(value, (int, float)) and np.isfinite(value)


def encoded(value):
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    if value is None:
        return ''
    return value


def write_csv(path, rows):
    """Retain unrounded numeric values and nested raw fields as JSON cells."""
    rows = list(rows)
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns or ['status'])
        writer.writeheader()
        for row in rows:
            writer.writerow({key: encoded(value) for key, value in row.items()})
    return len(rows)


def styling():
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 10,
        'axes.titlesize': 11, 'axes.labelsize': 10,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.edgecolor': '#8993A0', 'axes.linewidth': .7,
        'grid.color': '#DCE2E8', 'grid.linewidth': .65,
        'legend.frameon': False, 'legend.fontsize': 8.5,
        'savefig.facecolor': 'white', 'figure.facecolor': 'white',
        'svg.fonttype': 'none',
    })


def axis_grid(axis, log=False):
    axis.set_axisbelow(True)
    axis.grid(True, which='major', axis='y', alpha=.8)
    if log:
        axis.set_yscale('log')


def footer(fig, text):
    fig.text(.055, .025, text, ha='left', va='bottom', fontsize=8.2,
             color='#505C6A', linespacing=1.5)


class Inputs:
    def __init__(self, root, output, max_world_points):
        self.root = root
        self.output = output
        self.max_world_points = max_world_points
        self.cache = {}
        self.raw = []
        self.notes = []
        self.sources = []
        self.failures = []
        self.config = {}
        config = root/'docs/P0_FROZEN_CONFIG.json'
        try:
            self.config = json.loads(config.read_text())
        except (OSError, ValueError) as error:
            self.note('protocol', f'Frozen config unavailable: {error}. '
                      'The requested 0.1 line is shown; stored metrics are not rescaled.')
        self.threshold = self.config.get('numerics', {}).get('mesh_whitened_tolerance', .1)

    def relative(self, path):
        try:
            return str(path.relative_to(self.root))
        except ValueError:
            return str(path)

    def note(self, scope, text):
        self.notes.append(dict(scope=scope, note=text))

    def rows(self, directory):
        if directory in self.cache:
            return self.cache[directory]
        path = self.root/'results'/directory/'RESULTS.json'
        source = self.relative(path)
        if not path.exists():
            self.note(directory, f'Missing {source}; unavailable rows are skipped.')
            self.cache[directory] = []
            return []
        try:
            data = json.loads(path.read_text())
            if not isinstance(data, list):
                raise ValueError('expected a JSON list of result rows')
        except (OSError, ValueError) as error:
            self.note(directory, f'Partial/unreadable {source}: {error}; file skipped.')
            self.cache[directory] = []
            return []
        rows = []
        for index, row in enumerate(data):
            if not isinstance(row, dict):
                self.note(directory, f'Non-object result row {index} skipped.')
                continue
            record = dict(row, source_file=source, source_row=index)
            self.raw.append(record)
            rows.append(record)
        self.sources.append(dict(path=source, rows=len(rows),
                                 modified_utc=datetime.fromtimestamp(
                                     path.stat().st_mtime, timezone.utc).isoformat()))
        self.cache[directory] = rows
        return rows

    def read_failure_ledger(self):
        path = self.root/'results/FAILURE_LEDGER.jsonl'
        if path.exists():
            for index, line in enumerate(path.read_text().splitlines()):
                try:
                    row = json.loads(line)
                    self.failures.append(dict(row, source_file=self.relative(path),
                                              source_row=index))
                except ValueError as error:
                    self.note('failure ledger', f'Unreadable row {index}: {error}')
        for job in ('fitted-cylinder-002', 'wall-qualification-001', 'wall-pml-001'):
            path = self.root/'results/jobs'/job/'RECEIPT.json'
            if not path.exists():
                self.note(job, 'No job receipt at this snapshot; completion is not inferred.')
                continue
            try:
                receipt = json.loads(path.read_text())
                if receipt.get('status') != 'COMPLETE' or receipt.get('exit_code') != 0:
                    self.note(job, f"Job status={receipt.get('status')}, "
                              f"exit_code={receipt.get('exit_code')}. Saved successes may be "
                              'plotted individually; later configurations are not invented.')
            except (OSError, ValueError) as error:
                self.note(job, f'Partial/unreadable job receipt: {error}')


def cylinder(ctx):
    attribution = ctx.rows('attribution')
    galerkin = ctx.rows('galerkin-cylinder')
    unfitted = ctx.rows('fem-cylinder_fem-cylinder-002')
    fitted = ctx.rows('fitted-cylinder_fitted-cylinder-002')
    independent = ctx.rows('independent')
    circle_controls = ctx.rows('circle-controls')
    records = []
    fig, axes = plt.subplots(1, 3, figsize=(16.2, 6.3),
                             gridspec_kw={'width_ratios': [1.05, 1.05, 1.7]})
    fig.subplots_adjust(left=.055, right=.97, bottom=.19, top=.82, wspace=.4)
    fig.suptitle('Cylinder accuracy: controlled changes and saved FEM configurations',
                 x=.055, y=.96, ha='left', fontsize=16)

    ax = axes[0]
    controls = []
    for sampling, order in (('center', 4), ('fraction', 4), ('fraction', 8), ('fraction', 16)):
        found = [r for r in attribution if r.get('case') == 'cylinder'
                 and r.get('nx') == 96 and r.get('sampling') == sampling
                 and r.get('q') == order and finite(r.get('relative_error'))]
        if found:
            controls.append(found[-1])
        else:
            ctx.note('01 cylinder', f'Missing nx96 {sampling} q{order}; not plotted.')
    if controls:
        x = np.arange(len(controls))
        values = [100*r['relative_error'] for r in controls]
        bars = ax.bar(x, values, color=[GRAY if r['sampling'] == 'center' else BLUE
                                       for r in controls], width=.68)
        ax.set_xticks(x, [('Center' if r['sampling'] == 'center' else 'Exact area')
                         + f"\nq={r['q']}" for r in controls], fontsize=8.3)
        for bar, value in zip(bars, values):
            ax.annotate(f'{value:.6f}%', (bar.get_x()+bar.get_width()/2, value),
                        xytext=(0, 5), textcoords='offset points', ha='center',
                        fontsize=7.5, rotation=90, va='bottom')
        ax.set_ylim(0, max(values)*1.48)
        records.extend(dict(r, panel='A', plotted_percent=100*r['relative_error'])
                       for r in controls)
    else:
        ax.text(.5, .5, 'No saved nx96 controls', transform=ax.transAxes, ha='center')
    ax.set_title('A  Same nx=96: geometry and kernel order', loc='left')
    ax.set_ylabel('Relative plane-wave Mie error (%)')
    axis_grid(ax)

    ax = axes[1]
    for rows, label, color, marker, linestyle in (
            (attribution, 'VIE: exact fractions, q4', BLUE, 'o', '-'),
            (galerkin, 'VIE: double average, q4', ORANGE, 'D', '--')):
        selected = [r for r in rows if r.get('case') == 'cylinder'
                    and finite(r.get('relative_error')) and finite(r.get('nx'))
                    and (r.get('route') == 'double_average_VIE'
                         or (r.get('sampling') == 'fraction' and r.get('q') == 4))]
        selected.sort(key=lambda r: r['nx'])
        if selected:
            ax.plot([r['nx'] for r in selected], [100*r['relative_error'] for r in selected],
                    color=color, marker=marker, linestyle=linestyle, label=label, lw=1.8)
            records.extend(dict(r, panel='B', series=label,
                                plotted_percent=100*r['relative_error']) for r in selected)
    ax.set_title('B  Constant-cell field variants', loc='left')
    ax.set_xlabel('VIE cells along x (nx)')
    ax.set_xticks([64, 96, 128])
    ax.legend(loc='upper right', fontsize=7.6)
    axis_grid(ax, log=True)

    configs = []
    for rows, kind in ((unfitted, 'Unfitted'), (fitted, 'Fitted')):
        for row in rows:
            if row.get('case') != 'cylinder' or not finite(row.get('relative_error')):
                if row.get('status'):
                    ctx.note('01 cylinder', f"Unplotted saved status: {row.get('status')}; "
                              f"h={row.get('h')}, order={row.get('order')}")
                continue
            label = f"{kind} P{row['order']} | h={row['h']:g}"
            if row.get('polygon_nodes'):
                label += f" | n={row['polygon_nodes']}"
            configs.append((label, row['relative_error'],
                            dict(row, receiver_extraction='FE-point probe'), kind, row['order']))
    for row in independent:
        if row.get('case') == 'cylinder' and finite(row.get('plane_relative')):
            label = f"Fitted P3 | h={row['h']:g} | n={row['n']}"
            configs.append((label, row['plane_relative'],
                            dict(row, receiver_extraction='FE-point probe'), 'Fitted', 3))
    for row in circle_controls:
        if (finite(row.get('plane_volume_relative')) and finite(row.get('h'))
                and finite(row.get('degree')) and finite(row.get('n'))):
            label = (f"Fitted P{row['degree']} volume\n"
                     f"h={row['h']:g} | n={row['n']}")
            configs.append((label, row['plane_volume_relative'],
                            dict(row, receiver_extraction='volume'), 'Fitted', row['degree']))
        else:
            ctx.note('01 cylinder', 'Incomplete circle-control row omitted from panel C; '
                      'its original fields remain in the separate control CSV.')
    ax = axes[2]
    for position, (label, value, row, kind, order) in enumerate(configs):
        color = ORANGE if kind == 'Unfitted' else (GREEN if order == 2 else BLUE)
        marker = ('s' if row.get('receiver_extraction') == 'volume'
                  else 'D' if 'plane_relative' in row else 'o')
        ax.scatter(100*value, position, c=color, s=42, marker=marker, zorder=3)
        ax.annotate(f'{100*value:.4g}%', (100*value, position), xytext=(6, 0),
                    textcoords='offset points', va='center', fontsize=8)
        records.append(dict(row, panel='C', configuration=label,
                            plotted_percent=100*value, plotted_reference='plane-wave Mie'))
    if configs:
        ax.set_yticks(range(len(configs)), [item[0] for item in configs], fontsize=8)
        ax.invert_yaxis()
        ax.set_xscale('log')
        values = [100*item[1] for item in configs]
        ax.set_xlim(min(values)*.65, max(values)*2)
        ax.grid(True, axis='x', which='major')
    else:
        ax.text(.5, .5, 'No saved FEM successes', transform=ax.transAxes, ha='center')
    ax.set_title('C  FEM: measured configurations only', loc='left')
    ax.set_xlabel('Relative plane-wave Mie error (%)')
    footer(fig, 'All panels use saved plane-wave cylinder errors; lower is better. '
           'h is a FEM mesh-construction label; n is interface polygon resolution.\n'
           'The fitted-cylinder job stopped after saved successes. Missing later cases are omitted; '
           'double averaging retains its worse observed errors.\n'
           'Square markers labeled volume are saved circle controls; earlier FEM entries '
           'use FE-point extraction. These are measured configurations, with no fitted convergence law.')
    return fig, records, 'Cylinder controls and saved FEM successes, including negative VIE results.'


def walls(ctx):
    combined = ctx.rows('wall-qualification')+ctx.rows('wall-qualification-resume')
    complete = {}
    for row in combined:
        if row.get('case') == 'wall_complete' and finite(row.get('h')):
            complete[(row.get('seed'), row['h'])] = row
    seeds = sorted({seed for seed, _ in complete})
    if not seeds:
        ctx.note('02 wall', 'No saved wall_complete rows; figure skipped.')
        return None
    fig, axes = plt.subplots(1, len(seeds), figsize=(6.5*len(seeds), 5.5),
                             sharey=True, squeeze=False)
    fig.subplots_adjust(left=.075, right=.98, bottom=.2, top=.8, wspace=.16)
    fig.suptitle('Full-data wall convergence in all measured contexts',
                 x=.075, y=.96, ha='left', fontsize=16)
    records = []
    styles = {'wide_cal': (BLUE, 'o', '-'), 'wide_task': (ORANGE, 'D', '-'),
              'narrow_cal': (BLUE, 's', '--'), 'narrow_task': (ORANGE, '^', '--')}
    all_values = []
    series_count = 0
    for column, seed in enumerate(seeds):
        axis = axes[0, column]
        rows = sorted([row for (s, _), row in complete.items() if s == seed],
                      key=lambda r: -r['h'])
        prior_h = None
        points = {name: [] for name in styles}
        pair_labels = []
        for row in rows:
            metrics = row.get('mesh_white', {})
            if metrics:
                pair = (f'{prior_h:g} → {row["h"]:g}' if prior_h is not None
                        else f'previous unavailable → {row["h"]:g}')
                pair_index = len(pair_labels)
                pair_labels.append(pair)
                for name in styles:
                    if finite(metrics.get(name)):
                        value = metrics[name]
                        points[name].append((pair_index, value))
                        all_values.append(value)
                        records.append(dict(source_file=row['source_file'],
                                            source_row=row['source_row'], seed=seed,
                                            band=name.split('_')[0], label=name.split('_')[1],
                                            previous_h=prior_h, h=row['h'], mesh_white=value,
                                            geometry_tolerance=row.get('tolerance'),
                                            geometry_resolution=row.get('geometry_resolution'),
                                            status=row.get('status')))
                    else:
                        ctx.note('02 wall', f'Missing {name} at seed={seed}, h={row["h"]}.')
            prior_h = row['h']
        for name, values in points.items():
            if values:
                color, marker, line = styles[name]
                axis.plot([v[0] for v in values], [v[1] for v in values],
                          c=color, marker=marker, linestyle=line, lw=1.7,
                          label=name.replace('_', ' / '))
                series_count += 1
        axis.axhline(ctx.threshold, color='#303A45', ls=':', lw=1.4,
                     label=f'Protocol line {ctx.threshold:g}')
        axis.set_title(f'Seed {seed}', loc='left')
        axis.set_xticks(range(len(pair_labels)), pair_labels)
        axis.set_xlabel('Saved adjacent refinement pair (h labels)')
        axis_grid(axis, log=True)
        axis.legend(loc='upper right')
    if not all_values:
        plt.close(fig)
        ctx.note('02 wall', 'Only baseline wall rows exist; no convergence differences to plot.')
        return None
    axes[0, 0].set_ylabel('Full complex-data difference / 0.01')
    axes[0, 0].set_ylim(min(all_values)*.65, max(ctx.threshold*1.6, max(all_values)*1.3))
    if series_count != 8:
        ctx.note('02 wall', f'Partial snapshot: {series_count}/8 expected context curves available.')
    footer(fig, f'{series_count}/8 saved context curves. Each band contains all five frequencies, '
           'four sources and twelve receivers; no averaging over contexts.\n'
           'The refinement changes field mesh and interface geometry together. '
           'These are consistency differences, not continuum PDE error bounds.')
    return fig, records, 'Eight full-data wall context curves with the frozen 0.1 line.'


def kernels(ctx):
    rows = ctx.rows('kernels')
    rows = [r for r in rows if finite(r.get('relative_error')) and finite(r.get('q'))]
    if not rows:
        ctx.note('03 kernel', 'No saved kernel integration rows; figure skipped.')
        return None
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.6))
    fig.subplots_adjust(left=.085, right=.98, bottom=.2, top=.8, wspace=.25)
    fig.suptitle('Source-cell kernel integration: near cells and singular self term',
                 x=.085, y=.96, ha='left', fontsize=15)
    records = []
    for cell, color, marker in (('nearest', BLUE, 'o'), ('diagonal', GREEN, 's'),
                                 ('far', PURPLE, '^'), ('self', ORANGE, 'D')):
        selected = sorted([r for r in rows if r.get('cell') == cell], key=lambda r: r['q'])
        axis = axes[1] if cell == 'self' else axes[0]
        positive = [r for r in selected if r['relative_error'] > 0]
        if positive:
            axis.plot([r['q'] for r in positive], [r['relative_error'] for r in positive],
                      color=color, marker=marker, lw=1.7,
                      ls='--' if cell == 'self' else '-',
                      label='Raw self GL: diagnostic, NOT USED' if cell == 'self' else cell)
        if len(positive) != len(selected):
            ctx.note('03 kernel', f'{cell}: zero errors retained in CSV but omitted on log axis.')
        for row in selected:
            records.append(dict(row, production_role=(
                'raw GL diagnostic; production self uses analytic square integral'
                if cell == 'self' else 'nonself source integration quadrature study'),
                baseline_uses_this_entry=cell != 'self' and row['q'] == 4))
    for axis in axes:
        axis.set_xscale('log', base=2)
        axis.set_xticks([4, 8, 16, 32, 64], ['4', '8', '16', '32', '64'])
        axis.set_xlabel('Tensor Gauss–Legendre order per axis')
        axis.set_ylabel('Relative integrated-kernel difference')
        axis.axvline(4, c=GRAY, ls=':', lw=.9)
        axis_grid(axis, log=True)
        axis.legend(loc='upper right', fontsize=8)
    axes[0].set_title('A  Nonself entries (reference: GL96)', loc='left')
    axes[1].set_title('B  Raw self GL versus analytic square reference', loc='left')
    axes[1].text(.05, .09, 'Production self term: analytic actual-square integral.\n'
                 'The dashed self curve was never used by the old VIE.',
                 transform=axes[1].transAxes, fontsize=8.6,
                 bbox=dict(boxstyle='round,pad=.5', fc='#FFF4E8', ec='none'))
    footer(fig, 'Saved kernel study at h=4/96 and f=1. Near-cell GL96 is a higher-order '
           'numerical integration reference.\nThe raw self-GL error is shown only to explain '
           'the singular quadrature diagnostic; it is not the error of the implemented self term.')
    return fig, records, 'Kernel integration errors with production analytic self explicitly distinguished.'


def world_candidates(ctx):
    candidates = {}
    pattern = re.compile(r'WORLD_(\d+)_h([0-9.]+)\.npz$')
    for priority, directory in enumerate(('wall-qualification', 'wall-qualification-resume')):
        for path in (ctx.root/'results'/directory).glob('WORLD_*.npz'):
            match = pattern.match(path.name)
            if match:
                seed, h = int(match.group(1)), float(match.group(2))
                candidates.setdefault(seed, []).append((h, -priority, path))
    selected = []
    required = {'cal_points', 'cal_epsilon', 'task_points', 'task_epsilon', 'task_field'}
    for seed, entries in sorted(candidates.items()):
        for h, _, path in sorted(entries):
            try:
                with np.load(path, allow_pickle=False) as archive:
                    missing = sorted(required-set(archive.files))
                if missing:
                    ctx.note('04 world', f'{ctx.relative(path)} is partial; missing {missing}; skipped.')
                    continue
                selected.append((seed, h, path))
                break
            except (OSError, ValueError, EOFError, zipfile.BadZipFile) as error:
                ctx.note('04 world', f'Unreadable {ctx.relative(path)}: {error}; skipped.')
    return selected


def sample_world(archive, label, max_points, phase=False):
    points = np.asarray(archive[f'{label}_points'])
    if points.ndim != 2:
        raise ValueError('saved material points must be two-dimensional')
    if points.shape[1] != 2 and points.shape[0] == 2:
        points = points.T
    epsilon = np.asarray(archive[f'{label}_epsilon']).reshape(-1)
    if points.shape != (len(epsilon), 2):
        raise ValueError('saved material points and epsilon lengths differ')
    indices = np.linspace(0, len(points)-1, min(len(points), max_points), dtype=int)
    selected_points = points[indices].copy()
    selected_epsilon = epsilon[indices].copy()
    total_field = None
    if phase:
        values = np.asarray(archive[f'{label}_field'])
        if values.ndim == 1 and len(values) == len(points):
            total_field = values[indices].copy()
        elif values.ndim == 2 and values.shape[0] == len(points):
            total_field = values[indices, 0].copy()
        elif values.ndim == 2 and values.shape[1] == len(points):
            total_field = values[0, indices].copy()
        else:
            raise ValueError('saved total-field shape does not match quadrature points')
        del values
    count = len(points)
    del points, epsilon
    return selected_points, selected_epsilon, total_field, indices, count


def worlds(ctx):
    selected = world_candidates(ctx)
    if not selected:
        ctx.note('04 world', 'No complete saved WORLD quadrature fields; figure skipped.')
        return None
    fig, axes = plt.subplots(len(selected), 3, figsize=(13.4, 4.1*len(selected)),
                             squeeze=False)
    fig.subplots_adjust(left=.055, right=.92, bottom=.14, top=.88, wspace=.18, hspace=.3)
    fig.suptitle('Known material and its saved FEM field — qualification worlds',
                 x=.055, y=.965, ha='left', fontsize=16)
    records = []
    material_artists = []
    phase_artists = []
    maximum_epsilon = 1.
    query = ctx.config.get('query', [[-1.2, 1.2], [1.6, 2.8]])
    for row_index, (seed, h, path) in enumerate(selected):
        ctx.sources.append(dict(path=ctx.relative(path), role='saved f1 quadrature fields',
                                 seed=seed, h=h))
        # One archive, and at most one label's full field member, is in memory.
        with np.load(path, allow_pickle=False) as archive:
            for label, columns in (('cal', (0,)), ('task', (1, 2))):
                points, epsilon, total_field, indices, original_count = sample_world(
                    archive, label, ctx.max_world_points, phase=label == 'task')
                valid = np.isfinite(points).all(axis=1) & np.isfinite(epsilon)
                if total_field is not None:
                    valid &= np.isfinite(total_field)
                if not valid.all():
                    ctx.note('04 world', f'{path.name} {label}: '
                              f'{int((~valid).sum())} nonfinite sampled entries skipped.')
                points, epsilon, indices = points[valid], epsilon[valid], indices[valid]
                if total_field is not None:
                    total_field = total_field[valid]
                if not len(points):
                    raise ValueError(f'{path.name} {label}: no finite saved quadrature samples')
                maximum_epsilon = max(maximum_epsilon, float(epsilon.real.max()))
                for column in columns:
                    axis = axes[row_index, column]
                    phase = column == 2
                    values = np.angle(total_field) if phase else epsilon.real
                    artist = axis.scatter(points[:, 0], points[:, 1], c=values,
                                          cmap='twilight_shifted' if phase else 'viridis',
                                          vmin=-np.pi if phase else 1.,
                                          vmax=np.pi if phase else max(2., maximum_epsilon),
                                          s=3.2, marker='s', linewidths=0, rasterized=True)
                    (phase_artists if phase else material_artists).append(artist)
                    axis.add_patch(Rectangle((query[0][0], query[1][0]),
                                             query[0][1]-query[0][0],
                                             query[1][1]-query[1][0],
                                             fill=False, ec='#D8DDE3' if phase else '#111820',
                                             ls='--', lw=.8))
                    axis.set(xlim=(-2, 2), ylim=(0, 3), aspect='equal',
                             xlabel='x / λ₀', ylabel='y / λ₀')
                    heading = ('Known wall: cal Re ε' if column == 0 else
                               'Known wall + probe: task Re ε' if column == 1 else
                               'Saved task total-field phase, source 1')
                    axis.set_title(f'{heading}\nseed {seed}; h label={h:g}', fontsize=10)
                    if column == 2:
                        axis.text(.02, .025, f'f/f₀=1; {len(points):,}/{original_count:,} q-points',
                                  transform=axis.transAxes, fontsize=7.3, color='#20252B',
                                  bbox=dict(fc='white', ec='none', alpha=.85, pad=2))
                for offset, (xy, eps, index) in enumerate(zip(points, epsilon, indices)):
                    field = total_field[offset] if total_field is not None else None
                    records.append(dict(source_file=ctx.relative(path), seed=seed,
                                        h_label=h, material_label=label, frequency=1.,
                                        saved_qpoint_index=int(index),
                                        saved_qpoint_count=original_count,
                                        sampled_qpoint_count=len(points), x=float(xy[0]), y=float(xy[1]),
                                        epsilon_real=float(eps.real), epsilon_imag=float(eps.imag),
                                        field_source_index=0 if field is not None else None,
                                        total_field_real=float(field.real) if field is not None else None,
                                        total_field_imag=float(field.imag) if field is not None else None,
                                        total_field_phase_rad=float(np.angle(field)) if field is not None else None,
                                        sampling='equally-spaced saved quadrature indices; no interpolation'))
                del points, epsilon, total_field, indices
                gc.collect()
    norm = Normalize(1., maximum_epsilon if maximum_epsilon > 1 else 2.)
    for artist in material_artists:
        artist.set_norm(norm)
    material_bar = fig.add_axes([.935, .53, .012, .28])
    fig.colorbar(material_artists[0], cax=material_bar, label='Known Re ε')
    phase_bar = fig.add_axes([.935, .18, .012, .25])
    bar = fig.colorbar(phase_artists[0], cax=phase_bar, label='arg E total (rad)')
    bar.set_ticks([-np.pi, -np.pi/2, 0., np.pi/2, np.pi],
                  labels=['−π', '−π/2', '0', 'π/2', 'π'])
    footer(fig, 'These are the known qualification materials and already-solved f=1 FEM total fields, '
           'not reconstructions. Dashed rectangle: fixed query region.\n'
           f'Display uses at most {ctx.max_world_points:,} saved quadrature points per label, '
           'deterministically subsampled; no nodal approximation or phase interpolation.')
    return fig, records, 'Known cal/task materials and saved task total-field phase; bounded quadrature scatter.'


def pml(ctx):
    rows = ctx.rows('wall-pml')
    bands = [r for r in rows if r.get('case') == 'wall_PML_complete']
    thickness = [r for r in rows if r.get('case') == 'wall_PML_thickness'
                 and finite(r.get('white_difference'))]
    if not bands and not thickness:
        ctx.note('05 PML', 'No completed PML band/thickness summaries; figure skipped.')
        return None
    order = {'wide_cal': 0, 'wide_task': 1, 'narrow_cal': 2, 'narrow_task': 3}
    bands.sort(key=lambda r: (r.get('seed'), order.get(f"{r.get('band')}_{r.get('label')}", 99)))
    fig, axes = plt.subplots(1, 2, figsize=(14.4, 5.8),
                             gridspec_kw={'width_ratios': [1.4, 1]})
    fig.subplots_adjust(left=.07, right=.98, bottom=.23, top=.8, wspace=.22)
    fig.suptitle('Saved PML and receiver-extraction comparisons',
                 x=.07, y=.96, ha='left', fontsize=16)
    records = []
    values = []
    axis = axes[0]
    series = [('PML_strength_white', 'PML strength: 4 → 6 (probe)', BLUE),
              ('volume_probe_white', 'Volume vs FE-point extraction', ORANGE),
              ('combined_white', 'Volume @ strength 6 vs baseline', PURPLE)]
    x = np.arange(len(bands))
    for index, (key, label, color) in enumerate(series):
        xs, ys = [], []
        for position, row in enumerate(bands):
            if finite(row.get(key)) and row[key] > 0:
                xs.append(position+(index-1)*.25)
                ys.append(row[key])
                values.append(row[key])
            else:
                ctx.note('05 PML', f'Missing/nonpositive {key} for '
                          f"seed={row.get('seed')} {row.get('band')} {row.get('label')}.")
        axis.bar(xs, ys, width=.23, color=color, label=label)
    axis.set_xticks(x, [f"S{str(r['seed'])[-2:]}\n{r['band']} {r['label']}" for r in bands], fontsize=8)
    axis.set_title(f'A  Full-band comparisons ({len(bands)}/8 saved groups)', loc='left')
    axis.set_ylabel('Full complex-data difference / 0.01')
    axis.legend(loc='upper left', fontsize=8)
    records.extend(dict(r, panel='A') for r in bands)
    axis_grid(axis, log=True)
    axis = axes[1]
    seeds = sorted({r['seed'] for r in thickness})
    for seed_index, seed in enumerate(seeds):
        for label, marker, line in (('cal', 'o', '-'), ('task', 'D', '--')):
            selected = sorted([r for r in thickness if r['seed'] == seed and r['label'] == label],
                              key=lambda r: r['frequency'])
            if selected:
                color = BLUE if seed_index == 0 else ORANGE
                axis.plot([r['frequency'] for r in selected], [r['white_difference'] for r in selected],
                          c=color, marker=marker, ls=line, lw=1.7,
                          label=f'S{str(seed)[-2:]} / {label}')
                values.extend(r['white_difference'] for r in selected if r['white_difference'] > 0)
                records.extend(dict(r, panel='B') for r in selected)
    axis.set_title('B  Thickness 1.0 → 1.3, remeshed exterior', loc='left')
    axis.set_xlabel('Frequency f/f₀ (saved edge frequencies)')
    axis.set_xticks(sorted({r['frequency'] for r in thickness}))
    axis.legend(loc='upper left')
    axis_grid(axis, log=True)
    for axis in axes:
        axis.axhline(ctx.threshold, c='#303A45', ls=':', lw=1.3)
        axis.text(.985, ctx.threshold, f' {ctx.threshold:g} protocol line',
                  transform=axis.get_yaxis_transform(), ha='right', va='bottom', fontsize=8)
        if values:
            axis.set_ylim(min(values)*.55, max(ctx.threshold*1.55, max(values)*1.4))
    footer(fig, 'Panel A uses h label=0.145 and complete five-frequency bands. '
           'Panel B uses h label=0.18; the expanded exterior is remeshed.\n'
           'Thickness differences therefore combine exterior/PML and mesh changes; '
           'the three panel-A metrics are separate differences, not additive error bounds.')
    if len(bands) != 8:
        ctx.note('05 PML', f'Partial snapshot: {len(bands)}/8 complete band/state groups.')
    return fig, records, 'Full-band PML strength/extraction differences and remeshed exterior thickness checks.'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--outdir', type=Path, help='default: ROOT/figures')
    parser.add_argument('--max-world-points', type=int, default=12000)
    parser.add_argument('--dpi', type=int, default=190)
    args = parser.parse_args()
    if not 1 <= args.max_world_points <= 50000:
        parser.error('--max-world-points must be between 1 and 50000')
    if not 72 <= args.dpi <= 400:
        parser.error('--dpi must be between 72 and 400')
    root = args.root.resolve()
    output = (args.outdir or root/'figures').resolve()
    output.mkdir(parents=True, exist_ok=True)
    styling()
    ctx = Inputs(root, output, args.max_world_points)
    ctx.read_failure_ledger()
    manifest = dict(generated_utc=datetime.now(timezone.utc).isoformat(),
                    root=str(root), output=str(output), figures=[],
                    physical_solves=0, physical_rhs=0, physical_lu=0,
                    world_max_display_points_per_label=args.max_world_points,
                    world_scatter_rasterized_in_svg=True)
    builders = [('01_cylinder_accuracy', cylinder), ('02_wall_convergence', walls),
                ('03_kernel_quadrature', kernels), ('04_known_world_and_field', worlds),
                ('05_pml_extraction', pml)]
    failures = 0
    for stem, build in builders:
        item = dict(stem=stem, status='SKIPPED', outputs=[])
        try:
            result = build(ctx)
            if result is None:
                manifest['figures'].append(item)
                continue
            fig, rows, caption = result
            data_path = output/f'{stem}.csv'
            item.update(caption=caption, csv_rows=write_csv(data_path, rows))
            item['outputs'].append(data_path.name)
            if stem == '01_cylinder_accuracy':
                control_path = output/'01_circle_controls.csv'
                item['circle_control_csv_rows'] = write_csv(control_path, ctx.rows('circle-controls'))
                item['outputs'].append(control_path.name)
            for extension in ('png', 'svg'):
                path = output/f'{stem}.{extension}'
                fig.savefig(path, dpi=args.dpi, bbox_inches='tight')
                item['outputs'].append(path.name)
            item['status'] = 'GENERATED_NOT_YET_VISUALLY_REVIEWED'
            plt.close(fig)
            print(f'{stem}: generated PNG, SVG and {item["csv_rows"]} CSV rows', flush=True)
        except Exception as error:
            failures += 1
            item.update(status='FAILED', error=str(error))
            record = dict(utc=datetime.now(timezone.utc).isoformat(), figure=stem,
                          error_type=type(error).__name__, error=str(error),
                          traceback=traceback.format_exc())
            with (output/'PLOT_FAILURES.jsonl').open('a') as stream:
                stream.write(json.dumps(record, ensure_ascii=False)+'\n')
            ctx.note(stem, f'Plot failed and is preserved in PLOT_FAILURES.jsonl: {error}')
            plt.close('all')
            print(f'{stem}: FAILED: {type(error).__name__}: {error}', flush=True)
        manifest['figures'].append(item)
        gc.collect()
    write_csv(output/'raw_result_rows.csv', ctx.raw)
    write_csv(output/'raw_run_failures.csv', ctx.failures)
    generated_count = sum(item['status'] == 'GENERATED_NOT_YET_VISUALLY_REVIEWED'
                          for item in manifest['figures'])
    skipped_count = sum(item['status'] == 'SKIPPED' for item in manifest['figures'])
    status = ('PARTIAL_PLOT_FAILURE' if failures else
              'NO_FIGURES_GENERATED' if not generated_count else
              'PARTIAL_INPUTS_PENDING_VISUAL_REVIEW' if skipped_count else
              'GENERATED_PENDING_VISUAL_REVIEW')
    manifest.update(sources=ctx.sources, notes=ctx.notes,
                    input_failed_runs=len(ctx.failures), plot_failures=failures,
                    generated_figures=generated_count, skipped_figures=skipped_count,
                    raw_result_rows=len(ctx.raw),
                    status=status)
    (output/'MANIFEST.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False,
                                                 allow_nan=False)+'\n')
    notes = ['# WQ figures: source and completeness notes', '',
             f"Generated: {manifest['generated_utc']}", '',
             'All figures read saved arrays/rows. No Maxwell, LU or RHS action is performed. '
             'Generation does not certify scientific qualification or visual acceptance.', '',
             'Only files listed as generated in the current MANIFEST belong to this snapshot. '
             'Old files can remain after a skipped/failed rebuild.', '',
             '## Outputs', '']
    notes.extend(f"- {item['stem']}: {item['status']}" for item in manifest['figures'])
    notes.extend(['', '## Input failures preserved', '',
                  'raw_run_failures.csv retains the actual failed-run ledger. '
                  'A failed job may still have individually saved successful rows.', '',
                  '## Snapshot notes', ''])
    notes.extend(f"- [{row['scope']}] {row['note']}" for row in ctx.notes)
    if not ctx.notes:
        notes.append('- No missing or partial input was observed in this snapshot.')
    (output/'PLOT_NOTES.md').write_text('\n'.join(notes)+'\n')
    print(f'Manifest: {output/"MANIFEST.json"}', flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
