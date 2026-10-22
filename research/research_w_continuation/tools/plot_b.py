#!/usr/bin/env python3
"""Plot corrected, saved P0-B evidence. No search or Maxwell action is run.

Only corrected online calibration distances select displayed candidates.
Continuous truth is accessed in the offline image routine, after selection.
Standalone RHS bars are offline attributions, never actual-ledger totals.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import gc
import json
from pathlib import Path
import traceback

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Patch, Rectangle
import numpy as np

from plot_q import (BLUE, ORANGE, PURPLE, GREEN, GRAY, finite, styling,
                    write_csv, footer)


METHODS = ('FULL_J', 'RANDOM_MATCHED', 'MULTISTART_FWI', 'TASK_AWARE')
METHOD_LABELS = ('Full J', 'Random\nmatched', 'Multistart\nFWI', 'Task\naware')
METHOD_COLORS = (BLUE, GREEN, PURPLE, ORANGE)
DISPLAY_NX, DISPLAY_NY = 256, 192
EXCLUSION = ('Earlier results/p0b/search (p0b-search-001) is invalid for method '
             'comparison and excluded: its corrector continued from a rejected '
             'proposal instead of the accepted best/no-op anchor. Its paid costs '
             'remain in the original ledger; these plots do not reassign them.')


def positive(value):
    return finite(value) and value > 0


def context_name(context):
    seed, band = context
    return f'Seed {seed} / {band}'


class Evidence:
    def __init__(self, root, output):
        self.root, self.output = root, output
        self.notes = [dict(scope='excluded run', note=EXCLUSION)]
        self.sources = []
        self.raw = {}
        self.config = {}
        try:
            self.config = json.loads((root/'docs/P0_FROZEN_CONFIG.json').read_text())
        except (OSError, ValueError) as error:
            self.note('config', f'Frozen config unavailable: {error}; fixed display domain '
                      'and requested protocol lines are retained, metrics are not rescaled.')
        self.scored = self.load('results/p0b/SCORED_RESULTS.json')
        self.search = self.load('results/p0b/search-corrected/RESULTS.json')
        self.endpoints = self.load('results/p0b/task-endpoints/RESULTS.json')
        self.scored_by_key = {self.key(r): r for r in self.scored if self.key(r) is not None}
        self.endpoint_by_key = {self.key(r): r for r in self.endpoints if self.key(r) is not None}
        self.candidates = {}
        self.initial = {}
        self.not_run = []
        for row in self.search:
            if row.get('case') == 'calibration_candidate' and self.key(row) is not None:
                if positive(row.get('d_cal')) and finite(row.get('method_rhs')):
                    self.candidates[self.key(row)] = row
                else:
                    self.note('corrected search', f'Incomplete candidate row {row["source_row"]} skipped.')
            elif row.get('case') == 'context_initial':
                self.initial[(row.get('seed'), row.get('band'))] = row
            elif row.get('case') == 'context_not_run':
                self.not_run.append(row)
        self.contexts = sorted({key[:2] for key in self.candidates},
                               key=lambda c: (c[0], 0 if c[1] == 'wide' else 1))
        for key, row in self.candidates.items():
            score = self.scored_by_key.get(key)
            if score is not None and score.get('d_cal') != row['d_cal']:
                self.note('cross-check', f'{key}: scored and corrected calibration values differ; '
                          'online plots/selection use the corrected online record.')
        for row in self.not_run:
            self.note('NOT_RUN', f"{context_name((row.get('seed'), row.get('band')))}: "
                      f"{row.get('reason')}; remaining_rhs={row.get('remaining_rhs')}. "
                      'No missing method results are filled in.')

    @staticmethod
    def key(row):
        if row.get('method') not in METHODS or not finite(row.get('seed')):
            return None
        if row.get('band') not in ('wide', 'narrow'):
            return None
        return row['seed'], row['band'], row['method']

    def relative(self, path):
        try:
            return str(path.relative_to(self.root))
        except ValueError:
            return str(path)

    def note(self, scope, message):
        self.notes.append(dict(scope=scope, note=message))

    def load(self, relative):
        path = self.root/relative
        try:
            data = json.loads(path.read_text())
            if not isinstance(data, list):
                raise ValueError('expected a JSON row list')
        except (OSError, ValueError) as error:
            self.note(relative, f'Missing/partial input skipped: {error}')
            self.raw[relative] = []
            return []
        rows = []
        for index, row in enumerate(data):
            if isinstance(row, dict):
                rows.append(dict(row, source_file=relative, source_row=index))
            else:
                self.note(relative, f'Non-object row {index} skipped.')
        self.raw[relative] = rows
        self.sources.append(dict(path=relative, rows=len(rows),
                                 modified_utc=datetime.fromtimestamp(
                                     path.stat().st_mtime, timezone.utc).isoformat()))
        return rows

    def online_rows(self, context):
        return [self.candidates[context+(method,)] for method in METHODS
                if context+(method,) in self.candidates]

    def best_online(self, context):
        rows = self.online_rows(context)
        if not rows:
            return None, []
        chosen = min(rows, key=lambda r: (r['d_cal'], r['method']))
        tied = [r['method'] for r in rows if r['d_cal'] == chosen['d_cal']]
        return chosen, tied

    def not_run_caption(self):
        if not self.not_run:
            return 'Only recorded contexts are displayed; no missing context is imputed.'
        allowance = self.config.get('search', {}).get('per_context_rhs', 800)
        return '; '.join(f"seed {r['seed']} / {r['band']}: NOT_RUN "
                         f"({r.get('remaining_rhs')} RHS remained; full context allowance {allowance})"
                         for r in self.not_run)


def calibration(ctx, stem):
    if not ctx.contexts:
        ctx.note(stem, 'No corrected calibration candidates; figure skipped.')
        return None
    columns = len(ctx.contexts)
    fig, axes = plt.subplots(2, columns, figsize=(4.6*columns, 7.5), squeeze=False,
                             gridspec_kw={'height_ratios': [1.2, 1.]})
    fig.subplots_adjust(left=.065, right=.98, bottom=.21, top=.86, hspace=.36, wspace=.22)
    fig.suptitle('P0-B: calibration candidates and standalone RHS attribution',
                 x=.065, y=.97, ha='left', fontsize=16)
    records = []
    all_distances = [r['d_cal'] for r in ctx.candidates.values()]
    upper = max(70., max(all_distances)*1.6)
    for column, context in enumerate(ctx.contexts):
        initial = ctx.initial.get(context, {})
        top, bottom = axes[0, column], axes[1, column]
        attributed = []
        for index, method in enumerate(METHODS):
            row = ctx.candidates.get(context+(method,))
            if row is None:
                ctx.note(stem, f'{context_name(context)} {method}: missing; no bar inserted.')
                continue
            distance = row['d_cal']
            top.bar(index, distance, color=METHOD_COLORS[index], width=.65, zorder=2)
            top.annotate(f'{distance:.3f}', (index, distance), xytext=(0, 5),
                         textcoords='offset points', ha='center', fontsize=8)
            setup = 160 if method == 'TASK_AWARE' else 80
            standalone = row['method_rhs']+setup
            attributed.append(standalone)
            bottom.bar(index, row['method_rhs'], color=METHOD_COLORS[index], width=.65)
            bottom.bar(index, setup, bottom=row['method_rhs'], color='#E2E6EB',
                       edgecolor='#99A3AF', hatch='///', lw=.5, width=.65)
            bottom.annotate(str(standalone), (index, standalone), xytext=(0, 4),
                            textcoords='offset points', ha='center', fontsize=8)
            score = ctx.scored_by_key.get(context+(method,), {})
            if (finite(score.get('standalone_proposal_rhs'))
                    and score['standalone_proposal_rhs'] != standalone):
                ctx.note(stem, f'{context_name(context)} {method}: saved standalone attribution '
                          'differs from method_rhs+80/+160; both values remain in CSV.')
            records.append(dict(row, online_calibration_distance=distance,
                                actual_method_rhs=row['method_rhs'],
                                offline_standalone_setup_rhs=setup,
                                offline_standalone_rhs=standalone,
                                scored_standalone_proposal_rhs=score.get('standalone_proposal_rhs'),
                                actual_shared_setup_rhs_for_context=initial.get('shared_rhs'),
                                rhs_interpretation='offline standalone attribution; not a ledger sum'))
        top.axhline(1., color='#303A45', ls=':', lw=1.5, label='Calibration line: 1')
        if positive(initial.get('d_cal')):
            top.axhline(initial['d_cal'], color=GRAY, ls='--', lw=1.0, label='Saved BP start')
        top.set_yscale('log')
        top.set_ylim(.65, upper)
        top.set_title(context_name(context), loc='left', fontsize=10.5)
        top.set_ylabel('Calibration distance d_cal')
        top.legend(loc='lower left', fontsize=7.7)
        top.grid(axis='y', which='major', alpha=.5)
        bottom.set_ylim(0., max(attributed or [1])*1.3)
        bottom.set_ylabel('Standalone proposal RHS (attributed)')
        paid = initial.get('shared_rhs')
        bottom.set_title(f'Actual shared setup paid once: {paid} RHS' if paid is not None
                         else 'Actual shared setup unavailable', loc='left', fontsize=9)
        bottom.grid(axis='y', alpha=.5)
        for axis in (top, bottom):
            axis.set_axisbelow(True)
            axis.set_xticks(range(4), METHOD_LABELS, fontsize=8)
        if len(ctx.online_rows(context)) != 4:
            ctx.note(stem, f'{context_name(context)} is partial: '
                      f'{len(ctx.online_rows(context))}/4 methods recorded.')
    axes[1, 0].legend(handles=[Patch(facecolor=BLUE, label='Method RHS: recorded'),
                              Patch(facecolor='#E2E6EB', edgecolor='#99A3AF', hatch='///',
                                    label='Standalone setup: offline attribution')],
                       loc='upper left', fontsize=7)
    footer(fig, f'No observed d_cal ≤ 1; minimum saved online d_cal={min(all_distances):.6f}. '
           'Standalone = method_rhs +80 for baselines / +160 for TASK_AWARE.\n'
           'Shared work was actually paid once per context; standalone bars must not be summed '
           'as actual ledger usage or interpreted as a speedup.\n'+ctx.not_run_caption())
    count = write_csv(ctx.output/f'{stem}.csv', records)
    return fig, count, 'Corrected online calibration distances and separately labeled offline standalone RHS attribution.', []


def offline_images(ctx, stem):
    if not ctx.contexts:
        ctx.note(stem, 'No corrected online candidates; offline image figure skipped.')
        return None
    path = ctx.root/'results/p0b/search-corrected/ARRAYS.npz'
    if not path.exists():
        ctx.note(stem, 'Corrected theta archive missing; no candidate image manufactured.')
        return None
    # These helpers only evaluate existing coefficient coordinates. No FEM
    # constructor or solve is called; truth is confined to this offline image.
    from wq.inverse_support import PixelGrid, pixel_material
    from p0.geometry import wall_values, in_query
    grid = PixelGrid()
    x = -2+(np.arange(DISPLAY_NX)+.5)*4/DISPLAY_NX
    y = (np.arange(DISPLAY_NY)+.5)*3/DISPLAY_NY
    xx, yy = np.meshgrid(x, y)
    points = np.c_[xx.ravel(), yy.ravel()]
    fig, axes = plt.subplots(len(ctx.contexts), 3, figsize=(12.8, 3.05*len(ctx.contexts)),
                             squeeze=False)
    fig.subplots_adjust(left=.065, right=.9, bottom=.14, top=.9, wspace=.17, hspace=.38)
    fig.suptitle('Offline image inspection: known wall, saved BP and online-selected candidate',
                 x=.065, y=.975, ha='left', fontsize=15)
    query = ctx.config.get('query', [[-1.2, 1.2], [1.6, 2.8]])
    norm = Normalize(0., 8.)
    artists = []
    selected_rows = []
    display_path = ctx.output/f'{stem}.csv'
    theta_path = ctx.output/'07_p0b_selected_theta.csv'
    display_columns = ('seed', 'band', 'selected_method', 'selected_online_d_cal',
                       'selection_tied_methods', 'display_x', 'display_y',
                       'offline_known_contrast', 'saved_BP_contrast', 'selected_contrast',
                       'offline_known_conductivity', 'saved_BP_conductivity', 'selected_conductivity',
                       'theta_source_file', 'theta_candidate_key', 'theta_BP_key',
                       'reference_role')
    theta_columns = ('seed', 'band', 'role', 'method', 'array_key', 'block', 'iy', 'ix',
                     'theta_value', 'source_file')
    count = theta_count = 0
    with np.load(path, allow_pickle=False) as archive, display_path.open('w', newline='') as stream, \
            theta_path.open('w', newline='') as theta_stream:
        writer = csv.DictWriter(stream, fieldnames=display_columns)
        writer.writeheader()
        theta_writer = csv.DictWriter(theta_stream, fieldnames=theta_columns)
        theta_writer.writeheader()
        for context_index, context in enumerate(ctx.contexts):
            chosen, tied = ctx.best_online(context)
            seed, band = context
            key = f'{seed}_{band}_{chosen["method"]}'
            common_key = f'{seed}_{band}_COMMON'
            if key not in archive.files or common_key not in archive.files:
                ctx.note(stem, f'{context_name(context)}: chosen/common theta absent; row skipped.')
                for axis in axes[context_index]:
                    axis.text(.5, .5, 'Saved theta unavailable', transform=axis.transAxes, ha='center')
                    axis.set_axis_off()
                continue
            theta = np.asarray(archive[key]).copy()
            common = np.asarray(archive[common_key]).copy()
            candidate = pixel_material(theta, 1., grid)(points)
            bp = pixel_material(common, 1., grid)(points)
            er, sigma = wall_values(points, seed)
            query_mask = in_query(points)
            known = np.where(query_mask, 0., er-1)
            known_sigma = np.where(query_mask, 0., sigma)
            initial = ctx.initial.get(context, {})
            tie_note = f' ({len(tied)}-way tie)' if len(tied) > 1 else ''
            labels = ('Offline known wall (image reference only)',
                      f'Saved BP start; d_cal={initial.get("d_cal", float("nan")):.3f}',
                      f'{chosen["method"]}{tie_note}; d_cal={chosen["d_cal"]:.3f}')
            fields = (known, bp.real-1, candidate.real-1)
            for column, (label, field) in enumerate(zip(labels, fields)):
                axis = axes[context_index, column]
                artists.append(axis.imshow(field.reshape(DISPLAY_NY, DISPLAY_NX),
                                           extent=(-2, 2, 0, 3), origin='lower',
                                           cmap='viridis', norm=norm, interpolation='nearest',
                                           aspect='equal', rasterized=True))
                axis.add_patch(Rectangle((query[0][0], query[1][0]),
                                         query[0][1]-query[0][0], query[1][1]-query[1][0],
                                         fill=False, ec='#E2E6EB', ls='--', lw=.8))
                axis.set_title(label, fontsize=9)
                axis.set_xlabel('x / λ₀')
                axis.set_ylabel(f'{context_name(context)}\ny / λ₀' if column == 0 else 'y / λ₀',
                                fontsize=9)
            for index, xy in enumerate(points):
                writer.writerow(dict(seed=seed, band=band, selected_method=chosen['method'],
                                     selected_online_d_cal=chosen['d_cal'],
                                     selection_tied_methods=';'.join(tied),
                                     display_x=float(xy[0]), display_y=float(xy[1]),
                                     offline_known_contrast=float(known[index]),
                                     saved_BP_contrast=float(bp[index].real-1),
                                     selected_contrast=float(candidate[index].real-1),
                                     offline_known_conductivity=float(known_sigma[index]),
                                     saved_BP_conductivity=float(bp[index].imag),
                                     selected_conductivity=float(candidate[index].imag),
                                     theta_source_file=ctx.relative(path), theta_candidate_key=key,
                                     theta_BP_key=common_key,
                                     reference_role='OFFLINE_IMAGE_ONLY_NOT_SELECTION_OR_RESCORING'))
                count += 1
            for role, method, array_key, values in (
                    ('BP', 'COMMON', common_key, common),
                    ('selected', chosen['method'], key, theta)):
                for block_index, block in enumerate(('contrast', 'conductivity')):
                    for pixel in range(grid.n_pixels):
                        iy, ix = divmod(pixel, grid.nx)
                        theta_writer.writerow(dict(seed=seed, band=band, role=role, method=method,
                                                   array_key=array_key, block=block, iy=iy, ix=ix,
                                                   theta_value=float(values[block_index*grid.n_pixels+pixel]),
                                                   source_file=ctx.relative(path)))
                        theta_count += 1
            selected_rows.append(dict(seed=seed, band=band, selected_method=chosen['method'],
                                      online_d_cal=chosen['d_cal'], tied_methods=tied,
                                      selection='minimum corrected online d_cal; lexical method-name tie break',
                                      theta_key=key, common_theta_key=common_key,
                                      truth_use='offline image only; no metric recomputed'))
            del theta, common, candidate, bp, er, sigma, known, known_sigma, fields
            gc.collect()
    if not artists:
        plt.close(fig)
        ctx.note(stem, 'No selected theta arrays could be displayed; figure skipped.')
        return None
    bar_axis = fig.add_axes([.925, .25, .014, .48])
    fig.colorbar(artists[0], cax=bar_axis, label='Dielectric contrast Re ε − 1')
    footer(fig, 'Candidates are chosen solely by saved corrected online d_cal; exact ties use '
           'lexical method names. Continuous known truth enters only the offline first column.\n'
           'BP and candidate use the saved full 48×32 contrast/conductivity grids; display sampling '
           'is 256×192, with calibration-query air imposed pointwise. No reconstruction-success claim.')
    selected_path = ctx.output/'07_p0b_image_selection.csv'
    write_csv(selected_path, selected_rows)
    ctx.sources.append(dict(path=ctx.relative(path), role='saved full-grid theta only',
                            display_grid=[DISPLAY_NX, DISPLAY_NY],
                            inverse_grid=[grid.nx, grid.ny]))
    ctx.note(stem, f'Exported {theta_count} original theta entries separately. '
              'The coarse display does not recompute the original offline material scores.')
    return fig, count, 'Offline continuous known wall versus saved BP and candidate selected by online d_cal.', \
        [theta_path.name, selected_path.name]


def export_probe_arrays(ctx):
    calibration_path = ctx.root/'results/p0b/search-corrected/ARRAYS.npz'
    task_path = ctx.root/'results/p0b/task-endpoints/ARRAYS.npz'
    reference_path = ctx.root/'results/wall-pml/ARRAYS.npz'
    missing = [ctx.relative(p) for p in (calibration_path, task_path, reference_path) if not p.exists()]
    if missing:
        ctx.note('08 raw arrays', f'Missing archives {missing}; complex-array CSV skipped.')
        return []
    path = ctx.output/'08_p0b_fixed_probe8_complex_data.csv'
    columns = ('seed', 'band', 'method', 'frequency_index', 'frequency', 'source_index',
               'receiver_index', 'candidate_cal_real', 'candidate_cal_imag',
               'candidate_task_real', 'candidate_task_imag', 'reference_cal_real',
               'reference_cal_imag', 'reference_task_real', 'reference_task_imag',
               'candidate_cal_key', 'candidate_task_key', 'reference_cal_key', 'reference_task_key')
    count = 0
    with np.load(calibration_path, allow_pickle=False) as cal, \
            np.load(task_path, allow_pickle=False) as task, \
            np.load(reference_path, allow_pickle=False) as ref, path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for context in ctx.contexts:
            seed, band = context
            for method in METHODS:
                key = f'{seed}_{band}_{method}'
                keys = (key+'_data', key+'_task_data',
                        f'seed{seed}_{band}_cal_volume', f'seed{seed}_{band}_task_volume')
                if keys[0] not in cal.files or keys[1] not in task.files \
                        or keys[2] not in ref.files or keys[3] not in ref.files:
                    ctx.note('08 raw arrays', f'{context_name(context)} {method}: missing array; skipped.')
                    continue
                arrays = (cal[keys[0]], task[keys[1]], ref[keys[2]], ref[keys[3]])
                shape = arrays[0].shape
                if len(shape) != 3 or any(a.shape != shape for a in arrays):
                    ctx.note('08 raw arrays', f'{key}: shape mismatch; no reshaping or fabricated data.')
                    continue
                frequencies = ctx.config.get('frequencies', {}).get(band, [])
                for fi, si, ri in np.ndindex(shape):
                    values = [a[fi, si, ri] for a in arrays]
                    writer.writerow(dict(seed=seed, band=band, method=method,
                                         frequency_index=fi,
                                         frequency=frequencies[fi] if fi < len(frequencies) else '',
                                         source_index=si, receiver_index=ri,
                                         candidate_cal_real=float(values[0].real),
                                         candidate_cal_imag=float(values[0].imag),
                                         candidate_task_real=float(values[1].real),
                                         candidate_task_imag=float(values[1].imag),
                                         reference_cal_real=float(values[2].real),
                                         reference_cal_imag=float(values[2].imag),
                                         reference_task_real=float(values[3].real),
                                         reference_task_imag=float(values[3].imag),
                                         candidate_cal_key=keys[0], candidate_task_key=keys[1],
                                         reference_cal_key=keys[2], reference_task_key=keys[3]))
                    count += 1
                del arrays
    for archive_path in (calibration_path, task_path, reference_path):
        ctx.sources.append(dict(path=ctx.relative(archive_path), role='existing complex endpoint data'))
    ctx.note('08 raw arrays', f'Exported {count} existing frequency/source/receiver rows. '
              'Plotted distances remain exactly the saved SCORED values, not recomputed values.')
    return [path.name]


def fixed_probe(ctx, stem):
    contexts = [c for c in ctx.contexts if any(c+(m,) in ctx.scored_by_key for m in METHODS)]
    if not contexts:
        ctx.note(stem, 'No scored fixed-probe endpoint rows; figure skipped.')
        return None
    fig, axes = plt.subplots(1, len(contexts), figsize=(4.7*len(contexts), 5.5),
                             sharey=True, squeeze=False)
    fig.subplots_adjust(left=.065, right=.98, bottom=.25, top=.78, wspace=.18)
    fig.suptitle('Unfitted-grid probe 8 diagnostics (1/24; not a qualified witness)',
                 x=.065, y=.96, ha='left', fontsize=16)
    records = []
    all_values = []
    qdiffs = []
    metrics = (('fixed_probe8_total_distance', 'Total task data', BLUE, -.18),
               ('fixed_probe8_induced_distance', 'Target-induced difference', ORANGE, .18))
    for column, context in enumerate(contexts):
        axis = axes[0, column]
        for index, method in enumerate(METHODS):
            score = ctx.scored_by_key.get(context+(method,))
            if score is None:
                ctx.note(stem, f'{context_name(context)} {method}: missing scored values; skipped.')
                continue
            endpoint = ctx.endpoint_by_key.get(context+(method,), {})
            qdiff = endpoint.get('joint_centroid_difference')
            if finite(qdiff):
                qdiffs.append(qdiff)
            if (finite(endpoint.get('task8_response_distance'))
                    and endpoint['task8_response_distance'] != score.get('fixed_probe8_total_distance')):
                ctx.note(stem, f'{context_name(context)} {method}: endpoint/scored total differ; '
                          'plotted value remains SCORED_RESULTS.json.')
            for key, label, color, offset in metrics:
                value = score.get(key)
                if not finite(value):
                    ctx.note(stem, f'{context_name(context)} {method}: missing {key}; no bar inserted.')
                    continue
                axis.bar(index+offset, value, width=.34, color=color,
                         label=label if index == 0 else None)
                axis.annotate(f'{value:.2f}', (index+offset, value), xytext=(0, 4),
                              textcoords='offset points', ha='center', fontsize=7.2, rotation=90)
                all_values.append(value)
            records.append(dict(score, recorded_same_target_centroid_difference=qdiff,
                                endpoint_status=endpoint.get('status'),
                                endpoint_source_file=endpoint.get('source_file'),
                                plotted_metrics_source='SCORED_RESULTS.json values without recomputation'))
        axis.set_title(context_name(context), loc='left', fontsize=10.5)
        axis.set_xticks(range(4), METHOD_LABELS, fontsize=8)
        axis.set_axisbelow(True)
        axis.grid(axis='y', alpha=.5)
    if not all_values:
        plt.close(fig)
        ctx.note(stem, 'No finite scored fixed-probe distances; figure skipped.')
        return None
    axes[0, 0].set_ylim(0., max(all_values)*1.32)
    axes[0, 0].set_ylabel('Saved full complex-data distance / 0.01')
    axes[0, 0].legend(loc='upper left', fontsize=8)
    qtext = ('All recorded same-target centroid differences are 0.'
             if qdiffs and all(value == 0 for value in qdiffs)
             else 'Recorded centroid differences are retained in CSV; missing values are not filled.')
    footer(fig, 'Total: candidate task data versus reference task data. Induced: difference '
           'between the two task-minus-calibration responses. Both are saved SCORED values.\n'
           f'Same fixed target in both worlds; {qtext} Joint compensation was NOT_RUN. '
           'No candidate met d_cal ≤ 1.\n'
           'Unfitted target geometry failed the 0.1 refinement check; three selected endpoints '
           'have a separate circle-fitted audit. No full supremum, posterior or speedup claim.')
    count = write_csv(ctx.output/f'{stem}.csv', records)
    extra = export_probe_arrays(ctx)
    return fig, count, 'Saved same-target probe8 distances only: 1/24 coverage, no joint compensation.', extra


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--outdir', type=Path, help='default: ROOT/figures')
    parser.add_argument('--dpi', type=int, default=190)
    args = parser.parse_args()
    if not 72 <= args.dpi <= 400:
        parser.error('--dpi must be between 72 and 400')
    root = args.root.resolve()
    output = (args.outdir or root/'figures').resolve()
    output.mkdir(parents=True, exist_ok=True)
    styling()
    ctx = Evidence(root, output)
    manifest = dict(generated_utc=datetime.now(timezone.utc).isoformat(),
                    output=str(output), figures=[], physical_lu=0, physical_rhs=0,
                    physical_solves=0, display_grid=[DISPLAY_NX, DISPLAY_NY],
                    selection_rule='corrected online d_cal, lexical method name for exact ties',
                    excluded_run='results/p0b/search; p0b-search-001',
                    standalone_rhs_role='offline attribution, not actual-ledger sum')
    builders = [('06_p0b_calibration_and_rhs', calibration),
                ('07_p0b_offline_material_comparison', offline_images),
                ('08_p0b_fixed_probe8', fixed_probe)]
    failures = 0
    for stem, build in builders:
        item = dict(stem=stem, status='SKIPPED', outputs=[])
        try:
            result = build(ctx, stem)
            if result is None:
                manifest['figures'].append(item)
                continue
            fig, count, caption, extra = result
            item.update(csv_rows=count, caption=caption, outputs=[f'{stem}.csv']+extra)
            for extension in ('png', 'svg'):
                path = output/f'{stem}.{extension}'
                fig.savefig(path, dpi=args.dpi, bbox_inches='tight')
                item['outputs'].append(path.name)
            plt.close(fig)
            item['status'] = 'GENERATED_NOT_YET_VISUALLY_REVIEWED'
            print(f'{stem}: generated PNG/SVG; {count} plot CSV rows', flush=True)
        except Exception as error:
            failures += 1
            item.update(status='FAILED', error=str(error))
            failure = dict(utc=datetime.now(timezone.utc).isoformat(), figure=stem,
                           error_type=type(error).__name__, error=str(error),
                           traceback=traceback.format_exc())
            with (output/'P0B_PLOT_FAILURES.jsonl').open('a') as stream:
                stream.write(json.dumps(failure, ensure_ascii=False)+'\n')
            ctx.note(stem, f'Plot failure preserved in P0B_PLOT_FAILURES.jsonl: {error}')
            plt.close('all')
            print(f'{stem}: FAILED: {type(error).__name__}: {error}', flush=True)
        manifest['figures'].append(item)
        gc.collect()
    raw_outputs = []
    for source, stem in (('results/p0b/SCORED_RESULTS.json', 'p0b_raw_scored_results'),
                         ('results/p0b/search-corrected/RESULTS.json', 'p0b_raw_corrected_search'),
                         ('results/p0b/task-endpoints/RESULTS.json', 'p0b_raw_task_endpoints')):
        name = f'{stem}.csv'
        write_csv(output/name, ctx.raw.get(source, []))
        raw_outputs.append(name)
    generated = sum(item['status'] == 'GENERATED_NOT_YET_VISUALLY_REVIEWED'
                    for item in manifest['figures'])
    manifest.update(sources=ctx.sources, notes=ctx.notes, raw_csv_outputs=raw_outputs,
                    recorded_contexts=len(ctx.contexts),
                    completed_contexts=sum(len(ctx.online_rows(c)) == 4 for c in ctx.contexts),
                    recorded_online_method_rows=len(ctx.candidates),
                    not_run_contexts=ctx.not_run, generated_figures=generated,
                    plot_failures=failures,
                    status='PARTIAL_PLOT_FAILURE' if failures else
                    'GENERATED_PENDING_VISUAL_REVIEW' if generated == 3 else
                    'PARTIAL_INPUTS_PENDING_VISUAL_REVIEW')
    (output/'P0B_PLOT_MANIFEST.json').write_text(json.dumps(manifest, indent=2,
                                                          ensure_ascii=False, allow_nan=False)+'\n')
    notes = ['# P0-B figure sources and limits', '',
             f"Generated: {manifest['generated_utc']}", '',
             'These are saved corrected results, not new searches or field solves. '
             'Raw result files and the original ledgers are read-only.', '',
             '## Current outputs', '']
    notes.extend(f"- {item['stem']}: {item['status']}" for item in manifest['figures'])
    notes.extend(['', 'Only outputs listed as generated in the current manifest belong to this '
                  'snapshot; earlier files may remain after a skipped/failed rebuild.', '',
                  '## Interpretation and completeness', '',
                  '- Images use a 256×192 display grid; original theta coordinates are exported '
                  'separately. Continuous known truth is image-only and never selects a candidate.',
                  '- Standalone RHS uses method_rhs+80 for baselines and +160 for TASK_AWARE. '
                  'These are offline attributions, not sums of actual shared ledger work.',
                  '- Fixed probe8 covers 1/24 tasks and uses the same target in both worlds; '
                  'the recorded centroid differences are not a joint-compensation experiment.', ''])
    notes.extend(f"- [{row['scope']}] {row['note']}" for row in ctx.notes)
    (output/'P0B_PLOT_NOTES.md').write_text('\n'.join(notes)+'\n')
    print(f'Manifest: {output/"P0B_PLOT_MANIFEST.json"}', flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
