#!/usr/bin/env python3
"""Cached-only tail plots. Reads saved JSON/NPZ; never runs physics or rescoring."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gc
import json
from pathlib import Path
import traceback
import zipfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np

from plot_q import (BLUE, GREEN, PURPLE, ORANGE, GRAY, finite, styling,
                    write_csv, footer, sample_world)


METHODS = ('FULL_J', 'RANDOM_MATCHED', 'MULTISTART_FWI', 'TASK_AWARE')
LABELS = ('Full J', 'Random matched', 'Multistart FWI', 'Task aware')
COLORS = (BLUE, GREEN, PURPLE, ORANGE)
MARKERS = ('o', 's', 'D', '^')
SIZES = (55, 90, 150, 70)


def key(row):
    return row.get('seed'), row.get('band'), row.get('method')


def context_order(context):
    return context[0], 0 if context[1] == 'wide' else 1


class Cache:
    def __init__(self, root, output, max_field_points):
        self.root, self.output = root, output
        self.max_field_points = max_field_points
        self.notes = [dict(scope='scope', note='Only corrected scored search results are used. '
                           'Earlier results/p0b/search (p0b-search-001) is invalid for method '
                           'comparison and excluded; its costs remain in the original ledger.')]
        self.sources = []
        self.raw = {}
        self.scores = self.load('results/p0b/SCORED_RESULTS.json')
        self.unfitted = self.load('results/p0b/finer-endpoints/RESULTS.json')
        self.fitted = self.load('results/p0b/fitted-task-endpoints/RESULTS.json')
        self.mesh = self.load('results/p0b/fitted-task-endpoints/MESH_DIAGNOSTICS.json')
        for row in self.mesh:
            self.note('skinny fitted triangles',
                      f"circle n={row.get('circle_polygon_n')}: "
                      f"min angle={row.get('min_triangle_angle_degrees')} degrees; "
                      f"max edge/min altitude={row.get('max_edge_to_min_altitude_ratio')}; "
                      f"min triangle area={row.get('min_triangle_area')}; "
                      f"P3 DOFs={row.get('p3_dofs')}; "
                      f"requested_min_angle={row.get('requested_min_angle')}; "
                      f"requested_max_area={row.get('requested_max_area')}; "
                      f"Triangle options={row.get('triangle_options')}. "
                      'These mesh diagnostics are retained explicitly; circle-refinement '
                      'consistency is not a continuum error bound.')
        if not self.mesh:
            self.note('skinny fitted triangles', 'Fitted rows flag skinny_mesh=True, but mesh '
                      'diagnostics are unavailable in this snapshot; no values are invented.')

    def note(self, scope, message):
        self.notes.append(dict(scope=scope, note=message))

    def load(self, relative):
        path = self.root/relative
        try:
            data = json.loads(path.read_text())
            if not isinstance(data, list):
                raise ValueError('expected a JSON list')
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


def cost_scatter(ctx, stem):
    rows = [r for r in ctx.scores if r.get('method') in METHODS
            and finite(r.get('seed')) and r.get('band') in ('wide', 'narrow')
            and finite(r.get('standalone_proposal_rhs')) and r['standalone_proposal_rhs'] > 0
            and finite(r.get('d_cal')) and r['d_cal'] > 0]
    if not rows:
        ctx.note(stem, 'No finite scored calibration/cost pairs; figure skipped.')
        return None
    contexts = sorted({(r['seed'], r['band']) for r in rows}, key=context_order)
    seeds = sorted({r['seed'] for r in rows})
    fig, axes = plt.subplots(1, len(contexts), figsize=(4.25*len(contexts), 4.9),
                             sharex=True, sharey=True, squeeze=False)
    fig.subplots_adjust(left=.07, right=.98, bottom=.24, top=.71, wspace=.16)
    fig.suptitle(f'Calibration / cost scatter: {len(contexts)} contexts, {len(seeds)} independent seeds',
                 x=.07, y=.96, ha='left', fontsize=15)
    records = []
    all_x = [r['standalone_proposal_rhs'] for r in rows]
    all_y = [r['d_cal'] for r in rows]
    for column, context in enumerate(contexts):
        axis = axes[0, column]
        for index, method in enumerate(METHODS):
            selected = [r for r in rows if (r['seed'], r['band']) == context and r['method'] == method]
            if not selected:
                ctx.note(stem, f'{context} {method}: missing pair; no point inserted.')
                continue
            row = selected[-1]
            # Open, nested markers preserve exact coincident coordinates.
            axis.scatter(row['standalone_proposal_rhs'], row['d_cal'],
                         marker=MARKERS[index], s=SIZES[index],
                         facecolors=COLORS[index] if index in (0, 3) else 'none',
                         edgecolors=COLORS[index], linewidths=1.5,
                         zorder=5-index*.3)
            records.append(dict(row, x_metric='standalone_proposal_rhs', y_metric='d_cal',
                                x_interpretation='offline standalone attribution, not actual-ledger total',
                                display_jitter=0, frontier_fit=False))
        axis.axhline(1., color='#303A45', ls=':', lw=1.4)
        axis.text(.97, 1., 'd_cal = 1', transform=axis.get_yaxis_transform(),
                  ha='right', va='bottom', fontsize=8)
        axis.set_yscale('log')
        axis.set_ylim(.7, max(all_y)*1.4)
        axis.set_xlim(min(all_x)-20, max(all_x)+20)
        axis.set_title(f'Seed {context[0]} / {context[1]}', loc='left', fontsize=10.5)
        axis.set_xlabel('Standalone proposal RHS (attributed)')
        axis.set_axisbelow(True)
        axis.grid(which='major', alpha=.5)
    axes[0, 0].set_ylabel('Saved calibration distance d_cal')
    handles = [Line2D([], [], ls='', marker=MARKERS[i], color=COLORS[i],
                      markerfacecolor=COLORS[i] if i in (0, 3) else 'none',
                      markersize=7, label=LABELS[i]) for i in range(4)]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.54, .875), ncol=4,
               fontsize=9)
    footer(fig, 'Standalone x values are saved offline cost attributions; actual shared setup '
           'was paid once in the ledger. Coincident points use nested markers, without jitter.\n'
           f'Minimum shown d_cal={min(all_y):.6f} > 1. Scatter only: no fitted frontier, '
           'qualified-witness or end-to-end speedup claim; two seeds do not establish algorithmic advantage.')
    count = write_csv(ctx.output/f'{stem}.csv', records)
    return fig, count, 'Saved per-context calibration distance versus offline standalone RHS attribution.'


def refinement(ctx, stem):
    unfitted = {key(r): r for r in ctx.unfitted if r.get('state') == 'task'
                and finite(r.get('mesh_white')) and r['mesh_white'] > 0}
    fitted = {key(r): r for r in ctx.fitted if r.get('polygon_n') == 512
              and finite(r.get('fitted_refinement_white')) and r['fitted_refinement_white'] > 0}
    matched = sorted(set(unfitted) & set(fitted), key=lambda k: context_order(k[:2]))
    for missing in sorted(set(unfitted) ^ set(fitted)):
        ctx.note(stem, f'{missing}: one comparison is missing; unmatched candidate skipped.')
    if not matched:
        ctx.note(stem, 'No matched selected-endpoint refinement pairs; figure skipped.')
        return None
    fig, axis = plt.subplots(figsize=(10.6, 5.4))
    fig.subplots_adjust(left=.1, right=.98, bottom=.26, top=.72)
    fig.suptitle('Selected task endpoints: two empirical consistency checks',
                 x=.1, y=.96, ha='left', fontsize=15)
    records = []
    all_values = []
    for index, selected in enumerate(matched):
        old, fit = unfitted[selected], fitted[selected]
        for offset, row, metric, color, label in (
                (-.18, old, 'mesh_white', ORANGE, 'Original unfitted P3 → subdivided P2'),
                (.18, fit, 'fitted_refinement_white', BLUE, 'Fitted P3 circle: 256 → 512')):
            value = row[metric]
            axis.bar(index+offset, value, width=.34, color=color,
                     label=label if index == 0 else None)
            axis.annotate(f'{value:.6g}', (index+offset, value), xytext=(0, 5),
                          textcoords='offset points', ha='center', fontsize=8)
            all_values.append(value)
            records.append(dict(row, comparison=label, plotted_metric=metric,
                                plotted_white_difference=value,
                                interpretation='empirical discrete consistency, not continuum error bound'))
    axis.axhline(.1, color='#303A45', ls=':', lw=1.5, label='Original consistency line: 0.1')
    axis.set_yscale('log')
    axis.set_ylim(min(all_values)*.55, max(all_values)*2)
    axis.set_xticks(range(len(matched)), [f'Seed {k[0]} / {k[1]}\n{k[2]}' for k in matched], fontsize=9)
    axis.set_ylabel('Full task-data refinement difference / 0.01')
    axis.set_axisbelow(True)
    axis.grid(axis='y', which='major', alpha=.5)
    fig.legend(*axis.get_legend_handles_labels(), loc='upper center',
               bbox_to_anchor=(.54, .88), ncol=1, fontsize=8.5)
    angles = [r['min_triangle_angle_degrees'] for r in ctx.mesh
              if finite(r.get('min_triangle_angle_degrees'))]
    ratios = [r['max_edge_to_min_altitude_ratio'] for r in ctx.mesh
              if finite(r.get('max_edge_to_min_altitude_ratio'))]
    diagnostic = (f'Fitted meshes contain skinny triangles: min angle {min(angles):.6g}°; '
                  f'max edge/min altitude {max(ratios):.6g}.' if angles and ratios else
                  'Fitted rows flag skinny triangles; available geometry diagnostics are retained in notes/CSV.')
    footer(fig, 'Same three online-selected candidates and fixed probe 8 (1/24). '
           'The two series apply different discrete refinements; neither is a continuum bound.\n'
           +diagnostic+'\nThese checks do not qualify all twelve candidates or establish a '
           'calibration witness; the original unfitted task inconsistency remains recorded.')
    count = write_csv(ctx.output/f'{stem}.csv', records)
    return fig, count, 'Original unfitted task refinement and fitted-circle 256-to-512 consistency for matched selected candidates.'


def saved_fields(ctx, stem):
    paths = []
    for seed in (2026102201, 2026102202):
        path = ctx.root/f'results/wall-qualification-resume/WORLD_{seed}_h0.145.npz'
        if path.exists():
            paths.append((seed, path))
        else:
            ctx.note(stem, f'Missing saved WORLD for original seed {seed}; no field is generated.')
    if not paths:
        ctx.note(stem, 'No cached original-seed task fields; figure skipped.')
        return None
    fig, axes = plt.subplots(len(paths), 2, figsize=(10.4, 3.8*len(paths)), squeeze=False)
    fig.subplots_adjust(left=.07, right=.86, bottom=.14, top=.88, wspace=.18, hspace=.32)
    fig.suptitle('Saved task total field: amplitude and phase, f/f₀=1, source 1',
                 x=.07, y=.97, ha='left', fontsize=14.5)
    amplitude_artists, phase_artists, records = [], [], []
    maximum_amplitude = 0.
    for row_index, (seed, path) in enumerate(paths):
        try:
            # The helper reads one archive member at a time and selects exact
            # stored quadrature values. It does not interpolate or solve.
            with np.load(path, allow_pickle=False) as archive:
                points, epsilon, values, indices, original_count = sample_world(
                    archive, 'task', ctx.max_field_points, phase=True)
            del epsilon
            valid = np.isfinite(points).all(axis=1) & np.isfinite(values)
            if not valid.all():
                ctx.note(stem, f'{path.name}: {int((~valid).sum())} nonfinite display samples skipped.')
            points, values, indices = points[valid], values[valid], indices[valid]
            if not len(values):
                raise ValueError('no finite saved task-field samples')
            amplitude, phase = np.abs(values), np.angle(values)
            maximum_amplitude = max(maximum_amplitude, float(amplitude.max()))
            for column, (data, cmap) in enumerate(((amplitude, 'viridis'),
                                                  (phase, 'twilight_shifted'))):
                axis = axes[row_index, column]
                artist = axis.scatter(points[:, 0], points[:, 1], c=data, s=3.2,
                                      marker='s', linewidths=0, rasterized=True,
                                      cmap=cmap, vmin=0. if column == 0 else -np.pi,
                                      vmax=max(maximum_amplitude, 1e-15) if column == 0 else np.pi)
                (amplitude_artists if column == 0 else phase_artists).append(artist)
                axis.set(xlim=(-2, 2), ylim=(0, 3), aspect='equal',
                         xlabel='x / λ₀', ylabel='y / λ₀')
                axis.set_title(f'Seed {seed}: '+('amplitude |E_total|' if column == 0
                                                else 'phase arg E_total'), fontsize=10)
                axis.text(.02, .025, f'{len(values):,}/{original_count:,} saved q-points',
                          transform=axis.transAxes, fontsize=7.2,
                          bbox=dict(fc='white', ec='none', alpha=.85, pad=2))
            relative = str(path.relative_to(ctx.root))
            for point, value, index in zip(points, values, indices):
                records.append(dict(source_file=relative, seed=seed, h_label=.145,
                                    material_state='known task world', frequency=1.,
                                    source_index=0, source_display_label=1,
                                    saved_qpoint_index=int(index), original_qpoint_count=original_count,
                                    displayed_qpoint_count=len(values),
                                    x=float(point[0]), y=float(point[1]),
                                    total_field_real=float(value.real),
                                    total_field_imag=float(value.imag),
                                    total_field_amplitude=float(abs(value)),
                                    total_field_phase_rad=float(np.angle(value)),
                                    sampling='deterministic equally-spaced saved qpoint indices; no interpolation'))
            ctx.sources.append(dict(path=relative, role='existing f1/source1 task total field',
                                     display_samples=len(values), saved_qpoints=original_count))
            del points, values, indices, amplitude, phase
            gc.collect()
        except (OSError, ValueError, KeyError, EOFError, zipfile.BadZipFile) as error:
            ctx.note(stem, f'Partial/unreadable {path.name} skipped: {error}')
            for axis in axes[row_index]:
                axis.text(.5, .5, f'Seed {seed}: saved field unavailable',
                          transform=axis.transAxes, ha='center', fontsize=9)
                axis.set_axis_off()
    if not amplitude_artists:
        plt.close(fig)
        ctx.note(stem, 'No readable saved task fields; figure skipped.')
        return None
    norm = Normalize(0., maximum_amplitude if maximum_amplitude > 0 else 1.)
    for artist in amplitude_artists:
        artist.set_norm(norm)
    amplitude_bar = fig.add_axes([.9, .55, .014, .25])
    fig.colorbar(amplitude_artists[0], cax=amplitude_bar, label='Total-field amplitude')
    phase_bar = fig.add_axes([.9, .2, .014, .25])
    bar = fig.colorbar(phase_artists[0], cax=phase_bar, label='Total-field phase (rad)')
    bar.set_ticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
                  labels=['−π', '−π/2', '0', 'π/2', 'π'])
    footer(fig, f'Actual saved FEM task-world total fields from {len(amplitude_artists)}/2 original seeds; '
           'these are known-world fields, not inverse reconstructions or new solves.\n'
           f'At most {ctx.max_field_points:,} stored quadrature values per seed are shown '
           'with deterministic subsampling and no interpolation. Amplitude uses one common sampled-value scale.')
    count = write_csv(ctx.output/f'{stem}.csv', records)
    if len(amplitude_artists) != 2:
        ctx.note(stem, f'Partial field snapshot: {len(amplitude_artists)}/2 original seeds displayed.')
    return fig, count, 'Actual saved source1/f1 task total-field amplitude and phase from original seed WORLD files.'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--outdir', type=Path, help='default: ROOT/figures')
    parser.add_argument('--dpi', type=int, default=190)
    parser.add_argument('--max-field-points', type=int, default=12000)
    args = parser.parse_args()
    if not 72 <= args.dpi <= 400:
        parser.error('--dpi must be between 72 and 400')
    if not 1 <= args.max_field_points <= 50000:
        parser.error('--max-field-points must be between 1 and 50000')
    root = args.root.resolve()
    output = (args.outdir or root/'figures').resolve()
    output.mkdir(parents=True, exist_ok=True)
    styling()
    ctx = Cache(root, output, args.max_field_points)
    manifest = dict(generated_utc=datetime.now(timezone.utc).isoformat(),
                    figures=[], physical_rhs=0, physical_lu=0, physical_solves=0,
                    rescoring_run=False, source_mode='cached JSON and saved WORLD NPZ only',
                    max_field_display_points_per_seed=args.max_field_points)
    failures = 0
    for stem, build in (('09_p0b_calibration_cost_scatter', cost_scatter),
                        ('10_p0b_task_refinement_comparison', refinement),
                        ('11_saved_task_total_field', saved_fields)):
        item = dict(stem=stem, status='SKIPPED', outputs=[])
        try:
            result = build(ctx, stem)
            if result is not None:
                fig, count, caption = result
                item.update(csv_rows=count, caption=caption, outputs=[f'{stem}.csv'])
                for extension in ('png', 'svg'):
                    path = output/f'{stem}.{extension}'
                    fig.savefig(path, dpi=args.dpi, bbox_inches='tight')
                    item['outputs'].append(path.name)
                plt.close(fig)
                item['status'] = 'GENERATED_NOT_YET_VISUALLY_REVIEWED'
                print(f'{stem}: generated PNG/SVG and {count} CSV rows', flush=True)
        except Exception as error:
            failures += 1
            item.update(status='FAILED', error=str(error))
            row = dict(utc=datetime.now(timezone.utc).isoformat(), figure=stem,
                       error_type=type(error).__name__, error=str(error),
                       traceback=traceback.format_exc())
            with (output/'TAIL_PLOT_FAILURES.jsonl').open('a') as stream:
                stream.write(json.dumps(row, ensure_ascii=False)+'\n')
            ctx.note(stem, f'Plot failure retained in TAIL_PLOT_FAILURES.jsonl: {error}')
            plt.close('all')
            print(f'{stem}: FAILED: {error}', flush=True)
        manifest['figures'].append(item)
        gc.collect()
    raw_names = []
    for source, name in (
            ('results/p0b/SCORED_RESULTS.json', 'tail_raw_scores.csv'),
            ('results/p0b/finer-endpoints/RESULTS.json', 'tail_raw_finer_endpoints.csv'),
            ('results/p0b/fitted-task-endpoints/RESULTS.json', 'tail_raw_fitted_endpoints.csv'),
            ('results/p0b/fitted-task-endpoints/MESH_DIAGNOSTICS.json', 'tail_raw_fitted_mesh_diagnostics.csv')):
        write_csv(output/name, ctx.raw.get(source, []))
        raw_names.append(name)
    generated = sum(item['status'] == 'GENERATED_NOT_YET_VISUALLY_REVIEWED'
                    for item in manifest['figures'])
    manifest.update(sources=ctx.sources, notes=ctx.notes, raw_csv_outputs=raw_names,
                    plot_failures=failures, generated_figures=generated,
                    status='PARTIAL_PLOT_FAILURE' if failures else
                    'GENERATED_PENDING_VISUAL_REVIEW' if generated == 3 else
                    'PARTIAL_INPUTS_PENDING_VISUAL_REVIEW')
    (output/'TAIL_PLOT_MANIFEST.json').write_text(json.dumps(manifest, indent=2,
                                                           ensure_ascii=False, allow_nan=False)+'\n')
    notes = ['# Cached P0-B tail figures', '',
             f"Generated: {manifest['generated_utc']}", '',
             'No physics, optimization, rescoring, or raw-result mutation is performed. '
             'Only outputs marked generated in the current manifest belong to this snapshot; '
             'older files may remain after skipped/failed rebuilds.', '',
             '## Figure scope', '',
             '- 09: saved per-context calibration/cost scatter; standalone RHS is offline '
             'attribution. Actual shared work remains counted once in the ledger. '
             'Only two independent seeds; no fitted frontier, witness or speedup claim.',
             '- 10: saved unfitted task refinement versus fitted P3 circle 256→512, '
             'for the same three online-selected endpoints. They are different empirical '
             'consistency tests, not continuum error bounds or twelve-candidate qualification.',
             '- 11: existing original-seed task total-field amplitude and phase at f1/source1. '
             'Deterministically sampled saved quadrature values, with no interpolation, '
             'new solve or reconstruction. Amplitude scale is common over displayed samples.', '',
             '## Geometry and input notes', '']
    notes.extend(f"- [{row['scope']}] {row['note']}" for row in ctx.notes)
    (output/'TAIL_PLOT_NOTES.md').write_text('\n'.join(notes)+'\n')
    print(f'Manifest: {output/"TAIL_PLOT_MANIFEST.json"}', flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
