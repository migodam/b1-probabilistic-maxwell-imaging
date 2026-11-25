"""Draw publication vectors and tabulate final spatial classifications."""
import csv
import json

from run import HERE, RESULTS, get_freeze, write

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
import numpy as np


COLORS = ('#ffffff', '#21757a', '#c35826', '#d8dfe5')
CODE = {'EXCLUDED': 1, 'WITNESS_FOUND': 2, 'UNRESOLVED': 3}
FIGURES = HERE/'figures'
STYLE = dict(font_size=8.8, dark='#182831')


def totals(counts):
    return f"{counts['EXCLUDED']} excluded  ·  {counts['WITNESS_FOUND']} witnessed  ·  {counts['UNRESOLVED']} unresolved"


def headers():
    return [Patch(facecolor=c, edgecolor='none', label=label) for c, label in zip(
        COLORS[1:], ('Uniformly excluded', 'Feasible witness', 'Unresolved'), strict=True)]


def grid(ax, s, rows, coordinate, win=None):
    g = np.zeros((s, s), dtype=int)
    for row in rows:
        x, y = row['center']
        g[y, x] = CODE[row['status']]
    if win:
        x0, x1, y0, y1 = win
    else:
        x0, x1, y0, y1 = 0, s-1, 0, s-1
    ax.pcolormesh(np.arange(x0, x1+2)-.5, np.arange(y0, y1+2)-.5,
                  g[y0:y1+1, x0:x1+1], cmap=ListedColormap(COLORS),
                  norm=BoundaryNorm(np.arange(-.5, 4.5), 4), edgecolors='white',
                  linewidth=.32, rasterized=False)
    ax.set_aspect('equal')
    ax.set(xlim=(x0-.5, x1+.5), ylim=(y0-.5, y1+.5),
           xticks=np.linspace(x0, x1, 3, dtype=int), yticks=np.linspace(y0, y1, 3, dtype=int))
    ax.set_xticks(np.arange(x0, x1+1)-.5, minor=True)
    ax.set_yticks(np.arange(y0, y1+1)-.5, minor=True)
    ax.grid(which='minor', color='white', linewidth=.32)
    ax.tick_params(which='minor', length=0)
    ax.tick_params(which='major', length=2, labelsize=7)
    ax.set_xlabel(coordinate, labelpad=2)
    for spine in ax.spines.values():
        spine.set_color('#9aaab3')
        spine.set_linewidth(.5)


def amplitude_bars(ax, trees):
    y = np.arange(len(trees))[::-1]
    left = np.zeros(len(trees))
    counts = [r['status_counts'] for r in trees]
    for status, color in zip(CODE, COLORS[1:], strict=True):
        vals = np.array([c[status] for c in counts])
        ax.barh(y, vals, left=left, color=color, height=.62)
        left += vals
    ax.set_yticks(y, [f"{r['low']:g}–{r['cap']:g}" for r in trees])
    ax.set(xlim=(0, left.max()), ylim=(-.55, len(trees)-.45))
    ax.set_xticks([0, int(left.max()/2), int(left.max())])
    ax.set_xlabel('Number of centers', labelpad=2)
    ax.tick_params(which='major', length=2, labelsize=7)
    for name, spine in ax.spines.items():
        spine.set_visible(name == 'bottom')
        spine.set_color('#9aaab3')
        spine.set_linewidth(.5)


def write_artifacts(fig, name):
    fig.savefig(FIGURES/f'{name}.pdf', bbox_inches='tight', pad_inches=.12)
    fig.savefig(FIGURES/f'{name}.svg', bbox_inches='tight', pad_inches=.12)
    fig.savefig(FIGURES/f'{name}.png', dpi=150, bbox_inches='tight', pad_inches=.12)
    plt.close(fig)


def map_figure(cases):
    fig, axes = plt.subplots(3, 4, figsize=(10.3, 8.5),
                             gridspec_kw={'width_ratios': (1, 1, 1, 1.15)})
    fig.subplots_adjust(left=.09, right=.98, bottom=.16, top=.86, hspace=.61, wspace=.32)
    titles = ('Unknown center, vacuum', 'Unknown center, weak background',
              'Any high-valued cell', 'Fixed shape, amplitude classes')
    for j, title in enumerate(titles):
        axes[0, j].set_title(title, loc='center', y=1.09, fontsize=8.7, fontweight='bold')
    for i, (name, r) in enumerate(cases):
        for j, m in enumerate(r['center_trees']):
            grid(axes[i, j], 12, m['final_centers'], 'Candidate center (cell index)')
        grid(axes[i, 2], 12, r['arbitrary_high']['rows'], 'Tested cell index')
        amplitude_bars(axes[i, 3], r['amplitude_trees'])
        axes[i, 0].set_ylabel(name, fontsize=9.1, fontweight='bold', labelpad=11)
    fig.text(.085, .955, 'What the same finite-aperture data can rule out', fontsize=17,
             color=STYLE['dark'], fontweight='bold')
    fig.text(.086, .920, '12 × 12 cells  /  7 incident fields × 9 receiver channels  /  3% stored-data noise',
             fontsize=9.5, color='#506070')
    fig.legend(handles=headers(), loc='lower left', ncol=3, bbox_to_anchor=(.088, .054), frameon=False)
    fig.text(.097, .024, 'White cells are outside the declared center window. Right-hand intervals scale the prescribed profile.',
             fontsize=8, color='#506070')
    write_artifacts(fig, 'spatial-map')


def native_figure(case, r):
    fig, ax = plt.subplots(figsize=(5.6, 4.25))
    fig.subplots_adjust(left=.14, right=.92, bottom=.159, top=.796)
    tree = r['center_trees'][0]
    grid(ax, 96, tree['final_centers'], 'Center x coordinate (cell index)', case['center_window'])
    ax.set_ylabel('Center y coordinate (cell index)')
    fig.text(.15, .946, 'Native 96 × 96 disk cover', fontsize=15, color=STYLE['dark'], fontweight='bold')
    fig.text(.15, .897, '7 incident fields × 9 receiver channels, vacuum exterior', fontsize=8.9, color='#506070')
    fig.text(.15, .855, totals(tree['status_counts']), fontsize=9, color='#506070')
    fig.legend(handles=headers(), loc='lower center', ncol=3, bbox_to_anchor=(.48, -.009), frameon=False)
    write_artifacts(fig, 'native-center-map')


def time_figure(runs):
    fig, (upper, lower) = plt.subplots(2, 1, figsize=(7.4, 5.0), sharex=True,
                                      gridspec_kw={'height_ratios': (1, 1)})
    fig.subplots_adjust(left=.14, right=.97, top=.81, bottom=.115, hspace=.43)
    x = [m['count'] for m in runs]
    shared = ('common_geometry_setup_with_persistence', 'common_reading_and_model_load',
              'candidate_generation_for_largest_batch')
    full = ('full_float_cold_setup', 'full_float_all_assembly_factorizations')
    fixed = ('fixed_float_cold_setup', 'fixed_float_all_reduced_examinations')
    strict = ('fixed_arb_cold_setup', 'fixed_arb_all_nodes')
    for fields, label, color, marker in ((full, 'Full state', '#4a5871', 's'),
                                          (fixed, 'Fixed traces', '#21757a', 'o')):
        costs = [sum(m[k]['wall_seconds'] for k in (*shared, *fields)) for m in runs]
        upper.plot(x, costs, marker=marker, color=color, lw=1.8, label=label)
    costs = [sum(m[k]['wall_seconds'] for k in (*shared, *fixed, *strict)) for m in runs]
    lower.plot(x, costs, marker='o', color='#21757a', lw=1.8, label='Fixed traces, all boxes')
    sampled = [m['independent_direct_arb_replay']['wall_seconds'] for m in runs]
    lower.plot(x, sampled, marker='s', color='#4a5871', lw=1.5, linestyle='--',
               label='Full state, independent sample')
    for i, m in enumerate(runs):
        lower.annotate(f"{m['independent_replayed_boxes']} boxes",
                       (x[i], sampled[i]), textcoords='offset points', xytext=(5, -9), fontsize=7)
    for ax, caption in ((upper, 'Binary64 comparison, all boxes'),
                        (lower, '128-bit independent verification costs')):
        ax.grid(axis='y', color='#e8ecef')
        ax.spines[['top', 'right']].set_visible(False)
        ax.set_ylabel('Wall time (s)', labelpad=4)
        ax.set_xscale('log')
        ax.set_xticks(x, [str(v) for v in x])
        ax.minorticks_off()
        ax.set_title(caption, fontsize=9.5, loc='left', fontweight='bold')
        ax.legend(loc='upper left', frameon=False, fontsize=8.1)
    lower.set_xlabel('Distinct candidate boxes')
    fig.text(.14, .934, 'Fixed traces amortize identical border tests', fontsize=15,
             color=STYLE['dark'], fontweight='bold')
    fig.text(.14, .878, '151 full and 25 fixed-trace bordered coordinates per candidate',
             fontsize=9, color='#506070')
    write_artifacts(fig, 'map-amortization')


def export_rows(cases):
    with (RESULTS/'map_centers.csv').open('w', newline='') as file:
        out = csv.writer(file)
        out.writerow(['observation_id', 'object_id', 'task', 'x_index', 'y_index', 'status', 'proof_node'])
        for r in cases:
            for tree in r['center_trees']+r['amplitude_trees']:
                for center in tree['final_centers']:
                    out.writerow([r['id'], r['object_id'], tree['id'], *center['center'],
                                  center['status'], center['backing_node']])
            for cell in r['arbitrary_high'].get('rows', []):
                out.writerow([r['id'], r['object_id'], 'arbitrary_high', *cell['center'],
                              cell['status'], cell['id']])


def summary(cases, benchmark):
    offline = json.loads((RESULTS/'offline_score.json').read_text())
    replay = json.loads((RESULTS/'independent_replay.json').read_text())
    maps = []
    for r in cases:
        tree = []
        for t in r['center_trees']+r['amplitude_trees']:
            evidence = [node['certificate'] for node in t['nodes'] if node['certificate']]
            witness_ratios = [node['witness']['residual_over_radius_upper'] for node in t['nodes']
                              if node.get('witness')]
            tree.append(dict(contract=t['id'], background=t['b'], low=t['low'], cap=t['cap'],
                             candidate_centers=t['total_centers'], counts=t['status_counts'],
                             distinct_boxes=t['distinct_boxes'], nodes=len(t['nodes']),
                             surviving_centers_exclusive_bounds=t['total_centers']-t['status_counts']['EXCLUDED'],
                             min_certified_residual_lower_over_radius=min((q['residual_over_radius_lower']
                                 for q in evidence), default=None), witness_residual_upper_over_radius=witness_ratios,
                             node_budget_deferrals=sum(node.get('reason', '').startswith('retained vacuum') for node in t['nodes']),
                             seconds=t['elapsed']))
        surviving_bins = [[t['low'], t['cap']] for t in r['amplitude_trees']
                          if t['status_counts']['EXCLUDED'] < t['total_centers']]
        peak = r['arbitrary_high']
        peak_low = peak.get('global_low_peak_box', {}).get('certificate')
        maps.append(dict(id=r['id'], object_id=r['object_id'], grid=r['grid'],
                          center_and_amplitude_maps=tree, surviving_amplitude_bins=surviving_bins,
                          amplitude_surviving_hull=[surviving_bins[0][0], surviving_bins[-1][1]] if surviving_bins else None,
                          amplitude_hull_width=surviving_bins[-1][1]-surviving_bins[0][0] if surviving_bins else None,
                          full_arbitrary_high_cover=peak['union_complete'],
                          arbitrary_high=peak.get('status_counts'),
                          global_low_peak_cap=peak.get('global_low_peak_box', {}).get('cap'),
                          global_low_peak_excluded=peak_low is not None,
                          global_low_peak_residual_lower_over_radius=peak_low['residual_over_radius_lower'] if peak_low else None,
                          generator_status=offline[r['id']], seconds_online=r['total_online'],
                          native_model_replay_cost=r['native_arb_idft_setup'],
                          independent_replay={k:v for k,v in replay[r['id']].items()
                                              if k not in ('fixed_rows', 'full_rows', 'witness_rows')}))
    out = dict(scope='four distinct synthetic objects, all fixed-aperture conditions paired within object',
               all_counts_are_task_labels_or_boxes=True,
               decisions=['SPATIAL_SELECTIVITY_FOR_DECLARED_DISCRETE_CONTRACT',
                          'NO_UNIFORM_SINGLE_CELL_HIGH_EXCLUSION_IN_COMPLETE_SMALL_GRID_COVER',
                          'FIXED_TRACE_INERTIA_EQUIVALENCE_WITH_AMORTIZED_RUNTIME'],
               maps=maps, amortization=[{k:v for k,v in r.items() if k != 'boxes'} for r in benchmark])
    write(RESULTS/'summary.json', out)


if __name__ == '__main__':
    FIGURES.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': STYLE['font_size'],
                         'text.color': STYLE['dark'], 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    freeze = get_freeze()
    all_cases = [json.loads((RESULTS/f"{cfg['id']}.json").read_text()) for cfg in freeze['cases']]
    cases = list(zip(('Single disk', 'Two disks', 'Annulus'), all_cases[:3], strict=True))
    map_figure(cases)
    native_figure(freeze['cases'][-1], all_cases[-1])
    bench = [json.loads((RESULTS/f'benchmark_{count}.json').read_text())
             for count in freeze['benchmark_box_counts']]
    time_figure(bench)
    export_rows(all_cases)
    summary(all_cases, bench)
