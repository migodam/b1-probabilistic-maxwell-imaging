"""Redraw publication figures from the saved scientific results."""

import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BLUE, TEAL, ORANGE = '#22629B', '#00866A', '#BA5E19'
INK, PALE, GRID = '#162A3A', '#5D6F7A', '#DCE5E8'
plt.rcParams.update({'font.size': 10, 'font.family': 'DejaVu Sans', 'axes.labelcolor': INK,
                     'text.color': INK, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.edgecolor': '#AAB8C0', 'xtick.color': INK, 'ytick.color': INK,
                     'pdf.fonttype': 42, 'ps.fonttype': 42})


def read(path):
    return json.loads((ROOT / path).read_text(encoding='utf-8'))


def output(fig, name):
    fig.savefig(HERE / f'{name}.png', dpi=310, facecolor='white')
    fig.savefig(HERE / f'{name}.pdf', facecolor='white')
    plt.close(fig)


def heading(fig, title, detail):
    fig.text(.055, .952, title, fontsize=20, weight='bold', va='top')
    fig.text(.056, .892, detail, fontsize=10, color=PALE, va='top')


def log_axis(ax, ylabel, title):
    ax.set_title(title, loc='left', fontweight='bold', fontsize=11, pad=12)
    ax.set(xscale='log', yscale='log', xlabel=r'Adjacent component spacing $R/\lambda$',
           ylabel=ylabel, xlim=(15, 300000))
    ax.set_xticks([20, 200, 2000, 20000, 200000])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x/1000:g}k' if x >= 1000 else f'{x:g}'))
    ax.grid(which='major', color=GRID, linewidth=.7)
    ax.set_axisbelow(True)


def sci_lower_bound(value):
    exponent = int(np.floor(np.log10(value)))
    mantissa = np.floor(value / 10.**exponent * 1000) / 1000
    return '$' + rf'{mantissa:.3f}\times10^{{{exponent}}}$'


def root_table(fig, table, runs, verified):
    table.set_axis_off()
    table.set(xlim=(0, 1), ylim=(0, 1))
    table.text(0, 1., 'c   Interval-certified roots at one common local readout',
               weight='bold', fontsize=11, va='top')
    xs = (0., .175, .30, .43, .595, .76, .88)
    labels = ('Components', r'$R/\lambda$', 'Roots', 'Min. material gap', 'Min. cross-port gap', 'Search LU', 'Ball solves')
    for x, label in zip(xs, labels, strict=True):
        table.text(x, .71, label, fontsize=9, color=PALE, va='center')
    for i, cert in enumerate(verified['results']):
        pairwise = cert['pairwise_separation_certificates']
        low_m = min(pair['material_l2'][0] for pair in pairwise)
        low_data = min(pair['cross_ports_l2'][0] for pair in pairwise)
        shown_m = np.floor(low_m * 1e5) / 1e5
        entries = (str(cert['components']), f'{cert["spacing"]:,}', f'≥ {cert["certified_root_count"]}',
                   f'≥ {shown_m:.5f}', '≥ ' + sci_lower_bound(low_data),
                   str(runs[cert['id']]['cost']['full_state_lu']),
                   str(cert['state_enclosure_count'] + cert['preparation']['known_host_full_ball_solve']))
        y = .47 - i*.235
        table.axhline(y+.114, color=GRID, linewidth=.75)
        for j, (x, entry) in enumerate(zip(xs, entries, strict=True)):
            table.text(x, y, entry, fontsize=11.3, color=ORANGE if j == 2 else INK,
                       weight='bold' if j == 2 else 'normal', va='center')
    full_lu = sum(row['cost']['full_state_lu'] for row in runs.values())
    table.text(0, -.045, f'All tested spacings used {full_lu:,} search LU factorizations. '
               f'Interval certification, including model preparation, took {verified["cpu_seconds"]:.2f} CPU s.',
               fontsize=9, color=PALE, va='top')
    fig.text(.095, .053, 'Orange outlines mark certified settings. Unsuccessful six-slot searches at 20–20,000 λ give no bound on root count.',
             fontsize=8.6, color=PALE)
    reduced = read('multiplicity/conditioning/results.json')['coupled_certificates_at_first_complete_spacing']
    if [row['id'] for row in reduced] == ['L2_R30', 'L3_R30']:
        counts = [row['number_certified'] for row in reduced]
        fig.text(.095, .025, f'Separate two-slot control at 30 λ: ≥ {counts[0]} and ≥ {counts[1]} roots. '
                 'All distances describe mathematical test geometries.', fontsize=8.6, color=PALE)


def multiplicity():
    nominal = read('multiplicity/results/numeric_ladder.json')
    extra = read('multiplicity/results/intermediate/numeric_ladder.json')
    verified = read('multiplicity/results/krawczyk.json')
    runs = {row['id']: row for row in nominal['rows'] + extra['rows']}
    certified = {row['id'] for row in verified['results']}
    fig = plt.figure(figsize=(11.3, 7.05))
    heading(fig, 'Finite spacing preserves multiple material roots',
            'Non-symmetric six-slot components in 2D scalar TM. Each excitation records its local complex echo.')
    grid = fig.add_gridspec(2, 2, height_ratios=[2.9, 1.35], left=.095, right=.963,
                           bottom=.108, top=.815, hspace=.52, wspace=.38)
    gax, eax = [fig.add_subplot(grid[0, i]) for i in (0, 1)]
    log_axis(gax, r'Off-diagonal block norm $\|G_{\rm off}\|_F$', 'a   Cross-component propagation')
    log_axis(eax, r'Local data mismatch $\|F_R(m_{\rm seed})-y_*\|_2$', 'b   Correction to the local echo')
    for count, color, marker in ((2, BLUE, 'o'), (3, TEAL, 's')):
        rows = sorted((row for row in runs.values() if row['components'] == count),
                      key=lambda row: row['spacing_reference_wavelengths'])
        distance = [row['spacing_reference_wavelengths'] for row in rows]
        propagation = [row['off_block_green_frobenius'] for row in rows]
        branches = [[v['initial_local_residual'] for v in row['root_attempts']] for row in rows]
        median = [np.median(v) for v in branches]
        style = dict(color=color, marker=marker, markersize=4.7, linewidth=1.65, label=f'{count} components')
        gax.plot(distance, propagation, **style)
        eax.plot(distance, median, **style)
        eax.fill_between(distance, [min(v) for v in branches], [max(v) for v in branches], color=color, alpha=.22)
        for i, row in enumerate(rows):
            if row['id'] in certified:
                for ax, values in ((gax, propagation), (eax, median)):
                    ax.scatter(distance[i], values[i], facecolors='none', edgecolors=ORANGE,
                               s=100, linewidth=1.4, zorder=5)

    rr = np.geomspace(20, 200000, 50)
    gax.plot(rr, .49*(rr/20)**(-.5), '--', color=PALE, linewidth=1., label=r'Guide $R^{-1/2}$')
    eax.plot(rr, .00095*(rr/20)**(-1), '--', color=PALE, linewidth=1., label=r'Guide $R^{-1}$')
    for ax in (gax, eax):
        ax.legend(loc='upper right', frameon=False, fontsize=8.5)
    root_table(fig, fig.add_subplot(grid[1, :]), runs, verified)
    output(fig, 'multiplicity_scaling')


def contents_grid(ax, faces, caps):
    arms = ('base', 'angle', 'frequency', 'rich')
    names = [(f'disks{n}_{s}', f'{int(n)} slots, {s}') for n in ('08', '12', '16') for s in ('low', 'high')]
    subset = {(row['case'], row['arm']): row for row in faces['rows'] if row['width'] == .75}
    exact = caps['certified_fixed_reference_task_diameter_intervals']
    for r, (case, label) in enumerate(names):
        for c, arm in enumerate(arms):
            if subset[case, arm]['center_ball_both']:
                value, shade = '100%', '#D7EEE7'
            else:
                lower, upper = np.asarray(exact[arm]) / caps['prior_task_diameter'] * 100
                value = f'{np.floor(lower*10)/10:.1f}–{np.ceil(upper*10-1e-10)/10:.1f}%'
                shade = '#F7E5D4'
            ax.add_patch(Rectangle((c+.02, r+.04), .96, .92, color=shade, linewidth=0))
            ax.text(c+.5, r+.5, value, fontsize=8.9, ha='center', va='center', weight='bold')
    ax.set_xticks(np.arange(4)+.5, ['Base\n3 complex', 'Extra angles\n6 complex', 'Changed freq.\n6 complex', 'Rich\n24 complex'])
    ax.set_yticks(np.arange(6)+.5, [label for _, label in names])
    ax.set(xlim=(0, 4), ylim=(6, 0))
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0, pad=8)
    ax.xaxis.tick_top()
    return subset


def margin_points(ax, subset):
    rows = [('angle', 'roi_only', 'Angles, ROI only', BLUE), ('angle', 'linear', 'Angles, outside fit', BLUE),
            ('frequency', 'roi_only', 'Frequency, ROI only', ORANGE), ('frequency', 'linear', 'Frequency, outside fit', ORANGE)]
    ax.axvspan(.72, 1, color='#EAF4EE', zorder=0)
    ax.axvline(1, color=PALE, linestyle='--', linewidth=1.15, zorder=1)
    for i, (arm, proposal, label, color) in enumerate(rows):
        row = subset['disks08_low', arm]
        delta, record = row['delta'], row[proposal]
        distances = sorted(np.mean(bounds) / delta for bounds in record['distances_to_center'])
        separation = np.mean(record['end_to_end_distance']) / (2*delta)
        ax.plot(distances, [i+.075]*2, color=color, marker='o', markersize=5.4, linewidth=1.55)
        ax.scatter(separation, i-.15, color=color, marker='D', s=27, clip_on=False)
        ax.text(distances[-1]+.015, i+.075, f'{distances[-1]:.3f}', color=color, fontsize=8, va='center')
    ax.set(xlim=(.72, 1.56), ylim=(3.52, -.52), xlabel='Distance relative to its noise allowance')
    ax.set_yticks(range(4), [row[2] for row in rows])
    ax.tick_params(axis='y', length=0, pad=7)
    ax.axhline(1.52, color=GRID, linewidth=.65)
    ax.spines['left'].set_visible(False)
    ax.grid(axis='x', color=GRID, linewidth=.65)
    ax.set_axisbelow(True)
    symbols = [Line2D([0, 1], [0, 0], linewidth=1.5, color=PALE, marker='o', markersize=4.8,
                      label=r'Distances of $F_-$ and $F_+$ from $F_0$, divided by $\delta$'),
               Line2D([0], [0], linestyle='none', color=PALE, marker='D', markersize=4.8,
                      label=r'Distance between ends, divided by $2\delta$')]
    ax.legend(handles=symbols, loc='upper left', bbox_to_anchor=(-.05, -.17),
              frameon=False, fontsize=8.1, handlelength=2.6, borderaxespad=0)


def specification():
    faces = read('faces/results.json')
    caps = read('bounds/refined_caps.json')
    fig = plt.figure(figsize=(12.1, 5.85))
    heading(fig, 'How much of the material specification remains possible?',
            r'Left-region mean contrast at the fixed reference $F_0$, with 10% noise and independently bounded materials.')
    seq = fig.add_gridspec(1, 2, left=.126, right=.955, bottom=.225, top=.719, wspace=.49, width_ratios=[1.2, 1])
    heatmap, points = [fig.add_subplot(seq[0, i]) for i in (0, 1)]
    margin_points(points, contents_grid(heatmap, faces, caps))
    fig.text(.126, .845, 'a   Certified task diameter as a fraction of its prior width', weight='bold', fontsize=10.8)
    fig.text(.126, .805, 'Material half-width 0.75. Ranges give certified lower and upper limits.', color=PALE, fontsize=8.8)
    fig.text(.655, .845, 'b   Why material outside the region matters', weight='bold', fontsize=10.8)
    fig.text(.655, .805, '8-slot low-contrast case. Both arms use 6 complex readings', color=PALE, fontsize=8.8)
    fig.text(.655, .777, 'and equal incident exposure in the homogeneous calibration.', color=PALE, fontsize=8.8)
    center_ok = faces['summary']['center_ball_both']
    complete_small = sum(row['center_ball_both'] for row in faces['rows'] if row['width'] < .75)
    fig.text(.126, .139, f'Full prior diameter in {center_ok}/72 settings, including all {complete_small} tests',
             color=INK, fontsize=9.4, weight='bold')
    fig.text(.126, .107, 'with the two smaller half-widths, 0.05 and 0.25.', color=PALE, fontsize=9.2)
    fig.text(.126, .055, 'Below 1 at both circle markers means both materials fit the same reference. '
             'A diamond above 1 rules out any common reading for that particular pair.', color=PALE, fontsize=8.8)
    output(fig, 'specification_range')


if __name__ == '__main__':
    multiplicity()
    specification()
