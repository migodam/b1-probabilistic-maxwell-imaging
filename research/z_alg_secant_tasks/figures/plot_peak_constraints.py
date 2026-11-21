"""Plot the archived peak-cap bounds and detector-geometry transition."""

import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
import numpy as np


HERE = Path(__file__).resolve().parent
DATA = HERE.parent / 'highdim'
RING = 'r10_twod_0001'
BLUE, TEAL, VIOLET, RUST = '#1e6091', '#007c71', '#7054a5', '#b44726'


def read_json(path):
    return json.loads((DATA/path).read_text(encoding='utf-8'))


def peak_bounds(ax, geometry):
    positive = read_json('sector/results.json')['cases'][RING]
    inertia = read_json('sector_inertia/results.json')['cases'][RING]['grid']
    extension = read_json('sector_inertia/ring_cap_extension.json')['cases'][RING]['grid']
    delta = positive['archived_radii']['frozen_observation_ball']
    peak_grid = positive['grid']

    def plot(rows, read_distance, color, marker, label):
        distance = np.array([read_distance(row) for row in rows], dtype=float)/delta
        ax.plot([row['u'] for row in rows], distance, color=color, marker=marker,
                markersize=5.8, markerfacecolor='white', markeredgewidth=1.45,
                linewidth=1.55, label=label, zorder=3)

    plot(peak_grid, lambda r: r['controls']['direction_free_sector_amplitude'].get('distance_lower_using_receiver_amplitude'),
         TEAL, 'o', 'Received amplitude only')
    plot(peak_grid, lambda r: r['best_grid_result']['distance_lower_bound'],
         BLUE, 's', 'PD: trace')
    plot(peak_grid, lambda r: r['best_single_coherent_result']['best_single_coherent_distance_lower_bound'],
         VIOLET, '^', 'PD: best coherent mode')
    plot(inertia+extension,
         lambda row: (row['best_indefinite_with_source_directions'] or {}).get('lower_bound_numeric'),
         RUST, 'D', 'Indefinite inertia bound')

    cert = positive['independent_certificate']
    if cert['strictly_excludes']['frozen_observation_ball']:
        ax.scatter(cert['u'], cert['distance_lower_bound']/delta, color=BLUE, marker='s', s=29, zorder=5)
    cert = positive['positive_gram_projection_certificate']
    rank_one = cert['single_coherent_projection']
    if rank_one['strictly_excludes']['frozen_observation_ball']:
        ax.scatter(cert['u'], rank_one['distance_lower_bound']/delta, color=VIOLET, marker='^', s=31, zorder=5)
    for name in ('certificate.json', 'certificate_u5.json'):
        path = DATA/'sector_inertia'/name
        if path.exists():
            record = read_json(path.relative_to(DATA))
            if record['status'] == 'Arb_verified' and record['strictly_excludes']['frozen_observation_ball']:
                ax.scatter(record['u'], record['certified_frobenius_residual_lower']/delta,
                           color=RUST, marker='D', s=29, zorder=5)

    limit = geometry[RING]['receiver_layouts'][-1]['all_measured_rows_projected']['critical_material_u_if_positive']
    ax.axhline(1., color='.46', linewidth=1.1, linestyle='--', zorder=1)
    ax.axvline(limit, color='#41484c', linewidth=1.1, linestyle=':', zorder=1)
    ax.annotate(f'Global-sector limit\n$u_*={limit:.5f}$', (limit, 65.), (limit-.23, 65.),
                fontsize=8, ha='right', va='center', color='#41484c',
                arrowprops=dict(arrowstyle='-', color='#41484c', linewidth=.7))
    ax.text(.055, 1.27, 'Exclusion threshold', fontsize=8, color='.36')

    residuals = {row['u']: row['best_indefinite_with_source_directions']['lower_bound_numeric']/delta
                 for row in extension if row['best_indefinite_with_source_directions']}
    ax.annotate('5.5 and 6:\nthreshold not met', (5.5, residuals[5.5]), (5.05, .34), fontsize=8,
                color=RUST, ha='center', va='bottom',
                arrowprops=dict(arrowstyle='-', linewidth=.7, color=RUST))

    ax.set_title('a   Peak-cap exclusion on the 24² ring', loc='left', fontweight='bold')
    ax.set(xlabel='Tested peak cap $u$', ylabel=r'Distance lower bound / original noise radius $\delta$',
           xlim=(-.04, 6.22), ylim=(-.08, 135), yscale='symlog')
    ax.set_yscale('symlog', linthresh=1., linscale=.85)
    ax.yaxis.set_major_locator(FixedLocator([0, .5, 1, 3, 10, 30, 100]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:g}'))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xticks(range(7))
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Line2D([], [], color='.31', marker='o', linewidth=0, markersize=5.4))
    labels.append('Filled symbol: Arb verified')
    ax.legend(handles, labels, loc='upper right', bbox_to_anchor=(1., 1.01), frameon=False,
              handlelength=1.65, fontsize=8.1, labelspacing=.38, borderpad=.45)


def receiver_geometry(ax, geometry):
    plotted = (
        (RING, BLUE, 'o', '-', 'R10 ring', (8, 0)),
        ('austria_eps4_400MHz', RUST, 's', '-', 'Austria', (8, 10)),
        ('s2d_007__pw360_10x20', TEAL, 'D', '--', 's007 full', (8, -13)),
        ('s2d_007__pw180_10x20', VIOLET, '^', ':', 's007/s010 half', (8, 0)),
    )
    for case, color, marker, style, label, offset in plotted:
        row = geometry[case]
        values = [row['unconstrained_receiver_bound_u0']['critical_material_u_if_positive']]
        values.extend(record['all_measured_rows_projected']['critical_material_u_if_positive']
                      for record in row['receiver_layouts'])
        ax.plot(range(4), values, color=color, linestyle=style, marker=marker, linewidth=1.55,
                markersize=5.8, markerfacecolor='white', markeredgewidth=1.3, zorder=3)
        count = row['receiver_layouts'][-1]['receivers']
        ax.annotate(f'{label}\n{values[-1]:.3f} · {count} Rx', (3, values[-1]), offset,
                    textcoords='offset points', ha='left', va='center', color=color, fontsize=8)

    ax.set_title('b   Geometry limit of the global sector', loc='left', fontweight='bold')
    ax.set(xlabel='Archived receiver rows retained', ylabel='Critical cap $u_*$',
           xlim=(-.14, 4.17), ylim=(.29, 5.23), yscale='log')
    ax.set_xticks(range(4), ['None', 'First', 'First half', 'All'])
    ax.yaxis.set_major_locator(FixedLocator([.3, .5, 1, 2, 4]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:g}'))
    ax.yaxis.set_minor_formatter(NullFormatter())


def main():
    geometry = {row['id']: row for row in read_json('sector_geometry/results.json')['rows']}
    with matplotlib.rc_context({'font.size': 9.6, 'font.family': 'DejaVu Sans', 'pdf.fonttype': 42,
                                'axes.spines.top': False, 'axes.spines.right': False, 'axes.axisbelow': True}):
        fig, ax = plt.subplots(1, 2, figsize=(12.4, 4.65), gridspec_kw={'width_ratios': [1.36, 1.]})
        fig.subplots_adjust(left=.087, right=.983, top=.915, bottom=.19, wspace=.26)
        peak_bounds(ax[0], geometry)
        receiver_geometry(ax[1], geometry)
        for panel in ax:
            panel.grid(axis='y', color='.89', linewidth=.68)
        fig.text(.09, .065, 'a  Zero denotes a zero lower bound; missing points denote unavailable bounds.',
                 color='.35', fontsize=8.15)
        fig.text(.09, .035, r'     R10 original $\delta=0.01051718$; stored $\mu$ grid. PD = positive definite.',
                 color='.35', fontsize=8.15)
        fig.text(.62, .065, 'b  Rx = receivers. Both SRAP half-ring limits overlap.',
                 color='.35', fontsize=8.15)
        fig.text(.62, .035, '     Frozen arrays and archived receive order; dimensionless caps.',
                 color='.35', fontsize=8.15)
        fig.savefig(HERE/'peak_constraints.png', dpi=250)
        fig.savefig(HERE/'peak_constraints.pdf')
        plt.close(fig)


if __name__ == '__main__':
    main()
