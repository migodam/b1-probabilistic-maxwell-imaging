"""Plot all six fine-grid objectives against their inflated observation radii."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'results/theory'
LABELS = ('01  Slotted block', '02  Four lobes', '03  Two chevrons',
          '04  Crescent & plate', '05  Offset channels', '06  Graded pinwheel')
VIEWS = (('180', .5, 'Half aperture · contrast cap 0.50'),
         ('360', .75, 'Full aperture · contrast cap 0.75'))
NOISE = (('02', '#087773', 'o', -.15), ('10', '#b04d27', 's', .15))
SPATIAL_COLOR = '#7e3095'


def read(name):
    return json.loads((DATA / f'{name}.json').read_text(encoding='utf-8'))


def group_max(rows, cap_name, value_name):
    values = {}
    for row in rows:
        key = (row['id'], row[cap_name])
        values[key] = max(values.get(key, -np.inf), row[value_name])
    return values


def main():
    numerical = group_max(read('fine_scan')['rows'], 'u', 'd')
    strict = group_max(read('fine_verified')['cases'], 'cap', 'distance_lower')
    radius = {row['id']: row['inflated_radius'] for row in read('fine_feasibility')['observations']}
    spatial = read('fine_spatial_verified')
    spatial_case = spatial['case']
    spatial_ratio = spatial['outward_distance_lower'] / radius[spatial_case]
    if not spatial['excluded_at_padded_radius'] or spatial_ratio <= 1.:
        raise ValueError('The spatial receipt does not exclude its cap at the inflated radius.')

    fig, axes = plt.subplots(1, 2, figsize=(11.1, 4.65), sharex=True, sharey=True,
                             layout='constrained')
    y = np.arange(len(LABELS))
    for ax, (aperture, cap, title) in zip(axes, VIEWS, strict=True):
        ax.set_axisbelow(True)
        ax.grid(axis='x', alpha=.2)
        for noise, color, marker, offset in NOISE:
            for i in range(len(LABELS)):
                case = f'optfresh_{i + 1:02d}__pw{aperture}_10x20__n{noise}'
                key = (case, cap)
                verified = key in strict
                scalar = strict[key] if verified else numerical[key]
                ratio = scalar / radius[case]
                pos = y[i] + offset
                ax.plot([0, ratio], [pos, pos], color='#d8dddb', lw=1.3, zorder=1)
                ax.scatter(ratio, pos, marker=marker, s=55, edgecolors=color,
                           facecolors=color if verified else 'white', linewidths=1.5, zorder=3)
                if case == spatial_case and cap == spatial['cap']:
                    ax.plot([ratio, spatial_ratio], [pos, pos], color=SPATIAL_COLOR,
                            lw=1.15, linestyle=':', zorder=2)
                    ax.scatter(spatial_ratio, pos, marker='*', s=180, color=SPATIAL_COLOR,
                               edgecolors='white', linewidths=.75, zorder=4)
        ax.axvline(1, color='#30383a', linestyle='--', lw=1)
        ax.set_yticks(y, LABELS)
        ax.set_xlim(left=-.13)
        ax.set_title(title, fontsize=11.5, loc='left', fontweight='bold')
        ax.set_xlabel('Residual lower bound ÷ inflated radius')
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].invert_yaxis()
    legends = [Line2D([], [], linestyle='none', marker=marker, markersize=6.5,
                     markerfacecolor=color, markeredgecolor=color, label=f'Scalar · {int(noise)}% noise')
               for noise, color, marker, _ in NOISE]
    legends += [Line2D([], [], linestyle='none', marker='o', markersize=6.5, markerfacecolor='white',
                       markeredgecolor='#30383a', label='Hollow = numerical'),
                Line2D([], [], linestyle='none', marker='*', markersize=11, color=SPATIAL_COLOR,
                       label='Spatial certificate (#03, 10%)'),
                Line2D([], [], ls='--', color='#30383a', label='Exclusion threshold')]
    fig.legend(handles=legends, loc='outside lower center', ncol=5, frameon=False, fontsize=8.7,
               columnspacing=1.3, handletextpad=.5)
    fig.suptitle('Peak contrast on the 96 × 96 forward grid', fontsize=14, fontweight='bold')
    for extension in ('png', 'pdf'):
        fig.savefig(ROOT / 'figures' / f'fine_peaks.{extension}', dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    main()
