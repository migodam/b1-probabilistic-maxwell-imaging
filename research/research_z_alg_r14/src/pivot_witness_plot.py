"""Plot certified diameter lower bounds and their split-band upper comparison."""
import json

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from data14 import ROOT
from pivot_witnesses import ARMS

GROUPS = ('spectral', 'natural', 'aperture')
LABELS = ('R12 conditional (36)', 'R10 natural (18)', 'R10 finite aperture (8)')
COLORS = ('#255d9b', '#cc673c', '#278875')
ARM_LABELS = ('Repeat', '−2%', '+2%', '±2%', '±8%')


def main():
    base = ROOT / 'results'
    cases = json.loads((base / 'pivot_diameter_witnesses.json').read_text())['cases']
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10.5, 4.1), sharey=True,
                                 gridspec_kw=dict(width_ratios=[.86, 1.14]))
    rng = np.random.default_rng(1448)
    next_case = 0
    boundaries, intervals = [], []
    for j, (cohort, color) in enumerate(zip(GROUPS, COLORS)):
        group = [(name, item) for name, item in cases.items() if item['cohort'] == cohort]
        for i, act in enumerate(ARMS):
            found = [q['actions'][act]['witness']['rms_lower'] for _, q in group
                     if q['actions'][act]['witness'] is not None]
            pos = i + .23 * (j - 1)
            ax.scatter(pos + rng.uniform(-.08, .08, len(found)), found,
                       color=color, s=11, alpha=.42, linewidths=0, zorder=2)
            ax.scatter(pos, np.median(found), facecolor=color, edgecolors='white',
                       linewidths=.60, marker='D', s=27, zorder=3)
        arm = 'split_0.02'
        limits = []
        for name, item in group:
            witness = item['actions'][arm]['witness']
            if witness is None:
                continue
            upper = json.loads((base / f'frequency_covers/{name}.json').read_text())['actions'][arm]['D_upper']
            limits.append((witness['rms_lower'], upper))
        lo, hi = np.array(sorted(limits)).T
        x = np.arange(next_case, next_case + len(lo))
        bx.vlines(x, lo, hi, color=color, linewidth=.66, alpha=.68)
        bx.scatter(x, lo, s=11, color=color, alpha=.96, linewidths=0, zorder=3)
        bx.scatter(x, hi, s=10, facecolors='white', edgecolors=color, linewidths=.6, zorder=3)
        intervals.append(float(x.mean()))
        boundaries.append(float(x[-1] + 1.5))
        next_case = int(x[-1] + 4)
    for panel in (ax, bx):
        panel.set_yscale('log')
        panel.set_ylim(.0017, 35.)
        panel.grid(axis='y', color='#e0e5e9', lw=.55)
        panel.spines[['top', 'right']].set_visible(False)
        panel.tick_params(labelsize=8.6)
    ax.set_title('A   Certified diameter lower bounds', loc='left', fontsize=10.5, pad=9)
    ax.set_ylabel('Material difference (RMS)', fontsize=9.5)
    ax.set_xticks(range(len(ARM_LABELS)), ARM_LABELS)
    ax.set_xlabel('New measurement, equal total energy', fontsize=9.2)
    bx.set_title('B   Diameter bounds for the ±2% split', loc='left', fontsize=10.5, pad=9)
    for position in boundaries[:-1]:
        bx.axvline(position, color='#d9dfe3', lw=.7, linestyle=':')
    bx.set_xticks(intervals, ['Conditional', 'Natural', 'Finite aperture'])
    bx.set_xlabel('Cases ordered by lower bound within each group', fontsize=9.2)
    handles = [Line2D([0], [0], color=c, marker='D', markersize=4.7, linestyle='none', label=label)
               for c, label in zip(COLORS, LABELS)]
    handles.extend([Line2D([0], [0], marker='o', markersize=5, linestyle='none', color='#364856', label='Lower bound'),
                    Line2D([0], [0], marker='o', markersize=5, linestyle='none', markerfacecolor='white',
                           color='#364856', label='Upper bound')])
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .016), ncol=5,
               fontsize=8.1, frameon=False)
    fig.subplots_adjust(left=.077, right=.98, top=.87, bottom=.21, wspace=.12)
    for suffix in ('pdf', 'png'):
        fig.savefig(ROOT / 'figures' / f'pivot_diameter_witnesses.{suffix}', dpi=220, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
