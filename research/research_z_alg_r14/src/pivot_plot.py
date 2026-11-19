"""Plot equal-cost material-region contraction and matched-quality CPU times."""
from collections import Counter
import json

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
import numpy as np

from data14 import ROOT
from pivot_experiment import METHODS

LABEL = ('R13 pivot', 'R10 scalar', 'R10 shared', 'Direct 6D')
COLORS = ('#26867f', '#255d9b', '#8965a9', '#d16a35')
COHORTS = ('spectral', 'natural', 'aperture')
HEADINGS = ('R12 conditional old reads', 'R10 natural six-cell', 'R10 finite aperture')


def main():
    data = json.loads((ROOT / 'results/pivot_old.json').read_text())
    cost = json.loads((ROOT / 'results/pivot_cost.json').read_text())
    fig, axs = plt.subplots(2, 3, figsize=(10.6, 5.9), sharex=True,
                             gridspec_kw=dict(wspace=.18, hspace=.28))
    rng = np.random.default_rng(71)
    for column, (group, title) in enumerate(zip(COHORTS, HEADINGS)):
        cases = [x['methods'] for x in data['cases'].values() if x['cohort'] == group]
        axs[0, column].set_title(f'{title}  ·  n={len(cases)}', fontsize=9.5, loc='left', pad=8)
        for row, metric in enumerate(('unresolved_relative_volume', 'mean_relative_hull_width')):
            ax = axs[row, column]
            if row == 0:
                ax.set_yscale('log')
                ax.set_ylim(1e-11, 2)
                ax.set_yticks([1e-11, 1e-9, 1e-7, 1e-5, 1e-3, 1e-1, 1.])
            else:
                ax.set_ylim(0., 1.10)
                ax.set_yticks([0, .25, .5, .75, 1.])
                ax.yaxis.set_major_formatter(PercentFormatter(1.))
            for i, (method, color) in enumerate(zip(METHODS, COLORS)):
                value = [x[method][metric] for x in cases]
                xx = rng.uniform(-.15, .15, len(value)) + i
                ax.scatter(xx, value, s=12, alpha=.46, color=color, edgecolors='none', zorder=2)
                q1, med, q3 = np.quantile(value, [.25, .5, .75])
                ax.plot([i, i], [q1, q3], color=color, lw=4., solid_capstyle='round', zorder=3)
                ax.plot([i-.17, i+.17], [med, med], color='#0b222f', lw=1.6, zorder=4)
            ax.spines[['right', 'top']].set_visible(False)
            ax.grid(axis='y', color='#dfe5ea', linewidth=.55)
            ax.set_axisbelow(True)
            ax.set_xlim(-.55, 3.55)
            ax.set_xticks(range(len(LABEL)), LABEL, fontsize=8.2)
            ax.tick_params(axis='x', length=0, pad=5)
            if column:
                ax.tick_params(axis='y', labelleft=False)
            else:
                ax.tick_params(axis='y', labelsize=8)
    axs[0, 0].set_ylabel('Relative volume retained', fontsize=9.5)
    axs[1, 0].set_ylabel('Mean global interval width', fontsize=9.5)
    fig.subplots_adjust(left=.090, right=.98, top=.91, bottom=.12)
    fig.text(.093, .025, f"Starting prior: six axes [0.5, 30]. Equal CPU allowance: {data['cpu_cap_per_case_method_s']:.2f} s/case/method. Bar = median; stem = middle 50%.",
             color='#3a4d5e', fontsize=8.4)
    for suffix in ('pdf', 'png'):
        fig.savefig(ROOT / 'figures' / f'pivot_shrink.{suffix}', dpi=220, bbox_inches='tight')
    plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(10.6, 3.2), sharey=True)
    thresholds = (.5, .1, .01, .001)
    cap = cost['cpu_cap_per_case_method_s']
    for column, (group, title) in enumerate(zip(COHORTS, HEADINGS)):
        ax = axs[column]
        cases = [x['methods'] for x in cost['cases'].values() if x['cohort'] == group]
        for i, (method, color) in enumerate(zip(METHODS, COLORS)):
            for j, rows in enumerate(cases):
                attained = rows[method]['cpu_to_unresolved_volume_upper_bound']
                for index, target in enumerate(thresholds):
                    hit = attained.get(str(target))
                    xx = index + (i - 1.5) * .16 + (j - .5) * .05
                    ax.scatter(xx, hit['cpu_s'] if hit else cap, s=19, color=color,
                               marker='o' if hit else 'v', alpha=.87, linewidths=0, zorder=3)
        ax.set_xticks(range(4), ['50%', '10%', '1%', '0.1%'])
        ax.set_xlabel('Relative volume target', fontsize=9)
        ax.set_title(title, fontsize=9.2, loc='left')
        ax.spines[['right', 'top']].set_visible(False)
        ax.axhline(cap, color='#d4d9df', lw=.8)
        ax.grid(axis='y', color='#e2e6e9', linewidth=.55)
        ax.set_ylim(0., cap * 1.103)
        ax.tick_params(labelsize=8.3)
    axs[0].set_ylabel('CPU seconds to target', fontsize=9)
    handles = [Line2D([0], [0], color=c, marker='o', markersize=5, linestyle='none', label=l)
               for l, c in zip(LABEL, COLORS)]
    handles.append(Line2D([0], [0], color='#3a4d5e', marker='v', markersize=5, linestyle='none', label='Budget cap, target not reached'))
    fig.legend(handles=handles, loc='lower center', ncol=5, fontsize=8.4, frameon=False)
    fig.subplots_adjust(left=.07, right=.98, top=.86, bottom=.23, wspace=.13)
    for suffix in ('pdf', 'png'):
        fig.savefig(ROOT / 'figures' / f'pivot_cpu_to_quality.{suffix}', dpi=220, bbox_inches='tight')
    plt.close(fig)

    for group in COHORTS:
        cases = [x['methods'] for x in data['cases'].values() if x['cohort'] == group]
        print(group, len(cases))
        for method in METHODS:
            arrays = [x[method] for x in cases]
            stats = [float(np.median([x[metric] for x in arrays]))
                     for metric in ('unresolved_relative_volume', 'mean_relative_hull_width', 'cpu_s',
                                    'setup_cpu_s', 'evaluated_boxes', 'retained_boxes')]
            events = Counter()
            for x in arrays:
                events.update(x['decisions'])
            print(method, stats, dict(events))


if __name__ == '__main__':
    main()
