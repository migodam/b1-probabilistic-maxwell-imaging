"""Draw the fixed-layout and fixed-receiver comparisons from the saved results."""

import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


HERE = Path(__file__).resolve().parent
LAYOUT_COLORS = ('#236f94', '#b35436')
FREQUENCY_COLORS = ('#676eaa', '#236f94', '#13776a', '#b35436')
ANGLES = ('$x$', '45°', '$y$')


def layout_panels(radius, cross, cases):
    for count, color, offset in zip((2, 3), LAYOUT_COLORS, (-.13, .13), strict=True):
        rows = [c['readings'] for c in cases if c['components'] == count]
        x = np.arange(3)+offset
        bounds = 100*np.array([row['covering_radius_fraction_shared_x_reference_bounds'] for row in rows])
        near = 100*np.array([row['nearest_pair_radius_fraction_shared_x_reference'] for row in rows])
        radius.errorbar(x, bounds.mean(axis=1), yerr=(bounds[:, 1]-bounds[:, 0])/2, fmt='o',
                        markersize=4, color=color, capsize=5, elinewidth=2, zorder=3)
        radius.scatter(x, near, marker='v', facecolor='white', edgecolor=color, linewidth=1.2, zorder=4)
        norms = np.array([row['cross_reference_L2'] for row in rows])
        cross.plot(x-offset, 1000*norms, 'o-', color=color, linewidth=1.55, label=f'{count} components')
        for pos, row in zip(x-offset, rows, strict=True):
            values = list(row['cross_L2_by_branch'].values())
            cross.vlines(pos, 1000*min(values), 1000*max(values), color=color, alpha=.28, linewidth=4)
    radius.set_title('a   Full-data radius for the coupled roots', loc='left', fontweight='bold')
    radius.set_ylabel('Radius / shared $x$-layout signal (%)')
    radius.axhline(10, color='.46', linestyle=':', linewidth=1.15)
    radius.text(2.15, 11.4, '10%', ha='right', color='.37')
    radius.set_yscale('log')
    radius.set_ylim(.08, 50)
    radius.legend(handles=[Line2D([], [], marker='o', color='.25', lw=0, label='Cover all roots (bounds)'),
                           Line2D([], [], marker='v', color='.25', markerfacecolor='white', lw=0,
                                  label='Cover nearest pair')], loc='lower left', frameon=False)
    cross.set_title('b   Cross-port signal across the layouts', loc='left', fontweight='bold')
    cross.set_ylabel(r'Cross-port norm $\times\,10^3$')
    cross.legend(loc='upper left', frameon=False)
    cross.text(.975, .028, 'Points: reference    Bars: range across roots', transform=cross.transAxes,
               ha='right', color='.34', fontsize=8.3)
    for ax in (radius, cross):
        ax.set_xlim(-.3, 2.3)
        ax.set_xticks(range(3), ANGLES)
        ax.set_xlabel(r'Translation direction, spacing $30\lambda_0$')


def aperture_panels(single, symmetry, cases):
    h = np.array([.01, .02, .04, .10, .20])
    for ratio, color in zip((.5, 1., 2., 4.), FREQUENCY_COLORS, strict=True):
        rows = [row for c in cases if c['frequency_ratio'] == ratio for row in c['pairs']]
        single_values = np.concatenate([np.array([row['joint_half_gap_over_fixed_10_percent_radius'][1:6],
                                                 row['joint_half_gap_over_fixed_10_percent_radius'][6:]])
                                        for row in rows])
        symmetry_values = [[v['min_even_odd_over_clean_radius_for_unknown_gain']
                            for v in row['symmetric_receiver_pair']] for row in rows]
        for ax, values in ((single, single_values), (symmetry, symmetry_values)):
            lo, hi = np.min(values, axis=0), np.max(values, axis=0)
            ax.fill_between(h, lo, hi, color=color, alpha=.18)
            ax.plot(h, (lo+hi)/2, 'o-', color=color, markersize=3.6, linewidth=1.5, label=f'{ratio:g}')
    calibrated = [row['symmetric_receiver_pair'][-1]['odd_over_clean_radius_for_calibrated_gain'] for row in rows]
    symmetry.scatter(h[-1], np.mean(calibrated), marker='D', facecolor='white', edgecolor=FREQUENCY_COLORS[-1],
                     s=32, zorder=5)
    single.set_title('c   Move one receiving port, fixed scene', loc='left', fontweight='bold')
    single.set_ylabel(r'Mirror-pair gap / $(2\delta)$')
    symmetry.set_title('d   Symmetric receiver pair, unknown gain', loc='left', fontweight='bold')
    symmetry.set_ylabel(r'$\min(\Vert z_{\mathrm{even}}\Vert,\,\Vert z_{\mathrm{odd}}\Vert)/\delta$')
    single.text(.975, .026, r'On axis: field gap $<3\times10^{-17}$ at all four frequencies',
                transform=single.transAxes, ha='right', fontsize=8.3, color='.34')
    symmetry.text(.975, .026, 'Open diamond: calibrated gain at the last 4× point',
                  transform=symmetry.transAxes, ha='right', fontsize=8.3, color='.34')
    for ax in (single, symmetry):
        ax.set(xscale='log', yscale='log', xlim=(.0085, .235), ylim=(.025, 9.5))
        ax.set_xticks(h, [f'{value:.2f}' for value in h])
        ax.set_xlabel(r'Lateral offset magnitude $h/\lambda_0$')
        ax.axhline(1., color='.46', linestyle=':', linewidth=1.15)
        ax.axvline(.01, color='.92', linewidth=.6, zorder=0)


def main():
    data = json.loads((HERE/'results.json').read_text())
    with matplotlib.rc_context({'font.family': 'DejaVu Sans', 'font.size': 9.5,
                                'axes.spines.right': False, 'axes.spines.top': False,
                                'pdf.fonttype': 42, 'axes.axisbelow': True}):
        fig, ax = plt.subplots(2, 2, figsize=(10.85, 7.45))
        fig.subplots_adjust(left=.092, right=.974, bottom=.134, top=.964, wspace=.28, hspace=.40)
        layout_panels(*ax[0], data['layout_trials'][-1]['cases'])
        aperture_panels(*ax[1], data['mirror_trials'][-1]['cases'])
        for panel in ax.flat:
            panel.grid(axis='y', color='.9', linewidth=.65)
        handles, labels = ax[1, 1].get_legend_handles_labels()
        fig.legend(handles, labels, title='Frequency ratio $k/k_0$  (bands across mirror pairs)', ncol=4,
                   loc='lower center', bbox_to_anchor=(.5, .009), frameon=False)
        fig.savefig(HERE/'finite_aperture.png', dpi=230)
        fig.savefig(HERE/'finite_aperture.pdf')
        plt.close(fig)


if __name__ == '__main__':
    main()
