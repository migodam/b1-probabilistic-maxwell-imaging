#!/usr/bin/env python3
"""Plot the sealed reconstruction path and interval-audited readout."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from run import CASES, HERE, RING, read_json

STYLES = (
    ('Archived ring · 24 × 24', .96),
    ('Native slotted block · 96 × 96', 2.),
)
GN = '#295681'
PROOF = '#97462d'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r10-reference', type=Path, required=True,
                        help='Original R10 synthetic reference, available only in an authorized workspace.')
    args = parser.parse_args()
    metrics = read_json(HERE / 'interval_and_offline.json')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.4,
                         'mathtext.fontset': 'dejavusans', 'pdf.fonttype': 42,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axs = plt.subplots(2, 3, figsize=(6.02, 4.65),
                            gridspec_kw={'width_ratios': [1., 1., 2.12]})
    fig.subplots_adjust(left=.045, right=.98, top=.885, bottom=.15, wspace=.26, hspace=.75)

    for idx, (case, (title, width)) in enumerate(zip(CASES, STYLES, strict=True)):
        ax_ref, ax_gn, ax_path = axs[idx]
        case_result = metrics[case]
        online = read_json(HERE / f'{case}_online.json')
        last = case_result['endpoint_idx']
        with np.load(HERE / f'{case}_path.npz', allow_pickle=False) as seal:
            endpoint = seal['materials'][last]
        reference_file = args.r10_reference if case == RING else HERE / f'{case}_offline_reference.npz'
        with np.load(reference_file, allow_pickle=False) as reference:
            truth = reference['m']
        n = int(np.sqrt(truth.size))
        extent = [-width/2, width/2, -width/2, width/2]
        for ax, arr, label in ((ax_ref, truth, 'Reference'), (ax_gn, endpoint, 'Terminal iterate')):
            img = ax.imshow(arr.reshape(n, n), extent=extent, origin='lower',
                            interpolation='nearest', cmap='cividis', vmin=0, vmax=truth.max())
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(label, fontsize=8.4, pad=3)
            for spine in ax.spines.values():
                spine.set_visible(False)
        pos1, pos2 = ax_ref.get_position(), ax_gn.get_position()
        fig.text(pos1.x0, pos1.y1 + .058, f'({chr(97 + 2*idx)})  {title}', fontsize=9.3, weight='bold')
        cax = fig.add_axes([pos1.x0, pos1.y0 - .028, pos2.x1 - pos1.x0, .012])
        cb = fig.colorbar(img, cax=cax, orientation='horizontal',
                          ticks=[0., truth.max()/2, truth.max()])
        cb.outline.set_visible(False)
        cb.ax.tick_params(labelsize=7.3, length=0, pad=2)
        cb.ax.set_xticklabels(['0', f'{truth.max()/2:.2g}', f'{truth.max():.2g}'])
        fig.text((pos1.x0+pos2.x1)/2, pos1.y0-.07, r'Contrast $\chi$', ha='center', fontsize=7.6)

        accepted = [online['trials'][i] for i in online['accepted_indices']]
        radius, cap = (case_result[key] for key in ('radius', 'cap'))
        floor = case_result['cap_wide_residual_floor'] / radius
        peaks = [r['peak'] for r in accepted]
        ratios = [r['residual_frobenius']/radius for r in accepted]
        ax_path.set_yscale('log')
        ax_path.set_xlim(-.018*max(cap, peaks[-1]), 1.09*max(cap, peaks[-1]))
        ax_path.set_ylim(.73, 1.55 * max(ratios))
        ax_path.set_axisbelow(True)
        ax_path.axvline(cap, color='.42', linestyle=':', lw=.85)
        ax_path.axhline(1, color='.35', linestyle='--', lw=.80)
        ax_path.plot([0, cap], [floor, floor], color=PROOF, lw=2)
        ax_path.fill_between([0, cap], [.73, .73], [floor, floor], color=PROOF, alpha=.11, linewidth=0.)
        ax_path.plot(peaks, ratios, color=GN, marker='o', lw=1.2, ms=3.3)
        ax_path.scatter(peaks[-1], ratios[-1], marker='D', edgecolors='white',
                        linewidths=.45, color=GN, s=38, zorder=4)
        ax_path.set_title(f'({chr(98 + 2*idx)})  Full-wave residual / noise radius', loc='left',
                          fontsize=8.6, weight='bold')
        ax_path.set_xlabel(r'Peak contrast $\max_i\,\chi_i$')
        if idx:
            ax_path.annotate('Update 2', (peaks[2], ratios[2]), xytext=(-11, -18), textcoords='offset points',
                             ha='right', fontsize=7.2, color=GN,
                             arrowprops={'arrowstyle': '-', 'color': GN, 'lw': .55})
        ax_path.annotate(f'cap {cap:g}', (cap, 1.85), xytext=(3 if idx else -3, 0),
                         textcoords='offset points', ha='left' if idx else 'right',
                         va='bottom', fontsize=7.1, color='.30')

    legend = [
        Line2D([], [], marker='o', lw=1.2, ms=3.3, color=GN, label='Ridge–GN path'),
        Line2D([], [], marker='D', lw=0, ms=4.9, color=GN, label='Terminal iterate'),
        Line2D([], [], lw=2, color=PROOF, label='Cap-wide lower bound'),
        Line2D([], [], lw=.8, ls='--', color='.35', label='Noise-ball boundary'),
    ]
    fig.legend(handles=legend, loc='lower center', bbox_to_anchor=(.5, .012), ncol=4,
               frameon=False, fontsize=8, columnspacing=1.2, handlelength=2.)
    for fmt in ('pdf', 'png'):
        options = {'dpi': 300}
        if fmt == 'pdf':
            options['metadata'] = {'CreationDate': None, 'Author': ''}
        fig.savefig(HERE / f'figure_gn_diagnostic.{fmt}', **options)
    plt.close(fig)


if __name__ == '__main__':
    main()
