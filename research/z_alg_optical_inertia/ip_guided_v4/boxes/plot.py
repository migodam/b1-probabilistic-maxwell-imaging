"""Plot only the saved numerical scans and their independently verified points."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from run import CONTROL_ID, HERE, LEGACY_ID, NATIVE_IDS, SMALL_ID, read

COLORS = {'positive': '#245c80', 'rank': '#af5038', 'neumann': '#3a806b'}


def split(rows):
    out = {}
    for kind, predicate in [('positive', lambda n: n == 0), ('rank', lambda n: n > 0)]:
        out[kind] = [max((r['distance_over_radius_numeric'] or 0. for r in b['by_mu']
                          if predicate(r['nu'])), default=np.nan) for b in rows]
    return out


def widths(path, out):
    scan = read(path / 'scan.json')['cases']
    proofs = read(path / 'verified.json')['cases']
    disk = read(path / 'followup_disk.json')['cases'][CONTROL_ID]
    host = read(path / 'followup_host.json')['cases']
    names = {LEGACY_ID: '24² · annulus', SMALL_ID: '32² · independent same-grid',
             NATIVE_IDS[0]: '96² · optfresh 01', NATIVE_IDS[1]: '96² · optfresh 02'}
    fig, axs = plt.subplots(2, 3, figsize=(13.0, 7.4), layout='constrained')
    legend = [Line2D([0], [0], color=COLORS[k], label=lab) for k, lab in
              [('positive', 'Positive-definite scan'), ('rank', 'Indefinite scan'),
               ('neumann', '2-norm of centered iteration')]]
    legend += [Line2D([0], [0], color='#282822', marker='o', markersize=5, linestyle='None',
                      label='Strict lower endpoint')]

    def axes(ax, title, rows, selected):
        tx = ax.twinx()
        x = [b.get('width', b.get('parameters', {}).get('width')) for b in rows]
        for k, y in split(rows).items():
            ax.plot(x, y, color=COLORS[k], linewidth=1.8, marker='.', markersize=4)
        tx.plot(x, [b['norm_T_numeric'] for b in rows], color=COLORS['neumann'],
                linewidth=1.3, linestyle='--')
        tx.axhline(1., color=COLORS['neumann'], linewidth=.8, alpha=.45)
        for row, kind in selected:
            if row.get('status') != 'outward_verified':
                continue
            w = row.get('width', row.get('parameters', {}).get('width'))
            ax.scatter(w, row['distance_over_radius_lower'], s=32, marker='o',
                       color=COLORS[kind], edgecolor='white', linewidth=.8, zorder=8)
        ax.axhline(1., color='#222222', linewidth=.75, linestyle=':', alpha=.88)
        ax.set(title=title, xlabel='Box half-width', ylabel='Residual lower endpoint / radius')
        ax.set_yscale('symlog', linthresh=1)
        ax.set_ylim(bottom=0.)
        tx.set(ylabel='Centered 2-norm', ylim=(0., max(1.25, 1.08*max(b['norm_T_numeric'] for b in rows))))
        tx.spines['right'].set_color(COLORS['neumann'])
        tx.tick_params(axis='y', colors=COLORS['neumann'])
        tx.yaxis.label.set_color(COLORS['neumann'])
        ax.spines[['top']].set_visible(False)
        tx.spines[['top']].set_visible(False)

    for ax, (key, title) in zip((axs[0,0], axs[0,1], axs[0,2], axs[1,0]), names.items(), strict=True):
        p = proofs[key]['tests']
        anchor = p['rank_tail']
        rows = [b for b in scan[key]['boxes'] if b['placement'] == anchor['placement'] and
                b['radius_fraction'] == anchor['radius_fraction'] and b['peak'] == anchor['peak']]
        selected = [(v, 'rank' if label == 'rank_tail' else 'positive') for label, v in p.items()]
        axes(ax, title, rows, selected)
    disk_selected = [(value, 'rank') for value in disk['strict_checks'].values()]
    axes(axs[1,1], '96² · disk data, high box on same support', disk['proposals'], disk_selected)

    ax = axs[1,2]
    x = np.arange(len(NATIVE_IDS))
    for j, (label, text) in enumerate([('rank_tail_identity', 'Identity weight'),
                                       ('rank_tail_inverse_center_illumination', 'Center field weight')]):
        strict = [host[key]['strict_checks'][label]['distance_over_radius_lower'] for key in NATIVE_IDS]
        ax.bar(x + (.17 if j else -.17), strict, width=.30, color=['#af5038', '#d89570'][j], label=text)
    ax.axhline(1., color='#222222', linewidth=.75, linestyle=':', alpha=.88)
    ax.set(xticks=x, xticklabels=['optfresh 01', 'optfresh 02'], ylabel='Residual lower endpoint / radius',
           title='96² · paired nonzero fixed host hypotheses', ylim=(0, max(35, *strict)))
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(loc='upper right', frameon=False, fontsize=8)
    fig.legend(handles=legend, loc='outside upper center', ncol=4, frameon=False)
    fig.savefig(out / 'box-widths.png', dpi=190)
    fig.savefig(out / 'box-widths.pdf')
    plt.close(fig)


def maps(path, out):
    proof = read(path / 'verified.json')['cases']
    keys = ((NATIVE_IDS[0], 'rank_tail'), (NATIVE_IDS[1], 'rank_tail'))
    fig, axs = plt.subplots(2, 4, figsize=(12, 6.6), layout='constrained')
    with np.load(path / 'verified_witnesses.npz', allow_pickle=False) as states,             np.load(path / 'scan_maps.npz', allow_pickle=False) as obs:
        for i, (key, label) in enumerate(keys):
            prefix = key + '__' + label + '__'
            j = proof[key]['side']
            region = states[prefix + 'active']
            y = obs[key + '__Y']
            center_port = states[prefix + 'center_ports']
            img = [np.zeros(j*j) for _ in range(2)]
            for k, field in enumerate(('lower', 'upper')):
                img[k][region] = states[prefix + field]
                m = axs[i,k].imshow(img[k].reshape(j,j), vmin=0., vmax=6., cmap='magma',
                                    extent=(-1,1,1,-1))
                axs[i,k].set(xlabel='x (m)', ylabel='y (m)', title=f'optfresh {i+1:02d} · {field} map')
            fig.colorbar(m, ax=list(axs[i,:2]), fraction=.04, label='Contrast')
            v = max(np.max(abs(y)), np.max(abs(center_port)))
            for k, (matrix, title) in enumerate(((y, 'Observed'), (center_port, 'Center forward'))):
                im = axs[i,k+2].imshow(abs(matrix), vmin=0., vmax=v, cmap='viridis', aspect='auto')
                axs[i,k+2].set(xlabel='Source', ylabel='Receiver', title=f'{title} ports')
            fig.colorbar(im, ax=list(axs[i,2:]), fraction=.04, label='Magnitude')
    fig.savefig(out / 'native-boxes.png', dpi=190)
    fig.savefig(out / 'native-boxes.pdf')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--results', type=Path, default=HERE / 'results')
    ap.add_argument('--out', type=Path, default=HERE / 'figures')
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    with plt.rc_context({'font.family': 'DejaVu Sans', 'font.size': 10}):
        widths(args.results, args.out)
        maps(args.results, args.out)


if __name__ == '__main__':
    main()
