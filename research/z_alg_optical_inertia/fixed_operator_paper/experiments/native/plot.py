"""Single-page figures from original material, source-index and scaling records."""
from pathlib import Path
import csv
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[2]
RING = STUDY/'ip_final_v5/strong_control'
V4 = STUDY/'ip_guided_v4'
LOW_ID = 'control_disk075__pw360_32x32__pair2pct'
DARK, MUTED, GRID = '#172733', '#445867', '#d8e3e8'
BLUE, ORANGE, GREEN = '#17638c', '#ac502c', '#1e776b'
LIGHT = '#dae2e6'


def read(path):
    return json.loads(path.read_text())


def save(fig, name):
    fig.savefig(HERE/f'{name}.pdf', bbox_inches='tight', pad_inches=.05)
    fig.savefig(HERE/f'{name}.svg', bbox_inches='tight', pad_inches=.05)
    fig.savefig(HERE/f'{name}.png', dpi=170, bbox_inches='tight', pad_inches=.05)
    plt.close(fig)


def native_materials():
    with np.load(RING/'strong_ring_arrays.npz', allow_pickle=False) as z:
        ids, active = z['ids'], z['chi_active']
    directory = STUDY/'ip_unified/numerics/feasibility/strict_endpoints'
    with np.load(directory/f'{LOW_ID}.npz', allow_pickle=False) as z:
        disk = z['material']
    ring = np.zeros_like(disk)
    mask = np.zeros_like(disk, dtype=bool)
    ring[ids], mask[ids] = active, True
    if not np.all(mask[np.flatnonzero(disk)]) or np.count_nonzero(mask) != 1153:
        raise AssertionError('materials do not share the declared permitted region')
    return disk.reshape(96, 96), ring.reshape(96, 96), mask.reshape(96, 96)


def material_panel(ax, material, permitted, title, x_offset, y_offset, first):
    xy = np.indices(permitted.shape)
    ii = xy[0, permitted]
    jj = xy[1, permitted]
    i0, i1 = ii.min()-4, ii.max()+5
    j0, j1 = jj.min()-4, jj.max()+5
    x = (np.arange(j0, j1+1)-48)/48-x_offset
    y = (np.arange(i0, i1+1)-48)/48-y_offset
    im = ax.pcolormesh(x, y, material[i0:i1, j0:j1], cmap='Blues', vmin=0, vmax=3.5,
                       shading='flat', rasterized=False)
    ax.contour((x[:-1]+x[1:])/2, (y[:-1]+y[1:])/2, permitted[i0:i1, j0:j1],
               levels=[.5], colors=[MUTED], linewidths=.62, linestyles='--')
    ax.set(xlim=(-.49, .49), ylim=(-.49, .49), aspect='equal',
           xticks=[-.4, 0, .4], yticks=[-.4, 0, .4], xlabel='Horizontal offset (m)')
    ax.set_title(title, fontsize=8.5, pad=5)
    if first:
        ax.set_ylabel('Vertical offset (m)')
    else:
        ax.tick_params(labelleft=False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color(LIGHT)
    return im


def native_plot():
    data = read(HERE/'source_index.json')['rows']
    disk, ring, permitted = native_materials()
    yy, xx = np.nonzero(permitted)
    center = ((xx.mean()+.5-48)/48, (yy.mean()+.5-48)/48)
    fig = plt.figure(figsize=(7.25, 4.90))
    fig.subplots_adjust(top=.87, bottom=.11, left=.09, right=.98, hspace=.53, wspace=.30)
    layout = fig.add_gridspec(2, 3, height_ratios=(1.1, .85))
    ax0, ax1, info = [fig.add_subplot(layout[0, k]) for k in range(3)]
    material_panel(ax0, disk, permitted, 'Low-peak disk\n416 nonzero cells', *center, first=True)
    im = material_panel(ax1, ring, permitted, 'Strong hollow ring\n768 nonzero cells', *center, first=False)
    info.set_axis_off()
    info.plot([.0, .16], [.91, .91], color=MUTED, lw=.85, linestyle='--',
              transform=info.transAxes, clip_on=False)
    info.text(.20, .91, '1,153 permitted cells', ha='left', va='center', transform=info.transAxes)
    cax = info.inset_axes([0., .67, .90, .09])
    color = fig.colorbar(im, cax=cax, orientation='horizontal', ticks=[0., .625, 2., 3.5])
    color.ax.set_xticklabels(['0', '0.625', '2', '3.5'])
    color.outline.set_visible(False)
    color.ax.set_title('True dielectric contrast', loc='left', fontsize=8, pad=4)
    handles = [
        Line2D([], [], color=BLUE, marker='o', lw=0, label='Low-peak disk'),
        Line2D([], [], color=ORANGE, marker='s', lw=0, label='Strong hollow ring'),
        Line2D([], [], color=MUTED, lw=1.35, label='Receiver deficit'),
        Line2D([], [], color=MUTED, lw=1.35, linestyle='--', label='Source-only deficit'),
    ]
    info.legend(handles=handles, loc='lower left', bbox_to_anchor=(-.04, -.01),
                frameon=False, fontsize=7.7, handlelength=1.8)
    for k, (key, title) in enumerate((('solid_half_1.5', 'Solid box, half-width 1.5'),
                                       ('solid_half_2.5', 'Solid box, half-width 2.5'),
                                       ('cap_4', 'Local cap 4'))):
        ax = fig.add_subplot(layout[1, k])
        for obj, color, marker in (('disk', BLUE, 'o'), ('hollow ring', ORANGE, 's')):
            observed = sorted((r for r in data if r['object'] == obj and r['box'] == key), key=lambda r:r['s'])
            x = np.log2([r['s'] for r in observed])
            actual = [r['receiver_directions'] - r['direct_receiver_inertia'][0]
                      if r['direct_receiver_inertia'] is not None else np.nan for r in observed]
            source = [r.get('deficit_source_only', np.nan) for r in observed]
            ax.plot(x, source, color=color, lw=1.2, alpha=.73, linestyle='--')
            # A thin white outline keeps coincident colored markers legible.
            ax.plot(x, actual, color=color, lw=1.35, marker=marker, ms=4.2,
                    markeredgewidth=.43, markeredgecolor='white')
        ax.set(title=title, xlim=(-.16, 5.16), ylim=(-.38, 13.7), yticks=[0, 4, 8, 12])
        ax.set_xticks(range(6), [1, 2, 4, 8, 16, 32])
        ax.grid(axis='y', lw=.55, color=GRID)
        ax.axhline(0, color=LIGHT, lw=.75)
        if k == 0:
            ax.set_ylabel('Missing positive directions')
        if k == 1:
            ax.set_xlabel('Number of old illuminations')
    fig.text(.09, .967, 'Native TM source index and conditional receiving directions',
             fontsize=10.5, fontweight='bold', color=DARK)
    save(fig, 'native-source-index')


def scaling_plot():
    with (V4/'routing/principal_28.csv').open() as f:
        large = list(csv.DictReader(f))
    compact = read(STUDY/'validation/max_struct/summary.json')['frozen_r10_profile']['cases']
    if len(large) != 28 or len(compact) != 7:
        raise AssertionError('historical scale cohorts changed')
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.25, 2.95), gridspec_kw={'width_ratios':[1.67, 1]})
    fig.subplots_adjust(left=.11, right=.98, bottom=.17, top=.76, wspace=.25)
    xs = np.arange(1, len(large)+1)
    for col, label, color, marker in [
        ('l1_lower_over_radius', 'First Born, 28/28', BLUE, 'o'),
        ('l2_cap05_distance_over_radius', 'Cap 0.50, 28/28', GREEN, '^'),
        ('l2_cap075_distance_over_radius', 'Cap 0.75, 0/28', ORANGE, 'v'),
    ]:
        ratio = [float(r[col]) for r in large]
        expected = int(label.split()[-1].split('/')[0])
        if np.count_nonzero(np.array(ratio)>1) != expected:
            raise AssertionError('archived large-scale counts differ')
        ax.scatter(xs, ratio, s=13, c=color, marker=marker, label=label, zorder=3)
    ax.legend(loc='upper right', frameon=False, fontsize=6.9, handletextpad=.35)
    ax.set(ylabel='Distance lower bound / noise radius', xlabel='Material',
           xticks=[1, 6, 11, 16, 21, 28], xlim=(.2, 28.8), ylim=(-.25, 9.5))
    ratios = np.array([r['joint_interval'][0]/r['radius'] for r in compact.values()])
    if np.count_nonzero(ratios>1) != 5:
        raise AssertionError('archived compact-scale count differs')
    indices = np.arange(1, 8)
    bx.bar(indices, ratios, width=.55, color=GREEN, zorder=3)
    record_of_witness = list(compact).index('r14_new_0018')+1
    remaining = [i for i, v in zip(indices, ratios, strict=True) if v == 0 and i != record_of_witness]
    bx.scatter(remaining, [0]*len(remaining), s=19, edgecolors=MUTED, facecolors='white', zorder=4)
    bx.scatter([record_of_witness], [0], color=ORANGE, marker='D', s=17, zorder=4)
    bx.legend(handles=[
        Line2D([], [], color=GREEN, marker='s', lw=0, label='Cap 5, 5/7'),
        Line2D([], [], color=ORANGE, marker='D', lw=0, label='Feasible peak 4.282'),
    ], loc='upper right', frameon=False, fontsize=6.9, handletextpad=.35)
    bx.set(xlabel='Material', xticks=indices, ylim=(-.28, 15))
    headings = ((ax, 'Large domain, 28 materials', r'$ka=11.85$ · $10\times20$ ports · 10% noise'),
                (bx, 'Compact annuli, 7 materials', r'$ka=4.27$ · $12\times48$ ports · 1% noise'))
    for panel, title, detail in headings:
        x = panel.get_position().x0
        fig.text(x, .866, title, color=DARK, fontsize=8.5)
        fig.text(x, .805, detail, color=MUTED, fontsize=7)
        panel.grid(axis='y', lw=.52, color=GRID)
        panel.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
        panel.axhline(1, color=MUTED, lw=.8, linestyle='--')
    fig.text(.11,.977,'Global cap exclusions across two electrical scales',
             color=DARK, fontsize=10.5, fontweight='bold')
    save(fig, 'historical-scaling')


if __name__ == '__main__':
    plt.rcParams.update({
        'font.family':'DejaVu Sans', 'font.size':7.5, 'axes.titlesize':8.5,
        'axes.labelsize':7.6, 'axes.spines.top':False, 'axes.spines.right':False,
        'axes.edgecolor':MUTED, 'axes.labelcolor':DARK, 'text.color':DARK,
        'xtick.color':MUTED, 'ytick.color':MUTED, 'pdf.fonttype':42, 'svg.fonttype':'none',
    })
    native_plot()
    scaling_plot()
