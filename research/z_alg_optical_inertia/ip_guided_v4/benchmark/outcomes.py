#!/usr/bin/env python3
"""Render the branch, candidate-pool and physical-decision figures."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch
from matplotlib.ticker import NullFormatter, NullLocator
import numpy as np

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
OUT = STUDY / 'figures'
BOXES = STUDY / 'boxes/results'
FAMILIES = STUDY / 'families'
CONFIRMATION = STUDY / 'threshold/confirmation'
FLOW = STUDY / 'pool/flow_results'

INK, SOFT, GRID = '#202b33', '#53636c', '#dce4e8'
BLUE, GREEN, RED, VIOLET = '#245e7d', '#18765b', '#ad4f37', '#725c97'
OPEN, CONTROL = '#e5ebed', '#f0f4f6'


def read(path):
    return json.loads(path.read_text())


def figure_style():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.1,
                         'axes.titlesize': 9.5, 'axes.titlepad': 7,
                         'axes.labelsize': 8.4, 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.edgecolor': SOFT,
                         'text.color': INK, 'axes.labelcolor': INK,
                         'xtick.color': SOFT, 'ytick.color': SOFT,
                         'pdf.fonttype': 42, 'pdf.use14corefonts': False})


def save(fig, stem):
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / f'{stem}.pdf')
    fig.savefig(OUT / f'{stem}.png', dpi=270)
    plt.close(fig)


def grid_y(ax):
    ax.grid(axis='y', color=GRID, linewidth=.55)
    ax.set_axisbelow(True)


def disk_panel(geom, ax):
    disk = read(BOXES / 'disk_radius.json')
    obs = disk['observation']
    width = disk['width_bracket']
    verified = STUDY.parent / 'ip_unified/numerics/feasibility/strict_endpoints'
    with np.load(verified / f'{obs}.npz', allow_pickle=False) as a:
        material = a['material'].reshape(96, 96)
    prefix = obs + '__rank_tail_disk__'
    with np.load(BOXES / 'followup_disk_witnesses.npz', allow_pickle=False) as a:
        active = a[prefix + 'active']
    assert len(active) == width['support_cells']
    mask = np.zeros(96 * 96)
    mask[active] = 1.
    scale = np.linspace(-1 + 1 / 96, 1 - 1 / 96, 96)
    im = geom.imshow(material, origin='lower', extent=(-1, 1, -1, 1), vmin=0.,
                     vmax=width['witnessed_peak'], cmap='Blues', interpolation='nearest')
    geom.contour(scale, scale, mask.reshape(96, 96), levels=[.5],
                 colors=[RED], linewidths=1.1, linestyles='--')
    geom.set(title='A  Verified fitting disk', xlabel='Horizontal position (m)',
             ylabel='Vertical position (m)', xticks=[-1, 0, 1], yticks=[-1, 0, 1])
    geom.spines[['top', 'right']].set_visible(True)
    bar = geom.figure.colorbar(im, ax=geom, shrink=.77, fraction=.07, pad=.04)
    bar.set_ticks([0., width['witnessed_peak']], labels=['0', f"{width['witnessed_peak']:.3f}"])
    bar.set_label('Material contrast', fontsize=7.5)
    geom.text(-.90, .80, 'Dashed outline\nShared support', fontsize=6.9, color=RED,
              ha='left', va='top', bbox=dict(boxstyle='square,pad=.24', fc='white', ec='none'))

    lower_w = width['exclusion_width_safe_from_saved_endpoints'][0]
    fit_w = width['width_witness_entry'][1]
    contraction = disk['norm_threshold']['sufficient_2norm_half_width_bracket'][0]
    lower_d = width['exclusion_residual_over_radius_lower']
    upper_d = width['witness_residual_over_radius_upper']
    ax.axvspan(0., contraction, facecolor='#e9e3f0', edgecolor='none')
    ax.axvline(contraction, ls=':', lw=1.0, color=VIOLET)
    ax.axhline(1., color=SOFT, ls='--', lw=.85)
    ax.plot([lower_w, fit_w], [.29, .29], color=SOFT, lw=1., zorder=2)
    ax.plot([lower_w, fit_w], [.29, .29], marker='|', linestyle='none',
            color=SOFT, markersize=7, zorder=2)
    ax.text((lower_w+fit_w)/2, .16, 'First-fit width bracket', ha='center', va='center',
            fontsize=7.25, color=SOFT)
    ax.scatter(lower_w, lower_d, marker='^', s=40, c=RED, zorder=5)
    ax.scatter(fit_w, upper_d, marker='v', s=37, c=GREEN, zorder=5)
    ax.text(.14, 1.10, 'Neumann norm cutoff', color=VIOLET, fontsize=7.2, rotation=90)
    ax.text(lower_w-.07, lower_d+.13, f'Entire box excluded\nLower bound {lower_d:.3f}',
            color=RED, fontsize=7.6, ha='right')
    ax.text(fit_w+.02, upper_d+.14, f'A verified fit\nUpper bound {upper_d:.5f}',
            color=GREEN, fontsize=7.45, ha='right')
    ax.text(.66, .87, 'Noise threshold', fontsize=7.0, color=SOFT)
    rho = disk['spectral_check']['ideal_exclusion_width_operator_rho_lower']
    ax.text(.61, .51, rf'At the excluded width, $\rho(T)\geq {rho:.3f}$',
            fontsize=7.3, color=SOFT)
    ax.set(title='B  How wide can the same-support box be?',
           xlim=(0., 4.10), ylim=(0., 2.25), xlabel='Half-width of box centred at peak 4',
           ylabel='Residual bound / noise radius',
           xticks=[0, contraction, lower_w, fit_w],
           xticklabels=['0', f'{contraction:.3f}', f'{lower_w:.2f}', f'{fit_w:.3f}'],
           yticks=[0, .5, 1, 1.5, 2])
    grid_y(ax)


def nested_panels(axes):
    original = read(FAMILIES / 'results.json')
    nested = read(FAMILIES / 'nested.json')
    ids = original['source_cases'] + [row['id'] for row in original['controls']]
    # The same parent multiplier is used at both sites and transported to each child.
    mu, bg = 2**-10, .125
    candidates = [(None,'Patch range [0, 6]', SOFT, 'x'),
                  (2., 'Patch range [2, 6]', BLUE, 'o'),
                  (1., 'Patch range [3, 5]', GREEN, 's'),
                  (.5, 'Patch range [3.5, 4.5]', RED, '^')]
    for i, ax in enumerate(axes, start=1):
        broad = next(row for row in original['runs'] if (row['site'], row['b'], row['mu'],
                     row['center'], row['halfwidth']) == (i, bg, mu, 3, 3))
        conditions = {None: {row['id']:row['ratio_lower'] for row in broad['readings']}}
        for a in nested['runs']:
            if (a['site'],a['b'],a['mu']) == (i,bg,mu):
                conditions[a['width']] = a['results']
        assert set(conditions) == {k for k,_,_,_ in candidates}
        x = np.arange(8)
        ax.axvspan(5.53, 7.49, color=CONTROL, lw=0, zorder=0)
        ax.axhline(1., ls='--', color=SOFT, lw=.8)
        for j,(w,_,color,marker) in enumerate(candidates):
            values = [conditions[w][name] for name in ids]
            ax.scatter(x+(j-1.5)*.17, values, marker=marker, s=21, linewidth=.9,
                       c=color, zorder=3)
        ax.text(6.51, 6.28, 'Site controls', color=SOFT, fontsize=6.8, ha='center')
        ax.set(title=f'{"C" if i==1 else "D"}  Fixed site {i}, background contrast 0–⅛',
               xticks=x, xticklabels=['1','2','3','4','5','6','A','B'],
               xlim=(-.49, 7.49), ylim=(0., 6.7), xlabel='Material',
               ylabel='Lower bound / noise radius' if i==1 else '')
        if i==2:
            ax.tick_params(axis='y', labelleft=False)
        grid_y(ax)
    return [Line2D([],[], ls='none', marker=m, markerfacecolor=c, markeredgecolor=c,
                   color=c, markersize=4.5, label=lab) for _,lab,c,m in candidates]


def branches():
    f = plt.figure(figsize=(7.3, 6.15), layout='constrained')
    outer = f.add_gridspec(3, 1, height_ratios=(1.07,1.,.10))
    top = outer[0].subgridspec(1, 2, width_ratios=(1., 1.60))
    disk_panel(f.add_subplot(top[0]), f.add_subplot(top[1]))
    lower = outer[1].subgridspec(1, 2)
    left = f.add_subplot(lower[0])
    legend = nested_panels((left, f.add_subplot(lower[1], sharey=left)))
    footer = f.add_subplot(outer[2])
    footer.axis('off')
    footer.legend(handles=legend, loc='center', ncol=4, frameon=False, fontsize=7.1,
                  columnspacing=1.2, handletextpad=.4)
    save(f,'branches')


def flow_counts():
    root = read(FLOW / 'native_l4.json')['cases']
    point = read(FLOW / 'strict_point_witnesses.json')['cases']
    counts = []
    for key in root:
        boxes = root[key]['boxes']
        excluded = {b for b, facts in boxes.items() if facts['status']=='STRICTLY_EXCLUDED'}
        witness = {b for b, facts in point[key]['boxes'].items() if facts['strictly_feasible']}
        assert not excluded & witness and excluded | witness <= set(boxes)
        counts.append((len(excluded),len(witness),len(set(boxes)-excluded-witness)))
    assert tuple(np.sum(counts, axis=0)) == (2, 15, 19)
    return np.asarray(counts)


def pools():
    frozen = read(CONFIRMATION / 'public/cases.json')
    online = read(CONFIRMATION / 'online.json')
    offline = read(CONFIRMATION / 'offline.json')
    names = ['Gaussian','Ring','Two lobes','Step','Crescent','Ridges','Rectangle','Shoulders']
    ids = [v['id'] for v in frozen]
    n = len(ids)
    np.testing.assert_equal(n, len(names))
    x = np.arange(n)
    rejected = np.array([len(offline[k]['excluded']) for k in ids])
    covered = np.array([len(offline[k]['covered']) for k in ids])
    box_counts = np.array([len(online['cases'][k]['methods']['fast_boxes']['decisions']) for k in ids])
    for key in ids:
        assert not set(offline[key]['excluded']) & set(offline[key]['covered'])

    fig = plt.figure(figsize=(7.3, 5.8), layout='constrained')
    axes = fig.add_gridspec(2, 2, height_ratios=(.89,1.14), width_ratios=(1.32,1.))
    top = fig.add_subplot(axes[0,:])
    top.bar(x,rejected,color=RED,label='Entire box excluded',width=.63,zorder=2)
    top.bar(x,box_counts-rejected,bottom=rejected,color=OPEN,label='Other retained boxes',width=.63,zorder=2)
    top.bar(x,covered,bottom=box_counts-covered,color=GREEN,label='Retained box containing the generating material',
            width=.63,zorder=3)
    for i,num in enumerate(rejected):
        top.text(i,num-.52,str(num),color='white',fontweight='bold',ha='center',va='top',fontsize=7.7)
    top.set(title='A  Eight prospective materials, 16 candidate boxes each',
            ylim=(0,18.55),yticks=[0,4,8,12,16],
            ylabel='Candidate boxes',xticks=x,xticklabels=names,xlim=(-.60, n-.40))
    grid_y(top)
    top.legend(frameon=False,loc='upper center',fontsize=7.1,ncol=3,
               columnspacing=1.2,handlelength=1.1,handletextpad=.4)

    timing = fig.add_subplot(axes[1,0])
    tools = [('all_cluster','Refine\nall',BLUE),('top_early','Point ranked,\nstop at fit',GREEN),
             ('fast_boxes','Screen and\nrefine',VIOLET),('fast_boxes','Screen and\nrefine + replay',RED)]
    for i,(engine,label,c) in enumerate(tools):
        vals = np.array([online['cases'][k]['methods'][engine]['wall_seconds'] for k in ids])
        if engine=='fast_boxes':
            vals += online['setup']['wall_seconds']/n
        if i==3:
            vals += [offline[k]['independent_ball_recheck']['wall_seconds'] for k in ids]
        offsets = np.linspace(-.13,.13,n)
        timing.scatter(i+offsets,vals,c=c,s=16,alpha=.84,zorder=3)
        mdn = np.median(vals)
        timing.plot([i-.24,i+.24],[mdn,mdn],color=INK,lw=1.4,zorder=4)
        timing.text(i,mdn*1.20,f'{mdn:.3f}',fontsize=6.6,ha='center',va='bottom',
                    bbox=dict(fc='white',ec='none',pad=.55),zorder=5)
    timing.set(title='B  Paid end-to-end wall time',
               ylabel='Wall time (s, log scale)',xticks=range(4),
               xticklabels=[label for _,label,_ in tools],xlim=(-.40,3.40),yscale='log',
               ylim=(.024,.66))
    timing.tick_params(axis='x', labelsize=7.2)
    timing.set_yticks([.03,.06,.12,.24,.48],labels=['.03','.06','.12','.24','.48'])
    timing.yaxis.set_minor_locator(NullLocator())
    timing.yaxis.set_minor_formatter(NullFormatter())
    grid_y(timing)

    flow = fig.add_subplot(axes[1,1])
    values = flow_counts()
    sx = np.arange(len(values))
    bottom = np.zeros(len(values))
    cats = [('Excluded',RED),('Contains a verified fit',GREEN),('Unresolved',OPEN)]
    for i,(label,c) in enumerate(cats):
        col=values[:,i]
        flow.bar(sx,col,bottom=bottom,color=c,width=.58,zorder=4,label=label)
        for j,total in enumerate(col):
            if total:
                flow.text(j,bottom[j]+total/2,str(total),ha='center',va='center',
                          color=INK if i==2 else 'white',fontsize=8.3,fontweight='bold',zorder=5)
        bottom += col
    flow.set(title='C  Signed-contrast flow proposals',
             ylim=(0,15.9),yticks=[0,3,6,9,12],xticks=sx,xticklabels=['1','2','3'],
             xlabel='Material',ylabel='Candidate boxes')
    grid_y(flow)
    handles=[Patch(facecolor=c,label=f'{label} ({values[:,i].sum()})')
             for i,(label,c) in enumerate(cats)]
    flow.legend(handles=handles,frameon=False,loc='upper center',fontsize=7.0,
                ncol=1,handlelength=1,labelspacing=.30)
    save(fig,'pools')


def workflow():
    fig,ax = plt.subplots(figsize=(7.3,4.15))
    ax.set(xlim=(0,1),ylim=(0,1))
    ax.axis('off')
    fig.subplots_adjust(left=.012,right=.988,bottom=.025,top=.975)

    def box(x,y,w,h,color,heading,body=None):
        style = 'round,pad=.010,rounding_size=.013'
        bg = FancyBboxPatch((x,y),w,h,boxstyle=style,edgecolor=color,
                            linewidth=.9,facecolor='white',zorder=2)
        ax.add_patch(bg)
        pad=.019 if body else w/2
        ax.text(x+pad,y+h-(.035 if body else h*.5),heading,ha='left' if body else 'center',
                va='top' if body else 'center',fontsize=10.1 if body else 10.5,
                fontweight='bold',color=color,zorder=3)
        if body:
            ax.text(x+pad,y+h-.10,body,ha='left',va='top',fontsize=8.0,linespacing=1.36,
                    color=INK,zorder=3)

    def arrow(p,q):
        ax.add_patch(FancyArrowPatch(p,q,arrowstyle='-|>',
                                    mutation_scale=9,lw=.80,color='#788994',zorder=1))

    group = FancyBboxPatch((.012,.230),.976,.513,boxstyle='round,pad=.004,rounding_size=.013',
                           edgecolor=GRID,linewidth=.8,facecolor='#f4f7f8',zorder=0)
    ax.add_patch(group)
    box(.238,.856,.524,.083,BLUE,'Measured complex fields and error radius')
    arrow((.5,.843),(.5,.754))
    ax.text(.5,.714,'Apply the tests that match the physical assumptions',ha='center',va='center',
            fontsize=8.6,color=SOFT)
    box(.042,.475,.433,.186,BLUE,'Test the first Born model',
        'Can a nonnegative first-Born image fit?\nReject it if the bound exceeds the error radius.')
    box(.525,.475,.433,.186,GREEN,'Test a peak or loss constraint',
        'Can any full-wave material in this class fit?\nApply the same threshold to the whole class.')
    box(.042,.260,.433,.177,VIOLET,'Check a current representation',
        'Can the representation carry the required current?\nKeep the necessary component for the tested class.')
    box(.525,.260,.433,.177,RED,'Test a proposed material interval',
        'Can any material in the interval fit?\nDrop a branch only when the whole interval fails.')
    arrow((.5,.227),(.5,.158))
    box(.147,.061,.706,.085,SOFT,'Refine survivors and check the full-wave predictions')
    save(fig,'workflow')


def material_alternatives():
    observations = read(FLOW / 'online.json')['cases']
    branches = read(FLOW / 'continuation.json')['current_pass']
    materials, ratios = [], []
    with np.load(FLOW / 'offline_geometry_scores.npz') as scores, np.load(FLOW / 'continuation_arrays.npz') as fits:
        for name, row in observations.items():
            with np.load(FLOW / f'{name}.npz') as points:
                best_point = int(np.argmin(row['point_residuals']))
                tested = branches[name]['refined_box_results']
                best_box = min(tested, key=lambda key: tested[key]['residual'])
                chosen = fits[f'{name}__box{int(best_box):02d}__material']
                materials.append([scores[f'{name}__truth'], points['samples'][best_point], chosen])
                ratios.append([None, row['point_residuals'][best_point]/row['radius'],
                               tested[best_box]['residual']/row['radius']])
    limit = np.ceil(4 * max(np.max(np.abs(x)) for group in materials for x in group))/4
    fig, axes = plt.subplots(3, 3, figsize=(6.45, 5.05), layout='constrained')
    labels = ['Generating material', 'Best sampled point', 'Refined physical fit']
    for column, (group, residuals) in enumerate(zip(materials, ratios)):
        for row, (field, ratio) in enumerate(zip(group, residuals)):
            ax = axes[row, column]
            image = ax.imshow(field.reshape(32, 32), vmin=-limit, vmax=limit,
                              cmap='RdBu_r', interpolation='nearest')
            ax.set(xticks=[], yticks=[])
            if column == 0:
                ax.set_ylabel(labels[row])
            if row == 0:
                ax.set_title(f'Material {column+1}')
            if ratio is not None:
                ax.set_xlabel(f'Residual / radius   {ratio:.2f}')
    fig.colorbar(image, ax=axes, fraction=.042, pad=.025, shrink=.82, label='Dielectric contrast')
    save(fig, 'material-alternatives')


def main():
    figure_style()
    branches()
    pools()
    workflow()
    material_alternatives()


if __name__ == '__main__':
    main()
