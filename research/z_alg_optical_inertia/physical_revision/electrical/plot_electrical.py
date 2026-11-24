"""Plot the angular comparison and chamber state counts from the frozen tables."""

import csv
from pathlib import Path

import matplotlib.font_manager as fonts
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator, PercentFormatter
import numpy as np

HERE = Path(__file__).resolve().parent
NAMES = {
    'en': ('Slotted block', 'Four lobes', 'Double chevron',
           'Crescent and plate', 'Staggered grooves', 'Graded pinwheel'),
    'zh': ('开槽块', '四叶交叠', '双折带', '月牙和板', '错列沟槽', '梯度风车'),
}
COLORS = ('#146994', '#d2741d', '#96448b', '#087b6d', '#765803', '#527ad2')
MARKERS = ('o', '^', 'D', 'v', 's', 'P')
INK, MID, GRID = '#1b2b36', '#4e5966', '#e5eaee'


def read_csv(name):
    with (HERE/name).open(newline='') as source:
        return list(csv.DictReader(source))


def load_angular():
    rows = read_csv('angular_interpolation.csv')
    values = {}
    for row in rows:
        key = tuple(int(row[entry]) for entry in ('shape', 'receivers', 'from_sources'))
        if key in values or int(row['to_sources']) != 32:
            raise ValueError('Angular table requires unique rows and a common 32-angle reference.')
        values[key] = float(row['clean_ports_rel_error'])
    planned = {(shape, p, s) for shape in range(1,7) for p in (20,32) for s in (10,16,24)}
    if set(values) != planned or any(v <= 0 for v in values.values()):
        raise ValueError('Angular conditions are incomplete.')
    return values


def load_chamber():
    g = sorted((r for r in read_csv('geometry.csv') if r['dataset'].startswith('Manitoba ')),
               key=lambda r: float(r['frequency_GHz']))
    count = {(r['dataset'], float(r['cap'])): r for r in read_csv('state_counts.csv')}
    caps = (.5, .75, 1.)
    sources = {int(r['s']) for r in g}
    if len(g) != 7 or len(sources) != 1:
        raise ValueError('Expected seven Manitoba frequencies with a common per-group source count.')
    source_count = sources.pop()
    for r in g:
        if any(int(count[r['dataset'], u]['source_columns']) != source_count for u in caps):
            raise ValueError('The geometry and state-count tables disagree on source columns.')
    f = np.array([float(r['frequency_GHz']) for r in g])
    size = np.array([float(r['two_ka']) for r in g])
    data = {u: [int(count[r['dataset'], u]['hidden']) for r in g] for u in caps}
    return f, size, source_count, data


def configure(lang):
    cjk = Path('/System/Library/Fonts/STHeiti Light.ttc')
    if lang == 'zh' and cjk.is_file():
        fonts.fontManager.addfont(str(cjk))
        cjk_name = fonts.FontProperties(fname=cjk).get_name()
    else:
        cjk_name = 'Noto Sans CJK SC'
    plt.rcParams.update({
        'font.family': ['DejaVu Sans'] if lang == 'en' else [cjk_name, 'DejaVu Sans'],
        'font.size': 9.2, 'axes.titlesize': 10.5, 'axes.labelsize': 9.2,
        'axes.spines.top': False, 'axes.spines.right': False,
        'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': '#839099',
        'xtick.color': MID, 'ytick.color': MID, 'grid.color': GRID,
        'grid.linewidth': .62, 'savefig.facecolor': 'white', 'pdf.fonttype': 42,
    })


def panel(ax, tag, title, *, y=None):
    ax.set_title(f'({tag})  {title}', loc='left', pad=7, y=y)
    ax.set_axisbelow(True)
    ax.grid(axis='y')
    ax.tick_params(axis='both', length=2.4, pad=3.2)


def angular(ax, lang, values):
    en = lang == 'en'
    prep = np.array([10, 16, 24])
    for i, (color, marker) in enumerate(zip(COLORS, MARKERS, strict=True), 1):
        y = [values[i,20,s] for s in prep]
        y32 = [values[i,32,s] for s in prep]
        ax.plot(prep, y, marker=marker, markersize=4.7, color=color, linewidth=1.14,
                label=f'{i:02d}  {NAMES[lang][i-1]}')
        ax.plot(prep, y32, marker=marker, markersize=6.6, linestyle='none',
                color=color, markerfacecolor='none', markeredgewidth=.7)
    ax.set_yscale('log')
    ax.set_xticks(prep)
    ax.set_xlim(8.25, 25.75)
    ax.set_ylim(0.000055,1.23)
    ax.set_yticks([.0001, .001, .01, .1, 1],
                  labels=['0.01%', '0.1%', '1%', '10%', '100%'])
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_xlabel('Incident angles used' if en else '用于插值的入射角数')
    ax.set_ylabel('Relative field error vs 32 angles' if en else '相对 32 个入射角的复场误差')
    panel(ax,'a','Full-wave angular interpolation' if en else '全波复场的角向插值', y=1.23)
    ax.text(.5, 1.129, 'Equal total incident fluence' if en else '总入射能量通量相同',
            transform=ax.transAxes, ha='center', va='bottom', color=MID, fontsize=8.5)
    return ax.get_legend_handles_labels()


def chamber(ax, lang, f, size, source_count, data):
    en = lang == 'en'
    for (cap, values), color, marker in zip(data.items(),
                                            ('#216c99','#087b6d','#b45c1b'),
                                            ('o','s','^'), strict=True):
        ax.plot(size, values, color=color, marker=marker, markersize=4.2,
                linewidth=1.45, label=f'$u={cap:.2f}$')
    ax.axhline(source_count, color=MID, linewidth=.85, linestyle=(0,(4,3)),
               label=f'{source_count} inputs per group' if en else f'每组 {source_count} 束')
    ratio = size[0]/f[0]
    if not np.allclose(size, ratio*f, rtol=1e-12, atol=1e-14):
        raise ValueError('Frequency and electrical size are not proportional.')
    top = ax.secondary_xaxis('top', functions=(lambda x: x/ratio, lambda x: x*ratio))
    top.set_xticks([3,4,5,6], labels=['3', '4', '5', '6'])
    top.tick_params(length=2.4, pad=1.5, colors=MID)
    top.set_xlabel('Frequency (GHz)' if en else '频率（GHz）', labelpad=3)
    ax.set_xlabel('Electrical size $2ka$' if en else '电尺寸 $2ka$')
    ax.set_ylabel('Number of negative directions' if en else '负方向数')
    ax.set_xlim(size[0]-.4, size[-1]+.45)
    ax.set_xticks(size[[0,2,4,6]], labels=[f'{x:.2f}' for x in size[[0,2,4,6]]])
    ax.set_ylim(-.9, 21.4)
    ax.set_yticks([0,4,8,12,16,20])
    panel(ax,'b','Manitoba receiver-dark states' if en else 'Manitoba 接收盲区的状态', y=1.23)
    ax.legend(loc='upper left', frameon=False, ncols=2, handletextpad=.42,
              columnspacing=.65, fontsize=7.5, borderaxespad=.35)


def save(fig, stem):
    fig.savefig(HERE/(stem+'.pdf'), bbox_inches='tight')
    fig.savefig(HERE/(stem+'.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)


def main_figure(lang, values, chamber_data):
    en = lang == 'en'
    fig, axes = plt.subplots(1,2,figsize=(7.25,3.5))
    fig.subplots_adjust(left=.105,right=.968,top=.778,bottom=.34,wspace=.35)
    handles, labels = angular(axes[0], lang, values)
    chamber(axes[1], lang, *chamber_data)
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.49,.075),
               frameon=False,ncols=3,fontsize=7.7,handletextpad=.45,columnspacing=1.35)
    fig.text(.49,.044,
             'Filled markers, 20 receivers    Open outlines, 32 receivers' if en else
             '实心标记为 20 个接收点    空心轮廓为 32 个接收点',
             color=MID,ha='center',fontsize=7.5)
    save(fig, 'electrical_scaling_'+lang)


def chamber_cap(lang):
    en = lang == 'en'
    rows = sorted((r for r in read_csv('chamber_by_frequency.csv') if r['object'] == 'nylon_large'),
                  key=lambda r: float(r['frequency_GHz']))
    if len(rows) != 7:
        raise ValueError('Expected the seven frequency rows for the large nylon object.')
    f = [float(r['frequency_GHz']) for r in rows]
    data = {
        'lower': [float(r['u1_relative_distance']) for r in rows],
        'first': [float(r['first_calibration']) for r in rows],
        'second': [float(r['second_calibration']) for r in rows]
    }
    for row, d, a, b in zip(rows, *data.values(), strict=True):
        if not np.isclose(float(row['u1_signed_margin']), d-max(a,b), rtol=1e-10, atol=1e-13):
            raise ValueError('The margin does not match the two calibration scales.')
    fig, ax = plt.subplots(figsize=(4.25,2.55),layout='constrained')
    ax.plot(f,data['lower'], '-o',color='#216c99',markersize=4,linewidth=1.55,
            label='Lower bound at $u=1$' if en else '$u=1$ 时的距离下界')
    ax.plot(f,data['first'], '--s',color='#8f4386',markersize=3.6,linewidth=1.1,
            label='Two metal references' if en else '两件金属标定参照')
    ax.plot(f,data['second'], ':^',color='#b7671d',markersize=4,linewidth=1.3,
            label='Empty scene + held-out metal' if en else '空场与留出金属参考')
    ax.set_title('Manitoba large nylon object' if en else 'Manitoba 大尼龙目标',loc='left')
    ax.set_xlabel('Frequency (GHz)' if en else '频率（GHz）')
    ax.set_ylabel('Relative field distance' if en else '相对复场距离')
    ax.set_xticks(f, labels=[f'{x:g}' for x in f])
    ax.yaxis.set_major_formatter(PercentFormatter(1,decimals=0))
    ax.set_ylim(-.025,.95)
    ax.set_yticks([0,.2,.4,.6,.8])
    ax.set_axisbelow(True)
    ax.grid(axis='y')
    ax.legend(frameon=False,loc='upper left',fontsize=7.2,handletextpad=.4,columnspacing=.55)
    save(fig,'manitoba_large_nylon_cap_'+lang)


def main():
    values = load_angular()
    chamber_data = load_chamber()
    for lang in ('en','zh'):
        configure(lang)
        main_figure(lang,values,chamber_data)
        chamber_cap(lang)


if __name__ == '__main__':
    main()
