"""Draw bilingual program and material figures from sealed result arrays."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.font_manager as fonts
from matplotlib.ticker import PercentFormatter
import numpy as np

from gate import HERE

COMPACT = HERE/'compact_synthetic'
OUT = HERE/'figures'
BLUE, GREEN, GREY, ORANGE, PURPLE = '#216b9d','#087d77','#687380','#d17b31','#76358c'


def read(p):
    return json.loads(p.read_text())


def finish(fig, stem):
    fig.savefig(OUT/(stem+'.png'),dpi=220,bbox_inches='tight')
    fig.savefig(OUT/(stem+'.pdf'),bbox_inches='tight')
    plt.close(fig)


def originals(en, suffix):
    g = read(HERE/'original_gate.json')['decisions']
    labels = [f'{i:02d}' for i in range(1,7)]
    by = {(r['id'],r['cap']):r for r in g}
    fig,axes = plt.subplots(1, 2,figsize=(9.8,3.45),layout='constrained',sharey=True)
    times = [('pair1pct','1%',PURPLE,'o'),('pair2pct','2%',BLUE,'s'),('n10','10%',GREEN,'D')]
    for ax, cap in zip(axes,[.5,.75],strict=True):
        for shift,(key,level,color,mark) in zip([-.15,0,.15],times,strict=True):
            y = [by[(f'optfresh_{short}__pw360_10x20__{key}',cap)]['strict_ratio'] for short in labels]
            ax.scatter(np.arange(6)+shift,y,marker=mark,color=color,label=level,s=32,zorder=3)
        ax.axhline(1,color=GREY,linestyle='--',linewidth=1)
        ax.axhline(2.5,color=GREY,linestyle=':',linewidth=1)
        ax.set_yscale('symlog',linthresh=1)
        ax.set_yticks([0,1,2.5,10,25],labels=['0','1','2.5','10','25'])
        ax.set_ylim(-.10,33)
        ax.set_xticks(np.arange(6),labels)
        ax.set_title((f'Tested upper limit {cap:g}' if en else f'检验的对比度上限 {cap:g}'),loc='left',weight='bold')
        ax.set_xlabel('Object' if en else '物体编号')
        ax.grid(axis='y',color='#dee6ea',linewidth=.65)
    axes[0].set_ylabel('Strict lower bound / noise radius' if en else '严格距离下界 / 噪声半径')
    axes[1].legend(frameon=False,title='Paired 1%, 2%; independent 10%' if en else '配对 1%、2%，独立 10%',loc='upper right')
    finish(fig,'original_class_exclusion_'+suffix)


def compact(en, suffix):
    gen = read(COMPACT/'offline_truth/generation.json')['records']
    criteria = read(COMPACT/'gate.json')
    score = read(COMPACT/'offline_scores.json')
    more = read(COMPACT/'exposed_followup/offline_scores.json')
    ids = [case['id'] for case in gen]
    noise = np.array([case['relative_noise'] for case in gen])
    low = [score[case]['fixed_low'] for case in ids]
    high = [score[case]['high_from_start'] for case in ids]
    cont = [more[case]['class_guided_continuation'] for case in ids]
    exact = [more[case]['exact_discrepancy'] for case in ids]
    fig,axes = plt.subplots(1,3,figsize=(12.5,3.55),layout='constrained')
    titles = ['Fit to the measured field','Prediction at unseen angles','Material error'] if en else [
              '观测复场的拟合', '未见角度的复场预测', '材料图误差']
    for ax,title in zip(axes,titles,strict=True):
        ax.set_title(title,loc='left',weight='bold')
        ax.set_xlabel('Relative noise' if en else '相对噪声')
        ax.set_xticks(noise)
        ax.xaxis.set_major_formatter(PercentFormatter(1.,decimals=1))
        ax.grid(axis='y',color='#dee6ea',linewidth=.65)
    ax = axes[0]
    ax.axhspan(1.,2.5,color='#e8eef3',zorder=0)
    ax.axhline(1,color=GREY,linestyle='--',linewidth=1)
    ax.set_ylim(.84,2.65)
    ax.set_ylabel('Distance / noise radius' if en else '距离 / 噪声半径')
    ax.plot(noise,[r['strict_ratio'] for r in criteria],linestyle=':',marker='o',
            color=PURPLE,label='Global lower bound' if en else '整类严格下界')
    for name,rows,res,mat,mark,color,lw in [
        ('Fixed 0.5 and 2.5r policy' if en else '固定 0.5 与 2.5r 残差',
         low,'train_residual_ratio_upper','relative_l2','s',GREY,1.9),
        ('Open to 2 immediately' if en else '立即放宽至 2',
         high,'train_residual_ratio_upper','relative_l2','^',BLUE,1.8),
        ('Open to 2 at exact r' if en else '按精确半径放宽至 2',
         exact,'residual_over_radius_upper','full_relative_l2','v',ORANGE,1.1),
        ('Increase cap in steps' if en else '逐级放宽',
         cont,'residual_over_radius_upper','full_relative_l2','D',GREEN,1.8),
    ]:
        axes[0].plot(noise,[r[res] for r in rows],marker=mark,color=color,lw=lw,label=name,ms=4)
        axes[1].plot(noise,[r['withheld_relative_l2'] for r in rows],marker=mark,color=color,lw=lw,ms=4)
        axes[2].plot(noise,[r[mat] for r in rows],marker=mark,color=color,lw=lw,ms=4)
    axes[1].set_ylabel('Relative complex-field error' if en else '复场相对误差')
    axes[1].yaxis.set_major_formatter(PercentFormatter(1.,decimals=0))
    axes[1].set_ylim(0.,.31)
    axes[2].set_ylabel('Material relative L2' if en else '材料相对 L2 误差')
    axes[2].set_ylim(.23,1.02)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncols=3,frameon=False)
    finish(fig,'compact_fit_vs_material_'+suffix)


def maps(en,suffix):
    key = 'compact_tm_g060'
    with np.load(COMPACT/'offline_truth/reference.npz',allow_pickle=False) as ref:
        truth = ref['material']
    with np.load(COMPACT/'fits'/(key+'.npz'),allow_pickle=False) as initial:
        low, high = initial['fixed_low_material'], initial['high_from_start_material']
    with np.load(COMPACT/'exposed_followup'/(key+'.npz'),allow_pickle=False) as replay:
        exact, path = replay['exact_discrepancy_material'], replay['class_guided_continuation_material']
    maps = (truth,low,high,exact,path)
    names = ('Generating material','Fixed at 0.5','Immediate release','Exact discrepancy','Cap continuation') if en else (
             '生成材料','固定在 0.5','立即放宽上限','精确半径判定','逐级放宽上限')
    fig,axes = plt.subplots(1,5,figsize=(11.3,2.75),layout='constrained')
    fig.suptitle('Compact dielectric patch at 16.5% noise' if en else '电小介质块，噪声 16.5%',x=.015,ha='left',weight='bold')
    for ax,name,data in zip(axes,names,maps,strict=True):
        im = ax.imshow(data.reshape(12,12),origin='lower',interpolation='nearest',cmap='viridis',vmin=0.,vmax=2.)
        ax.set_title(name,loc='left')
        ax.text(.03,.045,(f'Peak {data.max():.3g}' if en else f'峰值 {data.max():.3g}'),
                transform=ax.transAxes,color='white',fontsize=9,
                bbox=dict(facecolor='black',alpha=.57,edgecolor='none',pad=2.1))
        ax.set_xticks([])
        ax.set_yticks([])
    fig.colorbar(im,ax=axes,aspect=22,label='Dielectric contrast' if en else '介电对比度')
    finish(fig,'compact_material_maps_'+suffix)


def main():
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({'savefig.facecolor':'white','axes.spines.top':False,'axes.spines.right':False,
                         'font.size':10,'axes.labelcolor':'#283340','text.color':'#182633'})
    zh = Path('/System/Library/Fonts/STHeiti Light.ttc')
    if zh.exists():
        fonts.fontManager.addfont(str(zh))
    zh_font = fonts.FontProperties(fname=zh).get_name() if zh.exists() else 'Noto Sans CJK SC'
    for en,suffix in ((True,'en'),(False,'zh')):
        plt.rcParams['font.family']=['DejaVu Sans'] if en else [zh_font,'DejaVu Sans']
        originals(en,suffix)
        compact(en,suffix)
        maps(en,suffix)


if __name__ == '__main__':
    main()
