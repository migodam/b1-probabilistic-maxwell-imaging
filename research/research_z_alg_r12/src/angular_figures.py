"""Render the Chinese angle-power and equal-radiated-power comparison."""
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.font_manager import FontProperties, fontManager
from matplotlib.lines import Line2D
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def chinese_font():
    selected = os.environ.get('ZALG_CJK_FONT')
    candidates = [Path(selected)] if selected else [Path('/System/Library/Fonts/Hiragino Sans GB.ttc'),
                                                   Path('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')]
    font = next((p for p in candidates if p.is_file()), None)
    if font is None:
        raise FileNotFoundError('Set ZALG_CJK_FONT to a Chinese font file')
    fontManager.addfont(font)
    return FontProperties(fname=font).get_name()


def main():
    plt.rcParams.update({'font.family':chinese_font(), 'axes.unicode_minus':False,
                         'font.size':11, 'axes.spines.top':False, 'axes.spines.right':False,
                         'pdf.fonttype':42, 'mathtext.fontset':'dejavusans'})
    data = json.loads((ROOT/'results/angular_gap.json').read_text())['cases']
    cases = list(data)
    with np.load(ROOT/'results/angular_gap.npz', allow_pickle=False) as z:
        orders = z['angular_orders']
        c = np.stack([z[case+'_modal_coefficients'] for case in cases])
    cap = 28
    power = ((c*c.conj()).real).sum(axis=2)
    fraction = np.stack([power[:,:,np.abs(orders)==n].sum(axis=2) for n in range(cap+1)],axis=2)
    fraction /= power.sum(axis=2)[:,:,None]
    fig = plt.figure(figsize=(12.0,8.1),facecolor='white',layout='constrained')
    grid = fig.add_gridspec(2, 3, height_ratios=[1.12,1.],width_ratios=[1,1,.039])
    tops = []
    for i, title in enumerate(('a  外圈半径 1.60','b  近圈半径 0.72')):
        ax = fig.add_subplot(grid[0,i])
        im = ax.imshow(np.log10(np.maximum(fraction[:,i],1e-8)), aspect='auto', cmap='magma',
                       norm=Normalize(vmin=-5,vmax=0), interpolation='nearest')
        ax.set_title(title,loc='left',fontweight='bold',fontsize=13,pad=7)
        ax.axvline(5.5,color='white',lw=1.2,ls=(0,(4,3)),alpha=.9)
        ax.set_xticks([0,3,6,10,14,18,22,28])
        ax.set_yticks(np.arange(len(cases)), [key.split('_')[1] for key in cases] if i==0 else [])
        ax.tick_params(axis='both',length=2)
        ax.set_xlabel(r'角谐波阶数 $|n|$')
        if i==0:
            ax.set_ylabel('历史材料对')
        tops.append(ax)
    cbax = fig.add_subplot(grid[0,2])
    cbar = fig.colorbar(im,cax=cbax,ticks=np.arange(-5,1))
    cbar.ax.set_yticklabels([r'$10^{%d}$'%i for i in range(-5,1)])
    cbar.set_label('占各材料对总角谱功率的比例')
    ax = fig.add_subplot(grid[1,:2])
    ax.set_title('c  同一发射功率下，一次 48 接收点曝光',loc='left',fontweight='bold',fontsize=13,pad=7)
    left, right = '#4463a1','#d56728'
    labels = []
    for index,(size,method,label) in enumerate((size,method,label) for size in (16,24) for method,label in
                                                [('single','单源'),('coherent_l2_2','馈电上限 2'),('coherent_l2_4','馈电上限 4')]):
        values = np.array([[row['arrays'][f'{R}_{size}']['designs'][method]['actual_gap_over_noise']
                            for R in ('1.60','0.72')] for row in data.values()])
        x = index + np.array([-.17,.17])
        for pair in values:
            ax.plot(x,pair,color='#898c94',lw=.8,alpha=.3,zorder=2)
        for i,color in enumerate((left,right)):
            ax.scatter(np.full(len(values),x[i]), values[:,i], color=color,s=19,alpha=.64,zorder=3,lw=0)
            m = np.median(values[:,i])
            ax.plot([x[i]-.10,x[i]+.10],[m,m],color=color,lw=3,solid_capstyle='round',zorder=4)
        labels.append(label)
    ax.axvline(2.5,color='#d7d7df',lw=1)
    ax.axhline(2.,color='#686b72',ls=(0,(4,3)),lw=.95,zorder=1)
    ax.set_yscale('log')
    ax.set_ylim(.087,19.5)
    ax.set_yticks([.1,.2,.5,1,2,5,10],labels=['0.1','0.2','0.5','1','2','5','10'])
    ax.set_ylabel('候选分离量 ÷ 单次接收噪声')
    ax.set_xticks(range(6),labels)
    ax.tick_params(axis='x',length=0,pad=6)
    ax.grid(axis='y',alpha=.14,which='major')
    top = ax.secondary_xaxis('top')
    top.set_xticks([1,4],['16 发射点','24 发射点'])
    top.tick_params(axis='x',length=0,pad=4,labelsize=10)
    top.spines['top'].set_visible(False)
    ax.legend(handles=[Line2D([0],[0],marker='o',ls='',color=left,label='外圈'),
                       Line2D([0],[0],marker='o',ls='',color=right,label='近圈'),
                       Line2D([0],[0],color='#686b72',ls=(0,(4,3)),label='2 倍噪声参考')],
              loc='upper left',ncol=3,frameon=False,fontsize=9)
    fig.savefig(ROOT/'figures/angular_modes.png',dpi=225,bbox_inches='tight')
    fig.savefig(ROOT/'figures/angular_modes.pdf',bbox_inches='tight')


if __name__=='__main__':
    main()
