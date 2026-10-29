"""Plot frozen Z-ALG evidence only; does not run a Maxwell solve or quadrature.

Usage: python make_two_cell_structure.py [--e4-json PATH] [--output-dir DIR]
The two input JSON files are existing experiment records. This script only reads
their saved values and draws the exact sum/product curves used in the proof.
"""
from pathlib import Path
import argparse
import json
import os
import tempfile

os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'zalg_matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np


HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--e4-json',type=Path,default=ROOT/'zalgo_probability/e4_two_cell_posterior.json')
parser.add_argument('--witness-json',type=Path,default=ROOT/'zalgo_maxwell/two_cell_witness.json')
parser.add_argument('--output-dir',type=Path,default=HERE)
args=parser.parse_args()
e4=json.loads(args.e4_json.read_text(encoding='utf-8'))
witness=json.loads(args.witness_json.read_text(encoding='utf-8'))
latest=e4['results']['64']
event_values=[latest['original_box']['event_x1_gt_x2'],
              latest['ordered_quotient']['event_x1_gt_x2'],
              latest['ordered_quotient']['uniform_branch_wrong_control']['event_x1_gt_x2']]

plt.rcParams.update({
    'font.family':'DejaVu Sans',
    'font.size':11,
    'axes.titlesize':12,
    'axes.labelsize':11,
    'xtick.labelsize':10,
    'ytick.labelsize':10,
    'legend.fontsize':9.5,
    'axes.spines.top':False,
    'axes.spines.right':False,
    'axes.linewidth':0.8,
    'savefig.facecolor':'white',
    'pdf.fonttype':42,
    'ps.fonttype':42,
})
navy='#234F78'
teal='#00857D'
amber='#D07A20'
purple='#76528D'
grey='#6B7280'
light='#D9DEE4'

fig=plt.figure(figsize=(7.2,6.2))
grid=fig.add_gridspec(2,2,width_ratios=[.85,1.15],height_ratios=[1.2,1.0],
                     left=.105,right=.97,bottom=.115,top=.93,hspace=.62,wspace=.48)
axes=[fig.add_subplot(grid[0,0]),fig.add_subplot(grid[0,1]),fig.add_subplot(grid[1,:])]

# A. Preserve the actual physical geometry and equal x/y length scale.
ax=axes[0]
ax.set_title('(A) Fixed labelled geometry',pad=12,fontsize=11)
ax.axvline(0,color=light,lw=1,ls=(0,(3,3)),zorder=0)
for index,(center,color) in enumerate(zip(witness['centers'],[navy,purple]),1):
    ax.add_patch(Circle(center,witness['radius'],facecolor=color,edgecolor=color,lw=1))
    offset=-0.085 if index==1 else 0.085
    ax.annotate(rf'$x_{index}$',xy=center,xytext=(center[0]+offset,0.17),
                ha='center',va='bottom',color=color,fontsize=12,
                arrowprops={'arrowstyle':'-','lw':0.7,'color':color,'shrinkB':5})
for index,(x,y) in enumerate(witness['sources'],1):
    ax.scatter(x,y,marker='^',s=62,color=amber,edgecolors='white',linewidths=0.6,zorder=3)
    ax.annotate(f'S{index}',(x,y),xytext=(10,0),textcoords='offset points',ha='left',va='center',color=amber)
for index,(x,y) in enumerate(witness['receivers'],1):
    ax.scatter(x,y,marker='s',s=48,color=teal,edgecolors='white',linewidths=0.6,zorder=3)
    ax.annotate(f'R{index}',(x,y),xytext=(10,0),textcoords='offset points',ha='left',va='center',color=teal)
ax.text(-0.47,-0.16,'Two disks\n'+r'$a=0.04\lambda$',ha='left',va='top',fontsize=9,color=grey)
ax.set_xlim(-0.52,0.52)
ax.set_ylim(-1.27,1.37)
ax.set_aspect('equal',adjustable='box')
ax.set_xticks([-0.4,0,0.4])
ax.set_yticks([-1,-0.5,0,0.5,1])
ax.set_xlabel(r'$x/\lambda$')
ax.set_ylabel(r'$y/\lambda$',labelpad=3)
ax.tick_params(length=3)

# B. Algebraic curves and the two real material states; no field evaluation.
ax=axes[1]
ax.set_title('(B) Exact material fiber',pad=12,fontsize=11)
x=np.linspace(0.5,4,700)
ax.plot(x,4-x,color=navy,lw=1.9,label=r'$x_1+x_2=4$')
ax.plot(x,3/x,color=teal,lw=1.9,label=r'$x_1x_2=3$')
states=np.asarray([[1,3],[3,1]])
ax.scatter(states[:,0],states[:,1],s=60,color=navy,edgecolors='white',lw=1,zorder=4)
ax.annotate('(1, 3)',(1,3),xytext=(7,6),textcoords='offset points',ha='left',va='bottom',fontsize=10,color=navy)
ax.annotate('(3, 1)',(3,1),xytext=(7,7),textcoords='offset points',ha='left',va='bottom',fontsize=10,color=navy)
mean=np.asarray(e4['exact_conditional_sanity_at_pair_1_3']['conditional_mean'])
ax.scatter(*mean,marker='X',s=78,color=amber,edgecolors='white',lw=0.6,zorder=5)
ax.annotate('Conditional mean\n(2.5, 1.5)\nNot on the fiber',xy=mean,xytext=(3.02,2.38),
            fontsize=9,color=amber,ha='center',va='center',
            bbox={'facecolor':'white','edgecolor':'none','alpha':.92,'pad':1.5},
            arrowprops={'arrowstyle':'-','color':amber,'lw':0.8,'shrinkB':6})
ax.set_xlim(.5,4)
ax.set_ylim(.5,4)
ax.set_aspect('equal',adjustable='box')
ax.set_xticks([.5,1,2,3,4])
ax.set_yticks([.5,1,2,3,4])
ax.set_xlabel(r'$x_1$')
ax.set_ylabel(r'$x_2$',labelpad=3)
ax.grid(color=light,lw=.6,alpha=.55)
ax.set_axisbelow(True)
ax.legend(loc='upper right',frameon=False,handlelength=1.5,borderaxespad=.3,labelspacing=.5)

# C. Frozen n=64 quadrature values; correct posterior and quotient agree.
ax=axes[2]
ax.set_title('(C) Finite-noise event probability',pad=12,fontsize=11,loc='left')
bars=ax.bar(np.arange(3),event_values,width=.62,color=[navy,teal,amber],zorder=3)
for bar,value in zip(bars,event_values):
    ax.text(bar.get_x()+bar.get_width()/2,value+.018,f'{value:.6f}',ha='center',va='bottom',fontsize=10)
ax.set_ylim(0,.88)
ax.set_yticks([0,.25,.5,.75])
ax.set_xticks(np.arange(3),['Original\nposterior','Correct\nquotient','Uniform\nbranches'])
ax.tick_params(axis='x',length=0,pad=8)
ax.set_ylabel(r'$P(x_1>x_2\mid y)$',labelpad=4)
ax.yaxis.grid(color=light,lw=.7,alpha=.75,zorder=0)
ax.text(.02,.98,r'$\rho(x)\propto x_1,\quad\sigma=10^{-3}$',transform=ax.transAxes,
        ha='left',va='top',fontsize=10,color=grey)

args.output_dir.mkdir(parents=True,exist_ok=True)
fig.savefig(args.output_dir/'two_cell_structure.pdf',bbox_inches='tight',pad_inches=.08)
fig.savefig(args.output_dir/'two_cell_structure.png',dpi=240,bbox_inches='tight',pad_inches=.08)
plt.close(fig)
print(json.dumps({'png':str(args.output_dir/'two_cell_structure.png'),
                  'pdf':str(args.output_dir/'two_cell_structure.pdf'),
                  'event_values':event_values},indent=2))
