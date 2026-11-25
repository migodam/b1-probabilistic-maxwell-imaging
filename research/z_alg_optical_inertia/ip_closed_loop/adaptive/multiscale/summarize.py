"""Summarize and render the completed adaptive centre search."""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE=Path(__file__).resolve().parent
IDS=('M32-D01','M32-R01','M32-P01','M32-A01','M96-R01','M96-P01')
COLOURS=('#faf9f5','#dad7cc','#89c4ba','#469489','#185a59')
DESCRIPTIONS=('Unresolved','Cached singleton exclusion','Two-centre parent',
              '2 × 2 parent','3 × 3 parent')


def read(path):
    return json.loads(path.read_text())


def plot(rows):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                         'pdf.fonttype':42,'ps.fonttype':42})
    fig,axes=plt.subplots(2,3,figsize=(10.5,7.7),layout='constrained')
    for ax,key in zip(axes.flat,IDS,strict=True):
        val=rows[key]
        side=val['side']
        arr=np.load(HERE/f"{val['oid']}.npz",allow_pickle=False)
        img=np.zeros((side,side),dtype=int)
        img[arr['excluded'].reshape(side,side)]=1
        for x,y,w,h,depth,code,margin,n,scale,weight in arr['records']:
            if code!=3:
                continue
            x,y,w,h=map(int,(x,y,w,h))
            k=2 if w*h==2 else 3 if w*h==4 else 4
            img[y:y+h,x:x+w]=k
        ax.imshow(img,origin='lower',interpolation='nearest',cmap=ListedColormap(COLOURS),vmin=0,vmax=4)
        x,y=val['old_physical_witness']['centre']
        ax.scatter([x],[y],s=48,facecolors='#eb7436',edgecolors='#5d392b',linewidths=.7,zorder=3)
        ax.set_xticks((0,side-1))
        ax.set_yticks((0,side-1))
        ax.set_title(f"{key} · {val['pruned_by_positive_parent']:,} centres rejected by parents",fontsize=9)
    handles=[Patch(facecolor=c,edgecolor='#bcb6aa',label=l) for c,l in zip(COLOURS,DESCRIPTIONS,strict=True)]
    handles.append(Line2D([],[],marker='o',color='none',markerfacecolor='#eb7436',markeredgecolor='#5d392b',
                          label='Retained physical witness'))
    fig.get_layout_engine().set(rect=(0,0.075,1,0.925))
    fig.legend(handles=handles,loc='lower center',ncol=3,frameon=False)
    fig.savefig(HERE/'multiscale_cover.pdf')
    fig.savefig(HERE/'multiscale_cover.png',dpi=155)
    plt.close(fig)


def main():
    papers=[read(HERE/f'panel_{side}.json') for side in (32,96)]
    rows={obs['oid'].split('__')[1]:obs for p in papers for obs in p['rows'].values()}
    overall={}
    for side in (32,96):
        v=[r for r in rows.values() if r['side']==side]
        strict=sum(r['strict_parent_certificates'] for r in v)
        indefinite=sum(node.get('strict_negative_inertia_parents',0)
                       for r in v for node in r['counts_by_tile'].values())
        overall[str(side)]=dict(objects=len(v),
                                possible_centres=sum(r['total_centres'] for r in v),
                                strict_parents=strict,strict_indefinite_parents=indefinite,
                                centres_pruned=sum(r['pruned_by_positive_parent'] for r in v),
                                cached_singletons_consulted=sum(r['fallback_singletons'] for r in v),
                                archived_singleton_exclusions=sum(r['archived_exhaustive_exclusions'] for r in v),
                                final_exclusions=sum(r['final_exclusions'] for r in v),
                                sum_wall_s=sum(r['wall_s'] for r in v),sum_cpu_s=sum(r['cpu_s'] for r in v))
    output=dict(overall=overall,objects=rows,
                interpretive_scope='All centres covered; group results verified prospectively on public old readings, with strict archived leaves only at unresolved branches. Cached fallback avoids recomputation and cannot support an end-to-end speedup claim against original batched runs.')
    (HERE/'summary.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    plot(rows)


if __name__=='__main__':main()
