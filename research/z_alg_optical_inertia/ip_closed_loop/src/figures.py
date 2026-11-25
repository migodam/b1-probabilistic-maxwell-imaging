"""Manuscript graphics. All marks are derived from saved fields or decisions."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy import stats
P=Path(__file__).resolve().parents[1];D=P/'maps/results';OUT=P/'figures'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,
 'axes.labelsize':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,
 'axes.linewidth':.65,'xtick.major.size':3,'ytick.major.size':3,'figure.dpi':120})
blue='#236a86';orange='#c57324';red='#a33741';green='#168375';pale='#e5eced'
def js(path):return json.loads(path.read_text())
def save(f,name):
    f.savefig(OUT/(name+'.pdf'),bbox_inches='tight',pad_inches=.04)
    f.savefig(OUT/(name+'.png'),dpi=165,bbox_inches='tight',pad_inches=.04)
    plt.close(f)
def title(ax,txt):ax.set_title(txt,loc='left',weight='semibold')
def spatial():
    fig,ax=plt.subplots(2,3,figsize=(7.1,4.8),layout='constrained')
    for row,(id,desc) in enumerate([('M96-R01','Annulus'),('M96-P01','Clipped pair')]):
        oid=f'ipcl-fov-f1__{id}__weak__7tx9rx__n03'
        t=np.load(P/'maps/offline'/(oid+'.npz'))
        x=t['truth'].reshape(96,96); ctr=np.load(D/(oid+'__b0.03125__witness.npz'))['centre']
        grid=np.linspace(-1,1,97)
        cover=np.load(D/(oid+'__b0.03125__primary_sector.npz'))['margin_sector'].reshape(96,96)
        color=ax[row,0].pcolormesh(grid,grid,x,cmap='viridis',vmin=0,vmax=4.25,edgecolors='none',rasterized=False)
        ax[row,1].pcolormesh(grid,grid,cover<=0,cmap=ListedColormap([pale,orange]),vmin=0,vmax=1,rasterized=False)
        ax[row,1].scatter((ctr[0]+.5)/48-1,(ctr[1]+.5)/48-1,marker='*',color=green,edgecolor='black',lw=.4,s=43)
        item=P/'adaptive'/oid
        if (item/'update.npz').exists():
            u=np.load(item/'update.npz');mask=np.zeros(96**2);mask[u['ids'][~(u['margins'][:,3]>0)]]=1
            ax[row,2].pcolormesh(grid,grid,mask.reshape(96,96),cmap=ListedColormap([pale,orange]),vmin=0,vmax=1,rasterized=False)
            if (P/'adaptive/updated_witnesses/summary.json').is_file():
                ax[row,2].scatter((ctr[0]+.5)/48-1,(ctr[1]+.5)/48-1,marker='*',color=green,edgecolor='black',lw=.4,s=43)
            ax[row,2].text(.03,.96,f'{int(mask.sum())} possible',va='top',transform=ax[row,2].transAxes)
        ax[row,1].text(.03,.96,f'{int((cover<=0).sum())} possible',va='top',transform=ax[row,1].transAxes)
        for a in ax[row]:
            a.set_aspect('equal');a.set_xticks([-1,0,1]);a.set_yticks([-1,0,1]);a.set_xlabel('x (m)')
        ax[row,0].set_ylabel(desc+'\ny (m)')
    for a,txt in zip(ax[0],['(a) Material contrast','(b) After seven illuminations','(c) After one additional beam']):title(a,txt)
    fig.colorbar(color,ax=ax[:,0],location='bottom',shrink=.85,aspect=25,label='Dielectric contrast')
    fig.legend(handles=[Patch(color=pale,label='Excluded centre'),Patch(color=orange,label='Retained centre'),
                         Line2D([],[],marker='*',ls='',color=green,label='Verified material')],loc='outside lower right',ncol=1,frameon=False,fontsize=7)
    save(fig,'location')
def population():
    f,(a,b)=plt.subplots(1,2,figsize=(7.0,2.85),layout='constrained')
    q=js(D/'summary.json');rows=sorted([r for r in q['rows'] if r['grid']==32 and r['actual_scene']=='weak'], key=lambda z:('solid ring double austria'.split().index(z['shape']),z['obs_id']))
    xx=np.arange(len(rows))
    for j,r in enumerate(rows):
        k=r['obs_id'];up=P/'adaptive'/k/'update.json'
        old=1024-r['sector_excluded']
        new=js(up)['after']['coherent'] if up.is_file() else old
        a.plot([j,j],[new,old],color='#83949b',lw=1)
        a.scatter(j,old,c='white',s=24,edgecolors=blue,zorder=3)
        a.scatter(j,new,c=red,marker='x',s=19,zorder=4)
    for v in (3.5,7.5,11.5):a.axvline(v,color='#e5e5e5')
    a.set_xticks([1.5,5.5,9.5,13.5],['Disk','Ring','Pair','Ring +\ndisks']);a.set_ylim(bottom=0)
    a.set_ylabel('Centres remaining out of 1,024')
    a.legend(handles=[Line2D([],[],marker='o',mfc='white',color=blue,ls='',label='Original cover'),Line2D([],[],marker='x',color=red,ls='',label='One extra coherent beam')],fontsize=7,frameon=False)
    title(a,'(a) Independent 32 by 32 specimens')
    wide=[]
    for j in range(1,5):
        oid=f'ipcl-fov-f1__M32-R{j:02d}__weak__7tx9rx__n03'
        counts=[]
        for v in ('primary_sector','broad_sector','broad_psd_sector'):
            ar=np.load(D/(oid+'__b0.03125__'+v+'.npz'))
            counts.append(int(np.sum(ar['margin_sector']>0)))
        wide.append(counts)
    z=np.transpose(wide)
    for j,(name,col) in enumerate([('Specified layer',blue),('Uncertain layer',orange),('Uncertain, positive only',red)]):
        b.plot(1+np.arange(4),100*z[j]/1024,'o-',color=col,ms=4,label=name)
    b.set_xticks(range(1,5));b.set_xlabel('Independent annulus');b.set_ylim(-4,105);b.set_ylabel('Excluded centres (%)')
    b.legend(frameon=False,fontsize=7,loc='center right',bbox_to_anchor=(1.,.3))
    title(b,'(b) Effect of layer uncertainty')
    save(f,'population')
def power():
    fig,(a,b)=plt.subplots(1,2,figsize=(7,3.1),layout='constrained')
    vals=[]; counts=[]
    roll=js(P/'beams/beam_summary.json')['secant_arm_trials']
    for n in (32,96):
        cuts=js(P/f'beams/paired_replay_power_{n}.json')['rows']
        ang=js(P/f'beams/angle_bound_power_{n}.json')['rows']
        proposals=js(P/f'pair/power_{n}/proposal.json')['rows']
        for i,v in cuts.items():
            lo=roll[f'confirm_{n}']['necessary_fluence_rows'][i]['necessary_fluence_lower']
            up=v['outcomes']['indefinite']['certified']['indefinite_quadratic']['sufficient_fluence_for_fixed_direction_upper']
            pos=v['outcomes']['positive']['certified']['positive_quadratic']['sufficient_fluence_for_fixed_direction_upper']
            vals.append((int(i),n,lo,pos,up,proposals[i]['fluence']))
            counts.append((ang[i]['maximum_midpoint_residual_over_new_radius_upper'],
                           v['outcomes']['indefinite']['certified']['indefinite_quadratic']['lower_distance_over_radius']/2))
    for j,(idx,n,lo,pos,up,e) in enumerate(vals):
        co=blue if n==32 else red
        a.plot([j,j],[lo,up],color=co,lw=1.6)
        a.scatter(j,pos,facecolor='white',edgecolor=co,s=20,zorder=3)
        a.scatter(j,e,c=co,s=14)
    a.axhline(.014,color='#555555',ls=':',lw=.85);a.set_yscale('log')
    a.set_xticks([0,3,7,11,12,13],['1','4','8','12','1','2'])
    a.axvline(11.5,color='#d9dede');a.set_xlabel('Twelve 32 by 32 cases, then two native cases')
    a.set_ylabel('Incident energy in old-shot units')
    a.legend(handles=[Line2D([],[],c=blue,label='Witness to indefinite bound'),Line2D([],[],marker='o',ls='',mfc='white',color=blue,label='Positive-state upper bound'),Line2D([],[],marker='o',ls='',color=blue,label='Chosen energy')],frameon=False,fontsize=6.7,loc='upper left')
    title(a,'(a) Necessary and sufficient energy')
    angle,coherent=np.array(counts).T
    b.scatter(np.arange(14),angle,color=orange,s=22,marker='v',label='Largest single-angle ambiguity radius')
    b.scatter(np.arange(14),coherent,color=blue,s=19,marker='o',label='Coherent whole-class half-gap')
    b.axhline(1,color='#454545',ls='--',lw=.9);b.axvline(11.5,color='#d9dede');b.set_xticks([0,3,7,11,12,13],['1','4','8','12','1','2'])
    b.set_xlabel('Independent cases, as in (a)'); b.set_ylabel('In units of new receiver error')
    b.legend(frameon=False,loc='lower left',fontsize=6.6);b.set_ylim(bottom=min(.3,angle.min()*.95))
    title(b,'(b) Can the two classes share a reading?')
    save(fig,'energy')
def lab():
    rows=list(csv.DictReader((P/'lab/held_ports.csv').open()))
    fig,ax=plt.subplots(1,2,figsize=(7,2.6),layout='constrained')
    for a,facility in zip(ax,['Fresnel','Manitoba']):
        groups=(['own_empty','common_empty'] if facility=='Fresnel' else ['cal_metal_3p5in_empty_first','held_metal_2in_empty_first','held_metal_2in_empty_last'])
        annot=['Own reference','Shared reference'] if facility=='Fresnel' else ['Large metal, first empty','Small metal, first empty','Small metal, last empty']
        for idx,(cal,lbl) in enumerate(zip(groups,annot)):
            values=sorted(float(r['held_ratio']) for r in rows if r['facility']==facility and r['calibration']==cal)
            a.step(values,np.arange(1,len(values)+1)/len(values),where='post',label=lbl,color=[blue,orange,green][idx],linestyle=['-','--',':'][idx],lw=1.5)
        a.axvline(1,color='#aaa',lw=.9);a.set_xlabel('Coherent / best individual (held arc)');a.set_ylim(0,1.02);a.legend(frameon=False,fontsize=6.7)
        a.set_ylabel('Fraction of paired conditions');title(a,facility)
    save(fig,'measured')
if __name__=='__main__':
    OUT.mkdir(exist_ok=True)
    spatial();population();power();lab()
