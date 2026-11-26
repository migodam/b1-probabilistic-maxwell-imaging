"""Rebuild all overview figures from recorded arrays, including unresolved trials."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
P=Path(__file__).resolve().parents[1];O=P/'figures';O.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'figure.dpi':135,'savefig.dpi':180,'pdf.fonttype':42,'axes.titlelocation':'left','axes.spines.right':False,'axes.spines.top':False})
def records(n):return list(csv.DictReader((P/'results'/n).open()))
def save(fig,name):
    for ext in ('png','pdf'):fig.savefig(O/f'{name}.{ext}',bbox_inches='tight',facecolor='white')
    plt.close(fig)
def masks():
    ids=['n32_c01','n96_c01'];fig,ax=plt.subplots(2,5,figsize=(12.4,5.1),layout='constrained')
    cols=['True contrast','Blind CIE-I','Same data: ROI zero','Same data: bright patch','What the data allow']
    rows=records('region_intervals.csv')
    for i,c in enumerate(ids):
        a=np.load(P/'imaging/results'/c/'stationary/baseline.npz');rr=a['rois'];x=a['chi'].reshape(rr.shape[1:]);truth=np.load(P/'ports/private'/(c+'.npz'))['material']
        e=json.loads((P/'results/extended'/(c+'.json')).read_text());hi=np.load(P/'results/extended'/(c+'.npz'));lo=np.load((P/e['zero_file']).with_suffix('.npz'))
        rs=rr[e['roi_index']-1];patch=np.zeros_like(rs);patch.ravel()[hi['patch']]=True
        for j,im in enumerate([truth,x,lo['material'].reshape(x.shape),hi['material'].reshape(x.shape)]):
            h=ax[i,j].imshow(im,origin='lower',cmap='magma',vmin=0,vmax=4.5,interpolation='nearest')
            ax[i,j].contour(rs,levels=[.5],origin='lower',colors='#43e4d8',linewidths=.7,extent=(-.5,x.shape[0]-.5,-.5,x.shape[0]-.5))
            if j==3:ax[i,j].contour(patch,levels=[.5],origin='lower',colors='white',linewidths=.65,extent=(-.5,x.shape[0]-.5,-.5,x.shape[0]-.5))
        state=np.zeros(x.shape,dtype=int)
        for rec in [v for v in rows if v['case']==c]:
            state[rr[int(rec['roi'])-1]]=2 if rec['width_exact']=='True' else 1
        ax[i,4].imshow(state,origin='lower',interpolation='nearest',vmin=0,vmax=2,cmap=ListedColormap(['#f4f3ed','#c8943a','#297e84']))
        for j in range(5):
            ax[i,j].set(xticks=[],yticks=[]);ax[i,j].spines[:].set_visible(False)
            if i==0:ax[i,j].set_title(cols[j],fontsize=10)
        ax[i,0].set_ylabel('32 x 32' if i==0 else '96 x 96')
        ax[i,2].set_xlabel('Regional mean = 0')
        ax[i,3].set_xlabel(f"{e['changed_patch_cells']} cells at contrast >= 2")
    fig.colorbar(h,ax=ax[:,:4],location='bottom',fraction=.05,pad=.03,aspect=65,label='Dielectric contrast (all eight image panels)')
    ax[1,4].legend(handles=[Patch(color='#297e84',label='Both 0 and 4.5 fit'),Patch(color='#c8943a',label='Remaining uncertainty'),Patch(color='#f4f3ed',label='Outside selected ROIs')],loc='upper center',bbox_to_anchor=(.5,-.10),fontsize=7.3,frameon=False)
    save(fig,'roi_witnesses')
def widths():
    a=records('region_intervals.csv');fig,axes=plt.subplots(1,2,figsize=(10.4,6.7),gridspec_kw={'width_ratios':[1.9,1]},layout='constrained')
    for ax,rows in zip(axes,([r for r in a if r['case'].startswith('n32')],[r for r in a if r['case'].startswith('n96')])):
        y=np.arange(len(rows))[::-1]
        for j,v in zip(y,rows):
            ax.plot([0,4.5],[j,j],lw=4,color='#d9d9d4',solid_capstyle='butt',zorder=1)
            if v['known_attained_min']:
                l=float(v['known_attained_min']);h=float(v['known_attained_max']);ax.plot([l,h],[j,j],lw=1.2,ls=':',color='#127980');ax.plot([l,h],[j,j],'o',ms=3,color='#127980')
            if v['continuous_from'] and float(v['continuous_to'])>float(v['continuous_from']):
                ax.plot([float(v['continuous_from']),float(v['continuous_to'])],[j,j],lw=1.8,color='#127980',zorder=2)
            ax.plot(float(v['true_peak']),j+.16,marker='v',ms=3.7,color='#963e2b')
        ax.set_yticks(y,[v['case'].split('_c')[1]+'.'+v['roi'] for v in rows],fontsize=7.6)
        ax.set_xlim(-.15,4.68);ax.set_xlabel('Regional peak contrast');ax.set_ylabel('Object.region');ax.set_title('32 x 32' if len(rows)>6 else 'Native 96 x 96')
        ax.grid(axis='x',alpha=.18)
    fig.suptitle('Which peak values still explain the original measurements?',fontsize=14)
    fig.supxlabel('Grey: guaranteed outer interval.  Teal dots: attainable endpoints.  Solid teal: verified continuous path.   Triangle: true peak.',fontsize=8.8)
    save(fig,'intervals')
def transfer():
    rows=records('grid_transfer.csv');fig,axes=plt.subplots(1,2,figsize=(10.4,3.4),layout='constrained')
    for i,(key,lab) in enumerate([('scoring_material','True materials'),('zero','ROI zero'),('extended','Bright patch')]):
        vals=[float(z['other_residual_over_r']) for z in rows if z['role']==key];xx=i+np.linspace(-.12,.12,len(vals))
        axes[0].scatter(xx,vals,color=['#242f43','#287a84','#b07120'][i],s=19)
    axes[0].set(xticks=[0,1,2],xticklabels=['True material','ROI zero','Bright patch'],ylabel='Residual / original noise radius',title='Changing grid in the same geometry');axes[0].axhline(1,color='#8d3c2f',ls='--',lw=1)
    chamber=records('manitoba_empty_repeat.csv')
    for i,q in enumerate([6,8,12,16]):
        vv=[float(r['min_radius_as_fraction_of_contrast'])*100 for r in chamber if int(r['q'])==q]
        axes[1].scatter(i+np.linspace(-.14,.14,len(vv)),sorted(vv),s=9,c='#27828b',alpha=.75)
    axes[1].axhline(3,color='#8d3c2f',ls='--',lw=1);axes[1].set(xticks=range(4),xticklabels=['6','8','12','16'],xlabel='Subsampled ports per side',ylabel='Minimum common empty radius (%)',title='Manitoba repeated blank measurements')
    save(fig,'transfer')
def mibox():
    root=P/'ports/results';a=json.loads((root/'mechanism_class_separation.json').read_text())
    b=json.loads((root/'mechanism_witness_arb.json').read_text())['conditions']
    fig,axes=plt.subplots(1,2,figsize=(10,3.6),layout='constrained')
    for freq,col in [(400,'#1e7d86'),(460,'#a96526')]:
        qs=[6,8,12,16]; vals=[float(b[f'{freq}_{q}']['separation_over_2r'][0]) for q in qs]
        axes[0].plot(qs,vals,'o-',label=f'{freq} MHz',color=col)
    axes[0].axhline(1,color='#585858',lw=.8,ls='--');axes[0].set(xticks=qs,title='Same two centered materials',xlabel='Transmitters and receivers per arc',ylabel='Full-wave distance / twice noise radius');axes[0].legend(frameon=False)
    for key,lab,c in [('nominal_margin_lower','Known impedance','#1e7d86'),('directional_calibration_margins_lower','0.5% impedance uncertainty','#a96526')]:
        xx=sorted({float(row['width']) for row in a['rows']});yy=[max(row[key] for row in a['rows'] if float(row['width'])==x and row[key] is not None) for x in xx]
        axes[1].semilogx(xx,1000*np.array(yy),'o-',color=c,label=lab)
    axes[1].axhline(0,color='#585858',ls='--',lw=.8);axes[1].set(xlabel='Simultaneous half-width on all 1,024 pixels',ylabel='Directional margin (10$^{-3}$ $\\sqrt{W}$)',title='8 x 8: from points to material boxes')
    axes[1].legend(frameon=False,fontsize=8);save(fig,'class_width')
if __name__=='__main__':
    masks();widths();transfer();mibox()
