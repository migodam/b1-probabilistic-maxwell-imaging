"""F figures from stored disclosed-development arrays only; no new solves."""
import json,time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[1]


def materials(out):
    real=np.load(ROOT/'results/families/F_analytic_real_h0.03.npz')
    pas=np.load(ROOT/'results/families/F_passive_complex_h0.03.npz')
    p=real['points'];a=real['axis'];n=len(a);iz=n//2
    fig=plt.figure(figsize=(12,7.3),constrained_layout=True)
    vals=[real['epsilon'][:,2,2].real,real['branch'][:,2,2].real]
    mask=np.linalg.norm(real['V'].reshape(-1,9),axis=1)>0
    vmax=max(x.max() for x in vals);axs=[]
    for i,(v,title) in enumerate(zip(vals,['Analytic real anchor: epsilon_zz','Analytic real branch: epsilon_zz'])):
        ax=fig.add_subplot(2,3,i+1,projection='3d');axs.append(ax)
        im=ax.scatter(*p[mask].T,c=v[mask],s=20,cmap='viridis',vmin=1,vmax=vmax)
        ax.set(xlabel='x / lambda',ylabel='y / lambda',zlabel='z / lambda',title=title)
        ax.set_box_aspect((1,1,1));ax.view_init(22,35)
    fig.colorbar(im,ax=axs,shrink=.6,label='Physical relative epsilon_zz')
    arrays=[(real['V'][:,2,2].real,'Real branch change: Delta epsilon_zz','magma',False),
            (pas['V'][:,0,2].real,'Passive branch: Re Delta epsilon_xz','RdBu_r',True),
            (pas['V'][:,2,2].imag,'Passive branch: Im Delta epsilon_zz','RdBu_r',True),
            (np.linalg.norm(pas['V'],axis=(1,2)),'Passive branch: tensor Frobenius change','magma',False)]
    # iz is a real saved center plane, not interpolation through z=0.
    for idx,(arr,title,cmap,symmetric) in zip([3,4,5,6],arrays):
        ax=fig.add_subplot(2,3,idx)
        sl=arr.reshape(n,n,n)[:,:,iz].T
        lim=max(abs(sl.min()),abs(sl.max()))
        kw=dict(vmin=-lim,vmax=lim) if symmetric else dict(vmin=0,vmax=sl.max())
        im=ax.imshow(sl,origin='lower',extent=[-.12,.12,-.12,.12],cmap=cmap,**kw)
        ax.set(xlabel='x / lambda',ylabel='y / lambda',title=title+'\nz='+f'{a[iz]:.3f}'+' lambda')
        fig.colorbar(im,ax=ax,shrink=.8)
    fig.suptitle('Disclosed material branches, h=.03 lambda: cell centers / one saved 2D plane\n'
                 'Same old-source fields in finite DDA; materials shown are not reconstructed images',fontsize=12)
    fig.savefig(out/'F_MATERIAL_BRANCHES.png',dpi=170)
    fig.savefig(out/'F_MATERIAL_BRANCHES.pdf');plt.close(fig)


def gaps(out):
    rows={}
    for kind in ['analytic_real','passive_complex']:
        rows[kind]=[json.loads((ROOT/f'results/families/F_{kind}_h{h:.2f}.json').read_text()) for h in [.06,.03]]
    specs=[('Old-source ports',72,lambda r:r['old_source_gap']['port_whitened']),
           ('A: added-receiver data only',144,lambda r:r['A_extra_receivers']['whitened_gap']),
           ('B: unitary old-source mix',72,lambda r:r['B_source_mix']['whitened_gap']),
           ('C: new-source prediction',36,lambda r:r['new_source_gap']['best_coherent_prediction'])]
    fig,axs=plt.subplots(2,4,figsize=(14,6.8),constrained_layout=True)
    metrics=[]
    for ri,kind in enumerate(['analytic_real','passive_complex']):
        for ci,(title,m,get) in enumerate(specs):
            ax=axs[ri,ci];ys=[float(get(r)) for r in rows[kind]];r=np.sqrt(m)
            ax.set_yscale('log');ax.set_ylim(1e-15,120);ax.set_xlim(-.4,1.4)
            ax.axhline(r,color='#737373',ls='--',lw=1)
            ax.axhline(2*r,color='#c44e52',ls=':',lw=1.5)
            positive=[i for i,y in enumerate(ys) if y>0]
            if positive:
                ax.plot(positive,[ys[i] for i in positive],color='#277da1',marker='o',lw=1.3)
            for i,y in enumerate(ys):
                if y==0:
                    # Annotation in axes coordinates: no fabricated positive data floor.
                    ax.text(i,.045,'exact 0',transform=ax.get_xaxis_transform(),ha='center',fontsize=9)
                else:
                    ax.annotate(f'{y:.3g}',(i,y),xytext=(0,7),textcoords='offset points',ha='center',fontsize=9)
            label='Analytic real' if ri==0 else 'Passive complex'
            ax.set(title=title+f'\nm={m}; r={r:.3g}; 2r={2*r:.3g}',xticks=[0,1],xticklabels=['h=.06','h=.03'])
            ax.set_ylabel(label+'\nwhitened port gap' if ci==0 else 'whitened port gap')
            ax.grid(axis='y',alpha=.2)
            metrics.append(dict(family=kind,experiment=title,complex_channels=m,noise_radius=float(r),pair_disjoint_threshold=float(2*r),gaps=ys))
    handles=[Line2D([],[],color='#277da1',marker='o',label='Stored nonzero gap'),
             Line2D([],[],color='#737373',ls='--',label='Noise-ball radius r = sqrt(m)'),
             Line2D([],[],color='#c44e52',ls=':',label='Two equal noise balls disjoint if gap > 2r')]
    fig.legend(handles=handles,loc='outside lower center',ncol=3,fontsize=9)
    fig.suptitle('A/B/C disclosed branch diagnostics: ||Delta E_ports|| / sigma, sigma=.001\n'
                 'Deterministic complex Frobenius balls; C is prediction only, unit ideal-channel coefficient power\n'
                 'Exact zeros are text annotations, omitted from logarithmic data; no display-floor values are inserted',fontsize=11)
    fig.savefig(out/'F_SOURCE_RECEIVER_GAPS.png',dpi=170)
    fig.savefig(out/'F_SOURCE_RECEIVER_GAPS.pdf');plt.close(fig)
    return metrics


def main():
    start=time.perf_counter();cpu=time.process_time();out=ROOT/'figures';out.mkdir(exist_ok=True)
    r=dict(status='RUNNING',LU=0,RHS=0,new_physical_states=0,
           inputs='results/families/F_{analytic_real,passive_complex}_h{0.06,0.03}.{json,npz}',
           scope='stored development-branch diagnostics only; no reconstruction or new acquisition',
           zero_plot_semantics='actual zeros omitted from log data and annotated in axes coordinates')
    prior=out/'F_PLOT_RECEIPT.json'
    if prior.exists():
        old=json.loads(prior.read_text())
        r['previous_plot_runs']=old.get('previous_plot_runs',[])+[
            {k:old.get(k) for k in ['status','wall_s','cpu_s','LU','RHS']}]
    try:
        materials(out);r['stored_gap_metrics']=gaps(out)
        r.update(status='COMPLETED',outputs=['F_MATERIAL_BRANCHES.png','F_MATERIAL_BRANCHES.pdf',
                                          'F_SOURCE_RECEIVER_GAPS.png','F_SOURCE_RECEIVER_GAPS.pdf'])
    except Exception as exc:r.update(status='FAILED',error=repr(exc));raise
    finally:
        r.update(wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu)
        (out/'F_PLOT_RECEIPT.json').write_text(json.dumps(r,indent=2,allow_nan=False))

if __name__=='__main__':main()
