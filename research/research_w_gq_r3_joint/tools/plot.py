"""Post-seal plots of actual saved fields, branches and prediction residuals."""
from pathlib import Path
import argparse,json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from generate import ROOT,maps,write_json
def case_label(case):return 'boundary-connected air defect [blind_cavity]' if case=='blind_cavity' else case
plt.rcParams.update({'font.size':10,'axes.titlesize':11,'figure.dpi':140,'savefig.dpi':160})
def image(ax,field,x,y,title,vmin,vmax,cmap='viridis',roi=None):
    p=ax.imshow(field,origin='lower',extent=[x[0],x[-1],y[0],y[-1]],vmin=vmin,vmax=vmax,cmap=cmap,interpolation='nearest',aspect='equal');ax.set_title(title);ax.set_xlabel(r'$x/\lambda_0$');ax.set_ylabel(r'$y/\lambda_0$')
    if roi:ax.set_xlim(roi[:2]);ax.set_ylim(roi[2:])
    return p
def plot_one(out):
    out=Path(out).resolve()
    scores=json.loads((out/'SCORES.json').read_text());case=scores['case'];repeat=scores['repeat'];truth=np.load(ROOT/'inputs/private'/case/'truth.npz');e=np.load(out/'endpoint.npz');actual=maps(truth['z']);estimated=maps(e['z']);x=actual['pixel_x'];y=actual['pixel_y'];roi=[-1,1,1.6,2.65]
    fig,ax=plt.subplots(2,2,figsize=(9,6.5),layout='constrained');vmax=max(actual['target_material'].real.max(),estimated['target_material'].real.max());lmax=max(actual['target_material'].imag.max(),estimated['target_material'].imag.max(),.01)
    for a,d,t,vm in [(ax[0,0],actual['target_material'].real,'Offline true target real permittivity',vmax),(ax[0,1],estimated['target_material'].real,'Measurement-selected target real permittivity',vmax),(ax[1,0],actual['target_material'].imag,'Offline true target loss at f=1',lmax),(ax[1,1],estimated['target_material'].imag,'Measurement-selected target loss at f=1',lmax)]:
        p=image(a,d,x,y,t,1 if 'permittivity' in t else 0,vm,roi=roi);fig.colorbar(p,ax=a,shrink=.85)
    s=scores['selected'];fig.suptitle(f"{case_label(case)} repeat {repeat}: position error {s['position_error']:.4f}, IoU {s['shape_IoU']:.3f}\nContinuous off-grid ellipse; private truth read after online seal",fontsize=12)
    fig.savefig(out/'target_material.png');plt.close(fig)
    fig,ax=plt.subplots(2,3,figsize=(13,4.6),layout='constrained');wroi=[-1.5,1.5,.5,1.25]
    for row,key,label,vm in [(0,'wall_eps','Wall real permittivity',6),(1,'wall_loss','Wall loss coefficient',.25)]:
        delta=estimated[key]-actual[key];dm=max(float(np.max(abs(delta))),.001)
        for j,d,title,cmap,lo,hi in [(0,actual[key],'Offline truth','viridis',1 if row==0 else 0,vm),(1,estimated[key],'Measurement-selected estimate','viridis',1 if row==0 else 0,vm),(2,delta,'Estimate minus truth','coolwarm',-dm,dm)]:
            p=image(ax[row,j],d,x,y,title,lo,hi,cmap,wroi);fig.colorbar(p,ax=ax[row,j],shrink=.85)
        ax[row,0].set_ylabel(label+'\n'+r'$y/\lambda_0$')
    fig.suptitle(f'{case_label(case)} r{repeat}: all six wall cells unknown online; no true support mask');fig.savefig(out/'estimated_wall.png');plt.close(fig)
    obs=np.load(ROOT/'inputs/online'/f'{case}_r{repeat}'/'observation.npz')['y'];pred=e['data'];clean=truth['clean'];fig,ax=plt.subplots(1,2,figsize=(11,3.5),layout='constrained')
    q=(pred-obs).ravel()/.01;qc=(pred-clean).ravel()/.01
    ax[0].plot(q.real,lw=.8,label='real');ax[0].plot(q.imag,lw=.8,label='imaginary');ax[0].axhline(0,color='k',lw=.5);ax[0].set_title('Predicted minus observed / sigma');ax[0].set_xlabel('Frequency/source/receiver channel');ax[0].legend()
    ax[1].plot(qc.real,lw=.8,label='real');ax[1].plot(qc.imag,lw=.8,label='imaginary');ax[1].axhline(0,color='k',lw=.5);ax[1].set_title('Predicted minus offline noiseless truth / sigma');ax[1].set_xlabel('Frequency/source/receiver channel');ax[1].legend();fig.suptitle(f'{case_label(case)} r{repeat}: actual saved full-wave predictions');fig.savefig(out/'predictions.png');plt.close(fig)
    branches=scores['branches']
    if branches:
        fig,ax=plt.subplots(figsize=(7,4),layout='constrained');ax.scatter([b['observation_white_RMS'] for b in branches],[b['position_error'] for b in branches],label='All retained branches');ax.scatter(s['observation_white_RMS'],s['position_error'],s=120,marker='*',color='red',label='Online selected');ax.axhline(.1,color='k',ls='--',lw=.8);ax.set_xlabel('Observation residual, real/imag RMS / sigma');ax.set_ylabel(r'Position error / $\lambda_0$');ax.legend();ax.set_title('Post-seal branch diagnosis; errors did not select endpoint');fig.savefig(out/'branches.png');plt.close(fig)
    return {k:s[k] for k in ['position_error','shape_IoU','epsilon_real_relative_error','loss_absolute_error','wall_epsilon_RMSE','wall_loss_RMSE','observation_white_RMS','clean_prediction_white_RMS']}|dict(case=case,repeat=repeat,method=scores.get('method',out.name),cohort=scores.get('cohort','DEVELOPMENT_EXPOSED' if case=='dev' else 'INDEPENDENT_FROZEN_BLIND'),joint_imaging_gate=s.get('joint_imaging_gate'),final_exact_projected_gradient=s.get('final_exact_projected_gradient'),stop=s.get('stop'),converged=s.get('converged'),objective=s.get('objective'),output=str(out.relative_to(ROOT)))
def overview(rows):
    methods=['FULL_GN','FIXED_CURRENT_CHART','ADAPTIVE_CURRENT_CHART'];cases=['blind_normal','blind_cavity','blind_strong'];lookup={(r['case'],r['method']):r for r in rows if r['repeat']==0 and r['case'] in cases}
    if not all((c,m) in lookup for c in cases for m in methods):return
    fig,axes=plt.subplots(3,4,figsize=(14,8),layout='constrained')
    labels=['Offline truth','Ordinary joint GN','Fixed forcing-current\ntangent chart (rank 16)','Adaptive forcing-current\ntangent chart (rank 16)']
    for i,c in enumerate(cases):
        actual=np.load(ROOT/'inputs/private'/c/'truth.npz');tm=maps(actual['z']);x=tm['pixel_x'];y=tm['pixel_y']
        image(axes[i,0],tm['target_material'].real,x,y,labels[0] if i==0 else case_label(c).replace('blind_',''),1,9,roi=[-1,1,1.6,2.65])
        for j,m in enumerate(methods,1):
            r=lookup[c,m];z=np.load(ROOT/r['output']/'endpoint.npz')['z'];em=maps(z);p=image(axes[i,j],em['target_material'].real,x,y,(labels[j]+'\n' if i==0 else '')+f"position {r['position_error']:.3f}; IoU {r['shape_IoU']:.2f}",1,9,roi=[-1,1,1.6,2.65])
        axes[i,0].set_ylabel(('boundary-connected air defect' if c=='blind_cavity' else c.replace('blind_','').replace('_',' '))+'\n'+r'$y/\lambda_0$')
    fig.colorbar(p,ax=axes.ravel().tolist(),shrink=.75,label='Target real permittivity at f=1');fig.suptitle('Three independent off-grid blind worlds, noise repeat 0\nAll six wall cells unknown; nonlinear forward F is full Maxwell for every method\nCharts approximate material tangent J; endpoints selected before offline scoring',fontsize=13);fig.savefig(ROOT/'figures/BLIND_TARGET_OVERVIEW.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4.5),layout='constrained')
    marker={'FULL_GN':'o','FIXED_CURRENT_CHART':'s','ADAPTIVE_CURRENT_CHART':'^','WALL_FIRST_FIXED':'D','KNOWN_WALL_ORACLE':'*','PRIOR_ONLY':'x'}
    for method in sorted(set(r['method'] for r in rows if r['case'] in cases)):
        rr=[r for r in rows if r['case'] in cases and r['method']==method];ax.scatter([r['position_error'] for r in rr],[r['shape_IoU'] for r in rr],marker=marker.get(method,'o'),s=65,label=method)
    ax.axvline(.1,color='k',ls='--',lw=.8);ax.axhline(.7,color='k',ls='--',lw=.8);ax.set_xlabel(r'Position error / $\lambda_0$');ax.set_ylabel('Shape IoU');ax.set_title('Post-seal imaging gates; three parent worlds\nNoise repeats within each parent are not independent worlds');ax.legend(fontsize=8,loc='best');fig.savefig(ROOT/'figures/BLIND_IMAGING_GATES.png');plt.close(fig)

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,action='append',default=[]);p.add_argument('--all',action='store_true');a=p.parse_args()
    outs=a.output
    if a.all:outs=sorted({p.parent for p in (ROOT/'results').glob('**/SCORES.json') if not any('INVALID' in q or 'REPAIR' in q or 'SUPERSEDED' in q or 'UNDERRESOLVED' in q or 'RANK8' in q for q in p.parts)})
    rows=[plot_one(o) for o in outs]
    if rows and a.all:
        with (ROOT/'results/IMAGING_SCORES.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        write_json(ROOT/'results/IMAGING_SCORES.json',rows)
        overview(rows)
    print(json.dumps(dict(plotted=len(rows),figures_per_endpoint='target_material,estimated_wall,predictions,branches if retained')))
if __name__=='__main__':main()
