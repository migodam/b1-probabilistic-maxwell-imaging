"""Offline-only evaluation and publication plots from saved Research F arrays.

No solver imports, no LU reads, no physics calls. Truth and N3 are used only
after F1 is closed, and never enter initialization, optimization, or gating.
"""
from pathlib import Path
import csv, json, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'research/research_f_intrinsic_transport'
FIG=OUT/'figures'
N3FILE=ROOT.parents[1]/'A26/a26-om-feedback-prior-completion/results/a26/evaluation/261009/PREDICTIONS.npz'


def pack(z): return np.concatenate((z.real.T,z.imag.T),axis=1).ravel()
def norm(z): return float(np.linalg.norm(z))
def load(path,key):
    with np.load(path) as z: return z[key]
def save(fig,name):
    fig.savefig(FIG/f'{name}.png',dpi=180,bbox_inches='tight')
    fig.savefig(FIG/f'{name}.pdf',bbox_inches='tight')
    plt.close(fig)


def main():
    start=time.perf_counter(); FIG.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
        'axes.spines.right':False,'figure.facecolor':'#fafaf7','axes.facecolor':'#fafaf7',
        'axes.titleweight':'bold','savefig.facecolor':'#fafaf7'})
    summary=json.loads((OUT/'results/F1/SUMMARY.json').read_text())
    records=[]; bundles=[]; allfields=[]
    with np.load(N3FILE) as z:
        n3ids=z['case_id']; n3=z['N3']; n3truth=z['truth']
    for s in summary:
        case=s['case']; dest=OUT/'results/F1'/case
        hist=json.loads((dest/'HISTORY.json').read_text())
        truth=load(OUT/'inputs/offline_only'/f'{case}.npz','truth')
        initial=load(ROOT/hist[0]['array'],'chi')
        final=load(ROOT/s['final_array'],'chi')
        observation=load(OUT/'inputs/online'/f'{case}.npz','observation')
        sigma=float(load(OUT/'inputs/online'/f'{case}.npz','sigma'))
        truthprediction=load(OUT/'inputs/offline_only'/f'{case}.npz','truth_prediction')
        idx=np.flatnonzero(n3ids==case)
        nn=n3[idx[0]] if len(idx)==1 else None
        match=float(np.max(abs(n3truth[idx[0]]-truth))) if nn is not None else None
        if match is not None and match>1e-12: raise RuntimeError('same-case N3 truth mismatch')
        fields=[initial,final,nn,truth]
        rec=dict(case=case,exposed_validation=True,same_parent_across_backgrounds=True,
            conditioning=s['known_conditional_geometry'],
            stationarity_pass=s['qualified'],final_optimality=s['optimality'],
            final_residual=s['final_residual'],noise_threshold=s['tau'],
            noise_ball_pass=s['final_residual']<=s['tau'],
            truth_cached_residual=norm(pack(truthprediction-observation))*np.sqrt(2)/sigma,
            N3_pixel_cache=str(N3FILE),N3_row=int(idx[0]) if len(idx)==1 else None,
            N3_cached_truth_max_abs_error=match,N3_full_forward_residual='MISSING_NOT_RECOMPUTED',
            metrics={},posterior_samples=False,F2_status='NOT_RUN_GATE2_FAILED')
        for method,x in zip(['initial','physics_final','N3','offline_truth'],fields):
            if x is None:rec['metrics'][method]='MISSING';continue
            err=x-truth
            rec['metrics'][method]=dict(complex_MAE=float(np.mean(abs(err))),complex_RMS=norm(err)/24,
                real_MAE=float(np.mean(abs(err.real))),imaginary_MAE=float(np.mean(abs(err.imag))))
            allfields.append(x)
        records.append(rec);bundles.append((s,hist,fields))
    limits=[(min(float(x.real.min()) for x in allfields),max(float(x.real.max()) for x in allfields)),
            (min(float(x.imag.min()) for x in allfields),max(float(x.imag.max()) for x in allfields))]
    titles=['Prior initial','Physics final','Cached N3','Offline truth']
    for s,hist,fields in bundles:
        fig,axs=plt.subplots(2,4,figsize=(12.4,6.2),layout='constrained')
        for row in range(2):
            for col,(title,x) in enumerate(zip(titles,fields)):
                ax=axs[row,col]
                if x is None:ax.text(.5,.5,'MISSING',ha='center',va='center');continue
                vals=(x.real if row==0 else x.imag).reshape(24,24)
                im=ax.imshow(vals,origin='lower',cmap='viridis',vmin=limits[row][0],vmax=limits[row][1],interpolation='nearest')
                ax.set_xticks([]);ax.set_yticks([])
                if row==0:ax.set_title(title)
                if col==0:ax.set_ylabel('Re(chi)' if row==0 else 'Im(chi)')
            fig.colorbar(im,ax=axs[row,:],shrink=.78,pad=.01,label='Re(chi)' if row==0 else 'Im(chi)')
        bg=s['case'].split('_')[1]
        fig.suptitle(f'{bg.upper()} | exposed conditional reconstruction | stationarity gate FAILED',fontsize=14)
        fig.supxlabel('Shared colors across methods and cases. Truth/N3 offline only; conditioning differs.\nN3 forward residual missing. No posterior samples; F2 not run.',fontsize=10)
        save(fig,f'01_material_{bg}')
    colors=['#2a6f97','#d97732']
    fig,axs=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    for (s,hist,fields),color in zip(bundles,colors):
        bg=s['case'].split('_')[1]; ids=np.arange(1,len(hist)+1)
        axs[0].semilogy(ids,[h['residual'] for h in hist],'-o',ms=3,color=color,label=bg)
        axs[1].semilogy(ids,[h['likelihood_cost'] for h in hist],'-o',ms=3,color=color,label=bg)
    axs[0].axhline(summary[0]['tau'],color='#5a5a5a',ls='--',lw=1,label='Frozen noise threshold')
    axs[0].set(title='Data fit improves',xlabel='Paid state evaluation',ylabel='Whitened residual')
    axs[1].set(title='Likelihood objective',xlabel='Paid state evaluation',ylabel='Phi = residual² / 2')
    axs[0].legend(fontsize=8)
    vals=[s['optimality'] for s in summary]
    axs[2].bar([s['case'].split('_')[1] for s in summary],vals,color=colors,width=.55)
    axs[2].set_yscale('log');axs[2].axhline(.001,color='#c13f3f',ls='--',label='Required 0.001')
    axs[2].set(title='Final stationarity: both fail',ylabel='Bound-scaled first-order optimality')
    for i,v in enumerate(vals):axs[2].text(i,v*1.15,f'{v:.5g}',ha='center')
    axs[2].legend(fontsize=8);axs[2].set_ylim(.0005,max(vals)*2)
    fig.suptitle('F1: residual-qualified endpoints, without a qualified stationary anchor',fontsize=14)
    fig.supxlabel('24 new states / 960 complex RHS total. Only final optimality is available; no intermediate KKT is invented.',fontsize=10)
    save(fig,'02_F1_convergence_and_gate')
    with (OUT/'CURVATURE_DECOMPOSITION.csv').open() as f:
        geometry=[r for r in csv.DictReader(f) if r['derivative']=='analytic_fixed_max']
    fig,axs=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
    x=np.arange(2); width=.23
    for j,(column,label,color) in enumerate([
        ('C_coordinate_prior_unit','Coordinate','#2a6f97'),
        ('C_prior_prior_unit','Prior normal','#d97732'),
        ('historical_Maxwell_norm_prior_unit','Maxwell norm only*','#979797')]):
        bars=axs[0].bar(x+(j-1)*width,[float(r[column]) for r in geometry],width,color=color,label=label)
        if j==2:
            for bar in bars:bar.set_hatch('///')
    axs[0].set_yscale('log');axs[0].set_xticks(x,['TRAIN weak','TRAIN ring'])
    axs[0].set(title='F0: three second-order terms',ylabel='Whitened response / prior-coordinate speed²')
    axs[0].legend(fontsize=8)
    for i,r in enumerate(geometry):
        lo=float(r['historical_norm_only_intrinsic_cokernel_lower']);hi=float(r['historical_norm_only_intrinsic_cokernel_upper'])
        mid=float(r['Qout_prior_prior_unit'])
        axs[1].errorbar(i,mid,yerr=[[mid-lo],[hi-mid]],fmt='o',capsize=7,color=colors[1-i],lw=2)
        axs[1].text(i+.06,mid,f'[{lo:.4g}, {hi:.4g}]',fontsize=8,va='center')
    axs[1].set_yscale('log');axs[1].set_xticks(x,['TRAIN weak','TRAIN ring']);axs[1].set_xlim(-.35,1.7)
    axs[1].set(title='Intrinsic cokernel norm: bounded',ylabel='||Qout (C_Maxwell + C_prior)||')
    for i,r in enumerate(records):
        names=['Initial','Physics','N3']; metrics=r['metrics']
        ys=[metrics[k]['complex_RMS'] for k in ['initial','physics_final','N3']]
        axs[2].plot(names,ys,'-o',color=colors[i],label=r['case'].split('_')[1])
    axs[2].set(title='F1: offline image risk',ylabel='Complex material RMS to offline truth');axs[2].legend(fontsize=8)
    fig.suptitle('Coordinate acceleration dominates one TRAIN chart; a real normal response remains',fontsize=13)
    fig.supxlabel('*Historical Maxwell scalar norm only: its vector is MISSING. Bounds use triangle inequalities.\nBoth F0 tangent and restricted-response ranks are 7; no nonzero exact prior tangent kernel. F0 TRAIN and F1 validation cases differ.\nOffline risk uses two exposed backgrounds of one parent; N3 and physics conditioning differ.',fontsize=9)
    save(fig,'03_geometry_and_offline_risk')
    result=dict(scope='OFFLINE_ONLY_AFTER_F1_NO_SOLVER_IMPORTS',records=records,
        shared_color_limits=dict(real=limits[0],imaginary=limits[1]),
        new_full_states=0,new_complex_rhs=0,new_lu_factorizations=0,
        caveats=['Same material parent observed under two backgrounds; not two independent parents.',
            'Physics uses known family and global geometry; N3 has a different learned-input/conditioning regime.',
            'Two exposed cases do not establish average risk superiority.',
            'F1 final endpoints pass the noise threshold but fail the frozen stationarity gate.',
            'F2 was not run. No posterior mass or sampling claim.'],
        wall_seconds=time.perf_counter()-start,
        figures=[p.name for p in sorted(FIG.glob('*'))])
    (OUT/'OFFLINE_EVALUATION.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status='COMPLETE',figures=len(result['figures']),offline_metrics=str(OUT/'OFFLINE_EVALUATION.json'),new_states=0,new_RHS=0)))


if __name__=='__main__':main()
