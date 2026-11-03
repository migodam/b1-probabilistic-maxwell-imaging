"""Review sealed saved arrays only. No physics or optimization imports."""
from pathlib import Path
from datetime import datetime, timezone
import json, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];R2=ROOT.parent/'research_max_struct_r2';CASE='obs_b20d5a55677d4eca'

def alternative_otsu(m):
    histogram,edges=np.histogram(m,bins=256);centers=(edges[:-1]+edges[1:])/2
    best=-np.inf;cut=centers[0]
    for k in range(256):
        nl=histogram[:k+1].sum();nr=histogram[k+1:].sum()
        if not nl or not nr:continue
        ml=np.dot(histogram[:k+1],centers[:k+1])/nl;mr=np.dot(histogram[k+1:],centers[k+1:])/nr
        score=nl*nr*(ml-mr)**2
        if score>best:best,cut=score,centers[k]
    return m>cut,float(cut)

def contract_otsu(m):
    histogram,edges=np.histogram(m,bins=256);centers=(edges[:-1]+edges[1:])/2
    count=np.cumsum(histogram);moment=np.cumsum(histogram*centers)
    product=count*(count[-1]-count);variance=np.full(256,-np.inf);nonempty=product>0
    variance[nonempty]=(moment[-1]*count[nonempty]-moment[nonempty]*count[-1])**2/product[nonempty]
    cut=centers[np.argmax(variance)]
    return m>cut,float(cut)

def main():
    tick,cpu=time.perf_counter(),time.process_time()
    scored=json.loads((ROOT/'results/DEVELOPMENT_SCORES.json').read_text());seal=json.loads((ROOT/'audit/DEVELOPMENT_SEAL.json').read_text())
    assert seal['count']==len(scored['rows'])==11
    with np.load(R2/'sealed_truth'/f'{CASE}.npz',allow_pickle=False) as f:truth=f['m'];mask=f['mask'].astype(bool).ravel();clean=f['Y_clean']
    with np.load(R2/'data_public'/f'{CASE}.npz',allow_pickle=False) as f:Y=f['Y'];sigma=float(f['sigma']);radius=float(f['noise_radius'])
    checks=[];images={}
    for row in scored['rows']:
        rec=json.loads((ROOT/row['receipt']).read_text());assert rec['status']=='SEALED_DEVELOPMENT_PREDICTION'
        assert rec['finished_utc'] < seal['checked_utc'] <= scored['scored_utc']
        with np.load(ROOT/row['prediction'],allow_pickle=False) as f:
            m=f['m'];yp=f['Y_prediction'];g=f['gradient'] if 'gradient' in f else None;images[row['arm']]=m.reshape(24,24).copy()
        segment,cut=contract_otsu(m);alternative,alternative_cut=alternative_otsu(m);alternative_iou=np.count_nonzero(alternative&mask)/np.count_nonzero(alternative|mask);iou=np.count_nonzero(segment&mask)/np.count_nonzero(segment|mask)
        n=np.sqrt(np.mean((m-truth)**2))/12;res=np.linalg.norm(yp-Y)/np.linalg.norm(Y);white=np.linalg.norm(yp-Y)/sigma
        cleanres=np.linalg.norm(yp-clean)/np.linalg.norm(clean);targetres=cleanres if row['arm']=='clean_native' else res
        mm=m.reshape(24,24);d=np.concatenate((np.diff(mm,axis=0).ravel(),np.diff(mm,axis=1).ravel()))/19.5
        phi=.5*targetres**2+1e-4*np.mean(d*d)+1e-4*np.mean(np.sqrt(d*d+.01**2)-.01)
        pg=float(np.max(np.abs(m-np.minimum(20,np.maximum(.02,m-g))))) if g is not None else None
        delta=dict(N=abs(n-row['normalized_material_error']),IoU=abs(iou-row['otsu_iou']),
            noisy_residual=abs(res-row['original_noisy_relative_residual']),target_residual=abs(targetres-row['target_relative_residual']),
            Phi_receipt=abs(phi-rec['physical_objective']),projected_gradient=abs(pg-rec['projected_gradient_inf']) if pg is not None else None)
        assert max(v for v in delta.values() if v is not None)<1e-12,(row['arm'],delta)
        assert row['original_noisy_joint_gate']==bool(white<=radius and n<=.15 and iou>=.7 and m.min()>=.02 and m.max()<=20)
        checks.append(dict(arm=row['arm'],absolute_deltas=delta,all_available_array_checks_passed=True,projected_gradient_check='PASS' if pg is not None else 'NOT_RUN_SAVED_GRADIENT_ABSENT',otsu_algebraic_sensitivity=dict(contract_threshold=cut,alternative_threshold=alternative_cut,alternative_iou=float(alternative_iou),iou_change=float(alternative_iou-iou))))
    readback=dict(status='PASS',classification='ARITHMETIC_ONLY_SEALED_ARRAY_READBACK',checked_utc=datetime.now(timezone.utc).isoformat(),
        arms=checks,sealed_count=seal['count'],joint_passes=sum(r['original_noisy_joint_gate'] for r in scored['rows']),
        otsu_note='Inherited histogram-center threshold and exact cumulative implementation retained. Equivalent between-class arithmetic can select a different tied-bin index and slightly alter IoU; all material and joint dispositions remain FAIL. No GT threshold tuning.',no_forward_model_imported=True,new_physical_LU=0,new_complex_RHS=0,hash_checks=False,
        timing_scope='Saved-array/receipt/truth read and independent arithmetic; plot work recorded separately',
        wall_seconds=time.perf_counter()-tick,cpu_seconds=time.process_time()-cpu)
    (ROOT/'audit/FINAL_ARITHMETIC_READBACK.json').write_text(json.dumps(readback,indent=2)+'\n')
    tick,cpu=time.perf_counter(),time.process_time()
    old=json.loads((R2/'results/BLIND_SCORES.json').read_text());base=next(r for r in old['rows'] if r['observation_id']==CASE and r['arm']=='native')
    with np.load(R2/base['prediction'],allow_pickle=False) as f:native=f['m'].reshape(24,24)
    byarm={r['arm']:r for r in scored['rows']}
    panels=[('Published truth',truth.reshape(24,24),None),('R2 native CIE',native,dict(N=base['amplitude_normalized_nrmse'],IoU=base['otsu_iou'],res=base['relative_residual'])),
        ('Strong-mode release',images['noisy_strong_release'],byarm['noisy_strong_release']),
        ('Clean-data CIE control',images['clean_native'],byarm['clean_native']),
        ('Passive-loss path',images['passive_path'],byarm['passive_path']),
        ('Full-wave paired rays',images['dressed_neighbor_rays'],byarm['dressed_neighbor_rays'])]
    fig,axes=plt.subplots(2,3,figsize=(11.8,8.7),constrained_layout=True)
    for ax,(title,m,row) in zip(axes.flat,panels):
        im=ax.imshow(m,vmin=.02,vmax=20,cmap='viridis',interpolation='nearest',origin='lower')
        ax.set_title(title,fontsize=12,pad=9);ax.set_xticks([]);ax.set_yticks([])
        if row:
            n=row.get('N',row.get('normalized_material_error'));iou=row.get('IoU',row.get('otsu_iou'));res=row.get('res',row.get('target_relative_residual'))
            target='clean target' if title=='Clean-data CIE control' else 'original noisy target'
            ax.set_xlabel(f'N = {n:.4f}   IoU = {iou:.4f}\nResidual = {res:.5f} ({target})',fontsize=9,labelpad=7)
        else:ax.set_xlabel('Already public R2 case; evaluation only',fontsize=9,labelpad=7)
    cbar=fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.86,pad=.025);cbar.set_label('Real material contrast m (same scale in every panel)',fontsize=10)
    fig.suptitle('R3 development: lower residual does not establish material recovery\n24 x 24 scalar-TM curved bar | all 11 development controls fail the inherited joint gate',fontsize=14)
    out=ROOT/'figures';out.mkdir(exist_ok=True);dest=out/'DEVELOPMENT_RECONSTRUCTIONS.png';fig.savefig(dest,dpi=170);plt.close(fig)
    receipt=dict(status='RENDERED_PENDING_VISUAL_REVIEW',classification='SAVED_ARRAY_PLOT_NO_PHYSICS',file=str(dest.relative_to(ROOT)),
        scale=[.02,20],interpolation='nearest',origin='lower',panels=[p[0] for p in panels],
        target_distinction='Clean-data control residual uses clean Y; all other panels use original noisy Y.',
        wall_seconds=time.perf_counter()-tick,cpu_seconds=time.process_time()-cpu,new_physical_LU=0,new_complex_RHS=0)
    (ROOT/'audit/RECONSTRUCTION_FIGURE_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(arithmetic=readback['status'],arms=len(checks),joint_passes=readback['joint_passes'],figure=str(dest))))
if __name__=='__main__':main()
