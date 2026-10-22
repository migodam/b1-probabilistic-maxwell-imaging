"""Offline attribution of frozen material-risk gains; no model or physics work."""
from pathlib import Path
import json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r2.core import OUT,write,save
from p1_material.prior import MaterialPrior


def main():
    for path in ('FRESH_CHECKPOINT.json','L/evaluation/COMPLETE.json','L/trajectory_evaluation/COMPLETE.json'):
        if json.loads((OUT/path).read_text()).get('status')!='COMPLETE':raise RuntimeError('Frozen first-cohort online work and scoring must close first')
    dest=OUT/'L/saturation_risk';dest.mkdir(parents=True,exist_ok=False)
    cfg=json.loads((OUT/'configs/TARGET.json').read_text())
    geometry=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())
    prior=MaterialPrior(cfg,geometry['points']);roi=np.flatnonzero(prior.mask)
    parents=[f'R2_locked_{i:05d}' for i in range(256)]
    contributions=[];correction_energy=[];secant_energy=[];phase_counts=[];phase_gains=[];ranks=[]
    for parent in parents:
        with np.load(OUT/'L/contexts'/f'{parent}.npz') as z:q=z['q0'].copy();point=z['chi0'].copy()
        with np.load(OUT/'L/estimates'/f'{parent}.npz') as z:learned=z['fused'].copy()
        with np.load(OUT/'inputs/offline_only'/f'{parent}.npz') as z:truth=z['chi'].copy()
        D=prior.derivative(q);M=np.concatenate((D.real,D.imag),axis=0)/24
        U,s,_=np.linalg.svd(M,full_matrices=False);rank=int((s>max(s[0]*1e-8,1e-12)).sum());Q=U[:,:rank]
        error=truth-point;delta=learned-point
        e=np.r_[error.real,error.imag]/24;d=np.r_[delta.real,delta.imag]/24
        et=Q@(Q.T@e);dt=Q@(Q.T@d);en=e-et;dn=d-dt
        # Positive entries are material MSE improvements, adding exactly because
        # tangent and normal are orthogonal in the declared full-grid metric.
        gains=np.array([2*et@dt-dt@dt,2*en@dn-dn@dn])
        actual=float(np.mean(abs(error)**2)-np.mean(abs(truth-learned)**2))
        if not np.isclose(gains.sum(),actual,atol=1e-12,rtol=1e-10):raise ValueError('orthogonal risk identity failed')
        contributions.append(gains);correction_energy.append([dt@dt,dn@dn]);secant_energy.append([et@et,en@en]);ranks.append(rank)
        field=prior.decode(q)[1]['field'][roi];labels=np.where(field<=-prior.temperature,0,np.where(field>=prior.temperature,2,1))
        phase_counts.append(np.bincount(labels,minlength=3))
        cellgain=abs(error[roi])**2-abs((truth-learned)[roi])**2
        phase_gains.append([cellgain[labels==k].sum()/576 for k in range(3)])
    contributions=np.asarray(contributions);phase_gains=np.asarray(phase_gains);phase_counts=np.asarray(phase_counts)
    correction_energy=np.asarray(correction_energy);secant_energy=np.asarray(secant_energy)
    indices=np.random.default_rng(2026102292).integers(0,256,(2000,256))
    def summary(values):
        values=np.asarray(values);return dict(mean=float(values.mean()),parent_CI95=np.quantile(values[indices].mean(1),[.025,.975]).tolist())
    result=dict(status='COMPLETE',scope='Post-scoring exploratory attribution only; all original estimators remain frozen',
        parents=parents,new_RHS=0,GT_used_online=False,material_metric='complex full-grid MSE, equivalently [Re;Im]/24',
        tangent_rank_rule='s>max(smax*1e-8,1e-12)',rank_min=int(min(ranks)),rank_max=int(max(ranks)),
        mean_risk_gain=summary(contributions.sum(1)),
        tangent_risk_gain=summary(contributions[:,0]),outside_tangent_risk_gain=summary(contributions[:,1]),
        correction_energy_tangent=float(correction_energy[:,0].mean()),correction_energy_normal=float(correction_energy[:,1].mean()),
        correction_normal_energy_fraction=float(correction_energy[:,1].sum()/correction_energy.sum()),
        GT_secant_normal_energy_fraction=float(secant_energy[:,1].sum()/secant_energy.sum()),
        point_decoder_phase={name:dict(cells=int(phase_counts[:,k].sum()),parents_with_cells=int((phase_counts[:,k]>0).sum()),
            fullgrid_risk_gain_contribution=summary(phase_gains[:,k])) for k,name in enumerate(('host_plateau','transition','inclusion_plateau'))},
        interpretation='Normal-to-local-prior-tangent risk correction can be a valid unconstrained risk action. It is not a legal P1 sample, finite reachable path, posterior mode, or exact-fiber certificate.',
        phase_definition='Anchor decoder field, not GT-selected phase; host is original complex background, never air')
    save(dest/'RAW.npz',parent=np.array(parents),risk_gain_tangent_normal=contributions,correction_energy_tangent_normal=correction_energy,
        GT_secant_energy_tangent_normal=secant_energy,point_phase_counts=phase_counts,point_phase_gain_contributions=phase_gains,ranks=np.array(ranks))
    write(dest/'RESULTS.json',result)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10,4.3),constrained_layout=True)
    means=contributions.mean(0);ci=np.quantile(contributions[indices].mean(1),[.025,.975],axis=0)
    axes[0].bar(['Prior tangent','Outside local tangent'],means,yerr=np.array([means-ci[0],ci[1]-means]),capsize=5)
    axes[0].axhline(0,color='black',lw=1);axes[0].set(ylabel='Positive = full-grid material MSE improvement',title='Orthogonal risk attribution; all256')
    axes[1].bar(['Point host plateau','Point transition','Point inclusion plateau'],phase_gains.mean(0))
    axes[1].axhline(0,color='black',lw=1);axes[1].set(ylabel='Full-grid MSE gain contribution',title='Anchor decoder phase, not GT selection')
    axes[1].tick_params(axis='x',rotation=15)
    for ext in ('png','pdf'):fig.savefig(dest/('SATURATION_RISK.'+ext),dpi=170)
    plt.close(fig);write(dest/'COMPLETE.json',dict(status='COMPLETE',new_RHS=0,post_scoring_only=True))


if __name__=='__main__':main()
