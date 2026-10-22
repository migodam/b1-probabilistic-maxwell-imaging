"""Read-only old SMC energy/ancestry/decoder audit; zero new physics.

Uses only paid particles, original target decoder, paid DEV030 anchor and its
source-isolated online observation. Never imports Book/Target or opens GT.
"""
from pathlib import Path
import argparse,csv,json,sys
import numpy as np
from scipy.special import ndtr,logsumexp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p1_material.prior import MaterialPrior,saturation
from p3l.tempered import choose_beta
from p3r2.core import OUT,write,save


def stats(values,w):
    a=np.asarray(values,float);m=w@a
    return dict(mean=float(m),sd=float(np.sqrt(w@(a-m)**2)),minimum=float(a.min()),maximum=float(a.max()))


def correlation(a,b,w):
    a=np.asarray(a,float);b=np.asarray(b,float);a=a-w@a;b=b-w@b
    den=np.sqrt((w@a**2)*(w@b**2))
    return float(w@(a*b)/den) if den>0 else None


def decoder_metrics(prior,x,w):
    u=ndtr(x[:,64:]);field=prior.config['field_mean']+x[:,:64]@prior.A.T
    local=field[:,prior.mask];temperature=prior.temperature
    host=(local<=-temperature).mean(axis=1);pure=(local>=temperature).mean(axis=1)
    transition=1-host-pure
    fraction=saturation(field/temperature)*prior.mask[None,:]
    delta=(prior.real_interval[0]+np.diff(prior.real_interval)[0]*u[:,0])+1j*(prior.imag_interval[0]+np.diff(prior.imag_interval)[0]*u[:,1])
    material=prior.host+delta[:,None]*fraction;mean=w@material
    variance=np.einsum('i,ij->j',w,np.abs(material-mean)**2)
    centered=x[:,:64]-w@x[:,:64];eigen=np.linalg.eigvalsh(centered.T@(w[:,None]*centered))
    trace=float(eigen.sum());rank=float(trace**2/(eigen@eigen)) if eigen@eigen>0 else None
    return dict(native_property_real=stats(u[:,0],w),native_property_imag=stats(u[:,1],w),
        property_eta_real=stats(x[:,64],w),property_eta_imag=stats(x[:,65],w),
        property_real_tail_mass_01=float(w@((u[:,0]<.01)|(u[:,0]>.99))),
        property_imag_tail_mass_01=float(w@((u[:,1]<.01)|(u[:,1]>.99))),
        host_plateau_fraction_inside_ROI=stats(host,w),pure_plateau_fraction_inside_ROI=stats(pure,w),
        transition_fraction_inside_ROI=stats(transition,w),field_latent_covariance_effective_rank=rank,
        field_latent_covariance_trace=trace,material_total_variance_mean=float(variance.mean()),
        material_RMS_spread=float(np.sqrt(variance.mean()))),material,host,pure,transition


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--old',type=Path,default=ROOT.parent/'b1-p3l-local-posterior/research/p3l_local_posterior_imaging')
    parser.add_argument('--output',type=Path,default=OUT/'B/old_smc_mechanism')
    args=parser.parse_args();old=args.old.resolve();dest=args.output.resolve()
    if dest.exists(): raise RuntimeError('audit destination already exists; use a new explicit output')
    destination=old/'A/jobs/smc_DEV_030_mala_r0_memory_resume1'
    report=json.loads((destination/'RESULT.json').read_text());history={r['stage']:r for r in report['history']}
    cfg=json.loads((old/'configs/TARGET.json').read_text())
    points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
    prior=MaterialPrior(cfg,points)
    observation_path=old/'inputs/online/DEV_030.npz'
    with np.load(observation_path,allow_pickle=False) as z:
        if set(z.files)!={'observation','target_id','sigma'} or str(z['target_id'])!=cfg['target_id'] or float(z['sigma'])!=cfg['noise']['sigma']:
            raise ValueError('old source-isolated observation identity')
        y=z['observation'].copy();sigma=float(z['sigma'])
    anchorpath=old/'A/anchors/DEV_030.npz'
    with np.load(anchorpath,allow_pickle=False) as z:
        anchorq=z['q'].copy();anchorchi=z['chi'].copy();prediction=z['prediction'].copy()
    decode_error=float(np.linalg.norm(prior.decode(anchorq)[0]-anchorchi)/24)
    if decode_error>1e-12: raise ValueError('old anchor target decoder changed')
    anchorPhi=float(np.vdot(prediction-y,prediction-y).real/sigma**2)
    anchor_sourcePhi=np.sum(abs(prediction-y)**2,axis=0)/sigma**2
    paths=sorted((destination/'particles').glob('STAGE_*.npz'))
    rows=[];flat=[];material_means=[];stagearrays=[]
    for path in paths:
        stage=int(path.stem.split('_')[-1]);z=np.load(path,allow_pickle=False)
        x=z['x'];lw=z['logweights'];lw=lw-logsumexp(lw);w=np.exp(lw);phi=-z['loglikelihood'];roots=z['ancestor_ids'];beta=float(z['beta'])
        unique,counts=np.unique(roots,return_counts=True);masses=np.array([w[roots==r].sum() for r in unique]);entropy=-float(np.sum(masses*np.log(masses)))
        metrics,material,host,pure,transition=decoder_metrics(prior,x,w)
        hr=history.get(stage,{});mutation=hr.get('mutation',{})
        row=dict(stage=stage,beta=beta,particles=len(x),weight_ESS=float(1/(w@w)),
            distinct_prior_ancestors=len(unique),ancestor_entropy_effective_count=float(np.exp(entropy)),
            largest_ancestor_weight_mass=float(masses.max()),largest_ancestor_particle_fraction=float(counts.max()/len(x)),
            energy_Phi=stats(phi,w),energy_excess_over_paid_anchor=stats(phi-anchorPhi,w),
            mutation_acceptance=mutation.get('acceptance'),log_evidence=hr.get('log_evidence'),
            energy_property_real_correlation=correlation(phi,ndtr(x[:,64]),w),
            energy_transition_fraction_correlation=correlation(phi,transition,w),
            root_ids=unique.tolist(),**metrics)
        rows.append(row);material_means.append(w@material)
        flat.append(dict(stage=stage,beta=beta,weight_ESS=row['weight_ESS'],ancestors=len(unique),
            ancestor_entropy_effective_count=row['ancestor_entropy_effective_count'],
            Phi_mean=row['energy_Phi']['mean'],Phi_min=float(phi.min()),Phi_max=float(phi.max()),
            anchor_Phi=anchorPhi,mutation_acceptance=row['mutation_acceptance'],
            property_real_mean=metrics['native_property_real']['mean'],property_imag_mean=metrics['native_property_imag']['mean'],
            host_plateau_fraction=metrics['host_plateau_fraction_inside_ROI']['mean'],
            pure_plateau_fraction=metrics['pure_plateau_fraction_inside_ROI']['mean'],
            transition_fraction=metrics['transition_fraction_inside_ROI']['mean'],
            material_RMS_spread=metrics['material_RMS_spread']))
        stagearrays.append(dict(stage=stage,x=x.copy(),phi=phi.copy(),logweights=lw.copy(),ancestor_ids=roots.copy()))
    first=rows[0];last=rows[-1];collapse=next((r for r in rows if r['distinct_prior_ancestors']==1),None)
    initial=stagearrays[0];first_beta=choose_beta(0.,initial['logweights'],-initial['phi'],target_fraction=.8)
    finalmaterial=material_means[-1]
    conclusions=dict(
        prior_to_likelihood_concentration='Initial safe beta is set by likelihood-energy dispersion, before any decoder or proposal fit',
        genealogy='Resampling loses prior lineages irreversibly; finite weight ESS can recover without restoring lost lineages',
        energy_path='Terminal exact-likelihood energy is compared with the already-paid anchor for this identical observation; a much higher energy is direct poor-exploration evidence',
        decoder='Plateau and transition fractions expose C2 decoder saturation; correlations are descriptive and do not establish saturation as the sole cause',
        property='Native property and probit paths distinguish boundary concentration from an algorithmic change of prior',
        unsupported='No mathematical impossibility, global multimodality proof, calibrated posterior or causal attribution from one run')
    result=dict(status='READ_ONLY_MECHANISM_DIAGNOSTIC_COMPLETE',old_campaign_ledger_modified=False,new_full_states=0,new_RHS=0,GT_used=False,
        target_id=cfg['target_id'],online_observation=str(observation_path),paid_anchor=str(anchorpath),old_SMC_result=str(destination/'RESULT.json'),
        all_particle_snapshots=[str(p) for p in paths],paid_anchor_Phi=anchorPhi,paid_anchor_source_Phi=anchor_sourcePhi,
        paid_anchor_native_properties=anchorq[64:],anchor_decoder_RMS_error=decode_error,
        computed_initial_safe_tempering_step=first_beta,initial=first,terminal=last,
        first_single_ancestor_stage=None if collapse is None else dict(stage=collapse['stage'],beta=collapse['beta'],energy_Phi=collapse['energy_Phi']),
        terminal_mean_material_RMS_from_paid_anchor=float(np.linalg.norm(finalmaterial-anchorchi)/24),
        original_report_terminal_weight_ESS=report['weight_ess'],original_report_status=report['status'],
        mechanism_interpretation=conclusions,rows=rows,
        reference_qualification=False,caution='Historical exposed descriptive evidence; no new target evaluations or posterior reference qualification')
    dest.mkdir(parents=True)
    write(dest/'RESULT.json',result)
    with (dest/'ENERGY_ANCESTRY_DECODER_PATH.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(flat)
    save(dest/'MATERIAL_PATH.npz',stage=np.array([r['stage'] for r in rows]),beta=np.array([r['beta'] for r in rows]),
        mean=np.asarray(material_means),paid_anchor=anchorchi)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,7.5),constrained_layout=True)
    b=np.array([r['beta'] for r in rows]);a=axes.ravel()
    a[0].plot(b,[r['energy_Phi']['mean'] for r in rows],label='SMC weighted mean');a[0].plot(b,[r['energy_Phi']['minimum'] for r in rows],label='SMC minimum');a[0].axhline(anchorPhi,color='black',ls='--',label='Paid anchor');a[0].set_yscale('log');a[0].set_ylabel('Exact likelihood energy Phi');a[0].legend(fontsize=8)
    a[1].plot(b,[r['distinct_prior_ancestors'] for r in rows],label='Distinct prior roots');a[1].plot(b,[r['weight_ESS'] for r in rows],label='Weight ESS');a[1].set_ylabel('Count');a[1].legend(fontsize=8)
    a[2].plot(b,[r['native_property_real']['mean'] for r in rows],label='Native real property');a[2].plot(b,[r['native_property_imag']['mean'] for r in rows],label='Native imaginary property');a[2].set_ylabel('Weighted native property mean');a[2].legend(fontsize=8)
    for key,label in [('host_plateau_fraction_inside_ROI','Host plateau'),('pure_plateau_fraction_inside_ROI','Pure plateau'),('transition_fraction_inside_ROI','Transition band')]: a[3].plot(b,[r[key]['mean'] for r in rows],label=label)
    a[3].set_ylabel('Decoder cell fraction inside ROI');a[3].legend(fontsize=8)
    for axis in a: axis.set_xscale('symlog',linthresh=1e-5);axis.set_xlim(0,1);axis.set_xlabel('Likelihood power beta');axis.grid(alpha=.2)
    fig.suptitle('Old DEV030 SMC: historical energy and genealogy failure diagnostic')
    fig.savefig(dest/'MECHANISM_PATH.png',dpi=180);fig.savefig(dest/'MECHANISM_PATH.pdf');plt.close(fig)
    summary={k:result[k] for k in ('status','new_RHS','paid_anchor_Phi','computed_initial_safe_tempering_step','first_single_ancestor_stage')}
    summary.update(terminal_Phi_mean=last['energy_Phi']['mean'],terminal_weight_ESS=last['weight_ESS'],terminal_ancestors=last['distinct_prior_ancestors'],output=str(dest))
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__': main()
