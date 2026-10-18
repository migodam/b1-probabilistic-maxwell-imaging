"""Execute the new P1 material campaign. No electromagnetic backend imports."""
from pathlib import Path
import json,sys,time,subprocess,collections
import numpy as np
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p1_material.prior import MaterialPrior,topology,structural_change
from p1_material.proposals import normalized_direction,frozen_latent_step,material_retraction
from p1_material.target import MaxwellTarget,BudgetDenied,write_json
OUT=ROOT/'research/p1_material_target_20261018'


def wilson(k,n):
    z=1.959963984540054;p=k/n;den=1+z*z/n;mid=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [float(mid-half),float(mid+half)]


def main():
    if (OUT/'MATERIAL_STAGE_CLOSED.json').exists():raise RuntimeError('preserve material-stage results; choose a new protocol version instead of overwrite')
    start=time.perf_counter();cfg=json.loads((OUT/'configs/PROTOCOL.json').read_text());points=json.loads((OUT/'inputs/geometry.json').read_text())['points'];P=MaterialPrior(cfg,points)
    commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True,cwd=ROOT).strip();metadata=[];all_q=[];all_x=[];all_f=[];groups={}
    existing=(OUT/'results/MATERIAL_MANIFEST.json').exists()
    for group,seed,n in [('calibration',cfg['calibration_seed'],128),('dev',cfg['development_seed'],32),('material_check',cfg['check_seed'],32)]:
        cache=OUT/f'inputs/offline_materials/{group}.npz'
        if existing:
            with np.load(cache) as z:q=z['latent'].copy()
        else:q=P.draw(np.random.default_rng(seed),n)
        xx=[];ff=[]
        for i,latent in enumerate(q):
            x,a=P.decode(latent);xx.append(x);ff.append(a['fraction'])
            events={str(t):topology(a['fraction'],t) for t in cfg['topology_thresholds']}
            metadata.append(dict(parent_id=f'P1_{group}_{i:03}',group=group,seed=seed,index=i,log_latent_prior=P.logprior(latent),latent_squared_norm=float(np.dot(latent[:64],latent[:64])),inclusion_increment=[a['inclusion_increment'].real,a['inclusion_increment'].imag],fraction_mean=float(np.mean(a['fraction'][P.mask])),host_cells=int(np.sum((a['fraction']==0)&P.mask)),pure_inclusion_cells=int(np.sum(a['fraction']==1)),mixed_cells=int(np.sum((a['fraction']>0)&(a['fraction']<1))),events=events,Maxwell_evaluated=False,posterior_mass=None))
        xx=np.array(xx);ff=np.array(ff);groups[group]=(q,xx,ff)
        if not existing:np.savez_compressed(cache,latent=q,chi=xx,fraction=ff,points=P.points)
        all_q.extend(q);all_x.extend(xx);all_f.extend(ff)
    all_q=np.array(all_q);all_x=np.array(all_x);all_f=np.array(all_f)
    # All 32 development noise vectors fixed before any forward calculation.
    rng=np.random.default_rng(cfg['noise_seed']);sigma=cfg['noise']['sigma'];noise=(rng.normal(size=(32,32,8))+1j*rng.normal(size=(32,32,8)))*sigma/np.sqrt(2)
    if not existing:
        np.savez_compressed(OUT/'inputs/noise_draws.npz',noise=noise,sigma=sigma)
        write_json(OUT/'results/MATERIAL_MANIFEST.json',dict(frozen_code_commit=commit,target_id=cfg['target_id'],rows=metadata,all_material_only=True,all_inspected_development=True,not_Maxwell_validation=True))
    # Material calculus, independent of any Maxwell Jacobian or truth fit.
    derivative=[];spectrum=[];rng=np.random.default_rng(2026101811)
    for i in range(6):
        q=groups['calibration'][0][i].copy();q[64:]=.5;v=rng.normal(size=66);v/=np.linalg.norm(v);D=P.derivative(q);exact=D@v;errors=[]
        for h in [1e-3,3e-4,1e-4]:
            fd=(P.decode(q+h*v)[0]-P.decode(q-h*v)[0])/(2*h)
            errors.append(dict(h=h,relative_error=float(np.linalg.norm(fd-exact)/max(np.linalg.norm(exact),1e-300))))
        sv=np.linalg.svd(np.r_[D.real,D.imag]/24,compute_uv=False);rank=int(sum(sv>sv[0]*1e-8)) if sv[0]>0 else 0
        derivative.append(dict(calibration_index=i,property_coordinates_fixed_at_midpoint_for_calculus=True,steps=errors));spectrum.append(dict(index=i,rank=rank,singular_values=sv,zero_fraction_cells=int(sum(P.decode(q)[1]['fraction']==0))))
    # Actual material-distance comparison on fixed prior draws, not inverse anchors.
    rng=np.random.default_rng(2026101812);retractions=[];retraction_arrays={}
    for i in range(6):
        q=groups['calibration'][0][i]
        try:v,w=normalized_direction(P,q,rng)
        except ValueError:
            for target in [.02,.04]:
                for method in ['frozen_latent','material_native']:
                    retractions.append(dict(index=i,target_RMS=target,method=method,actual_RMS=0.,material_goal_error=target,log_latent_prior=P.logprior(q),target_matched=False,failure='DEGENERATE_MATERIAL_TANGENT; no redraw',candidate_only=True,Maxwell_validated=False))
            continue
        for target in [.02,.04]:
            frozen,matched=frozen_latent_step(P,q,v,target);native,rec=material_retraction(P,q,v,w,target)
            for method,new,info in [('frozen_latent',frozen,dict(target_matched=matched)),('material_native',native,rec)]:
                x=P.decode(new)[0];base=P.decode(q)[0];desired=base+target*w
                row=dict(index=i,target_RMS=target,method=method,actual_RMS=P.material_rms(x,base),material_goal_error=P.material_rms(x,desired),log_latent_prior=P.logprior(new),structural_change=structural_change(P,q,new),candidate_only=True,Maxwell_validated=False);row.update(info)
                retractions.append(row);key=f'{i}_{target}_{method}';retraction_arrays[key+'_q']=new;retraction_arrays[key+'_chi']=x
                write_json(OUT/'results/RETRACTION_RESULTS.json',retractions)
                np.savez_compressed(OUT/'results/MATERIAL_RETRACTIONS.npz',**retraction_arrays)
    np.savez_compressed(OUT/'results/MATERIAL_RETRACTIONS.npz',**retraction_arrays)
    write_json(OUT/'results/RETRACTION_RESULTS.json',retractions)
    # Prior-reversible proposal exercise is a prior chain, never posterior output.
    q=groups['calibration'][0][0].copy();rng=np.random.default_rng(2026101813);prior_chain=[]
    for _ in range(128):q=P.local_prior_reversible(q,rng);prior_chain.append(q.copy())
    np.savez_compressed(OUT/'results/PRIOR_ONLY_PROPOSALS.npz',latent=np.array(prior_chain))
    # Confirm the zero allocation stops before loading any physical backend.
    target=MaxwellTarget();denied=None
    try:target.forward(np.r_[np.zeros(64),.5,.5],purpose='zero-allocation-enforcement')
    except BudgetDenied as e:denied=str(e)
    imported=[m for m in sys.modules if m=='a30' or m.startswith(('a30.','a24.','a25.'))]
    event_counts=collections.Counter(r['events']['0.5']['event'] for r in metadata if r['group']=='calibration')
    event_prior={k:dict(count=v,n=128,estimated_prior_mass=v/128,binomial_Wilson95=wilson(v,128),posterior_mass=None) for k,v in event_counts.items()}
    # Pair selection is a fixed adjacency matching, no image-quality selection.
    xx=groups['calibration'][1];pair_dist=[P.material_rms(xx[2*i],xx[2*i+1]) for i in range(64)]
    dx=(all_x-P.host).reshape(-1,576);physical_ranks=np.sum(abs(dx)>0,axis=1)
    metrics=dict(target_id=cfg['target_id'],frozen_code_commit=commit,material_samples=len(metadata),material_calibration=128,development=32,independent_material_check=32,independent_Maxwell_validation=False,ROI_cells=int(sum(P.mask)),latent_dimension=66,field_basis_rank=int(np.linalg.matrix_rank(P.A)),field_basis_orthogonality_offdiag=float(np.max(abs(P.A.T@P.A/576-np.diag(np.diag(P.A.T@P.A/576))))),material_prior_bounds=dict(Re_chi=[float(all_x.real.min()),float(all_x.real.max())],Im_chi=[float(all_x.imag.min()),float(all_x.imag.max())],fraction=[float(all_f.min()),float(all_f.max())]),physical_DeltaX_rank_range=[int(min(physical_ranks)),int(max(physical_ranks))],pair_material_RMS_quantiles=np.quantile(pair_dist,[0,.1,.5,.9,1]),material_derivative_checks=derivative,material_tangent_spectra=spectrum,prior_event_mass_estimates=event_prior,prior_empty_material_draws=int(sum(np.max(all_f,axis=1)==0)),prior_reversible_proposals=128,zero_physics_guard=dict(denied=denied,backend_modules_loaded=imported),material_retraction_summary={method:dict(n=sum(r['method']==method for r in retractions),matched=sum(r['method']==method and r['target_matched'] for r in retractions),mean_goal_error=float(np.mean([r['material_goal_error'] for r in retractions if r['method']==method]))) for method in ['frozen_latent','material_native']},noise_norms_first_two=np.linalg.norm(noise[:2].reshape(2,-1),axis=1)*np.sqrt(2)/sigma,noise_ball_tau=float(np.sqrt(chi2.ppf(.95,512))),noise_sigma=sigma,new_Maxwell_states=0,new_solver_RHS=0,new_physical_LU=0,NN_training=0,posterior_mass=None,observations='NOT_RUN',physics_only_anchors='NOT_RUN',A17_competition='NOT_RUN')
    write_json(OUT/'METRICS.json',metrics)
    receipt=dict(event='MATERIAL_CONSTRUCTION_COMPLETED',wall_seconds=time.perf_counter()-start,frozen_code_commit=commit,raw_material_samples=192,derivative_cases=6,retraction_endpoints=24,prior_only_proposal_steps=128,new_Maxwell_states=0,new_solver_RHS=0,new_physical_LU=0,GPU_training=0)
    with (OUT/'RESOURCE_LEDGER.jsonl').open('a') as f:f.write(json.dumps(receipt)+'\n')
    write_json(OUT/'MATERIAL_STAGE_CLOSED.json',dict(closed=True,scope='prior_material_and_model_construction_only',model_status='IMPLEMENTED_AND_MATERIAL_AUDITED',observation_status='NOT_RUN',anchor_status='NOT_RUN',physics_authorization=0))
    print(json.dumps(dict(receipt=receipt,material_prior_event_counts=dict(event_counts),ROI_cells=metrics['ROI_cells'],retraction=metrics['material_retraction_summary'],backend_modules=imported,noise_norms_first_two=metrics['noise_norms_first_two'].tolist()),ensure_ascii=False))

if __name__=='__main__':main()
