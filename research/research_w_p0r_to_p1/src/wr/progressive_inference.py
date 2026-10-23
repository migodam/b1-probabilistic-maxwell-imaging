"""Controlled-family progressive inference on a paid-public-bank surrogate.

No function runs at import. No Book/FEM/Family.point or truth input is used.
main reads observation-only P1B data after the P1A reference gate. Every target
here is an explicitly approximate finite quadrature of W_FIXED_SUPPORT_JOINT_V1,
not a general unknown-wall posterior. Point profiling is ONLY an ML/MAP baseline.
"""
from __future__ import annotations
import json
import math
from pathlib import Path
import time
import numpy as np
from scipy.special import logsumexp
from .prob_reference import interpolate,marginal,posterior,cdf_coherence
from .sequential_particles import smc,make_proposal_library,exact_finite_reference

ROOT=Path(__file__).resolve().parents[2]
FAMILY_ID='W_FIXED_SUPPORT_JOINT_V1'
SIGMA=.01
FREQUENCIES=[.9975,.99875,1.,1.00125,1.0025]
CENTERS=np.array([[-.9,2.15],[-.8,2.15]])
MAX_ARRAY_BYTES=1_200_000_000
AMPLITUDE=(.99,1.01)
PHASE=(-.01,.01)


def _encode(value):
    if isinstance(value,Path):return str(value)
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)


def _write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,default=_encode,indent=2,allow_nan=False)+'\n')
    temp.replace(path)


def _load_rule(path):
    with np.load(path,allow_pickle=False) as archive:r={k:archive[k].copy() for k in archive.files}
    n=len(r['weights'])
    if r['F'].shape!=(n,5,4,12) or not np.iscomplexobj(r['F']) or not np.isfinite(r['F']).all():
        raise ValueError('Invalid complete public F bank:'+str(path))
    for k in ('branches','u','contrast','x'):
        if r[k].shape!=(n,) or not np.isfinite(r[k]).all():raise ValueError('Invalid rule coordinates:'+str(path))
    if not np.isfinite(r['weights']).all() or np.any(r['weights']<=0):raise ValueError('Invalid rule prior weights')
    if np.any((r['branches']!=0)&(r['branches']!=1)):raise ValueError('Invalid rule branches')
    return r


def material_grid(nu=97,nc=9):
    """Separable positive quadrature, C order(branch,iu,ic)."""
    zu,wu=np.polynomial.legendre.leggauss(nu);zc,wc=np.polynomial.legendre.leggauss(nc)
    u=.01*zu;c=3+.06*zc
    branch,uu,cc=np.meshgrid(np.arange(2),u,c,indexing='ij')
    weights=.5*np.broadcast_to((wu/2)[None,:,None]*(wc/2)[None,None,:],(2,nu,nc))
    return dict(branches=branch.ravel(),u=uu.ravel(),contrast=cc.ravel(),
        x=(CENTERS[branch.ravel(),0]+uu.ravel()),weights=weights.ravel(),
        u_axis=u,contrast_axis=c,u_weights=wu/2,contrast_weights=wc/2,shape=(2,nu,nc))


def gamma_rule(n):
    z,w=np.polynomial.legendre.leggauss(n);aa,pp=np.meshgrid(1+.01*z,.01*z,indexing='ij')
    return dict(amplitude=aa.ravel(),phase=pp.ravel(),weights=(w[:,None]*w[None,:]/4).ravel(),
        amplitude_axis=1+.01*z,phase_axis=.01*z,axis_weights=w/2,n=int(n))


def source_statistics(F,y):
    """One F per material node; sums frequency/receiver, preserves FOUR sources."""
    f=np.asarray(F);y=np.asarray(y)
    if f.ndim!=4 or f.shape[1:]!=(5,4,12) or y.shape!=(5,4,12) or not np.isfinite(f).all() or not np.isfinite(y).all():
        raise ValueError('Expected finite F(nmaterial,5,4,12),y(5,4,12)')
    fr,fi=f.real,f.imag;yr,yi=y.real,y.imag
    return dict(Fn=np.sum(fr*fr+fi*fi,axis=(1,3)),
        C=np.sum((fr*yr+fi*yi)+1j*(fi*yr-fr*yi),axis=(1,3)),
        Y=np.sum(yr*yr+yi*yi,axis=(0,2)),source_complex_count=60)


def _squares(stats,amplitude,phase,sources):
    """Residual norm square for all material rows at one gamma, on a prefix."""
    ids=np.asarray(sources,dtype=int)
    fn=np.sum(stats['Fn'][:,ids],axis=1);cross=np.sum(stats['C'][:,ids],axis=1);yn=float(np.sum(stats['Y'][ids]))
    square=amplitude**2*fn-2*amplitude*np.real(np.exp(1j*phase)*cross)+yn
    scale=64*np.finfo(float).eps*(amplitude**2*fn+2*abs(amplitude*cross)+yn)
    if np.any(square < -scale):raise FloatingPointError('Material residual square below roundoff scale')
    return np.maximum(square,0.)


def _likelihood_panel(stats,gamma,lo,hi):
    """(material chunk,gamma node,source), never repeats any240-entry F."""
    a=gamma['amplitude'][None,:,None];p=gamma['phase'][None,:,None]
    fn=stats['Fn'][lo:hi,None,:];cross=stats['C'][lo:hi,None,:];yn=stats['Y'][None,None,:]
    projected=np.cos(p)*cross.real-np.sin(p)*cross.imag
    square=a*a*fn-2*a*projected+yn
    scale=64*np.finfo(float).eps*(a*a*fn+2*abs(a*projected)+yn)
    if np.any(square < -scale):raise FloatingPointError('Block residual below numerical roundoff scale')
    constant=-2*stats['source_complex_count']*np.log(SIGMA*np.sqrt(2*np.pi))
    return constant-np.maximum(square,0.)/(2*SIGMA**2)


def nuisance_log_marginals(stats,n,chunk=32):
    gamma=gamma_rule(n);m=len(stats['Fn']);result=np.empty(m)
    for lo in range(0,m,chunk):
        hi=min(m,lo+chunk);blocks=_likelihood_panel(stats,gamma,lo,hi)
        total=np.sum(blocks,axis=-1);anchor=np.max(total,axis=1)
        result[lo:hi]=anchor+logsumexp((total-anchor[:,None])+np.log(gamma['weights'])[None,:],axis=1)
    return result


def select_gamma(stats,tolerance=2e-4):
    """Use21 if checked by41; else41 if checked by65; else65 UNQUALIFIED."""
    values={};history=[];previous=None
    for n in (21,41,65):
        current=nuisance_log_marginals(stats,n);values[n]=current
        error=None if previous is None else float(np.max(abs(current-values[previous])))
        history.append(dict(order=n,max_material_logZ_change=error))
        if error is not None and error<=tolerance:
            return previous,True,history,values[previous]
        previous=n
    return 65,False,history,values[65]


def joint_finite_target(stats,material,n,*,max_array_bytes=MAX_ARRAY_BYTES,chunk=32):
    gamma=gamma_rule(n);m=len(material['weights']);k=n*n;N=m*k
    # ll_blocks plus prior and short-lived material panels; no N*240 F array.
    estimated=32*N+8*N+64*chunk*k*4
    if estimated>max_array_bytes:raise MemoryError('Finite-target likelihood array allocation exceeds explicit cap')
    blocks=np.empty((N,4));view=blocks.reshape(m,k,4)
    for lo in range(0,m,chunk):
        hi=min(m,lo+chunk);view[lo:hi]=_likelihood_panel(stats,gamma,lo,hi)
    prior=(material['weights'][:,None]*gamma['weights'][None,:]).ravel()
    if abs(float(prior.sum())-1)>1e-12:raise ValueError('Joint branch/u/c/a/phi quadrature prior not normalized')
    return blocks,prior,gamma


def profile_shared_gamma(stats,sources=(0,1,2,3)):
    """Exact conditional box LS for each material, ONLY for point baselines.

    For each C, maximize Re(exp(i phi)C) over phi in[-.01,.01], then clip
    a=projection/Fn into[.99,1.01]. No nuisance marginal/evidence uses profiles.
    """
    ids=np.asarray(sources,dtype=int);fn=np.sum(stats['Fn'][:,ids],axis=1)
    cross=np.sum(stats['C'][:,ids],axis=1);yn=float(np.sum(stats['Y'][ids]))
    phases=np.stack((np.full(len(fn),PHASE[0]),np.full(len(fn),PHASE[1]),
        np.clip(-np.angle(cross),*PHASE)),axis=1)
    projected=np.real(np.exp(1j*phases)*cross[:,None])
    amplitudes=np.clip(np.divide(projected,fn[:,None],out=np.ones_like(projected),where=fn[:,None]>0),*AMPLITUDE)
    square=amplitudes**2*fn[:,None]-2*amplitudes*projected+yn
    scale=64*np.finfo(float).eps*(amplitudes**2*fn[:,None]+2*abs(amplitudes*projected)+yn)
    if np.any(square < -scale):raise FloatingPointError('Profile residual below numerical scale')
    square=np.maximum(square,0.);choice=np.argmin(square,axis=1);row=np.arange(len(fn))
    a=amplitudes[row,choice];p=phases[row,choice];phi=square[row,choice]/(2*SIGMA**2)
    zero=fn==0;a[zero]=1.;p[zero]=0.
    return dict(amplitude=a,phase=p,Phi=phi)


def point_baselines(stats,material):
    direct_cpu,direct_wall=time.process_time(),time.perf_counter()
    direct_profile=profile_shared_gamma(stats);j=int(np.argmin(direct_profile['Phi']))
    direct=dict(material_index=j,amplitude=float(direct_profile['amplitude'][j]),
        phase=float(direct_profile['phase'][j]),Phi=float(direct_profile['Phi'][j]),
        method='direct joint conditional-profile point',cpu_seconds=time.process_time()-direct_cpu,
        wall_seconds=time.perf_counter()-direct_wall,new_Maxwell_RHS=0)
    initial=np.flatnonzero((material['branches']==0)&(material['u']==0)&(material['contrast']==3.))
    if len(initial)!=1:raise ValueError('Odd material rules must contain the branch0,u0,c3 initial node')
    selected=int(initial[0]);a,p=1.,0.;history=[]
    progressive_cpu,progressive_wall=time.process_time(),time.perf_counter()
    for source in range(4):
        sources=tuple(range(source+1));cycles=[]
        for cycle in range(2):
            phi=_squares(stats,a,p,sources)/(2*SIGMA**2);selected=int(np.argmin(phi))
            prof=profile_shared_gamma(stats,sources)
            a=float(prof['amplitude'][selected]);p=float(prof['phase'][selected])
            cycles.append(dict(cycle=cycle+1,material_index=selected,amplitude=a,phase=p,
                prefix_Phi=float(prof['Phi'][selected])))
        history.append(dict(source_prefix=list(sources),cycles=cycles))
    progressive=dict(material_index=selected,amplitude=a,phase=p,
        Phi=float(_squares(stats,a,p,(0,1,2,3))[selected]/(2*SIGMA**2)),
        method='progressive source-prefix alternating point',stage_history=history,
        cpu_seconds=time.process_time()-progressive_cpu,wall_seconds=time.perf_counter()-progressive_wall,new_Maxwell_RHS=0)
    for row in (direct,progressive):
        idx=row['material_index'];row.update(branch=int(material['branches'][idx]),
            u=float(material['u'][idx]),contrast=float(material['contrast'][idx]),
            x=float(material['x'][idx]),target_y=2.15,
            validation_status='NEEDS_FULL_VALIDATION',response='branch-specific interpolated rawF approximation',
            objective='Phi only; shared-gamma profile is a point ML/MAP baseline, never posterior/evidence',
            fallback='NOT_RUN; paid full validation/fallback belongs to orchestrator')
    return {'direct joint point':direct,'progressive block point':progressive}


def continuous_reference(y,rules,protocol):
    """Fresh y: paid full5/full9 plus paid direct truncated integrals only."""
    orders=protocol['reference_orders'];cuts=protocol['cdf_checkpoints'];full={};nuisance={};qualified=True
    for order in orders:
        key=f'full_{order}';r=rules[key]
        if abs(float(r['weights'].sum())-1)>1e-12:raise ValueError('Full prior rule not normalized')
        answer,history,ok=marginal(r['F'],y);qualified &=ok;nuisance[key]=history
        full[order]=posterior(r,np.asarray(answer['logZ'],float))
    cdf={}
    for i,cut in enumerate(cuts):
        values=[]
        for order in orders:
            key=f'cdf_{i}_{order}';rule=rules[key]
            answer,history,ok=marginal(rule['F'],y);qualified &=ok;nuisance[key]=history
            ll=np.asarray(answer['logZ'],float);anchor=float(np.max(ll))
            lognumerator=anchor+float(logsumexp((ll-anchor)+np.log(rule['weights'])))
            values.append(float(np.exp(lognumerator-full[order]['logZ'])))
        cdf[str(cut)]={'order5':values[0],'order9':values[1],'error':abs(values[1]-values[0])}
    coherent=all(cdf_coherence([cdf[str(c)][f'order{n}'] for c in cuts]) for n in orders)
    event_error=abs(full[9]['branch1']-full[5]['branch1']);cdf_error=max(r['error'] for r in cdf.values())
    qualified &=coherent and event_error<=protocol['posterior_event_error_max'] and cdf_error<=protocol['posterior_cdf_error_max']
    return dict(qualified=bool(qualified),status='QUALIFIED_EMPIRICAL_QUADRATURE' if qualified else 'UNQUALIFIED',
        branch_mass=[1-full[9]['branch1'],full[9]['branch1']],x_mean=full[9]['x_mean'],x_sd=full[9]['x_sd'],
        logZ=full[9]['logZ'],CDF=cdf,CDF_coherent=coherent,event_order_error=event_error,CDF_order_max_error=cdf_error,
        nuisance_convergence=nuisance,meaning='continuous controlled-prior integral approximated by independent positive/truncated quadrature; not continuum exact')


def _weighted_material_summary(weights,material,cuts):
    w=np.asarray(weights,float)
    if w.shape!=(len(material['x']),) or not np.isfinite(w).all() or np.any(w<0) or abs(float(w.sum())-1)>1e-10:
        raise ValueError('Invalid normalized material posterior masses')
    x=material['x'];mean=float(w@x);branch=[float(w[material['branches']==b].sum()) for b in (0,1)]
    return dict(branch_mass=branch,x_mean=mean,x_sd=float(np.sqrt(w@((x-mean)**2))),
        CDF={str(c):float(w[x<=c].sum()) for c in cuts})


def finite_reference(blocks,prior,material,gamma,continuous,protocol,nuisance_ok):
    exact=exact_finite_reference(blocks,prior);w=exact['posterior_weights']
    wm=w.reshape(len(material['x']),gamma['n']**2).sum(axis=1)
    result=_weighted_material_summary(wm,material,protocol['cdf_checkpoints'])
    cdf_error=max(abs(v-continuous['CDF'][k]['order9']) for k,v in result['CDF'].items())
    branch_error=max(abs(a-b) for a,b in zip(result['branch_mass'],continuous['branch_mass']))
    coherent=cdf_coherence([result['CDF'][str(c)] for c in protocol['cdf_checkpoints']])
    qualified=bool(continuous['qualified'] and nuisance_ok and coherent and cdf_error<=protocol['posterior_cdf_error_max'] and branch_error<=protocol['posterior_event_error_max'])
    result.update(logZ=exact['log_evidence'],logZ_raw_prior_quadrature=exact['log_evidence_raw_prior_quadrature'],
        prior_quadrature_mass=exact['prior_input_mass'],CDF_coherent=coherent,CDF_error_vs_continuous=cdf_error,
        branch_error_vs_continuous=branch_error,nuisance_qualified=bool(nuisance_ok),qualified=qualified,
        status='QUALIFIED_CONTROLLED_FINITE_APPROXIMATION' if qualified else 'UNQUALIFIED',
        grid_shape=[*material['shape'],gamma['n'],gamma['n']],node_count=len(prior),
        physical_prior_changed=False,
        response_surrogate='branch-specific full9 interpolation; no new Maxwell material nodes',
        interpretation='exact finite enumeration of surrogate/quadrature target; finite atomic CDF checked at declared cuts against independently truncated continuous reference')
    return result,w


def _particle_summary(result,material,gamma,finite,continuous,cuts):
    idx=result['particle_indexes'];w=result['weights'];material_index=idx//(gamma['n']**2)
    x=material['x'][material_index];b=material['branches'][material_index];mean=float(w@x)
    masses=[float(w[b==branch].sum()) for branch in (0,1)]
    cdf={str(c):float(w[x<=c].sum()) for c in cuts}
    eligible=[mass>.05 for mass in finite['branch_mass']]
    missing=[bool(eligible[branch] and not np.any(b==branch)) for branch in (0,1)]
    return dict(branch_mass=masses,x_mean=mean,x_sd=float(np.sqrt(w@((x-mean)**2))),CDF=cdf,
        log_evidence_estimate=float(result['log_evidence_estimate']),
        log_evidence_error=float(result['log_evidence_estimate']-finite['logZ']),
        branch_mass_max_error_vs_finite=max(abs(a-c) for a,c in zip(masses,finite['branch_mass'])),
        CDF_max_error_vs_finite=max(abs(cdf[k]-finite['CDF'][k]) for k in cdf),
        CDF_max_error_vs_continuous=max(abs(cdf[k]-continuous['CDF'][k]['order9']) for k in cdf),
        x_mean_error_vs_finite=mean-finite['x_mean'],x_sd_error_vs_finite=float(np.sqrt(w@((x-mean)**2)))-finite['x_sd'],
        branch_miss_eligible=eligible,branch_missing=missing,branch_MISS=any(missing),
        branch_MISS_reference='same-observation exact finite enumeration; criterion mass>.05',
        missing_branches=sum(missing),eligible_branches=sum(eligible),
        surviving_initial_ancestors=result['surviving_initial_ancestors'],
        probability_claim='FINITE_TARGET_ESTIMATE_REFERENCE_QUALIFIED_METHOD_ACCURACY_NOT_ASSUMED' if finite['qualified'] else 'UNQUALIFIED_NO_PROBABILITY_SUCCESS_CLAIM',
        method_success_claim='NOT_ASSESSED; all finite-reference errors and missing-branch outcomes retained')


def _ensemble(rows,finite,continuous,cuts):
    means=np.array([r['x_mean'] for r in rows]);sds=np.array([r['x_sd'] for r in rows])
    masses=np.mean([r['branch_mass'] for r in rows],axis=0)
    cdf={str(c):float(np.mean([r['CDF'][str(c)] for r in rows])) for c in cuts}
    mean=float(means.mean());mixture_sd=float(np.sqrt(max(0.,np.mean(sds*sds+(means-mean)**2))))
    logz=float(logsumexp([r['log_evidence_estimate'] for r in rows])-np.log(len(rows)))
    eligible=sum(r['eligible_branches'] for r in rows)
    return dict(seed_count=len(rows),average_branch_mass=masses.tolist(),average_x_mean=mean,
        average_individual_x_sd=float(sds.mean()),ensemble_mixture_x_sd=mixture_sd,average_CDF=cdf,
        log_mean_evidence_estimate=logz,mean_log_evidence_estimate=float(np.mean([r['log_evidence_estimate'] for r in rows])),
        log_mean_evidence_error=logz-finite['logZ'],
        branch_mass_max_error_vs_finite=float(np.max(abs(masses-finite['branch_mass']))),
        CDF_max_error_vs_finite=max(abs(cdf[k]-finite['CDF'][k]) for k in cdf),
        CDF_max_error_vs_continuous=max(abs(cdf[k]-continuous['CDF'][k]['order9']) for k in cdf),
        branch_MISS_run_rate=float(np.mean([r['branch_MISS'] for r in rows])),
        branch_MISS_eligible_rate=(sum(r['missing_branches'] for r in rows)/eligible if eligible else None),
        average_surviving_initial_ancestors=float(np.mean([r['surviving_initial_ancestors'] for r in rows])),
        probability_claim='FINITE_TARGET_ESTIMATE_REFERENCE_QUALIFIED_METHOD_ACCURACY_NOT_ASSUMED' if finite['qualified'] else 'UNQUALIFIED_NO_PROBABILITY_SUCCESS_CLAIM',
        method_success_claim='NOT_ASSESSED; all finite-reference errors and missing-branch outcomes retained')


def _save_particles(path,result):
    keys=('particle_indexes','weights','ancestor_indexes','trajectory_indexes','trajectory_weights',
        'trajectory_stage','trajectory_sweep','mutation_accepted','mutation_global','resampling_indexes','resampling_performed','final_beta')
    np.savez_compressed(path,**{k:result[k] for k in keys},trajectory_kind=np.asarray(result['trajectory_kind'],dtype=str))


def main():
    gate_path=ROOT/'results/P1A/reference/GATE.json';observations=ROOT/'results/P1B/data/OBSERVATIONS.npz'
    if not gate_path.exists() or json.loads(gate_path.read_text()).get('probability_reference_qualified') is not True:
        raise RuntimeError('P1B inference blocked: P1A probability_reference_qualified must beTrue')
    if not observations.exists():raise RuntimeError('P1B observation-only file missing')
    protocol=json.loads((ROOT/'configs/P1_PROBABILITY_PROTOCOL.json').read_text())
    prior=json.loads((ROOT/'configs/W_CONTROLLED_PRIOR.json').read_text())
    if prior.get('family_id')!=FAMILY_ID:raise RuntimeError('Only the explicit controlled family is supported')
    expected_tx=np.c_[np.linspace(-1.4,1.4,4),np.full(4,-.3)]
    expected_rx=np.c_[np.linspace(-1.5,1.5,12),np.full(12,-.3)]
    if (protocol['reference_orders']!=[5,9] or prior['acquisition']['sigma_each_real_quadrature']!=SIGMA
        or prior['acquisition']['frequencies']!=FREQUENCIES
        or not np.array_equal(prior['acquisition']['tx'],expected_tx)
        or not np.array_equal(prior['acquisition']['rx'],expected_rx)
        or prior['wall_support']['wall_probability']!=[.5,.5] or prior['wall_support']['seed']!=2026102202
        or prior['wall_support']['cavities']!=[]
        or prior['target']['branch_centers']!=CENTERS.tolist() or prior['target']['radius']!=.12
        or prior['target']['conductivity']!=.02 or prior['target']['offset_x']['bounds']!=[-.01,.01]
        or prior['target']['contrast']['bounds']!=[2.94,3.06]
        or prior['nuisance']['amplitude']['bounds']!=list(AMPLITUDE)
        or prior['nuisance']['phase_radians']['bounds']!=list(PHASE)
        or any(item['distribution']!='uniform' for item in (prior['target']['offset_x'],prior['target']['contrast'],prior['nuisance']['amplitude'],prior['nuisance']['phase_radians']))):
        raise RuntimeError('Unexpected frozen reference/noise contract')
    settings=protocol['P1B'];seeds=settings['SMC_seeds'];P=settings['SMC_particles'];S=settings['rejuvenation_sweeps_each_bridge']
    if (len(seeds)!=4 or P!=2048 or S!=4 or settings['resampling_ESS_fraction']!=.5
        or settings['direct_temperatures']!=[0,.25,.5,.75,1]
        or settings['progressive_sources']!=[0,1,2,3]):
        raise RuntimeError('Frozen four-seed/2048-particle/four-sweep/ESS contract changed')
    directory=ROOT/'results/P1A/rules';needed=['full_5','full_9']+[f'cdf_{i}_{n}' for i in range(len(protocol['cdf_checkpoints'])) for n in (5,9)]
    rules={name:_load_rule(directory/(name+'.npz')) for name in needed}
    with np.load(observations,allow_pickle=False) as archive:ys=archive['y'].copy();parents=archive['parent_ids'].copy()
    if ys.ndim!=4 or ys.shape[1:]!=(5,4,12) or not np.iscomplexobj(ys) or not np.isfinite(ys).all() or parents.shape!=(len(ys),):
        raise ValueError('Malformed observation-only data')
    out=ROOT/'results/P1B/inference';out.mkdir(parents=True,exist_ok=True)
    contract=dict(family_id=FAMILY_ID,protocol=protocol,public_prior=prior,material_rule=[2,97,9],
        nuisance_check_orders=[21,41,65],nuisance_selected_rule='lower checked order;65unqualified if no agreement',
        global_proposal_probability=.25,max_array_bytes=MAX_ARRAY_BYTES,
        original_minimum_feature='UNRESOLVED',general_unknown_wall_online_claim=False,
        physical_prior_changed=False,finite_rule='fixed GL97u/GL9contrast quadrature approximation to the unchanged continuous controlled prior')
    contract_path=out/'MATH_CONTRACT.json'
    if contract_path.exists() and json.loads(contract_path.read_text())!=contract:raise RuntimeError('Existing inference math contract differs; refusing mixed restart')
    if not contract_path.exists():_write(contract_path,contract)
    setup_cpu,setup_wall=time.process_time(),time.perf_counter()
    material=material_grid();Fbase=interpolate(rules['full_9'],material['branches'],material['u'],material['contrast'])
    if Fbase.shape!=(1746,5,4,12) or not np.isfinite(Fbase).all():raise RuntimeError('Invalid public-bank surrogate')
    if (out/'PUBLIC_MATERIAL_GRID.npz').exists():
        with np.load(out/'PUBLIC_MATERIAL_GRID.npz',allow_pickle=False) as archived:
            if not np.array_equal(archived['F_approx'],Fbase):
                raise RuntimeError('Public bank/surrogate response changed; refusing stale inference restart')
            for key in ('branches','u','contrast','weights'):
                if not np.array_equal(archived[key],material[key]):
                    raise RuntimeError('Public material quadrature changed; refusing mixed restart')
    else:
        np.savez_compressed(out/'PUBLIC_MATERIAL_GRID.npz',F_approx=Fbase,**{k:v for k,v in material.items()})
    _write(out/'SETUP_COST.json',dict(cpu_seconds=time.process_time()-setup_cpu,wall_seconds=time.perf_counter()-setup_wall,
        new_Maxwell_RHS=0,public_bank_setup='paid P1A full9; setup/validation/fallback is not free deployment',
        response_label='branch-specific degree8x8 interpolation, an approximation'))
    rows=[];run_cpu,run_wall=time.process_time(),time.perf_counter()
    for observation,y in enumerate(ys):
        obsroot=out/f'observation_{observation:03d}';obsroot.mkdir(exist_ok=True)
        completion=obsroot/'COMPLETE.json'
        if completion.exists():
            row=json.loads(completion.read_text());attempt=obsroot/row['attempt']
            with np.load(attempt/'OBSERVATION.npz',allow_pickle=False) as archived:
                if not np.array_equal(archived['y'],y):raise RuntimeError('Restart observation differs from saved observation')
            rows.append(row);_write(out/'RESULTS.json',rows);continue
        attempt_name=f'attempt_{len(list(obsroot.glob("attempt_*"))):03d}'
        attempt=obsroot/attempt_name;attempt.mkdir(exist_ok=False)
        np.savez_compressed(attempt/'OBSERVATION.npz',y=y,parent=parents[observation])
        cpu_start,wall_start=time.process_time(),time.perf_counter()
        try:
            ref_cpu,ref_wall=time.process_time(),time.perf_counter()
            continuous=continuous_reference(y,rules,protocol)
            _write(attempt/'CONTINUOUS_REFERENCE.json',continuous)
            stats=source_statistics(Fbase,y)
            n,nuisance_ok,history,material_logz=select_gamma(stats,protocol['nuisance_log_evidence_difference_max'])
            blocks,node_prior,gamma=joint_finite_target(stats,material,n)
            finite,weights=finite_reference(blocks,node_prior,material,gamma,continuous,protocol,nuisance_ok)
            finite['nuisance_convergence']=history;finite['gamma_order']=n
            finite['material_rule_order']=[97,9]
            np.savez_compressed(attempt/'FINITE_REFERENCE.npz',posterior_weights=weights,
                material_log_marginal_likelihood=material_logz,material_branches=material['branches'],material_u=material['u'],
                material_contrast=material['contrast'],material_x=material['x'],material_prior_weights=material['weights'],
                amplitude_axis=gamma['amplitude_axis'],phase_axis=gamma['phase_axis'],gamma_axis_weights=gamma['axis_weights'],
                grid_shape=np.asarray(finite['grid_shape']),logZ=finite['logZ'])
            _write(attempt/'FINITE_REFERENCE.json',finite)
            reference_cost=dict(cpu_seconds=time.process_time()-ref_cpu,wall_seconds=time.perf_counter()-ref_wall,
                new_Maxwell_RHS=0,teacher_cost='paid bank plus pure inference reference integration; not free method setup')
            point_cpu,point_wall=time.process_time(),time.perf_counter();points=point_baselines(stats,material)
            for name,point in points.items():
                label='direct_joint_point' if name.startswith('direct') else 'progressive_block_point'
                idx=point['material_index'];np.savez_compressed(attempt/(label+'.npz'),
                    params=[point[k] for k in ('branch','u','contrast','amplitude','phase')],F_approx=Fbase[idx],Phi=point['Phi'],
                    parameter_order=np.asarray(['branch','u','contrast','amplitude','phase']))
                _write(attempt/(label+'.json'),point)
            point_cost=dict(cpu_seconds=time.process_time()-point_cpu,wall_seconds=time.perf_counter()-point_wall,new_Maxwell_RHS=0)
            runs={'direct':[],'progressive':[]};proposal_costs=[]
            for seed in seeds:
                lib_cpu,lib_wall=time.process_time(),time.perf_counter()
                library=make_proposal_library(node_prior,finite['grid_shape'],rngseed=seed,nparticles=P,sweeps=S,max_array_bytes=MAX_ARRAY_BYTES)
                proposal_costs.append(dict(seed=seed,cpu_seconds=time.process_time()-lib_cpu,wall_seconds=time.perf_counter()-lib_wall))
                np.savez_compressed(attempt/f'COMMON_PROPOSALS_seed{seed}.npz',
                    **{k:library[k] for k in ('initial_nodes','global_nodes','global_uniforms','coordinates','signs','acceptance_uniforms','resampling_uniforms')},
                    grid_shape=np.asarray(finite['grid_shape']),prior='separable factors in FINITE_REFERENCE.npz')
                for schedule in ('direct','progressive'):
                    sample_cpu,sample_wall=time.process_time(),time.perf_counter()
                    result=smc(blocks,node_prior,finite['grid_shape'],schedule=schedule,nparticles=P,sweeps=S,
                        global_prob=.25,proposal_library=library,max_array_bytes=MAX_ARRAY_BYTES)
                    summary=_particle_summary(result,material,gamma,finite,continuous,protocol['cdf_checkpoints'])
                    summary.update(seed=seed,schedule=schedule,cpu_seconds=time.process_time()-sample_cpu,
                        wall_seconds=time.perf_counter()-sample_wall,new_Maxwell_RHS=0,
                        trajectory_frames=len(result['trajectory_kind']),estimated_array_bytes=result['estimated_working_array_bytes'])
                    _save_particles(attempt/f'{schedule}_seed{seed}_TRAJECTORY.npz',result)
                    _write(attempt/f'{schedule}_seed{seed}_STAGES.json',dict(summary=summary,stages=result['stages'],
                        final_beta=result['final_beta'],trajectory_interpretation='includes all rejects/repeated states/resampling; never pool different beta frames as finalposterior'))
                    runs[schedule].append(summary)
                    _write(attempt/'SMC_PROGRESS.json',runs)
                    del result
                del library
            ensembles={name:_ensemble(records,finite,continuous,protocol['cdf_checkpoints']) for name,records in runs.items()}
            row=dict(observation=observation,parent=int(parents[observation]),attempt=attempt_name,
                family_id=FAMILY_ID,scope='new restricted public fixed-support family, not general unknown-wall online',
                original_minimum_feature='UNRESOLVED',continuous_ref=continuous,finite_ref=finite,
                methods=points,SMC_runs=runs,ensemble=ensembles,reference_cost=reference_cost,
                point_cost=point_cost,shared_proposal_library_cost=proposal_costs,
                cpu_seconds=time.process_time()-cpu_start,wall_seconds=time.perf_counter()-wall_start,
                new_Maxwell_RHS=0,point_physics_validation='NEEDS_FULL_VALIDATION',
                numerical_probability_status='QUALIFIED_CONTROLLED_APPROXIMATION' if finite['qualified'] else 'UNQUALIFIED',
                outside_family_mass=None)
            _write(attempt/'RESULTS.json',row);_write(completion,row)
            rows.append(row);_write(out/'RESULTS.json',rows)
            print(json.dumps(dict(observation=observation,finite_qualified=finite['qualified'],gamma_order=n,
                direct_branch_error=ensembles['direct']['branch_mass_max_error_vs_finite'],
                progressive_branch_error=ensembles['progressive']['branch_mass_max_error_vs_finite'])),flush=True)
            del blocks,node_prior,weights
        except BaseException as exc:
            _write(attempt/'FAILURE.json',dict(observation=observation,error_type=type(exc).__name__,error=str(exc),
                cpu_seconds=time.process_time()-cpu_start,wall_seconds=time.perf_counter()-wall_start,
                no_physical_recovery_attempted=True))
            raise
    _write(out/'COMPLETE.json',dict(observations=len(rows),independent_parents=len(np.unique(parents)),
        all_finite_references_qualified=all(r['finite_ref']['qualified'] for r in rows),
        points='NEEDS_FULL_VALIDATION; not validated by this algebra/surrogate process',
        original_minimum_feature='UNRESOLVED',outside_family_mass=None,
        cpu_seconds=time.process_time()-run_cpu,wall_seconds=time.perf_counter()-run_wall,
        new_Maxwell_RHS=0,Book_events=0,CPU_billing='launcher alone bills process CPU once'))
    return rows


if __name__=='__main__':main()
