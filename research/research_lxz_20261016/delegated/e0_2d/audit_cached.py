"""LXZ E0: existing arrays only, no solver/adjoint/JVP imports or calls."""
from pathlib import Path
import csv,json,sys,time
import numpy as np
from scipy.stats import chi2

ROOT=Path(__file__).resolve().parents[4]
OUT=Path(__file__).resolve().parent
A30=ROOT.parents[1]/'A30/a30-maxwell-material-manifold'
RELEASE=A30/'inputs/a29_release'
PC=ROOT.parent/'b1-prior-constrained'
PCS=PC/'research/prior_constrained_20261015'
T=ROOT.parent/'b1-research-t/research/research_t_20261015'
sys.path.insert(0,str(ROOT/'src'))
from prior_constrained.prior import phase_audit,draw_train_row,verify_exact_endpoint

def read(path,keys):
    if 'locked' in str(path).lower():raise RuntimeError('LOCKED_FORBIDDEN')
    with np.load(path,allow_pickle=False)as z:return {k:z[k].copy()for k in keys}

def dump(path,x):path.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def csvwrite(path,rows):
    with path.open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)

def main():
    start=time.perf_counter();tau=float(np.sqrt(chi2.ppf(.95,512)))
    pooldir=RELEASE/'results/a29/pools';cases=sorted(p.stem for p in pooldir.glob('p*_a1_snr20.npz'))
    online=[];offline=[];sigma_map={};obs_map={};bg_map={};n3cache={}
    for case in cases:
        a=read(pooldir/f'{case}.npz',['n3','observation','sigma','background'])
        p=RELEASE/'inputs/a28/results/a28/exact_states'/f'{case}.npz'
        c=read(p,['chi','data'])
        if np.max(abs(c['chi'][0]-a['n3']))>1e-12:raise RuntimeError('CACHE_N3_MISMATCH')
        sig=float(a['sigma']);y=a['observation'];bg=a['background'];prediction=c['data'][0]
        sigma_map[case]=sig;obs_map[case]=y;bg_map[case]=bg;n3cache[case]=prediction
        residual=float(np.sqrt(2)*np.linalg.norm(prediction-y)/sig)
        online.append(dict(case_id=case,whitened_residual=residual,ell=.5*residual**2,tau=tau,
            noise_ball=residual<=tau,**phase_audit(a['n3'],bg),
            prior_exact=False,eligible_original_P0_anchor=False,source=str(p)))
    # Independent offline diagnostic: only known exposed validation clean data;
    # material truth/geometry are not needed and are not loaded here.
    valp=RELEASE/'inputs/a25/validation.npz'
    val=read(valp,['case_id','clean','data','sigma','background_id','amplitude','noise_db'])
    for row in online:
        case=row['case_id'];ii=int(np.flatnonzero(val['case_id']==case)[0]);clean=val['clean'][ii]
        sig=sigma_map[case];y=obs_map[case]
        if np.max(abs(y-val['data'][ii]))>1e-14 or abs(sig-float(val['sigma'][ii]))>1e-14:raise RuntimeError('VALIDATION_OBSERVATION_MISMATCH')
        distance=float(np.sqrt(2)*np.linalg.norm(n3cache[case]-clean)/sig)
        truthnoise=float(np.sqrt(2)*np.linalg.norm(clean-y)/sig)
        offline.append(dict(case_id=case,N3_to_true_clean_prediction_distance=distance,truth_noise_residual=truthnoise,
            truth_noise_inside95=truthnoise<=tau,N3_noise_compatible=row['noise_ball'],
            N3_prior_phase_compatible=row['phase_valid'],use='OFFLINE_DIAGNOSTIC_ONLY_NOT_ANCHOR_SELECTION'))

    # Score the public TRAIN material bank against development observations.
    # This can be an online prior-reservoir calculation; no development chi,
    # family, geometry, truth, or clean reference guides the choices.
    trainp=RELEASE/'inputs/a25/train.npz'
    tr=read(trainp,['parent_id','case_id','background_id','amplitude','noise_db','chi','data','clean','sigma','current'])
    reservoir=[]
    for case in cases:
        bgid='ring' if '_ring_'in case else'weak'
        ids=np.flatnonzero((tr['background_id']==bgid)&np.isclose(tr['amplitude'],1.2)&np.isnan(tr['noise_db']))
        norms=np.sqrt(2)*np.linalg.norm(tr['clean'][ids]-obs_map[case],axis=(1,2))/sigma_map[case]
        loc=int(np.argmin(norms));idx=int(ids[loc]);minnorm=float(norms[loc])
        reservoir.append(dict(case_id=case,public_prior_candidates=len(ids),best_TRAIN_parent=int(tr['parent_id'][idx]),
            whitened_residual=minnorm,ell=.5*minnorm**2,noise_ball=minnorm<=tau,
            original_P0_exact=True,selection='MIN_RESIDUAL_PUBLIC_TRAIN_ONLY_NO_DEV_TRUTH',
            source=str(trainp),source_index=idx))

    # Two reproducible known TRAIN calibration anchors, not inverse-discovered
    # answers to the development y and not hidden prior-stage witness arrays.
    calibration=[]
    for parent,bgid,fi,local in [(250033,'weak',2,0),(250034,'ring',2,1)]:
        ii=int(np.flatnonzero((tr['parent_id']==parent)&(tr['background_id']==bgid)&np.isclose(tr['amplitude'],1.2)&np.isclose(tr['noise_db'],20))[0])
        rr=draw_train_row(251008,fi,local,parent_id=parent)
        cert=verify_exact_endpoint(tr['chi'][ii],rr,background_id=bgid)
        residual=float(np.sqrt(2)*np.linalg.norm(tr['clean'][ii]-tr['data'][ii])/float(tr['sigma'][ii]))
        if cert['prior_status']!='P0_EXACT':raise RuntimeError('TRAIN_REPLAY_FAILED')
        if residual>tau:continue  # Do not modify tau or choose a favorable seed.
        name=f'PUBLIC_TRAIN_CALIBRATION_{parent}_{bgid}.npz'
        backgroundp=RELEASE/'inputs/a25/background'/f'{bgid}.npz'
        cached=read(backgroundp,['G','S','incident','points','volume'])
        field=cached['incident']+cached['G']@tr['current'][ii]
        constitutive=float(np.linalg.norm(tr['current'][ii]-tr['chi'][ii,:,None]*field)/np.linalg.norm(tr['current'][ii]))
        data_defect=float(np.linalg.norm(cached['S']@tr['current'][ii]-tr['clean'][ii])/np.linalg.norm(tr['clean'][ii]))
        np.savez_compressed(OUT/name,chi=tr['chi'][ii],prediction=tr['clean'][ii],observation=tr['data'][ii],
            sigma=tr['sigma'][ii],current=tr['current'][ii],field=field,points=cached['points'],volume=cached['volume'],
            background=bg_map[f'p250065_{bgid}_a1_snr20'])
        gen_path=OUT/name.replace('.npz','.json')
        dump(gen_path,dict(generator=rr,background_id=bgid,amplitude=1.2,prior_certificate=cert,
            role='PUBLIC_TRAIN_CALIBRATION_KNOWN_MATERIAL_NOT_INVERSE_DISCOVERY',
            source=str(trainp),source_index=ii,cached_equation_checks=dict(constitutive_relative=constitutive,prediction_relative=data_defect),
            target_measure=dict(identity=f'LXZ_P0_CONDITIONAL_SMOOTH_{parent}_{bgid}_A1_SNR20',
                parent_measure='Original A25 P0 conditioned on family=smooth, amplitude=1.2 and known background; not the original four-family mixture target',
                original_family_probability=.25,family_condition='smooth',
                parameter_prior='Original independent uniform TRAIN generator parameters, original support/rasterizer/normalization; inactive parameters marginalized, no field-uniform or geometry-radius-uniform substitution',
                likelihood='exp(-.5*||sqrt(2)/sigma*(F(chi)-y)||^2);512 real channels;gain=1;proper fixed sigma',
                posterior_mass_estimation=False),
            normalization='parent_material computes first two Gaussians in transformed coordinates then divides by actual raster maximum; chi=chi_b+1.2*p/max(p)*(1+.2i)/sqrt(1.04)',
            continuous_parameter_cautions=['Gaussian common weight scaling cancels under max normalization','Global coordinate transforms and Gaussian local coordinates can be redundant','Unused shell/rectangle/third-Gaussian parameters do not generate material changes; do not call them physical null directions'],
            no_new_physics=True,cached_matrix_products=dict(G_source_columns=8,S_source_columns=8)))
        calibration.append(dict(anchor_id=f'TRAIN_{parent}_{bgid}',source=str(trainp),source_index=ii,
            source_fields=['chi','clean','data','sigma','current'],copied_arrays=str(OUT/name),
            generator_json=str(gen_path),family='smooth',
            prior_certificate=cert,whitened_residual=residual,tau=tau,noise_ball_compatible=True,
            allowed_use='PUBLIC_SUPERVISED_TRAIN_CALIBRATION_ONLY',
            online_inverse_discovery=False,development_observation_compatible=False,
            solver_factorization_saved=False,exact_anchor_J_saved=False,
            leakage_warning='Known public TRAIN material/observation; do not call this an unknown-object reconstruction or use offline hidden E1 witnesses.'))

    audit=list(csv.DictReader((PCS/'delegated/cached_audit/E0_AUDIT.csv').open()))
    role_counts={}
    for role in sorted(set(r['role']for r in audit)):
        rs=[r for r in audit if r['role']==role]
        role_counts[role]=dict(records=len(rs),phase_pass=sum(r['phase']=='True'for r in rs),noise_pass=sum(r['noise_gate']=='True'for r in rs))
    g_inv=json.loads((PCS/'delegated/g_cache/inventory.json').read_text())
    geometry={};jassets=[]
    for bgid in ['weak','ring']:
        p=RELEASE/'inputs/a25/background'/f'{bgid}.npz'
        with np.load(p,allow_pickle=False)as z:
            # Shapes only; no derivative computation or multiplication.
            j=z['J'];chi_b=z['chi_b']
            jassets.append(dict(path=str(p),keys=['J','chi_b'],J_shape=list(j.shape),J_dtype=str(j.dtype),
                anchor='KNOWN_BACKGROUND_ONLY',original_P0_anchor=False,
                saved_full_J_finite=bool(np.isfinite(j).all()),usable_for_P0_weak_direction_claim=False,
                missing='No persisted exact Jacobian at an original-P0 observation-compatible development reconstruction.'))
    # Distinct paid PC endpoints only: inspect full secants already paid, never
    # compute a new directional response or new rasterized material.
    pc_arrays=[]
    for p in sorted((PCS/'offline_only').rglob('*.npz')):
        meta=p.with_suffix('.json')
        if not meta.exists():continue
        m=json.loads(meta.read_text());a=read(p,['chi','prediction'])
        pc_arrays.append(dict(path=str(p),metadata=m,chi=a['chi'],prediction=a['prediction']))
    duplicates=[];secants=[];paid_secant_rows=[]
    for i,a in enumerate(pc_arrays):
        for b in pc_arrays[i+1:]:
            ba=a['metadata'].get('background_id');bb=b['metadata'].get('background_id')
            if ba!=bb:continue
            dm=float(np.sqrt(np.mean(abs(a['chi']-b['chi'])**2)))
            df=float(np.linalg.norm(a['prediction']-b['prediction']))
            if dm<=1e-14:
                duplicates.append(dict(a=a['path'],b=b['path'],material_rms=dm,response_norm=df,
                    status='SAME_MATERIAL_REPEATED_CACHE_NOT_PHYSICAL_WEAK_DIRECTION'))
            elif df<=1e-14:
                secants.append(dict(a=a['path'],b=b['path'],material_rms=dm,response_norm=df,
                    status='DISTINCT_MATERIAL_IDENTICAL_RESPONSE_REQUIRES_REVIEW'))
            else:
                sigma_fixed=float(np.median([s for c,s in sigma_map.items() if f'_{ba}_'in c]))
                paid_secant_rows.append(dict(a=a['path'],b=b['path'],background=ba,material_rms=dm,
                    response_norm=df,d_pred=float(np.sqrt(2)*df/sigma_fixed),sigma_fixed=sigma_fixed,
                    scale='ONE_MEDIAN_HISTORICAL_SIGMA_PER_BACKGROUND',
                    scope='EXISTING_PAID_NONLINEAR_SECANT_NOT_JVP_NOT_CONTINUUM_CERTIFICATE'))

    csvwrite(OUT/'N3_ONLINE_CACHE_AUDIT.csv',online)
    csvwrite(OUT/'OFFLINE_N3_VS_CLEAN.csv',offline)
    csvwrite(OUT/'PUBLIC_TRAIN_BANK_VS_DEV.csv',reservoir)
    if duplicates:csvwrite(OUT/'PAID_CACHE_EXACT_DUPLICATES.csv',duplicates)
    if paid_secant_rows:csvwrite(OUT/'PAID_P0_NONLINEAR_SECANTS.csv',paid_secant_rows)
    result=dict(schema='LXZ_E0_2D_ASSETS_V1',new_physics=dict(states=0,RHS=0,LU=0,adjoints=0,JVPs=0,GPU=0),
        no_new_material_generation=True,no_samplers=True,no_new_P1=True,locked_reads=0,
        existing_material_decoder_replays=2,cached_matrix_products=dict(G_source_columns=16,S_source_columns=16),
        original_development_P0_compatible_anchors_found=sum(r['noise_ball']for r in reservoir),
        N3_noise_pass=sum(r['noise_ball']for r in online),N3_original_phase_pass=sum(r['phase_valid']for r in online),
        N3_prior_exact=0,N3_pass_cases=[r['case_id']for r in online if r['noise_ball']],
        public_TRAIN_calibration_anchors=calibration,old_candidate_records=role_counts,
        offline_N3_to_truth_prediction_distance_range=[min(r['N3_to_true_clean_prediction_distance']for r in offline),max(r['N3_to_true_clean_prediction_distance']for r in offline)],
        observed_truth_noise_outliers=[r['case_id']for r in offline if not r['truth_noise_inside95']],
        nearest_public_TRAIN_residual_range=[min(r['whitened_residual']for r in reservoir),max(r['whitened_residual']for r in reservoir)],
        background_J_assets=jassets,existing_P0_paid_endpoint_records=len(pc_arrays),
        repeated_same_material_cache_pairs=len(duplicates),distinct_material_zero_secant_pairs=secants,
        closest_nonduplicate_paid_P0_pair=min(paid_secant_rows,key=lambda r:r['d_pred']) if paid_secant_rows else None,
        raster_duplicate_conclusion='Zero raster secant means same discrete material, not a nontrivial physical weak/null direction. No local parameter derivative was newly evaluated.',
        P0_weak_direction_status='NOT_ESTABLISHED_FROM_EXISTING_PERSISTED_J_AT_A_QUALIFIED_DEV_ANCHOR',
        old_corrector_status='T/G saved raw and corrected endpoints fail original increment phase/range/peak; G median edit projection ~0.006687 for valid trajectories signals collapse, not proved topology retention.',
        usable_backend=dict(path=str(A30/'src/a30/physics.py'),warning='Do not call cached_anchor in E0: it performs 8 newly charged LU-audit RHS even with a saved N3 LU.',
            exact_J='receiver_transfer needs32complex adjoint RHS at current state; subsequent analytic material products do not make preparation free.'),
        offline_hidden_P0_E1_anchor_status='REUSABLE_ONLY_AS_EXPLICIT_OFFLINE_BENCHMARK_NOT_ONLINE_RECONSTRUCTION',
        online_truth_leakage=False,truth_usage='validation clean only in separate OFFLINE_N3_VS_CLEAN diagnostic; no chi/family/geometry used for dev selection',
        sources=[str(PCS/'PRIOR_SUPPORT_AUDIT.md'),str(PCS/'SOURCE_MANIFEST.json'),str(PCS/'E0_E1_E2_DECISION.md'),str(PCS/'delegated/cached_audit/SUMMARY.json'),str(PCS/'delegated/g_cache/inventory.json'),str(trainp),str(valp)],
        arithmetic_wall_seconds=time.perf_counter()-start)
    dump(OUT/'assets.json',result)
    print(json.dumps({k:result[k]for k in ['original_development_P0_compatible_anchors_found','N3_noise_pass','N3_original_phase_pass','N3_pass_cases','nearest_public_TRAIN_residual_range','arithmetic_wall_seconds']},indent=2))

if __name__=='__main__':main()
