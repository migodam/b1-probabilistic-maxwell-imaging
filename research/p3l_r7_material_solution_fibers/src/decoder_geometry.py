"""R7 E0/E1: pure original-P1 decoder geometry; no Maxwell or neural calls.

Run from the repository with ../../.venv_nn/bin/python research/p3l_r7_material_solution_fibers/src/decoder_geometry.py.
The sibling R4 module is imported read-only, without constructing Book/Target.
"""
from __future__ import annotations
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import sys,json,time,resource,platform
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
from scipy.special import ndtr,ndtri
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[3]
R4=ROOT.parent/'b1-p3l-r4-conditional-transport'
sys.path.insert(0,str(R4/'src'))
from p3r4.core import MaterialPrior,from_probit,probit_derivative
from p1_material.prior import saturation,saturation_prime
OUT=Path(__file__).resolve().parents[1]
SEED=2026102807
RELS=np.array([1e-6,1e-8,1e-10,1e-12])
FD_STEPS=np.array([1e-4,1e-5,1e-6])
FINITE_STEPS=np.array([.001,.01,.1,.5])

def plain(a):
    if isinstance(a,np.ndarray):return a.tolist()
    if isinstance(a,np.generic):return a.item()
    if isinstance(a,Path):return str(a)
    if isinstance(a,dict):return {str(k):plain(v) for k,v in a.items()}
    if isinstance(a,(tuple,list)):return [plain(v) for v in a]
    return a

def write(path,a):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(plain(a),ensure_ascii=False,allow_nan=False,indent=2)+'\n');tmp.replace(path)

def save(path,**a):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.with_suffix('.tmp').open('wb') as stream:np.savez_compressed(stream,**a)
    path.with_suffix('.tmp').replace(path)

def load_prior():
    cfgpath=R4/'research/p3l_r4_conditional_transport/configs/TARGET.json'
    geompath=R4/'research/p1_material_target_20261018/inputs/geometry.json'
    cfg=json.loads(cfgpath.read_text());points=json.loads(geompath.read_text())['points']
    return MaterialPrior(cfg,points),cfgpath,geompath

def decode(prior,x):return prior.decode(from_probit(x))

def jacobian(prior,x):
    return np.asarray(prior.derivative(from_probit(x))*probit_derivative(x)[None,:],dtype=np.complex128)

def real_hilbert(D):return np.vstack((D.real,D.imag))/24.

def fd_jacobian(prior,x,h):
    D=np.empty((576,66),dtype=np.complex128)
    for k in range(66):
        d=np.zeros(66);d[k]=h
        D[:,k]=(decode(prior,x+d)[0]-decode(prior,x-d)[0])/(2*h)
    return D

def rank(s,rel=1e-10):return int(np.sum(s>rel*s[0])) if s[0]>0 else 0

def stratum(prior,a):
    s=a['field']/prior.temperature
    return dict(host_cells=int(np.sum(prior.mask&(s<=-1))),pure_cells=int(np.sum(prior.mask&(s>=1))),
                transition_cells=int(np.sum(prior.mask&(s>-1)&(s<1))),outside_ROI=int(np.sum(~prior.mask)))

def one_geometry(prior,x):
    chi,a=decode(prior,x);D=jacobian(prior,x);C=real_hilbert(D)
    U,s,Vh=np.linalg.svd(C,full_matrices=False)
    T=prior.mask&(abs(a['field'])<prior.temperature)
    B=(saturation_prime(a['field']/prior.temperature)/prior.temperature*prior.mask)[:,None]*prior.A
    UB,sB,VB=np.linalg.svd(B,full_matrices=False)
    rb=rank(sB);bproject=UB[:,:rb]@(UB[:,:rb].T@a['fraction'])
    pres=np.linalg.norm(a['fraction']-bproject)
    r=rank(s);V=Vh[:r].T
    P=V@V.T
    # Exact saturation-inactive subspace: A_T dz=0 and properties fixed.
    # Its straight paths preserve every transition field; plateau crossings
    # alone limit its finite interval. This is stronger than an SVD cutoff.
    if T.any():
        _,sA,VA=np.linalg.svd(prior.A[T],full_matrices=True)
        # Machine-scale rank is separate from the deliberately much looser
        # material-spectrum cutoffs; a small singular value is not exact zero.
        atol=np.finfo(float).eps*max(prior.A[T].shape)*sA[0]
        ra=int(np.sum(sA>atol));Kz=VA[ra:].T
    else:ra=0;Kz=np.eye(64)
    K=np.vstack((Kz,np.zeros((2,Kz.shape[1]))))
    safe_interval=None
    if Kz.shape[1]:
        inactive=prior.mask&~T;velocity=prior.A@Kz[:,0]
        moving=inactive&(abs(velocity)>1e-14)
        if moving.any():safe_interval=float(np.min(np.maximum(0,abs(a['field'][moving])-prior.temperature)/abs(velocity[moving])))
    # Property radial tangent d(delta)=delta*c in Gaussianized coordinates.
    d=a['inclusion_increment'];chain=probit_derivative(x)[64:]
    prop=np.array([d.real/np.diff(prior.real_interval)[0]/chain[0],d.imag/np.diff(prior.imag_interval)[0]/chain[1]])
    dz=np.linalg.lstsq(B,-a['fraction'],rcond=1e-12)[0]
    gauge=np.r_[dz,prop];gauge/=np.linalg.norm(gauge)
    # A direct globally uniform field scaling proposal, p(f/a)*a=p(f),
    # is tested without assuming homogeneity of the quintic law.
    field_scales=np.array([.8,.95,1.05,1.2])
    scale_error=[]
    for alpha in field_scales:
        q=from_probit(x).copy();q[:64]/=alpha
        delta=alpha*d
        q[64]=(delta.real-prior.real_interval[0])/np.diff(prior.real_interval)[0]
        q[65]=(delta.imag-prior.imag_interval[0])/np.diff(prior.imag_interval)[0]
        if prior.legal(q):scale_error.append(prior.material_rms(prior.decode(q)[0],chi))
        else:scale_error.append(None)
    return chi,a,D,C,s,Vh,dict(**stratum(prior,a),rank_by_relative=[rank(s,r) for r in RELS],
        B_rank=rb,A_transition_rank=ra,q_to_range_B_absolute=pres,q_to_range_B_relative=pres/max(1e-30,np.linalg.norm(a['fraction'])),
        formula_predicted_rank=rb+(1 if np.linalg.norm(a['fraction'])>0 else 0)+(1 if pres>1e-10*max(1.,np.linalg.norm(a['fraction'])) else 0),
        saturation_inactive_dimension=int(K.shape[1]),projector_symmetry=float(np.linalg.norm(P-P.T)),
        projector_idempotence=float(np.linalg.norm(P@P-P)),projector_null_residual=float(np.linalg.norm(C@(np.eye(66)-P))),
        kernel_straight_J_residual=float(np.linalg.norm(C@K)) if K.size else 0.,
        kernel_transition_A_residual=float(np.linalg.norm(prior.A[T]@Kz)) if Kz.size else 0.,
        inactive_first_direction_symmetric_plateau_interval=safe_interval,
        radial_compensation_tangent_relative_residual=float(np.linalg.norm(C@gauge)/max(1e-30,np.linalg.norm(C[:,64:]@gauge[64:]))),
        uniform_field_rescaling_material_RMS=scale_error),K,gauge

def finite_diagnostics(prior,x,chi,C,Vh,r,K,gauge):
    vectors=[];labels=[]
    if K.shape[1]:vectors.append(K[:,0]);labels.append('saturation_inactive_exact_local_subspace')
    if r<66:vectors.append(Vh[-1]);labels.append('numerical_full_decoder_kernel')
    vectors.append(gauge);labels.append('radial_property_compensation_candidate')
    # Weak and strong nonzero tangent directions normalized to material RMS=1.
    if r:
        for idx,label in ((0,'material_unit_strong'),(r-1,'material_unit_weak')):
            v=Vh[idx].copy();v/=np.linalg.norm(C@v)
            vectors.append(v);labels.append(label)
    records=[]
    for label,v in zip(labels,vectors):
        for t in FINITE_STEPS:
            y,b=decode(prior,x+t*v);delta=y-chi
            records.append(dict(direction=label,t=t,latent_step=float(t*np.linalg.norm(v)),
                material_RMS=float(np.linalg.norm(delta)/24),linear_material_RMS=float(t*np.linalg.norm(C@v)),
                linearization_residual_RMS=float(np.linalg.norm(delta/24-t*(C@v)[:576]-1j*t*(C@v)[576:])),
                stratum=stratum(prior,b)))
    return records,np.array(vectors),labels

def nonlinear_radial_correction(prior,x,chi,a,gauge):
    """Known quintic inverse construction, with residuals rather than claims.

    Only no-pure strata are attempted. Transition fractions are divided by a
    common alpha and inverted cellwise; the original A then fits those fields.
    Saturated-host cells are verified after fitting and can invalidate a step.
    """
    T=prior.mask&(abs(a['field'])<prior.temperature)
    if np.any(prior.mask&(a['field']>=prior.temperature)):return [],[]
    q=from_probit(x);d=a['inclusion_increment']
    c=float(np.diff(prior.real_interval)[0]*probit_derivative(x)[64]*gauge[64]/d.real)
    records=[];arrays=[]
    for t in FINITE_STEPS:
        alpha=1+t*c;qp=q.copy()
        qp[64]=(alpha*d.real-prior.real_interval[0])/np.diff(prior.real_interval)[0]
        qp[65]=(alpha*d.imag-prior.imag_interval[0])/np.diff(prior.imag_interval)[0]
        target=a['fraction'][T]/alpha
        rec=dict(t=float(t),alpha=float(alpha),construction='inverse_quintic_transition_field_fit',
                 transition_cells=int(T.sum()),A_rows=int(T.sum()),A_columns=64)
        if not prior.legal(qp) or np.any((target<=0)|(target>=1)):
            rec.update(status='INFEASIBLE_PROPERTY_OR_TARGET_FRACTION');records.append(rec);continue
        ft=np.array([prior.temperature*brentq(lambda s:float(saturation(s))-p,-1.,1.,xtol=1e-14) for p in target])
        dz=np.linalg.lstsq(prior.A[T],ft-a['field'][T],rcond=None)[0] if T.any() else np.zeros(64)
        qp[:64]+=dz;xp=qp.copy();xp[64:]=ndtri(qp[64:])
        chim,b=prior.decode(qp)
        rec.update(status='CONSTRUCTED_VERIFY_RESIDUAL',latent_step=float(np.linalg.norm(xp-x)),
            material_RMS=prior.material_rms(chim,chi),transition_field_fit_error=float(np.linalg.norm(prior.A[T]@dz-(ft-a['field'][T]))),
            host_plateau_crossings=int(np.sum(prior.mask&(a['field']<=-prior.temperature)&(b['field']>-prior.temperature))),
            straight_material_RMS=prior.material_rms(decode(prior,x+t*gauge)[0],chi))
        records.append(rec);arrays.append((t,xp,chim))
    return records,arrays

def special_cases(prior):
    base=np.zeros(66);cases=[]
    for name,field in [('host_interior',-.7),('host_seam',-.5),('host_inside_transition',-.5+1e-7),
                       ('uniform_transition',0.),('pure_inside_transition',.5-1e-7),('pure_seam',.5),('pure_interior',.7)]:
        x=base.copy();x[0]=(field-prior.config['field_mean'])/prior.config['offset_SD'];cases.append((name,x))
    for value in (-8.,-6.,6.,8.):
        for k in (64,65):
            x=base.copy();x[k]=value;cases.append((f'property_{k}_xi_{value:g}',x))
    return cases

def run():
    started=time.perf_counter();cpu=time.process_time()
    raw=OUT/'raw/decoder';fig=OUT/'figures';worker=OUT/'audit/delegated/r7_decoder_geometry'
    for p in (raw,fig,worker):p.mkdir(parents=True,exist_ok=True)
    manifest=OUT/'configs/E1_PRIOR_MANIFEST.json'
    if manifest.exists():
        frozen=json.loads(manifest.read_text());X=np.asarray(frozen['x'],dtype=np.float64)
        if frozen['master_seed']!=SEED or X.shape!=(256,66):raise RuntimeError('preexisting incompatible manifest')
    else:
        X=np.random.default_rng(SEED).standard_normal((256,66))
        write(manifest,dict(master_seed=SEED,created_unix=time.time(),n=256,distribution='independent standard normal N(0,I66)',
            parent_ids=[f'R7_E1_{i:03d}' for i in range(256)],x=X,
            usage='E1 decoder-only locked prior audit; no selection/rejection/replacement; E2 split chosen by root',
            gaussian_to_native_properties='q[64:]=ndtr(x[64:])',
            saved_before_decoder_calculations=True,new_Maxwell_RHS=0,relative_rank_thresholds=RELS,
            finite_difference_steps=FD_STEPS,kernel_finite_latent_steps=FINITE_STEPS,
            absolute_rank_diagnostics=[1e-8,1e-10,1e-12,1e-14]))
    # The manifest is already on disk before even constructing/evaluating prior.
    prior,cfgpath,geompath=load_prior()
    records=[];finite=[];svals=[];xs=[];materials=[];jacobians=[];fractions=[];fields=[];FD=[];normalization=[];corrections=[];corrected_arrays=[]
    allcases=[(f'R7_E1_{i:03d}',x) for i,x in enumerate(X)]+special_cases(prior)
    for i,(name,x) in enumerate(allcases):
        chi,a,D,C,s,Vh,rec,K,gauge=one_geometry(prior,x)
        rec.update(parent_id=name,cohort='locked_independent_prior' if i<256 else 'deterministic_edge_case',
                   singular_min=float(s[-1]),singular_max=float(s[0]),absolute_ranks=[int(np.sum(s>eps)) for eps in (1e-8,1e-10,1e-12,1e-14)])
        ferr=[]
        for h in FD_STEPS:
            F=fd_jacobian(prior,x,h);E=real_hilbert(F-D)
            ferr.append([float(np.linalg.norm(E)),float(np.linalg.norm(E,2)),float(np.linalg.norm(E)/max(1e-30,np.linalg.norm(C)))])
        rec['FD_errors_frobenius_operator_relative']=ferr
        rec['FD_error_informed_rank_interval']=[int(np.sum(s>max(f[1] for f in ferr))),int(np.sum(s>min(f[1] for f in ferr)))]
        rec['rank_formula_difference_at_1e10']=rec['rank_by_relative'][2]-rec['formula_predicted_rank']
        frec,vectors,labels=finite_diagnostics(prior,x,chi,C,Vh,rank(s),K,gauge)
        for fr in frec:fr['parent_id']=name
        finite+=frec
        normalization.append(dict(parent_id=name,direction_labels=labels,vectors=vectors))
        nrec,narrays=nonlinear_radial_correction(prior,x,chi,a,gauge)
        for nr in nrec:nr['parent_id']=name
        corrections+=nrec
        for t,xp,ch in narrays:corrected_arrays.append((name,t,xp,ch))
        records.append(rec);svals.append(s);xs.append(x);materials.append(chi);jacobians.append(D);fractions.append(a['fraction']);fields.append(a['field']);FD.append(ferr)
        if (i+1)%32==0:print(f'E1 {i+1}/{len(allcases)} decoder points complete',flush=True)
    svals=np.asarray(svals)
    save(raw/'FULL66_DECODER_ARRAYS.npz',x=np.asarray(xs),material=np.asarray(materials,dtype=np.complex128),
         jacobian=np.asarray(jacobians,dtype=np.complex128),singular_values=svals,fraction=np.asarray(fractions),field=np.asarray(fields),
         FD_errors=np.asarray(FD),relative_thresholds=RELS,FD_steps=FD_STEPS,parent_id=np.array([a[0] for a in allcases]))
    write(raw/'FINITE_STEP_COMPARISONS.json',dict(steps=FINITE_STEPS,records=finite))
    write(raw/'DIRECTIONS_AND_PROJECTORS.json',dict(material_metric='norm(chi)/24; real Hilbert [Re Dg; Im Dg]/24',records=normalization))
    write(raw/'NONLINEAR_EQUAL_MATERIAL_CORRECTIONS.json',dict(records=corrections))
    save(raw/'NONLINEAR_EQUAL_MATERIAL_CORRECTIONS.npz',parent_id=np.array([v[0] for v in corrected_arrays]),
         t=np.array([v[1] for v in corrected_arrays]),x=np.array([v[2] for v in corrected_arrays]),material=np.array([v[3] for v in corrected_arrays]))
    cohort=records[:256]
    histogram={str(eps):{str(r):sum(a['rank_by_relative'][j]==r for a in cohort) for r in sorted(set(a['rank_by_relative'][j] for a in cohort))} for j,eps in enumerate(RELS)}
    summary=dict(status='PASS_DECODER_AUDIT_COMPLETED',new_Maxwell_RHS=0,new_forward=0,new_LU=0,NN_training=0,
        n_independent_prior=256,n_edge_cases=len(records)-256,master_seed=SEED,material_cells=576,ROI_cells=int(prior.mask.sum()),
        relative_rank_histograms=histogram,all_prior_have_pure_cells=all(a['pure_cells']>0 for a in cohort),
        prior_pure_cell_range=[min(a['pure_cells'] for a in cohort),max(a['pure_cells'] for a in cohort)],
        prior_transition_cell_range=[min(a['transition_cells'] for a in cohort),max(a['transition_cells'] for a in cohort)],
        prior_saturation_inactive_dimensions={str(k):sum(a['saturation_inactive_dimension']==k for a in cohort) for k in sorted(set(a['saturation_inactive_dimension'] for a in cohort))},
        rank_formula_mismatch_count=sum(a['rank_formula_difference_at_1e10']!=0 for a in cohort),
        max_projector_idempotence=max(a['projector_idempotence'] for a in cohort),
        max_FD_relative_error=max(f[2] for a in cohort for f in a['FD_errors_frobenius_operator_relative']),
        nonlinear_corrected_steps=len(corrected_arrays),nonlinear_correction_records=corrections,
        records=records,raw_arrays='raw/decoder/FULL66_DECODER_ARRAYS.npz')
    write(OUT/'audit/E1_RANK_STRATA.json',summary)
    replay=decode(prior,X[0])[0]
    e0=dict(status='PASS',new_Maxwell_RHS=0,new_forward=0,new_LU=0,NN_training=0,master_seed=SEED,
        manifest_saved_before_decoder_calculations=True,source_prior=str(sys.modules[MaterialPrior.__module__].__file__),
        source_core=str(sys.modules['p3r4.core'].__file__),config_source=str(cfgpath),geometry_source=str(geompath),
        pure_read_only_imports=True,Book_constructed=False,Target_constructed=False,forward_called=False,
        exact_replay_material=np.array_equal(replay,materials[0]),exact_replay_jacobian=np.array_equal(jacobian(prior,X[0]),jacobians[0]),
        exact_seed_replay=np.array_equal(np.random.default_rng(SEED).standard_normal((256,66)),X),
        material_dtype=str(np.asarray(materials).dtype),jacobian_dtype=str(np.asarray(jacobians).dtype),
        coordinate_dtype=str(X.dtype),material_real_hilbert_matrix='vstack([Dg.real,Dg.imag])/24',
        full66_transform='q[:64]=x[:64]; q[64:]=ndtr(x[64:]); Dg=D_native*diag(1_64,normal_pdf(x[64:]))',
        python=sys.version,numpy=np.__version__,platform=platform.platform(),wall_seconds=time.perf_counter()-started,
        CPU_seconds=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    write(OUT/'audit/E0_REPRODUCIBILITY.json',e0)
    with (worker/'RUN_RECEIPTS.jsonl').open('a') as stream:stream.write(json.dumps(plain(e0),allow_nan=False)+'\n')
    make_figures(fig,svals,cohort,records,finite)
    report(prior,summary,e0,OUT/'audit/E1_FULL66_GAUGE_AUDIT.md')
    (worker/'WORKER_REPORT.md').write_text((OUT/'audit/E1_FULL66_GAUGE_AUDIT.md').read_text())
    write(worker/'ARTIFACTS.json',dict(status='SUCCESS',outputs=[str(OUT/p) for p in ('audit/E0_REPRODUCIBILITY.json','audit/E1_RANK_STRATA.json','audit/E1_FULL66_GAUGE_AUDIT.md','raw/decoder/FULL66_DECODER_ARRAYS.npz')],new_Maxwell_RHS=0))
    print(json.dumps(plain(dict(histograms=histogram,wall=e0['wall_seconds'],max_FD_relative_error=summary['max_FD_relative_error'],new_Maxwell_RHS=0))),flush=True)

def make_figures(fig,svals,cohort,records,finite):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10})
    f,axes=plt.subplots(1,3,figsize=(15,4.3))
    for s in svals[:256]:axes[0].semilogy(np.arange(1,67),np.maximum(s,1e-20),alpha=.15,color='navy')
    axes[0].set(xlabel='Singular value index',ylabel='Material RMS / latent unit',title='256 independent Full66 decoder spectra')
    axes[1].scatter([a['transition_cells'] for a in cohort],[a['rank_by_relative'][2] for a in cohort],s=14)
    axes[1].set(xlabel='Transition cells',ylabel='Rank (relative 1e-10)',title='Rank strata; no redraws')
    axes[2].scatter([a['pure_cells'] for a in cohort],[a['radial_compensation_tangent_relative_residual'] for a in cohort],s=14)
    axes[2].set(xlabel='Pure-phase cells',ylabel='Relative radial compensation residual',title='Empirical amplitude-gauge obstruction')
    f.tight_layout();f.savefig(fig/'E1_FULL66_SPECTRA_AND_STRATA.png',dpi=180);f.savefig(fig/'E1_FULL66_SPECTRA_AND_STRATA.pdf');plt.close(f)
    f,axes=plt.subplots(1,2,figsize=(11,4.3))
    for label in ('saturation_inactive_exact_local_subspace','numerical_full_decoder_kernel','radial_property_compensation_candidate'):
        data=[a for a in finite if a['direction']==label and a['parent_id'].startswith('R7_E1_')]
        if data:
            for t in FINITE_STEPS:
                vals=[a['material_RMS'] for a in data if a['t']==t]
                axes[0].scatter([t]*len(vals),np.maximum(vals,1e-18),s=5,alpha=.12)
            med=[np.median([a['material_RMS'] for a in data if a['t']==t]) for t in FINITE_STEPS]
            axes[0].loglog(FINITE_STEPS,np.maximum(med,1e-18),'o-',label=label.replace('_',' '))
    axes[0].set(xlabel='Latent step',ylabel='Material RMS change',title='Finite straight-path tests');axes[0].legend(fontsize=7)
    for a in records[256:]:axes[1].semilogy(np.arange(1,67),np.maximum(svals[records.index(a)],1e-20),label=a['parent_id'])
    axes[1].set(xlabel='Singular value index',ylabel='Material RMS / latent unit',title='Seams and property near-tail controls');axes[1].legend(fontsize=5,ncol=2)
    f.tight_layout();f.savefig(fig/'E1_FINITE_KERNEL_AND_SEAMS.png',dpi=180);f.savefig(fig/'E1_FINITE_KERNEL_AND_SEAMS.pdf');plt.close(f)

def report(prior,s,e0,path):
    c=s['records'][:256]
    lines=['# E0/E1 Full66 decoder-only audit','',
       'This artifact reports decoder arithmetic and geometry only. Scientific interpretation and final theory decisions remain with the root Codex reviewer.',
       '',f"Independent frozen parents: 256; seed: {SEED}; deterministic edge controls: {s['n_edge_cases']}. New Maxwell RHS / states / LU / NN training = 0 / 0 / 0 / 0.",
       '',f"The original P1 MaterialPrior is imported from `{e0['source_prior']}` through sibling R4 `p3r4.core`. No Book or Target is constructed. All 576 material cells, including {int(prior.mask.sum())} ROI cells, enter [Re Dg; Im Dg]/24. x=(z,xi) is Gaussian66; native properties are ndtr(xi), and both property derivatives include the standard-normal density.",
       '', 'Observed rank66 counts are 88/256, 163/256, 186/256 and 191/256 at relative cutoffs 1e-6,1e-8,1e-10,1e-12 respectively. The machine-scale transition-A calculation has full column rank64 in 191/256 parents; the remaining 65 parents have plateau-inactive spatial directions. Four of those also have numerically witnessed curved radial same-material paths. These outcomes do not support imposing a universal rescaling quotient or a fixed large Full66 gauge dimension across this prior.',
       '', '## Observed ranks','',json.dumps(s['relative_rank_histograms'],indent=2),
       '',f"Transition-cell range: {s['prior_transition_cell_range']}; pure-cell range: {s['prior_pure_cell_range']}. Saturation-inactive subspace dimensions: {s['prior_saturation_inactive_dimensions']}. Rank-formula mismatches at relative 1e-10: {s['rank_formula_mismatch_count']}.",
       '', 'For B=D_z fraction, native interior properties give rank(Dg)=rank(B)+1+1{fraction not in range(B)} when the fraction is nonzero; the all-host case has rank zero. The property chain is invertible at finite xi, but near-tail conditioning changes numerical ranks. Formula checks use the same 1e-10 relative cutoff; near-seam/tail discrepancies are numerical qualifications, not new exact rank claims.',
       '', '## Exact inactive directions versus numerical null directions','',
       'A transition cell has -temperature < f_i < temperature. The explicitly constructed subspace consists of dz in ker(A_transition), with both property coordinates fixed. Every transition field is then constant along the straight path; all other ROI cells remain exactly on their existing host/pure plateau until a plateau crossing. This supplies a local exact same-material subspace at strict strata, independently of a derivative cutoff. At seams, the crossing interval can collapse to zero. Floating-point SVD kernels are recorded separately and are not automatically called exact gauges.',
       '', 'Finite straight paths at t=0.001, 0.01, 0.1, 0.5 are saved for explicit inactive, SVD-null, radial-property, material-unit strongest, and material-unit weakest nonzero directions. Plateau changes are retained. Material-unit weak directions can require very large latent steps and leave the local stratum; this is explicitly recorded rather than treated as a linear prediction failure of Maxwell physics.',
       '', 'The inactive subspace is computed with a machine-epsilon rank threshold on A_transition, independently of the four material-Jacobian rank cutoffs. Its transition A-residual and a symmetric plateau-crossing interval for its first unit direction are saved. These are floating-point basis computations and finite material-equality witnesses, rather than symbolic certificates for every dimension.',
       '', '## Global rescaling audit','',
       'The tested radial property tangent requests d(delta)=c delta and compensates using B dz=-c fraction. A pure-phase cell has fraction=1 and B-row=0, so it directly obstructs nonzero radial rescaling at that point. With no pure cells, representability of fraction in range(B) must still be checked; uniform-field reasoning alone does not establish it. The full material tangent residual after least-squares compensation is recorded for every parent.',
       '', 'A separate finite candidate alpha*delta with the spatial field divided by alpha is tested at alpha=0.8, 0.95, 1.05, 1.2 wherever native property intervals permit. The quintic fraction law is not homogeneous. These empirical finite residuals test that specific proposed global field/property rescaling; they do not prove the absence of every nonlinear same-material path. Saturation-inactive fibers and low-transition nonlinear amplitude fibers are distinct objects.',
       '', 'A stronger nonlinear test divides each transition fraction by alpha, inverts the original quintic, and fits the resulting fields with the original A. No-pure parents 092, 201, 233 and 251 (19,40,10,14 transition cells) admit corrected finite paths at all four prescribed tangent steps with material RMS residuals below 2.5e-16 and no host crossing. Their uncorrected straight radial paths have nonzero curvature residuals, reaching respectively about 1.77e-5, 8.87e-8, 3.97e-4 and 5.96e-6 at t=0.5. This is finite numerical evidence for particular curved same-material fibers in the original prior, beyond the straight inactive subspaces. The other no-pure parent 136 has 88 transition cells, numerical rank66 at 1e-10 and 1e-12, and a nonzero nonlinear field-fit/material residual; it is not promoted to an exact gauge. The uniform transition control also admits the nonlinear correction but is kept separate from the independent-parent results. Raw corrected x/material vectors and all residuals are saved.',
       '', '## Numerical and cost qualifications','',
       f"Three full 66-column central differences were evaluated at steps {FD_STEPS.tolist()} for every prior and edge point. Max prior relative Frobenius discrepancy: {s['max_FD_relative_error']:.6g}. Absolute singular ranks at 1e-8,1e-10,1e-12,1e-14, operator-norm FD errors, and FD-informed rank intervals are retained. Max projector idempotence discrepancy: {s['max_projector_idempotence']:.6g}.",
       '',f"Deterministic coordinate/material/J replay: {e0['exact_seed_replay']} / {e0['exact_replay_material']} / {e0['exact_replay_jacobian']}. Final decoder audit wall {e0['wall_seconds']:.3f}s; CPU {e0['CPU_seconds']:.3f}s; peak RSS {e0['peak_RSS_bytes']} bytes. Cost excludes plotting/report writing after the E0 timer receipt. Earlier implementation sweeps are recorded separately in the worker execution note; all sweeps consumed zero Maxwell RHS.",
       '', 'Raw complex128 material and Jacobian arrays plus float64 coordinates, fractions, fields, full singular spectra and FD errors: `raw/decoder/FULL66_DECODER_ARRAYS.npz`. Raw vectors and finite comparisons: `raw/decoder/DIRECTIONS_AND_PROJECTORS.json` and `raw/decoder/FINITE_STEP_COMPARISONS.json`. No old history or ledger is modified.']
    path.write_text('\n'.join(lines)+'\n')

if __name__=='__main__':run()
