"""MAX-STRUCT-R1: disclosed development anchors, no blind inversion.

Finite CM/RR DDA same-field branches; one analytic continuum control.
Run from study root: python src/tensor_families.py
"""
from pathlib import Path
import json, time, resource
import numpy as np
from scipy import linalg as la
from scipy.interpolate import RegularGridInterpolator
from threadpoolctl import threadpool_limits
from vector_backend import TensorDDA, block_diagonal, plane_waves

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/families'
BASE='7317601b43b914525c21f09e5a60020552ea677b'
DIR=np.array([[0,0,1],[0,0,1],[1,0,0],[0,1,0],[1,0,0],[0,1,0]],float)
POL=np.array([[1,0,0],[0,1,0],[0,0,1],[0,0,1],[0,1,0],[1,0,0]],float)
SIGMA=0.001

def dump(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False))

def grid(h):
    a=-.12+h/2+np.arange(round(.24/h))*h
    return a,np.stack(np.meshgrid(a,a,a,indexing='ij'),-1).reshape(-1,3)

def bump(x):
    q=np.sum(x*x,axis=1)/.1**2
    b=np.zeros(len(x));ok=q<1
    b[ok]=np.exp(1-1/(1-q[ok]))
    return b

def sphere(n,r):
    z=1-2*(np.arange(n)+.5)/n
    phi=np.arange(n)*np.pi*(3-np.sqrt(5))
    return r*np.column_stack([np.sqrt(1-z*z)*np.cos(phi),np.sqrt(1-z*z)*np.sin(phi),z])

def rel(a,b):return float(la.norm(a-b)/max(la.norm(a),1e-300))

def legal(eps):
    re=(eps+eps.conj().transpose(0,2,1))/2
    im=(eps-eps.conj().transpose(0,2,1))/(2j)
    return dict(reciprocity_defect=float(la.norm(eps-eps.transpose(0,2,1))),
                real_min=float(np.linalg.eigvalsh(re).min()),real_max=float(np.linalg.eigvalsh(re).max()),
                loss_min=float(np.linalg.eigvalsh(im).min()),loss_max=float(np.linalg.eigvalsh(im).max()))

def rank_info(E,N):
    u=E.reshape(N,3,-1);v=np.concatenate([u.real,u.imag],axis=2)
    sc=np.linalg.svd(u,compute_uv=False);sr=np.linalg.svd(v,compute_uv=False)
    tol=64*np.finfo(float).eps
    rc=np.sum(sc>tol*sc[:,0,None],axis=1);rr=np.sum(sr>tol*sr[:,0,None],axis=1)
    # A single block tensor must annihilate all fields at all cells.
    uc=u.transpose(1,0,2).reshape(3,-1);vr=v.transpose(1,0,2).reshape(3,-1)
    bc=np.linalg.svd(uc,compute_uv=False);br=np.linalg.svd(vr,compute_uv=False)
    rows=np.linalg.norm(u,axis=2)
    return dict(complex_ranks={str(int(k)):int(np.sum(rc==k)) for k in np.unique(rc)},
                real_ranks={str(int(k)):int(np.sum(rr==k)) for k in np.unique(rr)},
                complex_dof_sum=int(np.sum((3-rc)*(4-rc)//2)),
                real_dof_sum=int(np.sum((3-rr)*(4-rr)//2)),
                block_constant_complex_rank=int(np.sum(bc>tol*bc[0])),
                block_constant_real_rank=int(np.sum(br>tol*br[0])),
                axis_diagonal_null_dof=int(np.sum(rows<=tol*np.max(rows,axis=1)[:,None])),
                relative_rank_threshold=tol,rank_threshold_semantics='roundoff diagnostic, not weak-singular-value truncation',
                smallest_returned_complex_singular_value=float(sc[:,-1].min()),
                block_complex_singular_values=bc.tolist(),block_real_singular_values=br.tolist())

def base_receipt(case,h,N):
    return dict(run_id=case,case_id=case,branch='F',git_sha=BASE,code_hash=None,config_hash=None,
       hash_reason='User forbids hash verification; src/tensor_families.py and configs/initial_plan.json retained',
       seed=20261102,data_id=case,material_class=None,
       physical_convention='exp(-i omega t); lambda0=epsilon_b=mu_b=Z0=1',
       mesh=dict(side=.24,h=h,N=N,discretization='CM+RR dyadic DDA',PML='not applicable; outgoing free-space kernel'),
       sources=dict(directions=DIR.tolist(),polarizations=POL.tolist(),old_columns=[0,1],new_columns=[2,3,4,5]),
       receivers=dict(old=sphere(12,.6).tolist(),added=sphere(24,.8).tolist(),type='external near-field xyz E'),
       frequencies=[1.],noise_model=dict(type='deterministic complex Frobenius ball',sigma=SIGMA,radius='sigma*sqrt(m_complex)',Gaussian_probability=None),
       whitening='W=I/sigma; complex norm (not sqrt(2) real Gaussian whitening)',prior='R=0 inside disclosed reciprocal tensor set',
       initialization='public analytic/development anchor',state_representation='all vector internal E; j=(eps-I)E; H via analytic curl dipole kernel',
       state_dof=6*N*2,material_dof=None,baseline_id='not an imaging comparison; frozen scalar CIE incompatible with tensor task',
       stage='OFFLINE_DIAGNOSTIC',termination_reason='PREDECLARED_FINITE_BRANCHES_COMPLETE',online_selection_rule='none; retain both branches',
       physical_residual=None,auxiliary_residual=None,constitutive_defect=None,constraint_violation=None,
       cost=None,candidate_ids=['anchor','branch'],truth_access='DISCLOSED_DEV_ANCHOR; no hidden GT',
       validation=dict(continuum='analytic real example only; passive branch finite DDA',independent_reference='Mie control in results/backend_qualification'))

def run_case(kind,h):
    tic=time.perf_counter();cpu=time.process_time();a,points=grid(h);N=len(points)
    case=f'F_{kind}_h{h:.2f}';rec=base_receipt(case,h,N);costs=[];extras=0
    mod=TensorDDA(points,h**3);b=bump(points)
    if kind=='analytic_real':
        eps=np.repeat(np.eye(3)[None,:,:],N,axis=0).astype(complex)
        eps[:,2,2]+=b
        V=np.zeros_like(eps);V[:,2,2]=2*b
        eps2=eps+V
        rec.update(material_class='real symmetric positive definite; axis diagonal',material_dof=6*N)
    else:
        T=np.array([[2.2,.12,.08],[.12,2.7,.09],[.08,.09,3.1]],complex)+.25j*np.eye(3)
        eps=np.repeat(T[None,:,:],N,axis=0)
        rec.update(material_class='complex symmetric strictly passive',material_dof=12*N)
    s0=mod.solve(eps,DIR,POL,OUT/f'{case}_anchor_solver.json');costs.append(s0.receipt)
    Eold=s0.E[:,:2].reshape(N,3,2)
    if kind=='passive_complex':
        q=np.cross(Eold[:,:,0],Eold[:,:,1])
        raw=b[:,None,None]*q[:,:,None]*q[:,None,:]  # transpose dyad, NOT Hermitian
        scale=.12/np.linalg.norm(raw,ord=2,axis=(1,2)).max()
        V=scale*raw;eps2=eps+V
    else:q=np.tile([0,0,1.],(N,1));scale=2.
    s1=mod.solve(eps2,DIR,POL,OUT/f'{case}_branch_solver.json');costs.append(s1.receipt)
    S=mod.external_electric_operator(sphere(12,.6));Sa=mod.external_electric_operator(sphere(24,.8))
    Y0=S@s0.j;Y1=S@s1.j;Ya0=Sa@s0.j;Ya1=Sa@s1.j
    Q=np.array([[1,1j],[1j,1]],complex)/np.sqrt(2)
    # Full retained factor solves with actually mixed incident columns; no new LU.
    inc=s0.incident[:,:2]@Q
    mix0=la.lu_solve(s0.factor,inc);mix1=la.lu_solve(s1.factor,inc);extras+=4
    Xm0=block_diagonal(eps-np.eye(3));Xm1=block_diagonal(eps2-np.eye(3))
    mixgap=float(la.norm(S@(Xm1@mix1-Xm0@mix0))/SIGMA)
    A=(Y1[:,2:]-Y0[:,2:])/SIGMA
    _,sv,vh=la.svd(A,full_matrices=False)
    c=vh[0].conj();pred=float(la.norm(A@c))
    m=Y0[:,:2].size;noise_radius=np.sqrt(m)
    rank=rank_info(s0.E[:,:2],N)
    leg=legal(eps2)
    support=dict(real_positive=leg['real_min']>0,passive=leg['loss_min']>=0,nontrivial=bool(la.norm(V)>0),algebraic_annihilator='exact by construction in real arithmetic; floating defects reported, not threshold-tuned')
    rec.update(status='PASS' if all(support[k] for k in ['real_positive','passive','nontrivial']) else 'FAIL',support_checks=support,scientific_scope='finite-model same-field construction, not reconstruction or continuum certification',
       physical_residual=max(x['linear_relative_residual'] for x in costs),
       constitutive_defect=max(x['constitutive_relative_residual'] for x in costs),
       same_field_annihilator_defect=float(la.norm(np.einsum('nij,njs->nis',V,Eold))),
       constraint_violation=0. if legal(eps2)['real_min']>=1-1e-12 and legal(eps2)['loss_min']>=-1e-12 else None,
       real_or_complex_rank=rank, tensor_legal_margin=dict(anchor=legal(eps),branch=legal(eps2)),
       branch_distance=dict(frobenius_rms=float(la.norm(V)/np.sqrt(N)),real_rms=float(la.norm(V.real)/np.sqrt(N)),loss_rms=float(la.norm(V.imag)/np.sqrt(N)),normalization_amplitude=1.),
       old_source_gap=dict(E_relative=rel(s0.E[:,:2],s1.E[:,:2]),H_relative=rel(s0.H[:,:2],s1.H[:,:2]),
           j_absolute=float(la.norm(s1.j[:,:2]-s0.j[:,:2])),port_whitened=float(la.norm(Y1[:,:2]-Y0[:,:2])/SIGMA)),
       A_extra_receivers=dict(whitened_gap=float(la.norm(Ya1[:,:2]-Ya0[:,:2])/SIGMA),added_complex_channels=int(Ya0[:,:2].size)),
       B_source_mix=dict(whitened_gap=mixgap,unitary=Q.tolist() if False else '[[1,i],[i,1]]/sqrt(2)',noise_covariance='sigma^2 I preserved',RHS=4),
       new_source_gap=dict(each_library_whitened=np.linalg.norm(A,axis=0).tolist(),best_coherent_prediction=pred,
           coefficient_real=c.real.tolist(),coefficient_imag=c.imag.tolist(),power_metric='I: normalized ideal external plane-wave channel coefficients; not a finite antenna power calibration',
           single_source_noise_radius=float(np.sqrt(Y0.shape[0])),two_ball_disjoint=bool(pred>2*np.sqrt(Y0.shape[0])),
           actual_new_measurement='NOT_RUN; prediction-only development demonstration',candidate_library_contains_truth='not asserted'),
       full_original_objective=dict(Y='clean disclosed anchor prediction',anchor=0.,branch=.5*float(la.norm(Y1[:,:2]-Y0[:,:2])/SIGMA)**2,
           discrepancy_radius=float(noise_radius),both_qualify=bool(la.norm(Y1[:,:2]-Y0[:,:2])/SIGMA<=noise_radius)),
       far_field_gap=float(la.norm(s1.far(sphere(20,1.))['E'][:,:2]-s0.far(sphere(20,1.))['E'][:,:2])),
       E_state_dof_all_six_solved=6*N*6,
       cost=dict(stage='validation',LU=sum(z['LU'] for z in costs),forward_RHS=sum(z['RHS'] for z in costs),additional_triangular_RHS=extras,
            FFT=0,elliptic_solves=0,local_rank_SVD=2*N+2,Green_kernel_assemblies=1,setup=mod.setup_receipt,
            wall_s=time.perf_counter()-tic,cpu_s=time.process_time()-cpu,memory_peak_process_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            GPU_s=0,currency=None,currency_reason='no billed-service usage meter available'))
    dump(OUT/f'{case}.json',rec)
    np.savez_compressed(OUT/f'{case}.npz',points=points,epsilon=eps,branch=eps2,V=V,E=s0.E,H=s0.H,j=s0.j,E_branch=s1.E,H_branch=s1.H,Y=Y0,Y_branch=Y1,q=q,scale=scale,axis=a)
    public=ROOT/'data_public/dev_anchors';public.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(public/f'{case}.npz',points=points,epsilon=eps,Y=Y0[:,:2],directions=DIR[:2],polarizations=POL[:2],receivers=sphere(12,.6))
    rec['cost']['wall_including_npz_io_s']=time.perf_counter()-tic
    rec['cost']['cpu_including_npz_io_s']=time.process_time()-cpu
    rec['cost']['wall_s_scope']='assembly+compute before array output; summary total includes NPZ I/O'
    dump(OUT/f'{case}.json',rec)
    return rec

def interpolation_audit():
    """Fixed coarse direction transported to refined model, not re-nulled."""
    start=time.perf_counter();cpu=time.process_time()
    c=np.load(OUT/'F_passive_complex_h0.06.npz');f=np.load(OUT/'F_passive_complex_h0.03.npz')
    ax=c['axis'];pf=f['points'];nf=len(pf)
    interp=RegularGridInterpolator((ax,ax,ax),c['q'].reshape(len(ax),len(ax),len(ax),3),bounds_error=False,fill_value=None)
    qc=interp(pf);bf=bump(pf)
    Vc=float(c['scale'])*bf[:,None,None]*qc[:,:,None]*qc[:,None,:]
    fineq=f['q'];Vref=float(c['scale'])*bf[:,None,None]*fineq[:,:,None]*fineq[:,None,:]
    mod=TensorDDA(pf,.03**3)
    s=mod.solve(f['epsilon']+Vc,DIR[:2],POL[:2],OUT/'F_interpolated_solver.json')
    S=mod.external_electric_operator(sphere(12,.6))
    gap=float(la.norm(S@s.j-f['Y'][:,:2])/SIGMA)
    r=base_receipt('F_fixed_coarse_branch_interpolation',.03,nf)
    r.update(status='INCONCLUSIVE',termination_reason='ONE_NECESSARY_REFINEMENT_COMPLETE_NO_CONTINUUM_RATE_CLAIM',
       material_class='complex symmetric strictly passive',material_dof=12*nf,
       physical_residual=s.receipt['linear_relative_residual'],constitutive_defect=s.receipt['constitutive_relative_residual'],
       same_field_annihilator_defect=float(la.norm(np.einsum('nij,njs->nis',Vc,f['E'][:,:2].reshape(nf,3,2)))),
       constraint_violation=None,tensor_legal_margin=legal(f['epsilon']+Vc),
       branch_distance=float(la.norm(Vc)/np.sqrt(nf)),old_source_gap=gap,new_source_gap=None,
       direction_relative_error_to_fine_field_construction=float(la.norm(Vc-Vref)/la.norm(Vref)),
       transport='coarse q linearly interpolated/extrapolated then analytic bump times q q^T, not direct tensor interpolation; extrapolated centers lie outside bump support',
       note='Coarse and fine re-nulled branches are different discretizations. This fixed interpolant is NOT exact on fine grid. Two grids do not establish convergence rate or continuum noise certification.',
       cost=dict(stage='validation',LU=1,forward_RHS=2,additional_triangular_RHS=0,setup=mod.setup_receipt,FFT=0,elliptic_solves=0,wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu))
    dump(OUT/'F_interpolation.json',r)
    return r

def finite_noise_ray():
    """One nonzero weak excitation ray; uniform Neumann bound, no grid sweep."""
    tic=time.perf_counter();cpu=time.process_time()
    _,p=grid(.06);N=len(p);mod=TensorDDA(p,.06**3)
    data=np.load(OUT/'F_passive_complex_h0.06.npz');eps=data['epsilon'];pol=POL[:3].copy();pol[2]*=.02
    s=mod.solve(eps,DIR[:3],pol,OUT/'F_noise_anchor_solver.json')
    u=s.E.reshape(N,3,3);P,sv,_=np.linalg.svd(u)
    weak=P[:,:,-1].conj();V=bump(p)[:,None,None]*weak[:,:,None]*weak[:,None,:]
    V/=la.norm(V);Vm=block_diagonal(V)
    # Reuse field LU to obtain current resolvent: (I-XG)^-1=I+X(I-GX)^-1G.
    X=block_diagonal(eps-np.eye(3));Rcur=np.eye(3*N)+X@la.lu_solve(s.factor,mod.G)
    S=mod.external_electric_operator(sphere(12,.6))/SIGMA
    a=float(la.norm(Vm@mod.G@Rcur,'fro'))
    tlegal=.10/np.linalg.norm(V,ord=2,axis=(1,2)).max()
    T=min(tlegal,.4/a);K=float(la.norm(S@Rcur,'fro')/(1-T*a))
    excite=float(la.norm(Vm@s.E));radius=np.sqrt(36*3);r0=.5*radius
    t=min(T,(radius-r0)/(K*excite))
    st=mod.solve(eps+t*V,DIR[:3],pol,OUT/'F_noise_branch_solver.json')
    ya=S@s.j;yt=S@st.j
    rng=np.random.default_rng(20261102);noise=rng.normal(size=ya.shape)+1j*rng.normal(size=ya.shape)
    noise*=r0/la.norm(noise);y=ya+noise
    # Local excitation spectrum only: it excludes receiver and resolvent.
    spectrum=np.column_stack([sv,np.sqrt((sv[:,0]**2+sv[:,1]**2)/2),np.sqrt((sv[:,0]**2+sv[:,2]**2)/2),np.sqrt((sv[:,1]**2+sv[:,2]**2)/2)])
    r=base_receipt('F_finite_noise_weak_ray',.06,N)
    supported=(la.norm(yt-ya)<=t*K*excite and la.norm(yt-y)<=radius and legal(eps+t*V)['loss_min']>0)
    r.update(status='PASS' if supported else 'FAIL',support_checks=dict(pair_bound_holds=bool(la.norm(yt-ya)<=t*K*excite),actual_observation_both_fit=bool(la.norm(yt-y)<=radius),strict_passivity=bool(legal(eps+t*V)['loss_min']>0)),scientific_scope='floating finite-DDA bound on one ray; no outward-rounded continuum certificate',material_class='complex symmetric strictly passive',material_dof=12*N,
       sources=dict(directions=DIR[:3].tolist(),polarizations=pol.tolist(),third_amplitude=.02),
       physical_residual=st.receipt['linear_relative_residual'],constitutive_defect=st.receipt['constitutive_relative_residual'],
       unit_ray_excitation=excite,state_dof=6*N*3,
       constraint_violation=0,tensor_legal_margin=legal(eps+t*V),real_or_complex_rank=rank_info(s.E,N),
       branch_distance=dict(frobenius=float(t),frobenius_rms=float(t/np.sqrt(N)),norm_V_mat=1.,norm_definition='unweighted global Frobenius'),
       old_source_gap=float(la.norm(yt-ya)),new_source_gap=None,
       noise_qualification=dict(K=K,Neumann_coefficient=a,interval=[0,float(T)],t_legal=float(tlegal),t_selected=float(t),
           excitation=excite,pair_bound=float(t*K*excite),pair_actual=float(la.norm(yt-ya)),noise_radius=float(radius),
           anchor_residual=float(la.norm(ya-y)),branch_residual=float(la.norm(yt-y)),
           sufficient_slack_bound=float(radius-r0),both_fit=bool(la.norm(yt-y)<=radius),
           interpretation='constructive ambiguity lower bound, not recovery upper bound; K only uniform on this ray; numerical norms not interval certified'),
       local_excitation_spectrum_min=np.min(spectrum,axis=0).tolist(),
       cost=dict(stage='validation',LU=2,forward_RHS=6,additional_triangular_RHS=3*N,local_SVD=N,setup=mod.setup_receipt,FFT=0,elliptic_solves=0,wall_s=time.perf_counter()-tic,cpu_s=time.process_time()-cpu))
    dump(OUT/'F_noise_ray.json',r);return r

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();cpu=time.process_time()
    with threadpool_limits(limits=1):
        rows=[run_case(kind,h) for kind in ['analytic_real','passive_complex'] for h in [.06,.03]]
        inter=interpolation_audit();noise=finite_noise_ray()
    summary=dict(status='COMPLETED',interpretation='no imaging advantage claim',rows=rows,interpolation=inter,finite_noise=noise,
                 wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu)
    dump(OUT/'SUMMARY.json',summary)
    print(json.dumps({'rows':[{k:r[k] for k in ['run_id','branch_distance','old_source_gap','new_source_gap']} for r in rows],
                      'interpolation_gap':inter['old_source_gap'],'noise':noise['noise_qualification'],'wall_s':summary['wall_s']}))

if __name__=='__main__':main()
