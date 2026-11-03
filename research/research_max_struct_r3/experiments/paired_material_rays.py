"""Bounded paired material rays on the fixed1104nearest-neighbor edges.

Both directions (1,1)/(1,-1); roots of exact rational data derivative, not a
t-grid. Full unchanged prior scores finite candidates. Numerical roots are
not certified, and data-stationary candidates need not optimize totalPhi.
"""
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
import json,sys,time,traceback,resource
import numpy as np
from scipy import linalg as la
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];R2=REPO/'research/research_max_struct_r2'
sys.path.insert(0,str(REPO));sys.path.insert(0,str(ROOT/'experiments'))
from finite_coordinate_paths import resolvent,fixed_edges
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_fbg_r3_experiments.src.full_maxwell_corrector import corrector
from research.research_max_struct_r2.experiments.run_blind import metric,dump,compact

CASE='obs_b20d5a55677d4eca';REAL_TOL=128*np.finfo(float).eps

def coefficients(C,A,BO,BQ,rnorm,s):
    """Norm polynomial of p*r+t*A1+t²*A2 using2×2 Gram contractions."""
    adj=np.array([[C[1,1],-C[0,1]],[-C[1,0],C[0,0]]])
    K1=np.diag([1.,s]);K2=-s*adj
    p=np.array([1.,-(C[0,0]+s*C[1,1]),s*la.det(C)],complex)
    a1=np.vdot(A,K1);a2=np.vdot(A,K2)
    b11=np.trace(K1.conj().T@BO@K1@BQ)
    b12=np.trace(K1.conj().T@BO@K2@BQ)
    b22=np.trace(K2.conj().T@BO@K2@BQ)
    g01=p[1]*rnorm+a1;g02=p[2]*rnorm+a2
    g11=abs(p[1])**2*rnorm+2*np.real(p[1].conj()*a1)+b11.real
    g12=p[1].conj()*p[2]*rnorm+p[1].conj()*a2+p[2]*a1.conj()+b12
    g22=abs(p[2])**2*rnorm+2*np.real(p[2].conj()*a2)+b22.real
    N=np.array([rnorm,2*g01.real,2*g02.real+g11,2*g12.real,g22],float)
    D=np.polynomial.polynomial.polymul(p,p.conj()).real
    return N,D,p,K1,K2

def scaled_derivative(N,D,lo,hi):
    center=(lo+hi)/2;half=(hi-lo)/2
    change=np.polynomial.Polynomial([center,half])
    nx=np.polynomial.Polynomial(N)(change).coef;dx=np.polynomial.Polynomial(D)(change).coef
    # Explicit cancellation k=l removes the degree7 term structurally.
    f=np.zeros(7)
    for k in range(len(nx)):
        for l in range(len(dx)):
            if k!=l:f[k+l-1]+=(k-l)*nx[k]*dx[l]
    scale=max(abs(f))
    if scale==0:return center,half,f,np.empty(0,complex),0
    f/=scale;trim=np.trim_zeros(f,'b');degree=len(trim)-1
    roots=np.roots(trim[::-1])if degree>0 else np.empty(0,complex)
    return center,half,f,roots,degree

def root_residual(x,f):
    value=abs(np.polynomial.polynomial.polyval(x,f))
    denominator=np.sum(abs(f)*abs(x)**np.arange(len(f)))
    return float(value),float(value/max(denominator,np.finfo(float).tiny))

def run(arm='paired_neighbor_rays',edge_selector=None):
    directory=ROOT/'results/development'/arm;directory.mkdir(parents=True,exist_ok=True)
    path=directory/'RECEIPT.json'
    if path.exists():raise FileExistsError('Preserve single pairedray disposition')
    start=time.perf_counter();cpu=time.process_time();model=None;cost=Counter()
    rec=dict(status='RUNNING',run_id=arm,case_id=CASE,
        classification='FORMULA_DRIVEN_POST_HOC_DEVELOPMENT',code='experiments/paired_material_rays.py',
        started_utc=datetime.now(timezone.utc).isoformat(),truth_access='NONE',seed='same scaled_original_restart endpoint',
        dictionary='1104fixednearestneighbor pairs ×signpatterns(1,1)/(1,-1), no t-grid'if edge_selector is None else 'one dressedfullwave interaction neighbor percell, unionunorderedpairs; both signs, no t-grid',
        objective='unchanged original normalizedphysicalPhi+fullprior',max_finite_moves=1,
        corrector='same originalphysical L-BFGS-B1000 only afteractual strictlyimproving finiteproposal',
        root_coordinates='eachlegalinterval mapped to[-1,1] before differentiation',
        real_root_classification_tolerance=REAL_TOL,
        root_qualification='numericalroots andresiduals only; no certifiedrootcompleteness, no totalPhi coordinateoptclaim',
        threads=1,device='CPU',hash_checks=False)
    dump(path,rec)
    try:
        with np.load(R2/'data_public'/(CASE+'.npz'),allow_pickle=False)as f:public={k:f[k]for k in f.files}
        with np.load(ROOT/'results/development/scaled_original_restart/prediction.npz',allow_pickle=False)as f:
            m=f['m'].copy();seed_gradient=f['gradient'].copy()
        model=MaxwellModel(public);Y=public['Y'];initial=metric(model,Y,m);state,C,O=resolvent(model,m,cost)
        r=state['Y']-Y;Q=state['Q'];rnorm=float(np.vdot(r,r).real);scale=float(np.vdot(Y,Y).real)
        prior=float(model.reg.evaluate(m)[0]);base=initial['objective']
        A=(O.conj().T@r)@Q.conj().T;BO=O.conj().T@O;BQ=Q@Q.conj().T
        cost.update(rank2_cross_Gram_builds=1,rank2_O_Gram_builds=1,rank2_Q_Gram_builds=1)
        edges=fixed_edges(model.shape)
        if edge_selector is not None:
            edges,metadata=edge_selector(model,C,O,Q,BO,BQ,public,cost,directory)
            rec['interaction_graph']=metadata
            if metadata.get('all_in_previous_nearest_dictionary'):
                rec['reused_previous_nearest_ray_decisions']=True
                rec['reuse_source']='results/development/paired_neighbor_rays; exact same savedseed'
                edges=[]
        rays=[];tables=[];polyN=[];polyD=[];polyF=[];root_arrays=[];best=None;all_root_res=[]
        for edge,(i,j)in enumerate(edges):
            cells=np.array([i,j]);ca=C[np.ix_(cells,cells)];ab=A[np.ix_(cells,cells)];ob=BO[np.ix_(cells,cells)];qb=BQ[np.ix_(cells,cells)]
            for s in (1,-1):
                cost['paired_direction_cases']+=1
                lo=max(.02-m[i],.02-m[j]if s==1 else m[j]-20.)
                hi=min(20.-m[i],20.-m[j]if s==1 else m[j]-.02)
                ray=dict(edge=edge,cells=[i,j],sign=s,legal_interval=[float(lo),float(hi)],candidate_count=0)
                if hi<=lo:
                    cost['zero_width_directions']+=1;ray['degree']=None;ray['stop']='ZERO_WIDTH_INTERVAL';rays.append(ray);continue
                N,D,p,K1,K2=coefficients(ca,ab,ob,qb,rnorm,s)
                center,half,f,roots,degree=scaled_derivative(N,D,lo,hi)
                cost['quartic_N_D_problems']+=1;cost['degree_at_most6_root_problems']+=1
                cand=[(float(lo),'lower_bound',None),(float(hi),'upper_bound',None),(0.,'zero',None)]
                rr=[]
                for x in roots:
                    residual,backward=root_residual(x,f)
                    numericallyreal=bool(abs(x.imag)<=REAL_TOL*max(1.,abs(x.real)))
                    rr.append(dict(x_real=float(x.real),x_imag=float(x.imag),absolute_residual=residual,
                        normalized_backward_residual=backward,numerically_real=numericallyreal))
                    if numericallyreal and -1<=x.real<=1:
                        t=float(center+half*x.real);cand.append((t,'data_stationary_root',root_residual(x.real,f)))
                    all_root_res.append(backward)
                ray.update(degree=degree,root_diagnostics=rr,center=float(center),half_width=float(half))
                rootpad=np.full(6,np.nan+1j*np.nan);rootpad[:len(roots)]=roots
                polyN.append(N);polyD.append(D);polyF.append(f);root_arrays.append(rootpad)
                ray['polynomial_row']=len(polyN)-1;seen=set();bestlocal=None
                for t,kind,rootres in cand:
                    if t in seen:continue
                    seen.add(t);cost['finite_ray_candidates']+=1;ray['candidate_count']+=1
                    trial=m.copy();trial[i]=np.clip(m[i]+t,.02,20.);trial[j]=np.clip(m[j]+s*t,.02,20.)
                    # Bounds only correct floating representation; this ray uses
                    # common t, not independently optimized cell changes.
                    den=np.polynomial.polynomial.polyval(t,p)
                    if den==0:cost['singular_candidates']+=1;continue
                    K=(t*K1+t*t*K2)/den;cost['closed_form2x2_transfer_evaluations']+=1
                    cross=2*np.vdot(ab,K).real;norm=np.trace(K.conj().T@ob@K@qb).real
                    pnew=float(model.reg.evaluate(trial)[0]);cost['full_prior_candidate_evaluations']+=1
                    val=base+.5*(cross+norm)/scale+pnew-prior
                    if not np.isfinite(val):cost['nonfinite_candidates']+=1;continue
                    tables.append((edge,i,j,s,t,val));move=float(np.linalg.norm(trial-m))
                    if move==0:continue
                    item=dict(edge=edge,coordinates=[i,j],sign=s,t=t,candidate_kind=kind,
                        predicted_Phi=float(val),material_distance_l2=move,
                        legal_interval=[float(lo),float(hi)],scaled_polynomial_degree=degree,
                        scaled_root_residual=rootres,denominator_squared=float(abs(den)**2),trial=trial)
                    if bestlocal is None or val<bestlocal['predicted_Phi']:bestlocal=item
                    if best is None or val<best['predicted_Phi']:best=item
                ray['best_nonzero_Phi']=bestlocal['predicted_Phi']if bestlocal else None;rays.append(ray)
        np.savez_compressed(directory/'RAY_POLYNOMIALS.npz',N_t=np.asarray(polyN),D_t=np.asarray(polyD),
            derivative_scaled_x=np.asarray(polyF),numerical_roots_x=np.asarray(root_arrays),candidate_table=np.asarray(tables),anchor_m=m)
        dump(directory/'RAY_ROOT_DIAGNOSTICS.json',rays)
        rec['max_numerical_root_backward_residual']=max(all_root_res)if all_root_res else None
        incumbents=[initial];accepted=False
        if best is not None:
            trial=best.pop('trial');actual=metric(model,Y,trial);cost['chosen_candidate_full_forward_validations']+=1
            rec.update(best_candidate=best,actual_candidate_Phi=actual['objective'],
                predicted_full_forward_Phi_difference=actual['objective']-best['predicted_Phi'])
            if actual['objective']<base:
                accepted=True;incumbents.append(actual);cost['accepted_finite_moves']+=1
                np.savez_compressed(directory/'finite_proposal.npz',m=trial)
                corr=corrector(model,Y,trial,dict(method='lbfgs',max_nfev=1000,box=(.02,20.)))
                dump(directory/'CORRECTOR.json',compact(corr));incumbents.append(metric(model,Y,corr['m']))
                rec.update(corrector_stop=corr['stop'],corrector_cost=corr['cost'])
        bestphys=min(incumbents,key=lambda c:c['objective']);state=model.state(bestphys['m'])
        if rec.get('reused_previous_nearest_ray_decisions'):
            value,g=initial['objective'],seed_gradient
        else:value,g=model.data_value_gradient(bestphys['m'],Y)
        pg=bestphys['m']-np.clip(bestphys['m']-g,.02,20.)
        np.savez_compressed(directory/'prediction.npz',m=bestphys['m'],Y_prediction=state['Y'],j=state['j'],Q=state['Q'],gradient=g,projected_gradient=pg)
        rec.update(status='SEALED_DEVELOPMENT_PREDICTION',finished_utc=datetime.now(timezone.utc).isoformat(),
            prediction=str((directory/'prediction.npz').relative_to(ROOT)),accepted_finite_move=accepted,
            physical_objective=float(value),initial_Phi=base,physical_residual=float(np.linalg.norm(state['Y']-Y)/np.linalg.norm(Y)),
            physical_whitened_residual=float(np.linalg.norm(state['Y']-Y)/float(public['sigma'])),
            projected_gradient_inf=float(np.linalg.norm(pg,np.inf)),material_distance_from_seed_l2=float(np.linalg.norm(bestphys['m']-m)),
            qualification='nearestpairedrays/data-stationaryfinitecandidates with fullprior scores; not certifiedalltotalPhiminima or globalmaterialsearch')
    except Exception as exc:rec.update(status='FAILED_DEVELOPMENT_DISPOSITION',error=repr(exc),traceback=traceback.format_exc())
    finally:
        peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        rec.update(wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,
            physical_model_cost=dict(model.cost)if model else {},auxiliary_cost=dict(cost),FFT=0,GPU=0,
            process_peak_memory_bytes=int(peak if sys.platform=='darwin' else peak*1024),
            timer_scope='load,onefullresolvent,allpolynomials/roots/fullpriors,onefullcandidatevalidation,conditionalcorrector,finalgradient,serialization;receiptwriteexcluded')
        dump(path,rec)
    print(json.dumps({k:rec.get(k)for k in ['run_id','status','accepted_finite_move','physical_residual','projected_gradient_inf','wall_seconds','physical_model_cost','auxiliary_cost','error']}),flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=1):run()
