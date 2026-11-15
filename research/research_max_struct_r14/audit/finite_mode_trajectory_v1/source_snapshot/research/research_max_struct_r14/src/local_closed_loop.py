"""Declared finite-root mechanism: ideal four metric maps vs one native PR/readout.

No inverse benchmark or truth input. The prescribed physical root is used only
for OFFLINE_DIAGNOSTIC directional comparisons, including nonadditivity.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import argparse
from collections import Counter
from pathlib import Path
import sys,time,json,csv,traceback

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
sys.path[:0]=[str(ROOT/'src'),str(REPO),str(REPO/'vendor')]
import numpy as np
from scipy.linalg import null_space,cho_factor,cho_solve,eigh
from threadpoolctl import threadpool_limits
from research.research_fbg_r3_experiments.src.physics import MaxwellModel,make_grid
from research.research_max_struct_r8.src.benchmark_solver import arm_config
from dualmetric_core import StatefulDualMetric,_save_state_npz


def plain(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,Path):return str(v)
    raise TypeError(type(v).__name__)


def dump(p,v):p.write_text(json.dumps(v,default=plain,indent=2,allow_nan=False)+'\n')


def norm(a):return float(np.linalg.norm(a))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,default=ROOT/'audit/local_closed_loop_v1')
    args=ap.parse_args();out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    (out/'states').mkdir();started_cpu=time.process_time();started_wall=time.perf_counter()
    costs=Counter();trials=[];tables=[];errors=[];summary={}
    try:
      with threadpool_limits(limits=1):
        setup_cpu=time.process_time();setup_wall=time.perf_counter()
        public=make_grid(8,sources=6,receivers=16);public['coords']=public['centers'].copy()
        costs['finite_grid_operator_setups']+=1
        model=MaxwellModel(public);model.reg.smooth=0.;model.reg.tv=0.
        radius=np.linalg.norm(public['centers'],axis=1)
        m=.4+8.*((radius>.19)&(radius<.34))
        state=model.state(m);J=state['j'].copy();Q=state['Q'].copy();Y=state['Y'].copy()
        cfg=arm_config('a_plus',dict(host=.1,material_box=[.02,20.],dpcs_rank=16),Y,model.E)
        cfg['stages']=[dict(mf=4,maxiter=1,beta_current=1.,beta_material=.2)]
        cfg.update(fft_geometry=dict(coords=public['coords'],radius=float(public['radius']),k=float(public['k'])),
            regularizer=lambda x:model.reg.evaluate(x)[:2],aplus_regularizer=model.reg,
            sourcewise_cg=False,aplus_mode='mobius',aplus_arc_candidates=True,
            aplus_inherited_gtol=True,aplus_sweeps=2,aplus_step_tol=1e-8)
        base=StatefulDualMetric(model,Y,m,cfg,pair='01')
        # Native initialization is cold; replace only coefficients with the
        # declared physical root for this offline diagnostic, using its full U.
        UU=base.U.apply(np.eye(base.U.k,dtype=complex));costs['root_native_U_columns']+=base.U.k
        A,_,rank,sv=np.linalg.lstsq(UU,J-base.JD,rcond=1e-12)
        costs['root_full_U_lstsq_SVD']+=1;costs['root_full_U_lstsq_RHS']+=model.ns
        base.A=A.copy();base.J=base.JD+base.U.apply(A)
        base.Q=model.E+base.GJD+base.green.apply(base.J-base.JD)
        base.old_grad=base.old_direction=None
        root_errors=dict(current_reconstruction_relative=norm(base.J-J)/norm(J),
            field_reconstruction_relative=norm(base.Q-Q)/norm(Q),
            constitutive_relative=norm(base.J-m[:,None]*base.Q)/norm(J),
            physical_relative=norm(state['L']@J-m[:,None]*model.E)/norm(m[:,None]*model.E),
            native_full_U_rank=int(rank),native_full_U_columns=base.U.k,
            receiver_rank=base.receiver_rank,receiver_smallest_relative=float(base.receiver_singular_values[-1]/base.receiver_singular_values[0]))
        base.checkpoint(out/'ROOT_FULL_STATE.npz')
        np.savez_compressed(out/'ROOT_PUBLIC_AND_PHYSICAL.npz',**public,m_anchor=m,J_anchor=J,Q_anchor=Q,Y_anchor=Y,
                            root_coefficients=A,root_U_singular_values=sv)
        setup=dict(CPU_seconds=time.process_time()-setup_cpu,wall_seconds=time.perf_counter()-setup_wall,
                   inherited_setup_counts=dict(base.cost),physical_cost=dict(model.cost),root_errors=root_errors)
        U=null_space(base.Vd.conj().T);costs['independent_current_nullspace_SVD']+=1
        L=np.eye(model.n)-m[:,None]*model.G;K=L@U;D=model.S@U
        P0=np.diag(1/(1+m));P1=base.blocks.apply(m,.2,np.eye(model.n,dtype=complex))
        costs['ideal_onehop_transform_columns']+=model.n
        metrics={'0':P0.conj().T@P0,'1':P1.conj().T@P1}
        ws,wd=base.ws.copy(),base.wd.copy();ideal={};B={};X={};ideal_rows=[]
        ideal_cpu=time.process_time();ideal_wall=time.perf_counter()
        current_residuals=[]
        for c in '01':
            M=metrics[c];Xs=[]
            for s in range(model.ns):
                H=ws[s]*(K.conj().T@M@K)+wd[s]*(D.conj().T@D)
                H=(H+H.conj().T)/2
                rhs=ws[s]*K.conj().T@M@np.diag(Q[:,s])
                fac=cho_factor(H,lower=True,check_finite=False);x=cho_solve(fac,rhs,check_finite=False)
                costs['ideal_current_normal_Cholesky']+=1;costs['ideal_current_normal_RHS']+=model.n
                Xs.append(x);current_residuals.append(norm(H@x-rhs)/norm(rhs))
            X[c]=np.array(Xs)
        for material in '01':
            M=metrics[material]
            b=sum(ws[s]*np.real(Q[:,s].conj()[:,None]*M*Q[:,s][None,:]) for s in range(model.ns))
            B[material]=(b+b.T)/2
        vals,vec=eigh(B['0']);Bc_sqrt=(vec*np.sqrt(vals))@vec.T;Bc_invsqrt=(vec/np.sqrt(vals))@vec.T
        for pair in ('00','01','10','11'):
            M=metrics[pair[1]];b=B[pair[1]]
            numerator=sum(ws[s]*np.real(Q[:,s].conj()[:,None]*(M@K@X[pair[0]][s])) for s in range(model.ns))
            fac=cho_factor(b,lower=True,check_finite=False);T=cho_solve(fac,numerator,check_finite=False)
            costs['ideal_material_normal_Cholesky']+=1;costs['ideal_material_normal_RHS']+=model.n
            ideal[pair]=T;ev=np.linalg.eigvals(T);costs['ideal_material_eigen_decompositions']+=1
            common=Bc_sqrt@T@Bc_invsqrt
            val=float(np.linalg.svd(common,compute_uv=False)[0]);costs['ideal_common_norm_SVD']+=1
            np.savez_compressed(out/f'IDEAL_{pair}.npz',T=T,B=b,eigenvalues=ev,X=X[pair[0]],Uind=U,K=K,D=D,Mc=metrics[pair[0]],Mm=M)
            ideal_rows.append(dict(pair=pair,spectral_radius=float(np.max(abs(ev))),same_B00_norm=val,
                min_material_gram=float(np.linalg.eigvalsh(b)[0]),material_normal_relative_residual=norm(b@T-numerator)/norm(numerator)))
        ideal_timing=dict(CPU_seconds=time.process_time()-ideal_cpu,wall_seconds=time.perf_counter()-ideal_wall,
                          max_current_normal_relative_residual=max(current_residuals),current_normal_relative_residuals=current_residuals,independent_current_dimensions=U.shape[1])
        _,_,vt=np.linalg.svd(ideal['01']-ideal['11']);costs['direction_difference_SVD']+=1
        h1=vt[0];h1=h1*(1 if h1[np.argmax(abs(h1))]>=0 else -1)
        rng=np.random.default_rng(1401501);h2=rng.normal(size=model.n);h2/=norm(h2)
        directions={'h1':h1,'h2':h2,'hsum':h1+h2}
        np.savez_compressed(out/'DIRECTIONS.npz',**directions,seed=np.array(1401501),hsum_normalized=np.array(False))
        root_responses={};responses={};step_costs=Counter();clone_costs=Counter()
        def run_trial(pair,name,epsilon,sign,h):
            label=f'{pair}__{name}__eps{epsilon:g}__sign{sign:+d}'
            start_cpu=time.process_time();start_wall=time.perf_counter()
            before_clone=dict(base.cost);sol=base.clone(pair=pair,reset_memory=True)
            clone_delta=Counter({k:v-before_clone.get(k,0) for k,v in sol.cost.items() if v-before_clone.get(k,0)})
            clone_costs.update(clone_delta)
            sol.m=m+sign*epsilon*h
            _save_state_npz(out/'states'/f'{label}__input.npz',sol.state_dict())
            before=dict(sol.cost);row=sol.step()
            increment={k:v-before.get(k,0) for k,v in sol.cost.items() if v-before.get(k,0)}
            step_costs.update(increment)
            _save_state_npz(out/'states'/f'{label}__output.npz',dict(state=sol.state_dict(),response=row))
            scalar={k:v for k,v in row.items() if not isinstance(v,np.ndarray) and k not in ['operation_counts']}
            info=dict(label=label,pair=pair,direction=name,epsilon=epsilon,sign=sign,
                CPU_seconds=time.process_time()-start_cpu,wall_seconds=time.perf_counter()-start_wall,
                cost_increment=increment,clone_cost_increment=dict(clone_delta),diagnostics=scalar)
            trials.append(info)
            return sol.m.copy(),sol.J.copy(),info
        for pair in ('00','01','10','11'):
            root_responses[pair]=run_trial(pair,'root',0.,1,np.zeros(model.n))
            for name,h in directions.items():
                for eps in (.01,.003,.001):
                    positive=run_trial(pair,name,eps,1,h);negative=run_trial(pair,name,eps,-1,h)
                    dm=(positive[0]-negative[0])/(2*eps);dj=(positive[1]-negative[1])/(2*eps)
                    predj=np.column_stack([U@X[pair[0]][s]@h for s in range(model.ns)])
                    ec=dj-predj
                    M=metrics[pair[1]];b=B[pair[1]]
                    def readout(z):
                        rhs=sum(ws[s]*np.real(Q[:,s].conj()*(M@L@z[:,s])) for s in range(model.ns))
                        costs['diagnostic_material_error_RHS']+=1
                        costs['diagnostic_material_error_dense_LU']+=1
                        return np.linalg.solve(b,rhs)
                    actual_br=readout(dj);material_current_error=readout(ec)
                    em=dm-actual_br;gap=dm-ideal[pair]@h
                    closure=gap-material_current_error-em
                    np.savez_compressed(out/f'RESPONSE_{pair}_{name}_{eps:g}.npz',h=h,actual_dm=dm,actual_dJ=dj,
                        ideal_dm=ideal[pair]@h,ideal_dJ=predj,current_error=ec,material_best_response_from_actual_current=actual_br,
                        current_error_material_response=material_current_error,material_error=em,idealgap=gap,decomposition_closure=closure)
                    info=dict(row_type='directional',pair=pair,direction=name,epsilon=eps,actual_dm_norm=norm(dm),ideal_dm_norm=norm(ideal[pair]@h),
                        idealgap_norm=norm(gap),current_error_norm=norm(ec),current_error_material_norm=norm(material_current_error),
                        material_error_norm=norm(em),decomposition_absolute_defect=norm(closure),
                        plus_sweeps=positive[2]['diagnostics'].get('aplus_sweeps'),minus_sweeps=negative[2]['diagnostics'].get('aplus_sweeps'),
                        plus_fallback=positive[2]['diagnostics'].get('material_fallback'),minus_fallback=negative[2]['diagnostics'].get('material_fallback'),
                        plus_material_movement=norm(positive[0]-m),minus_material_movement=norm(negative[0]-m),
                        even_material_response_over_epsilon=norm(positive[0]+negative[0]-2*root_responses[pair][0])/(2*eps),
                        bound_active_plus=int(np.count_nonzero((positive[0]<=.02+1e-8)|(positive[0]>=20-1e-8))),
                        bound_active_minus=int(np.count_nonzero((negative[0]<=.02+1e-8)|(negative[0]>=20-1e-8))))
                    tables.append(info);responses[pair,name,eps]=(dm,dj)
        for pair in ('00','01','10','11'):
            for eps in (.01,.003,.001):
                a=responses[pair,'hsum',eps];b=responses[pair,'h1',eps];c=responses[pair,'h2',eps]
                am=a[0]-b[0]-c[0];aj=a[1]-b[1]-c[1]
                np.savez_compressed(out/f'ADDITIVITY_{pair}_{eps:g}.npz',material_additivity_error=am,current_additivity_error=aj)
                tables.append(dict(row_type='additivity',pair=pair,epsilon=eps,material_additivity_absolute=norm(am),
                    material_additivity_relative=norm(am)/max(norm(b[0])+norm(c[0]),1e-300),
                    current_additivity_absolute=norm(aj),current_additivity_relative=norm(aj)/max(norm(b[1])+norm(c[1]),1e-300)))
        root_step=[dict(pair=k,material_absolute_gap=norm(v[0]-m),current_relative_gap=norm(v[1]-J)/norm(J),diagnostics=v[2]['diagnostics']) for k,v in root_responses.items()]
        summary.update(status='COMPLETE',scope='OFFLINE_DIAGNOSTIC_ONLY_FINITE_TM_NO_INVERSE_OR_CONTINUOUS_PDE_CLAIM',
            protocol=dict(side=8,sources=6,receivers=16,anchor='.4+8*(.19<r<.34)',dpcs_rank=16,mf=4,beta_current=1.,beta_material=.2,
                material='R7Arc2sweeps inheritedgtol',prior='zero',noise='zero',directions='h1=toprightsingular(T01-T11);h2seed1401501;hsumUNNORMALIZED',epsilons=[.01,.003,.001]),
            setup=setup,ideal_timing=ideal_timing,ideal_maps=ideal_rows,root_step_gaps=root_step,
            trials=len(trials),directional_rows=36,additivity_rows=12,
            actual_step_counts=dict(step_costs),actual_clone_counts=dict(clone_costs),mechanism_setup_and_ideal_counts=dict(costs),
            physical_counts=dict(model.cost),limitations=['directionalfiniteepsilonresponseisnotassumedFrechetJacobian',
                'decompositionusesrootlinearizedbestresponseanddefinesresidualterms;doesnotcertifyfiniteepsilonlinearizationaccuracy',
                'rootnumericalseedgap,candidaterootswitches,gtolguards,andfinite2sweepscaninfluenceresponse',
                'allstatesrootseededoffline;notdeployableinverseinitialization;noGTloaded'])
    except Exception as exc:
        errors.append(dict(error=repr(exc),traceback=traceback.format_exc()));summary['status']='FAILED';traceback.print_exc()
    finally:
        summary.update(errors=errors,CPU_seconds=time.process_time()-started_cpu,wall_seconds=time.perf_counter()-started_wall,
            process_CPU_seconds_including_imports=time.process_time())
        dump(out/'SUMMARY.json',summary);dump(out/'TRIAL_LOG.json',trials)
        fields=list(dict.fromkeys(k for r in tables for k in r))
        with (out/'TABLES.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(tables)
        (out/'ANALYSIS_NOTES.md').write_text('''# Offline root mechanism audit\n\nOne declared8x8 finite-TM anchor; no ground-truth/test access and no new recovery claims. ROOT_FULL_STATE saves the Fourier/DPCS physical-root embedding. Every actual trial saves complete input/output state and PR/readout diagnostics, including zero-perturbation seed responses. IDEAL arrays use an independent orthonormal current space rather than singular redundant Fourier coordinates.\n\nTABLES directional rows compare central responses with exact-current/exact-material ideal maps and split actual-minus-ideal into current-block and material-block error responses. The closure is an algebraic accounting check using defined residuals, not a certificate that finite epsilon is asymptotic or that the program has a Frechet derivative. ADDITIVITY arrays test the unnormalized h1+h2 response. Interpret plateau/nonadditivity together with readout guard, root selection, finite sweeps, root seed gaps and box events in TRIAL_LOG.\n\nCounts separate one physical forward/setup, ideal normal factors/RHS, actual clone increments and actual-step increments. Dense diagnostic matrix products are offline work; inherited counters retain action/readout detail. Total CPU/wall include this program from main (plus process_CPU including imports), while individual trial times include state serialization. No search or hyperparameter tuning is performed. Parent retains final interpretation.\n''')
        print(json.dumps({k:summary.get(k) for k in ('status','trials','CPU_seconds','wall_seconds','errors')},default=plain))
    return summary['status']=='COMPLETE'


if __name__=='__main__':raise SystemExit(0 if main() else 1)
