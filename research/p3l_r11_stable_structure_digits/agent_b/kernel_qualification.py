"""Actual R11 fixed total-field test on the 24 saved R10 DEV path shapes.
No truth, decoder fitting, or new 576-state solve. All solve attempts are charged.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[k]='1'
import sys,json,time,csv,traceback
from pathlib import Path
import numpy as np
from scipy.linalg import lu_factor,lu_solve
D=Path(__file__).resolve().parent;P=D.parent;sys.path.insert(0,str(P/'shared'))
from core import load_setup,R10,SIGMA,write,plain
from fixed_total_field import TotalFieldSolver,CostLedger,block_krylov,LocalDissipativeBound
RANKS=(32,64,96);METHODS=[f'TOTAL_ANCHOR_R{r}' for r in RANKS]+[f'TOTAL_GREEN_R{r}' for r in RANKS]+['WEIGHTED_ANCHOR_R64','WEIGHTED_REBUILD_R64']

def main():
    ledger=CostLedger(D/'KERNEL_COST_LEDGER.jsonl',30000);rows=[];setups=[];arrays={};status='FAILED';error=None
    ts=time.perf_counter();S=load_setup();C,O,E,Fh,roi=S['C'],S['O'],S['E'],S['Fh'],S['ROI'];loadtime=time.perf_counter()-ts
    try:
        # Both shared-independent shape-independent and parent-anchor fixed total-field spaces are genuine new builds.
        green=TotalFieldSolver(C,O,E,Fh,ranks=RANKS,basis_kind='green',ledger=ledger)
        for pi in range(8):
            parent=f'R10BI_FRESH_{pi:03d}';root=R10/'followup_shape_reuse'
            z=dict(np.load(root/f'{parent}_ANCHOR.npz'));a0=z['fraction'][roi]
            anchor=TotalFieldSolver(C,O,E,Fh,a_anchor=a0,ranks=RANKS,basis_kind='anchor',ledger=ledger)
            t=time.perf_counter();p0=np.sqrt(a0);B0=p0[:,None]*C*p0[None,:];R0=p0[:,None]*E
            V0,meta0=block_krylov(B0,R0,64,ledger=ledger);oldsetup=time.perf_counter()-t
            y=np.load(R10/'results/candidates_fresh10/inputs/online'/f'{parent}.npz')['observation']
            for ti in range(3):
                saved=dict(np.load(root/f'{parent}_t{ti}.npz'));a=saved['fraction'][roi];ds=saved['delta']
                t=time.perf_counter();p=np.sqrt(a);B=p[:,None]*C*p[None,:];L=O*p[None,:];R=p[:,None]*E;weightedform=time.perf_counter()-t
                t=time.perf_counter();Vnew,metanew=block_krylov(B,R,64,ledger=ledger);newbasistime=time.perf_counter()-t
                bound=LocalDissipativeBound(C,O);bm,betatime=bound.shape(a)
                fixedmodels={};projection={}
                for kind,solver in [('ANCHOR',anchor),('GREEN',green)]:
                    models,dt=solver._models_for(a);projection[kind]=dt
                    for r,mm in models.items():fixedmodels[f'TOTAL_{kind}_R{r}']=mm
                weightedmodels={};weightedtimes={}
                for name,V in [('WEIGHTED_ANCHOR_R64',V0),('WEIGHTED_REBUILD_R64',Vnew)]:
                    t=time.perf_counter()
                    if V.shape[1]>=64:
                        ledger.green_columns+=64;BV=B@V;weightedmodels[name]=(V,BV,L@V,V.conj().T@BV,V.conj().T@R)
                    weightedtimes[name]=time.perf_counter()-t
                setups.append(dict(parent=parent,t_index=ti,t=float(saved['t']),
                                   anchor_basis_seconds=anchor.setup['basis_seconds'],green_basis_seconds=green.setup['basis_seconds'],
                                   weighted_anchor_setup_seconds=oldsetup,weighted_rebuild_basis_seconds=newbasistime,
                                   weighted_formation_seconds=weightedform,beta_port_setup_seconds=betatime,
                                   weighted_anchor_rank=V0.shape[1],weighted_rebuild_rank=Vnew.shape[1],
                                   anchor_rank=anchor.V.shape[1],green_rank=green.V.shape[1],
                                   **{k.lower()+'_projection_seconds':v for k,v in projection.items()},
                                   **{k.lower()+'_projection_seconds':v for k,v in weightedtimes.items()}))
                pred={m:[] for m in METHODS};returned={m:[] for m in METHODS};exact=[]
                for ni,delta in enumerate(ds):
                    # Independent exact A reference; no reuse as method fallback.
                    f,_,_,ex=anchor._exact(a,delta,False,'REFERENCE_EXACT_A');exact.append(f)
                    # Rotate arm order over nodes to reduce deterministic cache/order favoring.
                    order=METHODS[(pi+ti+ni)%len(METHODS):]+METHODS[:(pi+ti+ni)%len(METHODS)]
                    for method in order:
                        model=fixedmodels.get(method) if method.startswith('TOTAL') else weightedmodels.get(method)
                        fa=Fh.copy();gate=dict(accepted=False,white_bound=None,potential_bound=None);online=0.;solve=bt=0.;resnorm=None
                        rank=int(method.rsplit('R',1)[1]);available=model is not None
                        if available:
                            V,CV,OV,H,RR=model;ledger.reserve(rank,8,method)
                            t=time.perf_counter();fac=lu_factor(np.eye(rank)-delta*H,check_finite=False);u=lu_solve(fac,RR,check_finite=False)
                            fa=Fh+delta*(OV@u);res=(E-V@u+delta*(CV@u)) if method.startswith('TOTAL') else (R-V@u+delta*(CV@u))
                            solve=time.perf_counter()-t;t=time.perf_counter()
                            if method.startswith('TOTAL'):gate=bound(a,delta,res,fa,y)
                            else:
                                b=float(delta.imag/abs(delta)**2+bm['beta']);rn=float(np.linalg.norm(res));wb=float(np.sqrt(2)/SIGMA*bm['port_norm']*rn/b) if b>0 else None
                                pb=float(np.sqrt(2)/SIGMA*np.linalg.norm(fa-y)*wb+.5*wb*wb) if wb is not None else None
                                gate=dict(accepted=bool(wb is not None and wb<=.01 and pb<=.1),white_bound=wb,potential_bound=pb,b_delta=b,residual_weighted_norm=rn)
                            bt=time.perf_counter()-t;online=solve+bt;resnorm=float(np.linalg.norm(res))
                        fallback=not gate['accepted'];fb=0.
                        if fallback:
                            fr,_,_,fx=anchor._exact(a,delta,False,method+'_INDEPENDENT_EXACT_FALLBACK');fb=fx['exact_seconds']
                        else:fr=fa
                        pred[method].append(fa);returned[method].append(fr)
                        ae=float(np.sqrt(2)/SIGMA*np.linalg.norm(fa-f)) if available else None
                        row=dict(parent=parent,t_index=ti,node=ni,method=method,rank=rank,basis_available=available,
                                 delta_real=float(delta.real),delta_imag=float(delta.imag),gate_pass=bool(gate['accepted']),fallback=fallback,
                                 white_bound=gate['white_bound'],potential_bound=gate['potential_bound'],full_residual_norm=resnorm,
                                 actual_approx_white_error=ae,returned_white_error=float(np.sqrt(2)/SIGMA*np.linalg.norm(fr-f)),
                                 bound_dominates_actual=bool(ae<=gate['white_bound']*(1+1e-8)+1e-9) if available and gate['white_bound'] is not None else None,
                                 solve_seconds=solve,bound_seconds=bt,online_seconds=online,fallback_seconds=fb,total_online_seconds=online+fb,
                                 reference_seconds=ex['exact_seconds'])
                        rows.append(row);ledger.event('NODE',**row)
                np.savez_compressed(D/f'{parent}_t{ti}_RESULT.npz',a=a,delta=ds,F_exact=np.stack(exact),
                                    **{m:np.stack(pred[m]) for m in METHODS},**{m+'_returned':np.stack(returned[m]) for m in METHODS})
        summary={}
        for method in METHODS:
            rr=[r for r in rows if r['method']==method];computed=[r for r in rr if r['basis_available']]
            totalonline=sum(r['total_online_seconds'] for r in rr);pureonline=sum(r['online_seconds'] for r in rr)
            if method.startswith('TOTAL_ANCHOR'):
                basis=sum(s['anchor_basis_seconds'] for s in setups[::3]);shape=sum(s['anchor_projection_seconds']+s['beta_port_setup_seconds'] for s in setups)
            elif method.startswith('TOTAL_GREEN'):
                basis=green.setup['basis_seconds'];shape=sum(s['green_projection_seconds']+s['beta_port_setup_seconds'] for s in setups)
            elif method=='WEIGHTED_ANCHOR_R64':
                basis=sum(s['weighted_anchor_setup_seconds'] for s in setups[::3]);shape=sum(s['weighted_formation_seconds']+s['beta_port_setup_seconds']+s['weighted_anchor_r64_projection_seconds'] for s in setups)
            else:
                basis=0.;shape=sum(s['weighted_formation_seconds']+s['beta_port_setup_seconds']+s['weighted_rebuild_basis_seconds']+s['weighted_rebuild_r64_projection_seconds'] for s in setups)
            summary[method]=dict(nodes=len(rr),computed_nodes=len(computed),qualified=sum(r['gate_pass'] for r in rr),fallback_nodes=sum(r['fallback'] for r in rr),
                max_actual_approx_white_error=max(r['actual_approx_white_error'] for r in computed),
                median_actual_approx_white_error=float(np.median([r['actual_approx_white_error'] for r in computed])),
                max_returned_white_error=max(r['returned_white_error'] for r in rr),sampled_bounds_all_hold=all(r['bound_dominates_actual'] for r in computed),
                fixed_basis_setup_seconds=basis,current_shape_setup_seconds=shape,pure_online_seconds=pureonline,
                online_including_fallback_seconds=totalonline,pure_cold_task_seconds=basis+shape+pureonline,
                cold_task_including_fallback_seconds=basis+shape+totalonline)
        # Reference counted once per matched node, not once per method.
        direct=sum(r['reference_seconds'] for r in rows if r['method']==METHODS[0])
        for v in summary.values():v['saving_cold_vs_exact_A']=1-v['cold_task_including_fallback_seconds']/direct
        write(D/'KERNEL_SUMMARY.json',dict(status='COMPLETE',scope='R10 exposed DEV 24 saved path shapes x8 stored nodes; no GT',
              methods=summary,exact_A=dict(nodes=192,cold_task_seconds=direct),shape_setups=setups,
              cold_accounting='Each total-field rank charged whole build-to96 plus whole current projection-to96, independently. GREEN fixed basis built once for all24shapes; ANCHOR onceperparent. Each method pays full necessary current beta/portnorm. No reference reused as fallback.',
              setup_load_seconds=loadtime,TSOM='NOT_RUN',claim_scope='Fixed total-field Galerkin forward qualification only, not TSOM inversion or continuous/PDE claim'))
        status='COMPLETE'
    except BaseException as e:
        error=repr(e);write(D/'KERNEL_ERROR.json',dict(error=error,traceback=traceback.format_exc()));raise
    finally:
        if rows:
            keys=list(dict.fromkeys(k for row in rows for k in row))
            with (D/'KERNEL_NODE_RESULTS.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
        write(D/'KERNEL_RECEIPT.json',ledger.close(status,error=error))
if __name__=='__main__':main()
