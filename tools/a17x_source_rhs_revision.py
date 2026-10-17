"""One synthetic revision after the real-cache dense-port obstruction. No physics."""
from pathlib import Path
import sys,json,time,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from a17x.finite import finite_state_exchange,galerkin,q_operator
OUT=ROOT/'research/a17x_f3_competition_audit'

def main():
    path=OUT/'checks/SOURCE_RHS_REVISION.json'
    if path.exists():raise RuntimeError('preserve revision evidence')
    start=time.perf_counter();rng=np.random.default_rng(20261019)
    def cm(shape):return rng.normal(size=shape)+1j*rng.normal(size=shape)
    n,p,m=18,3,4;g=cm((n,n));g*=.4/np.linalg.norm(g,2);sensor=cm((m,n))/n
    x0=np.full(n,.4+.2j);l0=np.eye(n)-x0[:,None]*g;field=cm((n,p));P=np.linalg.solve(l0,np.eye(n));H=g@P;O=sensor@P
    basis=np.linalg.qr(cm((n,9)))[0];V=basis[:,:5];C=basis[:,5:8];W=np.c_[V,C];rows=[];arrays={}
    for i in range(12):
        d=.01*cm((n,));l=l0-d[:,None]*g;rhs=d[:,None]*field;kfull=np.linalg.solve(l,rhs)
        kv,dk,*_=finite_state_exchange(l,V,C,rhs);ka=kv+dk;kg=galerkin(l,W)@rhs
        exact=O@q_operator(d,H)@field;reference=sensor@kfull;rho=rhs-l@ka
        alpha=float(np.linalg.svd(l,compute_uv=False)[-1]);bound=float(np.linalg.norm(sensor,2)*np.linalg.norm(rho)/alpha);error=float(np.linalg.norm(sensor@ka-reference))
        rows.append(dict(index=i,physical_DeltaX_rank=int(sum(abs(d)>0)),finite_RHS_columns=p,exact_port_reference_error=float(np.linalg.norm(exact-reference)/np.linalg.norm(reference)),same_space_goal_vs_exchange_error=float(np.linalg.norm(kg-ka)/np.linalg.norm(kg)),direct_state_error=float(np.linalg.norm(sensor@ka-reference)),numeric_residual_output_bound=bound,bound_covers=error<=bound*(1+1e-10),certified=False,proposal_material_changed_by_evaluator=False))
        arrays[f'd_{i:02}']=d;arrays[f'rhs_{i:02}']=rhs;arrays[f'current_full_{i:02}']=kfull;arrays[f'current_reduced_{i:02}']=ka
    receipt=dict(scope='SYNTHETIC_ONLY_ALGORITHM_REVISION',seed=20261019,reason='575-rank cached P0 increments invalidate an unqualified seven-port setup',revision='fixed-material direct finite-state exchange, ordinary goal-ROM mandatory',matrix_dimension=n,examples=len(rows),rows=rows,wall_seconds=time.perf_counter()-start,new_Maxwell_states=0,new_complex_physical_RHS=0,new_physical_LU=0,real_speed_or_discovery_increment='NOT_ESTABLISHED',source_prior='abstract dense passive contrast; not P0 or P1')
    path.write_text(json.dumps(receipt,indent=2)+'\n');np.savez_compressed(OUT/'results/SYNTHETIC_SOURCE_RHS_REVISION.npz',**arrays,G=g,S=sensor,V=V,C=C,field=field)
    with (OUT/'RESOURCE_LEDGER.jsonl').open('a') as f:f.write(json.dumps(dict(event='SYNTHETIC_REVISION_1',wall_seconds=receipt['wall_seconds'],new_full_states=0,new_complex_RHS=0,new_Maxwell_LU=0,examples=12))+'\n')
    print(json.dumps(dict(examples=12,exact_port_max_error=max(x['exact_port_reference_error'] for x in rows),same_basis_max_error=max(x['same_space_goal_vs_exchange_error'] for x in rows),numeric_bounds_cover=all(x['bound_covers'] for x in rows),new_physics=0,wall_seconds=receipt['wall_seconds'])))

if __name__=='__main__':main()
