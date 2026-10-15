"""Four explicitly paid states: actual Jacobian and independent solve replay.

This is a numerical physics qualification, not another search campaign.
"""
import json
from types import SimpleNamespace
import numpy as np
from .budget import Book,dump
from .target import OUT,Targets,validate,save_state,pack
from .corrector import linearize,material_real

def main():
    book=Book(OUT);folder=OUT/'actual_physics_qualification';folder.mkdir(parents=True,exist_ok=True)
    with book.job('actual_numerical_physics_qualification_v1'):
        case='p250065_ring_a1_snr20';book.context=dict(phase='actual_physics_qualification',case=case,truth_access=False)
        model,_,y,sigma=Targets(book).case(case)
        with np.load(OUT/'qualification'/case/'qualified_minchange_base.npz') as z:chi=z['chi'].copy()
        base=model.full_state(chi,purpose='ResearchT_actual_Jacobian_qualification_base');save_state(folder/'base.npz',base)
        jac=linearize(model,base,sigma)
        # A legal inward direction even at active passivity/box boundaries.
        direction=(1.75+.25j)-chi;direction/=np.max(abs(direction))
        analytic=jac.jacobian@material_real(direction);fd=[]
        for step in [1e-4,1e-5]:
            endpoint=model.full_state(chi+step*direction,purpose=f'ResearchT_actual_Jacobian_FD_{step:g}')
            numerical=pack(endpoint.data-base.data)*np.sqrt(2)/sigma/step
            fd.append(dict(step=step,relative_error=float(np.linalg.norm(numerical-analytic)/np.linalg.norm(analytic))))
            save_state(folder/f'fd_{step:g}.npz',endpoint)
        # A separately assembled direct NumPy solve, not a free hidden call.
        book.reserve_state(chi,resolution=24,purpose='ResearchT_independent_numpy_full_endpoint')
        book.reserve_rhs(8,'independent_numpy_full_forward')
        with book.span('independent_numpy_full_endpoint',independent_numpy_full_states=1,LU_factorizations=1,full_LU_factorizations=1,full_L_assemblies=1,full_forward_RHS=8,linear_solve_rhs=8):
            matrix=np.eye(model.n)-chi[:,None]*model.G;rhs=chi[:,None]*model.incident
            current=np.linalg.solve(matrix,rhs);field=model.incident+model.G@current;data=model.S@current
            backward=float(np.linalg.norm(matrix@current-rhs)/np.linalg.norm(rhs))
            independent=SimpleNamespace(chi=chi,current=current,field=field,data=data,residual=backward)
        save_state(folder/'independent_numpy.npz',independent)
        matching=dict(current_relative=float(np.linalg.norm(current-base.current)/np.linalg.norm(base.current)),
                      prediction_relative=float(np.linalg.norm(data-base.data)/np.linalg.norm(base.data)),
                      actual_validation=validate(independent,y,sigma))
        # Exact finite-edit injection identity; matrix actions only, no solves.
        finite=[]
        for version in [1,2,3]:
            rows=json.loads((OUT/f'mechanism_v{version}/SUMMARY.json').read_text())
            for row in rows:
                path=OUT/f'mechanism_v{version}'/row['case'];endpoint_path=path/f'edit_{row["edit_id"]:02}'
                with np.load(path/'base.npz') as z:b={k:z[k].copy() for k in z.files}
                actual=endpoint_path/'best_valid.npz'
                if not actual.exists():actual=endpoint_path/'final.npz'
                with np.load(actual) as z:e={k:z[k].copy() for k in z.files}
                d=e['chi']-b['chi'];L=np.eye(model.n)-e['chi'][:,None]*model.G
                injection=d[:,None]*b['field'];defect=L@(e['current']-b['current'])-injection
                finite.append(dict(version=version,case=row['case'],edit_id=row['edit_id'],
                                   relative_defect=float(np.linalg.norm(defect)/max(np.linalg.norm(injection),1e-300))))
        result=dict(status='PASS' if max(r['relative_error'] for r in fd)<1e-3 and matching['prediction_relative']<1e-10 and max(r['relative_defect'] for r in finite)<1e-10 else 'FAIL',
                    finite_difference=fd,independent_solve=matching,finite_secant_identities=finite,
                    scope='actual fixed24grid TM numerical contract, not continuum output enclosure',
                    new_full_states=4,new_RHS=72,probability_claim=None)
        dump(folder/'RESULT.json',result);print(json.dumps(result),flush=True)
    dump(OUT/'RESOURCE_RECEIPT.json',book.receipt())

if __name__=='__main__':main()
