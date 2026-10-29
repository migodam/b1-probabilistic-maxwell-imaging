"""Paid adjoint-only diagnostics using cached independently validated fields."""
import json
import numpy as np
from accounting import P,Book,write

def main():
    b=Book('POSTHOC_SAVED_ENDPOINT_RANK_REVIEW');rows=[]
    for parent in range(3):
        for method in ('A','B'):
            folder=P/f'results/E1/P{parent}/{method}'
            with np.load(folder/'REFERENCE.npz') as z:G=z['G'].copy();S=z['S'].copy();dep=z['dep'].copy()
            summary=json.loads((folder/'SUMMARY.json').read_text())
            for old in summary['rows']:
                f=folder/old['key']/'ENDPOINT.npz'
                if not f.exists():continue
                with np.load(f) as z:m=z['m'].copy();E=z['E_full'].copy()
                with b.op('CACHED_FIELD_ADJOINT_J_REVIEW',rhs=1,LU=1,adjoint_RHS=1,field_cache_hits=1,new_states=0) as action:
                    L=np.eye(6)-m[:,None]*G
                    T=np.linalg.solve(L.T,S.T).T
                    C=np.einsum('ri,is->sri',T,E).reshape(2,6)
                    J=np.stack([C.real,C.imag],axis=1).reshape(4,6)
                    sv=np.linalg.svd(J,compute_uv=False);sd=np.linalg.svd(J[:,dep],compute_uv=False)
                    b.save(action,m=m,J=J,response_singular_values=sv,dependent_singular_values=sd)
                rows.append(dict(parent=parent,method=method,key=old['key'],original_qualified=old['qualified'],
                    numerical_rank={str(c):int(np.sum(sv>sv[0]*c)) for c in (1e-8,1e-10,1e-12)},response_sigma_min=float(sv[-1]),dependent_condition=float(sd[0]/sd[-1]),
                    minimum_material_boundary_margin=float(min(np.min(m-.5),np.min(4-m))),node_field_squared=float(np.min(np.sum(abs(E)**2,axis=1)))))
    result=dict(status='COMPLETE_POSTHOC_DIAGNOSTIC',points=len(rows),new_adjoint_RHS=b.run_used,rows=rows,
        timing_scope='Not part of original online stopping; paid after-the-fact rank diagnostics, not a retroactive online guarantee',
        continuous_rank='No between-node regularity certificate; narrow interval charts separately certified')
    write(P/'audit/ENDPOINT_RANK_REVIEW.json',result);b.end();print(json.dumps(dict(points=len(rows),paid=b.run_used,min_sigma=min(r['response_sigma_min'] for r in rows),max_dependent_condition=max(r['dependent_condition'] for r in rows))))
if __name__=='__main__':main()
