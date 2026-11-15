"""Analytic whole-box invertibility bounds for the actual residual transforms.

Parameter qualification only, not a frozen-identity error-floor experiment.
No truth, forward, optimization or parameter search.
"""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import argparse,json,time,sys
from pathlib import Path


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--legacy-repo',type=Path,required=True)
    a=p.parse_args();wall=time.perf_counter();sys.path.insert(0,str(a.legacy_repo))
    import numpy as np
    from collections import Counter
    from research.research_max_struct_r4.src.cie_onehop import OneHopAbsorption
    r=a.root;rows=[];cost=Counter()
    with np.load(r/'data_public/exposed48/r13_fresh_0017.npz',allow_pickle=False) as z:
        G=z['G'];shape=tuple(z['grid_shape']);lo=float(z['lower']);hi=float(z['upper'])
    blocks=OneHopAbsorption(G,shape,cost)
    for beta in (1.,.2):
        q=blocks.whole_box_bound((lo,hi),beta)
        dmax=np.maximum(abs(1+(beta-blocks.diag_g)*lo),abs(1+(beta-blocks.diag_g)*hi))
        pmin=float(np.min(beta/dmax));bound=q['bound_inf']
        rows.append(dict(physics='scalar_TM_24x24',**q,block_component_dimension=4,
             uniform_P1_singular_value_lower_bound=pmin*(1-bound)/2 if bound<1 else None,
             kernel_scope='Only this frozen inverse grid/kernel and declared material box'))
    with np.load(r/'data_public/vector4/mixed_crossed_bracket.npz',allow_pickle=False) as z:
        G=z['G'];shape=tuple(z['grid_shape']);lo=float(z['lower']);hi=float(z['upper'])
    ijk=np.indices(shape).reshape(3,-1).T;cell_block=ijk//2
    same=np.all(cell_block[:,None,:]==cell_block[None,:,:],axis=-1)
    goff=G.copy();goff[np.diag_indices(len(G))]=0;goff*=np.repeat(np.repeat(same,3,axis=0),3,axis=1)
    diag=np.diag(G);row_sum=np.sum(abs(goff),axis=1)
    for beta in (6.,1.):
        b=beta-diag
        critical=np.divide(-1.,b.real,out=np.full_like(b.real,np.nan),where=b.real!=0)
        cand=np.stack([np.full(len(G),lo),np.full(len(G),hi),critical])
        legal=np.isfinite(cand)&(cand>=lo)&(cand<=hi)
        safe=np.where(legal,cand,lo);ratio=np.where(legal,safe/abs(1+b*safe),-np.inf)
        bound=float(np.max(np.max(ratio,axis=0)*row_sum))
        closest=np.clip(-b.real/abs(b)**2,lo,hi);dmin=float(np.min(abs(1+b*closest)))
        pmin=float(np.min(beta/np.maximum(abs(1+b*lo),abs(1+b*hi))))
        rows.append(dict(physics='vector_CM_RR_6cubed',beta=beta,material_box=[lo,hi],bound_inf=bound,
            diagonal_denominator_min_abs=dmin,sufficient_invertibility_condition=bool(bound<1 and dmin>0),
            block_component_dimension=24,uniform_P1_singular_value_lower_bound=pmin*(1-bound)/np.sqrt(24) if bound<1 else None,
            scope='Actual disjoint 2x2x2 cells,24 coupled xyz rows; shared real isotropic chi'))
    out=r/'audit/METRIC_QUALIFICATION.json';out.write_text(json.dumps(dict(status='COMPLETE_PARAMETER_QUALIFICATION',
         rows=rows,operation_counts=dict(cost),CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-wall,
         new_forward_LU=0,new_complex_RHS=0,truth_access='NONE',
         interpretation='bound<1 is sufficient for no transformed spurious zeros; bound>=1 would be inconclusive, not failure proof'),indent=2)+'\n')
    print(json.dumps(rows))


if __name__=='__main__':main()
